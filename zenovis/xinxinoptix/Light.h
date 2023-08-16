#pragma once
#include "Sampling.h"
#include "LightTree.h"

#include "TraceStuff.h"
#include "optixPathTracer.h"
#include "zeno/types/LightObject.h"
#include "zxxglslvec.h"

#include "DisneyBRDF.h"
#include "DisneyBSDF.h"

#include "Shape.h"
#include <cuda/random.h>
#include <cuda/helpers.h>
#include <sutil/vec_math.h>

#define _DefaultSkyLightProb_ 0.5f

// static __inline__ __device__
// int GetLightIndex(float p, GenericLight* lightP, int n)
// {
//     int s = 0, e = n-1;
//     while( s < e )
//     {
//         int j = (s+e)/2;
//         float pc = lightP[j].CDF/lightP[n-1].CDF;
//         if(pc<p)
//         {
//             s = j+1;
//         }
//         else
//         {
//             e = j;
//         }
//     }
//     return e;
// }

static __inline__ __device__
vec3 ImportanceSampleEnv(float* env_cdf, int* env_start, int nx, int ny, float p, float &pdf)
{
    if(nx*ny == 0)
    {
        pdf = 1.0f;
        return vec3(0);
    }
    int start = 0; int end = nx*ny-1;
    while(start<end-1)
    {
        int mid = (start + end)/2;
        if(env_cdf[mid]<p)
        {
            start = mid;
        }
        else
        {
            end = mid;
        }
    }
    pdf = 1.0f;
    start = env_start[start];
    int i = start%nx;
    int j = start/nx;
    float theta = ((float)i + 0.5f)/(float) nx * 2.0f * 3.1415926f - 3.1415926f;
    float phi = ((float)j + 0.5f)/(float) ny * 3.1415926f;
    float twoPi2sinTheta = 2.0f * M_PIf * M_PIf * sin(phi);
    //pdf = env_cdf[start + nx*ny] / twoPi2sinTheta;
    vec3 dir = normalize(vec3(cos(theta), sin(phi - 0.5f * 3.1415926f), sin(theta)));
    dir = dir.rotY(to_radians(-params.sky_rot))
             .rotZ(to_radians(-params.sky_rot_z))
             .rotX(to_radians(-params.sky_rot_x))
             .rotY(to_radians(-params.sky_rot_y));
    return dir;
}

static __inline__ __device__ float sampleIES(const float* iesProfile, float h_angle, float v_angle) {

    if (iesProfile == nullptr) { return 0.0f; }

    int h_num = *(int*)iesProfile; ++iesProfile;
    int v_num = *(int*)iesProfile; ++iesProfile;

    const float* h_angles = iesProfile; iesProfile+=h_num;
    const float* v_angles = iesProfile; iesProfile+=v_num;
    const float* intensity = iesProfile;

    auto h_angle_min = h_angles[0];
    auto h_angle_max = h_angles[h_num-1];

    if (h_angle > h_angle_max || h_angle < h_angle_min) { return 0.0f; }

    auto v_angle_min = v_angles[0];
    auto v_angle_max = v_angles[v_num-1];

    if (v_angle > v_angle_max || v_angle < v_angle_min) { return 0.0f; }

    auto lambda = [](float angle, const float* angles, uint num) -> uint { 

        auto start = 0u, end = num-1u;
        auto _idx_ = start;

        while (start<end) {
            _idx_ = (start + end) / 2;

            if(angles[_idx_] > angle) {
                end = _idx_; continue;
            }

            if(angles[_idx_+1] < angle) {
                start = _idx_+1; continue;
            }

            break;
        }
        return _idx_;
    };

    auto v_idx = lambda(v_angle, v_angles, v_num);
    auto h_idx = lambda(h_angle, h_angles, h_num);

    auto _a_ = intensity[h_idx * v_num + v_idx];
    auto _b_ = intensity[h_idx * v_num + v_idx+1];

    auto _c_ = intensity[(h_idx+1) * v_num + v_idx];
    auto _d_ = intensity[(h_idx+1) * v_num + v_idx+1];

    auto v_ratio = (v_angle-v_angles[v_idx]) / (v_angles[v_idx+1]-v_angles[v_idx]);
    auto h_ratio = (h_angle-h_angles[h_idx]) / (h_angles[h_idx+1]-h_angles[h_idx]);

    auto _ab_ = mix(_a_, _b_, v_ratio);
    auto _cd_ = mix(_c_, _d_, v_ratio);

    return mix(_ab_, _cd_, h_ratio);
}

static __inline__ __device__ void sampleSphereIES(LightSampleRecord& lsr, const float2& uu, const float3& shadingP, const float3& center, float radius) {

    float3 vector = center - shadingP;
    float dist2 = dot(vector, vector);
    float dist = sqrtf(dist2);

    if (dist < radius) {
        lsr.PDF = 0.0f;
        return; 
    }

    lsr.dist = dist - radius;
    lsr.dir = vector/dist;
    lsr.n = -lsr.dir;
    lsr.NoL = 1.0f;
    lsr.p = center + radius*lsr.n;
    lsr.PDF = 1.0f;

    lsr.intensity = 1.0f / dist2;
}

namespace detail {
    template <typename T> struct is_void {
        static constexpr bool value = false;
    };
    template <> struct is_void<void> {
        static constexpr bool value = true;
    };
}

template<bool _MIS_, typename TypeEvalBxDF, typename TypeAux = void>
static __inline__ __device__
void DirectLighting(RadiancePRD *prd, RadiancePRD& shadow_prd, const float3& shadingP, const float3& ray_dir, TypeEvalBxDF& evalBxDF, TypeAux* taskAux=nullptr) {

    const float3 wo = normalize(-ray_dir); 
    float3 light_attenuation = vec3(1.0f);

    const float _SKY_PROB_ = params.num_lights>0? _DefaultSkyLightProb_ : 1.0f;

    float scatterPDF = 1.0f;
    float UF = prd->rndf();

    if(UF > _SKY_PROB_) {

        const uint3 idx = optixGetLaunchIndex();

        float lightPickProb = 1.0f - _SKY_PROB_;
        UF = (UF - _SKY_PROB_) / lightPickProb;

        auto lightTree = reinterpret_cast<pbrt::LightTreeSampler*>(params.lightTreeSampler);

        const Vector3f& SP = reinterpret_cast<const Vector3f&>(shadingP);
        const Vector3f& SN = reinterpret_cast<const Vector3f&>(prd->geometryNormal);

        auto pick = lightTree->sample(UF, SP, SN);

        if (pick.prob <= 0.0f) { return; }

        assert(pick.prob >= 0.0f && pick.prob <= 1.0f);

        uint lighIdx = min(pick.lightIdx, params.num_lights-1);
        auto& light = params.lights[lighIdx];

        lightPickProb *= pick.prob;

        LightSampleRecord lsr{};
        float3 emission = light.emission;

        const float* iesProfile = reinterpret_cast<const float*>(light.ies);

        if (nullptr != iesProfile) {
        
            sampleSphereIES(lsr, {}, shadingP, light._cone_.p, _FLT_EPL_);

            if (lsr.PDF <= 0.0f) return; 

            auto v_angle = acos(dot(-lsr.dir, light.N));
            auto h_angle = acos(dot(-lsr.dir, light.T));

            auto intensity = sampleIES(iesProfile, h_angle, v_angle);
            emission *= intensity * lsr.intensity;

        } else if (light.type == zeno::LightType::Direction) {

            bool valid = light.rect.hitAsLight(&lsr, shadingP, -light.N);
            if (!valid) { emission = {}; lsr.PDF = 0.0f; return; }

        } else {

            float2 uu = {prd->rndf(), prd->rndf()};

            switch (light.shape) {
                case zeno::LightShape::Plane:
                    light.rect.sample(&lsr, uu, shadingP); break;
                case zeno::LightShape::Sphere:
                    light.sphere.sample(&lsr, uu, shadingP); break;
                default: break;
            }
        }

        lsr.PDF *= lightPickProb;
        //lsr.p = rtgems::offset_ray(lsr.p, lsr.n);

            if (light.config & zeno::LightConfigDoubleside) {
                lsr.NoL = abs(lsr.NoL);
            }

            if (lsr.NoL > _FLT_EPL_ && lsr.PDF > __FLT_DENORM_MIN__) {

                traceOcclusion(params.handle, shadingP, lsr.dir, 0, lsr.dist-1e-5f, &shadow_prd);
                
                light_attenuation = shadow_prd.shadowAttanuation;

                if (length(light_attenuation) > 0.0f) {
                    
                    auto bxdf_value = evalBxDF(lsr.dir, wo, scatterPDF);
                    float misWeight = 1.0f;

                    if (!light.isDeltaLight()) {
                        misWeight = BRDFBasics::PowerHeuristic(lsr.PDF, scatterPDF);
                    }

                    prd->radiance = light_attenuation * emission * bxdf_value;
                    prd->radiance *= misWeight / lsr.PDF;
                }
            }

    } else {

        float env_weight_sum = 1e-8f;
        int NSamples = prd->depth<=2?1:1;//16 / pow(4.0f, (float)prd->depth-1);
        for(int samples=0;samples<NSamples;samples++) {

            bool hasenv = params.skynx | params.skyny;
            hasenv = params.usingHdrSky && hasenv;
            float envpdf = 1.0f;

            vec3 sunLightDir = hasenv? ImportanceSampleEnv(params.skycdf, params.sky_start,
                                                            params.skynx, params.skyny, rnd(prd->seed), envpdf)
                                    : vec3(params.sunLightDirX, params.sunLightDirY, params.sunLightDirZ);
            auto sun_dir = BRDFBasics::halfPlaneSample(prd->seed, sunLightDir,
                                                    params.sunSoftness * 0.0f); //perturb the sun to have some softness
            sun_dir = hasenv ? normalize(sunLightDir):normalize(sun_dir);

            float tmpPdf;
            auto illum = float3(envSky(sun_dir, sunLightDir, make_float3(0., 0., 1.),
                                        40, // be careful
                                        .45, 15., 1.030725f * 0.3f, params.elapsedTime, tmpPdf));

            auto LP = shadingP;
            auto Ldir = sun_dir;

            if (envpdf < __FLT_DENORM_MIN__) {
                return;
            }

            //LP = rtgems::offset_ray(LP, sun_dir);
            traceOcclusion(params.handle, LP, sun_dir,
                        1e-5f, // tmin
                        1e16f, // tmax,
                        &shadow_prd);

            light_attenuation = shadow_prd.shadowAttanuation;

            auto inverseProb = 1.0f/_SKY_PROB_;
            auto bxdf_value = evalBxDF(sun_dir, wo, scatterPDF, illum);

            vec3 tmp(1.0f);

            if constexpr(_MIS_) {
                float misWeight = BRDFBasics::PowerHeuristic(tmpPdf, scatterPDF);
                misWeight = misWeight>0.0f?misWeight:1.0f;
                misWeight = scatterPDF>1e-5f?misWeight:0.0f;
                misWeight = tmpPdf>1e-5?misWeight:0.0f;

                tmp = (1.0f / NSamples) * misWeight * inverseProb * light_attenuation  / tmpPdf;
            } else {
                tmp = (1.0f / NSamples) * inverseProb * light_attenuation / tmpPdf;
            }

            prd->radiance += (float3)(tmp) * bxdf_value;

            if constexpr (!detail::is_void<TypeAux>::value) {
                if (taskAux != nullptr) {
                    (*taskAux)(tmp);
                }
            }// TypeAux
        }
    }
};
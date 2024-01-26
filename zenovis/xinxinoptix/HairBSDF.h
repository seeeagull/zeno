#pragma once
#include "zxxglslvec.h"
#include "TraceStuff.h"
#include "IOMat.h" 
#include "DisneyBRDF.h"
namespace HairBSDF{
    /*
    * Define as fllow:
    * Hair as a curve r(s), r is a vector from origin to point on curve, s is the length of curve.
    * We have tanget vector t = (dr/ds) / |dr/ds|
    * Normal vector n = (dt/ds) / |dt/ds|
    * Plane p is perpendicular to t
    * Project w_i,w_o to p, get w_i' and w_o'
    * Theta_i is the angle from w_i' to wi [0,pi/2]
    * Theta_o is the angle from w_o' to wo [0,pi/2]
    * Phi_i is the angle from n to w_i' [0,2pi]
    * Phi_o is the angle from n to w_o' [0,2pi]
    * Phi = Phi_i - Phi_o
    * 
    * 
    * 
    * 
    * 
    */

    static __inline__ __device__ float
    sinh(float x)
    {
        float a = exp(2*x) - 1;
        float b = 2 * exp(x);
        return a / b;
    }
    static __inline__ __device__ float
    csch(float x)
    {
        float a = exp(2*x) - 1;
        float b = 2 * exp(x);
        return b / a;
    }
    static __inline__ __device__ float
    I_0(float x)
    {
        float sum = 1.0f;
        float P = 1.0f;
        for(int i=1;i<11;i++){
            P *= x * x / 4 / i / i;
            sum +=P;
        }
        return sum;
    }

    /*From weta, can be repalce */
    static __inline__ __device__ float 
    M_p_Weta(float sinTheta_i, float cosTheta_i, float sinTheta_o, float cosTheta_o, float beta)
    {
        float v = beta * beta;
        float term1 = csch(1/v) / (2*v);
        float term2 = exp(-sinTheta_i * sinTheta_o / v);
        float term3 = cosTheta_i * cosTheta_o / v;
        float term4 = I_0(term3);
        return term1 * term2 * term4;

    }
    static __inline__ __device__ float
    M_p_UE(float sinTheta_i, float sinTheta_o, float beta, float alpha)
    {
        float term1 = 1 / beta / 2.5066283f; // sqrt(2*pi)
        float term2 = sinTheta_i + sinTheta_o - alpha ;
        float term3 = exp(-term2 / 2 / beta/beta);
        return term1 * term3;
    }
    static __inline__ __device__ vec3 
    N_r(float Phi, float h,float ior)
    {
        float term1 = 0.25f * sqrtf((1.0f+cosf(Phi))/2.0f);
        float term2 = BRDFBasics::SchlickDielectic(h,ior);
        return vec3(term1 * term2);
    }
    static __inline__ __device__ vec3
    Ap(float cosTheta_o, float ior, float h,int p, vec3 color)
    {
        float cosTheta = cosTheta_o * sqrtf(1-h*h);
        float f = BRDFBasics::SchlickDielectic(cosTheta,ior);
    }
    static __inline__ __device__ vec3 
    EvaluteHair(float sinTheta_i, float cosTheta_i,
                float sinTheta_o, float cosTheta_o, 
                float Phi, 
                float h, float ior, vec3 sigma_a,
                float beat_m, float beta_n, float alpha)
    {
        vec3 R = M_p_Weta(sinTheta_i,cosTheta_i,sinTheta_o,cosTheta_o,beat_m) * N_r(Phi,h,ior) / abs(cosTheta_o);
        return R;
        //return vec3(M_p_Weta(sinTheta_i,cosTheta_i,sinTheta_o,cosTheta_o,beat_m));
    }
}
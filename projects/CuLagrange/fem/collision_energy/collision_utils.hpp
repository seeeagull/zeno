#pragma once

#include "Structures.hpp"
#include "Utils.hpp"

#include "zensim/math/VecInterface.hpp"


namespace zeno {
namespace COLLISION_UTILS {

    // using namespace std;

    using REAL = float;
    using VECTOR12 = typename zs::vec<REAL,12>;
    using VECTOR4 = typename zs::vec<REAL,4>;
    using VECTOR3 = typename zs::vec<REAL,3>;
    using VECTOR2 = typename zs::vec<REAL,2>;
    using MATRIX3x12 = typename zs::vec<REAL,3,12>;
    using MATRIX12 = typename zs::vec<REAL,12,12>;

    ///////////////////////////////////////////////////////////////////////
    // should we reverse the direction of the force?
    ///////////////////////////////////////////////////////////////////////
    bool reverse(const std::vector<VECTOR3>& v,const std::vector<VECTOR3>& e)
    {
        // get the normal
        VECTOR3 n = e[2].cross(e[0]);
        n = n / n.norm();

        // e[1] is already the collision vertex recentered to the origin
        // (v[0] - v[2])
        const REAL dotted = n.dot(e[1]);
        return (dotted < 0) ? true : false;
    }


    VECTOR12 flatten(const std::vector<VECTOR3>& v) {
        auto res = VECTOR12::zeros();
        for(size_t i = 0;i < 4;++i)
            for(size_t j = 0;j < 3;++j)
                res[i * 3 + j] = v[i][j];
        return res;
    }

    void setCol(MATRIX3x12& m,int col,const VECTOR3& v) {
        for(int i = 0;i < 3;++i)
            m[i][col] = v[i];
    }

    VECTOR3 getCol(MATRIX3x12& m,int col) {
        VECTOR3 res{0};
        for(int i = 0;i < 3;++i)
            res[i] = m[i][col];
        return res;
    }


    ///////////////////////////////////////////////////////////////////////
    // partial of (va - vb)
    ///////////////////////////////////////////////////////////////////////
    MATRIX3x12 vDiffPartial(const VECTOR2& a, const VECTOR2& b)
    {
        auto tPartial = MATRIX3x12::zeros();
        tPartial(0,0) = tPartial(1,1)  = tPartial(2,2) = -a[0];
        tPartial(0,3) = tPartial(1,4)  = tPartial(2,5) = -a[1];
        tPartial(0,6) = tPartial(1,7)  = tPartial(2,8) = b[0];
        tPartial(0,9) = tPartial(1,10) = tPartial(2,11) = b[1];

        return tPartial;
    }

 ///////////////////////////////////////////////////////////////////////
    // gradient of the cross product used to compute the normal,
    // edge-edge case
    ///////////////////////////////////////////////////////////////////////
    MATRIX3x12 crossGradientEE(const std::vector<VECTOR3>& e)
    {
        MATRIX3x12 crossMatrix;

        const REAL e0x = e[0][0];
        const REAL e0y = e[0][1];
        const REAL e0z = e[0][2];

        const REAL e1x = e[1][0];
        const REAL e1y = e[1][1];
        const REAL e1z = e[1][2];

        setCol(crossMatrix,0,VECTOR3(0, -e1z, e1y));
        setCol(crossMatrix,1,VECTOR3(e1z, 0, -e1x));
        setCol(crossMatrix,2,VECTOR3(-e1y, e1x, 0));

        setCol(crossMatrix,3,VECTOR3(0, e1z, -e1y));
        setCol(crossMatrix,4,VECTOR3(-e1z, 0, e1x));
        setCol(crossMatrix,5,VECTOR3(e1y, -e1x, 0));

        setCol(crossMatrix,6,VECTOR3(0, e0z, -e0y));
        setCol(crossMatrix,7,VECTOR3(-e0z, 0, e0x));
        setCol(crossMatrix,8,VECTOR3(e0y, -e0x, 0));

        setCol(crossMatrix,9,VECTOR3(0, -e0z, e0y));
        setCol(crossMatrix,10,VECTOR3(e0z, 0, -e0x));
        setCol(crossMatrix,11,VECTOR3(-e0y, e0x, 0));

        return crossMatrix;
    }


    ///////////////////////////////////////////////////////////////////////
    // gradient of the normal, edge-edge case
    ///////////////////////////////////////////////////////////////////////
    MATRIX3x12 normalGradientEE(const std::vector<VECTOR3>& e)
    {
        VECTOR3 crossed = e[1].cross(e[0]);
        const REAL crossNorm = crossed.norm();
        const REAL crossNormInv = (crossNorm > 1e-8) ? 1.0 / crossed.norm() : 0.0;
        const REAL crossNormCubedInv = (crossNorm > 1e-8) ? 1.0 / pow(crossed.dot(crossed), 1.5) : 0.0;
        MATRIX3x12 crossMatrix = crossGradientEE(e);

        MATRIX3x12 result;
        for (int i = 0; i < 12; i++)
        {
            VECTOR3 crossColumn = getCol(crossMatrix,i);
            auto col_vec = crossNormInv * crossColumn - 
                            ((crossed.dot(crossColumn)) * crossNormCubedInv) * crossed;
            setCol(result,i,col_vec);
        }
        return result;
    }

    ///////////////////////////////////////////////////////////////////////
    // one entry of the rank-3 hessian of the cross product used to compute 
    // the triangle normal, edge-edge case
    ///////////////////////////////////////////////////////////////////////
    VECTOR3 crossHessianEE(const int iIn, const int jIn)
    {
        int i = iIn;
        int j = jIn;

        if (i > j)
        {
            int temp = j;
            j = i;
            i = temp;
        }

        if ((i == 1 && j == 11)  || (i == 2 && j == 7) || (i == 4 && j == 8) || (i == 5 && j == 10))
            return VECTOR3(1, 0, 0);

        if ((i == 0 && j == 8) || (i == 2 && j == 9) || (i == 3 && j == 11) || (i == 5 && j == 6))
            return VECTOR3(0, 1, 0);

        if ((i == 0 && j == 10)  || (i == 1 && j == 6) || (i == 3 && j == 7) || (i == 4 && j == 9))
            return VECTOR3(0, 0, 1);

        if ((i == 1 && j == 8) || (i == 2 && j == 10) || (i == 4 && j == 11) || (i == 5 && j == 7))
            return VECTOR3(-1, 0, 0);

        if ((i == 0 && j == 11) || (i == 2 && j == 6) || (i == 3 && j == 8) || (i == 5 && j == 9))
            return VECTOR3(0, -1, 0);

        if ((i == 0 && j == 7) || (i == 1 && j == 9) || (i == 3 && j == 10) || (i == 4 && j == 6))
            return VECTOR3(0, 0, -1);

        return VECTOR3(0, 0, 0);
    }

    ///////////////////////////////////////////////////////////////////////
    // hessian of the triangle normal, edge-edge case
    ///////////////////////////////////////////////////////////////////////
    std::vector<MATRIX12> normalHessianEE(const std::vector<VECTOR3>& e)
    {
        using namespace std;

        std::vector<MATRIX12> H(3);
        for (int i = 0; i < 3; i++)
            H[i] = MATRIX12::zeros();

        VECTOR3 crossed = e[1].cross(e[0]);
        MATRIX3x12 crossGrad = crossGradientEE(e);
        const VECTOR3& z = crossed;
        
        //denom15 = (z' * z) ^ (1.5);
        REAL denom15 = pow(crossed.dot(crossed), 1.5);
        REAL denom25 = pow(crossed.dot(crossed), 2.5);

        for (int j = 0; j < 12; j++)
            for (int i = 0; i < 12; i++)
            {
                VECTOR3 zGradi = getCol(crossGrad,i);
                VECTOR3 zGradj = getCol(crossGrad,j);
                VECTOR3 zHessianij = crossHessianEE(i,j);

                // z = cross(e2, e0);
                // zGrad = crossGradientVF(:,i);
                // alpha= (z' * zGrad) / (z' * z) ^ (1.5);
                REAL a = z.dot(zGradi) / denom15;

                // final = (zGradj' * zGradi) / denom15 + (z' * cross_hessian(i,j)) / denom15;
                // final = final - 3 * ((z' * zGradi) / denom25) * (zGradj' * z);
                REAL aGrad = (zGradj.dot(zGradi)) / denom15 + 
                            z.dot(crossHessianEE(i,j)) / denom15;
                aGrad -= 3.0 * (z.dot(zGradi) / denom25) * zGradj.dot(z);
                
                //entry = -((zGradj' * z) / denom15) * zGradi + 
                //          1 / norm(z) * zHessianij - 
                //          alpha * zGradj - alphaGradj * z;
                VECTOR3 entry = -((zGradj.dot(z)) / denom15) * zGradi + 
                                    (REAL)1.0 / z.norm() * zHessianij - 
                                    a * zGradj - aGrad * z;

                H[0](i,j) = entry[0];
                H[1](i,j) = entry[1];
                H[2](i,j) = entry[2];
            }
        return H;
    }



    MATRIX3x12 crossGradientVF(const std::vector<VECTOR3>& e)
    {
        MATRIX3x12 crossMatrix;

        const REAL e0x = e[0][0];
        const REAL e0y = e[0][1];
        const REAL e0z = e[0][2];

        const REAL e2x = e[2][0];
        const REAL e2y = e[2][1];
        const REAL e2z = e[2][2];

        setCol(crossMatrix,0,VECTOR3(0,0,0)); 
        setCol(crossMatrix,1,VECTOR3(0,0,0)); 
        setCol(crossMatrix,2,VECTOR3(0,0,0)); 
        setCol(crossMatrix,3,VECTOR3(0, -e0z, e0y)); 
        setCol(crossMatrix,4,VECTOR3(e0z, 0, -e0x)); 
        setCol(crossMatrix,5,VECTOR3(-e0y, e0x, 0));
        setCol(crossMatrix,6,VECTOR3(0, (e0z - e2z), (-e0y + e2y)));
        setCol(crossMatrix,7,VECTOR3((-e0z + e2z), 0, (e0x - e2x)));
        setCol(crossMatrix,8,VECTOR3((e0y - e2y), (-e0x + e2x), 0));
        setCol(crossMatrix,9,VECTOR3(0, e2z, -e2y));
        setCol(crossMatrix,10,VECTOR3(-e2z, 0, e2x));
        setCol(crossMatrix,11,VECTOR3(e2y, -e2x, 0));

        return crossMatrix;
    }

    MATRIX3x12 normalGradientVF(const std::vector<VECTOR3>& e)
    {
        //crossed = cross(e2, e0);
        VECTOR3 crossed = e[2].cross(e[0]);
        REAL crossNorm = crossed.norm();
        const REAL crossNormCubedInv = 1.0 / pow(crossed.dot(crossed), 1.5);
        MATRIX3x12 crossMatrix = crossGradientVF(e);

        //final = zeros(3,12);
        //for i = 1:12
        //  crossColumn = crossMatrix(:,i);
        //  final(:,i) = (1 / crossNorm) * crossColumn - ((crossed' * crossColumn) / crossNormCubed) * crossed;
        //end
        MATRIX3x12 result;
        for (int i = 0; i < 12; i++)
        {
            auto crossColumn = VECTOR3(crossMatrix[0][i],
                crossMatrix[1][i],
                crossMatrix[2][i]);
            VECTOR3 resc = ((REAL)1. / crossNorm) * crossColumn - 
                            ((crossed.dot(crossColumn)) * crossNormCubedInv) * crossed;
            setCol(result,i,resc);
        }
        return result;
    }

    ///////////////////////////////////////////////////////////////////////
    // gradient of spring length, n' * (va - vb)
    ///////////////////////////////////////////////////////////////////////
    VECTOR12 springLengthGradient(const std::vector<VECTOR3>& e,
                                                const VECTOR3& n,
                                                const VECTOR3& diff,
                                                const VECTOR2& a,
                                                const VECTOR2& b){
        MATRIX3x12 nPartial = normalGradientEE(e);
        MATRIX3x12 tPartial = vDiffPartial(a,b);
        const REAL sign = (diff.dot(n) > (REAL)0.0) ? (REAL)-1.0 : (REAL)1.0;
        return sign * nPartial.transpose() * diff + tPartial.transpose() * (sign * n);
    }

    VECTOR12 springLengthGradient(const std::vector<VECTOR3>& v,const std::vector<VECTOR3>& e,const VECTOR3& n)
    {
        const MATRIX3x12 nPartial = normalGradientVF(e);
        const VECTOR3 tvf = v[0] - v[2];

        MATRIX3x12 tvfPartial{0};
        tvfPartial(0,0) = tvfPartial(1,1) = tvfPartial(2,2) = 1.0;
        tvfPartial(0,6) = tvfPartial(1,7) = tvfPartial(2,8) = -1.0;

        //f = nPartial' * (v2 - v0) + tvfPartial' * n;
        return nPartial.transpose() * tvf + tvfPartial.transpose() * n;
    }   
///////////////////////////////////////////////////////////////////////
// one entry of the rank-3 hessian of the cross product used to compute 
// the triangle normal, vertex-face case
///////////////////////////////////////////////////////////////////////
    VECTOR3 crossHessianVF(const int iIn, const int jIn)
    {
        int i = iIn;
        int j = jIn;

        if (i > j)
        {
            int temp = j;
            j = i;
            i = temp;
        }

        if ((i == 5 && j == 7)  || (i == 8 && j == 10) || (i == 4 && j == 11))
            return VECTOR3(1, 0, 0);

        if ((i == 6 && j == 11) || (i == 3 && j == 8) || (i == 5 && j == 9))
            return VECTOR3(0, 1, 0);

        if ((i == 4 && j == 6)  || (i == 7 && j == 9) || (i == 3 && j == 10))
            return VECTOR3(0, 0, 1);

        if ((i == 7 && j == 11) || (i == 4 && j == 8) || (i == 5 && j == 10))
            return VECTOR3(-1, 0, 0);

        if ((i == 5 && j == 6)  || (i == 8 && j == 9) || (i == 3 && j == 11))
            return VECTOR3(0, -1, 0);

        if ((i == 6 && j == 10) || (i == 3 && j == 7) || (i == 4 && j == 9))
            return VECTOR3(0, 0, -1);

        return VECTOR3(0, 0, 0);
    }

///////////////////////////////////////////////////////////////////////
// hessian of the triangle normal, vertex-face case
///////////////////////////////////////////////////////////////////////
    std::vector<MATRIX12> normalHessianVF(const std::vector<VECTOR3>& e)
    {
        using namespace std;

        std::vector<MATRIX12> H(3);
        for (int i = 0; i < 3; i++)
            H[i] = MATRIX12::zeros();

        //crossed = cross(e2, e0);
        //crossNorm = norm(crossed);
        //crossGradient = cross_gradient(x);
        VECTOR3 crossed = e[2].cross(e[0]);
        MATRIX3x12 crossGrad = crossGradientVF(e);
        const VECTOR3& z = crossed;
        
        //denom15 = (z' * z) ^ (1.5);
        REAL denom15 = zs::pow(crossed.dot(crossed), (REAL)1.5);
        REAL denom25 = zs::pow(crossed.dot(crossed), (REAL)2.5);

        for (int j = 0; j < 12; j++)
            for (int i = 0; i < 12; i++)
            {
            auto zGradi = VECTOR3(crossGrad[0][i],crossGrad[1][i],crossGrad[2][i]);
            auto zGradj = VECTOR3(crossGrad[0][j],crossGrad[1][j],crossGrad[2][j]);
            VECTOR3 zHessianij = crossHessianVF(i,j);

            // z = cross(e2, e0);
            // zGrad = crossGradientVF(:,i);
            // alpha= (z' * zGrad) / (z' * z) ^ (1.5);
            REAL a = z.dot(zGradi) / denom15;

            // final = (zGradj' * zGradi) / denom15 + (z' * cross_hessian(i,j)) / denom15;
            // final = final - 3 * ((z' * zGradi) / denom25) * (zGradj' * z);
            REAL aGrad = (zGradj.dot(zGradi)) / denom15 + z.dot(crossHessianVF(i,j)) / denom15;
            aGrad -= (REAL)3.0 * (z.dot(zGradi) / denom25) * zGradj.dot(z);
            
            //entry = -((zGradj' * z) / denom15) * zGradi + 1 / norm(z) * zHessianij - alpha * zGradj - alphaGradj * z;
            VECTOR3 entry = -((zGradj.dot(z)) / denom15) * zGradi + (REAL)1.0 / z.norm() * zHessianij - a * zGradj - aGrad * z;

            H[0](i,j) = entry[0];
            H[1](i,j) = entry[1];
            H[2](i,j) = entry[2];
            }
        return H;
    }


    ///////////////////////////////////////////////////////////////////////
    // hessian of spring length, n' * (v[2] - v[0])
    ///////////////////////////////////////////////////////////////////////
    MATRIX12 springLengthHessian(const std::vector<VECTOR3>& e,
                                                const VECTOR3& n,
                                                const VECTOR3& diff,
                                                const VECTOR2& a,
                                                const VECTOR2& b)
    {
        MATRIX3x12 tPartial = vDiffPartial(a,b);
        const REAL sign = (diff.dot(n) > (REAL)0.0) ? (REAL)-1.0 : (REAL)1.0;

        //% mode-3 contraction
        //[nx ny nz] = normal_hessian(x);
        //final = nx * delta(1) + ny * delta(2) + nz * delta(3);
        std::vector<MATRIX12> normalH = normalHessianEE(e);

        MATRIX12 contracted = diff[0] * normalH[0] + 
                                diff[1] * normalH[1] + 
                                diff[2] * normalH[2];
        contracted *= sign;
        
        //nGrad= normal_gradient(x);
        MATRIX3x12 nGrad = sign * normalGradientEE(e);

        //product = nGrad' * vGrad;
        //final = final + product + product';
        MATRIX12 product = nGrad.transpose() * tPartial;

        return contracted + product + product.transpose();
    }



    ///////////////////////////////////////////////////////////////////////
    // hessian of spring length, n' * (v[0] - v[2])
    ///////////////////////////////////////////////////////////////////////
    MATRIX12 springLengthHessian(const std::vector<VECTOR3>& v,
                                                        const std::vector<VECTOR3>& e,
                                                        const VECTOR3& n){
        const VECTOR3 tvf = v[0] - v[2];

        MATRIX3x12 tvfPartial{0};
        // tvfPartial.setZero();
        tvfPartial(0,0) = tvfPartial(1,1) = tvfPartial(2,2) = 1.0;
        tvfPartial(0,6) = tvfPartial(1,7) = tvfPartial(2,8) = -1.0;

        //% mode-3 contraction
        //[nx ny nz] = normal_hessian(x);
        //final = nx * tvf(1) + ny * tvf(2) + nz * tvf(3);
        const std::vector<MATRIX12> normalH = normalHessianVF(e);
        const MATRIX12 contracted = tvf[0] * normalH[0] + tvf[1] * normalH[1] + 
                                    tvf[2] * normalH[2];
        
        const MATRIX3x12 nGrad = normalGradientVF(e);

        //product = nGrad' * vGrad;
        const MATRIX12 product = nGrad.transpose() * tvfPartial;

        return contracted + product + product.transpose();
    }

    ///////////////////////////////////////////////////////////////////////
    // get the linear interpolation coordinates from v0 to the line segment
    // between v1 and v2
    ///////////////////////////////////////////////////////////////////////
    VECTOR2 getLerp(const VECTOR3 v0, const VECTOR3& v1, const VECTOR3& v2)
    {
        const VECTOR3 e0 = v0 - v1;
        const VECTOR3 e1 = v2 - v1;
        const VECTOR3 e1hat = e1 / e1.norm();
        const REAL projection = e0.dot(e1hat);

        if (projection < 0.0)
            return VECTOR2(1.0, 0.0);

        if (projection >= e1.norm())
            return VECTOR2(0.0, 1.0);

        const REAL ratio = projection / e1.norm();
        return VECTOR2(1.0 - ratio, ratio);
    }


    ///////////////////////////////////////////////////////////////////////
    // find the distance from a line segment (v1, v2) to a point (v0)
    ///////////////////////////////////////////////////////////////////////
    REAL pointLineDistance(const VECTOR3 v0, const VECTOR3& v1, const VECTOR3& v2)
    {
        const VECTOR3 e0 = v0 - v1;
        const VECTOR3 e1 = v2 - v1;
        const VECTOR3 e1hat = e1 / e1.norm();
        const REAL projection = e0.dot(e1hat);

        // if it projects onto the line segment, use that length
        if (projection > 0.0 && projection < e1.norm())
        {
            const VECTOR3 normal = e0 - projection * e1hat;
            return normal.norm();
        }

        // if it doesn't, find the point-point distances
        const REAL diff01 = (v0 - v1).norm();
        const REAL diff02 = (v0 - v2).norm();

        return (diff01 < diff02) ? diff01 : diff02;
    }


    ///////////////////////////////////////////////////////////////////////
    // get the barycentric coordinate of the projection of v[0] onto the triangle
    // formed by v[1], v[2], v[3]
    ///////////////////////////////////////////////////////////////////////
    VECTOR3 getBarycentricCoordinates(const std::vector<VECTOR3>& vertices)
    {
        const VECTOR3 v0 = vertices[1];
        const VECTOR3 v1 = vertices[2];
        const VECTOR3 v2 = vertices[3];
            
        const VECTOR3 e1 = v1 - v0;
        const VECTOR3 e2 = v2 - v0;
        const VECTOR3 n = e1.cross(e2);
        const VECTOR3 nHat = n / n.norm();
        const VECTOR3 v = vertices[0] - (nHat.dot(vertices[0] - v0)) * nHat;

        // get the barycentric coordinates
        const VECTOR3 na = (v2 - v1).cross(v - v1);
        const VECTOR3 nb = (v0 - v2).cross(v - v2);
        const VECTOR3 nc = (v1 - v0).cross(v - v0);
        const VECTOR3 barycentric(n.dot(na) / n.l2NormSqr(),
                                    n.dot(nb) / n.l2NormSqr(),
                                    n.dot(nc) / n.l2NormSqr());

        return barycentric;
    }


    ///////////////////////////////////////////////////////////////////////
    // get the barycentric coordinate of the projection of v[0] onto the triangle
    // formed by v[1], v[2], v[3]
    //
    // but, if the projection is actually outside, project to all of the
    // edges and find the closest point that's still inside the triangle
    ///////////////////////////////////////////////////////////////////////
    VECTOR3 getInsideBarycentricCoordinates(const std::vector<VECTOR3>& vertices)
    {
        VECTOR3 barycentric = getBarycentricCoordinates(vertices);

        // if it's already inside, we're all done
        if (barycentric[0] >= 0.0 &&
            barycentric[1] >= 0.0 &&
            barycentric[2] >= 0.0)
            return barycentric;

        // find distance to all the line segments
        //
        // there's lots of redundant computation between here and getLerp,
        // but let's get it working and see if it fixes the actual
        // artifact before optimizing
        REAL distance12 = pointLineDistance(vertices[0], vertices[1], vertices[2]);
        REAL distance23 = pointLineDistance(vertices[0], vertices[2], vertices[3]);
        REAL distance31 = pointLineDistance(vertices[0], vertices[3], vertices[1]);

        // less than or equal is important here, otherwise fallthrough breaks
        if (distance12 <= distance23 && distance12 <= distance31)
        {
            VECTOR2 lerp = getLerp(vertices[0], vertices[1], vertices[2]);
            barycentric[0] = lerp[0];
            barycentric[1] = lerp[1];
            barycentric[2] = 0.0;
            return barycentric;
        }
        
        // less than or equal is important here, otherwise fallthrough breaks
        if (distance23 <= distance12 && distance23 <= distance31)
        {
            VECTOR2 lerp = getLerp(vertices[0], vertices[2], vertices[3]);
            barycentric[0] = 0.0;
            barycentric[1] = lerp[0];
            barycentric[2] = lerp[1];
            return barycentric;
        }

        // else it must be the 31 case
        VECTOR2 lerp = getLerp(vertices[0], vertices[3], vertices[1]);
        barycentric[0] = lerp[1];
        barycentric[1] = 0.0;
        barycentric[2] = lerp[0];
        return barycentric;
    }

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    MATRIX3x12 tDiffPartial(const VECTOR3& bary)
    {
        MATRIX3x12 tPartial{0};
        tPartial(0,0) = tPartial(1,1)  = tPartial(2,2) = 1.0;
        tPartial(0,3) = tPartial(1,4)  = tPartial(2,5) = -bary[0];
        tPartial(0,6) = tPartial(1,7)  = tPartial(2,8) = -bary[1];
        tPartial(0,9) = tPartial(1,10) = tPartial(2,11) = -bary[2];

        return tPartial;
    }

    ///////////////////////////////////////////////////////////////////////
    // are the two edges nearly parallel?
    ///////////////////////////////////////////////////////////////////////
    bool nearlyParallel(const std::vector<VECTOR3> e){
        const VECTOR3 e0 = e[0].normalized();
        const VECTOR3 e1 = e[1].normalized();
        const REAL dotted = zs::abs(e0.dot(e1));

        // too conservative, still seeing some conditioning problems
        // in the simulation. If the mesh suddenly pops, it means
        // that the conditioning problem made the solve go haywire.
        //const REAL eps = 1e-4;
        
        // this is still quite conservative, with some popping visible
        // in the simulation
        //const REAL eps = 1e-3;
        
        // this seems too permissive, and ends up missing some collisions,
        // but is what we're using for now
        const REAL eps = 1e-2;

        return (dotted > (REAL)1.0 - eps);
    }

   


    ///////////////////////////////////////////////////////////////////////
    // does this face and edge intersect?
    ///////////////////////////////////////////////////////////////////////
    bool faceEdgeIntersection(const std::vector<VECTOR3>& triangleVertices, 
                            const std::vector<VECTOR3>& edgeVertices)
    {
        assert(triangleVertices.size() == 3);
        assert(edgeVertices.size() == 2);

        const VECTOR3& a = triangleVertices[0];
        const VECTOR3& b = triangleVertices[1];
        const VECTOR3& c = triangleVertices[2];

        const VECTOR3& origin = edgeVertices[0];
        const VECTOR3& edgeDiff = (edgeVertices[1] - edgeVertices[0]);
        const VECTOR3& direction = edgeDiff.normalized();

        const VECTOR3 geometricNormal = ((b - a).cross(c - a)).normalized();

        const VECTOR3 diff = a - origin;
        REAL denom = direction.dot(geometricNormal);
        if (fabs(denom) <= 0.0) return false;

        REAL t = diff.dot(geometricNormal) / denom;
        if (t < 0) return false;

        VECTOR3 h = origin + direction * t;

        VECTOR3 test = (b - a).cross(h - a);
        if (geometricNormal.dot(test) < 0) return false; 
        test = (c - b).cross(h - b);
        if (geometricNormal.dot(test) < 0) return false; 
        test = (a - c).cross(h - c);
        if (geometricNormal.dot(test) < 0) return false; 

        if (t < edgeDiff.norm())
            return true;

        return false;
    }
};
};
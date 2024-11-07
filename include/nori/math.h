#pragma once

#include <nori/common.h>
#include <nori/transform.h>
#include <nori/vector.h>
#include <nori/color.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class Math 
{
    public:

    static constexpr float EPSILON = 1e-6;

    static constexpr float PI = M_PI;
    static constexpr float ONE_MINUS_EPSILON = 1 - EPSILON;
    //static constexpr float INV_PI = 1 / M_PI;

    /// Compute the dot product of two vectors
    inline static float dot(const Vector3f &a, const Vector3f &b) { return a.dot(b); }
    inline static float dot(const Normal3f &a, const Normal3f &b) { return a.dot(b); }
    inline static float dot(const Vector3f &a, const Normal3f &b) { return a.dot(b); }
    inline static float dot(const Normal3f &a, const Vector3f &b) { return a.dot(b); }

    inline static float clamp(float v, float min, float max) { return std::max(min, std::min(max, v)); }
    inline static Color3f clamp(const Color3f &v, float min, float max) { return Color3f(clamp(v.x(), min, max), clamp(v.y(), min, max), clamp(v.z(), min, max)); }

    inline static float toRadians(float degrees) { return degrees * PI / 180.0f; }
    inline static float toDegrees(float radians) { return radians * 180.0f / PI; }

    /// Compute the cross product of two vectors
    inline static Vector3f cross(const Vector3f &a, const Vector3f &b) { return a.cross(b); }

    /// Compute the normalization of a vector
    inline static Vector3f normalize(const Vector3f &a) { return a.normalized(); }

    /// Compute the squared norm of a vector
    inline static float squaredNorm(const Vector3f &a) { return a.squaredNorm(); }

    /// Compute the norm of a vector
    inline static float norm(const Vector3f &a) { return a.norm(); }

    /// Compute the absolute value of a float
    inline static float abs(float a) { return std::abs(a); }
    inline static Color3f abs(const Color3f &a) { return Color3f(std::abs(a.x()), std::abs(a.y()), std::abs(a.z())); }
    inline static float absDot(const Vector3f &a, const Vector3f &b) { return std::abs(dot(a, b)); }

    /// Compute the square root of a float
    inline static float sqrt(float a) { return std::sqrt(a); }

    inline static float exp(float a) { return std::exp(a); }
    inline static Color3f sqrt(const Color3f &a) { return Color3f(std::sqrt(a.x()), std::sqrt(a.y()), std::sqrt(a.z())); }
    inline static Color3f exp(const Color3f &a) { return Color3f(std::exp(a.x()), std::exp(a.y()), std::exp(a.z())); }

    inline static float safeSqrt(float a) { return std::sqrt(std::max(a, 0.0f)); }
    
    inline static float pow(float a, float b) { return std::pow(a, b); }
    inline static float pow2(float a) { return a * a; }
    inline static float pow3(float a) { return a * a * a; }
    inline static float pow4(float a) { float a2 = a*a; return a2*a2; }
    inline static float pow5(float a) { float a2 = a*a; return a2*a2*a; }

    inline static Color3f pow(const Color3f &a, float b) { return Color3f(std::pow(a.x(), b), std::pow(a.y(), b), std::pow(a.z(), b)); }
    inline static Color3f pow2(const Color3f &a) { return Color3f(a.x() * a.x(), a.y() * a.y(), a.z() * a.z()); }
    inline static Color3f pow3(const Color3f &a) { return Color3f(a.x() * a.x() * a.x(), a.y() * a.y() * a.y(), a.z() * a.z() * a.z()); }
    inline static Color3f pow4(const Color3f &a) { return Color3f(a.x() * a.x() * a.x() * a.x(), a.y() * a.y() * a.y() * a.y(), a.z() * a.z() * a.z() * a.z()); }

    /// Compute the cosine of an angle
    inline static float cos(float a) { return std::cos(a); }
    inline static float cos(const Vector3f &a, const Vector3f &b) { return a.dot(b) / (a.norm() * b.norm()); }
    inline static float absCos(const Vector3f &a, const Vector3f &b) { return std::abs(cos(a, b)); }

    /// Compute the sine of an angle
    inline static float sin(float a) { return std::sin(a); }

    inline static float cosTheta(const Vector3f &w) { return w.z(); }
    inline static float cos2Theta(const Vector3f &w) { return w.z() * w.z(); }
    inline static float absCosTheta(const Vector3f &w) { return std::abs(w.z()); }

    inline static float max(float a, float b) { return std::max(a, b); }
    inline static float min(float a, float b) { return std::min(a, b); }
    
    inline static Point2f max(const Point2f &a, const Point2f &b) { return Point2f(std::max(a.x(), b.x()), std::max(a.y(), b.y())); }
    inline static Point2f min(const Point2f &a, const Point2f &b) { return Point2f(std::min(a.x(), b.x()), std::min(a.y(), b.y())); }

    inline static float asin(float a) { return std::asin(a); }

    inline static Color3f max (const Color3f &a, const Color3f &b) { return Color3f(std::max(a.x(), b.x()), std::max(a.y(), b.y()), std::max(a.z(), b.z())); }
    inline static Color3f min (const Color3f &a, const Color3f &b) { return Color3f(std::min(a.x(), b.x()), std::min(a.y(), b.y()), std::min(a.z(), b.z())); }

    inline static Color3f max (const Color3f &a, float b) { return Color3f(std::max(a.x(), b), std::max(a.y(), b), std::max(a.z(), b)); }
    inline static Color3f min (const Color3f &a, float b) { return Color3f(std::min(a.x(), b), std::min(a.y(), b), std::min(a.z(), b)); }

    inline static float lerp(float a, float b, float t) { return a + (b - a) * t; }
    inline static Color3f lerp(const Color3f &a, const Color3f &b, float t1, float t2) { return a * t1 + b * t2; } 


    // Compute sin of a vector in local coordinates
    inline static float sin2Theta(const Vector3f &w) {
        return std::max((float)0, (float)1 - cos2Theta(w));
    }
    inline static float sinTheta(const Vector3f &w) {
        return std::sqrt(sin2Theta(w));
    }

    inline static float cosPhi(const Vector3f &w) {
        float sinTh = sinTheta(w);
        return (sinTh == 0) ? 1 : clamp(w.x() / sinTh, -1.0f, 1.0f);
    }
    inline static float sinPhi(const Vector3f &w) {
        float sinTh = sinTheta(w);
        return (sinTh == 0) ? 0 : clamp(w.y() / sinTh, -1.0f, 1.0f);
    }

    inline static float cos2Phi(const Vector3f &w) {
        return cosPhi(w) * cosPhi(w);
    }
    inline static float sin2Phi(const Vector3f &w) {
        return sinPhi(w) * sinPhi(w);
    }

    inline static float tanTheta(const Vector3f &w) 
    {
        return sinTheta(w) / cosTheta(w);
    }

    inline static float tanTheta2(const Vector3f &w) 
    {
        return sin2Theta(w) / cos2Theta(w);
    }

    inline static float floor(float a) { return std::floor(a); }
    inline static float ceil(float a) { return std::ceil(a); }

    inline static bool sameSign(float a, float b) { return a * b > 0; }
    inline static int sign(float a) { return (a > 0) - (a < 0); }

    inline static bool sameDirection(const Vector3f &a, const Vector3f &b) { return dot(a, b) > 0; }
    inline static bool sameHemisphere(const Vector3f &a, const Vector3f &b) { return dot(a, b) > 0; }

    inline static Vector3f reflect(const Vector3f &wo, const Vector3f &n) { return -wo + 2 * dot(wo, n) * n; }
    inline static float distanceSquared(const Point3f &a, const Point3f &b) { return (a - b).squaredNorm(); }
    inline static float distance(const Point3f &a, const Point3f &b) { return (a - b).norm(); }

    // Round up to the next power of two
    inline static int roundUpPow2 (int x) 
    {
        x--;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x + 1;
    }

    inline static bool isInf(float v) 
    {
        return std::isinf(v);
    }

    inline static uint32_t floatToBits(float v) 
    {
        uint32_t ui;
        memcpy(&ui, &v, sizeof(float));
        return ui;
    }

    inline static float bitsToFloat(uint32_t ui) 
    {
        float v;
        memcpy(&v, &ui, sizeof(uint32_t));
        return v;
    }

    // Compute the next representable floating point value
    inline static float nextFloatUp(float v) 
    {
        // Infinity and zero
        if (isInf(v) && v > 0.f)
            return v;

        if (v == -0.f)
            v = 0.f;

        // Next representable value
        uint32_t ui = floatToBits(v);

        if (v >= 0) 
            ++ui;
        else    
            --ui;

        return bitsToFloat(ui);
    }

    // Compute the previous representable floating point value
    inline static float nextFloatDown(float v) 
    {
        // Infinity and zero
        if (isInf(v) && v < 0.f)
            return v;

        if (v == 0.f)
            v = -0.f;

        // Next representable value
        uint32_t ui = floatToBits(v);

        if (v > 0) 
            --ui;
        else    
            ++ui;

        return bitsToFloat(ui);
    }

    inline static Vector3f faceForward(const Vector3f &n, const Vector3f &v) 
    {
        return (dot(n, v) < 0.f) ? -n : n;
    }

    inline static float balanceHeuristic(float fPdf, float gPdf) 
    {
        return fPdf / (fPdf + gPdf);
    }

    inline static float balanceHeuristic(int nf, float fPdf, int ng, float gPdf) 
    {
        return (nf * fPdf) / (nf * fPdf + ng * gPdf);
    }

    inline static float powerHeuristic(int nf, float fPdf, int ng, float gPdf) 
    {
        float f = nf * fPdf, g = ng * gPdf;
        float term1 = (f * f);
        float term2 = (f * f + g * g);
        if (term2 > 1e-6)
            return term1 / term2;
        else
            return 0;
    }

    inline static bool isNaN(float v) 
    {
        return std::isnan(v);
    }

    // Computes the 3D gradient of slice 1,
    //  x and y are slice1 gradients and 
    //  and Z is the gradient going to slice 2 
    static Eigen::MatrixXf finiteDifferenceGradientModulus3D(const Eigen::MatrixXf& slice1, const Eigen::MatrixXf& slice2);
    
    // Function to compute the Sobel gradient modulus in 3D (two slices)
    static Eigen::MatrixXf sobelGradientModulus3D(const Eigen::MatrixXf& slice1, const Eigen::MatrixXf& slice2);

    static float matrixVariance(const Eigen::MatrixXf& block, float &mean);
    static float matrixStdDeviation(const Eigen::MatrixXf& block, float &mean);
    static float blockPercentualStdDeviation(const Eigen::MatrixXf& block);

    static Eigen::MatrixXf computeVariance3x3(const Eigen::MatrixXf& slice1, const Eigen::MatrixXf& slice2);
    static Eigen::MatrixXf luminanceMatrix(const Bitmap *image);

    static Eigen::MatrixXf equalize(const Eigen::MatrixXf &m);

    // Finds a root for the function f using the secant method.
    //  Input: Function f, initial x0
    static float findRoot(std::function<float(float)> f, float x0);

    static Color3f fresnel(float cosThetaI, const Color3f &extIOR, const Color3f &intIOR)
    {
        Color3f eta = extIOR / intIOR;

        float cosThetaI2 = cosThetaI * cosThetaI;
        float sinThetaI2 = std::max(0.0f, 1.0f - cosThetaI2);
        float sinThetaI4 = sinThetaI2 * sinThetaI2;

        // Calculate a^2 + b^2 term
        Color3f eta2 = eta * eta;
        Color3f term1 = eta2 - Color3f(sinThetaI2, sinThetaI2, sinThetaI2);
        Color3f a2b2 = (term1 * term1 + Color3f(4.0f, 4.0f, 4.0f) * (eta2)).sqrt();

        // r_perpendicular (r_perp)
        Color3f twoACosTheta = eta * (2.0f * cosThetaI);
        Color3f rPerpNumerator = a2b2 - twoACosTheta + Color3f(cosThetaI2, cosThetaI2, cosThetaI2);
        Color3f rPerpDenominator = a2b2 + twoACosTheta + Color3f(cosThetaI2, cosThetaI2, cosThetaI2);
        Color3f rPerp = rPerpNumerator / rPerpDenominator;

        // r_parallel (r_parallel)
        Color3f rParallelNumerator = Color3f(cosThetaI2, cosThetaI2, cosThetaI2) * a2b2 - twoACosTheta + Color3f(sinThetaI4, sinThetaI4, sinThetaI4);
        Color3f rParallelDenominator = Color3f(cosThetaI2, cosThetaI2, cosThetaI2) * a2b2 + twoACosTheta + Color3f(sinThetaI4, sinThetaI4, sinThetaI4);
        Color3f rParallel = rPerp * (rParallelNumerator / rParallelDenominator);

        // Average the perpendicular and parallel reflectance
        return (rPerp + rParallel) * 0.5f;
    }

    static Color3f fresnel(float cosThetaI, const Color3f &extIOR, const Color3f &intIOR, const Color3f &K)
    {
        Color3f eta = extIOR / intIOR;

        float cosThetaI2 = cosThetaI * cosThetaI;
        float sinThetaI2 = std::max(0.0f, 1.0f - cosThetaI2);
        float sinThetaI4 = sinThetaI2 * sinThetaI2;

        // Calculate a^2 + b^2 term
        Color3f eta2 = eta * eta;
        Color3f k2 = K * K;
        Color3f term1 = eta2 - k2 - Color3f(sinThetaI2, sinThetaI2, sinThetaI2);
        Color3f a2b2 = (term1 * term1 + Color3f(4.0f, 4.0f, 4.0f) * (eta2 * k2)).sqrt();

        // r_perpendicular (r_perp)
        Color3f twoACosTheta = eta * (2.0f * cosThetaI);
        Color3f rPerpNumerator = a2b2 - twoACosTheta + Color3f(cosThetaI2, cosThetaI2, cosThetaI2);
        Color3f rPerpDenominator = a2b2 + twoACosTheta + Color3f(cosThetaI2, cosThetaI2, cosThetaI2);
        Color3f rPerp = rPerpNumerator / rPerpDenominator;

        // r_parallel (r_parallel)
        Color3f rParallelNumerator = Color3f(cosThetaI2, cosThetaI2, cosThetaI2) * a2b2 - twoACosTheta + Color3f(sinThetaI4, sinThetaI4, sinThetaI4);
        Color3f rParallelDenominator = Color3f(cosThetaI2, cosThetaI2, cosThetaI2) * a2b2 + twoACosTheta + Color3f(sinThetaI4, sinThetaI4, sinThetaI4);
        Color3f rParallel = rPerp * (rParallelNumerator / rParallelDenominator);

        // Average the perpendicular and parallel reflectance
        return (rPerp + rParallel) * 0.5f;
    }

    static Vector3f refract (Vector3f wi, Vector3f n, float etaI, float etaT)
    {
        float eta = etaI / etaT;

        float cosThetaI = dot(wi, n);
        float sin2ThetaI = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
        float sin2ThetaT = eta * eta * sin2ThetaI;

        if (sin2ThetaT >= 1.0f)
            return Vector3f(0.0f, 0.0f, 0.0f);

        float cosThetaT = std::sqrt(1.0f - sin2ThetaT);
        return eta * (-wi) + (eta * cosThetaI - cosThetaT) * n;
    }

    static Eigen::MatrixXf weibull_stretched_exponential(const Eigen::MatrixXf &m, float alpha=0.5, float beta=0.2);
    inline static float weibull_stretched_exponential(float x, float alpha=0.5, float beta=0.2)
    {
        return 1 - std::exp(-std::pow(x/alpha, beta));   
    }
};

NORI_NAMESPACE_END
/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN


Vector3f Warp::to_cartesian(float r, float theta, float phi) {
    return Vector3f(r * std::sin(theta) * std::cos(phi),
                     r * std::sin(theta) * std::sin(phi),
                     r * std::cos(theta));
}

Point2f Warp::to_polar (const Vector3f &v) {
    float r = v.norm();
    float theta = std::acos(v.z() / r);
    float phi = std::atan2(v.y(), v.x());
    return Point2f(theta, phi);
}

int Warp::sampleDiscrete(const std::vector<float> &weights, float sample, float &pdf,
                   float &sampleRemapped) 
{
    if (weights.empty()) 
    {
        pdf = 0;
        return -1;
    }

    // Sum of weights
    float sumWeights = 0;
    for (float w : weights)
        sumWeights += w;

    // Remap the sample to [0, sumWeights)
    float up = sample * sumWeights;
    if (up == sumWeights)
        up = Math::nextFloatDown(up);

    // Find the offset of the sample
    int offset = 0;
    float sum = 0;
    while (sum + weights[offset] <= up)
        sum += weights[offset++];
  
    // Compute the PDF and remapped sample
    pdf = weights[offset] / sumWeights;
    sampleRemapped = std::min((up - sum) / weights[offset], 1 - Math::EPSILON);

    return offset;
}

int Warp::sampleDiscrete(float weights[],
                    int nWeights,
                    float sample, float &pdf,
                    float &sampleRemapped) 
{
    if (nWeights == 0) 
    {
        pdf = 0;
        return -1;
    }

    // Sum of weights
    float sumWeights = 0;
    for (int i = 0; i < nWeights; ++i)
        sumWeights += weights[i];

    // Remap the sample to [0, sumWeights)
    float up = sample * sumWeights;
    if (up == sumWeights)
        up = Math::nextFloatDown(up);

    // Find the offset of the sample
    int offset = 0;
    float sum = 0;
    while (sum + weights[offset] <= up)
        sum += weights[offset++];
  
    // Compute the PDF and remapped sample
    pdf = weights[offset] / sumWeights;
    sampleRemapped = std::min((up - sum) / weights[offset], 1 - Math::EPSILON);

    return offset;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformTriangle(const Point2f &sample) 
{
    float x = std::sqrt(sample.x());
    return Point2f(1 - x, x * sample.y());
}

float Warp::squareToUniformTrianglePdf(const Point2f &p) 
{
    return (p.x() >= 0.0f && p.x() <= 1.0f && p.y() >= 0.0f && p.y() <= 1.0f) ? 2.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) 
{
    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f &p) 
{
    throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) 
{
    float r = std::sqrt(sample.x());
    float theta = 2 * M_PI * sample.y();
    return Point2f(r * std::cos(theta), r * std::sin(theta));
}

Point2f Warp::SrToDisk(float sample, float r) 
{
    float theta = 2 * M_PI * sample;
    return Point2f(r * std::cos(theta), r * std::sin(theta));
}

float Warp::SrToDiskPdf(const Point2f &p) 
{
    float thetaPdf = 1.0f / (2 * M_PI);
    return thetaPdf;
}

float Warp::squareToUniformDiskPdf(const Point2f &p) 
{
    return (p.norm() <= 1.0f) ? 1.0f / M_PI : 0.0f;
}


Point2f Warp::concentricSampleDisk(const Point2f &sample) 
{
    Point2f p;
    float r = 0, theta=0;
    float sx = 2 * sample.x() - 1;
    float sy = 2 * sample.y() - 1;
    if (sx == 0.0f && sy == 0.0f) {
        p = Point2f(0.0f, 0.0f);
    } else if (sx >= -sy) {
        if (sx > sy) {
            r = sx;
            if (sy > 0.0f) {
                theta = sy / r;
            } else {
                theta = 8 + sy / r;
            }
        } else {
            r = sy;
            theta = 2 - sx / r;
        }
    } else {
        if (sx <= sy) {
            r = -sx;
            theta = 4 - sy / r;
        } else {
            r = -sy;
            theta = 6 + sx / r;
        }
    }

    theta *= M_PI / 4;
    p.x() = r * std::cos(theta);
    p.y() = r * std::sin(theta);
    
    return p;
}

float Warp::concentricSampleDiskPdf(const Point2f &v)
{
    return (v.norm() <= 1.0f) ? 1.0f / M_PI : 0.0f;
} 

Vector3f Warp::squareToUniformSphere(const Point2f &sample) 
{
    float z = 1 - 2 * sample.x();
    float r = std::sqrt(std::max(0.0f, 1 - z * z));
    float phi = 2 * M_PI * sample.y();
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) 
{
    return (v.squaredNorm() <= 1.0f) ? 1.0f / (4 * M_PI) : 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) 
{
    float z = sample[0];
    float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
    float phi = 2 * M_PI * sample[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) 
{
    return 1 / (2 * M_PI);
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) 
{
    Point2f p = squareToUniformDisk(sample);
    float z = std::sqrt(std::max(0.0f, 1 - p.x() * p.x() - p.y() * p.y()));
    return Vector3f(p.x(), p.y(), z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) 
{
    return (v.z() >= 0.0f) ? v.z() / M_PI : 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    // Step 1: Sample azimuthal angle phi
    float phi = 2.0f * M_PI * sample.x();

    // Step 2: Sample elevation angle theta using the inverse CDF
    float tanTheta2 = -alpha * alpha * std::log(1.0f - sample.y());
    float theta = std::atan(std::sqrt(tanTheta2));

    // Step 3: Convert spherical to Cartesian coordinates
    float sinTheta = std::sin(theta);
    float cosTheta = std::cos(theta);

    return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}


float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    // Extract theta from the z-component of the normalized vector
    float cosTheta = m.z();
    if (cosTheta <= 0.0f)
        return 0.0f;

    float tanTheta2 = (1.0f - cosTheta * cosTheta) / (cosTheta * cosTheta);
    float cosTheta3 = cosTheta * cosTheta * cosTheta;

    // Beckmann distribution PDF
    return std::exp(-tanTheta2 / (alpha * alpha)) / (M_PI * alpha * alpha * cosTheta3);
}

float Warp::squareToSrDecay(const float &sample, float sigmaT)
{
    return (-std::log(1 - sample) / sigmaT);
}

float Warp::squareToSrDecayPdf(const float &sample, float sigmaT)
{
    return sigmaT * std::exp(-sigmaT * sample);
}


NORI_NAMESPACE_END

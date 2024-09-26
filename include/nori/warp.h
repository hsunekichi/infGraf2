/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/common.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

/// A collection of useful warping functions for importance sampling
class Warp {
public:
    // Convert between polar and cartesian
    static Vector3f to_cartesian(float r, float theta, float phi);
    static Point2f to_polar(const Vector3f &v);

    /// Dummy warping function: takes uniformly distributed points in a square and just returns them
    static Point2f squareToUniformSquare(const Point2f &sample);

    /// Probability density of \ref squareToUniformSquare()
    static float squareToUniformSquarePdf(const Point2f &p);

    /// Sample a 2D tent distribution
    static Point2f squareToTent(const Point2f &sample);

    /// Probability density of \ref squareToTent()
    static float squareToTentPdf(const Point2f &p);

    /// Uniformly sample a vector on a 2D isosceles right triangle with area 1/2 based on its barycentric
    static Point2f squareToUniformTriangle(const Point2f& sample);
    
    /// Probability density of \ref squareToUniformTriangle()
    static float squareToUniformTrianglePdf(const Point2f& p);

    /// Uniformly sample a vector on a 2D disk with radius 1, centered around the origin
    static Point2f squareToUniformDisk(const Point2f &sample);

    /// Probability density of \ref squareToUniformDisk()
    static float squareToUniformDiskPdf(const Point2f &p);

    /// Sample a point on a unit disk with a distribution that is uniform with respect to area
    static Point2f concentricSampleDisk(const Point2f &sample);

    /// Probability density of \ref concentricSampleDisk()
    static float concentricSampleDiskPdf(const Point2f &v);

    /// Uniformly sample a vector on the unit sphere with respect to solid angles
    static Vector3f squareToUniformSphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformSphere()
    static float squareToUniformSpherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to solid angles
    static Vector3f squareToUniformHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformHemisphere()
    static float squareToUniformHemispherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to projected solid angles
    static Vector3f squareToCosineHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToCosineHemisphere()
    static float squareToCosineHemispherePdf(const Vector3f &v);

    /// Warp a uniformly distributed square sample to a Beckmann distribution * cosine for the given 'alpha' parameter
    static Vector3f squareToBeckmann(const Point2f &sample, float alpha);

    /// Probability density of \ref squareToBeckmann()
    static float squareToBeckmannPdf(const Vector3f &m, float alpha);

    static int sampleDiscrete(const std::vector<float> &weights, 
                    float sample, float &pdf,
                    float &sampleRemapped);

    static int sampleDiscrete(float weights[],
                    int nWeights,
                    float sample, float &pdf,
                    float &sampleRemapped);
};

NORI_NAMESPACE_END

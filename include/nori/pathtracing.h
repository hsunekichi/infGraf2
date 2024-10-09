#pragma once

#include <nori/pathtracing.h>

#include <nori/common.h>
#include <nori/vector.h>
#include <nori/scene.h>
#include <nori/kdtree.h>



NORI_NAMESPACE_BEGIN


struct PathState
{
    bool intersectionComputed = false;
    bool intersected = false;
    Intersection intersection;

    Ray3f ray;
    int depth = 0;

    Color3f radiance = Color3f(0.0f);
    Color3f scatteringFactor = Color3f(1.0f);
};


struct Pth
{

enum IntegrationType {EMITTER, DIFFUSE, SPECULAR, NONE};


static Color3f nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                PathState &state)
{
    size_t nSamplesNes;
    return nextEventEstimation(scene, sampler, state, nSamplesNes, false);
}

static Color3f nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state,
                size_t &nSamplesNes,
                bool MIS);

static void sampleBSDF(const Scene *scene, Sampler *sampler,
            PathState &state, float &bsdfPdf);

static IntegrationType getIntegrationType(const PathState &mesh); 

};

NORI_NAMESPACE_END
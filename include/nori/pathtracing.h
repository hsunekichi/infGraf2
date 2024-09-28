#pragma once

#include <nori/pathtracing.h>

#include <nori/common.h>
#include <nori/vector.h>
#include <nori/scene.h>



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

static Color3f sampleRandomEmitter(const Scene *scene, Sampler *sampler, 
            const Point3f &surfaceP,
            Emitter *&emitterMesh,
            Point3f &lightP,
            float &lightPdf);

static Color3f nextEventEstimation(const Scene *scene, Sampler *sampler,
                const PathState &state, 
                bool MIS = false, bool applyF = true);

static void sampleBSDF(const Scene *scene, Sampler *sampler,
            PathState &state, float &bsdfPdf);

};

NORI_NAMESPACE_END
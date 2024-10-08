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

enum IntegrationType {EMITTER, DIFFUSE, SPECULAR, SUBSURFACE, NONE};


static Color3f sampleRandomEmitter(const Scene *scene, Sampler *sampler, 
            const Point3f &surfaceP,
            Emitter *&emitterMesh,
            Point3f &lightP,
            float &lightPdf);

static Color3f nextEventEstimation(const Scene *scene, Sampler *sampler,
                const PathState &state, 
                bool MIS = false, bool applyF = true)
{
    Vector3f wi;
    return nextEventEstimation(scene, sampler, state, wi, MIS, applyF);
}

static std::vector<Photon> generateSubsurfaceSamples(const Scene *scene, Sampler *sampler);

static Color3f nextEventEstimation(const Scene *scene, Sampler *sampler,
                const PathState &state, 
                Vector3f &wi,
                bool MIS = false, bool applyF = true);

static Color3f nextEventEstimationBSSRDF(const Scene *scene, Sampler *sampler,
                const PathState &state, const Point3f &pi,
                Vector3f &wi,
                bool MIS = false, bool applyF = true);

static void sampleBSDF(const Scene *scene, Sampler *sampler,
            PathState &state, float &bsdfPdf);

static IntegrationType getIntegrationType(const PathState &mesh); 


static void integrateSubsurfacePhotons(const Scene *scene,
                const PhotonMap &photons, 
                Sampler *sampler, PathState &state);


static void sampleBSSRDF(const Scene *scene, 
                Sampler *sampler, PathState &state,
                float &pdf);

static Point3f sampleBSSRDFpoint(const Scene *scene,
                Sampler *sampler,
                PathState &state,
                float &pdf);



};

NORI_NAMESPACE_END
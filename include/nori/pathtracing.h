#pragma once

#include <nori/pathtracing.h>

#include <nori/common.h>
#include <nori/vector.h>
#include <nori/scene.h>
#include <nori/kdtree.h>



NORI_NAMESPACE_BEGIN


struct PathState
{
    bool previous_diffuse = false; 
    size_t previous_n_samples = 1;
    float bsdfPdf = 0.0f; Point3f prevP = Point3f(0.0f);    
    Intersection intersection;

    Ray3f ray;
    int depth = 0;

    Color3f radiance = Color3f(0.0f);
    Color3f scatteringFactor = Color3f(1.0f);
};


struct Pth
{

enum IntegrationType {EMITTER, DIFFUSE, SUBSURFACE, SPECULAR, VOLUME, NONE};


static Color3f nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state,
                BSDFQueryRecord &bsdfQuery)
{
    float lightPdf, bsdfPdf;
    return nextEventEstimation(scene, sampler, state, bsdfQuery, lightPdf, bsdfPdf);
}

static BSDFQueryRecord initBSDFQuery(const Scene *scene, Sampler *sampler, const PathState &state);

static Color3f estimateDirectLight(const Scene *scene, 
                Sampler *sampler,
                const Point3f &p,
                float &lightPdf,
                Vector3f &g_wo,
                Emitter *&emitterMesh);

static Color3f nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state,
                BSDFQueryRecord &bsdfQuery,
                float &lightPdf, float &bsdfPdf);

static Color3f sampleBSDF(
        PathState &state,
        Sampler *sampler, 
        BSDFQueryRecord &bsdfQuery, 
        float &pdf);

static IntegrationType getIntegrationType(const Intersection &its); 


static Color3f integrateBSDF(const BSDF *bsdf, Sampler *sampler);
static Color3f integrateSkinSpecular(const BSDF *bsdf, std::unique_ptr<Sampler> sampler, float cosTh, float specWeight);

};

NORI_NAMESPACE_END
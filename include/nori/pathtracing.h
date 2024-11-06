#pragma once

#include <nori/pathtracing.h>

#include <nori/common.h>
#include <nori/vector.cu>
#include <nori/scene.h>
#include <nori/kdtree.h>
#include <nori/integrator.h>



NORI_NAMESPACE_BEGIN

struct Pth
{

enum IntegrationType {EMITTER, DIFFUSE, SPECULAR, NONE};


static Color3f nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state,
                BSDFQueryRecord &bsdfQuery)
{
    float lightPdf, bsdfPdf;
    return nextEventEstimation(scene, sampler, state, bsdfQuery, lightPdf, bsdfPdf);
}

static BSDFQueryRecord initBSDFQuery(const Scene *scene, const PathState &state);

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

};

NORI_NAMESPACE_END
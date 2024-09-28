#include <nori/pathtracing.h>


#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

#include <nori/scene.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN




Color3f Pth::sampleRandomEmitter(const Scene *scene, Sampler *sampler, 
            const Point3f &surfaceP,
            Emitter *&emitterMesh,
            Point3f &lightP,
            float &lightPdf)
{
    // Sample a random emitter
    float randomLightPdf = 0.0f;
    emitterMesh = scene->sampleEmitter(sampler, randomLightPdf);

    if (!emitterMesh)
        return Color3f(0.0f);


    // Sample a point on the emitter
    EmitterQueryRecord emitterQuery (ESolidAngle);
    emitterQuery.surfaceP = surfaceP;
    Color3f Le = emitterMesh->sampleLi(sampler, emitterQuery) / randomLightPdf;
    
    lightP = emitterQuery.lightP;
    lightPdf = emitterQuery.pdf * randomLightPdf;
    
    return Le;
}


bool checkVisibility (const Scene *scene, 
            const PathState &state,
            Emitter *emitterMesh,
            Vector3f &g_wi)
{
    Ray3f shadowRay(state.intersection.p, g_wi.normalized(), Epsilon, g_wi.norm() + Epsilon);

    // Check visibility
    bool objectSeesEmitter = true; // surface_wiNormalized.z() > 0.0f;

    Intersection lightIntersection;
    bool intersects = scene->rayIntersect(shadowRay, lightIntersection);

    //*********************** Compute Le ******************************
    return (objectSeesEmitter && !intersects)
        ||
        ( 
            objectSeesEmitter && intersects 
            && 
            (   (lightIntersection.p - state.intersection.p).norm() > g_wi.norm()
                ||
                lightIntersection.mesh->getEmitter() == emitterMesh
            )
        );
}



Color3f Pth::nextEventEstimation(const Scene *scene, Sampler *sampler,
                const PathState &state, 
                bool MIS, bool applyF)
{
    float lightPdf = 0.0f;
    Point3f lightP;
    Emitter *emitterMesh = nullptr;


    // Sample a point on a random emitter
    Color3f Le = sampleRandomEmitter(scene, sampler, state.intersection.p, 
            emitterMesh, lightP, lightPdf);
    
    if (Le == Color3f(0.0f))
        return Color3f(0.0f);


    Vector3f g_wi = (lightP - state.intersection.p);
    Vector3f surface_wiNormalized = state.intersection.toLocal(g_wi).normalized();


    //*********************** Sample emitter ******************************
    if (checkVisibility(scene, state, emitterMesh, g_wi))
    {
        // Evaluate bsdf
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
            BSDFQueryRecord bsdfQuery(state.intersection.toLocal(-state.ray.d), 
                    surface_wiNormalized, ESolidAngle);
        
        Color3f f = applyF ? bsdf->eval(bsdfQuery) : Color3f(1.0f);
        float bsdfPdf = bsdf->pdf(bsdfQuery);


        // Compute the geometric term
        float cosThetaP = Math::absCosTheta(surface_wiNormalized);
        float G = cosThetaP / g_wi.squaredNorm();

        // Combine all terms
        if (MIS)
        {
            float weight = Math::powerHeuristic(1, lightPdf, 1, bsdfPdf);
            return Le * f * G * weight;
        }
        else {
            return Le * f * G;
        }
    }

    return Color3f(0.0f);
}


void Pth::sampleBSDF(const Scene *scene, Sampler *sampler, PathState &state, float &pdf)
{
    Vector3f wi = state.intersection.toLocal(-state.ray.d).normalized();

    // Render non diffuse BSDF
    BSDFQueryRecord bsdfQuery(wi);

    Color3f f = state.intersection.mesh->getBSDF()->sample(bsdfQuery, sampler);
    pdf = state.intersection.mesh->getBSDF()->pdf(bsdfQuery);

    // Create the new ray
    Ray3f newRay(state.intersection.p, state.intersection.toWorld(bsdfQuery.wo), Epsilon, INFINITY);

    // Apply scattering factor
    state.scatteringFactor *= f;
    state.ray = newRay;
}



NORI_NAMESPACE_END
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



Color3f Pth::nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state, 
                Vector3f &wi,
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
    wi = g_wi;

    //*********************** Sample emitter ******************************
    if (checkVisibility(scene, state, emitterMesh, g_wi))
    {
        // Evaluate bsdf
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        BSDFQueryRecord bsdfQuery(state.intersection.toLocal(-state.ray.d), 
                surface_wiNormalized, ESolidAngle);
        bsdfQuery.uv = state.intersection.uv;
        
        
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
    bsdfQuery.isCameraRay = state.ray.isCameraRay;

    Color3f f = state.intersection.mesh->getBSDF()->sample(bsdfQuery, sampler);
    pdf = state.intersection.mesh->getBSDF()->pdf(bsdfQuery);

    // Create the new ray
    Ray3f newRay(state.intersection.p, state.intersection.toWorld(bsdfQuery.wo), Epsilon, INFINITY);
    newRay.isCameraRay = bsdfQuery.isCameraRay;

    // Apply scattering factor
    state.scatteringFactor *= f;
    state.ray = newRay;
}

std::vector<Photon> Pth::generateSubsurfaceSamples(const Scene *scene, Sampler *sampler)
{
    std::vector<Photon> photons;

    auto SS_meshes = scene->getSSMeshes();

    for (auto mesh : SS_meshes)
    {
        u_int32_t nTriangles = mesh->nTriangles();

        for (u_int32_t i = 0; i < nTriangles; i++)
        {
            float pdf;
            Point3f p; Normal3f n; Point2f uv;
            u_int32_t triangle_id;
            
            mesh->samplePosition(sampler, p, n, uv, pdf, triangle_id);
            photons.push_back(Photon(p, pdf, mesh));
        }
    }

    return photons;
}

void Pth::integrateSubsurface(const Scene *scene,
                const std::vector<Photon> &photons, 
                Sampler *sampler,
                PathState &state)
{
    // Static cast to SubsurfaceScattering bsdf
    //const SubsurfaceScattering &SSS = *static_cast<const SubsurfaceScattering *>(state.intersection.mesh->getBSDF());
    const BSDF &SSS = *state.intersection.mesh->getBSDF();
    
    BSDFQueryRecord bsdfQuery(state.intersection.toLocal(-state.ray.d));
    bsdfQuery.po = state.intersection.p;
    bsdfQuery.wo = state.intersection.toLocal(-state.ray.d);
    bsdfQuery.measure = ESolidAngle;
    bsdfQuery.isCameraRay = state.ray.isCameraRay;

    Color3f contributions = Color3f(0.0f);
    
    int N_SAMPLES = 100;
    float photonPdf = 1.0f / photons.size();
    
    for (int i = 0; i < N_SAMPLES; i++)
    {   
        // Choose random photon
        int randomPhoton = sampler->next1D() * photons.size();
        auto photon = photons[randomPhoton % photons.size()];

        if (photon.radiance == Color3f(0.0f)
            || photon.mesh != state.intersection.mesh)
            continue;

        bsdfQuery.pi = photon.p;
        bsdfQuery.ni = state.intersection.toLocal(photon.n);
        bsdfQuery.wi = state.intersection.toLocal(photon.d);
        Color3f f = SSS.eval(bsdfQuery);
        Color3f radiance = photon.radiance * f;

        contributions += radiance;
    }

    contributions = contributions / (N_SAMPLES * photonPdf);            
    
    state.radiance += contributions;
}

Pth::IntegrationType Pth::getIntegrationType(const PathState &state)
{
    const Mesh *mesh = state.intersection.mesh;
    const Emitter *emitter = mesh->getEmitter();

    if (emitter) { // Render emitter 
        return EMITTER;
    }
    else if (mesh->getBSDF()->isDiffuse()
        && !mesh->hasSubsurfaceScattering()) 
    { 
        // Render diffuse surface
        return DIFFUSE;
    }
    else if (mesh->hasSubsurfaceScattering()
        && !state.ray.isCameraRay)
    {
        return DIFFUSE;
    }
    else if (mesh->hasSubsurfaceScattering()
        && state.ray.isCameraRay) // Render subsurface scattering
    {
        return SUBSURFACE;
    }
    else { // Render specular surface
        return SPECULAR;
    }
}






NORI_NAMESPACE_END
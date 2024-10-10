#include <nori/pathtracing.h>

#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/bssrdf.h>

#include <nori/scene.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN




Color3f sampleRandomEmitter(const Scene *scene, Sampler *sampler, 
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

    Intersection lightIntersection;
    bool intersects = scene->rayIntersect(shadowRay, lightIntersection);

    // Check visibility
    bool objectSeesEmitter = true; //state.intersection.toLocal(g_wi).z() > 0.0f;

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


bool Pth::sampleBSSRDFpoint(const Scene *scene,
                Sampler *sampler,
                const PathState &state,
                BSDFQueryRecord &bsdfQuery,
                float &pdf)
{
    // Static cast to SubsurfaceScattering bsdf
    const BSSRDF &SSS = *static_cast<const BSSRDF *>(state.intersection.mesh->getBSDF());
    
    bsdfQuery.wi = state.intersection.toLocal(-state.ray.d);
    bsdfQuery.pi = state.intersection.p;
    bsdfQuery.po = state.intersection.p;
    bsdfQuery.measure = ESolidAngle;
    bsdfQuery.isCameraRay = state.ray.isCameraRay;
    bsdfQuery.uv = state.intersection.uv;
    bsdfQuery.frame = state.intersection.shFrame;
    bsdfQuery.scene = scene;
    bsdfQuery.mesh = state.intersection.mesh;

    return SSS.samplePoint(bsdfQuery, sampler, pdf);
}



Color3f Pth::nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state,
                BSDFQueryRecord &bsdfQuery,
                float &lightPdf, float &bsdfPdf)
{
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
        bsdfQuery.wo = surface_wiNormalized;
        bsdfQuery.wi = state.intersection.toLocal(-state.ray.d);
        bsdfQuery.measure = ESolidAngle;
        
        Color3f f = bsdf->eval(bsdfQuery);
        bsdfPdf = bsdf->pdf(bsdfQuery);


        // Compute the geometric term
        float cosThetaP = Math::absCosTheta(surface_wiNormalized);
        float G = cosThetaP / g_wi.squaredNorm();
        
        return Le * f * G;
    }

    return Color3f(0.0f);
}


void Pth::sampleBSDF(const Scene *scene, Sampler *sampler, PathState &state, float &pdf)
{
    Vector3f wi = state.intersection.toLocal(-state.ray.d).normalized();

    // Render non diffuse BSDF
    BSDFQueryRecord bsdfQuery(wi);
    bsdfQuery.pi = state.intersection.p;
    bsdfQuery.measure = ESolidAngle;
    bsdfQuery.isCameraRay = state.ray.isCameraRay;
    bsdfQuery.uv = state.intersection.uv;
    bsdfQuery.frame = state.intersection.shFrame;

    Color3f f = state.intersection.mesh->getBSDF()->sample(bsdfQuery, sampler);
    pdf = state.intersection.mesh->getBSDF()->pdf(bsdfQuery);

    // Create the new ray
    Ray3f newRay(state.intersection.p, state.intersection.toWorld(bsdfQuery.wo), Epsilon, INFINITY);
    newRay.isCameraRay = bsdfQuery.isCameraRay;

    // Apply scattering factor
    state.scatteringFactor *= f;
    state.ray = newRay;
}


Pth::IntegrationType Pth::getIntegrationType(const Intersection &its)
{
    const Mesh *mesh = its.mesh;
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
    else if (mesh->hasSubsurfaceScattering()) { // Render subsurface scattering
        return SUBSURFACE;
    }
    else { // Render specular surface
        return SPECULAR;
    }
}


/*


std::vector<Photon> generateSubsurfaceSamples(const Scene *scene, Sampler *sampler)
{
    std::vector<Photon> photons;

    auto SS_meshes = scene->getSSMeshes();

    for (auto mesh : SS_meshes)
    {
        u_int32_t nTriangles = mesh->nTriangles();
        u_int32_t nSamples = nTriangles;

        std::cout << "Generating " << nSamples << " photons for mesh " << mesh->getName() << std::endl;

        for (u_int32_t i = 0; i < nSamples; i++)
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


void Pth::integrateSubsurfacePhotons(const Scene *scene,
                const PhotonMap &photons,
                Sampler *sampler,
                PathState &state)
{
    // Static cast to SubsurfaceScattering bsdf
    const BSSRDF &SSS = *static_cast<const BSSRDF *>(state.intersection.mesh->getBSDF());
    
    BSDFQueryRecord bsdfQuery(state.intersection.toLocal(-state.ray.d));
    bsdfQuery.po = state.intersection.p;
    bsdfQuery.wo = state.intersection.toLocal(-state.ray.d);
    bsdfQuery.measure = ESolidAngle;
    bsdfQuery.isCameraRay = state.ray.isCameraRay;

    Color3f contributions = Color3f(0.0f);

    std::vector<Photon> nearest = photons.nearest_neighbors(state.intersection.p,
                                         (ulong)-1, 0.001f, state.intersection.mesh);


    if (nearest.size() == 0)
        return;

    std::cout << "Nearest size: " << nearest.size() << std::endl;
    
    for (auto &photon : nearest)
    {   
        // Choose random photon
        //int randomPhoton = sampler->next1D() * photons.size();
        //auto photon = photons[randomPhoton % photons.size()];

        bsdfQuery.pi = photon.p;
        bsdfQuery.ni = state.intersection.toLocal(photon.n);
        bsdfQuery.wi = state.intersection.toLocal(photon.d);
        Color3f f = SSS.eval(bsdfQuery);
        Color3f radiance = photon.radiance;

        contributions += radiance * f;
    }

    contributions = contributions / nearest.size();       
    //contributions *= state.intersection.mesh->meshArea();

    state.radiance += contributions;
}
*/







NORI_NAMESPACE_END
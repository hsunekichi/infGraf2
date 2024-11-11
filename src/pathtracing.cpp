#include <nori/pathtracing.h>

#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/warp.h>

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
            const Point3f &p,
            Emitter *emitterMesh,
            Vector3f &g_wi)
{
    Ray3f shadowRay(p, g_wi.normalized(), Epsilon, g_wi.norm() + Epsilon);

    Intersection lightIntersection;
    bool intersects = scene->rayIntersect(shadowRay, lightIntersection);

    //*********************** Compute Le ******************************
    return (!intersects)
        ||
        ( 
            intersects 
            && 
            (   (lightIntersection.p - p).norm() > g_wi.norm()
                ||
                lightIntersection.mesh->getEmitter() == emitterMesh
            )
        );
}


BSDFQueryRecord Pth::initBSDFQuery(const Scene *scene, const PathState &state)
{
    BSDFQueryRecord bsdfQuery;
    bsdfQuery.wi = state.intersection.vtoLocal(-state.ray.d);
    bsdfQuery.pi = state.intersection.p;
    bsdfQuery.fri = state.intersection.shFrame;

    bsdfQuery.measure = ESolidAngle;
    bsdfQuery.isCameraRay = state.ray.isCameraRay;
    bsdfQuery.uv = state.intersection.uv;
    
    bsdfQuery.mesh = state.intersection.mesh;
    bsdfQuery.scene = scene;

    return bsdfQuery;
}

Color3f Pth::estimateDirectLight(const Scene *scene, 
                Sampler *sampler,
                const Point3f &p,
                float &lightPdf,
                Vector3f &g_wo,
                Emitter *&emitterMesh)
{
    Point3f lightP;

    // Sample a point on a random emitter
    Color3f Le = sampleRandomEmitter(scene, sampler, p, 
            emitterMesh, lightP, lightPdf);
    
    if (Le == Color3f(0.0f))
        return Color3f(0.0f);

    g_wo = (lightP - p);

    //*********************** Sample emitter ******************************
    if (checkVisibility(scene, p, emitterMesh, g_wo))
        return Le;
    else
        return Color3f(0.0f);
}


Color3f Pth::nextEventEstimation(const Scene *scene, 
                Sampler *sampler,
                const PathState &state,
                BSDFQueryRecord &bsdfQuery,
                float &lightPdf, float &bsdfPdf)
{
    Vector3f g_wo; Emitter *emitterMesh;
    Color3f Le = estimateDirectLight(scene, sampler, state.intersection.p, lightPdf, g_wo, emitterMesh);

    if (Le != Color3f(0.0f) && state.intersection.vtoLocal(g_wo).z() > 0.0f)
    {
        Vector3f l_wo = bsdfQuery.fro.vtoLocal(g_wo).normalized();

        // Evaluate bsdf
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        bsdfQuery.wo = l_wo;
        bsdfQuery.wi = bsdfQuery.fro.vtoLocal(-state.ray.d);
        bsdfQuery.measure = ESolidAngle;
        
        Color3f f = bsdf->eval(bsdfQuery);
        bsdfPdf = bsdf->pdf(bsdfQuery);

        // Compute the geometric term
        float cosThetaP = Math::absCosTheta(l_wo);
        float G = cosThetaP / g_wo.squaredNorm();

        if (emitterMesh->getEmitterType() == EmitterType::EMITTER_ENVIRONMENT)
        {
            return Le*f *cosThetaP;
        }else
        {
            return Le * f * G;
        }
    }

    return Color3f(0.0f);
}


Color3f Pth::sampleBSDF(
        PathState &state,
        Sampler *sampler, 
        BSDFQueryRecord &bsdfQuery, 
        float &pdf)
{
    const BSDF *bsdf = state.intersection.mesh->getBSDF();

    Vector3f w_wi = bsdfQuery.fri.vtoWorld(bsdfQuery.wi);
    bsdfQuery.wi = bsdfQuery.fro.vtoLocal(w_wi).normalized();

    Color3f f = bsdf->sample(bsdfQuery, sampler);
    pdf = bsdf->pdf(bsdfQuery);

    // Create the new ray
    Ray3f newRay(bsdfQuery.po, bsdfQuery.fro.vtoWorld(bsdfQuery.wo), 0.001, INFINITY);
    newRay.isCameraRay = bsdfQuery.isCameraRay;

    // Apply scattering factor
    state.ray = newRay;
    return f;
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



Color3f Pth::integrateBSDF(const BSDF *bsdf, Sampler *sampler)
{
    Color3f result(0.0f, 0.0f, 0.0f);
    int nSteps = 10000;

    for (int i = 0; i < nSteps; i++)
    {
        BSDFQueryRecord bsdfQuery;
        bsdfQuery.wi = Warp::squareToCosineHemisphere(sampler->next2D());

        Color3f f = bsdf->sample(bsdfQuery, sampler);
        //float pdf = bsdf->pdf(bsdfQuery);

        result += f;
    }

    return result / nSteps;
}


Color3f Pth::computeInScattering(const Scene *scene, 
        Sampler *sampler, 
        const Ray3f &ray, 
        const Point3f &p,
        float sigma_s,
        float sigma_t,
        float g) 
{    
    Vector3f g_wo; Emitter *emitterMesh; float lightPdf;
    Color3f Le = estimateDirectLight(scene, sampler, p, lightPdf, g_wo, emitterMesh);

    if (Le != Color3f(0.0f))
    {
        Color3f direct = Le / g_wo.squaredNorm();

        // Henyey-Greenstein phase function
        float cosTheta = ray.d.dot(g_wo.normalized());
        float phaseFunctionNormalization = 1.0f / (4 * M_PI);
        float phaseFunction = phaseFunctionNormalization * (1.0f - g * g) / std::pow(1.0f + g * g - 2.0f * g * cosTheta, 1.5f);
    
        return sigma_s * phaseFunction * direct;
    }

    return Color3f(0.0f); 
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
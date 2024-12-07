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
    Ray3f shadowRay(p, g_wi.normalized(), Epsilon, g_wi.norm() - Epsilon);

    Intersection lightIntersection;
    bool intersects = scene->rayIntersect(shadowRay, lightIntersection);
    
    return !intersects
            || (intersects && lightIntersection.mesh->getEmitter() == emitterMesh);
}


BSDFQueryRecord Pth::initBSDFQuery(const Scene *scene, Sampler *sampler, const PathState &state)
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
    bsdfQuery.sampler = sampler;

    return bsdfQuery;
}

Color3f Pth::estimateDirectLight(const Scene *scene, 
                Sampler *sampler,
                const Frame &frame,
                float &lightPdf,
                Vector3f &g_wo,
                Emitter *&emitterMesh)
{
    Point3f lightP;

    // Sample a point on a random emitter
    Color3f Le = sampleRandomEmitter(scene, sampler, frame.o, 
            emitterMesh, lightP, lightPdf);
    
    if (Le == Color3f(0.0f))
        return Color3f(0.0f);

    g_wo = (lightP - frame.o);

    //*********************** Sample emitter ******************************
    bool objectSeesEmitter = frame.vtoLocal(g_wo).z() > 0.0f;
    if (objectSeesEmitter && checkVisibility(scene, frame.o, emitterMesh, g_wo))
        return Le;
    else
        return Color3f(0.0f);
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
    Color3f Le = estimateDirectLight(scene, sampler, bsdfQuery.fro, 
                                    lightPdf, g_wo, emitterMesh);

    if (Le != Color3f(0.0f))
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
            if (scene->getIntegrator()->toString() == "Whitted[]" 
                || scene->getIntegrator()->toString() == "directWhitted[]") 
            {
                return Le * f * cosThetaP;
            } else {
                return Le * f * cosThetaP / M_PI;
            }
            
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
    Ray3f newRay(bsdfQuery.po, bsdfQuery.fro.vtoWorld(bsdfQuery.wo), Epsilon, INFINITY);
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
    int nSteps = 2048;

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

template <std::size_t N>
std::array<float, N> precompute_sin(float constant)
{
    std::array<float, N> array {};

    for(unsigned i = 0; i < N; i++){
        float th = float(i) / float(N - 1) * M_PI * constant;
        array[i] = std::sin(th);
    }

    return array;
}

template <std::size_t N>
std::array<float, N> precompute_cos(float constant)
{
    std::array<float, N> array {};

    for(unsigned i = 0; i < N; i++){
        float th = float(i) / float(N - 1) * M_PI * constant;
        array[i] = std::cos(th);
    }

    return array;
}

Color3f Pth::integrateSkinSpecular(const BSDF *bsdf, std::unique_ptr<Sampler> sampler, float costheta, float specWeight)
{
    Color3f result(0.0f, 0.0f, 0.0f);

    constexpr unsigned nSteps = 512;
    auto sinphi = precompute_sin<nSteps>(2.0f);
    auto cosphi = precompute_cos<nSteps>(2.0f);
    auto sinth =  precompute_sin<nSteps>(1.0f/2.0f);
    auto costh =  precompute_cos<nSteps>(1.0f/2.0f);

    Vector3f V = Vector3f(0.0, sqrt(1.0 - costheta * costheta), costheta);

    for (unsigned i = 0; i < nSteps; i++)
    {
        float cosp = cosphi[i];
        float sinp = sinphi[i];

        Color3f localsum = 0.0f;

        for (unsigned j = 0; j < nSteps; j++)
        {
            float sint = sinth[j];
            float cost = costh[j];
            Vector3f L = Vector3f(sinp * sint, cosp * sint, cost);

            BSDFQueryRecord bsdfQuery;
            bsdfQuery.wi = L;
            bsdfQuery.wo = V;

            localsum += bsdf->eval(bsdfQuery) * specWeight * sint;
        }

        result += localsum * (M_PI / 2.0) / float(nSteps);
    }

    return result * (2.0 * M_PI) / (float(nSteps));
}


NORI_NAMESPACE_END
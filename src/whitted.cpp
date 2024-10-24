

#include <nori/scene.h>
#include <nori/integrator.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>
#include <nori/kdtree.h>

NORI_NAMESPACE_BEGIN

class Whitted : public Integrator {
public:

    Whitted(const PropertyList &props) {}

    Color3f specularIntegration(const Scene *scene, 
            Sampler *sampler,
            const Ray3f &ray,
            Intersection &its,
            int depth) const
    {
        // Apply roussian roulette
        float roulettePdf = 1.0f;
        if (depth > 3)
        {
            roulettePdf = 0.95f; 
            if (sampler->next1D() > roulettePdf)
                return Color3f(0.0f);
        }

        PathState state;
        state.intersection = its;
        state.ray = ray;
        
        // Sample the BSDF
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        if (fp == Color3f(0.0f))
            return Color3f(0.0f);

        float pdf;
        Color3f f = Pth::sampleBSDF(state, sampler, query, pdf);

        pdf = roulettePdf;
        depth++;

        // Compute the contribution
        return f * fp * Li(scene, sampler, state.ray, depth);
    }

    Color3f integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {   
        PathState state;
        state.intersection = its;
        state.ray = ray;

        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        float lightPdf, bsdfPdf;
        Color3f result = fp * Pth::nextEventEstimation(scene, sampler, state, query, lightPdf, bsdfPdf);
        return result;
    }

    Color3f integrateSubsurface(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {
        int nSamples = 4;
        Color3f radiance = Color3f(0.0f);

        for (int i = 0; i < nSamples ; i++)
        {
            PathState state;
            state.intersection = its;
            state.ray = ray;

            const BSDF *bsdf = state.intersection.mesh->getBSDF();
            
            auto query = Pth::initBSDFQuery(scene, state);
            Color3f fp = bsdf->samplePoint(query, sampler);

            radiance += fp * Pth::nextEventEstimation(scene, sampler, state, query);
        }

        return radiance / nSamples;
    }

    Color3f Li (const Scene *scene, Sampler *sampler,
            const Ray3f &ray,
            int depth) const
    {
        Intersection its;
        Color3f radiance(0.0f);

        /* Find the surface that is visible in the requested direction */
        if (!scene->rayIntersect(ray, its)){
            if (scene->getEnvironmentalEmitter() != nullptr){
                EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
                emitterQuery.lightP = ray.d*1e15;
                radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery);
                return radiance;
            }else{
                return Color3f(0.0f);
            }

        }

        Pth::IntegrationType type = Pth::getIntegrationType(its);        

        switch(type)
        {
            case Pth::EMITTER:
            {
                // Render emitter
                EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
                emitterQuery.lightP = its.p;
                radiance = its.mesh->getEmitter()->eval(emitterQuery);
                break;
            }
            case Pth::DIFFUSE:
                // Render diffuse surface
                radiance = integrateDiffuse(scene, sampler, ray, its);
                break;

            case Pth::SUBSURFACE:
                radiance = integrateDiffuse(scene, sampler, ray, its);
                break;

            case Pth::SPECULAR:
                radiance = specularIntegration(scene, sampler, ray, its, depth);
                break;
            
            default:
                break;
        }

        return radiance;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {        
        return Li(scene, sampler, ray, 0);
    }

    std::string toString() const {
        return "Whitted[]";
    }

    protected:
        PhotonMap photonMap;
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END


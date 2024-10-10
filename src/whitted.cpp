

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
        float pointPdf; 
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        bool valid = bsdf->samplePoint(query, sampler, pointPdf);

        if (!valid)
            return Color3f(0.0f);

        float pdf;
        Color3f f = Pth::sampleBSDF(state, sampler, query, pdf);

        pdf = roulettePdf * pointPdf;
        depth++;

        // Compute the contribution
        return f * Li(scene, sampler, state.ray, depth) / pdf;
    }

    Color3f integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {   
        PathState state;
        state.intersection = its;
        state.ray = ray;

        float pointPdf; 
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        bool valid = bsdf->samplePoint(query, sampler, pointPdf);

        if (!valid)
            return Color3f(0.0f);

        return Pth::nextEventEstimation(scene, sampler, state, query) / pointPdf;
    }

    Color3f integrateSubsurface(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {
        int nSamples = 32;
        Color3f radiance = Color3f(0.0f);

        for (int i = 0; i < nSamples ; i++)
        {
            PathState state;
            state.intersection = its;
            state.ray = ray;

            float pointPdf; 
            const BSDF *bsdf = state.intersection.mesh->getBSDF();
            
            auto query = Pth::initBSDFQuery(scene, state);
            bool valid = bsdf->samplePoint(query, sampler, pointPdf);

            if (!valid)
                continue;

            radiance += Pth::nextEventEstimation(scene, sampler, state, query) / pointPdf;
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
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

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
                radiance = integrateSubsurface(scene, sampler, ray, its);
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


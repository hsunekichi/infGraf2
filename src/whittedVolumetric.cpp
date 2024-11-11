

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

class WhittedVolumetric : public Integrator {
public:

    WhittedVolumetric(const PropertyList &props) 
    {
        sigma_s = props.getFloat("sigma_s", 0.f);
        // sigma_s = props.getFloat("sigma_s", 0.001f);
        sigma_t = props.getFloat("sigma_t", 0.f);

        g = props.getFloat("g", 0.f);
        helios_coeff = props.getFloat("helios_coeff", 0.f);
    }

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

        return fp * Pth::nextEventEstimation(scene, sampler, state, query);
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
            // Render emitter
            EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
            emitterQuery.lightP = ray.d*1e15;

            if (scene->getEnvironmentalEmitter() != nullptr)
            {
                radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery);
                return radiance;
            } else{
                return Color3f(0.0f);
            }
        }

        /*********** Sample volume point intersection **********************/
        float volumeD = Math::abs(std::log(1 - sampler->next1D()) / sigma_t);
        if (volumeD < its.t)
        {
            if (sigma_s==0 || helios_coeff==0)
                return 0.0f;
            
            float volumePdf = sigma_t * std::exp(-sigma_t * volumeD);
            Point3f lp = ray(volumeD);

            Color3f inScatter = Pth::computeInScattering(scene, sampler, ray,
                    lp, sigma_s, sigma_t, g) * helios_coeff / volumePdf;
    
            
            // Dispersión volumétrica
            return inScatter;           
        }

        Pth::IntegrationType type = Pth::getIntegrationType(its);        

        switch(type)
        {
            case Pth::EMITTER:
            {
                // Render emitter
                EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
                emitterQuery.lightP = its.p;
                radiance += its.mesh->getEmitter()->eval(emitterQuery);
                break;
            }
            case Pth::DIFFUSE:
                // Render diffuse surface
                radiance += integrateDiffuse(scene, sampler, ray, its);
                break;

            case Pth::SUBSURFACE:
                radiance += integrateDiffuse(scene, sampler, ray, its);
                break;

            case Pth::SPECULAR:
                radiance += specularIntegration(scene, sampler, ray, its, depth);
                break;
            
            default:
                break;
        }

        return radiance * std::exp(-sigma_t * its.t);
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
    private:
        float sigma_s; // Scatter coefficient
        float sigma_t; // Extinction Coefficient
        float g; // Scatter Phase Coefficient (G) 
        float helios_coeff; // Helios Coefficient

};

NORI_REGISTER_CLASS(WhittedVolumetric, "whitted_volume");
NORI_NAMESPACE_END


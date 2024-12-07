

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

    Whitted(const PropertyList &props) 
    {
        sigma_s = props.getFloat("sigma_s", 0.f);
        // sigma_s = props.getFloat("sigma_s", 0.001f);
        sigma_t = props.getFloat("sigma_t", 0.001f);

        g = props.getFloat("g", -0.1f);
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
        
        auto query = Pth::initBSDFQuery(scene, sampler, state);
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
        
        auto query = Pth::initBSDFQuery(scene, sampler, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        if (query.po.y() > 1e-5)
        {
            int a = 0;
        }

        float lightPdf, bsdfPdf;
        Color3f result = fp * Pth::nextEventEstimation(scene, sampler, state, query, lightPdf, bsdfPdf);
        return result;
    }

    Color3f shadeEnvironment(const Scene *scene, const Ray3f &ray, Intersection &its) const
    {
        Color3f radiance = BLACK;

        if (scene->getEnvironmentalEmitter() != nullptr)
        {
            EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
            emitterQuery.lightP = ray.d*1e15;
            
            radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery);
        }

        return radiance;
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
            
            auto query = Pth::initBSDFQuery(scene, sampler, state);
            Color3f fp = bsdf->samplePoint(query, sampler);

            radiance += fp * Pth::nextEventEstimation(scene, sampler, state, query);
        }

        return radiance / nSamples;
    }

    Color3f integrateVolume(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {
        const auto phaseFunction = [g = g]
                            (const Vector3f &wo,
                            const Vector3f &wi) -> float
        {
            return (1.0f / (4 * M_PI)) * (1.0f - g * g) / std::pow(1.0f + g * g - 2.0f * g * std::abs(wo.dot(wi)), 1.5f);
        };

        Point3f lp = ray(its.t);
        float transmittance = std::exp(-sigma_t * its.t);
        float density = sigma_t * transmittance;
        float atmos_pdf = density; //(density.x + density.y + density.z) / 3.0f;
        Color3f mediumScattering = transmittance * sigma_s * helios_coeff / atmos_pdf;

        /******************* Direct *********************/
        Vector3f g_wo; Emitter *emitterMesh; float lightPdf;
        Color3f Le = Pth::estimateDirectLight(scene, sampler, lp, lightPdf, g_wo, emitterMesh);

        if (Le != Color3f(0.0f))
        {
            Color3f direct = Le / g_wo.squaredNorm();
            float phase = phaseFunction(-ray.d, g_wo.normalized());

            Le = mediumScattering * phase * direct;
        }

        return Le;
    }

    Color3f Li (const Scene *scene, Sampler *sampler,
            const Ray3f &ray,
            int depth, bool godrays=true) const
    {
        Intersection its;
        Color3f radiance(0.0f);

        /* Find the surface that is visible in the requested direction */
        if (!scene->rayIntersect(ray, its))
        {
            return shadeEnvironment(scene, ray, its);
        }

        Pth::IntegrationType type = Pth::getIntegrationType(its);   

        #ifndef DISABLE_VOLUME_INTEGRATION
        if (sigma_s != 0 || helios_coeff != 0)
        {
            float volume_t = Math::abs(std::log(1 - sampler->next1D()) / sigma_t);
            if (volume_t < its.t)
            {
                its.t = volume_t;
                type = Pth::VOLUME;
            }
        }
        #endif

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

            case Pth::VOLUME:
                radiance += integrateVolume(scene, sampler, ray, its);
                break;
            
            default:
                break;
        }

        #ifndef DISABLE_VOLUME_INTEGRATION
        if (type != Pth::VOLUME)
        {
            radiance *= std::exp(-sigma_t * its.t);
        }
        #endif

        return radiance;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {        
        return Li(scene, sampler, ray, 0);
    }

    std::string toString() const {
        return "Whitted[]";
    }

    private:
        float sigma_s; // Scatter coefficient
        float sigma_t; // Extinction Coefficient
        float g; // Scatter Phase Coefficient (G) 
        float helios_coeff; // Helios Coefficient

};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END


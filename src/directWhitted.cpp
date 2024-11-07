

#include <nori/scene.h>
#include <nori/integrator.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN

class directWhitted : public Integrator {
public:

    directWhitted(const PropertyList &props) {}

    // Integrate a diffuse surface
    Color3f integrateDiffuse(const Scene *scene, Sampler *sampler, const Ray3f &ray,
                Intersection &its) const
    {
        PathState state;
        state.intersection = its;
        state.ray = ray;

        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, sampler, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        float lightPdf, bsdfPdf;
        Color3f result = fp * Pth::nextEventEstimation(scene, sampler, state, query, lightPdf, bsdfPdf);
        return result;
    }

    Color3f Li (const Scene *scene, Sampler *sampler, const Ray3f &ray, int depth) const
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
                return BLACK;
                break;
            
            default:
                break;
        }

        return radiance;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        PathState state;
        state.ray = ray;
        
        return Li(scene, sampler, ray, 0);
    }

    std::string toString() const {
        return "directWhitted[]";
    }
};

NORI_REGISTER_CLASS(directWhitted, "direct_whitted");
NORI_NAMESPACE_END


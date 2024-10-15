

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
    Color3f diffuseIntegration(const Scene *scene, Sampler *sampler, const Ray3f &ray,
                Intersection &intersection) const
    {
        PathState state;
        state.ray = ray;
        state.intersection = intersection;

        BSDFQueryRecord bsdfQuery;
        Color3f incoming = Pth::nextEventEstimation(scene, sampler, state, bsdfQuery);
        return incoming;
    }

    Color3f Li (const Scene *scene, Sampler *sampler, const Ray3f &ray, int depth) const
    {
        /* Find the surface that is visible in the requested direction */
        Intersection intersection;
        if (!scene->rayIntersect(ray, intersection))
            return Color3f(0.0f);

        /* Retrieve the emitter associated with the surface */
        const Emitter *emitter = intersection.mesh->getEmitter();
        
        Color3f radiance = Color3f(0.0f);

        if (emitter != nullptr)
        {
            // Render emitter
            EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
            emitterQuery.lightP = intersection.p;
            radiance += emitter->eval(emitterQuery);
        } 
        else if (intersection.mesh->getBSDF()->isDiffuse()) 
        {
            // Render diffuse surface
            radiance += diffuseIntegration(scene, sampler, ray, intersection);
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




#include <nori/scene.h>
#include <nori/integrator.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN

class Whitted : public Integrator {
public:

    Whitted(const PropertyList &props) {}

    // Integrate a diffuse surface
    Color3f diffuseIntegration(const Scene *scene, Sampler *sampler, const Ray3f &ray,
                Intersection &intersection) const
    {
        PathState state;
        state.ray = ray;
        state.intersection = intersection;

        return Pth::nextEventEstimation(scene, sampler, state);
    }

    Color3f specularIntegration(const Scene *scene, Sampler *sampler, const Ray3f &ray,
                Intersection &intersection,
                int &depth) const
    {
        // Apply roussian roulette
        float roulettePdf = 1.0f;
        if (depth > 3)
        {
            roulettePdf = 0.95f; 
            if (sampler->next1D() > roulettePdf)
                return Color3f(0.0f);
        }

        Vector3f wi = intersection.toLocal(-ray.d).normalized();

        // Render non diffuse BSDF
        BSDFQueryRecord bsdfQuery(wi);
        Color3f f = intersection.mesh->getBSDF()->sample(bsdfQuery, sampler->next2D());

        // Compute the new ray
        Ray3f newRay(intersection.p, intersection.toWorld(bsdfQuery.wo), Epsilon, INFINITY);

        depth++;

        // Compute the contribution
        return f * Li(scene, sampler, newRay, depth) / roulettePdf;
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
        else // Render specular surface
        {
            radiance += specularIntegration(scene, sampler, ray, intersection, depth);
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
        return "Whitted[]";
    }
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END


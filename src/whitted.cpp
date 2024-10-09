

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

        Vector3f wi = its.toLocal(-ray.d).normalized();

        // Render non diffuse BSDF
        BSDFQueryRecord bsdfQuery(wi);
        Color3f f = its.mesh->getBSDF()->sample(bsdfQuery, sampler);

        // Compute the new ray
        Ray3f newRay(its.p, its.toWorld(bsdfQuery.wo), Epsilon, INFINITY);
        newRay.isCameraRay = ray.isCameraRay;

        depth++;

        // Compute the contribution
        return f * Li(scene, sampler, newRay, depth) / roulettePdf;
    }

    Color3f integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {   
        PathState state;
        state.intersection = its;
        state.ray = ray;

        return Pth::nextEventEstimation(scene, sampler, state);
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


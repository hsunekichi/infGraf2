

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

    Color3f subsurfaceIntegration(const Scene *scene, Sampler *sampler, const Ray3f &ray,
                Intersection &intersection, bool precomputing) const
    {
        PathState state;
        state.ray = ray;
        state.intersection = intersection;

        if (precomputing)
        {
            return Pth::nextEventEstimation(scene, sampler, state);
        }
        else
        {
            const auto &photons = scene->getSubsurfacePhotons();

            BSDFQueryRecord bsdfQuery(intersection.toLocal(-ray.d));
            bsdfQuery.surfaceP = intersection.p;
            bsdfQuery.wo = intersection.toLocal(-ray.d);
            bsdfQuery.measure = ESolidAngle;


            float minD = std::numeric_limits<float>::infinity();
            Color3f radiance = Color3f(0.0f);

            for (const auto &photon : photons)
            {
                float d = (photon.p - intersection.p).norm();
                
                if (d < minD) //  && d < 0.001f
                {
                    bsdfQuery.wi = intersection.toLocal(photon.d);

                    minD = d;
                    Color3f f = intersection.mesh->getBSDF()->eval(bsdfQuery);
                    radiance = f * photon.radiance;
                }
            }

            return radiance;
        }
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

    Color3f Li (const Scene *scene, Sampler *sampler, const Ray3f &ray, int depth, bool precomputing=false) const
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
        else if (intersection.mesh->getBSDF()->isDiffuse() && !intersection.mesh->hasSubsurfaceScattering()) 
        {
            // Render diffuse surface
            radiance += diffuseIntegration(scene, sampler, ray, intersection);
        }
        else if (intersection.mesh->hasSubsurfaceScattering())
        {
            // Render subsurface scattering
            radiance += subsurfaceIntegration(scene, sampler, ray, intersection, precomputing);
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

    void precomputeLi(const Scene *scene, Sampler *sampler,
            Photon &ph) const
    {
        Vector3f off (0.0f, 0.0f, 1.0f);
        Point3f o = ph.p + Epsilon*off;
        Ray3f ray(o, -off); 

        Intersection intersection;
        if (!scene->rayIntersect(ray, intersection))
            return;

        PathState state;
        state.ray = ray;
        state.intersection = intersection;
        
        Color3f direct = Pth::nextEventEstimation(scene, sampler, state, false, false);

        ph.d = (intersection.p - ph.p).normalized();
        ph.radiance = direct; 
    }

    std::string toString() const {
        return "Whitted[]";
    }
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END


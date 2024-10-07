

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

    void preprocess(const Scene *scene, Sampler *sampler) 
    {
        photons = Pth::generateSubsurfaceSamples(scene, sampler);

        for (auto &photon : photons)
        {
            precomputeLi(scene, sampler, photon);
        }
    }


    void specularIntegration(const Scene *scene, Sampler *sampler,
            PathState &state) const
    {
        // Apply roussian roulette
        float roulettePdf = 1.0f;
        if (state.depth > 3)
        {
            roulettePdf = 0.95f; 
            if (sampler->next1D() > roulettePdf)
                return;
        }

        Vector3f wi = state.intersection.toLocal(-state.ray.d).normalized();

        // Render non diffuse BSDF
        BSDFQueryRecord bsdfQuery(wi);
        Color3f f = state.intersection.mesh->getBSDF()->sample(bsdfQuery, sampler);

        // Compute the new ray
        Ray3f newRay(state.intersection.p, state.intersection.toWorld(bsdfQuery.wo), Epsilon, INFINITY);
        newRay.isCameraRay = state.ray.isCameraRay;

        state.depth++;
        Color3f scatteringFactor = f / roulettePdf;

        // Compute the contribution
        Li(scene, sampler, state);
        state.radiance *= scatteringFactor;
        state.depth--;
    }

    void Li (const Scene *scene, Sampler *sampler, PathState &state) const
    {
        /* Find the surface that is visible in the requested direction */
        if (!scene->rayIntersect(state.ray, state.intersection))
            return;

        Pth::IntegrationType type = Pth::getIntegrationType(state);        

        switch(type)
        {
            case Pth::EMITTER:
            {
                // Render emitter
                EmitterQueryRecord emitterQuery(-state.ray.d, EDiscrete);
                emitterQuery.lightP = state.intersection.p;
                state.radiance += state.intersection.mesh->getEmitter()->eval(emitterQuery);
                break;
            }
            case Pth::DIFFUSE:
                // Render diffuse surface
                state.radiance += Pth::nextEventEstimation(scene, sampler, state);
                break;

            case Pth::SUBSURFACE:
                Pth::integrateSubsurface(scene, photons, sampler, state);
                break;

            case Pth::SPECULAR:
                specularIntegration(scene, sampler, state);
                break;
            
            default:
                break;
        }
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        PathState state;
        state.ray = ray;
        
        Li(scene, sampler, state);
        return state.radiance;
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
        
        Vector3f wi;
        Color3f direct = Pth::nextEventEstimation(scene, sampler, state, wi, false, false);

        ph.d = wi.normalized();
        ph.radiance = direct; 
        ph.n = intersection.shFrame.n;
    }

    std::string toString() const {
        return "Whitted[]";
    }

    protected:
        std::vector<Photon> photons;
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END


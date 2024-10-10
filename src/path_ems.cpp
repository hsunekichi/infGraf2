

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>


NORI_NAMESPACE_BEGIN


class Path_ems : public Integrator {
public:

    Path_ems(const PropertyList &props) {}

    void integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {   
        // Sample the contribution of a random emitter
        BSDFQueryRecord query = Pth::initBSDFQuery(scene, state);
        Color3f directLight = Pth::nextEventEstimation(scene, sampler, state, query);
        state.radiance += state.scatteringFactor * directLight * 0.5f;

        // Sample the BSDF
        float bsdfPdf;
        Color3f f = Pth::sampleBSDF(state, sampler, query, bsdfPdf);   
        state.scatteringFactor *= f;


        Point3f surfaceP = state.intersection.p;
        
        state.intersected = scene->rayIntersect(state.ray, state.intersection);
        state.intersectionComputed = true;

        if (state.intersected && state.intersection.mesh->isEmitter())
        {
            const Emitter *emitter = state.intersection.mesh->getEmitter();

            EmitterQueryRecord emitterQuery(surfaceP, state.intersection.p, state.intersection.vtoLocal(-state.ray.d), EMeasure::EDiscrete);
            state.radiance += state.scatteringFactor * emitter->eval(emitterQuery) * 0.5f;
            
            // Terminate the path
            state.scatteringFactor = Color3f(0.0f);
        }
    }

    void integrateSpecular(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Sample the specular BSDF
        float bsdfPdf;
        BSDFQueryRecord query = Pth::initBSDFQuery(scene, state);
        Color3f f = Pth::sampleBSDF(state, sampler, query, bsdfPdf); 
        state.scatteringFactor *= f;    
    }

    void integrateEmitter(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Retrieve the emitter associated with the surface
        const Emitter *emitter = state.intersection.mesh->getEmitter();

        EmitterQueryRecord emitterQuery (state.intersection.vtoLocal(-state.ray.d), EMeasure::EDiscrete);
        emitterQuery.lightP = state.intersection.p;
        state.radiance += state.scatteringFactor * emitter->eval(emitterQuery);

        // Terminate the path
        state.scatteringFactor = Color3f(0.0f);
    }

    void sampleIntersection(const Scene *scene, Sampler *sampler, PathState &state) const
    {
        Pth::IntegrationType integrationType = Pth::getIntegrationType(state.intersection);
        
        switch (integrationType)
        {
            case Pth::EMITTER:
                integrateEmitter(scene, sampler, state);
                break;
            case Pth::DIFFUSE:
                integrateDiffuse(scene, sampler, state);
                break;
            case Pth::SPECULAR:
                integrateSpecular(scene, sampler, state);
                break;
            default:
                break;
        }
    }

    // Compute radiance over a full path
    void Li(const Scene *scene, Sampler *sampler, PathState &state) const 
    {
        while (state.scatteringFactor != Color3f(0.0f))
        {
            /* Find the surface that is visible in the requested direction */
            if ((state.intersectionComputed && !state.intersected)
                ||
                (!state.intersectionComputed && !scene->rayIntersect(state.ray, state.intersection)))
            {
                state.scatteringFactor = Color3f(0.0f);
                return;
            }
           
            state.intersectionComputed = false;

            if (state.depth > 3)
            {
                // Apply roussian roulette
                float roulettePdf = 0.95f;
                if (sampler->next1D() > roulettePdf)
                {
                    state.scatteringFactor = Color3f(0.0f);
                    return;
                }
                else { // Apply roulette pdf
                    state.scatteringFactor /= roulettePdf;
                }
            }

            // Integrate the current intersection
            sampleIntersection(scene, sampler, state);

            state.depth++;
        }
    }

    // Compute radiance over a given ray
    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        // Initialize path state
        PathState state;
        state.ray = ray;

        Li(scene, sampler, state);

        return state.radiance;
    }

    std::string toString() const {
        return "Path_ems[]";
    }
};

NORI_REGISTER_CLASS(Path_ems, "path_ems");
NORI_NAMESPACE_END


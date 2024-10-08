

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
        Color3f directLight = Pth::nextEventEstimation(scene, sampler, state);
        state.radiance += state.scatteringFactor * directLight * 0.5f;

        // Sample the BSDF
        float bsdfPdf;
        Pth::sampleBSDF(scene, sampler, state, bsdfPdf);

        Point3f surfaceP = state.intersection.p;
        
        state.intersected = scene->rayIntersect(state.ray, state.intersection);
        state.intersectionComputed = true;

        if (state.intersected && state.intersection.mesh->isEmitter())
        {
            const Emitter *emitter = state.intersection.mesh->getEmitter();

            EmitterQueryRecord emitterQuery(surfaceP, state.intersection.p, state.intersection.toLocal(-state.ray.d), EMeasure::EDiscrete);
            state.radiance += state.scatteringFactor * emitter->eval(emitterQuery) * 0.5f;
            
            // Terminate the path
            state.scatteringFactor = Color3f(0.0f);
        }
    }

    void integrateSubsurface(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Sample the subsurface scattering BSDF
        float pointPdf;
        Point3f pi = Pth::sampleBSSRDFpoint(scene, sampler, state, pointPdf);
        Vector3f wi;

        Color3f Li = Pth::nextEventEstimationBSSRDF(scene, sampler, state, pi, wi);
        state.radiance += state.scatteringFactor * Li * 0.5f / pointPdf;

        // Sample the subsurface scattering BSDF
        float bsdfPdf;
        Pth::sampleBSSRDF(scene, sampler, state, bsdfPdf);

        Point3f surfaceP = pi;
        
        state.intersected = scene->rayIntersect(state.ray, state.intersection);
        state.intersectionComputed = true;

        if (state.intersected && state.intersection.mesh->isEmitter())
        {
            const Emitter *emitter = state.intersection.mesh->getEmitter();

            EmitterQueryRecord emitterQuery(surfaceP, state.intersection.p, state.intersection.toLocal(-state.ray.d), EMeasure::EDiscrete);
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
        Pth::sampleBSDF(scene, sampler, state, bsdfPdf);        
    }

    void integrateEmitter(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Retrieve the emitter associated with the surface
        const Emitter *emitter = state.intersection.mesh->getEmitter();

        EmitterQueryRecord emitterQuery (state.intersection.toLocal(-state.ray.d), EMeasure::EDiscrete);
        emitterQuery.lightP = state.intersection.p;
        state.radiance += state.scatteringFactor * emitter->eval(emitterQuery);

        // Terminate the path
        state.scatteringFactor = Color3f(0.0f);
    }

    void sampleIntersection(const Scene *scene, Sampler *sampler, PathState &state) const
    {
        Pth::IntegrationType integrationType = Pth::getIntegrationType(state);
        
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
            case Pth::SUBSURFACE:
                integrateSubsurface(scene, sampler, state);
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


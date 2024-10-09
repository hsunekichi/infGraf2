

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/kdtree.h>

#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN





class Path_mis : public Integrator {
public:

    Path_mis(const PropertyList &props) {}

    void integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {   
        size_t nSamplesNES;

        // Sample the contribution of a random emitter
        Color3f directLight = Pth::nextEventEstimation(scene, sampler, state, nSamplesNES, true);
        state.radiance += state.scatteringFactor * directLight;

        // Sample the BSDF
        float bsdfPdf;
        Pth::sampleBSDF(scene, sampler, state, bsdfPdf);

        // Check if next bounce is to a light source
        Point3f prevSurfaceP = state.intersection.p;      
        state.intersected = scene->rayIntersect(state.ray, state.intersection);
        state.intersectionComputed = true;

        if (state.intersected && state.intersection.mesh->isEmitter())
        {
            const Emitter *emitter = state.intersection.mesh->getEmitter();
            EmitterQueryRecord emitterQuery(prevSurfaceP, state.intersection.p, 
                    state.intersection.toLocal(-state.ray.d), EMeasure::ESolidAngle);


            // Compute Le
            float emitterPdf = emitter->pdf(emitterQuery);
            float weight = Math::powerHeuristic(1, bsdfPdf, nSamplesNES, emitterPdf);
            state.radiance += state.scatteringFactor * weight * emitter->eval(emitterQuery);
            
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

    /*
    void integrateSubsurface(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Sample the subsurface scattering Next Event Estimation
        Color3f Li = Pth::nextEventEstimationBSSRDF(scene, sampler, state, true);
        state.radiance += state.scatteringFactor * Li;
    
        // Sample the subsurface scattering BSDF
        float bsdfPdf;
        Pth::sampleBSSRDF(scene, sampler, state, bsdfPdf);

        // Check if next bounce is to a light source
        Point3f prevSurfaceP = state.intersection.p;      
        state.intersected = scene->rayIntersect(state.ray, state.intersection);
        state.intersectionComputed = true;

        if (state.intersected && state.intersection.mesh->isEmitter())
        {
            const Emitter *emitter = state.intersection.mesh->getEmitter();
            EmitterQueryRecord emitterQuery(prevSurfaceP, state.intersection.p, 
                    state.intersection.toLocal(-state.ray.d), EMeasure::ESolidAngle);


            // Compute Le
            float emitterPdf = emitter->pdf(emitterQuery);
            float weight = Math::powerHeuristic(1, bsdfPdf, 1, emitterPdf);
            state.radiance += state.scatteringFactor * weight * emitter->eval(emitterQuery);
            
            // Terminate the path
            state.scatteringFactor = Color3f(0.0f);
        }
    }
    */

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

            if (state.depth > 3)     // Apply roussian roulette
            {
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
        return "Path_mis[]";
    }

    protected:
        PhotonMap photonMap;
};

NORI_REGISTER_CLASS(Path_mis, "path_mis");
NORI_NAMESPACE_END


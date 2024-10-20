

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
        /*************** Sample outgoing point *******************/
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);
        state.scatteringFactor *= fp;

        /**************** Compute next event estimation ******************/
        float lightPdf, bsdfPdf;
        Color3f directLight = Pth::nextEventEstimation(scene, sampler, state, query, lightPdf, bsdfPdf);
        float weightLight = Math::powerHeuristic(1, lightPdf, 1, bsdfPdf);
        state.radiance += state.scatteringFactor * directLight * weightLight;


        /************************* Sample indirect light *******************/
        Color3f f = Pth::sampleBSDF(state, sampler, query, bsdfPdf);   
        state.scatteringFactor *= f;

        // Check if next bounce is to a light source
        Point3f prevSurfaceP = state.intersection.p;      
        state.intersected = scene->rayIntersect(state.ray, state.intersection);
        state.intersectionComputed = true;

        /*********************** Check for double counting *********************/
        if (state.intersected && state.intersection.mesh->isEmitter())
        {
            const Emitter *emitter = state.intersection.mesh->getEmitter();
            EmitterQueryRecord emitterQuery(prevSurfaceP, state.intersection.p, 
                    state.intersection.vtoLocal(-state.ray.d), EMeasure::ESolidAngle);


            // Compute Le
            lightPdf = emitter->pdf(emitterQuery);
            float weightBSDF = Math::powerHeuristic(1, bsdfPdf, 1, lightPdf);
            state.radiance += state.scatteringFactor * weightBSDF * emitter->eval(emitterQuery);
            
            // Terminate the path
            state.scatteringFactor = Color3f(0.0f);
        }
    }


    void integrateSpecular(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Sample the specular BSDF
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        float bsdfPdf;
        Color3f f = Pth::sampleBSDF(state, sampler, query, bsdfPdf);   
        state.scatteringFactor *= (f * fp);
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
            case Pth::SUBSURFACE:
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
                // Render emitter
                EmitterQueryRecord emitterQuery(-state.ray.d, ESolidAngle);
                emitterQuery.lightP = state.ray.d*1e15;

                if (scene->getEnvironmentalEmitter() != nullptr)
                    state.radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery)*state.scatteringFactor;
                state.scatteringFactor = Color3f(0.0f);
                return;
            }
           
            state.intersectionComputed = false;

            if (state.depth > 3)     // Apply roussian roulette
            {
                constexpr float roulettePdf = 0.95f;
                if (sampler->next1D() < roulettePdf) 
                {
                    // Apply roulette pdf
                    state.scatteringFactor /= roulettePdf;
                }
                else 
                { 
                    state.scatteringFactor = Color3f(0.0f);
                    return;
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


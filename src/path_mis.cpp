

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN





class Path_mis : public Integrator {
public:

    Path_mis(const PropertyList &props) {}

    void preprocess(const Scene *scene, Sampler *sampler) 
    {
        photons = Pth::generateSubsurfaceSamples(scene, sampler);

        for (auto &photon : photons)
        {
            precomputeLi(scene, sampler, photon);
        }
    }

    void integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {   
        // Sample the contribution of a random emitter
        Color3f directLight = Pth::nextEventEstimation(scene, sampler, state, true);
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
            float weight = Math::powerHeuristic(1, bsdfPdf, 1, emitterPdf);
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
                Pth::integrateSubsurface(scene, photons, sampler, state);
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
        ph.n = intersection.shFrame.n;

        float nesPdf = 0.5f;
        if (sampler->next1D() < nesPdf)
        {
            Vector3f wi;
            Color3f direct = Pth::nextEventEstimation(scene, sampler, state, wi, false, false);

            ph.d = wi.normalized();
            ph.radiance = direct / (ph.pdf*nesPdf); 
        }
        else // Compute indirect light
        {

        }
    }

    std::string toString() const {
        return "Path_mis[]";
    }

    protected:
        std::vector<Photon> photons;
};

NORI_REGISTER_CLASS(Path_mis, "path_mis");
NORI_NAMESPACE_END


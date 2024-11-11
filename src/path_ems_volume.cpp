

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>


NORI_NAMESPACE_BEGIN


class Path_ems_volume : public Integrator {
public:

    float sigma_s; // Scatter coefficient
    float sigma_t; // Extinction Coefficient
    float g; // Scatter Phase Coefficient (G)
    float helios_coeff; // Helios Coefficient


    Path_ems_volume(const PropertyList &props) 
    {
        sigma_s = props.getFloat("sigma_s", 0.f);
        // sigma_s = props.getFloat("sigma_s", 0.001f);
        sigma_t = props.getFloat("sigma_t", 0.f);

        g = props.getFloat("g", 0.f);
        helios_coeff = props.getFloat("helios_coeff", 0.f);
    }

    void integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {   
        // Sample the contribution of a random emitter
        const BSDF *bsdf = state.intersection.mesh->getBSDF();

        BSDFQueryRecord query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);
        state.scatteringFactor *= fp;

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

            auto phaseFunction = [g = g]
                                    (const Vector3f &wo,
                                    const Vector3f &wi) -> float
                                    {
                                        return (1.0f / (4 * M_PI)) * (1.0f - g * g) / std::pow(1.0f + g * g - 2.0f * g * std::abs(wo.dot(wi)), 1.5f);
                                    };

            /*********** Sample volume point intersection **********************/
            float volumeD = Math::abs(std::log(1 - sampler->next1D()) / sigma_t);
            if (volumeD < state.intersection.t)
            {
                if (sigma_s==0 || helios_coeff==0) {
                    state.scatteringFactor = Color3f(0.0f);
                    return;
                }
                
                Point3f lp = state.ray(volumeD);
                float transmittance = std::exp(-sigma_t * volumeD);
                float density = sigma_t * transmittance;
                float atmos_pdf = density; //(density.x + density.y + density.z) / 3.0f;
                Color3f mediumScattering = transmittance * sigma_s * helios_coeff / atmos_pdf;

                /******************* Direct *********************/
                Vector3f g_wo; Emitter *emitterMesh; float lightPdf;
                Color3f Le = Pth::estimateDirectLight(scene, sampler, lp, lightPdf, g_wo, emitterMesh);

                if (Le != Color3f(0.0f))
                {
                    Color3f direct = Le / g_wo.squaredNorm();
                    float phase = phaseFunction(-state.ray.d, g_wo.normalized());

                    state.radiance += state.scatteringFactor * mediumScattering * phase * direct * 0.5f;
                }

                /******************* Indirect *******************/
                Vector3f newDir = Warp::squareToUniformSphere(sampler->next2D());

                state.ray = Ray3f(lp, newDir);
                float phase = phaseFunction(-state.ray.d, newDir);

                Color3f inScatterFactor = phase * mediumScattering * 0.5f;
                state.scatteringFactor *= inScatterFactor;   
                
                continue;      
            }
            else
            {
                // Multiply transmittance
                state.scatteringFactor *= std::exp(-sigma_t * state.intersection.t);
            }

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
        return "Path_ems_volume[]";
    }
};

NORI_REGISTER_CLASS(Path_ems_volume, "path_ems_volume");
NORI_NAMESPACE_END




#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN



class Path_mats : public Integrator 
{
private:
    float sigma_s; // Scatter coefficient
    float sigma_t; // Extinction Coefficient
    float g; // Scatter Phase Coefficient (G)
    float helios_coeff; // Helios Coefficient

public:

    Path_mats(const PropertyList &props)
    {
        sigma_s = props.getFloat("sigma_s", 0.f);
        // sigma_s = props.getFloat("sigma_s", 0.001f);
        sigma_t = props.getFloat("sigma_t", 0.001f);

        g = props.getFloat("g", -0.1f);
        helios_coeff = props.getFloat("helios_coeff", 100.f);
    }

    void integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {   
        const BSDF *bsdf = state.intersection.mesh->getBSDF();

        auto query = Pth::initBSDFQuery(scene, sampler, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        float bsdfPdf;
        Color3f f = Pth::sampleBSDF(state, sampler, query, bsdfPdf);   
        state.scatteringFactor *= (f * fp);
    }

    void integrateSpecular(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Sample the specular BSDF
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, sampler, state);
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

    void shadeEnvironment(const Scene *scene, PathState &state) const
    {
        if (scene->getEnvironmentalEmitter() != nullptr && state.depth < 2)
        {
            EmitterQueryRecord emitterQuery (-state.ray.d, EMeasure::EDiscrete);
            emitterQuery.lightP = state.ray.d*1e15;
            state.radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery) * state.scatteringFactor;
        }

        state.scatteringFactor = Color3f(0.0f);
    }

    void integrateVolume(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        const auto phaseFunction = [g = g]
                            (const Vector3f &wo,
                            const Vector3f &wi) -> float
        {
            return (1.0f / (4 * M_PI)) * (1.0f - g * g) / std::pow(1.0f + g * g - 2.0f * g * std::abs(wo.dot(wi)), 1.5f);
        };

        Point3f lp = state.ray(state.intersection.t);
        float transmittance = std::exp(-sigma_t * state.intersection.t);
        float density = sigma_t * transmittance;
        float atmos_pdf = density; //(density.x + density.y + density.z) / 3.0f;
        Color3f mediumScattering = transmittance * sigma_s * helios_coeff / atmos_pdf;

        /******************* Indirect *******************/
        Vector3f newDir = Warp::squareToUniformSphere(sampler->next2D());

        state.ray = Ray3f(lp, newDir);
        float phase = phaseFunction(-state.ray.d, newDir);

        Color3f inScatterFactor = phase * mediumScattering;
        state.scatteringFactor *= inScatterFactor;   
    }

    void sampleIntersection(const Scene *scene, Sampler *sampler, PathState &state) const
    {
        Pth::IntegrationType type = Pth::getIntegrationType(state.intersection);     
        
        #ifndef DISABLE_VOLUME_INTEGRATION
        if (sigma_s != 0 || helios_coeff != 0)
        {
            float volume_t = Math::abs(std::log(1 - sampler->next1D()) / sigma_t);
            if (volume_t < state.intersection.t)
            {
                state.intersection.t = volume_t;
                type = Pth::VOLUME;
            }
            else
            {
                // Apply extinction coefficient
                state.scatteringFactor *= std::exp(-sigma_t * state.intersection.t);
            }
        }
        #endif

        switch (type)
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
            case Pth::VOLUME:
                integrateVolume(scene, sampler, state);
                break;

            default:
                break;
        }
    }

    // Compute radiance over a full path
    void Li(const Scene *scene, Sampler *sampler, PathState &state) const 
    {
        while (state.scatteringFactor != Color3f(0.0f) && state.depth < 100000)
        {
            /* Find the surface that is visible in the requested direction */
            if (!scene->rayIntersect(state.ray, state.intersection))
            {
                shadeEnvironment(scene, state);
                return;
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
        return "Path_mats[]";
    }
};

NORI_REGISTER_CLASS(Path_mats, "path_mats");
NORI_NAMESPACE_END


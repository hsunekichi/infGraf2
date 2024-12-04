

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





class LedaIntegrator : public Integrator 
{
private:
    const int N_SSS_NES_SAMPLES = 8;

    float sigma_s; // Scatter coefficient
    float sigma_t; // Extinction Coefficient
    float g, g_2; // Scatter Phase Coefficient (G)
    float helios_coeff; // Helios Coefficient

public:

    LedaIntegrator(const PropertyList &props) 
    {
        sigma_s = props.getFloat("sigma_s", 0.f);
        sigma_t = props.getFloat("sigma_t", 0.001f);

        g = props.getFloat("g", -0.1f);
        g_2 = Math::pow2(g);
        helios_coeff = props.getFloat("helios_coeff", 0.f);
    }

    void integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {   
        /*************** Sample outgoing point *******************/
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, sampler, state);
        Color3f fp = bsdf->samplePoint(query, sampler);
        state.scatteringFactor *= fp;

        /**************** Compute next event estimation ******************/
        Color3f direct (0.0f);

        const size_t nSamples = bsdf->isSubsurfaceScattering() ? N_SSS_NES_SAMPLES : 1;
        state.previous_n_samples = nSamples;

        for (size_t i = 0; i < nSamples; i++)
        {
            float lightPdf;
            Color3f directLight = Pth::nextEventEstimation(scene, sampler, state, query, lightPdf, state.bsdfPdf);
            float weightLight = Math::powerHeuristic(nSamples, lightPdf, 1, state.bsdfPdf);
            direct += directLight * weightLight;
        }
        
        state.radiance += state.scatteringFactor * direct / nSamples;

        /************************* Sample indirect light *******************/
        Color3f f = Pth::sampleBSDF(state, sampler, query, state.bsdfPdf);   
        state.scatteringFactor *= f;

        state.previous_diffuse = true;
        state.prevP = state.intersection.p;
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
        state.previous_diffuse = false;
    }

    Color3f phaseFunction(const Vector3f &wo,
                        const Vector3f &wi) const
    {
        constexpr float normalization = 1.0f / (4 * M_PI);

        float cosTh = Math::absDot(wo, wi);
        float denom = Math::pow(1.0f + g_2 - 2.0f*g*cosTh, 1.5f);

        return normalization * (1.0f - g_2) / denom;
    }

    void integrateVolume(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        Point3f sampleP = state.ray(state.intersection.t);
        
        float transmittance = std::exp(-sigma_t * state.intersection.t);
        float density = sigma_t * transmittance;
        float atmos_pdf = density; //(density.x + density.y + density.z) / 3.0f;

        Color3f mediumScattering = transmittance * sigma_s * helios_coeff / atmos_pdf;

        /******************* Direct *********************/
        Vector3f g_wo; Emitter *emitterMesh; float lightPdf;
        Color3f Le = Pth::estimateDirectLight(scene, sampler, sampleP, lightPdf, g_wo, emitterMesh);

        if (Le != Color3f(0.0f))
        {
            Vector3f wo = g_wo.normalized();
            Color3f direct = Le / g_wo.squaredNorm();
            Color3f phase = phaseFunction(-state.ray.d, wo);

            float phasePdf = Warp::squareToUniformSpherePdf(wo);
            float weightLight = Math::powerHeuristic(1, lightPdf, 1, phasePdf);

            state.radiance += state.scatteringFactor * mediumScattering * phase * direct * weightLight;
        }

        /******************* Indirect *******************/
        Vector3f newDir = Warp::squareToUniformSphere(sampler->next2D());
        float dirPdf = Warp::squareToUniformSpherePdf(newDir);

        state.ray = Ray3f(sampleP, newDir);
        Color3f phase = phaseFunction(-state.ray.d, newDir);

        Color3f inScatterFactor = phase * mediumScattering;
        state.scatteringFactor *= inScatterFactor; 

        state.bsdfPdf = dirPdf;
        state.previous_diffuse = true;
    }

    void integrateEmitter(const Scene *scene, 
                Sampler *sampler,
                PathState &state) const
    {
        // Retrieve the emitter associated with the surface
        const Emitter *emitter = state.intersection.mesh->getEmitter();

        EmitterQueryRecord emitterQuery (state.intersection.vtoLocal(-state.ray.d), EMeasure::EDiscrete);
        emitterQuery.lightP = state.intersection.p;

        float weight = 1.0f;
        if (state.previous_diffuse)
        {
            emitterQuery.surfaceP = state.prevP;
            emitterQuery.measure = EMeasure::ESolidAngle;
            
            float lightPdf = emitter->pdf(emitterQuery);
            size_t nSamples = state.previous_n_samples;
            weight = Math::powerHeuristic(1, state.bsdfPdf, nSamples, lightPdf);
        }

        state.radiance += state.scatteringFactor * emitter->eval(emitterQuery) * weight;

        // Terminate the path
        state.scatteringFactor = Color3f(0.0f);
    }

    bool isLedaMesh(std::string name) const
    {
        // Name contains "leda"
        //if (name.find("noVol") != std::string::npos)
        //    std::cout << name << std::endl;
        
        return name.find("noVol") != std::string::npos;
    }


    void sampleIntersection(const Scene *scene, Sampler *sampler, PathState &state) const
    {
        Pth::IntegrationType integrationType = Pth::getIntegrationType(state.intersection);

        if (!state.intersection.mesh->getBSDF()->volumeRenderAllowed)
            state.disableVolume = true;

        #ifndef DISABLE_VOLUME_INTEGRATION
        if (sigma_s != 0 || helios_coeff != 0)
        {
            float volume_t = Math::abs(std::log(1 - sampler->next1D()) / sigma_t);
            if (volume_t < state.intersection.t && !state.disableVolume)
            {
                state.intersection.t = volume_t;
                integrationType = Pth::VOLUME;
            }
            else
            {
                // Apply extinction coefficient
                state.scatteringFactor *= std::exp(-sigma_t * state.intersection.t);
            }
        }
        #endif

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
            case Pth::VOLUME:
                integrateVolume(scene, sampler, state);
                break;

            default:
                break;
        }
    }

    void shadeEnvironment(const Scene *scene, PathState &state) const
    {
        if (scene->getEnvironmentalEmitter() != nullptr && state.depth < 2)
        {
            EmitterQueryRecord emitterQuery (-state.ray.d, EMeasure::EDiscrete);
            emitterQuery.lightP = state.ray.d*1e15;

            float weight = 1.0f;
            if (state.previous_diffuse)
            {
                emitterQuery.surfaceP = state.prevP;
                emitterQuery.measure = EMeasure::ESolidAngle;
                
                float lightPdf = scene->getEnvironmentalEmitter()->pdf(emitterQuery);
                size_t nSamples = state.previous_n_samples;
                weight = Math::powerHeuristic(1, state.bsdfPdf, nSamples, lightPdf);
            }

            state.radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery)*state.scatteringFactor * weight;
        }

        state.scatteringFactor = Color3f(0.0f);
    }

    // Compute radiance over a full path
    void Li(const Scene *scene, Sampler *sampler, PathState &state) const 
    {
        while (state.scatteringFactor != Color3f(0.0f))
        {
            /* Find the surface that is visible in the requested direction */
            if (!scene->rayIntersect(state.ray, state.intersection))
            {
                shadeEnvironment(scene, state);
                return;
            }
           
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
        return "LedaIntegrator[]";
    }

    protected:
        PhotonMap photonMap;
};

NORI_REGISTER_CLASS(LedaIntegrator, "leda");
NORI_NAMESPACE_END


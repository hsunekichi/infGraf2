

#include <nori/scene.h>
#include <nori/integrator.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>
#include <nori/kdtree.h>

NORI_NAMESPACE_BEGIN

class Whitted : public Integrator {
public:

    Whitted(const PropertyList &props) {
        sigma_s = props.getFloat("sigma_s", 0.f);
        // sigma_s = props.getFloat("sigma_s", 0.001f);
        sigma_t = props.getFloat("sigma_t", 0.001f);
        g = props.getFloat("g", -0.1f);
        helios_coeff = props.getFloat("helios_coeff", 100.f);
        }

    Color3f specularIntegration(const Scene *scene, 
            Sampler *sampler,
            const Ray3f &ray,
            Intersection &its,
            int depth) const
    {
        // Apply roussian roulette
        float roulettePdf = 1.0f;
        if (depth > 3)
        {
            roulettePdf = 0.95f; 
            if (sampler->next1D() > roulettePdf)
                return Color3f(0.0f);
        }

        PathState state;
        state.intersection = its;
        state.ray = ray;
        
        // Sample the BSDF
        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        if (fp == Color3f(0.0f))
            return Color3f(0.0f);

        float pdf;
        Color3f f = Pth::sampleBSDF(state, sampler, query, pdf);

        pdf = roulettePdf;
        depth++;

        // Compute the contribution
        return f * fp * Li(scene, sampler, state.ray, depth);
    }

    Color3f integrateDiffuse(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {   
        PathState state;
        state.intersection = its;
        state.ray = ray;

        const BSDF *bsdf = state.intersection.mesh->getBSDF();
        
        auto query = Pth::initBSDFQuery(scene, state);
        Color3f fp = bsdf->samplePoint(query, sampler);

        return fp * Pth::nextEventEstimation(scene, sampler, state, query);
    }

    Color3f integrateSubsurface(const Scene *scene, 
                Sampler *sampler,
                const Ray3f &ray,
                Intersection &its) const
    {
        int nSamples = 4;
        Color3f radiance = Color3f(0.0f);

        for (int i = 0; i < nSamples ; i++)
        {
            PathState state;
            state.intersection = its;
            state.ray = ray;

            const BSDF *bsdf = state.intersection.mesh->getBSDF();
            
            auto query = Pth::initBSDFQuery(scene, state);
            Color3f fp = bsdf->samplePoint(query, sampler);

            radiance += fp * Pth::nextEventEstimation(scene, sampler, state, query);
        }

        return radiance / nSamples;
    }

    bool atmos_scatter(Sampler *sampler, float isect_dist) const {
        // Probabilidad base de dispersión atmosférica (ajustable)
        float density = sigma_t * std::exp(-sigma_t * isect_dist);
        float scatter_prob = 1.0f - exp(-sigma_s * density * isect_dist);

        // Incrementar la probabilidad de dispersión según la distancia
        float distance_factor = 1.0f - exp(-isect_dist * 0.1f);  // Ejemplo: más dispersión cuanto mayor es la distancia

        // Generar un número aleatorio entre 0 y 1, usando la semilla
        float random_val = sampler->next1D();  // Esta función usa la semilla para generar un valor aleatorio

        // La dispersión ocurre si el valor aleatorio está por debajo de la probabilidad ajustada
        if (random_val < scatter_prob * distance_factor) {
            return true;  // Dispersión ocurre
        }

        return false;  // No ocurre dispersión
    }

    // Method to compute godrays based on the light source and scattering parameters
    Color3f computeGodrays(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection &its) const {
        if (sigma_s==0 || helios_coeff==0)
            return 0.0f;
        Color3f radiance(0.0f);

        // float transmittance = computeTransmittance(its.p, ray.d, its.t);
        float transmittance = std::exp(-sigma_t * its.t);

        // Dispersión volumétrica
        return transmittance * computeInScattering(scene, sampler, ray, its, sigma_s, g) * helios_coeff;
    }

    Color3f computeInScattering(const Scene *scene, Sampler *sampler, const Ray3f &ray, const Intersection &its, float sigma_s, float g) const {
        Color3f inScatter(0.0f);
        int numSamples_dist = 2; // Controla el número de distancias muestreadas
        int numSamples_dir = 2;  // Controla el número de direcciones muestreadas
        float distanceToLight = its.t; // Distancia a la luz
        float stepSize = distanceToLight / numSamples_dist; // Tamaño del paso de integración
        float phaseFunctionNormalization = 1.0f / (4 * M_PI);


        for (float d = 0.0f; d < distanceToLight; d += stepSize) {
            Vector3f point = its.p + ray.d * d;
            Color3f radiance(0.0f);

            for (int i = 0; i < numSamples_dir; ++i) {
                // Muestreo aleatorio de una dirección en la esfera
                Vector3f sampleDir = Warp::squareToUniformSphere(sampler->next2D());
                
                // Henyey-Greenstein phase function
                float cosTheta = ray.d.dot(sampleDir);
                float phaseFunction = phaseFunctionNormalization * (1.0f - g * g) / std::pow(1.0f + g * g - 2.0f * g * cosTheta, 1.5f);
                
                // Consultar la luz entrante en esa dirección
                Ray3f inRay(point, sampleDir);
                
                PathState state;
                state.intersection = its;
                state.ray = inRay;
                auto query = Pth::initBSDFQuery(scene, state);
                // Color3f directLight = Pth::nextEventEstimation(scene, sampler, state, query); // Luz directa en esa dirección
                Color3f directLight = Li(scene, sampler, inRay, 0, false); // Radiancia entrante
                // if (directLight.maxCoeff()>0)
                // {
                //     cerr << "directlight: " << directLight.toString() << endl;
                // }
                
                radiance += phaseFunction * directLight;
            }
            inScatter += radiance / numSamples_dist;
        }

        return sigma_s * inScatter / numSamples_dist;  // Promediamos el muestreo Monte Carlo
    }



    Color3f Li (const Scene *scene, Sampler *sampler,
            const Ray3f &ray,
            int depth, bool godrays=true) const
    {
        Intersection its;
        Color3f radiance(0.0f);

        /* Find the surface that is visible in the requested direction */
        if (!scene->rayIntersect(ray, its)){
            // Render emitter
            EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
            emitterQuery.lightP = ray.d*1e15;

            if (scene->getEnvironmentalEmitter() != nullptr){
                radiance += scene->getEnvironmentalEmitter()->eval(emitterQuery);
                return radiance;
            }else{
                return Color3f(0.0f);
            }

        }

        if (godrays && atmos_scatter(sampler, its.t))
        {
            Color3f godrays_radiance = computeGodrays(scene, sampler, ray, its);
            // cerr << "godrays_radiance: " << godrays_radiance.toString() << endl;
            return godrays_radiance;
        }

        Pth::IntegrationType type = Pth::getIntegrationType(its);        

        switch(type)
        {
            case Pth::EMITTER:
            {
                // Render emitter
                EmitterQueryRecord emitterQuery(-ray.d, EDiscrete);
                emitterQuery.lightP = its.p;
                radiance += its.mesh->getEmitter()->eval(emitterQuery);
                break;
            }
            case Pth::DIFFUSE:
                // Render diffuse surface
                radiance += integrateDiffuse(scene, sampler, ray, its);
                break;

            case Pth::SUBSURFACE:
                radiance += integrateDiffuse(scene, sampler, ray, its);
                break;

            case Pth::SPECULAR:
                radiance += specularIntegration(scene, sampler, ray, its, depth);
                break;
            
            default:
                break;
        }

        return radiance;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {        
        return Li(scene, sampler, ray, 0);
    }

    std::string toString() const {
        return "Whitted[]";
    }

    protected:
        PhotonMap photonMap;
    private:
        float sigma_s; // Scatter coefficient
        float sigma_t; // Extinction Coefficient
        float g; // Scatter Phase Coefficient (G) 
        float helios_coeff; // Helios Coefficient

};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END


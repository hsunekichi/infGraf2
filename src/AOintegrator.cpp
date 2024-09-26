

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <chrono>
#include <thread>

NORI_NAMESPACE_BEGIN

class AmbientOclusion : public Integrator {
public:

    AmbientOclusion(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        // Generate random ray
        Point2f sample = sampler->next2D();
        Vector3f direction = Warp::squareToCosineHemisphere(sample);
        Ray3f w_sampleRay = Ray3f(its.p, its.shFrame.toWorld(direction));
        
        // Intersect shadow ray
        float visibility = 1.0f;

        if (scene->rayIntersect(w_sampleRay))
            visibility = 0.0f;

        //std::this_thread::sleep_for(std::chrono::milliseconds(1));

        return Color3f(visibility) * Math::dot(Vector3f(its.shFrame.n), w_sampleRay.d) * INV_PI;
    }

    std::string toString() const {
        return "AmbientOclusion[]";
    }
};

NORI_REGISTER_CLASS(AmbientOclusion, "ao");
NORI_NAMESPACE_END


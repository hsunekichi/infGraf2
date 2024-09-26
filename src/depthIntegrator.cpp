

#include <nori/scene.h>
#include <nori/integrator.h>
#include <nori/math.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN

class depthIntegrator : public Integrator {
public:

    depthIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        /* Find the surface that is visible in the requested direction */
        Intersection intersection;
        if (!scene->rayIntersect(ray, intersection))
            return Color3f(0.0f);

        float depth = intersection.t;

        return Color3f(1 / depth);
    }

    std::string toString() const {
        return "depthIntegrator[]";
    }
};

NORI_REGISTER_CLASS(depthIntegrator, "depth");
NORI_NAMESPACE_END

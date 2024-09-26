

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/math.h>
#include <nori/lightProbeMipMap.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:

    Point3f m_position;
    Color3f m_energy;

    SimpleIntegrator(const PropertyList &props) 
    {
        // Get position and energy
        m_position = props.getPoint("position");
        m_energy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        float visibility = 1.0f;

        // Intersect shadow ray
        Ray3f shadowRay(its.p, m_position - its.p, Epsilon, (m_position - its.p).norm() - Epsilon);
        visibility = scene->rayIntersect(shadowRay) ? 0.0f : 1.0f;

        /* Compute the distance between the shading point and the light */
        Vector3f wi = m_position - its.p;
        float d = wi.norm();

        /* Compute the energy received from the light */
        float G = std::max(0.0f, Math::cos(wi, its.shFrame.n));

        return (m_energy * G / (4 * M_PI * M_PI * d * d)) * visibility;
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END


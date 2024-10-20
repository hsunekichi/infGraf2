
#include <nori/common.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PointLight : public Emitter 
{
protected:
    Point3f m_point;

public:
    PointLight(const PropertyList &props) 
    {
        m_type = EmitterType::EMITTER_POINT;
        m_radiance = props.getColor("radiance");
        m_point = props.getPoint("position");
    }

    Color3f sampleLi(Sampler *sampler, EmitterQueryRecord &query) const
    {
        query.pdf = 1.0f;
        query.meshId = -1;
        query.lightP = m_point;
        query.wo = (query.surfaceP - query.lightP).normalized();

        return m_radiance;
    }

    void sampleRay(Sampler *sampler, EmitterQueryRecord &query) const
    {
        query.lightP = m_point;
        query.meshId = -1;

        query.wo = Warp::squareToUniformSphere(sampler->next2D());
        query.pdf = Warp::squareToUniformSpherePdf(query.wo);
    }

    float pdf(const EmitterQueryRecord &query) const
    {
        return 0.0f;
    }

    // Evaluate the emitted radiance on a surface point
    Color3f eval(const EmitterQueryRecord &query) const 
    {
        return m_radiance; 
    }

    std::string toString() const 
    {
        return "PointLight[]";
    }
};

NORI_REGISTER_CLASS(PointLight, "pointlight");
NORI_NAMESPACE_END
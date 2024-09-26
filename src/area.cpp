
#include <nori/common.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter 
{
protected:
    Color3f m_radiance; 


    Point3f samplePoint (Sampler *sampler, Normal3f &n, int &index, float &pdf) const
    {
        if (m_mesh == nullptr)
            throw NoriException("AreaLight::sample(): mesh is nullptr");
        
        uint32_t uint_index;
        Point3f lightP; 
        Point2f uv;

        m_mesh->samplePosition(sampler, lightP, n, uv, pdf, uint_index);
        index = (int)uint_index;

        return lightP;
    }


public:
    AreaLight(const PropertyList &props) 
    {
        m_type = EmitterType::EMITTER_AREA;
        m_radiance = props.getColor("radiance");
    }

    // Returns radiance from a random lightP, over the provided surfaceP, and the direction wo
    Color3f sampleLi(Sampler *sampler, EmitterQueryRecord &query) const
    {
        Normal3f normal;
        Color3f Le = Color3f(0.0f);

        query.lightP = samplePoint(sampler, normal, query.meshId, query.pdf);
        query.wo = (query.surfaceP - query.lightP).normalized();
        Le = m_radiance / query.pdf;

        Le *= Math::absCosTheta(query.wo);

        // Change pdf units if necesary
        if (query.measure == EMeasure::ESolidAngle)
        {
            Normal3f n = m_mesh->getNormal(query.meshId, 
                    m_mesh->uvFrom3D(query.meshId, query.lightP));  

            if (Math::sameDirection(n, query.wo))
                n = -n; 

            
            query.pdf = Math::distanceSquared(query.surfaceP, query.lightP) /
                   (Math::absDot(n, query.wo) * m_mesh->surfaceArea(query.meshId));
        }

        return Le;
    }

    void sampleRay(Sampler *sampler, EmitterQueryRecord &query) const
    {
        throw NoriException("AreaLight::sampleRay(): Not implemented");
    }


    float pdf(const EmitterQueryRecord &query) const
    {
        float pdf = 1.0f;
        int index = query.meshId;

        if (index == -1)
            index = m_mesh->getTriangleIndex(query.lightP);
        if (index == -1)
            return 0.0f;
        
    
        if (query.measure == EMeasure::EDiscrete) {
            pdf = m_mesh->pdf(query.lightP, query.meshId);
        }
        else if (query.measure == EMeasure::ESolidAngle)
        {
            Normal3f n = m_mesh->getNormal(index, m_mesh->uvFrom3D(index, query.lightP));      

            if (Math::sameDirection(n, query.wo))
                n = -n;      

            pdf = Math::distanceSquared(query.surfaceP, query.lightP) /
                   (Math::absDot(n, query.wo) * m_mesh->surfaceArea(index));
        }

        return pdf;
    }

    // Evaluate the emitted radiance on a given direction
    Color3f eval(const EmitterQueryRecord &query) const 
    {
        return m_radiance;
    }

    std::string toString() const 
    {
        return "AreaLight[]";
    }
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
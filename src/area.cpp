
#include <nori/common.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>
#include <nori/sampler.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter 
{
protected:

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

    Texture *txt_radiance = nullptr;
    float scale = 1.0f;

public:
    AreaLight(const PropertyList &props) 
    {
        m_type = EmitterType::EMITTER_AREA;

        txt_radiance = new ConstantSpectrumTexture(props.getColor("radiance", Color3f(0.5f)));
        m_radiance = txt_radiance->eval(Point2f(0.0f));

        scale = props.getFloat("scale", 1.0f);
    }

    ~AreaLight() 
    {
        delete txt_radiance;
    }

    // Returns radiance from a random lightP, over the provided surfaceP, and the direction wo
    Color3f sampleLi(Sampler *sampler, EmitterQueryRecord &query) const
    {
        Normal3f normal;
        Color3f Le = Color3f(0.0f);

        query.lightP = samplePoint(sampler, normal, query.meshId, query.pdf);
        query.wo = (query.surfaceP - query.lightP).normalized();
        Le = txt_radiance->eval(query.uv) * scale / query.pdf;

        if (Math::cos(query.wo, normal) <= 0.0f)
            return Color3f(0.0f);

        Le *= Math::absDot(query.wo, normal);

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
        if (query.wo.z() <= 0.0f)
            return Color3f(0.0f);

        return txt_radiance->eval(query.uv) * scale;
    }

    std::string toString() const 
    {
        return "AreaLight[]";
    }

    void addChild(NoriObject* obj, const std::string& name = "none") 
    {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "radiance")
            {
                if (txt_radiance != nullptr)
                    delete txt_radiance;
                
                txt_radiance = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("Diffuse::addChild(<%s>,%s) is not supported!",
                classTypeName(obj->getClassType()), name);
            break;

        default:
            throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
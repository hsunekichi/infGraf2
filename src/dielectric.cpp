/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/common.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF 
{
    Color3f kd = Color3f(1.0f);
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &info, Sampler *sampler) const 
    {
        float reflectionFr = fresnel(Frame::cosTheta(info.wi), m_extIOR, m_intIOR);
        float refractionFr = 1 - reflectionFr;

        float reflectionPdf = reflectionFr / (reflectionFr + refractionFr);
        float refractionPdf = refractionFr / (reflectionFr + refractionFr);

        Point2f sample = sampler->next2D();

        // Sample proportional to the fresnel term contribution
        if (sample.x() < reflectionPdf)
        {
            info.wo = Vector3f(-info.wi.x(), -info.wi.y(), info.wi.z());
            return kd * reflectionFr / reflectionPdf;
        }
        else
        {
            float originRefrIndex = m_extIOR;
            float destinyRefrIndex = m_intIOR;

            // Check if the ray is entering or exiting the medium
            if (Frame::cosTheta(info.wi) < 0)
                std::swap(originRefrIndex, destinyRefrIndex);
            
            float refractionIndex = originRefrIndex / destinyRefrIndex;
            Vector3f normal = Math::faceForward(Vector3f(0, 0, 1), info.wi);

            // Compute the refraction direction vector with snells law
            float cosThetaOrigin = Math::cos(normal, info.wi);
            float sin2ThetaOrigin = std::max(0.f, 1.f - cosThetaOrigin * cosThetaOrigin);
            float sin2ThetaDestiny = refractionIndex * refractionIndex * sin2ThetaOrigin;

            // Check if the ray is totaly reflected
            if (sin2ThetaDestiny >= 1)
                return Color3f(0);

            float cosThetaT = sqrt(1 - sin2ThetaDestiny);
            info.wo = refractionIndex * -info.wi + 
                (refractionIndex * cosThetaOrigin - cosThetaT) * normal;

            /************ Apply fresnel term for transmited radiance *************/
            return kd * refractionFr / refractionPdf;
        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END

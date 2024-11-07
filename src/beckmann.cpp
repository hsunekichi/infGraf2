/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/math.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

class Beckmann : public BSDF {
public:
    Beckmann(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        float roughness = propList.getFloat("roughness", -1.0f);
        if (roughness >= 0.0f)
            m_alpha = roughness * roughness;

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        internIOR = propList.getColor("internIOR", Color3f(1.0f));
        K = propList.getColor("K", Color3f(0.0f));


        /* Exterior IOR (default: air) */
        m_extIOR = propList.getColor("extIOR", Color3f(1.000277f));
        ior = internIOR / m_extIOR;

        /* Albedo of the diffuse base material (a.k.a "kd") */
        Color3f kd = propList.getColor("albedo", Color3f(0.5f));
        
        if (kd.maxCoeff() > 1.0f)
            kd /= 255.0f;

        m_kd = new ConstantSpectrumTexture(kd);
    }

    float lambda(const Vector3f &rayDir) const
    {
        float alphax = m_alpha;
        float alphay = m_alpha;
        
        float absTanTheta = abs(Math::tanTheta(rayDir));
        
        if (isinf(absTanTheta)) 
            return 0;

        float sin2Ph = std::abs(Math::sin2Phi(rayDir));
        float cos2Ph = std::abs(Math::cos2Phi(rayDir));

        // Compute alpha direction for rayDir
        float alpha = sqrt(cos2Ph * alphax * alphax +
                                sin2Ph * alphay * alphay);

        float a = 1 / (alpha * absTanTheta);

        if (a >= 1.6)
            return 0;

        return (1 - 1.259 * a + 0.396 * a * a) /
            (3.535 * a + 2.181 * a * a);
    }

    float GeometryTerm(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const 
    {
        // Calculate the G1 term for wi and wo
        float lambda1 = lambda(wi);
        float lambda2 = lambda(wo);

        return 1 / (1 + lambda1 + lambda2);
    }

    Color3f l_fresnel(float cosTheta) const 
    {
        if (K == Color3f(0.0f))
            return Math::schlick(cosTheta, ior);
        else
            return Math::fresnel(cosTheta, m_extIOR, internIOR, K);
    }

    bool effectivelySmooth() const {
        return m_alpha < 0.001f;
    }


    Color3f eval(const BSDFQueryRecord &bRec) const 
    {
        // Handle perfect reflection
        if (effectivelySmooth())
            return Color3f(0.0f); 

        // Get the incident and outgoing directions
        const Vector3f wi = bRec.wi;
        const Vector3f wo = bRec.wo;

        // Cosines of the angles between directions and surface normal (z-axis)
        const float cosThetaI = std::abs(Frame::cosTheta(wi));
        const float cosThetaO = std::abs(Frame::cosTheta(wo));

        // Half-vector calculation
        const Vector3f wh = (wi + wo).normalized();

        // Beckmann normal distribution function D(wh)
        const float D = Warp::squareToBeckmannPdf(wh, m_alpha);

        // Fresnel term F(wh ⋅ wi)
        const Color3f F = l_fresnel(wi.dot(wh));

        // Geometry term G(wi, wo, wh)
        const float G = GeometryTerm(wi, wo, wh);

        // Specular term
        const Color3f Kd = m_kd->eval(bRec.uv);
        const Color3f result = (Kd * D * F * G) / (4.0f * cosThetaI * cosThetaO);

        // Return the sum of diffuse and specular terms
        return result;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {

        // Calculate the half-vector
        const Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Calculate the specular density using the Beckmann distribution
        const float D = Warp::squareToBeckmannPdf(wh, m_alpha);
        const float jacobian = 1.0f / (4.0f * std::abs(bRec.wo.dot(wh)));
        
        // Specular density
        const float PDF = D * jacobian;

        return PDF;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const 
    {
        if (Frame::cosTheta(bRec.wi) <= Epsilon)
            return Color3f(0.0f);

        if (effectivelySmooth()) 
        {
            bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
            Color3f F = l_fresnel(Frame::cosTheta(bRec.wo));

            return F;
        }
            
        bRec.wo = Vector3f(0.0f, 0.0f, -1.0f);

        // Rejection sampling to prevent the corner case 
        //  of the microfacet normal reflecting below the surface
        while (bRec.wo.z() <= 0)
        {
            const Point2f sample = sampler->next2D();

            // Sample a normal according to the Beckmann distribution
            const Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
            
            // Reflect the incident direction in the normal 
            //  to get the outgoing direction
            bRec.wo = Math::reflect(bRec.wi, wh);
        }

        const float PDF = pdf(bRec);
        if (std::abs(PDF) <= Epsilon)
            return Color3f(0.0f);

        return eval(bRec) * Frame::cosTheta(bRec.wo) / PDF;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "albedo")
            {
                delete m_kd;
                m_kd = static_cast<Texture*>(obj);
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

    std::string toString() const {
        return tfm::format(
            "Beckmann[\n"
            "  alpha = %f,\n"
            "  etaI = %f,\n"
            "  etaT = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "]",
            m_alpha,
            internIOR.toString(),
            K.toString(),
            m_extIOR
        );
    }
private:
    float m_alpha;
    Color3f internIOR, K, m_extIOR, ior;
    Texture *m_kd;
};

NORI_REGISTER_CLASS(Beckmann, "beckmann");
NORI_NAMESPACE_END

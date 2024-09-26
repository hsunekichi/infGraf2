/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN

class Torrance : public BSDF {
public:
    Torrance(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        float roughness = propList.getFloat("roughness", -1.0f);
        if (roughness > 0.0f)
            m_alpha = roughness * roughness;

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        internIOR = propList.getColor("internIOR", Color3f(1.0f));
        K = propList.getColor("K", Color3f(0.0f));


        /* Exterior IOR (default: air) */
        m_extIOR = propList.getColor("extIOR", Color3f(1.000277f));

        /* Albedo of the diffuse base material (a.k.a "kd") */
        Kd = propList.getColor("kd", Color3f(0.5f));
    }

    float GeometryTerm(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const {
        // Calculate the G1 term for wi and wo
        float G1_wi = GeometryG1(wi, wh);
        float G1_wo = GeometryG1(wo, wh);
        
        // The G term is the product of the two G1 terms
        return G1_wi * G1_wo;
    }

    float GeometryG1(const Vector3f &wv, const Vector3f &wh) const {
        // Calculate the cosine of the angle between wv and the normal (z-axis)
        float cosThetaV = Frame::cosTheta(wv);

        // Tangent of the angle theta_v
        float tanThetaV = std::sqrt(1.0f - cosThetaV * cosThetaV) / cosThetaV;

        // Calculate the 'b' parameter
        float b = 1 / (m_alpha * tanThetaV);

        // If b is less than 1.6, use the rational formula
        if (b < 1.6f) {
            return (3.535f * b + 2.181f * b * b) / (1.0f + 2.276f * b + 2.577f * b * b);
        }

        // Otherwise, return 1
        return 1.0f;
    }

    Color3f l_fresnel(float cosTheta) const 
    {
        if (K == Color3f(0.0f))
            return Math::fresnel(cosTheta, m_extIOR, internIOR);
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
        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;

        // Cosines of the angles between directions and surface normal (z-axis)
        float cosThetaI = std::abs(Frame::cosTheta(wi));
        float cosThetaO = std::abs(Frame::cosTheta(wo));

        // Half-vector calculation
        Vector3f wh = (wi + wo).normalized();

        // Beckmann normal distribution function D(wh)
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);

        // Fresnel term F(wh ⋅ wi)
        Color3f F = l_fresnel(wi.dot(wh));

        // Geometry term G(wi, wo, wh)
        float G = GeometryTerm(wi, wo, wh);

        // Specular term
        Color3f result = (Kd * D * F * G) / (4.0f * cosThetaI * cosThetaO);

        // Return the sum of diffuse and specular terms
        return result;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {

        // Calculate the half-vector
        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Calculate the specular density using the Beckmann distribution
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);
        float jacobian = 1.0f / (4.0f * std::abs(bRec.wo.dot(wh)));
        
        // Specular density
        float PDF = D * jacobian;

        return PDF;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const 
    {
        if (Frame::cosTheta(bRec.wi) <= Epsilon)
            return Color3f(0.0f);

        if (effectivelySmooth()) 
        {
            bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
            Color3f F = l_fresnel(Frame::cosTheta(bRec.wo));

            return F;
        }
            

        // Sample a normal according to the Beckmann distribution
        Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
        
        // Reflect the incident direction in the normal 
        //  to get the outgoing direction
        bRec.wo = -bRec.wi + 2 * bRec.wi.dot(wh) * wh;

        float PDF = pdf(bRec);
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

    std::string toString() const {
        return tfm::format(
            "Torrance[\n"
            "  alpha = %f,\n"
            "  etaI = %f,\n"
            "  etaT = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "]",
            m_alpha,
            internIOR.toString(),
            K.toString(),
            m_extIOR,
            Kd.toString()
        );
    }
private:
    float m_alpha;
    Color3f internIOR, K, m_extIOR;
    Color3f Kd;
};

NORI_REGISTER_CLASS(Torrance, "torrance");
NORI_NAMESPACE_END

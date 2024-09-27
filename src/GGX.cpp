/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN

class GGX : public BSDF {
public:
    GGX(const PropertyList &propList) 
    {
        m_alpha = propList.getFloat("alpha", 0.1f);

        float roughness = propList.getFloat("roughness", -1.0f);
        if (roughness >= 0.0f)
            m_alpha = roughness * roughness;


        /* Interior IOR (default: BK7 borosilicate optical glass) */
        internIOR = propList.getColor("internIOR", Color3f(1.0f));
        K = propList.getColor("K", Color3f(0.0f));

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getColor("extIOR", Color3f(1.000277f));

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));
    }

    float GeometryTerm(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const 
    {
        // Calculate the G1 term for wi and wo
        float G1_wi = GeometryG1(wi, wh);
        float G1_wo = GeometryG1(wo, wh);
        
        // The G term is the product of the two G1 terms
        return G1_wi * G1_wo;
    }

    float GeometryG1(const Vector3f &wv, const Vector3f &wh) const 
    {
        if (m_alpha <= 0)
            return 0.0f;

        float mask = wv.dot(wh) / Frame::cosTheta(wv);
        if (mask < 0)
            return 0.0f;

        float alpha2 = Math::pow2(m_alpha);
        float tanTheta2 = Math::tanTheta2(wv);

        float denom = 1.0f + Math::sqrt(1.0f + alpha2*tanTheta2);
        float G1 = 2.0f / denom;

        return G1;
    }

    float distribution(const Vector3f &wh) const 
    {
        if (m_alpha <= 0 || Math::cosTheta(wh) < 0)
            return 0.0f;

        float numerator = Math::pow2(m_alpha);

        float denom1 = M_PI * Math::pow4(Math::cosTheta(wh));
        float denom2 = Math::pow2(Math::pow2(m_alpha) + Math::tanTheta2(wh));

        return numerator / (denom1*denom2);
    }

    Color3f l_fresnel(float cosTheta) const 
    {
        if (K == Color3f(0.0f))
            return Math::fresnel(cosTheta, m_extIOR, internIOR);
        else
            return Math::fresnel(cosTheta, m_extIOR, internIOR, K);

        float eta_i = m_extIOR.maxCoeff();
        float eta_t = internIOR.maxCoeff();

        // Compute c
        float c = std::abs(cosTheta);

        // Compute g
        float g = std::sqrt(Math::pow2(eta_t) / Math::pow2(eta_i) - 1 + Math::pow2(c));

        // Compute the formula
        float term1 = Math::pow2(g - c) / Math::pow2(g + c);
        float term2 = (Math::pow2((c * (g + c) - 1))) / Math::pow2((c * (g - c) + 1));
        float F = 0.5 * term1 * (1 + term2);

        return F;
    }

    bool effectivelySmooth() const {
        return m_alpha < 0.001f;
    }

    Color3f diffuseTerm(const BSDFQueryRecord &bRec) const
    {
        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;

        // Cosines of the angles between directions and surface normal (z-axis)
        float cosThetaI = Math::absCosTheta(wi);
        float cosThetaO = Math::absCosTheta(wo);
        
        // Half-vector calculation
        Vector3f wh = (wi + wo).normalized(); //  * Math::sign(cosThetaI)


        float D = distribution(wh);
        Color3f F = l_fresnel(wi.dot(wh));
        float G = GeometryTerm(wi, wo, wh);

        // Reflect term
        Color3f diffuse = (m_kd * D * G * F) / (4.0f * cosThetaI * cosThetaO );
        
        return diffuse;
    }

    Color3f specularTerm(const BSDFQueryRecord &bRec) const
    {
        // Get the incident and outgoing directions
        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;

        // Cosines of the angles between directions and surface normal (z-axis)
        float cosThetaI = Math::abs(Frame::cosTheta(wi));
        float cosThetaO = Math::abs(Frame::cosTheta(wo));        

        float etaI = m_extIOR.maxCoeff();
        float etaT = internIOR.maxCoeff();
        float etaT2 = Math::pow2(etaT);

        /************ Specular term ***********/

        // Specular direction
        Vector3f ht = -(etaI * wi + etaT * wo);
        ht.normalize();

        float coshtI = Math::absDot(ht, wi);
        float coshtO = Math::absDot(ht, wo);

        Color3f F = l_fresnel(coshtI);
        float G = GeometryTerm(wi, wo, ht);
        float D = distribution(ht);

        float term1 = coshtI * coshtO / (cosThetaI * cosThetaO);
        Color3f term2Num = etaT2 * (1 - F) * G * D;
        Color3f term2Denom = Math::pow2(etaI * coshtI + etaT * coshtO);

        Color3f specular = m_kd * term1 * term2Num / term2Denom;

        // Return the sum of diffuse and specular terms
        return specular;
    }

    
    Color3f eval(const BSDFQueryRecord &bRec) const 
    {
        if (Frame::cosTheta(bRec.wi) < 0)
            return Color3f(0.0f);

        // Handle perfect reflection
        if (effectivelySmooth())
            return Color3f(0.0f); 
        
        
        return diffuseTerm(bRec);
    }


    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const 
    {
        if (Frame::cosTheta(bRec.wo) < 0)
            return 0.0f;

        if (effectivelySmooth())
            return 0.0f;

        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        float D = GGXpdf(wh);
        float jacobian = 1.0f / (4.0f * std::abs(bRec.wo.dot(wh)));

        return D * jacobian;
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
            Point2f sample = sampler->next2D();

            // Sample a normal according to the GGX distribution
            Vector3f wh = sampleGGX(bRec.wi, sample);
            bRec.wo = Math::reflect(bRec.wi, wh);
        }

        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    float GGXpdf (const Vector3f &wh) const
    {
        return distribution(wh) * Math::absCosTheta(wh);
    }


    Vector3f sampleGGX(const Vector3f &wi, const Point2f &sample) const
    {
        float theta = atan(m_alpha * sqrt(sample.x()) / sqrt(1.0 - sample.x())); // Distribution of theta
        float phi = 2.0 * M_PI * sample.y();                            // Uniform sampling of phi

        // Convert spherical coordinates to Cartesian coordinates
        float sin_theta = sin(theta);
        float cos_theta = cos(theta);

        float x = sin_theta * cos(phi);
        float y = sin_theta * sin(phi);
        float z = cos_theta;

        // Return the sampled microfacet normal
        return Vector3f(x, y, z).normalized();
    }


    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }



    std::string toString() const {
        return tfm::format(
            "GGX[\n"
            "  alpha = %f,\n"
            "  etaI = %f,\n"
            "  etaT = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            internIOR.toString(),
            K.toString(),
            m_extIOR,
            m_kd.toString()
            );
    }
private:
    float m_alpha, roughness;
    Color3f internIOR, K, m_extIOR;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(GGX, "GGX");
NORI_NAMESPACE_END

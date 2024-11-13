/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
	
	v1 - Dec 01 2020
    v2 - Oct 30 2021
	Copyright (c) 2021 by Adrian Jarabo

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/reflectance.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

#define KS_THRES 0.


class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
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

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("albedo", Color3f(0.5f));

        if (m_kd.maxCoeff() > 1.0f)
            m_kd /= 255.0f;

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
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
        {
            return Color3f(0.0f); 
        }

        // Get the incident and outgoing directions
        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;

        // Cosines of the angles between directions and surface normal (z-axis)
        float cosThetaI = Frame::cosTheta(wi);
        float cosThetaO = Frame::cosTheta(wo);

        // Half-vector calculation
        Vector3f wh = (wi + wo).normalized();
        float cosThetaH = Frame::cosTheta(wh);

        // Diffuse term
        Color3f diffuse = m_kd / M_PI;

        // Beckmann normal distribution function D(wh)
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);

        // Fresnel term F(wh ⋅ wi)
        Color3f F = l_fresnel(wi.dot(wh));

        // Geometry term G(wi, wo, wh)
        float G = GeometryTerm(wi, wo, wh);

        // Specular term
        Color3f specular = (m_ks * D * F * G) / (4.0f * cosThetaI * cosThetaO * cosThetaH);

        Color3f f = diffuse + specular;

        // Return the sum of diffuse and specular terms
        return f;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {

        // Calculate the half-vector
        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Calculate the specular density using the Beckmann distribution
        float D = Warp::squareToBeckmannPdf(wh, m_alpha);
        float jacobian = 1.0f / (4.0f * bRec.wo.dot(wh));
        
        // Specular density
        float s_pdf = m_ks * D * jacobian;

        // If diffuse, use the cosine-weighted density
        float d_pdf = (1.0f - m_ks) * Warp::squareToCosineHemispherePdf(bRec.wo);
        return s_pdf + d_pdf;
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
            //F /= absCosine(sampleDirection);

            return F;
        }

        Point2f sample = sampler->next2D();
            
        // Decide whether to sample diffuse or specular reflection
        if (sample.x() <= m_ks) 
        {
            // Specular case: adjust sample.x() for reuse
            Point2f newSample = Point2f(sample.x() / m_ks, sample.y());

            // Sample a normal according to the Beckmann distribution
            Vector3f wh = Warp::squareToBeckmann(newSample, m_alpha);
            
            // Reflect the incident direction in the normal to get the outgoing direction
            bRec.wo = Math::reflect(bRec.wi, wh);

            if (Frame::cosTheta(bRec.wo) <= 0)
                return Color3f(0.0f);
        } 
        else 
        {

            // Diffuse case: adjust sample.x() for cosine-weighted sampling
            Point2f newSample = Point2f((sample.x() - m_ks) / (1.0f - m_ks), sample.y());

            // Sample an outgoing direction using the cosine-weighted distribution
            bRec.wo = Warp::squareToCosineHemisphere(newSample);
        }

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
            "Microfacet[\n"
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
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    Color3f internIOR, K, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

class RoughConductor : public BSDF {
public:
    RoughConductor(const PropertyList& propList) {
        /* RMS surface roughness */
        m_alpha = new ConstantSpectrumTexture(propList.getFloat("alpha", 0.1f));

        /* Reflectance at direction of normal incidence.
           To be used when defining the Fresnel term using the Schlick's approximation*/
        m_R0 = new ConstantSpectrumTexture(propList.getColor("R0", Color3f(0.5f)));
    }

    float GeometryTerm(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const {
        float G1_wi = GeometryG1(wi, wh);
        float G1_wo = GeometryG1(wo, wh);
        return G1_wi * G1_wo;
    }

    float GeometryG1(const Vector3f &wv, const Vector3f &wh) const {
        float cosThetaV = Frame::cosTheta(wv);
        float tanThetaV = std::sqrt(1.0f - cosThetaV * cosThetaV) / cosThetaV;
        float b = 1 / (m_alpha->eval(Point2f(0.0f)).maxCoeff() * tanThetaV);

        if (b < 1.6f) {
            return (3.535f * b + 2.181f * b * b) / (1.0f + 2.276f * b + 2.577f * b * b);
        }
        return 1.0f;
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
        is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float D = Warp::squareToBeckmannPdf(wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        Color3f F = Reflectance::fresnel(Frame::cosTheta(bRec.wi), m_R0->eval(Point2f(0.0f)));
        float G = GeometryTerm(bRec.wi, bRec.wo, wh);

        return D * F * G / (4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo));
    
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
        is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // Calculate the half-vector
        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Calculate the specular density using the Beckmann distribution
        float D = Warp::squareToBeckmannPdf(wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        float jacobian = 1.0f / (4.0f * bRec.wo.dot(wh));

        return D * jacobian;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        if (Frame::cosTheta(bRec.wi) <= Epsilon)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
        {
            bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
            Color3f F = Reflectance::fresnel(Frame::cosTheta(bRec.wi), m_R0->eval(Point2f(0.0f)));
            //F /= absCosine(sampleDirection);

            return F;
        }

        Point2f sample = sampler->next2D();
        
        
        Point2f newSample = Point2f(sample.x(), sample.y());

        // Sample a normal according to the Beckmann distribution
        Vector3f wh = Warp::squareToBeckmann(newSample, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        
        // Reflect the incident direction in the normal to get the outgoing direction
        bRec.wo = Math::reflect(bRec.wi, wh);

        if (Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

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

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "R0")
            {
                delete m_R0;
                m_R0 = static_cast<Texture*>(obj);
            }
            else if (name == "alpha")
            {
                delete m_alpha;
                m_alpha = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("RoughConductor::addChild(<%s>,%s) is not supported!",
                    classTypeName(obj->getClassType()), name);
            break;
        default:
            throw NoriException("RoughConductor::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "RoughConductor[\n"
            "  alpha = %f,\n"
            "  R0 = %s,\n"
            "]",
            m_alpha->toString(),
            m_R0->toString()
        );
    }
private:
    Texture* m_alpha;
    Texture* m_R0;
};


class RoughDielectric : public BSDF {
public:
    RoughDielectric(const PropertyList& propList) {
        /* RMS surface roughness */
        m_alpha = new ConstantSpectrumTexture(propList.getFloat("alpha", 0.1f));

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Tint of the glass, modeling its color */
        m_ka = new ConstantSpectrumTexture(propList.getColor("ka", Color3f(1.f)));
    }


    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BSDF -- return zero if the measure is wrong */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);


        throw NoriException("RoughDielectric::eval() is not yet implemented!");
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const {
        /* This is a smooth BSDF -- return zero if the measure is wrong */
        if (bRec.measure != ESolidAngle)
            return 0.0f;

        throw NoriException("RoughDielectric::eval() is not yet implemented!");
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        bRec.measure = ESolidAngle;

        throw NoriException("RoughDielectric::sample() is not yet implemented!");
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
            if (name == "m_ka")
            {
                delete m_ka;
                m_ka = static_cast<Texture*>(obj);
            }
            else if (name == "alpha")
            {
                delete m_alpha;
                m_alpha = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("RoughDielectric::addChild(<%s>,%s) is not supported!",
                    classTypeName(obj->getClassType()), name);
            break;
        default:
            throw NoriException("RoughDielectric::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "RoughDielectric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  ka = %s,\n"
            "]",
            m_alpha->toString(),
            m_intIOR,
            m_extIOR,
            m_ka->toString()
        );
    }
private:
    float m_intIOR, m_extIOR;
    Texture* m_alpha;
    Texture* m_ka;
};

class RoughSubstrate : public BSDF {
public:
    RoughSubstrate(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = new ConstantSpectrumTexture(propList.getFloat("alpha", 0.1f));

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = new ConstantSpectrumTexture(propList.getColor("kd", Color3f(0.5f)));
    }


    float GeometryTerm(const Vector3f &wi, const Vector3f &wo, const Vector3f &wh) const {
        // Calculate the G1 term for wi and wo
        // float G1_wi = GeometryG1(wi, wh);
        // float G1_wo = GeometryG1(wo, wh);
        float G1_wi = Reflectance::G1(wi, wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        float G1_wo = Reflectance::G1(wo, wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        
        // The G term is the product of the two G1 terms
        return G1_wi * G1_wo;
    }


    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
        is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        // Get the incident and outgoing directions
        Vector3f wi = bRec.wi;
        Vector3f wo = bRec.wo;

        // Cosines of the angles between directions and surface normal (z-axis)
        float cosThetaI = Frame::cosTheta(wi);
        float cosThetaO = Frame::cosTheta(wo);

        // Half-vector calculation
        Vector3f wh = (wi + wo).normalized();
        float cosThetaH = Frame::cosTheta(wh);

        // Diffuse term
        // Color3f diffuse = m_kd->eval(Point2f(0.0f)) / M_PI;

        // Cálculo de la diferencia de índices de refracción
        float etaRatio = (m_extIOR - m_intIOR) / (m_extIOR + m_intIOR);
        float fresnelFactor = 1.0f - etaRatio * etaRatio;

        // Cálculo del término difuso
        float diffFactor = (1.0f - std::pow(1.0f - 0.5f * cosThetaI, 5)) * (1.0f - std::pow(1.0f - 0.5f * cosThetaO, 5));
        Color3f diffuse = (28.0f / (23.0f * M_PI)) * m_kd->eval(Point2f(0.0f)) * fresnelFactor * diffFactor;


        // Beckmann normal distribution function D(wh)
        float D = Warp::squareToBeckmannPdf(wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        // float D = Reflectance::BeckmannNDF(wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        

        // Fresnel term F(wh ⋅ wi)
        Color3f F = Reflectance::fresnel(wi.dot(wh), m_extIOR, m_intIOR);

        // Geometry term G(wi, wo, wh)
        float G = GeometryTerm(wi, wo, wh);

        float m_ks = 1 - m_kd->eval(Point2f(0.0f)).maxCoeff();
        // Specular term
        Color3f specular = (m_ks * D * F * G) / (4.0f * cosThetaI * cosThetaO);// * cosThetaH);

        Color3f f = diffuse + specular*M_PI*2;

        // Return the sum of diffuse and specular terms
        return f;

		throw NoriException("RoughSubstrate::eval() is not yet implemented!");
	}

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
       is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        // Calculate the half-vector
        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        // Calculate the specular density using the Beckmann distribution
        float D = Warp::squareToBeckmannPdf(wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        // float D = Reflectance::BeckmannNDF(wh, m_alpha->eval(Point2f(0.0f)).maxCoeff());
        
        float jacobian = 1.0f / (4.0f * bRec.wo.dot(wh));
        
        float m_ks = 1 - m_kd->eval(Point2f(0.0f)).maxCoeff();

        // Specular density
        float s_pdf = m_ks * D * jacobian;

        // If diffuse, use the cosine-weighted density
        float d_pdf = (1.0f - m_ks) * Warp::squareToCosineHemispherePdf(bRec.wo);
        return s_pdf + d_pdf;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const {
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        if (m_alpha->eval(Point2f(0.0f)).maxCoeff() < 0.001f)
        {
            bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
            Color3f F = Reflectance::fresnel(Frame::cosTheta(bRec.wo), m_extIOR, m_intIOR);
            // F /= Math::absCosTheta(bRec.wo);

            return F;
        }

        Point2f sample = sampler->next2D();

        float m_ks = 1 - m_kd->eval(Point2f(0.0f)).maxCoeff();
            
        // Decide whether to sample diffuse or specular reflection
        if (sample.x() <= m_ks) 
        {
            // Specular case: adjust sample.x() for reuse
            Point2f newSample = Point2f(sample.x() / m_ks, sample.y());

            // Sample a normal according to the Beckmann distribution
            Vector3f wh = Warp::squareToBeckmann(newSample, m_alpha->eval(Point2f(0.0f)).maxCoeff());
            
            // Reflect the incident direction in the normal to get the outgoing direction
            bRec.wo = Math::reflect(bRec.wi, wh);
            // bRec.wo = Reflectance::refract(bRec.wi, wh, m_extIOR, m_intIOR);


            if (Frame::cosTheta(bRec.wo) <= 0)
                return Color3f(0.0f);
        } 
        else 
        {
            // Diffuse case: adjust sample.x() for cosine-weighted sampling
            Point2f newSample = Point2f((sample.x() - m_ks) / (1.0f - m_ks), sample.y());

            // Sample an outgoing direction using the cosine-weighted distribution
            bRec.wo = Warp::squareToCosineHemisphere(newSample);
        }

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

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "kd")
            {
                delete m_kd;
                m_kd = static_cast<Texture*>(obj);
            }
            else if (name == "alpha")
            {
                delete m_alpha;
                m_alpha = static_cast<Texture*>(obj);
            }
            else 
                throw NoriException("RoughSubstrate::addChild(<%s>,%s) is not supported!",
                    classTypeName(obj->getClassType()), name);
            break;
        default:
            throw NoriException("RoughSubstrate::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    std::string toString() const {
        return tfm::format(
            "RoughSubstrate[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "]",
            m_alpha->toString(),
            m_intIOR,
            m_extIOR,
            m_kd->toString()
        );
    }
private:
    float m_intIOR, m_extIOR;
    Texture* m_alpha;
    Texture* m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_REGISTER_CLASS(RoughConductor, "roughconductor");
NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_REGISTER_CLASS(RoughSubstrate, "roughsubstrate");

NORI_NAMESPACE_END

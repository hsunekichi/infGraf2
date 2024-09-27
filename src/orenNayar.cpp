/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN

class OrenNayar : public BSDF {
public:
    OrenNayar(const PropertyList &propList) 
    {
        sigma = propList.getFloat("roughness", 0.1f);
        sigma = Math::toRadians(sigma);

        Kd = propList.getColor("kd", Color3f(0.5f));
        
        float sigma2 = sigma * sigma;
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }

    Color3f eval(const BSDFQueryRecord &bRec) const 
    {
        float sinThetaI = Math::sinTheta(bRec.wi);
        float sinThetaO = Math::sinTheta(bRec.wo);

        // Compute cosine term
        float maxCos = 0;
        if (sinThetaI > 1e-4 && sinThetaO > 1e-4) {
            float sinPhiI = Math::sinPhi(bRec.wi), cosPhiI = Math::cosPhi(bRec.wi);
            float sinPhiO = Math::sinPhi(bRec.wo), cosPhiO = Math::cosPhi(bRec.wo);
            float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
            maxCos = std::max((float)0, dCos);
        }

        // Compute sin and tan terms
        float sinAlpha, tanBeta;
        if (Math::absCosTheta(bRec.wi) > Math::absCosTheta(bRec.wo)) {
            sinAlpha = sinThetaO;
            tanBeta = sinThetaI / Math::absCosTheta(bRec.wi);
        } else {
            sinAlpha = sinThetaI;
            tanBeta = sinThetaO / Math::absCosTheta(bRec.wo);
        }

        return Kd * INV_PI * (A + B * maxCos * sinAlpha * tanBeta);
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const 
    {
        return Warp::squareToCosineHemispherePdf(bRec.wo);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const 
    {
        Point2f sample = sampler->next2D();
        bRec.wo = Warp::squareToCosineHemisphere(sample);
        return eval(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "OrenNayar[\n"
            "  sigma = %f,\n"
            "  A = %f,\n"
            "  B = %f,\n"
            "  Kd = %s\n"
            "]",
            sigma,
            A,
            B,
            Kd.toString()
        );
    }
private:
    float sigma, A, B;
    Color3f Kd;
};

NORI_REGISTER_CLASS(OrenNayar, "orenNayar");
NORI_NAMESPACE_END

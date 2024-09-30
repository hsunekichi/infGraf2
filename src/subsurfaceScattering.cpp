/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    v1 - Dec 01 2020
    Copyright (c) 2020 by Adrian Jarabo

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
#include <nori/texture.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief SubsurfaceScattering / Lambertian BRDF model
 */
class SubsurfaceScattering : public BSDF 
{
public:
    SubsurfaceScattering(const PropertyList &propList) 
    {
        m_albedo = new ConstantSpectrumTexture(propList.getColor("albedo", Color3f(0.5f)));

        Color3f kd = propList.getColor("kd", Color3f(-1.0f));
        if (kd != Color3f(-1.0f))
        {
            delete m_albedo;
            m_albedo = new ConstantSpectrumTexture(kd);
        }

        sigmaA = propList.getColor("sigmaA", Color3f(0.0f));
        sigmaS = propList.getColor("sigmaS", Color3f(0.0f));
        Color3f sigmaT = propList.getColor("sigmaT", Color3f(-1.0f));

        g = propList.getFloat("g", 0.0f);
        etaT = propList.getFloat("eta", 1.0f);

        if (sigmaT != Color3f(-1.0f))
        {
            sigmaS = sigmaT * (1 - g);
            sigmaA = sigmaT - sigmaS;
        }
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        if (!bRec.isCameraRay) {
            return m_albedo->eval(bRec.uv) * INV_PI;
        }
        else
        {
            return computeMultipleScattering(bRec);
        }
    }

  

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return 0.0f;


        /* Importance sampling density wrt. solid angles:
           cos(theta) / pi.

           Note that the directions in 'bRec' are in local coordinates,
           so Frame::cosTheta() actually just returns the 'z' component.
        */
        return INV_PI * Frame::cosTheta(bRec.wo);
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const 
    {
        bRec.measure = ESolidAngle;

        Point2f sample = sampler->next2D();

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        bRec.wo = Warp::squareToCosineHemisphere(sample);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return m_albedo->eval(bRec.uv);
    }

    bool isSubsurfaceScattering() const {
        return true;
    }

    bool isDiffuse() const {
        return true;
    }


    Color3f dipoleDiffusionAproximation(float r) const
    {
        float eta = etaT;

        // Compute isotropic phase function
        Color3f _sigmaS = (1 - g) * sigmaS;
        Color3f _sigmaT = _sigmaS + sigmaA;
        Color3f _alpha = _sigmaS / _sigmaT;

        // Effective transport coefficient
        Color3f sigmaTr = Math::sqrt(3 * sigmaA * _sigmaT);
        
        // Aproximation for the diffuse reflectance (fresnel)
        float Fdr = (-1.440 / Math::pow2(eta)) + (0.710 / eta) + 0.668 + 0.0636 * eta;
        float A = (1 + Fdr) / (1 - Fdr);    // Boundary condition for the change between refraction indexes

        Color3f lu = 1 / _sigmaT;
        Color3f zr = lu;
        Color3f zv = lu * (1.0f + 4.0f/3.0f * A);

        Color3f dr = Math::sqrt(Math::pow2(r) + Math::pow2(zr)); 
        Color3f dv = Math::sqrt(Math::pow2(r) + Math::pow2(zv)); 

        // Compute main formula
        Color3f C1 = zr * (1 + sigmaTr * dr) * Math::exp(-sigmaTr * dr) / Math::pow3(dr);
        Color3f C2 = zv * (1 + sigmaTr * dv) * Math::exp(-sigmaTr * dv) / Math::pow3(dv);

        Color3f result = _alpha * INV_FOURPI * (C1 - C2);

        return result;
    }

    Color3f impPaper_dipoleDiffusionAproximation(float r) const
    {
        float eta = etaT;

        // Compute isotropic phase function
        Color3f _sigmaS = (1 - g) * sigmaS;
        Color3f _sigmaT = _sigmaS + sigmaA;
        Color3f _alpha = _sigmaS / _sigmaT;

        // Effective transport coefficient
        Color3f sigmaTr = Math::sqrt(3 * sigmaA * _sigmaT);
        
        // Aproximation for the diffuse reflectance (fresnel)
        float Fdr = (-1.440 / Math::pow2(eta)) + (0.710 / eta) + 0.668 + 0.0636 * eta;
        float A = (1 + Fdr) / (1 - Fdr);    // Boundary condition for the change between refraction indexes

        Color3f zr = sqrt(3*(1-_alpha)) / sigmaTr;
        Color3f zv = zr * A;

        Color3f distanceR = sqrt(Math::pow2(r) + Math::pow2(zr)); 
        Color3f distanceV = sqrt(Math::pow2(r) + Math::pow2(zv)); 
        Color3f sigmaTrDr = sigmaTr * distanceR;
        Color3f sigmaTrDv = sigmaTr * distanceV;

        Color3f Rd = (sigmaTrDr+1) * exp(-sigmaTrDr) * zr/Math::pow3(distanceR)
                        +
                        (sigmaTrDv+1) * exp(-sigmaTrDv) * zv/Math::pow3(distanceV);

        Color3f C1 = zr * (sigmaTr + 1/distanceR);

        Color3f result = C1 * Rd * (1-Fdr) * _alpha;

        return result;
    }


    Color3f computeMultipleScattering(const BSDFQueryRecord &bRec) const
    {
        float r = (bRec.po - bRec.pi).norm();
        float eta = 1.0f / etaT;

        Color3f Rd = dipoleDiffusionAproximation(r);

        float cosWi = Math::absCos(bRec.wi, bRec.ni);
        float cosWo = Math::cosTheta(bRec.wo);

        double Ft_o = 1 - fresnel(cosWo, 1.0, eta);
        double Ft_i = 1 - fresnel(cosWi, 1.0, eta);

        // Compute the diffusion term
        Color3f Sd = INV_PI * Ft_i * Rd * Ft_o;

        return Sd;
    } 
    

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "SubsurfaceScattering[\n"
            "  albedo = %s\n"
            "]", m_albedo->toString());
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case ETexture:
            if (name == "albedo")
            {
                delete m_albedo;
                m_albedo = static_cast<Texture*>(obj);
            }
            else
                throw NoriException("SubsurfaceScattering::addChild(<%s>,%s) is not supported!",
                classTypeName(obj->getClassType()), name);
            break;

        default:
            throw NoriException("SubsurfaceScattering::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }


    EClassType getClassType() const { return EBSDF; }
private:
    Texture *m_albedo;
    Color3f sigmaA, sigmaS; 
    float g, etaT;
};

NORI_REGISTER_CLASS(SubsurfaceScattering, "subsurface");
NORI_NAMESPACE_END

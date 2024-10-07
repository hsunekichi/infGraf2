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

        sigmaA = propList.getColor("sigmaA", Color3f(0.0f));
        sigmaS = propList.getColor("sigmaS", Color3f(0.0f));
        Color3f sigmaT = propList.getColor("sigmaT", Color3f(-1.0f));
        
        scale = propList.getFloat("scale", 1.0f);
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

        if (!bRec.isCameraRay)
            return m_albedo->eval(bRec.uv) * INV_PI;
        else
            return evalSubsurface(bRec);
    }

    Color3f evalSubsurface(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        return computeMultipleScattering(bRec) 
                * m_albedo->eval(bRec.uv) 
                * Math::absCosTheta (bRec.wo)
                * scale;
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

        bRec.isCameraRay = false;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return eval(bRec);
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
        float r2 = Math::pow2(r);

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
        Color3f zv = -lu * (1.0f + (4.0f/3.0f) * A); // Negative because it is inside the surface

        Color3f dr = Math::sqrt(r2 + Math::pow2(zr)); 
        Color3f dv = Math::sqrt(r2 + Math::pow2(zv)); 

        Color3f sigmaTrDr = sigmaTr * dr;
        Color3f sigmaTrDv = sigmaTr * dv;

        // Compute main formula
        Color3f C1 = zr * (1 + sigmaTrDr) * Math::exp(-sigmaTrDr) / Math::pow3(dr);
        Color3f C2 = zv * (1 + sigmaTrDv) * Math::exp(-sigmaTrDv) / Math::pow3(dv);

        Color3f result = _alpha * INV_FOURPI * (C1 - C2);

        return result;
    }


    Color3f computeMultipleScattering(const BSDFQueryRecord &bRec) const
    {
        float r = (bRec.po - bRec.pi).norm();
        float eta = etaT;

        Color3f Rd = betterAlternative(r);

        float cosWi = Math::absCos(bRec.wi, bRec.ni);
        float cosWo = Math::cosTheta(bRec.wo);

        double Ft_o = 1 - fresnel(cosWo, 1.0, eta);
        double Ft_i = 1 - fresnel(cosWi, 1.0, eta);

        // Compute the diffusion term
        Color3f Sd = INV_PI * Ft_i * Rd * Ft_o;

        Sd = Math::max(Sd, Color3f(0.0f));

        return Sd;
    } 

    double Fdr() const {
        const double eta = etaT;
        const double eta2 = eta * eta;
        const double eta3 = eta2 * eta;
        if (eta >= 1.0) {
            return -1.4399 / eta2 + 0.7099 / eta + 0.6681 +
                    0.0636 * eta;
        } else {
            return -0.4399 + 0.7099 / eta - 0.3319 / eta2 +
                    0.0636 / eta3;
        }
    }

    double C1() const {
        const double eta = etaT;
        const double eta2 = eta * eta;
        const double eta3 = eta2 * eta;
        const double eta4 = eta3 * eta;
        const double eta5 = eta4 * eta;
        if (eta < 1.0) {
            return (0.919317 - 3.4793 * eta + 6.75335 * eta2
                    -7.80989 * eta3 + 4.98554 * eta4 - 1.36881 * eta5) / 2.0;
        } else {
            return (-9.23372 + 22.2272 * eta - 20.9292 * eta2 + 10.2291 * eta3
                    - 2.54396 * eta4 + 0.254913 * eta5) / 2.0;
        }
    }

    double C2() const {
        const double eta = etaT;
        const double eta2 = eta * eta;
        const double eta3 = eta2 * eta;
        const double eta4 = eta3 * eta;
        const double eta5 = eta4 * eta;
        if (eta < 1.0) {
            return (0.828421 - 2.62051 * eta + 3.36231 * eta2 - 1.95284 * eta3
                    + 0.236494 * eta4 + 0.145787 * eta5) / 3.0;
        } else {
            return (-1641.1 + 135.926 / eta3 - 656.175 / eta2 + 1376.53 / eta
                    + 1213.67 * eta - 568.556 * eta2 + 164.798 * eta3
                    - 27.0181 * eta4 + 1.91826 * eta5) / 3.0;
        }
    }

    Color3f betterAlternative (float r) const
    {
        Color3f sigmap_s = (1 - g) * sigmaS;
        Color3f sigmap_t = sigmap_s + sigmaA;
        Color3f alphap = sigmap_s / sigmap_t;
        float A = (1.0 + 3.0 * C2()) / (1.0 - 2.0 * C1());
        Color3f D = (2.0 * sigmaA + sigmap_s) / (3.0 * sigmap_t * sigmap_t);
        Color3f sigma_tr = Math::sqrt(sigmaA / D);
        Color3f zr = 1.0 / sigmap_t;
        Color3f zb = 2.0 * A * D;
        Color3f zv = -zr - 2.0 * zb;
        Color3f Cphi = 0.25 * (1.0 - 2.0 * C1());
        Color3f CE   = 0.50 * (1.0 - 3.0 * C2());

        const double r2 = r * r;
        const Color3f dr = Math::sqrt(r2 + zr * zr);
        const Color3f dv = Math::sqrt(r2 + zv * zv);
        const Color3f ar = zr * (sigma_tr * dr + 1.0) / (dr * dr);
        const Color3f av = zv * (sigma_tr * dv + 1.0) / (dv * dv);
        const Color3f br = Math::exp(-sigma_tr * dr) / dr;
        const Color3f bv = Math::exp(-sigma_tr * dv) / dv;
        return (alphap * alphap / (4.0 * M_PI)) *
               ((CE * ar + Cphi / D) * br -
                (CE * av + Cphi / D) * bv);
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
    float g, etaT, scale;
};

NORI_REGISTER_CLASS(SubsurfaceScattering, "subsurface");
NORI_NAMESPACE_END

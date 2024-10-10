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

#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN


/**
 * \brief SubsurfaceScattering / Lambertian BRDF model
 */
class subsurface : public BSDF
{
public:
    subsurface(const PropertyList &propList) 
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

        return evalSubsurface(bRec);

        //if (!bRec.isCameraRay)
        //    return m_albedo->eval(bRec.uv) * INV_PI;
        //else
        //    return evalSubsurface(bRec);
    }

    Color3f evalSubsurface(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        return computeMultipleScattering(bRec) 
                * m_albedo->eval(bRec.uv)
                / Math::absCos(bRec.fri.n, bRec.fro.n);
    }
  

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return 0.0f;

        return Warp::squareToCosineHemispherePdf(bRec.wo);
    }


    void samplePointOld(BSDFQueryRecord &bRec, Sampler *sampler, float &pdf) const
    {
        //Color3f sigmaT = sigmaA + sigmaS;

        // Select random channel
        //int channel = sampler->next1D() * 3;
        //float sigma = sigmaT[channel];
        //float channelPdf = 1.0f / 3.0f;

        // Sample an offset proportional to the sigmaT
        //Point2f sample = Warp::squareToSquaredDecayDisk(sampler->next2D(), sigma);

        // Clamp to account for highly curved surfaces
        //sample = Math::max(sample, Point2f(1.0f/sigma));
   
        // Sampled offset (in mm) to world
        //Vector3f sampled = Vector3f(sample.x(), sample.y(), 0.0f) / 1000.0f;
        //bRec.po = bRec.pi + bRec.frame.toWorld(sampled); 
        //bRec.po = projectToSurface(bRec, bRec.po);

        //pdf = Warp::squareToSquaredDecayDiskPdf(sample, sigma) * channelPdf;
    }

    bool samplePoint(BSDFQueryRecord &bRec, Sampler *sampler, float &pdf) const
    {
        Color3f sigmaT = sigmaA + sigmaS;

        // Select random channel
        int channel =  sampler->next1D() * 3;
        float sigma = sigmaT[channel];
        float channelPdf = 1.0f / 3.0f;

        // Select sphere radius
        float rMax = Warp::squareToSrDecay(sampler->next1D(), sigma);
        float rMin = 1.0f/sigma;
        float r = Warp::squareToSrDecay(sampler->next1D(), sigma);

        // Clamp to account for highly curved surfaces, change to meters
        rMax = Math::max(rMax, rMin) / 1000.0f;
        r = Math::clamp(r, rMin, rMax) / 1000.0f;

        // Sample a point on the disk
        Point2f diskSample = Warp::SrToDisk(sampler->next1D(), r);
        float lhalf = Math::sqrt(Math::pow2(rMax) - Math::pow2(r)); // * 2 / 2

        // Project to sphere
        Point3f base = Point3f (diskSample.x(), diskSample.y(), lhalf);
        Point3f target = Point3f (diskSample.x(), diskSample.y(), -lhalf);

        // Sphere to world
        Point3f w_base = bRec.fri.ptoWorld(base);
        Point3f w_target = bRec.fri.ptoWorld(target);
 

        /************ Project point to shape ****************/

        Vector3f d = (w_target - w_base).normalized();
        Ray3f ray = Ray3f(w_base, d, 0.0f, lhalf*2.0f);

        Intersection its;
        bool intersected = bRec.scene->rayIntersect(ray, bRec.mesh, its);

        if (!intersected)
        {
            return false;
        }

        bRec.fro = its.shFrame;
        bRec.po = ray(its.t);

        /*********************** Pdf *************************/

        float rMaxPdf = Warp::squareToSrDecayPdf(rMax*1000, sigma);
        float rPdf = Warp::squareToSrDecayPdf(r*1000, sigma);
        float diskPdf = Warp::SrToDiskPdf(diskSample);

        pdf = rMaxPdf * rPdf * diskPdf * channelPdf;

        return true;
    }


    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        bRec.measure = ESolidAngle;
        float pdf;

        bRec.wo = Warp::squareToCosineHemisphere(sampler->next2D());
        pdf = Warp::squareToCosineHemispherePdf(bRec.wo);

        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf;
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
        Color3f _sigmaS = sigmaS;
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
        float r = (bRec.po - bRec.pi).norm() * 1000.0f;
        float eta = etaT;

        Color3f Rd = betterAlternative(r);
        //std::cout << Rd.toString() + ", " + std::to_string(r) << std::endl;

        Vector3f w_wi = bRec.fro.vtoWorld(bRec.wi);
        float cosWi = Math::absCos(w_wi, bRec.fri.n);
        float cosWo = Math::absCosTheta(bRec.wo);

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
        Color3f sigmap_s = sigmaS;
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

NORI_REGISTER_CLASS(subsurface, "subsurface");
NORI_NAMESPACE_END

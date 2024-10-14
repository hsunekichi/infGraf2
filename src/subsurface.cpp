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
        
        scale = propList.getFloat("scale", 1000.0f);
        etaT = propList.getFloat("eta", 1.0f);
        float eta = 1.0f / etaT;

        Color3f sigmaT = (sigmaA + sigmaS);
        
        // Effective transport coefficient
        sigmaTr = Math::sqrt(3 * sigmaA * sigmaT);
        alpha_4pi = (sigmaS / sigmaT) * INV_FOURPI;

        // Constant used on the sampling method
        Color3f D = (2 * sigmaA) / (3.0f * Math::pow2(sigmaT));
        ld = 1.0f / Math::sqrt(sigmaA/D);

        // Aproximation for the diffuse reflectance (fresnel)
        float Fdr = (-1.440 / Math::pow2(eta)) + (0.710 / eta) + 0.668 + 0.0636 * eta;
        float A = (1 + Fdr) / (1 - Fdr);    // Boundary condition for the change between refraction indexes

        Color3f invSigmaT = 1 / sigmaT;
        zr = invSigmaT;
        zv = -invSigmaT * (1.0f + (4.0f/3.0f) * A); // Negative because it is inside the surface
    
        precomputeSrCdfInverse();
    }

    float SrCdf(float r) const
    {
        return 1.0 - 0.25 * Math::exp(-r) - 0.75 * Math::exp(-r/3);
    }

    void precomputeSrCdfInverse()
    {
        size_t steps = 4;
        float x0 = 0.0f;

        SrCdfInverse.resize(steps+1);
        XICdfInverse.resize(steps+1);

        for (size_t i = 0; i < steps; i++)
        {
            float xi = float(i) / steps;
            auto f = [this, xi](float r) { return SrCdf(r) - xi; };
            float r = Math::findRoot(f, x0);

            x0 = r;
            SrCdfInverse[i] = r;
            XICdfInverse[i] = xi;
        }

        // Compute a final step close to 1
        //  to allow interpolation on the last element
        float last_step = 1.0f - (1.0f / steps) / 8.0f;
        auto f = [this, last_step](float r) { return SrCdf(r) - last_step; };
        float r = Math::findRoot(f, x0);
        SrCdfInverse[steps] = r;
        XICdfInverse[steps] = last_step;
    }

    Color3f Sr_(float r) const
    {        
        // THIS FUNCTION HAS A SINGULARITY ON r=0????????

        r *= 1000.0f; // Convert to mm
        //Color3f s = 3.5f + 100.0f * Math::pow4(R - 0.33f);
        Color3f d = ld;

        Color3f num1 = Math::exp(-r / d);
        Color3f num2 = Math::exp(-r / (3 * d));
        Color3f den = 8.0f * M_PI * d * r;
        
        Color3f result = (num1 + num2) / den;
    
        return result;
    }

    float sampleSr(float rnd, int channel) const
    {
        size_t steps = SrCdfInverse.size() - 1;
        int i_xi = Math::floor(rnd * steps);

        float r1 = SrCdfInverse[i_xi];
        float r2 = SrCdfInverse[i_xi + 1];

        float x1 = XICdfInverse[i_xi];
        float x2 = XICdfInverse[i_xi + 1];

        float r = Math::lerp(r1, r2, (rnd - x1) / (x2 - x1));
        
        //Color3f s = 3.5f + 100.0f * Math::pow4(R - 0.33f);
        Color3f d = ld;
        r *= d[channel];
        
        return r / 1000.0f; // Convert to meters
    }

    float sampleSrPdf(float r, int channel) const
    {
        // We use dipole as pdf, although the sampling is for a different Sr function.
        //  They are modeling the same effect so they are pretty close, 
        //  so it is good enough for now
        return dipoleDiffusionAproximation(r)[channel];
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        return evalMS(bRec);
    }


    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return 0.0f;

        return Warp::squareToCosineHemispherePdf(bRec.wo);
    }

    /*
    Color3f samplePointOld(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        Color3f sigmaT = sigmaA + sigmaS;

        // Select random channel
        int channel = sampler->next1D() * 3;
        float sigma = sigmaT[channel];
        float channelPdf = 1.0f / 3.0f;

        // Sample an offset proportional to the sigmaT
        float r = Warp::squareToSrDecay(sampler->next1D(), sigma);
        Point2f sample = Warp::SrToDisk(sampler->next1D(), r);

        // Clamp to account for highly curved surfaces
        sample = Math::max(sample, Point2f(1.0f/sigma));
   
        // Sampled offset (in mm) to world
        Vector3f sampled = Vector3f(sample.x(), sample.y(), 0.0f) / 1000.0f;
        bRec.po = bRec.pi + bRec.fri.ptoWorld(sampled);
        bRec.fro = bRec.fri; 
        //bRec.po = projectToSurface(bRec, bRec.po);

        float pdf = Warp::SrToDiskPdf(sample) 
                * Warp::squareToSrDecayPdf(r, sigma)
                * channelPdf;

        return Color3f(1.0f / pdf);
    }
    */
    

    Color3f samplePoint(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        // Select random frame
        Vector3f vx, vy, vz;
        float rnd = sampler->next1D();
        if (rnd < 0.5f) 
        {
            vx = bRec.fri.s;
            vy = bRec.fri.t;
            vz = bRec.fri.n;
            rnd *= 2; // Adjust random sample
        } 
        else if (rnd < 0.75f) 
        {
            vx = bRec.fri.t;
            vy = bRec.fri.n;
            vz = bRec.fri.s;
            rnd = (rnd - 0.5f) * 4; // Adjust random sample
        } 
        else 
        {
            vx = bRec.fri.n;
            vy = bRec.fri.s;
            vz = bRec.fri.t;
            rnd = (rnd - 0.75f) * 4; // Adjust random sample
        }

        // Select random channel
        int channel =  rnd * 3;
        rnd = rnd * 3.0f - channel; // Adjust random sample

        // Select sphere radius
        float rMax = sampleSr(0.999f, channel);

        float r;
        do {
            r = sampleSr(sampler->next1D(), channel);
        }
        while (r > rMax);

        // Sample a point on the disk
        float l = 2.0f * Math::sqrt(Math::pow2(rMax) - Math::pow2(r));
        float th = 2 * M_PI * sampler->next1D();

        // Project to sphere
        Point3f w_center = bRec.pi + r * (vx*std::cos(th) + vy*std::sin(th));
        Point3f w_base = w_center - l*vz*0.5f;
        Point3f w_target = w_base + l*vz;
 

        /************ Project point to shape ****************/

        Vector3f d = (w_target - w_base).normalized();
        Ray3f ray = Ray3f(w_base, d, 0.0f, l);

        Intersection its;
        bool intersected = bRec.scene->rayIntersect(ray, bRec.mesh, its);

        if (!intersected)
        {
            // This is inefficient, 
            //  but it should not happen frequently 
            return samplePoint(bRec, sampler);
        }

        bRec.fro = its.shFrame;
        bRec.po = ray(its.t);

        /*********************** Pdf *************************/

        float pdf = pointPdf(bRec);

        return Color3f(1.0f / pdf);
    }

    float pointPdf(const BSDFQueryRecord &bRec) const
    {
        Vector3f d = bRec.po - bRec.pi;
        Vector3f l_d = bRec.fri.vtoLocal(d);
        Vector3f l_n = bRec.fri.vtoLocal(bRec.fro.n);

        // Compute the radius that has been sampled for each axis
        //  (d projected into each axis)
        float rProjected[3] = { std::sqrt(l_d.y()*l_d.y() + l_d.z()*l_d.z()),
                                std::sqrt(l_d.z()*l_d.z() + l_d.x()*l_d.x()),
                                std::sqrt(l_d.x()*l_d.x() + l_d.y()*l_d.y()) };

        float pdf = 0, axisPdf[3] = { .25f, .25f, .5f };
        float chProb = 1 / 3.0f; // 3 channels

        for (int axis = 0; axis < 3; axis++)
        {
            for (int ch = 0; ch < 3; ch++)
            {
                pdf += sampleSrPdf(rProjected[axis], ch) 
                        * std::abs(l_n[axis]) 
                        * chProb * axisPdf[axis];
            }
        }
        
        return pdf;
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
        r *= 1000.0f; // Convert to mm
        float r2 = Math::pow2(r);

        Color3f dr = Math::sqrt(r2 + Math::pow2(zr)); 
        Color3f dv = Math::sqrt(r2 + Math::pow2(zv)); 

        Color3f sigmaTrDr = sigmaTr * dr;
        Color3f sigmaTrDv = sigmaTr * dv;

        // Compute main formula
        Color3f C1 = zr * (1 + sigmaTrDr) * Math::exp(-sigmaTrDr) / Math::pow3(dr);
        Color3f C2 = zv * (1 + sigmaTrDv) * Math::exp(-sigmaTrDv) / Math::pow3(dv);

        Color3f result = alpha_4pi * (C1 - C2);

        return result;
    }


    Color3f evalMS(const BSDFQueryRecord &bRec) const
    {
        float r = (bRec.po - bRec.pi).norm();
        float eta = etaT;

        Color3f Rd = dipoleDiffusionAproximation(r);

        Vector3f w_wi = bRec.fro.vtoWorld(bRec.wi);
        float cosWi = Math::absCos(w_wi, bRec.fri.n);
        float cosWo = Math::absCosTheta(bRec.wo);

        float Ft_o = 1 - fresnel(cosWo, 1.0, eta);
        float Ft_i = 1 - fresnel(cosWi, 1.0, eta);

        // Compute the diffusion term
        Color3f Sd = INV_PI * Ft_i * Rd * Ft_o;

        Sd = Math::max(Sd, Color3f(0.0f));

        return Sd * m_albedo->eval(bRec.uv);
    } 

    float Fdr() const 
    {
        const float eta = etaT;
        const float eta2 = eta * eta;
        const float eta3 = eta2 * eta;

        if (eta >= 1.0) 
        {
            return -1.4399 / eta2 + 0.7099 / eta + 0.6681 +
                    0.0636 * eta;
        } 
        else 
        {
            return -0.4399 + 0.7099 / eta - 0.3319 / eta2 +
                    0.0636 / eta3;
        }
    }

    float C1() const 
    {
        const float eta = etaT;
        const float eta2 = eta * eta;
        const float eta3 = eta2 * eta;
        const float eta4 = eta3 * eta;
        const float eta5 = eta4 * eta;

        if (eta < 1.0) 
        {
            return (0.919317 - 3.4793 * eta + 6.75335 * eta2
                    -7.80989 * eta3 + 4.98554 * eta4 - 1.36881 * eta5) / 2.0;
        } 
        else 
        {
            return (-9.23372 + 22.2272 * eta - 20.9292 * eta2 + 10.2291 * eta3
                    - 2.54396 * eta4 + 0.254913 * eta5) / 2.0;
        }
    }

    float C2() const 
    {
        const float eta = etaT;
        const float eta2 = eta * eta;
        const float eta3 = eta2 * eta;
        const float eta4 = eta3 * eta;
        const float eta5 = eta4 * eta;

        if (eta < 1.0) 
        {
            return (0.828421 - 2.62051 * eta + 3.36231 * eta2 - 1.95284 * eta3
                    + 0.236494 * eta4 + 0.145787 * eta5) / 3.0;
        } 
        else {
            return (-1641.1 + 135.926 / eta3 - 656.175 / eta2 + 1376.53 / eta
                    + 1213.67 * eta - 568.556 * eta2 + 164.798 * eta3
                    - 27.0181 * eta4 + 1.91826 * eta5) / 3.0;
        }
    }

    Color3f betterAlternative (float r) const
    {
        r *= 1000.0f; // Convert to mm
        const float r2 = r * r;

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
    std::string toString() const 
    {
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
    Color3f sigmaA, sigmaS, 
        alpha_4pi, sigmaTr, ld, zr, zv; // Precomputed values 
    float etaT, scale;

    std::vector<float> SrCdfInverse, XICdfInverse;
};

NORI_REGISTER_CLASS(subsurface, "subsurface");
NORI_NAMESPACE_END

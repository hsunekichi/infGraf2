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
class subsurfacePbrt : public BSDF
{
public:
    subsurfacePbrt(const PropertyList &propList) 
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
        return beamDiffusionMs(r)[channel];
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        if (!bRec.isCameraRay)
            return INV_PI * m_albedo->eval(bRec.uv);

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
    

    Color3f samplePoint(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        if (!bRec.isCameraRay)
        {
            bRec.fro = bRec.fri;
            bRec.po = bRec.pi;
            return Color3f(1.0f);
        }

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

        std::vector<Intersection> its; its.reserve(4);
        bool intersected = bRec.scene->rayProbe(ray, bRec.mesh, its);

        if (!intersected)
        {
            // This is inefficient, 
            //  but it should not happen frequently 
            return samplePoint(bRec, sampler);
        }

        // Select random intersection
        int iIts = sampler->next1D() * its.size();
        float itsPdf = 1.0f / its.size();
        
        bRec.po = its[iIts].p;
        bRec.fro = its[iIts].shFrame;

        /*********************** Pdf *************************/

        float pdf = pointPdf(bRec) * itsPdf;

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

        Color3f f = eval(bRec) * Frame::cosTheta(bRec.wo) / pdf;
        bRec.isCameraRay = false;

        return f;
    }

    bool isSubsurfaceScattering() const {
        return true;
    }

    bool isDiffuse() const {
        return true;
    }

    float fresnelMoment1(float eta) const
    {
        float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
            eta5 = eta4 * eta;

        if (eta < 1) 
        {
            return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 +
                2.49277f * eta4 - 0.68441f * eta5;
        }
        else
        {
            return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 -
                1.27198f * eta4 + 0.12746f * eta5;
        }
    }

    float fresnelMoment2(float eta) const
    {
        float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
            eta5 = eta4 * eta;

        if (eta < 1) 
        {
            return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 +
                0.07883f * eta4 + 0.04860f * eta5;
        } 
        else 
        {
            float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
            return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 +
                458.843f * r_eta + 404.557f * eta - 189.519f * eta2 +
                54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
        }
    }

    Color3f beamDiffusionMs(float r) const
    {
        float eta = etaT;
        const int nSamples = 100;
        Color3f Ed = 0;

        Color3f _sigmaS = sigmaS;
        Color3f _sigmaT = sigmaA + _sigmaS;
        Color3f rhop = _sigmaS / _sigmaT;

        Color3f D_g = (2 * sigmaA + _sigmaS) / (3 * _sigmaT * _sigmaT);

        Color3f sigmaTR = Math::sqrt(sigmaA / D_g);

        Color3f fm1 = fresnelMoment1(eta), fm2 = fresnelMoment2(eta);
        Color3f ze = -2 * D_g * (1 + 3 * fm2) / (1 - 2 * fm1);

        Color3f cPhi = .25f * (1 - 2 * fm1), cE  = .5f  * (1 - 3 * fm2);

        for (int i = 0; i < nSamples; ++i) 
        {
            Color3f zr = -std::log(1 - (i + .5f) / nSamples) / _sigmaT;

            Color3f zv = -zr + 2 * ze;
            Color3f dr = Math::sqrt(r * r + zr * zr), dv = Math::sqrt(r * r + zv * zv);

            Color3f phiD = INV_FOURPI / D_g *
                (Math::exp(-sigmaTR * dr) / dr - Math::exp(-sigmaTR * dv) / dv);

            Color3f EDn = INV_FOURPI *
                (zr * (1 + sigmaTR * dr) * Math::exp(-sigmaTR * dr) / (dr*dr*dr) -
                zv * (1 + sigmaTR * dv) * Math::exp(-sigmaTR * dv) / (dv*dv*dv));

            Color3f E = phiD * cPhi + EDn * cE;
            Color3f kappa = 1 - Math::exp(-2 * _sigmaT * (dr + zr));
            Ed += kappa * rhop * rhop * E;
        }

        return Ed / nSamples;
    }

    Color3f evalMS(const BSDFQueryRecord &bRec) const
    {
        float r = (bRec.po - bRec.pi).norm();
        float eta = etaT;

        Color3f Rd = beamDiffusionMs(r);

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

NORI_REGISTER_CLASS(subsurfacePbrt, "subsurfacePbrt");
NORI_NAMESPACE_END

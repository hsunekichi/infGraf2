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
 * \brief normalizedDiffusion / Lambertian BRDF model
 */
class NormalizedDiffusion : public BSDF
{
public:
    NormalizedDiffusion(const PropertyList &propList) 
    {
        Color3f kd = propList.getColor("albedo", Color3f(0.5f));

        if (kd.maxCoeff() > 1.0f)
            kd /= 255.0f;

        m_albedo = new ConstantSpectrumTexture(kd);

        ld = propList.getColor("ld", Color3f(0.0f)); // ld is in mm
        
        scale = 1.0f/1000.0f; // Sigmas are in mm^-1

        etaT = propList.getFloat("eta", 1.0f);
        
        // Constant used on the sampling method
        //Color3f D = (sigmaT + sigmaA) / (3.0f * Math::pow2(sigmaT));
        //Color3f sigmaTr = Math::sqrt(sigmaA/D);
        //ld = 1.0f / sigmaTr;

        precomputeSrCdfInverse();
    }

    Color3f scalingFactor(Point2f uv) const
    {
        Color3f kd = m_albedo->eval(uv);

        Color3f scalingFactor = 3.5f + 100*Math::pow4(kd - 0.33f);
        return scalingFactor;
    }

    Color3f Rdisney(float scatterDistance) const
    {
        Color3f a = Color3f(1, 0.674, 0.372);
        //Color3f a2 = a * a;
        //Color3f a3 = a2 * a;

        //Color3f alpha = Color3f(1.0f) - Math::exp(-5.09406f * a + 2.61188f * a2 - 4.31805f * a3);
        Color3f s = Color3f(1.9f) - a + 3.5f * (a - Color3f(0.8f)) * (a - Color3f(0.8f));

        return 1.0f / (s * scatterDistance);
    }

    Color3f Sr_(float r, Point2f uv) const
    {
        r /= scale;

        // Avoid the singularity at r = 0
        r = Math::max(r, 1e-6f);

        //Color3f s = scalingFactor(uv);
        Color3f d = ld;

        Color3f num1 = Math::exp(-r/d);
        Color3f num2 = Math::exp(-r / (3.0f*d));
        Color3f den = 8.0f * M_PI * d * r;
        
        Color3f result = (num1 + num2) / den;

        //std::cout << "Sr: " + result.toString() << std::endl;

        return result;
    }

    Color3f Sr(float r, Point2f uv) const
    {
        Color3f result = Sr_(r, uv);    
    
        return result * m_albedo->eval(uv);
    }

    // Sr CDF^-1
    /*
    float sampleSrAnalitic(float rnd, int channel, Point2f uv, float &pdfInv) const
    {        
        float rcpS = 1 / scalingFactor(uv)[channel];
        rnd = 1 - rnd; // Convert CDF to CCDF; the resulting value of (u != 0)

        float g = 1 + (4 * rnd) * (2 * rnd + sqrt(1 + (4 * rnd) * rnd));
        float n = exp2(log2(g) * (-1.0/3.0));            // g^(-1/3)
        float p = (g * n) * n;                           // g^(+1/3)
        float c = 1 + p + n;                             // 1 + g^(+1/3) + g^(-1/3)
        float x = (3 / LOG2_E) * log2(c / (4 * rnd));    // 3 * Log[c / (4 * rnd)]

        float rcpExp = ((c * c) * c) / ((4 * rnd) * ((c * c) + (4 * rnd) * (4 * rnd)));

        float r = x * rcpS;
        pdfInv = (8 * M_PI * rcpS) * rcpExp; // (8 * Pi) / s / (Exp[-s * r / 3] + Exp[-s * r])

        return r;
    }
    */

    float sampleSr(float rnd, int channel, Point2f uv, float &pdfInv) const
    {
        size_t steps = SrCdfInverse.size() - 1;
        int i_xi = Math::floor(rnd * steps);

        float r1 = SrCdfInverse[i_xi];
        float r2 = SrCdfInverse[i_xi + 1];

        float x1 = XICdfInverse[i_xi];
        float x2 = XICdfInverse[i_xi + 1];

        float r = Math::lerp(r1, r2, (rnd - x1) / (x2 - x1));
        
        //Color3f R = dipoleDiffusionAproximation(r);
        //Color3f s = scalingFactor(uv);
        Color3f d = ld;
        r *= d[channel];
        
        return r * scale; // Convert to meters
    }


    double sampleSrPdf(float r, int channel, Point2f uv) const
    {
        r = Math::max(r, 1e-6f); // Avoid singularity at r == 0.
        return Sr_(r, uv)[channel];
    }


    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle)
            return Color3f(0.0f);

        //if (!bRec.isCameraRay)
        //    return INV_PI * m_albedo->eval(bRec.uv);

        // Compute the Fresnel term
        float Fr_o = 1 - fresnel(Math::absCosTheta(bRec.wo), 1.0f, etaT);
        //float c = 1.0f;

        return Color3f(Fr_o); 
    }

    
    Color3f evalMS(const BSDFQueryRecord &bRec) const
    {
        float Fr_i = 1 - fresnel(Math::absCosTheta(bRec.wi), 1.0f, etaT);
        float r = (bRec.po - bRec.pi).norm();
        
        Color3f Sp = Sr(r, bRec.uv);

        return Fr_i * Sp * INV_PI;
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
        //if (!bRec.isCameraRay)
        //{
        //    bRec.fro = bRec.fri;
        //    bRec.po = bRec.pi;
        //    return Color3f(1.0f);
        //}

        bool validSample = false;
        float itsPdf, srPdfInv;

        while (!validSample)
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
            float rMax = sampleSr(0.9999f, channel, bRec.uv, srPdfInv);

            float r;
            do {
                r = sampleSr(sampler->next1D(), channel, bRec.uv, srPdfInv);
            }
            while (r > rMax);

            r = Math::max(r, 1e-6f);

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

            std::vector<Intersection> its; its.reserve(2);
            bool intersected = bRec.scene->rayProbe(ray, bRec.mesh, its);

            if (!intersected) {
                continue;
            }

            // Select random intersection
            int iIts = sampler->next1D() * its.size();
            itsPdf = 1.0f / its.size();
            
            bRec.po = its[iIts].p;
            bRec.fro = its[iIts].shFrame;

            validSample = true;
        }

        /*********************** Pdf *************************/

        double pdf = pointPdf(bRec) * itsPdf;

        //std::cout << evalMS(bRec).toString() + ", " + std::to_string(pointPdf(bRec)) << std::endl;

        return evalMS(bRec) / pdf;
    }

    double pointPdf(const BSDFQueryRecord &bRec) const
    {
        Vector3f d = bRec.po - bRec.pi;
        Vector3f l_d = bRec.fri.vtoLocal(d);
        Vector3f l_n = bRec.fri.vtoLocal(bRec.fro.n);

        // Compute the radius that has been sampled for each axis
        //  (d projected into each axis)
        double rProjected[3] = { std::sqrt(l_d.y()*l_d.y() + l_d.z()*l_d.z()),
                                std::sqrt(l_d.z()*l_d.z() + l_d.x()*l_d.x()),
                                std::sqrt(l_d.x()*l_d.x() + l_d.y()*l_d.y()) };

        double pdf = 0, axisPdf[3] = { .25f, .25f, .5f };
        double chProb = 1 / 3.0f; // 3 channels

        for (int axis = 0; axis < 3; axis++)
        {
            for (int ch = 0; ch < 3; ch++)
            {
                pdf += sampleSrPdf(rProjected[axis], ch, bRec.uv)
                        * std::abs(l_n[axis]) 
                        * chProb * axisPdf[axis];
            }
        }
        
        return pdf;
    }


    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        bRec.measure = ESolidAngle;
        //float pdf;

        bRec.wo = Warp::squareToCosineHemisphere(sampler->next2D());
        //pdf = Warp::squareToCosineHemispherePdf(bRec.wo);

        Color3f f = eval(bRec) * Frame::cosTheta(bRec.wo);
        bRec.isCameraRay = false;

        return f;
    }

    float SrCdf(float r) const
    {
        return 1.0 - 0.25 * Math::exp(-r) - 0.75 * Math::exp(-r/3);
    }

    void precomputeSrCdfInverse()
    {
        size_t steps = 2024;
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

    bool isSubsurfaceScattering() const {
        return true;
    }

    bool isDiffuse() const {
        return true;
    }

    /// Return a human-readable summary
    std::string toString() const 
    {
        return tfm::format(
            "normalizedDiffusion[\n"
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
                throw NoriException("normalizedDiffusion::addChild(<%s>,%s) is not supported!",
                classTypeName(obj->getClassType()), name);
            break;

        default:
            throw NoriException("normalizedDiffusion::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }


    EClassType getClassType() const { return EBSDF; }
private:
    Texture *m_albedo;
    Color3f sigmaA, sigmaS, sigmaT, ld;
    float etaT, scale;

    std::vector<float> SrCdfInverse, XICdfInverse;
};

NORI_REGISTER_CLASS(NormalizedDiffusion, "normalizedDiffusion");
NORI_NAMESPACE_END

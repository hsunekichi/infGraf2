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

#include <future>

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <nori/dpdf.h>
#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Skin BSSRDF model
 */
class Skin : public BSDF {
public:
    Skin(const PropertyList &propList) 
    {
        
    }

    void preprocess(Sampler *sampler) 
    {
        if (m_bsdfs.size() != 2)
            throw NoriException("Skin::preprocess(): A skin BSDF must have 2 components!");

        for (size_t i = 0; i < m_bsdfs.size(); i++)
            m_bsdfs[i]->preprocess(sampler);

        m_pdf.clear();

        Color3f c1 = Pth::integrateBSDF(m_bsdfs[0], sampler);
        float w1 = c1.getLuminance();
        float ws = 1.0f - w1;

        m_pdf.append(w1);
        m_pdf.append(ws);
        m_pdf.normalize();


        int nSamples = 1024;
        precomputed_specular_integral.resize(nSamples);
        std::vector<std::future<Color3f>> futures(nSamples);

        //auto init = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < nSamples; i++)
        {
            float costheta = float(i) / float(nSamples - 1);
            futures[i] = std::async(std::launch::async, Pth::integrateSkinSpecular, m_bsdfs[0], std::move(sampler->clone()), costheta, specWeight);
        }

        for (int i = 0; i < nSamples; i++) {
            precomputed_specular_integral[i] = futures[i].get();
        }

        //auto end = std::chrono::high_resolution_clock::now();

        //std::chrono::duration<double> elapsed = end - init;
        //std::cout << "Elapsed time: " << elapsed.count() << std::endl;
    }

    Color3f interpolateSpecularWeight(float cosTheta) const
    {
        int index1 = Math::floor(cosTheta * (precomputed_specular_integral.size() - 1));
        int index2 = Math::ceil(cosTheta * (precomputed_specular_integral.size() - 1));

        float cosTheta1 = float(index1) / float(precomputed_specular_integral.size() - 1);
        float cosTheta2 = float(index2) / float(precomputed_specular_integral.size() - 1);

        return Math::lerp(precomputed_specular_integral[index1], precomputed_specular_integral[index2], cosTheta1, cosTheta2);
    }

    Color3f eval(const BSDFQueryRecord &bRec) const 
    {
        if (bRec.agregate_id == 0) {
            return m_bsdfs[bRec.agregate_id]->eval(bRec) * specWeight;
        }
        else
        { 
            Color3f w_wi = interpolateSpecularWeight(bRec.wi.z());
            Color3f w_wo = interpolateSpecularWeight(bRec.wo.z());

            Color3f w = (1.0f - w_wi) * (1.0f - w_wo);
            w = Math::clamp(w, 0.0f, 1.0f);

            return m_bsdfs[bRec.agregate_id]->eval(bRec) * w;
        }
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const 
    {
        return m_pdf[bRec.agregate_id] * m_bsdfs[bRec.agregate_id]->pdf(bRec);
    }


    Color3f samplePoint(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        int index = m_pdf.sample(sampler->next1D());
        bRec.agregate_id = index;

        return m_bsdfs[bRec.agregate_id]->samplePoint(bRec, sampler) / m_pdf[index];
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const 
    {
        return m_bsdfs[bRec.agregate_id]->sample(bRec, sampler);
    }

    bool isSubsurfaceScattering() const 
    {
        return true;
    }

    bool isDiffuse() const {
        return true;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "Diffuse[\n"
            "  albedo = \n");
    }

    void addChild(NoriObject* obj, const std::string& name = "none") {
        switch (obj->getClassType()) {
        case EBSDF:
        {
            BSDF *bsdf = static_cast<BSDF*>(obj);

            if (bsdf->isSubsurfaceScattering() && m_bsdfs.size() != 1)
                throw NoriException("Skin::addChild(): Second skin BSDF must be subsurface scattering!");

            m_bsdfs.push_back(bsdf);

            break;
        }

        default:
            throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }


    EClassType getClassType() const { return EBSDF; }
private:
    //Texture *m_albedo;
    std::vector<BSDF*> m_bsdfs;
    float specWeight = 0.0277778;

    std::vector<Color3f> precomputed_specular_integral;

    DiscretePDF m_pdf;
};

NORI_REGISTER_CLASS(Skin, "skin");
NORI_NAMESPACE_END

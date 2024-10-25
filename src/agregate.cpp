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
#include <nori/dpdf.h>
#include <nori/pathtracing.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Agregate BSSRDF model
 */
class Agregate : public BSDF {
public:
    Agregate(const PropertyList &propList) 
    {
        int i_weight = 0;
        
        while (true)
        {
            i_weight++;

            std::string wname = tfm::format("w%d", i_weight);
            float w = propList.getFloat(wname, -1.0f);

            if (w == -1.0f)
            {
                break;
            }
            else
            {
                m_weights.push_back(w);
            }
        }

        if (m_weights.empty())
            throw NoriException("Agregate::Agregate(): At least one weight must be defined!");

        // Check weights sum to 1
        float sum = 0;
        for (size_t i = 0; i < m_weights.size(); i++)
            sum += m_weights[i];

        if (sum != 1)
            throw NoriException("Agregate::Agregate(): Weights must sum to 1!");


        for (size_t i = 0; i < m_weights.size(); i++)
        {
            m_pdf.append(m_weights[i]);
        }

        m_pdf.normalize();
    }

    void activate()
    {
        is_sss = false;
        for (size_t i = 0; i < m_bsdfs.size(); i++)
            is_sss |= m_bsdfs[i]->isSubsurfaceScattering();

        if (m_weights.size() != m_bsdfs.size())
            throw NoriException("Agregate::activate(): Number of weights and BSDFs must be the same!");
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const 
    {
        return m_bsdfs[bRec.agregate_id]->eval(bRec); // * m_weights[bRec.agregate_id]
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const 
    {
        return m_pdf[bRec.agregate_id] * m_bsdfs[bRec.agregate_id]->pdf(bRec);
    }

    Color3f samplePoint(BSDFQueryRecord &bRec, Sampler *sampler) const
    {
        //Color3f intGGX = Pth::integrateBSDF(m_bsdfs[1], sampler);
        //std::cout << intGGX.toString() << std::endl;

        int index = m_pdf.sample(sampler->next1D());
        bRec.agregate_id = index;

        return m_bsdfs[index]->samplePoint(bRec, sampler);
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const 
    {
        // Agregate index is always defined on sample point
        return m_bsdfs[bRec.agregate_id]->sample(bRec, sampler);
    }

    bool isSubsurfaceScattering() const 
    {
        return is_sss;
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
        case ETexture:
            //if (name == "albedo")
            //{
            //    delete m_albedo;
            //    m_albedo = static_cast<Texture*>(obj);
            //}
            //else
            //    throw NoriException("Diffuse::addChild(<%s>,%s) is not supported!",
            //    classTypeName(obj->getClassType()), name);
            break;

        case EBSDF:
            m_bsdfs.push_back(static_cast<BSDF*>(obj));
            break;

        default:
            throw NoriException("Diffuse::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }


    EClassType getClassType() const { return EBSDF; }
private:
    //Texture *m_albedo;
    std::vector<BSDF*> m_bsdfs;
    std::vector<float> m_weights;
    DiscretePDF m_pdf;

    bool is_sss;
};

NORI_REGISTER_CLASS(Agregate, "agregate");
NORI_NAMESPACE_END

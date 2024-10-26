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
 * \brief Skin BSSRDF model
 */
class Skin : public BSDF {
public:
    Skin(const PropertyList &propList) 
    {
        
    }

    void preprocess(Sampler *sampler) 
    {
        if (m_bsdfs.size() != 3)
            throw NoriException("Skin::preprocess(): A skin BSDF must have 3 components!");

        for (size_t i = 0; i < m_bsdfs.size(); i++)
            m_bsdfs[i]->preprocess(sampler);

        m_pdf.clear();

        BSDF *spec1 = m_bsdfs[1];
        BSDF *spec2 = m_bsdfs[2];
        
        Color3f c1 = Pth::integrateBSDF(spec1, sampler);
        Color3f c2 = Pth::integrateBSDF(spec2, sampler);

        w1 = c1.getLuminance();
        w2 = c2.getLuminance();

        ws = 1.0f - w1 * (1.0f - w2);

        m_pdf.append(ws);
        m_pdf.append(w1);
        m_pdf.append(w2);
        m_pdf.normalize();
    }

    /// Evaluate the BRDF model
    Color3f evalIntegrated(const BSDFQueryRecord &bRec) const 
    {
        float w = 1.0f;

        if (bRec.agregate_id == 0)
            w = ws; 

        return m_bsdfs[bRec.agregate_id]->eval(bRec) * w;
    }

    Color3f eval(const BSDFQueryRecord &bRec) const 
    {
        if (bRec.agregate_id == 1) {
            return m_bsdfs[bRec.agregate_id]->eval(bRec) * w1;
        }
        else if (bRec.agregate_id == 2) {
            return m_bsdfs[bRec.agregate_id]->eval(bRec) * w2;
        }
        else
        {
            Color3f _w1 = m_bsdfs[1]->eval(bRec).getLuminance();
            Color3f _w2 = m_bsdfs[2]->eval(bRec).getLuminance();

            Color3f w = 1.0f - _w1 * (1.0f - _w2);
            w = Math::clamp(w, 0.0f, 1.0f);

            return m_bsdfs[0]->eval(bRec) * w;
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

        return m_bsdfs[bRec.agregate_id]->samplePoint(bRec, sampler);
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

            if (bsdf->isSubsurfaceScattering() && m_bsdfs.size() > 0)
                throw NoriException("Skin::addChild(): First skin BSDF must be subsurface scattering!");

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
    float w1, w2, ws;
    
    DiscretePDF m_pdf;
};

NORI_REGISTER_CLASS(Skin, "skin");
NORI_NAMESPACE_END

// Medium.h
#pragma once
#include <nori/phase.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Representation of a participating medium (volumetric scattering)
 */
class Medium {
public:
    Medium(float scatteringCoeff, float absorptionCoeff, PhaseFunction *phaseFunction)
        : m_scatteringCoeff(scatteringCoeff), m_absorptionCoeff(absorptionCoeff), m_phaseFunction(phaseFunction) { }

    /// Get the scattering coefficient
    float getScatteringCoeff() const { return m_scatteringCoeff; }

    /// Get the absorption coefficient
    float getAbsorptionCoeff() const { return m_absorptionCoeff; }

    /// Get the phase function for scattering in the medium
    PhaseFunction *getPhaseFunction() const { return m_phaseFunction; }

    /**
     * \brief Sample a scattering event within the medium
     * Returns true if a scattering event occurs, false otherwise
     */
    bool sampleScattering(float distance, float& scatteringDistance) const {
        // Example: Exponentially distribute scattering based on distance
        float probScatter = 1 - exp(-m_scatteringCoeff * distance);
        if (random() < probScatter) {
            scatteringDistance = -log(1 - random()) / m_scatteringCoeff;
            return true;
        }
        return false;
    }

private:
    float m_scatteringCoeff;   ///< Scattering coefficient
    float m_absorptionCoeff;   ///< Absorption coefficient
    PhaseFunction *m_phaseFunction;  ///< Phase function for light scattering

};

NORI_NAMESPACE_END

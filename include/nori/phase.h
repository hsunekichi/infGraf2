// PhaseFunction.h
#pragma once
#include <vector>
#include <nori/math.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Henyey-Greenstein phase function
 */
class PhaseFunction {
public:
    PhaseFunction(float g = 0.0f) : m_g(g) {}

    /**
     * \brief Compute the phase function value
     * \param wo The outgoing direction
     * \param wi The incoming direction
     * \return The phase function value
     */
    float eval(const Vector3f &wo, const Vector3f &wi) const {
        float cosTheta = wo.dot(wi);
        float denom = 1 + m_g * (m_g + 2 * cosTheta);
        return (1 - m_g * m_g) / (4 * M_PI * denom * std::sqrt(denom));
    }

    /**
     * \brief Sample the phase function to obtain a scattering direction
     */
    Vector3f sample(const Vector3f &wo) const {
        // Sampling based on Henyey-Greenstein distribution
        float cosTheta = (1 + m_g * m_g - pow((1 - random() + (1 - m_g * m_g)), 2)) / (2 * m_g);
        float sinTheta = std::sqrt(1 - cosTheta * cosTheta);
        float phi = 2 * M_PI * random();
        return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    }

private:
    float m_g;  ///< Anisotropy parameter (-1 <= g <= 1)
};

NORI_NAMESPACE_END

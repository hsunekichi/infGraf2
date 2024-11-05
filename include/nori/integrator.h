/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

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

#pragma once

#include <nori/object.h>
#include <nori/scene.h>
#include <nori/PathState.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Abstract integrator (i.e. a rendering technique)
 *
 * In Nori, the different rendering techniques are collectively referred to as 
 * integrators, since they perform integration over a high-dimensional
 * space. Each integrator represents a specific approach for solving
 * the light transport equation---usually favored in certain scenarios, but
 * at the same time affected by its own set of intrinsic limitations.
 */
class Integrator : public NoriObject {
public:
    /// Release all memory
    virtual ~Integrator() { }

    /// Perform an (optional) preprocess step
    virtual void preprocess(const Scene *scene, Sampler *sampler) { }

    /**
     * \brief Sample the incident radiance along a ray
     *
     * \param scene
     *    A pointer to the underlying scene
     * \param sampler
     *    A pointer to a sample generator
     * \param ray
     *    The ray in question
     * \return
     *    A (usually) unbiased estimate of the radiance in this direction
     */
    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const 
    {
        throw NoriException("Li() is not implemented!");
    }

    // Apply the scattering of the intersection and generate new direction
    //  This is the main integrator function
    virtual void shadeIntersection(const Scene *scene, Sampler *sampler, PathState &state) const = 0;
    
    // Shade the ray when intersecting with the environment (ray intersected on the infinite)
    //  This is zero by default
    virtual void shadeEnvironment(const Scene *scene, Sampler *sampler, PathState &state) const {}

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EIntegrator; }
};

NORI_NAMESPACE_END

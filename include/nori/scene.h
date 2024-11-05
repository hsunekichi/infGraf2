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

#pragma once

#include <nori/accel.cu>
#include <nori/common.h>
#include <nori/PathState.h>



NORI_NAMESPACE_BEGIN


struct Photon
{
    Point3f p;
    Vector3f d; Normal3f n;
    float pdf;
    
    Mesh *mesh = nullptr;
    Color3f radiance = Color3f(0.0f);

    Photon (const Point3f &p, float pdf, Mesh *mesh) : p(p), pdf(pdf), mesh(mesh) {}
};


/**
 * \brief Main scene data structure
 *
 * This class holds information on scene objects and is responsible for
 * coordinating rendering jobs. It also provides useful query routines that
 * are mostly used by the \ref Integrator implementations.
 */
class Scene : public NoriObject {
public:
    /// Construct a new scene object
    Scene(const PropertyList &);

    /// Release all memory
    virtual ~Scene();

    /// Return a pointer to the scene's kd-tree
    const Accel *getAccel() const { return m_accel; }

    /// Return a pointer to the scene's integrator
    const Integrator *getIntegrator() const { return m_integrator; }

    /// Return a pointer to the scene's integrator
    Integrator *getIntegrator() { return m_integrator; }

    /// Return a pointer to the scene's camera
    const Camera *getCamera() const { return m_camera; }

    /// Return a pointer to the scene's sample generator (const version)
    const Sampler *getSampler() const { return m_sampler; }

    /// Return a pointer to the scene's sample generator
    Sampler *getSampler() { return m_sampler; }

    /// Return a reference to an array containing all meshes
    const std::vector<Mesh *> &getMeshes() const { return m_meshes; }

	/// Return a reference to an array containing all lights
	const std::vector<Emitter *> &getLights() const { return m_emitters; }
    const std::vector<Mesh *> &getSSMeshes() const { return sss_meshes; }

	/// Return a the scene background
	Color3f getBackground(const Ray3f& ray) const;

	/// Sample emitter
	Emitter *sampleEmitter(Sampler* sampler, float &pdf) const;

    float pdfEmitter(const Emitter *em) const;

    void preprocess();

	/// Get enviromental emmiter
	const Emitter *getEnvironmentalEmitter() const
	{
		return m_enviromentalEmitter;
	}

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool isShadowRay=false) const {
        return m_accel->rayIntersect(ray, its, isShadowRay);
    }

    static uint64_t float_to_morton(float value, int max_bits = 10) 
    {
        uint64_t morton_code = 0;
        uint32_t value_scaled = static_cast<uint32_t>(value * ((1 << max_bits) - 1)); // Scale to the number of bits
        
        for (int i = 0; i < max_bits; ++i) {
            morton_code |= ((value_scaled >> i) & 1) << (3 * i); // 3D interleaving for Morton code
        }
        
        return morton_code;
    }


    // Function to compute the Morton code for a 3D point in [0, 1]
    static uint64_t morton3D(Point3f p, int bits = 10)
    {
        uint64_t x_morton = float_to_morton(p.x(), bits);
        uint64_t y_morton = float_to_morton(p.y(), bits);
        uint64_t z_morton = float_to_morton(p.z(), bits);  

        uint64_t morton_code = 0;
        for (int i = 0; i < bits; ++i) {
            morton_code |= ((x_morton >> (3 * i)) & 1) << (3 * i);
            morton_code |= ((y_morton >> (3 * i)) & 1) << (3 * i + 1);
            morton_code |= ((z_morton >> (3 * i)) & 1) << (3 * i + 2);
        }

        return morton_code;
    }

    void order_rays_by_morton(
            const std::vector<PathState> &states,
            std::vector<Ray3f> &rays,
            std::vector<uint32_t> &indices) const
    {
        BoundingBox3f m_bbox;
        for (const auto &state : states)
            m_bbox.expandBy(state.ray.o);

        Vector3f bbox_size = m_bbox.max - m_bbox.min;
        std::vector<uint64_t> morton_codes(states.size());

        for (uint32_t i = 0; i < states.size(); ++i)
        {
            Vector3f p = (states[i].ray.o - m_bbox.min);
            
            p = p.cwiseQuotient(bbox_size);
            morton_codes[i] = morton3D(p);
        }

        for (uint32_t i = 0; i < states.size(); ++i)
            indices[i] = i;

        std::sort(indices.begin(), indices.end(), [&morton_codes](uint32_t a, uint32_t b) {
            return morton_codes[a] < morton_codes[b];
        });

        for (uint32_t i = 0; i < states.size(); ++i)
            rays[i] = states[indices[i]].ray;
    }

    void rayIntersect(const std::vector<bool> &aliveMask, 
                        std::vector<PathState> &states,
                        std::vector<Ray3f> &rays,
                        std::vector<Intersection> its,
                        std::vector<bool> &b_its,
                        std::vector<uint32_t> &indices) const
    {
        order_rays_by_morton(states, rays, indices);
        
        m_accel->rayIntersect(aliveMask, rays, its, b_its);

        for (size_t i = 0; i < states.size(); i++)
        {
            states[indices[i]].rayIntersected = b_its[i];

            if (b_its[i])
                states[indices[i]].intersection = its[i];
        }
    }

    bool rayProbe(const Ray3f &ray, const Mesh *mesh, std::vector<Intersection> &its) const 
    {
        int id = -1;

        for (int i = 0; i < (int)sss_meshes.size(); i++)
        {
            if (sss_meshes[i] == mesh)
            {
                id = i;
                break;
            }
        }

        if (id == -1)
            throw NoriException("Requesting ray intersection for a non sss mesh");
    
        return sss_accelerators[id]->rayProbe(ray, its);
    }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and \a only determine whether or not there is an intersection.
     *
     * This method much faster than the other ray tracing function,
     * but the performance comes at the cost of not providing any
     * additional information about the detected intersection
     * (not even its position).
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray) const {
        Intersection its; /* Unused */
        return m_accel->rayIntersect(ray, its, true);
    }

    /// \brief Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const {
        return m_accel->getBoundingBox();
    }

    /**
     * \brief Inherited from \ref NoriObject::activate()
     *
     * Initializes the internal data structures (kd-tree,
     * emitter sampling data structures, etc.)
     */
    void activate();

    /// Add a child object to the scene (meshes, integrators etc.)
    void addChild(NoriObject *obj, const std::string& name = "none");

    /// Return a string summary of the scene (for debugging purposes)
    std::string toString() const;

    EClassType getClassType() const { return EScene; }
private:
    std::vector<Mesh *> m_meshes;
    std::vector<Mesh *> sss_meshes;
	std::vector<Emitter *> m_emitters;

	Emitter *m_enviromentalEmitter = nullptr;
	
    DiscretePDF m_emitterPDF;

    Integrator *m_integrator = nullptr;
    Sampler *m_sampler = nullptr;
    Camera *m_camera = nullptr;
    Accel *m_accel = nullptr;
    std::vector<Accel *> sss_accelerators;
};

NORI_NAMESPACE_END

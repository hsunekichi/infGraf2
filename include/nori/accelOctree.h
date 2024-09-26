/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/mesh.h>
#include <nori/bbox.h>


NORI_NAMESPACE_BEGIN


struct OctreeNode {
    BoundingBox3f bbox;          // Bounding box of the node
    std::vector<uint32_t> triangles; // List of triangle indices (for leaf nodes)
    OctreeNode* children[8];     // Children of the node (for internal nodes)

    OctreeNode() {
        // Initialize children to nullptr
        std::fill(std::begin(children), std::end(children), nullptr);
    }

    bool isLeaf() const {
        // A node is a leaf if it has no children
        return children[0] == nullptr && children[1] == nullptr && children[2] == nullptr && children[3] == nullptr && children[4] == nullptr && children[5] == nullptr && children[6] == nullptr && children[7] == nullptr;
    }
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;


private:
    std::vector<Mesh*> v_mesh; ///< Mesh
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene

    std::vector<OctreeNode*> v_octreeRoot; // Pointer to the root of the Octree

    bool traverseOctree(Ray3f &ray, Intersection &its, const OctreeNode *node, const int idx, bool shadowRay) const;
};

NORI_NAMESPACE_END

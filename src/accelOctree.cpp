/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/accelOctree.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN


OctreeNode* buildOctree(const Mesh *m_mesh, const BoundingBox3f &bbox, const std::vector<uint32_t> &triangles, int depth) {

    const int MAX_DEPTH = 10; // Max octree depth (10-15)

    if (triangles.size() == 0){
        return nullptr;
    }

    // Create a new internal node
    OctreeNode *node = new OctreeNode();
    node->bbox = bbox;

    // Stop subdividing if there are fewer than 10 triangles or max depth is reached
    if (triangles.size() <= 10 || depth >= MAX_DEPTH) {
        node->triangles = triangles;
        return node;
    }

    // Subdivide the current bounding box into 8 smaller boxes
    BoundingBox3f childrenBBoxes[8];
    
     // Get the center point of the bounding box
    Point3f center = bbox.getCenter();

    // Define the 8 sub-boxes using the corners and the center
    childrenBBoxes[0] = BoundingBox3f(bbox.getCorner(0), center); // Bottom-left (min.x, min.y, min.z)
    childrenBBoxes[1] = BoundingBox3f({ center.x(), bbox.getCorner(1).y(), bbox.getCorner(1).z() }, { bbox.getCorner(1).x(), center.y(), center.z() }); // Bottom-right
    childrenBBoxes[2] = BoundingBox3f({ bbox.getCorner(2).x(), center.y(), bbox.getCorner(2).z() }, { center.x(), bbox.getCorner(2).y(), center.z() }); // Top-left
    childrenBBoxes[3] = BoundingBox3f({ center.x(), center.y(), bbox.getCorner(3).z() }, { bbox.getCorner(3).x(), bbox.getCorner(3).y(), center.z() }); // Top-right
    childrenBBoxes[4] = BoundingBox3f({ bbox.getCorner(4).x(), bbox.getCorner(4).y(), center.z() }, { center.x(), center.y(), bbox.getCorner(4).z() }); // Bottom-left (above)
    childrenBBoxes[5] = BoundingBox3f({ center.x(), bbox.getCorner(5).y(), center.z() }, { bbox.getCorner(5).x(), center.y(), bbox.getCorner(5).z() }); // Bottom-right (above)
    childrenBBoxes[6] = BoundingBox3f({ bbox.getCorner(6).x(), center.y(), center.z() }, { center.x(), bbox.getCorner(6).y(), bbox.getCorner(6).z() }); // Top-left (above)
    childrenBBoxes[7] = BoundingBox3f(center, bbox.getCorner(7)); // Top-right (above)

    // List of triangles for each child node
    std::vector<uint32_t> childTriangles[8];

    // Classify triangles into child nodes
    for (uint32_t triIdx : triangles) {
        BoundingBox3f triBBox = m_mesh->getBoundingBox(triIdx);
        for (int i = 0; i < 8; ++i) {
            if (childrenBBoxes[i].overlaps(triBBox)) {
                childTriangles[i].push_back(triIdx);
            }
        }
    }

    // Recursively build child nodes
    for (int i = 0; i < 8; ++i) {
        if (!childTriangles[i].empty()) {
            node->children[i] = buildOctree(m_mesh, childrenBBoxes[i], childTriangles[i], depth + 1);
        }
    }

    return node;
}

void Accel::addMesh(Mesh *mesh) 
{
    if (v_mesh.size() == 0)
        m_bbox = mesh->getBoundingBox();
    else
        m_bbox.expandBy(mesh->getBoundingBox());

    v_mesh.push_back(mesh);
}


void Accel::build() {
    if (v_mesh.empty()) {
        throw NoriException("No mesh available for building the acceleration structure!");
    }

    // Obtain all the idx of the triangles
    std::vector<uint32_t> triangles;

    // Iterate over the meshes & build the respectives sub-Octrees
    for (const auto& m_mesh : v_mesh) {
        BoundingBox3f meshBBox = m_mesh->getBoundingBox();

        // Obtain all the triangles of the current mesh
        triangles.clear();
        for (uint32_t i = 0; i < m_mesh->getTriangleCount(); ++i) {
            triangles.push_back(i);
        }

        // Build the sub-Octree for the current mesh
        v_octreeRoot.push_back(buildOctree(m_mesh, meshBBox, triangles, 0));
    }

    std::cout << "Octree built successfully with bounding box: " << m_bbox.toString() << std::endl;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    //uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    // Call the recursive function that check the Octree
    for (size_t i = 0; i < v_octreeRoot.size(); i++)
    {
        foundIntersection |= traverseOctree(ray, its, v_octreeRoot[i], i, shadowRay);
    }

    return foundIntersection;
}

// TODO: hacer que solo lance el rayo de interseccion con el indice del triangulo y el mesh
bool Accel::traverseOctree(Ray3f &ray, Intersection &its, const OctreeNode *node, const int i_mesh, bool shadowRay) const {
    // If the node is null or the ray doesn't intersect with the bounding box, exit
    if (!node || !node->bbox.rayIntersect(ray)) {
        return false;
    }

    // If the node is a leaf, check the triangles
    if (node->isLeaf()) {
        bool foundIntersection = false;
        uint32_t f = (uint32_t) -1;

        for (uint32_t idx : node->triangles) {
            float u, v, t;
            if (v_mesh[i_mesh]->rayIntersect(idx, ray, u, v, t)) 
            {
                // If there is an intersection and it is a shadow ray, finish
                if (shadowRay)
                    return true;
                
                // If the intersection is closer, update the info
                if (t < ray.maxt) {
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = v_mesh[i_mesh];
                    f = idx;
                    foundIntersection = true;
                }
            }
        }
        
        if (foundIntersection) {
            /* At this point, we now know that there is an intersection,
            and we know the triangle index of the closest such intersection.

            The following computes a number of additional properties which
            characterize the intersection (normals, texture coordinates, etc..)
            */

            /* Find the barycentric coordinates */
            Vector3f bary;
            bary << 1-its.uv.sum(), its.uv;

            /* References to all relevant mesh buffers */
            const Mesh *mesh   = its.mesh;
            const MatrixXf &V  = mesh->getVertexPositions();
            const MatrixXf &N  = mesh->getVertexNormals();
            const MatrixXf &UV = mesh->getVertexTexCoords();
            const MatrixXu &F  = mesh->getIndices();

            /* Vertex indices of the triangle */
            uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

            Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

            /* Compute the intersection positon accurately
            using barycentric coordinates */
            its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

            /* Compute proper texture coordinates if provided by the mesh */
            if (UV.size() > 0)
                its.uv = bary.x() * UV.col(idx0) +
                    bary.y() * UV.col(idx1) +
                    bary.z() * UV.col(idx2);

            /* Compute the geometry frame */
            its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

            if (N.size() > 0) {
                /* Compute the shading frame. Note that for simplicity,
                the current implementation doesn't attempt to provide
                tangents that are continuous across the surface. That
                means that this code will need to be modified to be able
                use anisotropic BRDFs, which need tangent continuity */

                its.shFrame = Frame(
                    (bary.x() * N.col(idx0) +
                    bary.y() * N.col(idx1) +
                    bary.z() * N.col(idx2)).normalized());
            } else {
                its.shFrame = its.geoFrame;
            }
        }

        return foundIntersection;
    }

    // If the node isn't a leaf, iterate to the children
    bool found = false;
    for (int i = 0; i < 8; ++i) {
        if (node->children[i]) {
            found |= traverseOctree(ray, its, node->children[i], i_mesh, shadowRay);
            if (shadowRay && found) {
                return true;
            }
        }
    }

    return found;
}


NORI_NAMESPACE_END


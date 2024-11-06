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

#include <nori/mesh.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


NORI_NAMESPACE_BEGIN



struct BVHmesh
{
    std::vector<Mesh *> originalMeshes;
	size_t nMeshes = 0;
	MatrixXf      m_V;                   ///< Vertex positions
	MatrixXu      m_F;                   ///< Faces
	BoundingBox3f m_bbox;                ///< Bounding box

	const Mesh *getOriginalMesh(n_UINT index) const { return originalMeshes[index]; }

	void uploadToDevice(MatrixXf *&d_V, MatrixXu *&d_F) const
	{
		cudaMalloc(&d_V, m_V.size() * sizeof(float));
		cudaMalloc(&d_F, m_F.size() * sizeof(n_UINT));

		cudaMemcpy(d_V, m_V.data(), m_V.size() * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_F, m_F.data(), m_F.size() * sizeof(n_UINT), cudaMemcpyHostToDevice);
	}

	void clear() {
		m_V.resize(3, 0);
		m_F.resize(3, 0);
		originalMeshes.clear();
		m_bbox = BoundingBox3f();
	}

	__host__ __device__
	bool rayIntersect(n_UINT index, 
			const Ray3f &ray, float &u, float &v, float &t) const
	{
		n_UINT i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    	const Eigen::Vector3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

		/* Find vectors for two edges sharing v[0] */
		Eigen::Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

		/* Begin calculating determinant - also used to calculate U parameter */
		Eigen::Vector3f pvec = ray.d.cross(edge2);

		/* If determinant is near zero, ray lies in plane of triangle */
		float det = edge1.dot(pvec);

		if (det > -1e-8f && det < 1e-8f)
			return false;
		float inv_det = 1.0f / det;

		/* Calculate distance from v[0] to ray origin */
		Eigen::Vector3f tvec = ray.o - p0;

		/* Calculate U parameter and test bounds */
		u = tvec.dot(pvec) * inv_det;
		if (u < 0.0 || u > 1.0)
			return false;

		/* Prepare to test V parameter */
		Eigen::Vector3f qvec = tvec.cross(edge1);

		/* Calculate V parameter and test bounds */
		v = ray.d.dot(qvec) * inv_det;
		if (v < 0.0 || u + v > 1.0)
			return false;

		/* Ray intersects triangle -> compute t */
		t = edge2.dot(qvec) * inv_det;

		return t >= ray.mint && t <= ray.maxt;
	}

	BoundingBox3f getBoundingBox(n_UINT index) const {
		BoundingBox3f result(m_V.col(m_F(0, index)));
		result.expandBy(m_V.col(m_F(1, index)));
		result.expandBy(m_V.col(m_F(2, index)));
		return result;
	}

	Point3f getCentroid(n_UINT index) const {
		return (1.0f / 3.0f) *
			(m_V.col(m_F(0, index)) +
			m_V.col(m_F(1, index)) +
			m_V.col(m_F(2, index)));
	}

	size_t size() const { return m_F.cols(); }


	void append(Mesh *other, size_t meshId)
	{
		if (m_V.cols() == 0)
		{
			m_V = other->getVertexPositions();
			m_F = other->getIndices();
		}
		else
		{
			m_V.conservativeResize(3, m_V.cols() + other->getVertexPositions().cols());
			m_F.conservativeResize(3, m_F.cols() + other->getIndices().cols());

			m_V.rightCols(other->getVertexPositions().cols()) = other->getVertexPositions();
			m_F.rightCols(other->getIndices().cols()) = other->getIndices();
		}

		// Add the original mesh
		for (size_t i = 0; i < other->getIndices().cols(); ++i)
			originalMeshes.push_back(other);

		m_bbox.expandBy(other->getBoundingBox());
	}
};  

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
	friend class BVHBuildTask;
public:
	/// Create a new and empty BVH
	Accel() { m_meshOffset.push_back(0u); }

	/// Release all resources
	void clear();

	/// Release all resources
	virtual ~Accel() { clear(); };

	/**
	 * \brief Register a triangle mesh for inclusion in the BVH.
	 *
	 * This function can only be used before \ref build() is called
	 */
	void addMesh(Mesh *mesh);

	/// Build the BVH
	void build();

	/**
	 * \brief Intersect a ray against all triangle meshes registered
	 * with the BVH
	 *
	 * Detailed information about the intersection, if any, will be
	 * stored in the provided \ref Intersection data record.
	 *
	 * The <tt>shadowRay</tt> parameter specifies whether this detailed
	 * information is really needed. When set to \c true, the
	 * function just checks whether or not there is occlusion, but without
	 * providing any more detail (i.e. \c its will not be filled with
	 * contents). This is usually much faster.
	 *
	 * \return \c true If an intersection was found
	 */
	__host__ __device__
	bool rayIntersect(const Ray3f &ray, Intersection &its,
		bool shadowRay = false) const;

	void rayIntersect(const std::vector<bool> &mask, 
						const std::vector<Ray3f> &ray, 
						std::vector<Intersection> &its,
						std::vector<bool> &b_its) const;

	bool rayProbe(const Ray3f &_ray, 
		std::vector<Intersection> &its) const;

	/// Return the total number of meshes registered with the BVH
	n_UINT getMeshCount() const { return (n_UINT)m_meshes.size(); }

	/// Return the total number of internally represented triangles 
	n_UINT getTriangleCount() const { return globalMesh.size(); }

	/// Return one of the registered meshes
	Mesh *getMesh(n_UINT idx) { return m_meshes[idx]; }

	/// Return one of the registered meshes (const version)
	const Mesh *getMesh(n_UINT idx) const { return m_meshes[idx]; }

	//// Return an axis-aligned bounding box containing the entire tree
	const BoundingBox3f &getBoundingBox() const {
		return m_bbox;
	}

public:
	/**
	 * \brief Compute the mesh and triangle indices corresponding to
	 * a primitive index used by the underlying generic BVH implementation.
	 */
	n_UINT findMesh(n_UINT &idx) const {
		auto it = std::lower_bound(m_meshOffset.begin(), m_meshOffset.end(), idx + 1) - 1;
		idx -= *it;
		return (n_UINT)(it - m_meshOffset.begin());
	}

	//// Return an axis-aligned bounding box containing the given triangle
	BoundingBox3f getBoundingBox(n_UINT index) const {
		return globalMesh.getBoundingBox(index);
	}

	//// Return the centroid of the given triangle
	Point3f getCentroid(n_UINT index) const {
		return globalMesh.getCentroid(index);
	}

	void compactMeshes();

	/// Compute internal tree statistics
	std::pair<float, n_UINT> statistics(n_UINT index = 0) const;

	Intersection fillIntersection(Point2f uv, n_UINT f, float t) const;

	/* BVH node in 32 bytes */
	struct BVHNode {
		union {
			struct {
				unsigned flag : 1;
				uint32_t size : 31;
				n_UINT start;
			} leaf;

			struct {
				unsigned flag : 1;
				uint32_t axis : 31;
				n_UINT rightChild;
			} inner;

			uint64_t data;
		};
		BoundingBox3f bbox;

		__host__ __device__
		bool isLeaf() const {
			return leaf.flag == 1;
		}

		__host__ __device__
		bool isInner() const {
			return leaf.flag == 0;
		}

		__host__ __device__
		bool isUnused() const {
			return data == 0;
		}

		__host__ __device__
		n_UINT start() const {
			return leaf.start;
		}

		__host__ __device__
		n_UINT end() const {
			return leaf.start + leaf.size;
		}
	};
private:
	std::vector<Mesh *> m_meshes;       ///< List of meshes registered with the BVH
	std::vector<n_UINT> m_meshOffset; ///< Index of the first triangle for each shape
	thrust::host_vector<n_UINT> m_indices;    ///< Index references by BVH nodes
	thrust::device_vector<n_UINT> device_indices;    ///< Index references by BVH nodes


	thrust::host_vector<BVHNode> m_nodes;       ///< BVH nodes
	thrust::device_vector<BVHNode> device_nodes; ///< BVH nodes

	BVHmesh globalMesh;                 ///< Global mesh containing all triangles
	MatrixXf *deviceVertices;            ///< Global vertex positions
	MatrixXu *deviceFaces;             ///< Global triangle indices

	BoundingBox3f m_bbox;               ///< Bounding box of the entire BVH
};


NORI_NAMESPACE_END

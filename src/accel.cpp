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

#include <nori/accel.h>
#include <nori/timer.h>
#include <tbb/tbb.h>
#include <Eigen/Geometry>
#include <atomic>
#include <ranges>
#include <algorithm>
#include <execution>


NORI_NAMESPACE_BEGIN

/* Bin data structure for counting triangles and computing their bounding box */
struct Bins {
	static const int BIN_COUNT = 16;
	Bins() { memset(counts, 0, sizeof(n_UINT) * BIN_COUNT); }
	n_UINT counts[BIN_COUNT];
	BoundingBox3f bbox[BIN_COUNT];
};

/**
 * \brief Build task for parallel BVH construction
 *
 * This class uses the task scheduling system of Intel' Thread Building Blocks
 * to parallelize the divide and conquer BVH build at all levels.
 *
 * The used methodology is roughly that described in
 * "Fast and Parallel Construction of SAH-based Bounding Volume Hierarchies"
 * by Ingo Wald (Proc. IEEE/EG Symposium on Interactive Ray Tracing, 2007)
 */
class BVHBuildTask : public tbb::task {
private:
	Accel &bvh;
	n_UINT node_idx;
	n_UINT *start, *end, *temp;

public:
	/// Build-related parameters
	enum {
		/// Switch to a serial build when less than 32 triangles are left
		SERIAL_THRESHOLD = 32,

		/// Process triangles in batches of 1K for the purpose of parallelization
		GRAIN_SIZE = 1000,

		/// Heuristic cost value for traversal operations
		TRAVERSAL_COST = 1,

		/// Heuristic cost value for intersection operations
		INTERSECTION_COST = 1
	};

public:
	/**
	 * Create a new build task
	 *
	 * \param bvh
	 *    Reference to the underlying BVH
	 *
	 * \param node_idx
	 *    Index of the BVH node that should be built
	 *
	 * \param start
	 *    Start pointer into a list of triangle indices to be processed
	 *
	 * \param end
	 *    End pointer into a list of triangle indices to be processed
	 *
	 *  \param temp
	 *    Pointer into a temporary memory region that can be used for
	 *    construction purposes. The usable length is <tt>end-start</tt>
	 *    unsigned integers.
	 */
	BVHBuildTask(Accel &bvh, n_UINT node_idx, n_UINT *start, n_UINT *end, n_UINT *temp)
		: bvh(bvh), node_idx(node_idx), start(start), end(end), temp(temp) { }

	task *execute() {
		n_UINT size = (n_UINT)(end - start);
		Accel::BVHNode &node = bvh.m_nodes[node_idx];

		/* Switch to a serial build when less than SERIAL_THRESHOLD triangles are left */
		if (size < SERIAL_THRESHOLD) {
			execute_serially(bvh, node_idx, start, end, temp);
			return nullptr;
		}

		/* Always split along the largest axis */
		int axis = node.bbox.getLargestAxis();
		float min = node.bbox.min[axis], max = node.bbox.max[axis],
			inv_bin_size = Bins::BIN_COUNT / (max - min);

		/* Accumulate all triangles into bins */
		Bins bins = tbb::parallel_reduce(
			tbb::blocked_range<n_UINT>(0u, size, GRAIN_SIZE),
			Bins(),
			/* MAP: Bin a number of triangles and return the resulting 'Bins' data structure */
			[&](const tbb::blocked_range<n_UINT> &range, Bins result) {
			for (n_UINT i = range.begin(); i != range.end(); ++i) {
				n_UINT f = start[i];
				float centroid = bvh.getCentroid(f)[axis];

				int index = std::min(std::max(
					(int)((centroid - min) * inv_bin_size), 0),
					(Bins::BIN_COUNT - 1));

				result.counts[index]++;
				result.bbox[index].expandBy(bvh.getBoundingBox(f));
			}
			return result;
		},
			/* REDUCE: Combine two 'Bins' data structures */
			[](const Bins &b1, const Bins &b2) {
			Bins result;
			for (int i = 0; i < Bins::BIN_COUNT; ++i) {
				result.counts[i] = b1.counts[i] + b2.counts[i];
				result.bbox[i] = BoundingBox3f::merge(b1.bbox[i], b2.bbox[i]);
			}
			return result;
		}
		);

		/* Choose the best split plane based on the binned data */
		BoundingBox3f bbox_left[Bins::BIN_COUNT];
		bbox_left[0] = bins.bbox[0];
		for (int i = 1; i < Bins::BIN_COUNT; ++i) {
			bins.counts[i] += bins.counts[i - 1];
			bbox_left[i] = BoundingBox3f::merge(bbox_left[i - 1], bins.bbox[i]);
		}

		BoundingBox3f bbox_right = bins.bbox[Bins::BIN_COUNT - 1], best_bbox_right;
		int64_t best_index = -1;
		float best_cost = (float)INTERSECTION_COST * size;
		float tri_factor = (float)INTERSECTION_COST / node.bbox.getSurfaceArea();

		for (int i = Bins::BIN_COUNT - 2; i >= 0; --i) {
			n_UINT prims_left = bins.counts[i], prims_right = (n_UINT)(end - start) - bins.counts[i];
			float sah_cost = 2.0f * TRAVERSAL_COST +
				tri_factor * (prims_left * bbox_left[i].getSurfaceArea() +
					prims_right * bbox_right.getSurfaceArea());
			if (sah_cost < best_cost) {
				best_cost = sah_cost;
				best_index = i;
				best_bbox_right = bbox_right;
			}
			bbox_right = BoundingBox3f::merge(bbox_right, bins.bbox[i]);
		}

		if (best_index == -1) {
			/* Could not find a good split plane -- retry with
			   more careful serial code just to be sure.. */
			execute_serially(bvh, node_idx, start, end, temp);
			return nullptr;
		}

		n_UINT left_count = bins.counts[best_index];
		int node_idx_left = node_idx + 1;
		int node_idx_right = node_idx + 2 * left_count;

		bvh.m_nodes[node_idx_left].bbox = bbox_left[best_index];
		bvh.m_nodes[node_idx_right].bbox = best_bbox_right;
		node.inner.rightChild = node_idx_right;
		node.inner.axis = axis;
		node.inner.flag = 0;

		std::atomic<n_UINT> offset_left(0),
			offset_right(bins.counts[best_index]);

		tbb::parallel_for(
			tbb::blocked_range<n_UINT>(0u, size, GRAIN_SIZE),
			[&](const tbb::blocked_range<n_UINT> &range) {
			n_UINT count_left = 0, count_right = 0;
			for (n_UINT i = range.begin(); i != range.end(); ++i) {
				n_UINT f = start[i];
				float centroid = bvh.getCentroid(f)[axis];
				int index = (int)((centroid - min) * inv_bin_size);
				(index <= best_index ? count_left : count_right)++;
			}
			n_UINT idx_l = offset_left.fetch_add(count_left);
			n_UINT idx_r = offset_right.fetch_add(count_right);
			for (n_UINT i = range.begin(); i != range.end(); ++i) {
				n_UINT f = start[i];
				float centroid = bvh.getCentroid(f)[axis];
				int index = (int)((centroid - min) * inv_bin_size);
				if (index <= best_index)
					temp[idx_l++] = f;
				else
					temp[idx_r++] = f;
			}
		}
		);
		memcpy(start, temp, size * sizeof(n_UINT));
		assert(offset_left == left_count && offset_right == size);

		/* Create an empty parent task */
		tbb::task& c = *new (allocate_continuation()) tbb::empty_task;
		c.set_ref_count(2);

		/* Post right subtree to scheduler */
		BVHBuildTask &b = *new (c.allocate_child())
			BVHBuildTask(bvh, node_idx_right, start + left_count,
				end, temp + left_count);
		spawn(b);

		/* Directly start working on left subtree */
		recycle_as_child_of(c);
		node_idx = node_idx_left;
		end = start + left_count;

		return this;
	}

	/// Single-threaded build function
	static void execute_serially(Accel &bvh, n_UINT node_idx, n_UINT *start, n_UINT *end, n_UINT *temp) {
		Accel::BVHNode &node = bvh.m_nodes[node_idx];
		n_UINT size = (n_UINT)(end - start);
		float best_cost = (float)INTERSECTION_COST * size;
		int64_t best_index = -1, best_axis = -1;
		float *left_areas = (float *)temp;

		/* Try splitting along every axis */
		for (int axis = 0; axis < 3; ++axis) {
			/* Sort all triangles based on their centroid positions projected on the axis */
			std::sort(start, end, [&](n_UINT f1, n_UINT f2) {
				return bvh.getCentroid(f1)[axis] < bvh.getCentroid(f2)[axis];
			});

			BoundingBox3f bbox;
			for (n_UINT i = 0; i < size; ++i) {
				n_UINT f = *(start + i);
				bbox.expandBy(bvh.getBoundingBox(f));
				left_areas[i] = (float)bbox.getSurfaceArea();
			}
			if (axis == 0)
				node.bbox = bbox;

			bbox.reset();

			/* Choose the best split plane */
			float tri_factor = INTERSECTION_COST / node.bbox.getSurfaceArea();
			for (n_UINT i = size - 1; i >= 1; --i) {
				n_UINT f = *(start + i);
				bbox.expandBy(bvh.getBoundingBox(f));

				float left_area = left_areas[i - 1];
				float right_area = bbox.getSurfaceArea();
				n_UINT prims_left = i;
				n_UINT prims_right = size - i;

				float sah_cost = 2.0f * TRAVERSAL_COST +
					tri_factor * (prims_left * left_area +
						prims_right * right_area);

				if (sah_cost < best_cost) {
					best_cost = sah_cost;
					best_index = i;
					best_axis = axis;
				}
			}
		}

		if (best_index == -1) {
			/* Splitting does not reduce the cost, make a leaf */
			node.leaf.flag = 1;
			node.leaf.start = (n_UINT)(start - bvh.m_indices.data());
			node.leaf.size = size;
			return;
		}

		std::sort(start, end, [&](n_UINT f1, n_UINT f2) {
			return bvh.getCentroid(f1)[best_axis] < bvh.getCentroid(f2)[best_axis];
		});

		n_UINT left_count = (n_UINT)best_index;
		n_UINT node_idx_left = node_idx + 1;
		n_UINT node_idx_right = node_idx + 2 * left_count;
		node.inner.rightChild = node_idx_right;
		node.inner.axis = best_axis;
		node.inner.flag = 0;

		execute_serially(bvh, node_idx_left, start, start + left_count, temp);
		execute_serially(bvh, node_idx_right, start + left_count, end, temp + left_count);
	}
};

void Accel::addMesh(Mesh *mesh) {
	m_meshes.push_back(mesh);
	m_meshOffset.push_back(m_meshOffset.back() + mesh->getTriangleCount());
	m_bbox.expandBy(mesh->getBoundingBox());
}

void Accel::clear() {
	m_meshes.clear();
	m_meshOffset.clear();
	m_meshOffset.push_back(0u);
	m_nodes.clear();
	m_indices.clear();
	m_bbox.reset();
	m_nodes.shrink_to_fit();
	m_meshes.shrink_to_fit();
	m_meshOffset.shrink_to_fit();
	m_indices.shrink_to_fit();
}

uint64_t float_to_morton(float value, int max_bits = 10) 
{
	uint64_t morton_code = 0;
	uint32_t value_scaled = static_cast<uint32_t>(value * ((1 << max_bits) - 1)); // Scale to the number of bits
	
	for (int i = 0; i < max_bits; ++i) {
		morton_code |= ((value_scaled >> i) & 1) << (3 * i); // 3D interleaving for Morton code
	}
	
	return morton_code;
}


// Function to compute the Morton code for a 3D point in [0, 1]
uint64_t morton3D(Point3f p, int bits = 10) 
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


void order_faces_by_morton(BVHmesh &mesh)
{
	Vector3f bbox_size = mesh.m_bbox.max - mesh.m_bbox.min;
	for (int i = 0; i < 3; ++i)
	{
		if (bbox_size[i] == 0)
			bbox_size[i] = 1.0f;
	}

	std::vector<uint64_t> morton_codes(mesh.size());
	for (uint32_t i = 0; i < mesh.size(); ++i)
	{
		uint32_t index1 = mesh.m_F(0, i);
		uint32_t index2 = mesh.m_F(1, i);
		uint32_t index3 = mesh.m_F(2, i);

		Vector3f p0 = mesh.m_V.col(index1);
		Vector3f p1 = mesh.m_V.col(index2);
		Vector3f p2 = mesh.m_V.col(index3);

		Vector3f centroid = (p0 + p1 + p2) / 3.0f;
		Vector3f normalized_centroid = (centroid - mesh.m_bbox.min);
		
		normalized_centroid = normalized_centroid.cwiseQuotient(bbox_size);
		morton_codes[i] = morton3D(normalized_centroid);
	}

	std::vector<uint32_t> indices(mesh.size());
	for (uint32_t i = 0; i < mesh.size(); ++i)
		indices[i] = i;

	std::sort(indices.begin(), indices.end(), [&morton_codes](uint32_t a, uint32_t b) {
		return morton_codes[a] < morton_codes[b];
	});


	MatrixXu newF = MatrixXu(3, mesh.size());

	for (uint32_t i = 0; i < mesh.size(); ++i)
	{
		newF.col(i) = mesh.m_F.col(indices[i]);
	}

	mesh.m_F = newF;
}


void Accel::compactMeshes() 
{
	globalMesh.clear();

	for (size_t i = 0; i < m_meshes.size(); ++i) 
	{
		globalMesh.append(m_meshes[i], i);
	}

	//order_faces_by_morton(globalMesh); // Does not work
}


void Accel::build() 
{
	compactMeshes();

	n_UINT size = getTriangleCount();
	if (size == 0)
		return;
	cout << "Constructing a SAH BVH (" << m_meshes.size()
		<< (m_meshes.size() == 1 ? " mesh, " : " meshes, ")
		<< size << " triangles) .. ";
	cout.flush();
	Timer timer;

	/* Conservative estimate for the total number of nodes */
	m_nodes.resize(2 * size);
	memset(m_nodes.data(), 0, sizeof(BVHNode) * m_nodes.size());
	m_nodes[0].bbox = m_bbox;
	m_indices.resize(size);

	cout << "Size of each node is " << sizeof(BVHNode);

	if ((sizeof(n_UINT) == 4) && (sizeof(BVHNode) != 32))
		throw NoriException("BVH Node is not packed! Investigate compiler settings.");

	for (n_UINT i = 0; i < size; ++i)
		m_indices[i] = i;

	n_UINT *indices = m_indices.data(), *temp = new n_UINT[size];
	BVHBuildTask& task = *new(tbb::task::allocate_root())
		BVHBuildTask(*this, 0u, indices, indices + size, temp);
	tbb::task::spawn_root_and_wait(task);
	delete[] temp;
	std::pair<float, n_UINT> stats = statistics();

	/* The node array was allocated conservatively and now contains
	   many unused entries -- do a compactification pass. */
	std::vector<BVHNode> compactified(stats.second);
	std::vector<n_UINT> skipped_accum(m_nodes.size());

	for (int64_t i = stats.second - 1, j = m_nodes.size(), skipped = 0; i >= 0; --i) {
		while (m_nodes[--j].isUnused())
			skipped++;
		BVHNode &new_node = compactified[i];
		new_node = m_nodes[j];
		skipped_accum[j] = (n_UINT)skipped;

		if (new_node.isInner()) {
			new_node.inner.rightChild = (n_UINT)
				(i + new_node.inner.rightChild - j -
				(skipped - skipped_accum[new_node.inner.rightChild]));
		}
	}
	cout << "done (took " << timer.elapsedString() << " and "
		<< memString(sizeof(BVHNode) * m_nodes.size() + sizeof(n_UINT)*m_indices.size())
		<< ", SAH cost = " << stats.first
		<< ")." << endl;

	m_nodes = std::move(compactified);
}

std::pair<float, n_UINT> Accel::statistics(n_UINT node_idx) const {
	const BVHNode &node = m_nodes[node_idx];
	if (node.isLeaf()) {
		return std::make_pair((float)BVHBuildTask::INTERSECTION_COST * node.leaf.size, 1u);
	}
	else {
		std::pair<float, n_UINT> stats_left = statistics(node_idx + 1u);
		std::pair<float, n_UINT> stats_right = statistics(node.inner.rightChild);
		float saLeft = m_nodes[node_idx + 1u].bbox.getSurfaceArea();
		float saRight = m_nodes[node.inner.rightChild].bbox.getSurfaceArea();
		float saCur = node.bbox.getSurfaceArea();
		float sahCost =
			2 * BVHBuildTask::TRAVERSAL_COST +
			(saLeft * stats_left.first + saRight * stats_right.first) / saCur;
		return std::make_pair(
			sahCost,
			stats_left.second + stats_right.second + 1u
		);
	}
}


Intersection Accel::fillIntersection(Point2f uv, n_UINT f, float t) const
{
	const Mesh *mesh = globalMesh.getOriginalMesh(f);

	Intersection its;
	its.uv = uv;
	its.mesh = mesh;
	its.t = t;
	its.f = f;

	/* Find the barycentric coordinates */
	Vector3f bary;
	bary << 1 - its.uv.sum(), its.uv;

	/* References to all relevant mesh buffers */
	const MatrixXf &V = mesh->getVertexPositions();
	const MatrixXf &N = mesh->getVertexNormals();
	const MatrixXf &UV = mesh->getVertexTexCoords();
	const MatrixXu &F = mesh->getIndices();

	/* Vertex indices of the triangle */
	n_UINT idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

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
	its.geoFrame = Frame(its.p, (p1 - p0).cross(p2 - p0).normalized());

	if (N.size() > 0) {
		/* Compute the shading frame. Note that for simplicity,
			the current implementation doesn't attempt to provide
			tangents that are continuous across the surface. That
			means that this code will need to be modified to be able
			use anisotropic BRDFs, which need tangent continuity */

		its.shFrame = Frame(its.p,
			(bary.x() * N.col(idx0) +
				bary.y() * N.col(idx1) +
				bary.z() * N.col(idx2)).normalized());
	}
	else {
		its.shFrame = its.geoFrame;
	}

	return its;
}

/*
bool rayIntersect(const MatrixXu *F, const MatrixXf *V, const Accel::BVHNode *nodes, const n_UINT *indices, 
				const Ray3f &_ray, Intersection &its) const 
{
	n_UINT node_idx = 0, stack_idx = 0, stack[64];

	its.t = std::numeric_limits<float>::infinity();

	// Use an adaptive ray epsilon 
	Ray3f ray(_ray);
	if (ray.mint == Epsilon)
	{
		ray.mint = std::max(ray.mint, ray.mint * ray.o.array().abs().maxCoeff());
	}
	
	if (ray.maxt < ray.mint)
		return false;

	bool foundIntersection = false;
	//n_UINT f = 0;

	while (true) 
	{
		const Accel::BVHNode &node = nodes[node_idx];

		if (node.isInner()) 
		{
			const n_UINT right = node.inner.rightChild;
			const n_UINT left = node_idx + 1;
			float tRight, tLeft;
			bool right_intersected = nodes[right].bbox.rayIntersect(ray, tRight);
			bool left_intersected  = nodes[left].bbox.rayIntersect(ray, tLeft);

			if (left_intersected && right_intersected)
			{
				if (tRight > tLeft)
				{
					stack[stack_idx++] = right;
					node_idx = left;
				}
				else
				{
					stack[stack_idx++] = left;
					node_idx = right;
				}

				assert(stack_idx < 64);
			}
			else if (left_intersected)
			{
				node_idx = left;
			}
			else if (right_intersected)
			{
				node_idx = right;
			}
			else
			{
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
		else 
		{
			for (n_UINT i = node.start(), end = node.end(); i < end; ++i) 
			{
				n_UINT idx = indices[i];

				float u, v, t;
				if (rayIntersect(F, V, idx, ray, u, v, t)) 
				{
					ray.maxt = its.t = t;
					its.uv = Point2f(u, v);

					foundIntersection = true;
					its.f = idx; //= f;
				}
			}
			if (stack_idx == 0)
				break;
			node_idx = stack[--stack_idx];
			continue;
		}
	}

	return foundIntersection;
}

bool rayIntersect(const MatrixXu *F, const MatrixXf *V, 
				n_UINT index, const Ray3f &ray, float &u, float &v, float &t) const	
{
	n_UINT i0 = (*F)(0, index), i1 = (*F)(1, index), i2 = (*F)(2, index);
	const Eigen::Vector3f p0 = (*V).col(i0), p1 = (*V).col(i1), p2 = (*V).col(i2);

	// Find vectors for two edges sharing v[0] 
	Eigen::Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

	// Begin calculating determinant - also used to calculate U parameter 
	Eigen::Vector3f pvec = ray.d.cross(edge2);

	// If determinant is near zero, ray lies in plane of triangle 
	float det = edge1.dot(pvec);

	if (det > -1e-8f && det < 1e-8f)
		return false;
	float inv_det = 1.0f / det;

	// Calculate distance from v[0] to ray origin 
	Eigen::Vector3f tvec = ray.o - p0;

	// Calculate U parameter and test bounds 
	u = tvec.dot(pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	// Prepare to test V parameter 
	Eigen::Vector3f qvec = tvec.cross(edge1);

	// Calculate V parameter and test bounds 
	v = ray.d.dot(qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	// Ray intersects triangle -> compute t 
	t = edge2.dot(qvec) * inv_det;

	return t >= ray.mint && t <= ray.maxt;
}
};
*/


void Accel::rayIntersect(const std::vector<bool> &mask, 
						const std::vector<Ray3f> &ray, 
						std::vector<Intersection> &its,
						std::vector<bool> &hit) const
{
	/*
	for (n_UINT i = 0; i < ray.size(); ++i) 
	{
		if (mask[i] && rayIntersect(ray[i], its[i]))
		{
			its[i] = fillIntersection(its[i].uv, its[i].f, its[i].t);
			hit[i] = true;
		}
		else
			hit[i] = false;
	}
	return;
	*/

	// Indices
	auto indices = std::ranges::views::iota(size_t{0}, ray.size());

	// transform indices into bools
	std::transform(std::execution::par_unseq,
					mask.begin(), mask.end(), 
					indices.begin(), hit.begin(), 
					[&](bool b, int i) 
						{ 
							if (!b)
								return false;
							else
								//return true;
								return rayIntersectImpl(ray[i], its[i], false);	
						});

	for (n_UINT i = 0; i < its.size(); ++i)
	{
		if (hit[i])
			its[i] = fillIntersection(its[i].uv, its[i].f, its[i].t);
	}
}

bool Accel::rayIntersect(const Ray3f &_ray, Intersection &its, bool shadowRay) const 
{
	bool res = rayIntersectImpl(_ray, its, shadowRay);
	if (res)
		its = fillIntersection(its.uv, its.f, its.t);

	return res;
}

bool Accel::rayIntersectImpl(const Ray3f &_ray, Intersection &its, bool shadowRay) const 
{
	n_UINT node_idx = 0, stack_idx = 0, stack[64];
	bool foundIntersection = false;

	its.t = std::numeric_limits<float>::infinity();

	// Use an adaptive ray epsilon 
	Ray3f ray(_ray);
	
	if (ray.mint == Epsilon)
	{
		const auto absP = ray.o.array().abs();
		float maxCoeff = std::max(absP.x(), std::max(absP.y(), absP.z()));
		ray.mint = std::max(ray.mint, ray.mint * maxCoeff);
	}
	
	if (m_nodes.empty() || ray.maxt < ray.mint)
		return false;

	
	//n_UINT f = 0;

	while (true) {
		const BVHNode &node = m_nodes[node_idx];

		if (node.isInner()) 
		{
			const n_UINT right = node.inner.rightChild;
			const n_UINT left = node_idx + 1;
			float tRight, tLeft;
			bool right_intersected = m_nodes[right].bbox.rayIntersect(ray, tRight);
			bool left_intersected  = m_nodes[left].bbox.rayIntersect(ray, tLeft);

			if (left_intersected && right_intersected)
			{
				if (tRight > tLeft)
				{
					stack[stack_idx++] = right;
					node_idx = left;
				}
				else
				{
					stack[stack_idx++] = left;
					node_idx = right;
				}

				assert(stack_idx < 64);
			}
			else if (left_intersected)
			{
				node_idx = left;
			}
			else if (right_intersected)
			{
				node_idx = right;
			}
			else
			{
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
		else 
		{
			for (n_UINT i = node.start(), end = node.end(); i < end; ++i) 
			{
				n_UINT idx = m_indices[i];

				float u, v, t;
				if (globalMesh.rayIntersect(idx, ray, u, v, t)) 
				{
					ray.maxt = its.t = t;
					its.uv = Point2f(u, v);
					its.f = idx; //= f;

					if (shadowRay)
						return true;

					foundIntersection = true;
				}
			}
			if (stack_idx == 0)
				break;
			node_idx = stack[--stack_idx];
			continue;
		}
	}

	
	return foundIntersection;
}


bool Accel::rayProbe(const Ray3f &_ray, std::vector<Intersection> &its) const 
{
	n_UINT node_idx = 0, stack_idx = 0, stack[64];

	/* Use an adaptive ray epsilon */
	Ray3f ray(_ray);
	if (ray.mint == Epsilon)
		ray.mint = std::max(ray.mint, ray.mint * ray.o.array().abs().maxCoeff());

	if (m_nodes.empty() || ray.maxt < ray.mint)
		return false;

	n_UINT f = 0;

	while (true) 
	{
		const BVHNode &node = m_nodes[node_idx];

		if (node.isInner()) 
		{
			const n_UINT right = node.inner.rightChild;
			const n_UINT left = node_idx + 1;
			float tRight, tLeft;
			bool right_intersected = m_nodes[right].bbox.rayIntersect(ray, tRight);
			bool left_intersected  = m_nodes[left].bbox.rayIntersect(ray, tLeft);

			if (left_intersected && right_intersected)
			{
				if (tRight > tLeft)
				{
					stack[stack_idx++] = right;
					node_idx = left;
				}
				else
				{
					stack[stack_idx++] = left;
					node_idx = right;
				}

				assert(stack_idx < 64);
			}
			else if (left_intersected)
			{
				node_idx = left;
			}
			else if (right_intersected)
			{
				node_idx = right;
			}
			else
			{
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
		else 
		{
			for (n_UINT i = node.start(), end = node.end(); i < end; ++i) 
			{
				n_UINT idx = m_indices[i];
				const Mesh *mesh = globalMesh.getOriginalMesh(idx);

				float u, v, t; 
				if (mesh->rayIntersect(idx, ray, u, v, t)) 
				{
					its.push_back(fillIntersection(Point2f(u, v), f, t));
				}
			}

			if (stack_idx == 0)
				break;

			node_idx = stack[--stack_idx];
			continue;
		}
	}


	return its.size() > 0;
}



NORI_NAMESPACE_END


/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    v1 - Dec 2020
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


#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { }

Mesh::~Mesh() {
    m_pdf.clear();
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }

    meshArea = 0.0f;

    for (Eigen::Index i = 0; i < m_F.cols(); ++i)
        meshArea += surfaceArea(i);

    std::cout << "Mesh triangles: " << m_F.cols() << ", area: " << meshArea << std::endl;

    // Build triangle distribution
    m_pdf.clear();
    m_pdf.reserve(m_F.cols());
    for (uint32_t i = 0; i < m_F.cols(); ++i) {
        m_pdf.append(surfaceArea(i) / meshArea);
    }

    m_pdf.normalize();
}

float Mesh::surfaceArea(n_UINT index) const
{
    n_UINT i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}



int Mesh::sampleTriangle(Sampler *sampler, float &pdf)  const
{
    uint32_t triangleIndex = m_pdf.sample(sampler->next1D(), pdf);
    //pdf = 1.0f / pdf;

    return triangleIndex;
}

float Mesh::sampleTrianglePdf(uint32_t index) const 
{
    return m_pdf[index];
}

Normal3f Mesh::getNormal(n_UINT index, const Point2f &uv) const
{
    n_UINT i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    Vector3f n = (p1 - p0).cross(p2 - p0);
    n.normalize();

    if (m_N.cols() > 0)
    {
        const Normal3f n0 = m_N.col(i0), n1 = m_N.col(i1), n2 = m_N.col(i2);
        n = (1 - uv[0] - uv[1]) * n0 + uv[0] * n1 + uv[1] * n2;
    }

    return Normal3f(n);
}

void Mesh::samplePosition(Sampler *sampler, Point3f &p, 
        Normal3f &n, Point2f &uv, 
        float &pdf, n_UINT &triangleId) const
{
    // Choose a random triangle
    float trianglePdf;
    uint32_t triangleIndex = sampleTriangle(sampler, trianglePdf);
    triangleId = triangleIndex;

    float s1 = sampler->next1D();
    float s2 = sampler->next1D();

    float alpha = 1 - sqrt(1 - s1);
    float beta = s2 * sqrt(1 - s1);

    // Compute the normal at the sampled point
    n = getNormal(triangleIndex, Point2f(alpha, beta));

    // Compute the position
    p = uvTo3D(triangleIndex, Point2f(alpha, beta));

    // Compute the PDF
    pdf = 1.0f / meshArea;
}

void Mesh::samplePositions(Sampler *sampler, 
        std::vector<Point3f> &points, 
        std::vector<Normal3f> &normals,
        std::vector<Point2f> &uvs, 
        std::vector<float> &pdfs, 
        std::vector<n_UINT> &triangleIds, 
        size_t samplesPerTriangle) const
{
    long long n_samples = samplesPerTriangle * m_F.cols();
    points.resize(n_samples);
    normals.resize(n_samples);
    uvs.resize(n_samples);
    pdfs.resize(n_samples);
    triangleIds.resize(n_samples);

    for (Eigen::Index i = 0; i < m_F.cols(); ++i)
    {
        for (size_t j = 0; j < samplesPerTriangle; ++j)
        {
            Point2f s = sampler->next2D();

            float alpha = 1 - sqrt(1 - s.x());
            float beta = s.y() * sqrt(1 - s.x());

            // Compute the normal at the sampled point
            Normal3f n = getNormal(i, Point2f(alpha, beta));

            // Compute the position
            Point3f p = uvTo3D(i, Point2f(alpha, beta));

            // Compute the PDF
            float pdf = 1 / surfaceArea(i);

            size_t index = i * samplesPerTriangle + j;
            points[index] = p;
            normals[index] = n;
            uvs[index] = Point2f(alpha, beta);
            pdfs[index] = pdf;
            triangleIds[index] = i;
        }
    }
}

float Mesh::pdf(const Point3f &p, n_UINT triangleId) const
{
    return (1.0f / surfaceArea(triangleId)) * sampleTrianglePdf(triangleId);
}

float Mesh::pdf(const Point3f &p) const 
{
    uint32_t triangleIndex = getTriangleIndex(p);

    return pdf(p, triangleIndex);
}

Point3f Mesh::uvTo3D(n_UINT index, const Point2f &uv) const
{
    n_UINT i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return uv[0] * p0 + uv[1] * p1 + (1 - uv[0] - uv[1]) * p2;
}



Point2f Mesh::uvFrom3D(n_UINT index, const Point3f &p) const
{
    n_UINT i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    Vector3f n = (p1 - p0).cross(p2 - p0);
    float area = n.norm();

    Vector3f n0 = (p1 - p).cross(p2 - p);
    Vector3f n1 = (p2 - p).cross(p0 - p);
    //Vector3f n2 = (p0 - p).cross(p1 - p);

    float alpha = n0.norm() / area;
    float beta = n1.norm() / area;
    //float gamma = n2.norm() / area;

    return Point2f(alpha, beta);
}


int Mesh::getTriangleIndex(Point3f p) const
{
    int index = -1;

    for (uint32_t i = 0; i < m_F.cols(); ++i) 
    {
        //BoundingBox3f bbox = getBoundingBox(i);
        //if (bbox.contains(p)) 
        {
            uint32_t i0 = m_F(0, i), i1 = m_F(1, i), i2 = m_F(2, i);

            Vector3f v0 = m_V.col(i0).template head<3>();
            Vector3f v1 = m_V.col(i1).template head<3>();
            Vector3f v2 = m_V.col(i2).template head<3>();

            Vector3f n = (v1 - v0).cross(v2 - v0);
            float area = n.norm();

            Vector3f n0 = (v1 - p).cross(v2 - p);
            Vector3f n1 = (v2 - p).cross(v0 - p);
            Vector3f n2 = (v0 - p).cross(v1 - p);

            float alpha = n0.norm() / area;
            float beta = n1.norm() / area;
            float gamma = n2.norm() / area;

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                index = i;
                break;
            }
        }
    }

    return index;
}

bool Mesh::rayIntersect(n_UINT index, const Ray3f &ray, float &u, float &v, float &t) const {
    n_UINT i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(n_UINT index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(n_UINT index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}



void Mesh::addChild(NoriObject *obj, const std::string& name) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END

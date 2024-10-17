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

#include <nori/mesh.h>
#include <nori/timer.h>
#include <filesystem/resolver.h>
#include <unordered_map>
#include <fstream>

NORI_NAMESPACE_BEGIN

/**
 * \brief Loader for Wavefront OBJ triangle meshes
 */
class WavefrontOBJ : public Mesh 
{
        
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
    uint64_t morton3D(Point3f p, int bits = 10) {
        uint64_t x_morton = float_to_morton(p.x(), bits);
        uint64_t y_morton = float_to_morton(p.y(), bits);
        uint64_t z_morton = float_to_morton(p.z(), bits);  

        uint64_t morton_code = 0;
        for (int i = 0; i < bits; ++i) {
            morton_code |= ((x_morton >> i) & 1) << (3 * i);
            morton_code |= ((y_morton >> i) & 1) << (3 * i + 1);
            morton_code |= ((z_morton >> i) & 1) << (3 * i + 2);
        }

        return morton_code;
    }

    void order_faces_by_morton()
    {
        std::vector<uint64_t> morton_codes(m_F.cols());
        for (uint32_t i = 0; i < m_F.cols(); ++i)
        {
            uint32_t x = m_F(0, i);
            uint32_t y = m_F(1, i);
            uint32_t z = m_F(2, i);
            Point3f centroid = (m_V.col(x) + m_V.col(y) + m_V.col(z)) / 3.0f;
            
            Vector3f bbox_size = m_bbox.max - m_bbox.min;
            Point3f normalized_centroid = (centroid - m_bbox.min);
            
            normalized_centroid.x() /= bbox_size.x();
            normalized_centroid.y() /= bbox_size.y();
            normalized_centroid.z() /= bbox_size.z();

            morton_codes[i] = morton3D(normalized_centroid);
        }

        std::vector<uint32_t> indices(m_F.cols());
        for (uint32_t i = 0; i < m_F.cols(); ++i)
            indices[i] = i;

        std::sort(indices.begin(), indices.end(), [&morton_codes](uint32_t a, uint32_t b) {
            return morton_codes[a] < morton_codes[b];
        });

        MatrixXu F(3, m_F.cols());
        for (uint32_t i = 0; i < m_F.cols(); ++i)
            F.col(i) = m_F.col(indices[i]);

        m_F = F;
    }

public:
    WavefrontOBJ(const PropertyList &propList) {
        typedef std::unordered_map<OBJVertex, uint32_t, OBJVertexHash> VertexMap;

        filesystem::path filename =
            getFileResolver()->resolve(propList.getString("filename"));

        std::ifstream is(filename.str());
        if (is.fail())
            throw NoriException("Unable to open OBJ file \"%s\"!", filename);
        Transform trafo = propList.getTransform("toWorld", Transform());

        cout << "Loading \"" << filename << "\" .. ";
        cout.flush();
        Timer timer;

        std::vector<Vector3f>   positions;
        std::vector<Vector2f>   texcoords;
        std::vector<Vector3f>   normals;
        std::vector<uint32_t>   indices;
        std::vector<OBJVertex>  vertices;
        VertexMap vertexMap;

        std::string line_str;
        while (std::getline(is, line_str)) {
            std::istringstream line(line_str);

            std::string prefix;
            line >> prefix;

            if (prefix == "v") {
                Point3f p;
                line >> p.x() >> p.y() >> p.z();
                p = trafo * p;
                m_bbox.expandBy(p);
                positions.push_back(p);
            } else if (prefix == "vt") {
                Point2f tc;
                line >> tc.x() >> tc.y();
                texcoords.push_back(tc);
            } else if (prefix == "vn") {
                Normal3f n;
                line >> n.x() >> n.y() >> n.z();
                normals.push_back((trafo * n).normalized());
            } else if (prefix == "f") {
                std::string v1, v2, v3, v4;
                line >> v1 >> v2 >> v3 >> v4;
                OBJVertex verts[6];
                int nVertices = 3;

                verts[0] = OBJVertex(v1);
                verts[1] = OBJVertex(v2);
                verts[2] = OBJVertex(v3);

                if (!v4.empty()) {
                    /* This is a quad, split into two triangles */
                    verts[3] = OBJVertex(v4);
                    verts[4] = verts[0];
                    verts[5] = verts[2];
                    nVertices = 6;
                }
                /* Convert to an indexed vertex list */
                for (int i=0; i<nVertices; ++i) {
                    const OBJVertex &v = verts[i];
                    VertexMap::const_iterator it = vertexMap.find(v);
                    if (it == vertexMap.end()) {
                        vertexMap[v] = (uint32_t) vertices.size();
                        indices.push_back((uint32_t) vertices.size());
                        vertices.push_back(v);
                    } else {
                        indices.push_back(it->second);
                    }
                }
            }
        }
        cout << "End reading \"" << filename << "\" .. ";


        m_F.resize(3, indices.size()/3);
        memcpy(m_F.data(), indices.data(), sizeof(uint32_t)*indices.size());

        m_V.resize(3, vertices.size());
        for (uint32_t i = 0; i < vertices.size(); ++i)
        {
            m_V.col(i) = positions.at(vertices[i].p - 1);
        }

        cout << "Passed with vertices \"" << filename << "\" .. ";

        if (!normals.empty()) {
            m_N.resize(3, vertices.size());
            for (uint32_t i=0; i<vertices.size(); ++i)
                m_N.col(i) = normals.at(vertices[i].n-1);
        }

        cout << "Passed with normals \"" << filename << "\" .. ";

        if (!texcoords.empty()) {
            m_UV.resize(2, vertices.size());
            for (uint32_t i=0; i<vertices.size(); ++i)
                m_UV.col(i) = texcoords.at(vertices[i].uv-1);
        }

        cout << "Passed with UVs \"" << filename << "\" .. ";

        //order_faces_by_morton();

        m_name = filename.str();
        cout << "done. (V=" << m_V.cols() << ", F=" << m_F.cols() << ", took "
             << timer.elapsedString() << " and "
             << memString(m_F.size() * sizeof(uint32_t) +
                          sizeof(float) * (m_V.size() + m_N.size() + m_UV.size()))
             << ")" << endl;
    }

protected:
    /// Vertex indices used by the OBJ format
    struct OBJVertex {
        uint32_t p = (uint32_t) -1;
        uint32_t n = (uint32_t) -1;
        uint32_t uv = (uint32_t) -1;

        inline OBJVertex() { }

        inline OBJVertex(const std::string &string) {
            std::vector<std::string> tokens = tokenize(string, "/", true);

            if (tokens.size() < 1 || tokens.size() > 3)
                throw NoriException("Invalid vertex data: \"%s\"", string);

            p = toUInt(tokens[0]);

            if (tokens.size() >= 2 && !tokens[1].empty())
                uv = toUInt(tokens[1]);

            if (tokens.size() >= 3 && !tokens[2].empty())
                n = toUInt(tokens[2]);
        }

        inline bool operator==(const OBJVertex &v) const {
            return v.p == p && v.n == n && v.uv == uv;
        }
    };

    /// Hash function for OBJVertex
    struct OBJVertexHash : std::unary_function<OBJVertex, size_t> {
        std::size_t operator()(const OBJVertex &v) const {
            size_t hash = std::hash<uint32_t>()(v.p);
            hash = hash * 37 + std::hash<uint32_t>()(v.uv);
            hash = hash * 37 + std::hash<uint32_t>()(v.n);
            return hash;
        }
    };
};

NORI_REGISTER_CLASS(WavefrontOBJ, "obj");
NORI_NAMESPACE_END

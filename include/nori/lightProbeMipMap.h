#pragma once

#include <nori/common.h>
#include <nori/bbox.h>
#include <fstream>


NORI_NAMESPACE_BEGIN

class LightProbeMipMap
{
    int levels;

    std::vector<Eigen::MatrixXf> pyramid;
    BoundingBox2f domain;
    
    public:

    // From an eigen matrix of floats, create a MipMap
    LightProbeMipMap (const Eigen::MatrixXf &image, const BoundingBox2f &domain);
    LightProbeMipMap () {}

    int getDepth() const { return levels; }

    // Get the value of the MipMap at a given point
    float l_lookup (const Point2f &st, float level = 0.0f) const;
    float l_lookup (const Point2i &st, float level = 0.0f) const { return l_lookup(Point2f(st.x(), st.y()), level); }
    float l_lookup (float level, int x, int y) const;

    Point2i sampleDiscrete (Point2f sample, float &pdf) const;
    float discretePdf(Point2i p) const;

    static Point2f staticSample (Point2f sample);
    static float staticSamplePdf (Point2f p);

    Point2f sampleContinuous(Point2f sample, float &pdf) const;
    float continuousPdf(Point2f p) const;

    static LightProbeMipMap *loadMatrixFromFile(const std::string &filename)
    {
        std::ifstream file(filename);

        std::vector<double> values;
        int rows = 0;
        int cols = 0;
        double value;

        // Read data from file line by line
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            int temp_cols = 0;
            while (ss >> value) {
                values.push_back(value);
                ++temp_cols;
            }
            if (cols == 0) {
                cols = temp_cols;  // Set the number of columns based on the first row
            }
            ++rows;
        }
        file.close();

        // Convert the vector to an Eigen matrix
        Eigen::MatrixXd matrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix(i, j) = std::pow(values[i * cols + j], 2);
            }
        }

        std::cout << "Loaded matrix from file: " << filename << " with dimensions: " << matrix.rows() << "x" << matrix.cols() << std::endl;

        // Initialize
        return new LightProbeMipMap (matrix.cast<float>(), BoundingBox2f(Point2f(0), Point2f(1, 1)));
    }
};

NORI_NAMESPACE_END
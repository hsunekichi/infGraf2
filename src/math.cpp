#include <nori/math.h>


#include <nori/common.h>
#include <nori/transform.h>
#include <nori/bitmap.h>
#include <nori/vector.h>

#include <chrono>
#include <thread>

NORI_NAMESPACE_BEGIN

// Computes the 3D gradient of slice 1,
//  x and y are slice1 gradients and 
//  and Z is the gradient going to slice 2 
Eigen::MatrixXf Math::finiteDifferenceGradientModulus3D(const Eigen::MatrixXf& slice1, const Eigen::MatrixXf& slice2) 
{
    // Get the dimensions of the input matrices
    int rows = slice1.rows();
    int cols = slice1.cols();

    // Create matrices to store the gradients
    Eigen::MatrixXf gradX = Eigen::MatrixXf::Zero(rows, cols);
    Eigen::MatrixXf gradY = Eigen::MatrixXf::Zero(rows, cols);
    Eigen::MatrixXf gradZ = Eigen::MatrixXf::Zero(rows, cols);
    Eigen::MatrixXf gradientModulus = Eigen::MatrixXf::Zero(rows, cols);

    // Apply finite difference method for x, y, and z gradients
    for (int i = 1; i < rows - 1; ++i) 
    {
        for (int j = 1; j < cols - 1; ++j) 
        {
            // Compute the x-gradient using finite differences
            float Gx = (slice1(i+1, j) - slice1(i-1, j)) / 2.0;

            // Compute the y-gradient using finite differences
            float Gy = (slice1(i, j+1) - slice1(i, j-1)) / 2.0;

            // Compute the z-gradient as the difference between the two slices
            float Gz = (slice2(i, j) - slice1(i, j)) / 2.0;

            // Store the gradients
            gradX(i, j) = Gx;
            gradY(i, j) = Gy;
            gradZ(i, j) = Gz;

            // Compute the gradient modulus
            gradientModulus(i, j) = std::sqrt(Gx * Gx + Gy * Gy); //  + Gz * Gz
        }
    }

    return gradientModulus;
}

// Function to compute the Sobel gradient modulus in 3D (two slices)
Eigen::MatrixXf Math::sobelGradientModulus3D(const Eigen::MatrixXf& slice1, const Eigen::MatrixXf& slice2) 
{
    // Sobel kernels for x and y gradients
    Eigen::Matrix3f sobelX;
    sobelX << -1, 0, 1,
            -2, 0, 2,
            -1, 0, 1;
    
    Eigen::Matrix3f sobelY;
    sobelY << -1, -2, -1,
            0,  0,  0,
            1,  2,  1;

    // Get the dimensions of the input matrices
    int rows = slice1.rows();
    int cols = slice1.cols();

    // Create matrices to store the gradients
    //Eigen::MatrixXf gradX = Eigen::MatrixXf::Zero(rows, cols);
    //Eigen::MatrixXf gradY = Eigen::MatrixXf::Zero(rows, cols);
    //Eigen::MatrixXf gradZ = Eigen::MatrixXf::Zero(rows, cols);
    Eigen::MatrixXf gradientModulus = Eigen::MatrixXf::Zero(rows, cols);

    // Apply Sobel filter for x and y gradients
    for (int i = 1; i < rows - 1; ++i) 
    {
        for (int j = 1; j < cols - 1; ++j) 
        {
            // Extract the 3x3 region from the first slice (slice1)
            Eigen::MatrixXf region = slice1.block<3,3>(i-1, j-1);

            // Compute the x and y gradients using Sobel filters
            float Gx = (region.cwiseProduct(sobelX)).sum();
            float Gy = (region.cwiseProduct(sobelY)).sum();

            // Compute the z-gradient using central difference
            float Gz = (slice2(i, j) - slice1(i, j)) / 2.0;

            // Store the gradients
            //gradX(i, j) = Gx;
            //gradY(i, j) = Gy;
            //gradZ(i, j) = Gz;

            // Compute the gradient modulus
            gradientModulus(i, j) = std::sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
        }
    }

    return gradientModulus;
}


float Math::findRoot(std::function<float(float)> f, float x0) 
{
    // Initial guess
    float x1 = x0 + 0.01;

    // Compute the function value at the initial guess
    float f0 = f(x0);

    // Compute the function value at the next guess
    float f1 = f(x1);

    // Maximum number of iterations
    int maxIter = 100;

    // Tolerance
    float tol = 1e-6;

    // Iteration counter
    int iter = 0;

    // Loop until the function value is close to zero
    while (std::abs(f1) > tol && iter < maxIter) 
    {
        // Compute the next guess
        float x2 = (x0 * f1 - x1 * f0) / (f1 - f0);

        // Update the initial guess
        x0 = x1;
        f0 = f1;

        // Update the function value
        x1 = x2;
        f1 = f(x1);

        // Update the iteration counter
        iter++;
    }

    if (f1 > tol) 
    {
        throw NoriException("Root finding did not converge!");
    }

    return x1;
}

float Math::matrixVariance(const Eigen::MatrixXf& block, float &mean) 
{
    // Compute the mean of the block
    mean = block.mean();

    // Compute the variance of the block
    Eigen::ArrayXXf difference = block.array() - mean;
    float variance = difference.square().sum() / (block.size() - 1);

    return variance;
}

float Math::matrixStdDeviation(const Eigen::MatrixXf& block, float &mean) 
{
    float variance = matrixVariance(block, mean);
    return std::sqrt(variance);
}

float Math::blockPercentualStdDeviation(const Eigen::MatrixXf& block) 
{
    // Compute the mean of the block
    float mean;
    float variance = matrixVariance(block, mean);

    if (mean > 0)
        return (std::sqrt(variance) / mean);
    else
        return 0;
}

Eigen::MatrixXf Math::weibull_stretched_exponential(const Eigen::MatrixXf &m, float alpha, float beta)
{
    // Get the dimensions of the input matrix
    int rows = m.rows();
    int cols = m.cols();

    // Create matrix to store the equalized values
    Eigen::MatrixXf equalized = Eigen::MatrixXf::Zero(rows, cols);

    // Apply equalization
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            equalized(i, j) = weibull_stretched_exponential(m(i, j), alpha, beta);
        }
    }

    return equalized;
}


Eigen::MatrixXf Math::equalize(const Eigen::MatrixXf &m)
{
    // Get the dimensions of the input matrix
    int rows = m.rows();
    int cols = m.cols();

    // Create matrix to store the equalized values
    Eigen::MatrixXf equalized = Eigen::MatrixXf::Zero(rows, cols);

    // Compute maximum value
    float max = m.maxCoeff();

    // Apply equalization
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            equalized(i, j) = m(i, j) / max;
        }
    }

    return equalized;
}

Eigen::MatrixXf Math::computeVariance3x3(const Eigen::MatrixXf& slice1, const Eigen::MatrixXf& slice2) 
{
    // Get the dimensions of the input matrices
    int rows = slice1.rows();
    int cols = slice1.cols();

    // Create matrix to store the variances
    Eigen::MatrixXf variances = Eigen::MatrixXf::Zero(rows, cols);

    // Apply Sobel filter for x and y gradients
    for (int i = 1; i < rows - 1; ++i) 
    {
        for (int j = 1; j < cols - 1; ++j) 
        {
            // Extract the 3x3 region from the first slice (slice1)
            Eigen::MatrixXf region1 = slice1.block<3,3>(i-1, j-1);
            //Eigen::MatrixXf region2 = slice2.block<3,3>(i-1, j-1);

            // Compute variance of the 3x3 region
            float stdDev1 = blockPercentualStdDeviation(region1);
            //float stdDev2 = blockPercentualStdDeviation(region2);

            //std::cout << "stdDev1: " << stdDev1 << std::endl;
            variances(i, j) = std::abs(stdDev1);
        }
    }


    // Compute variance
    return weibull_stretched_exponential(variances);
}

Eigen::MatrixXf Math::luminanceMatrix(const Bitmap *_bmp_image)
{
    const Bitmap &bmp_image = *_bmp_image;

    // Get the dimensions of the image
    size_t imgX = bmp_image.cols();
    size_t imgY = bmp_image.rows();

    Eigen::MatrixXf luminance = Eigen::MatrixXf::Zero(imgX, imgY);
    for (size_t i = 0; i < imgX; ++i)
    {
        for (size_t j = 0; j < imgY; ++j)
        {
            luminance(i, j) = bmp_image(j, i).getLuminance();
        }
    }

    return luminance;
}

NORI_NAMESPACE_END
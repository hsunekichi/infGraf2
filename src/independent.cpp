/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/sampler.h>
#include <nori/block.h>
#include <pcg32.h>
#include <random>

NORI_NAMESPACE_BEGIN

/**
 * Independent sampling - returns independent uniformly distributed
 * random numbers on <tt>[0, 1)x[0, 1)</tt>.
 *
 * This class is essentially just a wrapper around the pcg32 pseudorandom
 * number generator. For more details on what sample generators do in
 * general, refer to the \ref Sampler class.
 */
class Independent : public Sampler {
public:
    Independent(const PropertyList &propList) 
    {
        m_sampleCount = propList.getInteger("sampleCount", 1);
        _nPasses = propList.getInteger("nPasses", 1);
    }

    virtual ~Independent() { }

    std::unique_ptr<Sampler> clone() const 
    {
        std::unique_ptr<Independent> cloned(new Independent());
        cloned->m_sampleCount = m_sampleCount;
        cloned->imgX = imgX;
        cloned->imgY = imgY;
        cloned->_weightedPassesCount = _weightedPassesCount;

        cloned->m_random = m_random;
        cloned->samplesPerPixel = samplesPerPixel;
        cloned->weightedPassTotalSamples = weightedPassTotalSamples;

        cloned->currentPass = currentPass;
        cloned->prevLuminance = prevLuminance;
        cloned->mipMap = mipMap;
        cloned->_nPasses = _nPasses;

        return cloned;
    }

    void _initialize(const Bitmap *initialImage) 
    {
        samplesPerPixel = Eigen::MatrixXf::Zero(imgX, imgY);
        currentPass = 0;
        prevLuminance = Math::luminanceMatrix(initialImage);
        setConstantSpp(m_sampleCount);
    }

    void prepare(const ImageBlock &block) 
    {
        // Optionally use a different entropy source for the second seed
        std::random_device rd;
        unsigned int seed1 = rd();  // Gets a random seed from the random device
        unsigned int seed2 = rd();  // Gets a random seed from the random device

        m_random.seed(
            seed1,
            seed2
        );
    }

    Point2f generatePixelOffset(Point2i pixelId) {
        return next2D();
    }

    void generate() { /* No-op for this sampler */ }
    void advance()  { /* No-op for this sampler */ }

    float next1D() {
        return m_random.nextFloat();
    }
    
    Point2f next2D() {
        return Point2f(
            m_random.nextFloat(),
            m_random.nextFloat()
        );
    }

    std::string toString() const {
        return tfm::format("Independent[sampleCount=%i]", m_sampleCount);
    }

    void initializeWeightedSamples (const Eigen::MatrixXf &image)
    {
        // Compute gradient of the image
        //Eigen::MatrixXf pixelWeights = Math::sobelGradientModulus3D(image, prevImage);
        Eigen::MatrixXf pixelWeights = Math::computeVariance3x3(image, prevLuminance);

        // Generate weighted mip map to sample
        mipMap = LightProbeMipMap(pixelWeights, BoundingBox2f(Point2f(0), Point2f(imgX, imgY)));

        auto start = std::chrono::high_resolution_clock::now();
        // Generate number of samples per pixel
        for (size_t i = 0; i < weightedPassTotalSamples; ++i)
        {
            float pdf;
            Point2i pixel = mipMap.sampleDiscrete(next2D(), pdf);

            samplesPerPixel(pixel.x(), pixel.y()) += 1;
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time to generate samples: " << duration.count() << "ms" << std::endl;

        matrixToFile(pixelWeights, "spp" + std::to_string(currentPass) + ".txt");
    }

    void setConstantSpp(size_t spp)
    {
        samplesPerPixel.setZero();

        for (size_t i = 0; i < imgX; ++i)
        {
            for (size_t j = 0; j < imgY; ++j)
            {
                samplesPerPixel(i, j) = spp;
            }
        }
    }


    void matrixToFile(const Eigen::MatrixXf &m, const std::string &filename)
    {
        std::ofstream file(filename, std::ios::out);

        file << m << std::endl;

        file.close();
    }

    void next_pass() override
    { 
        currentPass++; 

        /*
        samplesPerPixel.setZero();

        // Convert RGB to luminance
        auto image = Math::luminanceMatrix(renderedImage);
        initializeWeightedSamples(image);

        prevLuminance = image;
        */
    }

    /// Return the number of configured pixel samples
    size_t getSampleCount(Point2i pixel) const { return samplesPerPixel(pixel.x(), pixel.y()); }
    size_t getSampleCount() const 
    { 
        if (_nPasses == 1 || currentPass > 0)
        {
            return m_sampleCount;
        }
        else
        {
            return 1;
        }
    }

    float getSamplePdf(Point2i pixel) const
    {
        if (currentPass == 0)
        {
            return 1.0f;
        }
        else
        {
            return mipMap.discretePdf(pixel);
        }
    }


protected:
    Independent() { }

private:
    pcg32 m_random;

    size_t currentPass;
    size_t weightedPassTotalSamples;

    Eigen::MatrixXf samplesPerPixel;
    Eigen::MatrixXf prevLuminance;

    LightProbeMipMap mipMap;
};

NORI_REGISTER_CLASS(Independent, "independent");
NORI_NAMESPACE_END

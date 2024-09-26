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
 * Stratified sampling - returns stratified uniformly distributed
 * random numbers on <tt>[0, 1)x[0, 1)</tt>.
 */
class Stratified : public Sampler {
public:
    Stratified(const PropertyList &propList) 
    {
        m_sampleCount = (size_t) propList.getInteger("FirstPassSpp", 1);
        _weightedPassesCount = (size_t) propList.getInteger("WeightedPassesCount", 0);
        weightedPassTotalSamples = (size_t) propList.getInteger("WeightedPassTotalSamples", 0);
        nStrats = (size_t) propList.getInteger("stratsX", 2);        
    }

    void _initialize(const Bitmap *initialImage) override
    {
        m_lastStrat = MatrixXf(imgX, imgY);
        m_lastStrat.setZero();

        samplesPerPixel = Eigen::MatrixXf::Zero(imgX, imgY);
        prevLuminance = Math::luminanceMatrix(initialImage);
        currentPass = 0;
        setConstantSpp(m_sampleCount);
    }

    virtual ~Stratified() { }

    std::unique_ptr<Sampler> clone() const override
    {
        std::unique_ptr<Stratified> cloned(new Stratified());

        cloned->m_sampleCount = m_sampleCount;
        cloned->imgX = imgX;
        cloned->imgY = imgY;
        cloned->_weightedPassesCount = _weightedPassesCount;

        cloned->m_random = m_random;
        cloned->samplesPerPixel = samplesPerPixel;
        cloned->weightedPassTotalSamples = weightedPassTotalSamples;

        cloned->nStrats = nStrats;
        cloned->m_lastStrat = m_lastStrat;

        return cloned;
    }

    void prepare(const ImageBlock &block) override
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

    Point2f generatePixelOffset(Point2i pixel) override
    {
        int lastStrat = m_lastStrat(pixel.x(), pixel.y());

        int xStrat = lastStrat % nStrats;
        int yStrat = lastStrat / nStrats;

        float x = (float(xStrat) + m_random.nextFloat()) / nStrats;
        float y = (float(yStrat) + m_random.nextFloat()) / nStrats;

        m_lastStrat(pixel.x(), pixel.y()) = (lastStrat + 1) % (nStrats * nStrats);

        return Point2f(x, y);
    }

    void generate() override { /* No-op for this sampler */ } 
    void advance()  override { /* No-op for this sampler */ }

    float next1D() override {
        return m_random.nextFloat();
    }
    
    Point2f next2D () override {
        return Point2f(
            m_random.nextFloat(),
            m_random.nextFloat()
        );
    }

    std::string toString() const override {
        return tfm::format("Stratified[sampleCount=%i]", m_sampleCount);
    }


    void matrixToFile(const Eigen::MatrixXf &m, const std::string &filename)
    {
        std::ofstream file(filename, std::ios::out);

        file << m << std::endl;

        file.close();
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

    void next_pass(const Bitmap *renderedImage)  override
    { 
        currentPass++; 
        samplesPerPixel.setZero();

        // Convert RGB to luminance
        auto image = Math::luminanceMatrix(renderedImage);
        initializeWeightedSamples(image);

        prevLuminance = image;
    }

    /// Return the number of configured pixel samples
    size_t getSampleCount(Point2i pixel) const override { return samplesPerPixel(pixel.x(), pixel.y()); }
    size_t getSampleCount() const override { return m_sampleCount; }

    float getSamplePdf(Point2i pixel) const override
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
    Stratified() 
    { 
        m_sampleCount = 1;
        _weightedPassesCount = 0;
        weightedPassTotalSamples = 0;
        nStrats = 4;
    }

private:
    pcg32 m_random;
    int nStrats;

    MatrixXf m_lastStrat;
    size_t currentPass;
    size_t weightedPassTotalSamples;

    Eigen::MatrixXf samplesPerPixel;
    Eigen::MatrixXf prevLuminance;

    LightProbeMipMap mipMap;
};

NORI_REGISTER_CLASS(Stratified, "stratified");
NORI_NAMESPACE_END

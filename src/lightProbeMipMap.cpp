
#include <nori/lightProbeMipMap.h>
#include <nori/math.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/bbox.cu>

NORI_NAMESPACE_BEGIN



void printMatrix(const Eigen::MatrixXf &m)
{
    for (int i = 0; i < m.rows(); i++)
    {
        for (int j = 0; j < m.cols(); j++)
        {
            std::cout << m(i, j) << " ";
        }
        std::cout << std::endl;
    }
}


// Based on https://pharr.org/matt/blog/2019/06/05/visualizing-env-light-warpings
LightProbeMipMap::LightProbeMipMap (const Eigen::MatrixXf &image, const BoundingBox2f &_domain) : domain(_domain)
{
    int sizeX = image.rows();
    int sizeY = image.cols();

    // Initialize first level
    pyramid.push_back(image);    

    // Compute level sizes
    sizeX = Math::roundUpPow2(sizeX);
    sizeY = Math::roundUpPow2(sizeY);

    while (sizeX > 2 || sizeY > 2) 
    {
        // Reduce level size by half
        sizeX = std::max(1, sizeX / 2); 
        sizeY = std::max(1, sizeY / 2);
        size_t prev = pyramid.size() - 1;

        // Create next level
        pyramid.push_back(Eigen::MatrixXf::Zero(sizeX, sizeY));
        
        // Reduce previous level into this one. 
        //  No needed to average, since we only care about relative values
        for (int y = 0; y < sizeY; ++y)
        {
            for (int x = 0; x < sizeX; ++x)
            {
                pyramid.back()(x, y) = l_lookup(prev, 2*x, 2*y) 
                                        + l_lookup(prev, 2*x + 1, 2*y) 
                                        + l_lookup(prev, 2*x, 2*y + 1) 
                                        + l_lookup(prev, 2*x + 1, 2*y + 1);
            }
        }
    }

    std::reverse(pyramid.begin(), pyramid.end());
}



Point2f LightProbeMipMap::staticSample(Point2f sample)
{
    float g_densities[2][3] = { {40, 10, 10}, 
                                {10, 20, 10} };

    // Build engen matrix
    Eigen::MatrixXf m(2, 3);
    m << g_densities[0][0], g_densities[0][1], g_densities[0][2],
         g_densities[1][0], g_densities[1][1], g_densities[1][2];

    static LightProbeMipMap *map = nullptr;

    if (!map) {
        map = loadMatrixFromFile("samplesPerPixel2.txt");
    }
    
    // Sample
    float pdf;
    return map->sampleContinuous(sample, pdf);
}

float LightProbeMipMap::staticSamplePdf(Point2f p)
{
    float g_densities[2][3] = { {60, 20}, 
                                {5, 15} };
                            
    // Build engen matrix
    Eigen::MatrixXf m(2, 3);
    m << g_densities[0][0], g_densities[0][1], g_densities[0][2],
         g_densities[1][0], g_densities[1][1], g_densities[1][2];

    LightProbeMipMap map(m, BoundingBox2f(Point2f(0, 0), Point2f(1, 1)));

    return map.continuousPdf(p);
}



// Get the value of the MipMap at a given level, local coordinates
float LightProbeMipMap::l_lookup (const Point2f &st, float level) const
{
    if (st.x() >= pyramid[level].rows() || st.y() >= pyramid[level].cols())
        return 0;

    return pyramid[level](st.x(), st.y());
}

float LightProbeMipMap::l_lookup (float level, int x, int y) const
{
    if (x >= pyramid[level].rows() || y >= pyramid[level].cols())
        return 0;

    return pyramid[level](x, y);
}


// Returns a random point in the MipMap sampled hierarchicaly by luminance
Point2i LightProbeMipMap::sampleDiscrete (Point2f sample, float &pdf) const
{
    pdf = 1;
    Point2i p(0, 0);

    // For each level
    for (unsigned int i = 0; i < pyramid.size(); ++i) 
    {
        // Get quadrant
        if (i > 0 && pyramid[i].rows() > pyramid[i-1].rows()) 
            p.x() *= 2;
        
        if (i > 0 && pyramid[i].cols() > pyramid[i-1].cols())
            p.y() *= 2;

        // Get quadrants sum to warp distribution
        float wx[2];
        wx[0] = l_lookup(Point2i(p.x(), p.y()), i) + l_lookup(Point2i(p.x(), p.y() + 1), i);
        wx[1] = l_lookup(Point2i(p.x() + 1, p.y()), i) + l_lookup(Point2i(p.x() + 1, p.y() + 1), i);
    
        // Warp X axis
        float samplePdf;
        p.x() += Warp::sampleDiscrete(wx, 2, sample.x(), samplePdf, sample.x());
        pdf *= samplePdf;

        // Warp Y axis
        float wy[2];
        wy[0] = l_lookup(Point2i(p.x(), p.y()), i);
        wy[1] = l_lookup(Point2i(p.x(), p.y() + 1), i);
        
        p.y() += Warp::sampleDiscrete(wy, 2, sample.y(), samplePdf, sample.y());
        pdf *= samplePdf;
    }

    return p;
}

float LightProbeMipMap::discretePdf(Point2i p) const 
{
    float pdf = 1;

    for (int i = pyramid.size() - 1; i >= 0; --i) 
    {
        // Make coordinates even, rounding down
        Point2i pe(p.x() & ~1, p.y() & ~1);

        float value = l_lookup(Point2i(p.x(), p.y()), i);

        if (value == 0) 
            return 0;

        pdf *= value / (l_lookup(Point2i(pe.x(), pe.y()), i) +
                    l_lookup(Point2i(pe.x() + 1, pe.y()), i) +
                    l_lookup(Point2i(pe.x(), pe.y() + 1), i) +
                    l_lookup(Point2i(pe.x() + 1, pe.y() + 1), i));

        if (i > 0 && pyramid[i-1].rows() > pyramid[i].rows()) 
            p.x() /= 2;
        if (i > 0 && pyramid[i-1].cols() > pyramid[i].cols()) 
            p.y() /= 2;
    }

    return pdf;
}


Point2f LightProbeMipMap::sampleContinuous(Point2f sample, float &pdf) const 
{
    pdf = 1;
    Point2i p(0, 0);

    for (size_t i = 0; i < pyramid.size(); ++i) 
    {
        if (i > 0 && pyramid[i].rows() > pyramid[i-1].rows()) 
            p.x() *= 2;
        if (i > 0 && pyramid[i].cols() > pyramid[i-1].cols())
            p.y() *= 2;

        float wx[2];
        wx[0] = l_lookup(Point2i(p.x(), p.y()), i) + l_lookup(Point2i(p.x(), p.y() + 1), i);
        wx[1] = l_lookup(Point2i(p.x() + 1, p.y()), i) + l_lookup(Point2i(p.x() + 1, p.y() + 1), i);
        
        float sampPdf;
        p.x() += Warp::sampleDiscrete(wx, 2, sample.x(), sampPdf, sample.x());
        pdf *= sampPdf;

        float wy[2];
        wy[0] = l_lookup(Point2i(p.x(), p.y()), i);
        wy[1] = l_lookup(Point2i(p.x(), p.y() + 1), i);
        p.y() += Warp::sampleDiscrete(wy, 2, sample.y(), sampPdf, sample.y());
        pdf *= sampPdf;
    }

    pdf *= (pyramid.back().rows() * pyramid.back().cols()) / domain.getSurfaceArea();

    float x = (p.x() + sample.x()) / pyramid.back().rows();
    float y = (p.y() + sample.y()) / pyramid.back().cols();

    return domain.lerp(Point2f(x, y));
}

float LightProbeMipMap::continuousPdf(Point2f p) const 
{
    Vector2f o = domain.getUnitCubeCoordinates(p);

    return discretePdf(Point2i(o.x() * pyramid.back().rows(), o.y() * pyramid.back().cols())) *
        (pyramid.back().rows() * pyramid.back().cols()) / domain.getSurfaceArea();
}



NORI_NAMESPACE_END


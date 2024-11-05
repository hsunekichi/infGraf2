#pragma once

NORI_NAMESPACE_BEGIN

struct PathState
{
    bool previous_diffuse = false; bool previous_sss = false;
    float bsdfPdf = 0.0f; Point3f prevP = Point3f(0.0f);    
    Intersection intersection;
    Point2f pixel;

    Ray3f ray; bool rayIntersected;
    int depth = 0;

    Color3f radiance = Color3f(0.0f);
    Color3f scatteringFactor = Color3f(1.0f);
};

NORI_NAMESPACE_END
#include "extra.h"
#include "bvh.h"
#include "common.h"
#include "glm/gtx/string_cast.hpp"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <algorithm>
#include <cassert>
#include <framework/trackball.h>
#include <texture.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <iostream>
#include <limits>
#include <iostream>
#include <limits>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    // Slides Ray Tracing : https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3475537/View
    if (!features.extra.enableMotionBlur) {
        return;
    }

    // we use control points for the Bezier cubic curve
    std::vector<glm::vec3> controlPoints(4);

    controlPoints[0] = glm::vec3(-0.5, 1, 0);
    controlPoints[1] = glm::vec3(0, -1, -2);
    controlPoints[2] = glm::vec3(0.5, 0.5, 1);
    controlPoints[3] = glm::vec3(-1, 0, -1);

    const int numSamples = 30;
    std::vector<glm::vec3> samplePixels(screen.resolution().x * screen.resolution().y);
    samplePixels = screen.pixels();

    for (int s = 0; s < numSamples; ++s) {
        Scene local_scene = scene;

        for (Mesh& m : local_scene.meshes) {
            for (Vertex& v : m.vertices) {

                // t = how far in 'time' are we, based on the number of the current sample
                float t = s * 1.0f / numSamples;

                // we use a cubic Bezier curve for the movement pattern, the formula is on the website:
                // (cubic curves explicit formula section)
                // https://en.wikipedia.org/wiki/B%C3%A9zier_curve

                glm::vec3 point_on_curve = (float)glm::pow((1 - t), 3) * controlPoints[0]
                    + 3 * (float)glm::pow((1 - t), 2) * t * controlPoints[1]
                    + 3 * (1 - t) * (float)glm::pow(t, 2) * controlPoints[2]
                    + (float)glm::pow(t, 3) * controlPoints[3];

                // the division by 15 controls the displacement of the vertices (makes it visually less rough)
                v.position += (point_on_curve / 15.0f);
            }

            // new bvh data structure for the current scene
            BVH bvh(local_scene, features);

            // the next 3 lines of code are borrowed from the file render.cpp, line 27
            #ifdef NDEBUG // Enable multi threading in Release mode
            #pragma omp parallel for schedule(guided)
            #endif

            // we generate rays and notice their contribution on the pixels' color
            for (int i = 0; i < screen.resolution().y; ++i) {
                for (int j = 0; j < screen.resolution().x; ++j) {

                    RenderState state = { local_scene, features, bvh };
                    glm::vec2 pixel = glm::vec2(j, i);
                    auto rays = generatePixelRays(state, camera, pixel, screen.resolution());
                    samplePixels[screen.indexAt(j, i)] += renderRays(state, rays);
                }
            }
        }
    }

    // display changes
    for (int i = 0; i < screen.resolution().y; i++) {
        for (int j = 0; j < screen.resolution().x; j++)
            screen.setPixel(j, i, samplePixels[screen.indexAt(j, i)] * (1.0f / numSamples));
        
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float factorial(int n) {        // calculates factorial

    float result = 1;
    for (int i = 1; i <= n; ++i)
        result *= i;

    return result;
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // if filter size is odd, it's good
    int filter_size;
    if (features.filterSize % 2 == 1)
        filter_size = features.filterSize;
    else
        filter_size = features.filterSize - 1;


    float scale = features.scaleSize / 10.0f;
    std::vector<float> filter(filter_size + 1);

    for (int i = 0; i <= filter_size; ++i) {    // calculates the filter 'box' and normalizez (/2^n)
        filter[i] = factorial(filter_size) / (factorial(i) * factorial(filter_size - i)) / glm::pow(2, filter_size);
    }

    Screen temp = image;
    glm::vec3 black = glm::vec3(0, 0, 0);

    // make a new pictures that only keeps the light pixels and the rest turn black
    for (int i = 0; i < image.resolution().y; ++i) {
        for (int j = 0; j < image.resolution().x; ++j) {
            glm::vec3 pixel = image.pixels()[image.indexAt(i, j)];
            if (pixel.x > 0.5f && pixel.y > 0.5f && pixel.z > 0.5f)
                temp.setPixel(i, j, pixel);
            else
                temp.setPixel(i, j, black);
        }
    }

    Screen temp_horizontal = temp;
    // filter horizontally
    for (int i = 0; i < image.resolution().y; ++i) {
        for (int j = 0; j < image.resolution().x; ++j) {
            if (temp.pixels()[image.indexAt(i, j)] != black) {

                glm::vec3 value(0, 0, 0);
                for (int k = 0; k < filter_size; ++k) {

                    int index = k + j - filter_size / 2;
                    if (index >= 0 && index < image.resolution().y)     // in between boundaries
                        value = value + (filter[k] * temp.pixels()[image.indexAt(i, index)]);
                }

                temp_horizontal.setPixel(i, j, value);
            }
        }
    }

    Screen temp_vertical = temp_horizontal;
    // filter vertically
    for (int i = 0; i < image.resolution().x; ++i) {    // do the same as before, but the matrix is transposed
        for (int j = 0; j < image.resolution().y; ++j) {
            if (temp.pixels()[image.indexAt(i, j)] != black) {

                glm::vec3 value(0, 0, 0);
                for (int k = 0; k < filter_size; ++k) {

                    int index = k + j - filter_size / 2;
                    if (index >= 0 && index < image.resolution().x)
                        value = value + (filter[k] * temp.pixels()[image.indexAt(i, index)]);
                }

                temp_vertical.setPixel(i, j, value);
            }
        }
    }

    // add changes to the image
    for (int i = 0; i < image.resolution().y; ++i) {
        for (int j = 0; j < image.resolution().x; ++j) {
            image.pixels()[image.indexAt(i, j)] += temp_vertical.pixels()[image.indexAt(i, j)] * scale;
        }
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// I took the following method from assignment 1c. I modified it a bit so that it fits my code here.
// This method return a basis for the orthogonal plane (2 vectors)
std::array<glm::vec3, 2> plane_basis(glm::vec3 n)
{
    std::array<glm::vec3, 2> result;

    glm::vec3 v1 {};
    glm::vec3 v2 {};

    // first vector //
    if (n.x == 0 && n.y == 0) {
        v1.x = v1.y = 1.0f;
        v1.z = 0;
    } else if (n.x == 0) {
        v1.x = v1.z = 1.0f;
        v1.y = (1.0f * -n.z) / n.y;
    } else {
        v1.y = v1.z = 1.0f;
        v1.x = -(n.y + n.z) * 1.0f / n.x;
    }

    // second vector //
    // cross product //
    v2.x = n.y * v1.z - n.z * v1.y;
    v2.y = n.z * v1.x - n.x * v1.z;
    v2.z = n.x * v1.y - n.y * v1.x;

    result[0] = glm::normalize(v1);
    result[1] = glm::normalize(v2);

    return result;
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // https://blog.yiningkarlli.com/2012/12/blurred-glossy-reflections.html
    // Used this link to better grasp the concept.

     auto numSamples = state.features.extra.numGlossySamples;
     Ray specular = generateReflectionRay(ray, hitInfo);

     for (int i = 0; i < numSamples; ++i) {
        // r: 0 - 1
        float r = state.sampler.next_1d();
        // angle: 0 - 2pi
        float angle = state.sampler.next_1d() * 2.0f * glm::pi<float>();

        glm::vec3 origin = ray.origin + ray.direction * (ray.t - 0.001f);

        // we use the method from assignment 1c, provided above
        std::array<glm::vec3, 2> plane = plane_basis(hitInfo.normal);
        glm::vec3 u = plane[0];
        glm::vec3 v = plane[1];

        // sqrt(r) helps have samples on the disk
        float radius = 1.0f * hitInfo.material.shininess / 64 * sqrt(r);

        // https://stats.stackexchange.com/questions/481543/generating-random-points-uniformly-on-a-disk/481544#481544
        // I used this website for a uniform distribution on a disk.
        float x = radius * cos(angle); 
        float y = radius * sin(angle);

        glm::vec3 glossyReflection = glm::normalize(x * u + y * v + specular.direction);

        Ray glossyRay(origin, glossyReflection);
        glm::vec3 result = renderRay(state, glossyRay, rayDepth + 1);
        hitColor += result * hitInfo.material.ks;
     }
     hitColor = hitColor * (1.0f / numSamples);
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
/*
Sources: 
    https://mechref.engr.illinois.edu/dyn/rvs.html
    Computer Graphics 4th Edition, chapters 4, 11, and 13
    Slides Ray Tracing : https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3475537/View
*/
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if(!state.features.extra.enableEnvironmentMap || !state.scene.environmentMap)
    {
        return glm::vec3(0.f); // Return black if environment mapping is not enabled or not present
    }

    // Normalize the ray direction for the spherical mapping
    glm::vec3 dir = glm::normalize(ray.direction);

    // Convert direction to UV coordinates
    float u = (M_PI + atan2(dir.z, dir.x)) / (2 * M_PI);
    float v = acos(glm::clamp(dir.y, -1.0f, 1.0f)) / M_PI;

    // Sample the environment map image
    return sampleTexture(*state.scene.environmentMap, u, v);
    
}
/*
Sources:
    https://mechref.engr.illinois.edu/dyn/rvs.html
    Computer Graphics 4th Edition, chapters 4, 11, and 13
    Slides Ray Tracing : https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3475537/View
*/
glm::vec3 sampleTexture(const Image& image, float u, float v)
{
    // Clamp the coordinates
    u = std::clamp(u, 0.0f, 1.0f);
    v = std::clamp(v, 0.0f, 1.0f);

    // Map the UV coordinate to the image's dimensions
    unsigned int x = static_cast<unsigned int>(u * image.width);
    unsigned int y = static_cast<unsigned int>(v * image.height);

    // Clamp to the image's bounds
    x = std::min(x, static_cast<unsigned int>(image.width - 1));
    y = std::min(y, static_cast<unsigned int>(image.height - 1));

    // Fetch the color from the image's pixels
    return image.pixels[y * image.width + x];
}

struct Bin {
    // Initialize as "negative" AABB
    AxisAlignedBox aabb {
        .lower = glm::vec3(std::numeric_limits<float>::max()),
        .upper = glm::vec3(std::numeric_limits<float>::lowest())
    };
    uint32_t numTriangles { 0 };

    constexpr void grow(BVHInterface::Primitive p)
    {
        numTriangles++;
        aabb.lower = glm::min(glm::min(aabb.lower, p.v0.position), glm::min(p.v1.position, p.v2.position));
        aabb.upper = glm::max(glm::max(aabb.upper, p.v0.position), glm::max(p.v1.position, p.v2.position));
    }
};
// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
//
// Formula:
// $$
// c(A, B) = t_{trav} + p_A \sum_{i=1}^{N_A}(t_{isect}(a_i)) \sum_{i=1}^{N_B}(t_{isect}(b_i))
// $$
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;
    const int NUM_BINS = 16;

    if (primitives.size() <= 1)
        return 0;

    float min = std::numeric_limits<float>::infinity();
    float max = -std::numeric_limits<float>::infinity();

    for (auto prim : primitives) {
        min = glm::min(min, computePrimitiveCentroid(prim, axis));
        max = glm::max(max, computePrimitiveCentroid(prim, axis));
    }

    // std::cout << "min: " << min << "; max: " << max << std::endl;

    // formula to determine primitive bin taken from:
    // https://www.sci.utah.edu/~wald/Publications/2007/ParallelBVHBuild/fastbuild.pdf
    // Chapter 3.3, p. 4
    const float factor = (float)NUM_BINS / (max - min);
    Bin bins[NUM_BINS];

    int maxBinID = 0;
    for (auto prim : primitives) {
        // 1. determine primitive bin
        const float centroidOnAxis = computePrimitiveCentroid(prim, axis);
        const int binID = glm::clamp((int)(factor * (centroidOnAxis - min)), 0, NUM_BINS - 1);
        // 2. expand aabb of bin to account for it.
        bins[binID].grow(prim);
    }
    // std::cout << "max bin id " << maxBinID << "; bin sizes: ";
    // for (int i = 0; i < NUM_BINS; i++) {
    //     std::cout << bins[i].numTriangles << "; ";
    // }
    // std::cout << std::endl;
    // We need to accumulate the left and right aabb and bounds so we can then process the splits
    // We do this by sweeping across the array from the left and right, summing up the triangle counts and joining the AABBs
    // On the right-sweep we immediately calculate the SAH.
    // (As described in previously mentioned paper)
    float leftSurfaceArea[NUM_BINS - 1];
    uint32_t leftTriangles[NUM_BINS - 1];

    AxisAlignedBox cumAABB = AxisAlignedBox::negative();
    uint32_t cumTriangles = 0;
    for (int i = 0; i < NUM_BINS - 1; i++) {
        cumTriangles += bins[i].numTriangles;
        cumAABB.merge(bins[i].aabb);
        leftSurfaceArea[i] = glm::min(cumAABB.surfaceArea(), std::numeric_limits<float>::max());
        leftTriangles[i] = cumTriangles;
    }

    cumAABB = AxisAlignedBox::negative();
    cumTriangles = 0;

    float bestSAH = std::numeric_limits<float>::infinity();
    float bestSAHSplit = 0;
    int winningBin;
    for (int i = NUM_BINS - 1; i > 0; i--) {
        cumTriangles += bins[i].numTriangles;
        cumAABB.merge(bins[i].aabb);
        float rightSurfaceArea = glm::min(cumAABB.surfaceArea(), std::numeric_limits<float>::max());

        const float SAH = rightSurfaceArea * (float)cumTriangles + leftSurfaceArea[i - 1] * (float)leftTriangles[i - 1];
        // std::cout << "bin " << i << ": SAH=" << SAH << "; rightSurfaceArea=" << rightSurfaceArea << "; rightTriangles=" << cumTriangles << "; leftSurfaceArea=" << leftSurfaceArea[i - 1] << "; leftTriangles=" << leftTriangles[i - 1] << std::endl;

        if (SAH < bestSAH) {
            bestSAH = SAH;
            winningBin = i;
            const float splitFactor = (float)i * (1.f / (float)NUM_BINS);
            // LERP by splitFactor
            bestSAHSplit = max * splitFactor + min * (1 - splitFactor);
        }
    }
    // std::cout << "winning bin: " << winningBin << "; best SAH: " << bestSAH << "; axis " << axis << std::endl;

    // We can split the primitives by just partitioning them based on if they are smaller than the split index.
    auto split = std::partition(primitives.begin(), primitives.end(), [axis, bestSAHSplit](Primitive p) {
        return computePrimitiveCentroid(p, axis) <= bestSAHSplit;
    });
    if (split == primitives.begin() || split == primitives.end()) {
        // we ran into a degenerate case (triangles on the same plane).
        return splitPrimitivesByMedian(aabb, axis, primitives);
    } else {
        return (size_t)std::distance(primitives.begin(), split);
    }
}

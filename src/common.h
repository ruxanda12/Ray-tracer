#pragma once

#include <framework/disable_all_warnings.h>
#include <limits>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/mesh.h>

enum class DrawMode {
    Filled,
    Wireframe
};

enum class ShadingModel {
    Lambertian = 0,
    Phong = 1,
    BlinnPhong = 2,
    LinearGradient = 3,
};

struct HitInfo {
    glm::vec3 normal;
    glm::vec3 barycentricCoord;
    glm::vec2 texCoord;
    Material material;
};

struct Plane {
    float D = 0.0f;
    glm::vec3 normal { 0.0f, 1.0f, 0.0f };
};

struct AxisAlignedBox {
    glm::vec3 lower { 0.0f };
    glm::vec3 upper { 1.0f };

    static constexpr AxisAlignedBox negative()
    {
        return AxisAlignedBox {
            .lower = glm::vec3(std::numeric_limits<float>::max()),
            .upper = glm::vec3(std::numeric_limits<float>::lowest())
        };
    }

    constexpr void merge(AxisAlignedBox other)
    {
        lower = glm::min(lower, other.lower);
        upper = glm::max(upper, other.upper);
    }

    constexpr float surfaceArea() const
    {
        const glm::vec3 dims = upper - lower;
        return dims.x * dims.y + dims.y * dims.z + dims.x * dims.z;
    }
};

struct Sphere {
    glm::vec3 center { 0.0f };
    float radius = 1.0f;
    Material material;
};

struct PointLight {
    glm::vec3 position;
    glm::vec3 color;
};

struct SegmentLight {
    glm::vec3 endpoint0, endpoint1; // Positions of endpoints
    glm::vec3 color0, color1; // Color of endpoints
};

struct ParallelogramLight {
    // A parallelogram light (see figure 3.14 of chapter 13.4.2 of Fundamentals of CG 4th Edition)
    glm::vec3 v0; // v0
    glm::vec3 edge01, edge02; // edges from v0 to v1, and from v0 to v2
    glm::vec3 color0, color1, color2, color3;
};

struct ExtraFeatures {
    bool enableBvhSahBinning = false;
    bool enableBloomEffect = false;
    bool enableDepthOfField = false;
    bool enableEnvironmentMap = false;
    bool enableGlossyReflection = false;
    bool enableMipmapTextureFiltering = false;
    bool enableMotionBlur = false;

    // Parameters for glossy reflection
    uint32_t numGlossySamples = 1;
};

struct Features {
    // Feature toggles
    bool enableShading = false;
    bool enableReflections = false;
    bool enableShadows = false;
    bool enableNormalInterp = false;
    bool enableTextureMapping = false;
    bool enableAccelStructure = false;
    bool enableBilinearTextureFiltering = false;
    bool enableTransparency = false;
    bool enableJitteredSampling = false;

    // Feature-specific settings
    ShadingModel shadingModel = ShadingModel::Lambertian;
    uint32_t numPixelSamples = 1;
    uint32_t numShadowSamples = 4;
    uint32_t scaleSize = 5;
    uint32_t filterSize = 5;

    // Extras-specific settings
    ExtraFeatures extra = {};
};

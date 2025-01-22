#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()

float eps = 1e-4f;


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    position = (1 - sample) * light.endpoint0 + sample * light.endpoint1; 
    color = (1 - sample) * light.color0 + sample * light.color1;

    Ray debugRaySegment = { .origin = position, .direction = glm::normalize(glm::cross(light.endpoint0, light.endpoint1)), .t = 0.1f };
    drawRay(debugRaySegment, color);
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    position = light.v0 + sample.x * light.edge01 + sample.y * light.edge02;
    glm::vec3 color_edge01 = (1 - sample.x) * light.color0 + sample.x * light.color1;
    glm::vec3 color_edge02 = (1 - sample.x) * light.color2 + sample.x * light.color3;
    color = (1 - sample.y) * color_edge01 + sample.y * color_edge02;

    glm::vec3 normal = glm::normalize(glm::cross(light.edge01, light.edge02));
    Ray debugRay = { .origin = position, .direction = normal, .t = 0.1 };
    drawRay(debugRay, color);
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
/*
Sources:
    Computer Graphics 4th Edition : chapter 4
    Slides Ray Tracing : https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3475537/View
*/
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{

    //TODO: Assume that method needs to be extended to allow for light sources with volume(area lights) and reflected/transp surfaces
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer
        // TODO: implement this function; currently, the light simply passes through
        Ray shadowRay;
        shadowRay.origin = ray.origin + ray.t * ray.direction + hitInfo.normal * eps;
        shadowRay.direction = glm::normalize(lightPosition - shadowRay.origin);
        shadowRay.t = glm::length(lightPosition - shadowRay.origin);

        HitInfo shadowHit;

        bool result = false;

        if (!state.bvh.intersect(state, shadowRay, shadowHit)) {
            result = true;
        }

        glm::vec3 debugColor = (result) ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
        drawRay(shadowRay, debugColor);

        return result;

    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
/*
Sources:
    Computer Graphics 4th Edition : chapters 4 and 13
    Slides Ray Tracing : https://brightspace.tudelft.nl/d2l/le/content/595314/viewContent/3475537/View
*/
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: implement this function; currently, the light simply passes through
    if (!state.features.enableShadows) {
        return lightColor;
    } else {
        Ray shadowRay;
        shadowRay.origin = ray.origin + ray.t * ray.direction + hitInfo.normal * eps;
        shadowRay.direction = glm::normalize(lightPosition - shadowRay.origin);
        shadowRay.t = glm::length(lightPosition - shadowRay.origin);

        HitInfo shadowHit;

        if (state.bvh.intersect(state, shadowRay, shadowHit)) {
            if (shadowHit.material.transparency < 1) {
                drawRay(shadowRay, lightColor * shadowHit.material.kd * (1 - shadowHit.material.transparency));
                return lightColor * shadowHit.material.kd * (1 - shadowHit.material.transparency);
            } else {
                drawRay(shadowRay, glm::vec3(1.0f, 0.0f, 0.0f));
                return glm::vec3(0);
            }
        } else {
            drawRay(shadowRay, lightColor);
            return lightColor;
        }
    }
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: modify this function to incorporate visibility corerctly
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;

    glm::vec3 lightAct = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);

    if (lightAct != glm::vec3(0)) {
        return computeShading(state, v, l, light.color, hitInfo);
    }else {
        return glm::vec3(0);
    }
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 totalCont = glm::vec3(0.0f);

    for (int i = 0; i < numSamples; ++i) {
        glm::vec3 sampledPos, sampledCol;

        sampleSegmentLight(state.sampler.next_1d(), light, sampledPos, sampledCol);

        glm::vec3 lightContribution = visibilityOfLightSample(state, sampledPos, sampledCol, ray, hitInfo);

        if (lightContribution != glm::vec3(0)) {
            glm::vec3 l = glm::normalize(sampledPos - (ray.origin + ray.t * ray.direction));

            totalCont += computeShading(state, -ray.direction, l, sampledCol, hitInfo);
        
        }
    }

    return totalCont / static_cast<float>(numSamples);
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model

    glm::vec3 totalCont = glm::vec3(0.0f);

    for (int i = 0; i < numSamples; ++i) {
        glm::vec3 sampledPos, sampledCol;

        sampleParallelogramLight(state.sampler.next_2d(), light, sampledPos, sampledCol);

        glm::vec3 lightContribution = visibilityOfLightSample(state, sampledPos, sampledCol, ray, hitInfo);

        if (lightContribution != glm::vec3(0)) {
            glm::vec3 l = glm::normalize(sampledPos - (ray.origin + ray.t * ray.direction));

            totalCont += computeShading(state, -ray.direction, l, sampledCol, hitInfo);
        }
    }

    return totalCont / static_cast<float>(numSamples);
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}
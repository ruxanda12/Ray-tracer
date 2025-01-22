#pragma once

#include "common.h"
#include "fwd.h"
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>

/* Sampler functions */

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo);

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo);

// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and color on the segment light, and write these into the reference return values.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color);

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and color on the paralellogram light, and write these into the reference return values.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color);

/* Contribution sub-functions */

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& pointLight, const Ray& ray, const HitInfo& hitInfo);

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& segmentLight, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples);

// TODO: Standard feature
// Given a single paralellogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the paralellogram, taking `numSamples` samples.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& parallelogramLight, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples);

/* Provided functions */

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forwards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo);

// This function is provided as-is. You do not have to implement it.
// Given an incident ray and an intersection, accumulates the contribution of all lights in the scene
// with respect to this ray and intersection. The contributions of individual light types are computed
// in `computeContributionPointLight()`, computeContributionSegmentLight()`, and
// `computeContributionParallelogramLight()`, which you must implement yourself.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo);
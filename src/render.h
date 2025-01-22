#pragma once

#include "common.h"
#include "fwd.h"
#include "sampler.h"
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>

// The configurative state inside renderer; collects
// handles to e.g. the BVH and the scene, and holds
// a per-thread random sampler and other things you
// might want to access in many places.
struct RenderState {
    // Handles to global objects accessed throughout the renderer
    const Scene& scene; // The scene being rendered
    const Features& features; // The feature config that is active
    const BVHInterface& bvh; // The bvh generated over the current scene

    // Small per-thread objects kept alive throughout the renderer
    // You can add your own objects here ...
    Sampler sampler; // 1d/2d sampler on the range [0, 1)
};

/* Baseline render code; you do not have to implement the following methods */

// This function is provided as-is. You do not have to implement it.
// Given relevant objects (scene, bvh, camera, etc) and an output screen, multithreaded fills
// each of the pixels using one of the below `renderPixel*()` functions, dependent on scene
// configuration. By default, `renderPixelNaive()` is called.
void renderImage(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen);

// This function is provided as-is. You do not have to implement it.
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples for this pixel.
// This method forwards to `generatePixelRaysMultisampled` and `generatePixelRaysStratified` when necessary.
std::vector<Ray> generatePixelRays(RenderState &state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution);

/* Unfinished render code; you have to implement the following method */

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// uniformly throughout this pixel.
// For a description of the method's arguments, refer to 'render.cpp'
// This method is unit-tested, so do not change the function signature.
std::vector<Ray> generatePixelRaysMultisampled(RenderState &state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution);

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// using stratified/jittered sampling throughout this pixel.
// For a description of the method's arguments, refer to 'extra.cpp'
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
std::vector<Ray> generatePixelRaysStratified(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution);

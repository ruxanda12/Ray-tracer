#include "render.h"
#include "bvh_interface.h"
#include "draw.h"
#include "extra.h"
#include "light.h"
#include "recursive.h"
#include "sampler.h"
#include "screen.h"
#include "shading.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

// This function is provided as-is. You do not have to implement it.
// Given relevant objects (scene, bvh, camera, etc) and an output screen, multithreaded fills
// each of the pixels using one of the below `renderPixel*()` functions, dependent on scene
// configuration. By default, `renderPixelNaive()` is called.
void renderImage(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    // Either directly render the image, or pass through to extra.h methods
    if (features.extra.enableDepthOfField) {
        renderImageWithDepthOfField(scene, bvh, features, camera, screen);
    } else if (features.extra.enableMotionBlur) {
        renderImageWithMotionBlur(scene, bvh, features, camera, screen);
    } else {
#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
        for (int y = 0; y < screen.resolution().y; y++) {
            for (int x = 0; x != screen.resolution().x; x++) {
                // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
                // Note; we seed the sampler for consistenct behavior across frames
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };
                auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
                auto L = renderRays(state, rays);
                screen.setPixel(x, y, L);
            }
        }
    }

    // Pass through to extra.h for post processing
    if (features.extra.enableBloomEffect) {
        postprocessImageWithBloom(scene, features, camera, screen);
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples for this pixel.
// This method forwards to `generatePixelRaysMultisampled` and `generatePixelRaysStratified` when necessary.
std::vector<Ray> generatePixelRays(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    if (state.features.numPixelSamples > 1) {
        if (state.features.enableJitteredSampling) {
            return generatePixelRaysStratified(state, camera, pixel, screenResolution);
        } else {
            return generatePixelRaysMultisampled(state, camera, pixel, screenResolution);
        }
    } else {
        // Generate single camera ray placed at the pixel's center
        // Note: (-1, -1) at the bottom left of the screen,
        //       (+1, +1) at the top right of the screen.
        glm::vec2 position = (glm::vec2(pixel) + 0.5f) / glm::vec2(screenResolution) * 2.f - 1.f;
        return { camera.generateRay(position) };
    }
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// uniformly throughout this pixel.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is unit-tested, so do not change the function signature.
std::vector<Ray> generatePixelRaysMultisampled(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples camera rays uniformly distributed across the pixel. Use
    // Hint; use `state.sampler.next*d()` to generate random samples in [0, 1).
    auto numSamples = state.features.numPixelSamples;
    std::vector<Ray> rays;

    for (int i = 0; i < numSamples; ++i) {
        glm::vec2 screen_coord = glm::vec2(pixel) + state.sampler.next_2d();
        glm::vec2 position = screen_coord / glm::vec2(screenResolution);
        rays.push_back(camera.generateRay(position * 2.0f - 1.0f));
    }

    return rays;
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// using jittered sampling throughout this pixel. Given NxN cells across one pixel, each ray sample is randomly
// placed somewhere within a cell.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
std::vector<Ray> generatePixelRaysStratified(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples * numSamples camera rays as jittered samples across the pixel.
    // Hint; use `state.sampler.next*d()` to generate random samples in [0, 1).
    auto numSamples = static_cast<uint32_t>(std::round(std::sqrt(float(state.features.numPixelSamples))));
    std::vector<Ray> rays;

    for (int i = 0; i < numSamples; ++i) {
        for (int j = 0; j < numSamples; ++j) {
            glm::vec2 screen_coord = glm::vec2(pixel) + state.sampler.next_2d();
            glm::vec2 position = (screen_coord + glm::vec2(i, j)) / glm::vec2(screenResolution);
            rays.push_back(camera.generateRay(position * 2.0f - 1.0f));
        }
    }

    return rays;
}
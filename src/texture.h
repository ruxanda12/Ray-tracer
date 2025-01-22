#pragma once

#include "common.h"
#include "fwd.h"
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image &image, const glm::vec2 &texCoord);

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// For a description of the method's arguments, refer to 'light.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image &image, const glm::vec2 &texCoord);
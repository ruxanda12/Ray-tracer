#pragma once
#include "fwd.h"
#include "common.h"
#include <framework/ray.h>

/* Baseline render code; you do not have to implement the following methods */

// This function is provided as-is. You do not have to implement it.
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState &state, const HitInfo &hitInfo);

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState &state, const glm::vec3 &cameraDirection, const glm::vec3 &lightDirection, const glm::vec3 &lightColor, const HitInfo &hitInfo);

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState &state, const glm::vec3 &cameraDirection, const glm::vec3 &lightDirection, const glm::vec3 &lightColor, const HitInfo &hitInfo);

/* Unfinished render code; you have to implement the following methods */

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// For a description of the method's arguments, refer to 'shading.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState &state, const glm::vec3 &cameraDirection, const glm::vec3 &lightDirection, const glm::vec3 &lightColor, const HitInfo &hitInfo);

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// For a description of the method's arguments, refer to 'shading.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState &state, const glm::vec3 &cameraDirection, const glm::vec3 &lightDirection, const glm::vec3 &lightColor, const HitInfo &hitInfo);

struct LinearGradient {
  struct Component {
    float t;
    glm::vec3 color;
  };

  std::vector<Component> components;

  // TODO: Standard feature
  // Given a number ti between [-1, 1], sample from the gradient's components and return the
  // linearly interpolated color, for which ti lies in the interval between the t-values of two
  // components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
  // the nearest component must be sampled.
  // For a description of the method's arguments, refer to 'shading.cpp'
  // This method is unit-tested, so do not change the function signature.
  glm::vec3 sample(float ti) const;
};

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
// For a description of the method's arguments, refer to 'shading.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState &state, const glm::vec3 &cameraDirection, const glm::vec3 &lightDirection, const glm::vec3 &lightColor, const HitInfo &hitInfo, const LinearGradient &gradient);

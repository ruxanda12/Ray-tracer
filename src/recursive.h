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

/* Baseline render code; you do not have to implement the following methods */

// This function is provided as-is. You do not have to implement it.
// Helper; given a set of rays, render all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth = 0);

// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - `computeLightContribution()` and its submethods
// - `renderRaySpecularComponent()`, `renderRayTransparentComponent()`, `renderRayGlossyComponent()`
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth = 0);

/* Unfinished render code; you have to implement the following methods */

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// For a description of the method's arguments, refer to 'recursive.cpp'
// This method is unit-tested, so do not change the function signature.
Ray generateReflectionRay(Ray ray, HitInfo hitInfo);

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// For a description of the method's arguments, refer to 'recursive.cpp'
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo);

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, compute the contribution
// of a mirror ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// For a description of the method's arguments, refer to 'recursive.cpp'
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth);

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, compute the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// For a description of the method's arguments, refer to 'recursive.cpp'
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth);
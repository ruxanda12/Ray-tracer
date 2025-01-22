#pragma once

#include "fwd.h"
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// For a description of the method's arguments, refer to 'interpolate.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p);

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// For a description of the method's arguments, refer to 'interpolate.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc);

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// For a description of the method's arguments, refer to 'interpolate.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc);
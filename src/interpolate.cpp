#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.
    //return glm::vec3(0.0);
    glm::vec3 v0p = p - v0;
    glm::vec3 v01 = v1 - v0;
    glm::vec3 v02 = v2 - v0;

    float dot_01 = glm::dot(v0p, v01);
    float dot_02 = glm::dot(v0p, v02);
    float dot_11 = glm::dot(v01, v01);
    float dot_12 = glm::dot(v01, v02);
    float dot_22 = glm::dot(v02, v02);

    float denominator = dot_11 * dot_22 - dot_12 * dot_12;

    float beta = (dot_22 * dot_01 - dot_12 * dot_02) / denominator;
    float kapsalon = (dot_11 * dot_02 - dot_12 * dot_01) / denominator;
    float alpha = 1.0f - beta - kapsalon;

    glm::vec3 barycentric_coord = glm::vec3(alpha, beta, kapsalon);
    return barycentric_coord;

}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    // TODO: implement this function.

    return glm::vec3(n0.x * bc.x + n1.x * bc.y + n2.x * bc.z,
        n0.y * bc.x + n1.y * bc.y + n2.y * bc.z,
        n0.z * bc.x + n1.z * bc.y + n2.z * bc.z);
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
// TODO: implement this function.

    return glm::vec2(t0.x * bc.x + t1.x * bc.y + t2.x * bc.z,
        t0.y * bc.x + t1.y * bc.y + t2.y * bc.z);
}

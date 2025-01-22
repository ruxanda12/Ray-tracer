#include "recursive.h"
#include "draw.h"
#include "bvh_interface.h"
#include "intersect.h"
#include "extra.h"
#include "light.h"

// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth)
{
    glm::vec3 L { 0.f };
    for (const auto& ray : rays) {
        L += renderRay(state, ray, rayDepth);
    }
    return L / static_cast<float>(rays.size());
}

// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - `computeLightContribution()` and its submethods
// - `renderRaySpecularComponent()`, `renderRayTransparentComponent()`, `renderRayGlossyComponent()`
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth)
{
    // Trace the ray into the scene. If nothing was hit, return early
    HitInfo hitInfo;
    if (!state.bvh.intersect(state, ray, hitInfo)) {
        drawRay(ray, glm::vec3(1, 0, 0));
        return sampleEnvironmentMap(state, ray);
    }

    // Return value: the light along the ray
    // Given an intersection, estimate the contribution of scene lights at this intersection
    glm::vec3 Lo = computeLightContribution(state, ray, hitInfo);

    // Draw an example debug ray for the incident ray (feel free to modify this for yourself)   <---------------------
    // Ray newRay = { .origin = glm::vec3(1, 0, 2), .direction = glm::vec3(0, 0, 1), .t = 3.0 };
    drawRay(ray, glm::vec3(10.0f, 1.0f, 0.0f));

    // Given that recursive components are enabled, and we have not exceeded maximum depth,
    // estimate the contribution along these components
    if (rayDepth < 6) {
        bool isReflective = glm::any(glm::notEqual(hitInfo.material.ks, glm::vec3(0.0f)));
        bool isTransparent = hitInfo.material.transparency != 1.f;

        // Default, specular reflections
        if (state.features.enableReflections && !state.features.extra.enableGlossyReflection && isReflective) {
            renderRaySpecularComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Alternative, glossy reflections
        if (state.features.enableReflections && state.features.extra.enableGlossyReflection && isReflective) {
            renderRayGlossyComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Transparency passthrough
        if (state.features.enableTransparency && isTransparent) {
            renderRayTransparentComponent(state, ray, hitInfo, Lo, rayDepth);
        }
    }

    return Lo;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a reflected ray
// This method is unit-tested, so do not change the function signature.
Ray generateReflectionRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a mirrored ray
    //       if you use glm::reflect, you will not get points for this method!
    glm::vec3 normalizedDirection = glm::normalize(ray.direction);
    glm::vec3 reflectionDirection = normalizedDirection - 2 * glm::dot(normalizedDirection, hitInfo.normal) * hitInfo.normal;

    Ray reflected;
    reflected.origin = ray.origin + ray.t * ray.direction + hitInfo.normal * 0.001f;        //offset 0.001
    reflected.direction = reflectionDirection;
    reflected.t = std::numeric_limits<float>::max();

    return reflected;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a passthrough ray for transparency
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a passthrough ray

    Ray r;

    r.origin = ray.origin + ray.t * ray.direction + ray.origin * 0.001f; // offset
    r.direction = ray.direction;
    r.t = std::numeric_limits<float>::max();

    return r;
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a mirrored ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generateReflectionRay()
    Ray r = generateReflectionRay(ray, hitInfo);

    glm::vec3 colorRay = renderRay(state, r, rayDepth + 1);
    hitColor = hitColor + colorRay * hitInfo.material.ks;
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generatePassthroughRay()
    Ray r = generatePassthroughRay(ray, hitInfo);
    hitColor = hitInfo.material.transparency * hitColor + (1.0f - hitInfo.material.transparency) * renderRay(state, r, rayDepth + 1);

}
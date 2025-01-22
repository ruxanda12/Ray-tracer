#include "bvh.h"
#include "draw.h"
#include "extra.h"
#include "glm/common.hpp"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>
#include <limits>
#include <tuple>

// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex);

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    const auto min = glm::min(primitive.v0.position, glm::min(primitive.v1.position, primitive.v2.position));
    const auto max = glm::max(primitive.v0.position, glm::max(primitive.v1.position, primitive.v2.position));
    return { .lower = min, .upper = max };
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    glm::vec3 min = glm::vec3(std::numeric_limits<float>::max());
    glm::vec3 max = glm::vec3(std::numeric_limits<float>::min());
    for (auto prim : primitives) {
        const auto aabb = computePrimitiveAABB(prim);
        min = glm::min(aabb.lower, min);
        max = glm::max(aabb.upper, max);
    }
    return { .lower = glm::min(min, max), .upper = max };
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    return (1.f / 3.f) * (primitive.v0.position + primitive.v1.position + primitive.v2.position);
}

// More efficient method
// We only need to calculate the axis we actually use to compare
float computePrimitiveCentroid(const BVHInterface::Primitive primitive, const uint32_t axis)
{
    return (1.f / 3.f) * (primitive.v0.position[(int)axis] + primitive.v1.position[(int)axis] + primitive.v2.position[(int)axis]);
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    glm::vec3 axisLengths = glm::abs(aabb.upper - aabb.lower);
    uint32_t longestAxis = 0;
    if (axisLengths.y > axisLengths.x)
        longestAxis = 1;
    if (axisLengths.z > axisLengths[(int)longestAxis])
        longestAxis = 2;
    return longestAxis;
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    using Primitive = BVHInterface::Primitive;

    if (primitives.size() <= 1)
        return 0;

    std::sort(primitives.begin(), primitives.end(), [axis](Primitive a, Primitive b) {
        return computePrimitiveCentroid(a, axis) < computePrimitiveCentroid(b, axis);
    });

    // primitive.size() / 2, rounded up.
    return primitives.size() - (primitives.size() / 2);
}

[[nodiscard]] bool intersectRecursive(RenderState& state, Ray& ray, HitInfo& hitInfo, std::span<const BVHInterface::Node> nodes, std::span<const BVHInterface::Primitive> primitives, uint32_t node)
{
    // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
    // data is not easily extracted. Helper methods are available, however:
    // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
    // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
    // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
    //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
    //
    // In short, you will have to step down the bvh, node by node, and intersect your ray
    // with the node's AABB. If this intersection passes, you should:
    // - if the node is a leaf, intersect with the leaf's primitives
    // - if the node is not a leaf, test the left and right children as well!
    //
    // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
    // and it is likewise possible for a ray to hit both children of a node.
    float old_t = ray.t;
    bool hitAABB = intersectRayWithShape(nodes[node].aabb, ray);
    ray.t = old_t;

    if (!hitAABB)
        return false;
    if (nodes[node].isLeaf()) {
        auto prims = primitives.subspan(nodes[node].primitiveOffset(), nodes[node].primitiveCount());
        bool is_hit = false;
        for (const auto& prim : prims) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
        return is_hit;
    } else {
        auto a = intersectRecursive(state, ray, hitInfo, nodes, primitives, nodes[node].leftChild());
        auto b = intersectRecursive(state, ray, hitInfo, nodes, primitives, nodes[node].rightChild());
        return a || b;
    }
}

// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
        // TODO: implement here your (probably stack-based) BVH traversal.
        //
        is_hit |= intersectRecursive(state, ray, hitInfo, nodes, primitives, 0);

    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);

    return is_hit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    Node node;
    node.aabb = computeSpanAABB(primitives);
    node.data[0] = (uint32_t)m_primitives.size() | Node::LeafBit;
    node.data[1] = (uint32_t)primitives.size();

    // Copy the current set of primitives to the back of the primitives vector
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));

    return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    Node node;
    node.aabb = aabb;
    node.data[0] = leftChildIndex & (~Node::LeafBit);
    node.data[1] = rightChildIndex;
    return node;
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    AxisAlignedBox aabb = computeSpanAABB(primitives);

    // As a starting point, we provide an implementation which creates a single leaf, and stores
    // all triangles inside it. You should remove or comment this, and work on your own recursive
    // construction algorithm that implements the following steps. Make sure to reuse the methods
    // you have previously implemented to simplify this process.
    //
    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    //    (hint; use the `LeafSize` constant)
    // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
    // 3. If it is a node:
    //    3a. Split the range of triangles along the longest axis into left and right subspans,
    //        using either median or SAH-Binning based on the `Features` object
    //    3b. Allocate left/right child nodes
    //        (hint: use `nextNodeIdx()`)
    //    3c. Fill in the current node's data; aabb, left/right child indices
    //    3d. Recursively build left/right child nodes over their respective triangles
    //        (hint; use `std::span::subspan()` to split into left/right ranges)
    // Just configure the current node as a giant leaf for now
    if (primitives.size() > BVH::LeafSize) {
        uint32_t longestAxis = computeAABBLongestAxis(aabb);
        size_t splitIndex = (features.extra.enableBvhSahBinning ? splitPrimitivesBySAHBin(aabb, longestAxis, primitives)
                                                                : splitPrimitivesByMedian(aabb, longestAxis, primitives));
        uint32_t leftNode = nextNodeIdx();
        uint32_t rightNode = nextNodeIdx();
        m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftNode, rightNode);
        buildRecursive(scene, features, primitives.subspan(0, splitIndex), leftNode);
        buildRecursive(scene, features, primitives.subspan(splitIndex), rightNode);
    } else {
        m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
    }
}

uint32_t getNumLevels(std::span<BVHInterface::Node> nodes, uint32_t nodeIndex)
{
    auto current = nodes[nodeIndex];
    uint32_t currentLevel = current.isLeaf() ? 0 : glm::max(getNumLevels(nodes, current.leftChild()), getNumLevels(nodes, current.rightChild()));
    return currentLevel + 1;
}

// TODO: Standard feature, or part of it
// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
    m_numLevels = getNumLevels(m_nodes, 0);
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    uint32_t numLeaves = 0;
    for (auto node : m_nodes) {
        numLeaves += node.isLeaf();
    }
    m_numLeaves = numLeaves;
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.

void drawLevelRecursive(uint32_t level, uint32_t currentLevel, std::span<BVHInterface::Node> nodes, uint32_t nodeIndex)
{
    auto current = nodes[nodeIndex];
    if (current.isLeaf())
        return;
    if (currentLevel == level) {
        drawAABB(current.aabb, DrawMode::Wireframe, glm::vec3(1.f, 0.05f, 0.05f), 1.0f);
    } else {
        drawLevelRecursive(level, currentLevel + 1, nodes, current.leftChild());
        drawLevelRecursive(level, currentLevel + 1, nodes, current.rightChild());
    }
}

void BVH::debugDrawLevel(int level)
{
    if ((uint32_t)level >= m_numLevels || level < 0)
        return;

    drawLevelRecursive((uint32_t)level, 0, m_nodes, 0);
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
    // colors, transparencies, etc.

    // drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use drawTriangle (see `draw.h`) to draw the contained primitives
    uint32_t curLeafIndex = 0;
    uint32_t nodeIndex;
    for (nodeIndex = 0; nodeIndex < m_nodes.size(); nodeIndex++) {
        if (!m_nodes[nodeIndex].isLeaf())
            continue;
        curLeafIndex++;
        if (curLeafIndex == (uint32_t)leafIndex)
            break;
    }
    // We actually found a leaf node. This should always be true, but still sanity check.
    if (nodeIndex < m_nodes.size()) {
        drawAABB(m_nodes[nodeIndex].aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 1.f);
        uint32_t off = m_nodes[nodeIndex].primitiveOffset(), count = m_nodes[nodeIndex].primitiveCount();

        auto color = glm::vec3(0.05f, 1.0f, 0.05f);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glColor4f(color.r, color.g, color.b, 1.f);
        glPolygonMode(GL_FRONT, GL_FILL);
        glPolygonMode(GL_BACK, GL_FILL);
        for (uint32_t i = off; i < off + count; i++) {
            auto prim = m_primitives[i];
            const auto [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            drawTriangle(v0, v1, v2);
        }
        glPopAttrib();
    }
}

#pragma once
#include "bvh_interface.h"
#include <framework/ray.h>
#include <vector>

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive);

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives);

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive);

float computePrimitiveCentroid(const BVHInterface::Primitive primitive, const uint32_t axis);

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis as (x = 0, y = 1, z = 2)
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb);

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element
// of the list is at least as big as the other, and both differ at most by one element in size.
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives);

// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(), you must implement this method.
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo);

// The implementing class where you will put most of the BVH implementation; this class must conform
// to BVHInterface for grading purposes; see `bvh_interface.h` for details
struct BVH : public BVHInterface {
    // Constants used throughout the BVH
    static constexpr uint32_t LeafSize = 4; // Maximum nr. of primitives in a leaf
    static constexpr uint32_t RootIndex = 0; // Index of root node in `m_nodes` vector

    // Constructor. Receives the scene and starts the build process
    // NOTE: this constructor is used in tests, so do not change its function signature.
    BVH(const Scene& scene, const Features& features);

    // See BVHInterface::intersect(...) for argument descriptions
    bool intersect(RenderState& state, Ray& ray, HitInfo& hitInfo) const override
    {
        return intersectRayWithBVH(state, *this, ray, hitInfo);
    }

private: // Private members
    uint32_t m_numLevels;
    uint32_t m_numLeaves;
    std::vector<Node> m_nodes;
    std::vector<Primitive> m_primitives;

private: // Private methods
    // Helper method; simply allocates a new node, and returns its index
    uint32_t nextNodeIdx();

    // TODO: Standard feature
    // Helper functions to instantiate Node objects as either parent nodes with children, or as leaf nodes
    // For a description of either method's arguments, refer to 'bounding_volume_hierarchy.cpp'
    // Note: you are free to modify these function's signatures, as long as the constructor builds a BVH
    Node buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives);
    Node buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex);

    // TODO: Standard feature
    // Hierarchy construction routine; called by the BVH's constructor, you must implement this method.
    // For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
    // Note: you are free to modify this function's signature, as long as the constructor builds a BVH
    void buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex);

private: // Visual debug helpers
    // Compute the nr. of levels in your hierarchy after construction; useful for debugDrawLevel()
    // You are free to modify this function's signature, as long as the constructor builds a BVH
    void buildNumLevels();

    // Compute the nr. of leaves in your hierarchy after construction; useful for debugDrawLeaf()
    // You are free to modify this function's signature, as long as the constructor builds a BVH
    void buildNumLeaves();

public: // Visual debug
    // Draw the bounding boxes of the nodes at the selected level.
    // For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
    // You are free to modify this function's signature.
    void debugDrawLevel(int level);

    // Draw data of the leaf at the selected index.
    // For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
    // You are free to modify this function's signature.
    void debugDrawLeaf(int leafIndex);

public: // Public getters
    // Accessors to underlying data
    std::span<const Node> nodes() const override { return m_nodes; }
    std::span<Node> nodes() override { return m_nodes; }
    std::span<const Primitive> primitives() const override { return m_primitives; }
    std::span<Primitive> primitives() override { return m_primitives; }

    // Return how many levels/leaves there are in the tree
    uint32_t numLevels() const override { return m_numLevels; }
    uint32_t numLeaves() const override { return m_numLeaves; }
};

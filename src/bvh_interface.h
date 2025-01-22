#pragma once
#include "common.h"
#include "fwd.h"
#include <array>
#include <span>
#include <framework/ray.h>

//! ********** DO NOT MODIFY THIS FILE (IT EXISTS FOR GRADING PURPOSES)! **************** //
//! ********** YOU SHOULD PROBABLY READ THROUGH IT CAREFULLY THOUGH! ******************** //

// The base interface to which your BVH implementation must conform for tests
struct BVHInterface {
    // A primitive represents a triangle stored inside the BVH's leaf nodes
    struct Primitive {
        // Index of scene mesh from which the primitive's vertices are sourced
        uint32_t meshID;

        // Stored vertices, forming a triangle
        Vertex v0, v1, v2;

        // Default equality operator (not relevant for understanding the code)
        [[nodiscard]] constexpr bool operator==(const Primitive&) const noexcept = default;
    };

    // Packed BVH node; a node either has two children, or it is a leaf, in
    // which case it refers to one or more primitives (triangles). Read this
    // carefully, and make sure you understand how the data is packed with LeafBit!
    struct Node {
        // A flag bit used to distinguish nodes and leaves
        static constexpr uint32_t LeafBit = 1u << 31;

        // Bounding box around the node's contained primitives
        AxisAlignedBox aabb;

        // Node data; the first integer's most significant bit is used as a flag
        // to determine the node's purpose as a leaf or node.
        // In short, the layout is as follows:
        // - node: [[0, index of left child], [index of right child]]
        // - leaf: [[1, offset to primitive], [count of primitives]]
        std::array<uint32_t, 2> data;

        // Return if the node is a leaf node
        [[nodiscard]] inline constexpr bool isLeaf() const
        {
            return (data[0] & LeafBit) == LeafBit;
        };

    public: // Getters
        [[nodiscard]] inline constexpr uint32_t primitiveOffset() const { return data[0] & (~LeafBit); }
        [[nodiscard]] inline constexpr uint32_t primitiveCount() const { return data[1]; }
        [[nodiscard]] inline constexpr uint32_t leftChild() const { return data[0]; }
        [[nodiscard]] inline constexpr uint32_t rightChild() const { return data[1]; }
    };
    static_assert(sizeof(Node) == 32);

public:
    virtual ~BVHInterface() {};

    // Return true if something is hit, false otherwise. On a hit, the distance 't' on the
    // ray object is updated, and hit information is stored in the 'hitInfo' object.
    virtual bool intersect(RenderState& state, Ray& ray, HitInfo& hitInfo) const = 0;

    // Accessors to underlying data
    virtual std::span<const Node> nodes() const = 0;
    virtual std::span<Node> nodes() = 0;
    virtual std::span<const Primitive> primitives() const = 0;
    virtual std::span<Primitive> primitives() = 0;

    // Accessors to tree structure layout
    virtual uint32_t numLevels() const = 0;
    virtual uint32_t numLeaves() const = 0;
};

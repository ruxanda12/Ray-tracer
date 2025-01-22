#include "bvh.h"
#include "common.h"
#include "render.h"
#include "scene.h"
#include <concepts>
#include <type_traits>

//! ********** DO NOT MODIFY THIS FILE (IT EXISTS FOR MISTAKE_CATCHING PURPOSES)! ********* //
//! ********** IF YOU GET A COMPILATION ERROR IN THIS FILE, SCROLL DOWN FOR HELP ********** //

template <typename Ty>
concept is_scene = requires(Ty t) {
    {
        t.type
    } -> std::same_as<SceneType&>;
    {
        t.meshes
    } -> std::same_as<std::vector<Mesh>&>;
    {
        t.spheres
    } -> std::same_as<std::vector<Sphere>&>;
    {
        t.lights
    } -> std::same_as<std::vector<typename Ty::SceneLight>&>;
};

template <typename Ty>
concept is_render_state = requires(Ty t) {
    {
        t.scene
    } -> std::same_as<const Scene&>;
    {
        t.features
    } -> std::same_as<const Features&>;
    {
        t.bvh
    } -> std::same_as<const BVHInterface&>;
    {
        t.sampler
    } -> std::same_as<Sampler&>;
};

template <typename Ty>
concept is_features = requires(Ty t) {
    {
        t.enableShading
    } -> std::same_as<bool&>;
    {
        t.enableReflections
    } -> std::same_as<bool&>;
    {
        t.enableShadows
    } -> std::same_as<bool&>;
    {
        t.enableNormalInterp
    } -> std::same_as<bool&>;
    {
        t.enableTextureMapping
    } -> std::same_as<bool&>;
    {
        t.enableAccelStructure
    } -> std::same_as<bool&>;
    {
        t.enableBilinearTextureFiltering
    } -> std::same_as<bool&>;
    {
        t.enableTransparency
    } -> std::same_as<bool&>;
    {
        t.shadingModel
    } -> std::same_as<ShadingModel&>;
    {
        t.numPixelSamples
    } -> std::same_as<uint32_t&>;
    {
        t.numShadowSamples
    } -> std::same_as<uint32_t&>;
    {
        t.extra
    } -> std::same_as<ExtraFeatures&>;
};

template <typename Ty>
concept is_bvh_interface_node = requires(Ty t) {
    {
        t.aabb
    } -> std::same_as<AxisAlignedBox&>;
    {
        t.data
    } -> std::same_as<std::array<uint32_t, 2>&>;
};

template <typename Ty>
concept is_bvh_interface_primitive = requires(Ty t) {
    {
        t.meshID
    } -> std::same_as<uint32_t&>;
    {
        t.v0
    } -> std::same_as<Vertex&>;
};

template <typename Ty>
concept is_bvh_interface = requires(Ty t, RenderState& state, Ray& ray, HitInfo& hitInfo) {
    {
        t.nodes()
    } -> std::same_as<std::span<BVHInterface::Node>>;
    {
        t.primitives()
    } -> std::same_as<std::span<BVHInterface::Primitive>>;
    {
        t.numLevels()
    } -> std::same_as<uint32_t>;
    {
        t.numLeaves()
    } -> std::same_as<uint32_t>;
    {
        t.intersect(state, ray, hitInfo)
    } -> std::same_as<bool>;
};

// Be warned, traveler!
// You wandered into the cavern of the vile dragon "I-modified-the-BVH-interface-or-RenderState".
// Unfortunately, you rolled a nat. 1 on your stealth check, and the heinous dragon has noticed
// your trespassing. To escape the creature's fiery wrath, you must undo your modifications and
// make these assertions pass once moooooooore!

// For the people who do have a life: if any of these assertions fails, you have changed something you
// should leave alone, as it may break our grading tests. Please undo this change, or you may lose points.

static_assert(is_scene<Scene>,
    "You broke Scene. Unbreak it or fail grading tests!");
static_assert(is_render_state<RenderState>,
    "You broke RenderState. Unbreak it or fail grading tests!");
static_assert(is_features<Features>,
    "You broke Features. Unbreak it or fail grading tests!");
static_assert(is_bvh_interface<BVHInterface>,
    "You broke BVHInterface. Unbreak it or fail grading tests!");
static_assert(is_bvh_interface_node<BVHInterface::Node>,
    "You broke BVHInterface::Node. Unbreak it or fail grading tests!");
static_assert(is_bvh_interface_primitive<BVHInterface::Primitive>,
    "You broke BVHInterface::Primitive. Unbreak it or fail grading tests!");
static_assert(std::is_base_of_v<BVHInterface, BVH>,
    "You broke BVH's BVHInterface compliance. Unbreak it or fail grading tests!");

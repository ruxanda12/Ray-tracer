#include "scene.h"
#include <cmath>
#include <iostream>

Scene loadScenePrebuilt(SceneType type, const std::filesystem::path& dataDir)
{
    Scene scene;
    scene.type = type;
    switch (type) {
    case SingleTriangle: {
        // Load a 3D model with a single triangle
        std::cout << "begin tri" << std::endl;
        auto subMeshes = loadMesh(dataDir / "triangle.obj");
        subMeshes[0].material.kd = glm::vec3(1.0f);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        std::filesystem::path imagePath = "../../../data/small_empty_room_3_4k.jpeg";
        try {
            std::cout << std::filesystem::current_path() << std::endl;
            std::cout << "loaded img : " << imagePath << std::endl;
            Image img = Image(imagePath);
            scene.environmentMap.emplace(img);
        } catch (const std::exception& e) {
            std::cerr << "Failed to load environment map: " << e.what() << std::endl;
            // Handle the error, perhaps set a default environment map or abort
        }
    } break;
    case Cube: {
        // Load a 3D model of a cube with 12 triangles
        auto subMeshes = loadMesh(dataDir / "cube.obj");
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        // scene.lights.push_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        scene.lights.emplace_back(SegmentLight {
            .endpoint0 = glm::vec3(1.5f, 0.5f, -0.6f),
            .endpoint1 = glm::vec3(-1, 0.5f, -0.5f),
            .color0 = glm::vec3(0.9f, 0.2f, 0.1f), // Red-ish
            .color1 = glm::vec3(0.2f, 1, 0.3f) // Green-ish
        });
    } break;
    case CubeTextured: {
        auto subMeshes = loadMesh(dataDir / "cube-textured.obj");
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1.0, 1.5, -1.0), glm::vec3(1) });
    } break;
    case CornellBox: {
        // Load a 3D model of a Cornell Box
        auto subMeshes = loadMesh(dataDir / "CornellBox-Mirror-Rotated.obj", true);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(0, 0.58f, 0), glm::vec3(1) }); // Light at the top of the box
    } break;
    case CornellBoxTransparency: {
        // Load a 3D model of a Cornell Box
        auto subMeshes = loadMesh(dataDir / "CornellBox-Mirror-Rotated.obj", true);
        // for (auto &mesh : subMeshes)
        //     mesh.material.transparency = 0.5f;
        subMeshes[6].material = Material {
            .kd = glm::vec3(1, 0.25, 0.25),
            .ks = glm::vec3(0),
            .transparency = 0.5f
        };
        subMeshes[5].material = Material {
            .kd = glm::vec3(0.25, 1, 0.25),
            .ks = glm::vec3(0),
            .transparency = 0.5f
        };
        // subMeshes[6].material.transparency = 0.5;
        // subMeshes[7].material.transparency = 0.5;
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(0, 0.58f, 0), glm::vec3(1) }); // Light at the top of the box
    } break;
    case CornellBoxParallelogramLight: {
        // Load a 3D model of a Cornell Box
        auto subMeshes = loadMesh(dataDir / "CornellBox-Mirror-Rotated.obj", true);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        // Light at the top of the box.
        scene.lights.emplace_back(ParallelogramLight {
            .v0 = glm::vec3(-0.2f, 0.5f, 0),
            .edge01 = glm::vec3(0.4f, 0, 0),
            .edge02 = glm::vec3(0.0f, 0.0f, 0.4f),
            .color0 = glm::vec3(1, 0, 0), // Red
            .color1 = glm::vec3(0, 1, 0), // Green
            .color2 = glm::vec3(0, 0, 1), // Blue
            .color3 = glm::vec3(0, 1, 1), // Purple
        });
    } break;
    case Monkey: {
        // Load a 3D model of a Monkey
        auto subMeshes = loadMesh(dataDir / "monkey.obj", true);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        scene.lights.emplace_back(PointLight { glm::vec3(1, -1, -1), glm::vec3(1) });
    } break;
    case Teapot: {
        // Load a 3D model of a Teapot
        auto subMeshes = loadMesh(dataDir / "teapot.obj", true);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
    } break;
    case Dragon: {
        // Load a 3D model of a Dragon
        auto subMeshes = loadMesh(dataDir / "dragon.obj", true);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
    } break;
    case Spheres: {
        scene.spheres.push_back(Sphere { glm::vec3(3.0f, -2.0f, 10.2f), 1.0f, Material { glm::vec3(0.8f, 0.2f, 0.2f) } });
        scene.spheres.push_back(Sphere { glm::vec3(-2.0f, 2.0f, 4.0f), 2.0f, Material { glm::vec3(0.6f, 0.8f, 0.2f) } });
        scene.spheres.push_back(Sphere { glm::vec3(0.0f, 0.0f, 6.0f), 0.75f, Material { glm::vec3(0.2f, 0.2f, 0.8f) } });
        scene.lights.emplace_back(PointLight { glm::vec3(3, 0, 3), glm::vec3(1) });
    } break;
    case Custom: {
        // Load a 3D model of a Monkey
        auto subMeshes = loadMesh(dataDir / "monkey.obj", true);
        std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));
        scene.lights.emplace_back(PointLight { glm::vec3(-1, 1, -1), glm::vec3(1) });
        scene.lights.emplace_back(PointLight { glm::vec3(1, -1, -1), glm::vec3(1) });
    } break;
    };

    return scene;
}

Scene loadSceneFromFile(const std::filesystem::path& path, const std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight>>& lights)
{
    Scene scene;
    scene.lights = std::move(lights);

    auto subMeshes = loadMesh(path);
    std::move(std::begin(subMeshes), std::end(subMeshes), std::back_inserter(scene.meshes));

    return scene;
}

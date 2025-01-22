# Ray Tracing Project

Implemented as part of the **Computer Graphics** course at TU Delft, this ray tracer is capable of rendering realistic 3D scenes. Features include environment mapping for realistic reflections, motion blur, glossy reflections, bloom filters for enhanced lighting effects, and BVH optimizations for efficient rendering.

## Features

- **Shading Models**:
  - Implemented Phong and Blinn-Phong shading models for realistic lighting.
  - Introduced a custom linear gradient shading model based on light direction and material properties.

- **Texture Sampling**:
  - Nearest neighbor and bilinear interpolation for texture sampling.
  - Enhanced visual fidelity through efficient texture handling.

- **Performance Optimizations**:
  - Multi-threading support for faster rendering.
  - Bounding Volume Hierarchy (BVH) optimizations for efficient ray traversal.

- **Visual Effects**:
  - Glossy reflections using randomized sampling techniques.
  - Motion blur for dynamic scenes with moving objects.
  - Bloom filter for simulating light scattering effects.
  - Environment mapping for realistic background reflections.

## Technologies Used

- **Programming Language**: C++
- **Graphics API**: OpenGL
- **Mathematics Library**: GLM (OpenGL Mathematics)
- **Development Tools**: Visual Studio, CMake

## Contribution

This project involved enhancing existing source files to add new features and optimize performance. Specific contributions include:

- Developing advanced shading techniques in `shading.cpp`.
- Implementing texture sampling methods in `texture.cpp` and `texture.h`.
- Optimizing ray tracing performance with BVH.
- Adding visual effects like motion blur and glossy reflections.

## How to Build and Run

### Prerequisites

- C++17 or later
- CMake 3.15 or later
- OpenGL-compatible hardware

### Steps to Build

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/advanced-ray-tracing.git
   cd advanced-ray-tracing
   ```

2. Create a build directory and configure the project:
   ```bash
   mkdir build && cd build
   cmake ..
   ```

3. Build the project:
   ```bash
   cmake --build .
   ```

4. Run the application:
   ```bash
   ./raytracer
   ```

## Team Members

This project was collaboratively developed by:

- Maria Ruxanda Tudor
- Nathan Huisman
- Ísak Bieltvedt Jónsson


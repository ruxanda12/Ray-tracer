#include "image.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
DISABLE_WARNINGS_POP()
#include <cassert>
#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include <glm/glm.hpp>
#include <filesystem>
#include <algorithm>


Image::Image(const std::filesystem::path& filePath)
{
    if (!std::filesystem::exists(filePath)) {
        std::cerr << "Texture file " << filePath << " does not exists!" << std::endl;
        throw std::runtime_error("Texture file does not exist");
    }

    int numChannelsInSourceImage;
    stbi_uc* stbPixels = stbi_load(filePath.string().c_str(), &width, &height, &numChannelsInSourceImage, STBI_rgb);

    if (!stbPixels) {
        std::cerr << "Failed to read texture " << filePath << " using stb_image.h" << std::endl;
        throw std::runtime_error("Failed to load texture");
    }

    constexpr size_t numChannels = 3; // STBI_rgb == 3 channels
    pixels.reserve(width * height);
    for (size_t i = 0; i < width * height * numChannels; i += numChannels) {
        pixels.emplace_back(stbPixels[i] / 255.0f, stbPixels[i + 1] / 255.0f, stbPixels[i + 2] / 255.0f);
    }

    stbi_image_free(stbPixels);
}

void Image::writeBitmapToFile(const std::filesystem::path& filePath)
{
    std::vector<glm::u8vec4> textureData8Bits(pixels.size());

    std::transform(pixels.begin(), pixels.end(), textureData8Bits.begin(),
        [](const glm::vec3& color) {
            const glm::vec3 clampedColor = glm::clamp(color, 0.0f, 1.0f);
            return glm::u8vec4(glm::vec4(clampedColor, 1.0f) * 255.0f);
        });

    stbi_write_bmp(filePath.string().c_str(), width, height, 4, textureData8Bits.data());
}

#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image

    int x = int((texCoord.x * image.width) - 0.5f);
    int y = int((texCoord.y * image.height) - 0.5f);

    x = glm::clamp(x, 0, image.width - 1);     // make sure we respect the bounds
    y = glm::clamp(y, 0, image.height - 1);

    int index = y * image.width + x;

    return image.pixels[index];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image

    int x1 = int((texCoord.x * image.width) - 0.5f);
    int y1 = int((texCoord.y * image.height) - 0.5f);
    int x2 = x1 + 1;
    int y2 = y1 + 1;

    float fr_x = (texCoord.x * image.width) - 0.5f - x1;       //fractional part
    float fr_y = (texCoord.y * image.height) - 0.5f - y1;

    x1 = glm::clamp(x1, 0, image.width - 1);
    y1 = glm::clamp(y1, 0, image.height - 1);
    x2 = glm::clamp(x2, 0, image.width - 1);
    y2 = glm::clamp(y2, 0, image.height - 1);

    glm::vec3 texel1 = image.pixels[y1 * image.width + x1];
    glm::vec3 texel2 = image.pixels[y1 * image.width + x2];
    glm::vec3 texel3 = image.pixels[y2 * image.width + x1];
    glm::vec3 texel4 = image.pixels[y2 * image.width + x2];

    glm::vec3 result = texel1 * (1.0f - fr_x) * (1.0f - fr_y) +
        texel2 * fr_x * (1.0f - fr_y) +
        texel3 * (1.0f - fr_x) * fr_y +
        texel4 * fr_x * fr_y;

    return result;
}
#pragma once

#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
DISABLE_WARNINGS_POP()
#include <random>
#include <utility>

// Simple fast 1d/2d sampler for drawing uniformly distributed samples in [0, 1),
// based on the pcg integer hash; https://www.pcg-random.org/
// - Not thread-safe, do not share across threads.
class Sampler {
    uint32_t m_state;

    uint32_t pcg_hash(uint32_t& state)
    {
        state = state * 747796405u + 2891336453u;

        uint32_t v = state;
        v ^= v >> ((v >> 28u) + 4u);
        v *= 277803737u;
        v ^= v >> 22u;
        return v;
    }

public:
    // Seeded constructor, by default draws from std::random_device
    Sampler(uint32_t seed = std::random_device()())
        : m_state(seed)
    {
        // ...
    }

    // Draw a 1d sample in [a, b)
    float next_1d()
    {
        return static_cast<float>(pcg_hash(m_state)) / 4294967295.f;
    }

    // Draw a 2d sample in [a, b)
    glm::vec2 next_2d()
    {
        return { next_1d(), next_1d() };
    }
};
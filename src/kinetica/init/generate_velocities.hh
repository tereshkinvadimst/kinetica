#ifndef MF_KINETICA_INITIALIZE_GENERATE_VELOCITIES_H
#define MF_KINETICA_INITIALIZE_GENERATE_VELOCITIES_H
#pragma once
#include <concepts>
#include <random>
#include <ranges>
#include <algorithm>

namespace mf {

template <std::floating_point T, std::uniform_random_bit_generator RAND>
inline T generateMaxwellVelocity(T velocity_scale_factor, RAND& gen) {
    std::normal_distribution<T> dist(T(0), velocity_scale_factor);
    return dist(gen);
}

template <std::ranges::range R, std::uniform_random_bit_generator RAND>
    requires std::floating_point<std::ranges::range_value_t<R>>
void generateMaxwellVelocity(R&& out, std::ranges::range_value_t<R> velocity_scale_factor, RAND& gen) {
    using T = std::ranges::range_value_t<R>;

    std::normal_distribution<T> dist(T{}, velocity_scale_factor);

    std::ranges::generate(out, [&dist, &gen] { return dist(gen); });
}

}  // namespace mf

#endif  // MF_KINETICA_INITIALIZE_GENERATE_VELOCITIES_H
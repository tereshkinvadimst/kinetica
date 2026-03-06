#ifndef MF_KINETICA_DSMC_UTILS_H
#define MF_KINETICA_DSMC_UTILS_H
#pragma once

#include <cmath>
#include <concepts>
#include <numbers>
#include <random>

#include "kinetica/constants.hh"

namespace mf {

template <std::floating_point T>
[[nodiscard]] constexpr auto meanFreePath(T n_density, T d) noexcept -> T {
    const auto sigma = std::numbers::pi_v<T> * d * d;
    return T{1} / (std::sqrt(T{2}) * n_density * sigma);
}

template <std::floating_point T, std::uniform_random_bit_generator RNG>
[[nodiscard]] inline auto computeNParticles(T number_density, T volume, T weight, RNG& gen) -> std::size_t {
    const auto                             N_avg = number_density * volume / weight;

    std::poisson_distribution<std::size_t> poisson(static_cast<double>(N_avg));

    return poisson(gen);
}

template <std::floating_point T>
[[nodiscard]] constexpr auto computeVelocityScaleFactor(T temperature, T mass) noexcept -> T {
    return std::sqrt(K_B * temperature / mass);
}

template <std::floating_point T>
[[nodiscard]] constexpr auto computeTimeStep(T mean_free_path, T max_velocity, T scale_factor) noexcept -> T {
#ifndef NDEBUG
    assert(mean_free_path > T{});
    assert(max_velocity > T{});
    assert(scale_factor > T{} && scale_factor <= T{1});
#endif
    return scale_factor * mean_free_path / max_velocity;
}

template <std::floating_point T>
[[nodiscard]] constexpr auto computeCellSize(T mean_free_path, T scale_factor) noexcept -> T {
#ifndef NDEBUG
    assert(mean_free_path > T{});
    assert(scale_factor > T{} && scale_factor <= T{1});
#endif
    return scale_factor * mean_free_path;
}

}  // namespace mf

#endif  // MF_KINETICA_DSMC_UTILS_H
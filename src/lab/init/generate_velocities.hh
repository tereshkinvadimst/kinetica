#ifndef MF_KINETICA_INITIALIZE_GENERATE_VELOCITIES_H
#define MF_KINETICA_INITIALIZE_GENERATE_VELOCITIES_H
#pragma once
#include <concepts>
#include <random>

#include "lab/Random/xoshiro256.hh"

namespace mf {

template <std::floating_point T>
inline T generateMaxwellVelocity(T velocity_scale_factor, xoshiro256& gen) {
    std::normal_distribution<T> dist(T(0), velocity_scale_factor);
    return dist(gen);
}

template <std::ranges::output_range<double> R, std::floating_point T>
void generateMaxwellVelocities(R&& out, T velocity_scale_factor, xoshiro256& gen) {
    std::normal_distribution<T> dist(T(0), velocity_scale_factor);

    std::ranges::generate(out, [&] { return dist(gen); });
}

}  // namespace mf

#endif  // MF_KINETICA_INITIALIZE_GENERATE_VELOCITIES_H
#ifndef MF_KINETICA_INITIALIZE_GENERATE_POSITIONS_H
#define MF_KINETICA_INITIALIZE_GENERATE_POSITIONS_H
#pragma once
#include <concepts>
#include <random>
#include <ranges>

#include "lab/Random/xoshiro256.hh"

namespace mf {

template <std::ranges::output_range<double> R, std::floating_point T>
void generateUniformPositions(R&& out, T left, T right, xoshiro256& gen) {
#ifndef NDEBUG
    assert(left < right);
#endif

    std::uniform_real_distribution<T> dist(left, right);

    std::ranges::generate(out, [&gen, &dist] { return dist(gen); });
}

}  // namespace mf

#endif  // MF_KINETICA_INITIALIZE_GENERATE_POSITIONS_H
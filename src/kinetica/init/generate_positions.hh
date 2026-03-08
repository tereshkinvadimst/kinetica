#ifndef MF_KINETICA_INITIALIZE_GENERATE_POSITIONS_H
#define MF_KINETICA_INITIALIZE_GENERATE_POSITIONS_H
#pragma once
#include <concepts>
#include <random>
#include <ranges>

namespace mf {

template <std::ranges::range R, std::uniform_random_bit_generator RAND>
    requires std::floating_point<std::ranges::range_value_t<R>>
inline void generateUniformPositions(R&&                           out,
                                     std::ranges::range_value_t<R> left,
                                     std::ranges::range_value_t<R> right,
                                     RAND&                         gen) noexcept {
    using T = std::ranges::range_value_t<R>;
#ifndef NDEBUG
    assert(left <= right);
#endif

    std::uniform_real_distribution<T> dist(left, right);

    std::ranges::generate(out, [&gen, &dist]() noexcept { return dist(gen); });
}

}  // namespace mf

#endif  // MF_KINETICA_INITIALIZE_GENERATE_POSITIONS_H
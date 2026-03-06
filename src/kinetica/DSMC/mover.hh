#ifndef MF_KINETICA_DSMC_MOVER_H
#define MF_KINETICA_DSMC_MOVER_H
#pragma once
#include <algorithm>
#include <concepts>
#include <execution>
#include <ranges>

namespace mf {

template <std::ranges::range R1, std::ranges::range R2>
    requires std::floating_point<std::ranges::range_value_t<R1>>
inline void mover(R1&& position, const R2& velocity, std::ranges::range_value_t<R1> dt) noexcept {
    std::transform(
        std::execution::unseq, position.begin(), position.end(), velocity.begin(), position.begin(), [dt](auto x, auto v) {
            return x + v * dt;
        });
}

}  // namespace mf

#endif  // MF_KINETICA_DSMC_MOVER_H
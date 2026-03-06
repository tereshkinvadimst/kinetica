#ifndef MF_KINETICA_BOUNDARIES_H
#define MF_KINETICA_BOUNDARIES_H
#pragma once

#include <cassert>
#include <cmath>
#include <concepts>
#include <ranges>

#include "kinetica/Box/box.hh"
#include "kinetica/Particles/particles.hh"
#include "kinetica/Random/random.hh"

namespace mf {

template <std::floating_point T>
constexpr auto applyPeriodic(T x, T left, T right) noexcept -> T {
#ifndef NDEBUG
    assert(left < right);
#endif

    const auto L = right - left;

    T          r = std::fmod(x - left, L);
    if (r < T{}) r += L;
    return left + r;
}

template <std::ranges::range R>
    requires std::floating_point<std::ranges::range_value_t<R>>
void applyPeriodic(R&& xs, std::ranges::range_value_t<R> left, std::ranges::range_value_t<R> right) noexcept {
    for (auto& x : xs) {
        x = applyPeriodic(x, left, right);
    }
}

void applyPeriodic(Particles& particles, Box box, bool px, bool py, bool pz) noexcept;

auto scatterDiffuse(char axis, int sign, double Tw, double m, random& rng) -> std::tuple<double, double, double>;

}  // namespace mf

#endif  // MF_KINETICA_BOUNDARIES_H
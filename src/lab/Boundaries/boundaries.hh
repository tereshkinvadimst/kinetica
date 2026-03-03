#ifndef MF_KINETICA_BOUNDARIES_H
#define MF_KINETICA_BOUNDARIES_H
#include "lab/constants.hh"
#pragma once

#include <cassert>
#include <cmath>
#include <concepts>
#include <ranges>

#include "lab/Box/box.hh"
#include "lab/Particles/particles.hh"
#include "lab/Random/xoshiro256.hh"

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

class DiffuseWallBoundary final {
   public:
    using value_type      = double;
    using size_type       = std::size_t;

    DiffuseWallBoundary() = default;

    explicit DiffuseWallBoundary(Box wall_position, value_type Tw);

    auto scatter(char axis, int sign, double m, xoshiro256& rng) const -> std::tuple<value_type, value_type, value_type>;

    void killParticlesInBox(Particles& particles) const;

    void operator()(Particles& particles, value_type dt, xoshiro256& rng) const;

   private:
    Box        wall_position_;
    value_type Tw_;
};

}  // namespace mf

#endif  // MF_KINETICA_BOUNDARIES_H
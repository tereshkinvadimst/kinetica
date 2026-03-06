#ifndef MF_KINETICA_BOX_H
#define MF_KINETICA_BOX_H
#pragma once

namespace mf {

struct Box final {
    // -------------------------------------------------------------------------
    // Получить центр коробки
    [[nodiscard]] constexpr auto cx() const noexcept -> double { return x0 + 0.5 * Lx; }
    [[nodiscard]] constexpr auto cy() const noexcept -> double { return y0 + 0.5 * Ly; }
    [[nodiscard]] constexpr auto cz() const noexcept -> double { return z0 + 0.5 * Lz; }
    // -------------------------------------------------------------------------

    // Проверка, что точка лежит внутри коробки
    [[nodiscard]] constexpr auto contains(double x, double y, double z) const noexcept -> bool {
        return (x >= x0 && x < x0 + Lx) && (y >= y0 && y < y0 + Ly) && (z >= z0 && z < z0 + Lz);
    }

    [[nodiscard]] constexpr auto volume() const noexcept -> double { return Lx * Ly * Lz; }

    // Положение коробки
    double x0, y0, z0;
    // Размеры коробки
    double Lx, Ly, Lz;
};

}  // namespace mf

#endif  // MF_KINETICA_BOX_H
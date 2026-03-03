#ifndef MF_KINETICA_PARTICLE_H
#define MF_KINETICA_PARTICLE_H
#pragma once
#include <vector>

namespace mf {

struct Particles final {
   public:
    using value_type = double;
    using size_type  = size_t;

    Particles()      = default;

    explicit Particles(size_type n_particles, value_type m_, value_type W_, value_type molecule_size_);

    [[nodiscard]] auto getNParticles() const noexcept -> size_type;

    [[nodiscard]] auto isAlive(size_type i) const -> bool;
    void               killParticle(size_type i);

    auto               sigma(value_type gij) const noexcept -> value_type;

    /// Масса одной молекулы в частице
    double m;
    /// Вес одной частицы
    double W;
    // Размер молекулы
    double molecule_size;
    /// Положение частицы
    std::vector<value_type> x, y, z;
    /// Скорость частицы
    std::vector<value_type> ux, uy, uz;
    /// Указатель на ячейку, в которой находится частица
    std::vector<size_type> cell_id;
    /// Жива ли частица в расчётной области? 1 - жива, 0 - мертва
    std::vector<char> is_alive;
};

}  // namespace mf

#endif  // MF_KINETICA_PARTICLE_H
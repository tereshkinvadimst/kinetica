#include "kinetica/Particles/particles.hh"

#include <numbers>

#if _DEBUG
#include <cassert>
#endif

mf::Particles::Particles(size_type n_particles, value_type m_, value_type W_, value_type molecule_size_)
    : m{m_}
    , W{W_}
    , molecule_size{molecule_size_}
    , x(n_particles)
    , y(n_particles)
    , z(n_particles)
    , ux(n_particles)
    , uy(n_particles)
    , uz(n_particles)
    , cell_id(n_particles)
    , is_alive(n_particles, 1) {}

auto mf::Particles::getNParticles() const noexcept -> size_type {
    const size_type N = static_cast<size_type>(x.size());
#if _DEBUG
    // Проверяем, что массив ненулевой
    assert(N != 0);
    // Проверяем, что размерности всех векторов равны
    assert(N == y.size());
    assert(N == z.size());
    assert(N == ux.size());
    assert(N == uy.size());
    assert(N == uz.size());
#endif
    return N;
}

auto mf::Particles::isAlive(size_type i) const -> bool { return is_alive[i] != 0; }

void mf::Particles::killParticle(size_type i) { is_alive[i] = 1; }

auto mf::Particles::sigma([[maybe_unused]] value_type gij) const noexcept -> value_type {
    return std::numbers::pi * molecule_size * molecule_size;
}

void mf::Particles::addParticles(size_type n_particles) {
    const size_type n_old = getNParticles();
    const size_type n_new = n_old + n_particles;

    x.resize(n_new);
    y.resize(n_new);
    z.resize(n_new);
    ux.resize(n_new);
    uy.resize(n_new);
    uz.resize(n_new);
    cell_id.resize(n_new);
    is_alive.resize(n_new, 1);
}
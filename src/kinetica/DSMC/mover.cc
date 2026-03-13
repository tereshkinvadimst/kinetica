#include "kinetica/DSMC/mover.hh"

#include <algorithm>
#include <execution>

void mf::mover(Particles& particles, std::span<const Particles::size_type> particles_id, Particles::value_type dt) noexcept {
    std::for_each(std::execution::unseq, particles_id.begin(), particles_id.end(), [dt, &particles](Particles::size_type id) {
        if (!particles.isAlive(id)) return;
        particles.x[id] += particles.ux[id] * dt;
        particles.y[id] += particles.uy[id] * dt;
        particles.z[id] += particles.uz[id] * dt;
    });
}
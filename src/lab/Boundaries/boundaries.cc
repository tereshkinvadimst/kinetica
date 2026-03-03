#include "lab/Boundaries/boundaries.hh"

void mf::applyPeriodic(Particles& particles, Box box, bool px, bool py, bool pz) noexcept {
    if (px) {
        applyPeriodic(particles.x, box.x0, box.Lx);
    }
    if (py) {
        applyPeriodic(particles.y, box.y0, box.Ly);
    }
    if (pz) {
        applyPeriodic(particles.z, box.z0, box.Lz);
    }
}

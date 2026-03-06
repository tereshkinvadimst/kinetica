#include "kinetica/Boundaries/boundaries.hh"

#include <random>

#include "kinetica/constants.hh"

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

auto mf::scatterDiffuse(char axis, int sign, double Tw, double m, random& rng) -> std::tuple<double, double, double> {
    // axis: 'x', 'y', 'z'
    // sign: +1 если нормаль направлена в +ось, -1 если в -ось
    const double RTw      = K_B * Tw / m;

    double       v_normal = Rayleigh(std::sqrt(RTw))(rng);  // >0
    double       v_t1     = std::normal_distribution<double>(0., std::sqrt(RTw))(rng);
    double       v_t2     = std::normal_distribution<double>(0., std::sqrt(RTw))(rng);

    double       vx{};
    double       vy{};
    double       vz{};

    switch (axis) {
        case 'x':
            vx = sign * v_normal;
            vy = v_t1;
            vz = v_t2;
            break;
        case 'y':
            vy = sign * v_normal;
            vx = v_t1;
            vz = v_t2;
            break;
        case 'z':
            vz = sign * v_normal;
            vx = v_t1;
            vy = v_t2;
            break;
        default:
            throw std::runtime_error("Invalid axis");
    }
    return {vx, vy, vz};
}

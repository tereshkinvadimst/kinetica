#include "lab/Boundaries/boundaries.hh"

#include <random>

#include "lab/Random/generators.hh"

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

mf::DiffuseWallBoundary::DiffuseWallBoundary(Box wall_position, value_type Tw) : wall_position_(wall_position), Tw_(Tw) {}

auto mf::DiffuseWallBoundary::scatter(char axis, int sign, value_type m, xoshiro256& rng) const
    -> std::tuple<value_type, value_type, value_type> {
    // axis: 'x', 'y', 'z'
    // sign: +1 если нормаль направлена в +ось, -1 если в -ось
    const value_type RTw      = K_B * Tw_ / m;

    value_type       v_normal = Rayleigh(std::sqrt(RTw))(rng);  // >0
    value_type       v_t1     = std::normal_distribution<value_type>(0., RTw)(rng);
    value_type       v_t2     = std::normal_distribution<value_type>(0., RTw)(rng);

    value_type       vx{};
    value_type       vy{};
    value_type       vz{};

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

void mf::DiffuseWallBoundary::killParticlesInBox(Particles& particles) const {
    for (size_type p{}; p < particles.getNParticles(); ++p) {
        if (wall_position_.contains(particles.x[p], particles.y[p], particles.z[p])) particles.killParticle(p);
    }
}

void mf::DiffuseWallBoundary::operator()(Particles& particles, value_type dt, xoshiro256& rng) const {
    const double x_min = wall_position_.x0;
    const double x_max = wall_position_.x0 + wall_position_.Lx;
    const double y_min = wall_position_.y0;
    const double y_max = wall_position_.y0 + wall_position_.Ly;
    const double z_min = wall_position_.z0;
    const double z_max = wall_position_.z0 + wall_position_.Lz;

    for (size_type p = 0; p < particles.getNParticles(); ++p) {
        const auto new_x = particles.x[p] + particles.ux[p] * dt;
        const auto new_y = particles.y[p] + particles.uy[p] * dt;
        const auto new_z = particles.z[p] + particles.uz[p] * dt;

        if (new_x <= x_min) {         // левая стенка
            const auto [vx, vy, vz] = scatter('x', +1, particles.m, rng);
            particles.ux[p]         = vx;
            particles.uy[p]         = vy;
            particles.uz[p]         = vz;
        } else if (new_x >= x_max) {  // правая стенка
            const auto [vx, vy, vz] = scatter('x', -1, particles.m, rng);
            particles.ux[p]         = vx;
            particles.uy[p]         = vy;
            particles.uz[p]         = vz;
        }

        if (new_y <= y_min) {         // нижняя стенка
            const auto [vx, vy, vz] = scatter('y', +1, particles.m, rng);
            particles.ux[p]         = vx;
            particles.uy[p]         = vy;
            particles.uz[p]         = vz;
        } else if (new_y >= y_max) {  // верхняя стенка
            const auto [vx, vy, vz] = scatter('y', -1, particles.m, rng);
            particles.ux[p]         = vx;
            particles.uy[p]         = vy;
            particles.uz[p]         = vz;
        }

        if (new_z <= z_min) {         // передняя стенка
            const auto [vx, vy, vz] = scatter('z', +1, particles.m, rng);
            particles.ux[p]         = vx;
            particles.uy[p]         = vy;
            particles.uz[p]         = vz;
        } else if (new_z >= z_max) {  // задняя стенка
            const auto [vx, vy, vz] = scatter('z', -1, particles.m, rng);
            particles.ux[p]         = vx;
            particles.uy[p]         = vy;
            particles.uz[p]         = vz;
        }
    }
}
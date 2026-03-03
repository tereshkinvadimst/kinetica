#include "lab/Domain/domain.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

#include "lab/Boundaries/boundaries.hh"
#include "lab/Random/xoshiro256.hh"
#include "lab/constants.hh"
#include "lab/init/generate_positions.hh"
#include "lab/init/generate_velocities.hh"
#include "version.hh"

namespace {

[[nodiscard]] auto computeNParticles(double n_density, double V, double W, mf::xoshiro256& gen) -> std::size_t {
    const auto                             N_avg = n_density * V / W;
    std::uniform_real_distribution<double> uni_01(0., 1.);
    if (N_avg > 20) {
        auto Np = static_cast<size_t>(N_avg);
        if (uni_01(gen) < (N_avg - static_cast<double>(Np))) ++Np;
        return Np;
    }
    std::poisson_distribution<std::size_t> poisson(N_avg);

    return poisson(gen);
}

[[nodiscard]] auto computeVelocityScaleFactor(double T, double m) noexcept -> double {
    using mf::K_B;
    return K_B * T / m;
}

}  // namespace

mf::Domain::Domain(Box domain_box, value_type m, value_type W, value_type molecule_size)
    : gen_{std::random_device{}()}
    , domain_box_(domain_box)
    , particles_{}
    , n_cells_x_{}
    , n_cells_y_{}
    , n_cells_z_{}
    , hx_{}
    , hy_{}
    , hz_{}
    , avg_sigma_g_max_{}
    , cells_{}
    , sigma_g_max_{}
    , cell_list_{} {
    particles_.m             = m;
    particles_.W             = W;
    particles_.molecule_size = molecule_size;
}

void mf::Domain::generateParticles(value_type n_density, value_type T) {
    const auto velocity_scale_factor = computeVelocityScaleFactor(T, particles_.m);
    const auto Np                    = computeNParticles(n_density, domain_box_.volume(), particles_.W, gen_);
    particles_                       = Particles(Np, particles_.m, particles_.W, particles_.molecule_size);
    // Генерация положений
    generateUniformPositions(particles_.x, domain_box_.x0, domain_box_.Lx, gen_);
    generateUniformPositions(particles_.y, domain_box_.y0, domain_box_.Ly, gen_);
    generateUniformPositions(particles_.z, domain_box_.z0, domain_box_.Lz, gen_);
    // Генерация скоростей
    generateMaxwellVelocities(particles_.ux, velocity_scale_factor, gen_);
    generateMaxwellVelocities(particles_.uy, velocity_scale_factor, gen_);
    generateMaxwellVelocities(particles_.uz, velocity_scale_factor, gen_);

    const auto max_velocity = std::max_element(particles_.ux.begin(), particles_.ux.end());
    avg_sigma_g_max_        = particles_.sigma(2. * (*max_velocity)) * (*max_velocity);
}

void mf::Domain::generateMesh(value_type hx, value_type hy, value_type hz) {
    const auto start = std::chrono::high_resolution_clock::now();

    hx_              = hx;
    hy_              = hy;
    hz_              = hz;
    // Вычисляем число ячеек
    n_cells_x_ = static_cast<size_type>(std::ceil(domain_box_.Lx / hx));
    n_cells_y_ = static_cast<size_type>(std::ceil(domain_box_.Ly / hy));
    n_cells_z_ = static_cast<size_type>(std::ceil(domain_box_.Lz / hz));
    cells_     = std::vector<Box>(n_cells_x_ * n_cells_y_ * n_cells_z_);
    std::cout << "n cells = " << cells_.size() << '\n';
    sigma_g_max_ = std::vector<value_type>(cells_.size(), avg_sigma_g_max_);
    // Заполняем ячейки
    for (size_type k = 0; k < n_cells_z_; ++k) {
        for (size_type j = 0; j < n_cells_y_; ++j) {
            for (size_type i = 0; i < n_cells_x_; ++i) {
                const size_type id   = i + n_cells_x_ * (j + n_cells_y_ * k);
                const double    xmin = domain_box_.x0 + i * hx;
                const double    ymin = domain_box_.y0 + j * hy;
                const double    zmin = domain_box_.z0 + k * hz;
                const double    xmax = std::min(xmin + hx, domain_box_.x0 + domain_box_.Lx);
                const double    ymax = std::min(ymin + hy, domain_box_.y0 + domain_box_.Ly);
                const double    zmax = std::min(zmin + hz, domain_box_.z0 + domain_box_.Lz);
                cells_[id].x0        = xmin;
                cells_[id].y0        = ymin;
                cells_[id].z0        = zmin;
                cells_[id].Lx        = xmax;
                cells_[id].Ly        = ymax;
                cells_[id].Lz        = zmax;
            }
        }
    }

    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Generate mesh " << duration_cast<std::chrono::milliseconds>(stop - start) << " ms\n";
}

void mf::Domain::makeCellList() {
    cell_list_ = CellList(cells_.size(), particles_.getNParticles());
    updateCellList();
}

void mf::Domain::updateCellList() {
    const auto start = std::chrono::high_resolution_clock::now();

    // Определяем, в каких ячейках находятся частицы
    for (size_type p{}; p < particles_.getNParticles(); ++p) {
        particles_.cell_id[p] = cellIndex(particles_.x[p], particles_.y[p], particles_.z[p]);
    }
    // Конструируем список ячеек
    cell_list_.build(particles_.cell_id);

    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Update cell list " << duration_cast<std::chrono::milliseconds>(stop - start) << " ms\n";
}

auto mf::Domain::cellIndex(double x, double y, double z) -> size_type const {
    const double rx = (x - domain_box_.x0) / hx_;
    const double ry = (y - domain_box_.y0) / hy_;
    const double rz = (z - domain_box_.z0) / hz_;

    size_type    i  = static_cast<std::size_t>(rx);
    size_type    j  = static_cast<std::size_t>(ry);
    size_type    k  = static_cast<std::size_t>(rz);

    // защита от x == xmax и от накопленных ошибок double
    if (i >= n_cells_x_) i = n_cells_x_ - 1;
    if (j >= n_cells_y_) j = n_cells_y_ - 1;
    if (k >= n_cells_z_) k = n_cells_z_ - 1;

    return i + n_cells_x_ * (j + n_cells_y_ * k);
}

void mf::Domain::applyPeriodicBoundaries(bool px, bool py, bool pz) { applyPeriodic(particles_, domain_box_, px, py, pz); }

void mf::Domain::moveParticles(value_type dt) {
    const auto start = std::chrono::high_resolution_clock::now();
    for (size_type p{}; p < particles_.getNParticles(); ++p) {
        particles_.x[p] += particles_.ux[p] * dt;
        particles_.y[p] += particles_.uy[p] * dt;
        particles_.z[p] += particles_.uz[p] * dt;
    }
    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Moves particles " << duration_cast<std::chrono::milliseconds>(stop - start) << " ms\n";
}

void mf::Domain::collideParticles(value_type dt) {
    const auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Start computing collisions...\n";

    std::uniform_real_distribution<value_type> u01(0., 1.);
    std::normal_distribution<double>           normal(0., 1.);
    for (size_type c{}; c < cells_.size(); ++c) {
        auto       particles_in_cell = cell_list_.getParticlesInCell(c);
        const auto Np                = particles_in_cell.size();
        if (Np < 2) continue;  // пропускаем ячейки с <2 частицами
        const auto collision_rate = sigma_g_max_[c] / cells_[c].volume();
        // Вычисляем число сталкивающихся частиц
        const auto                               N_avg      = Np;
        const auto                               N_coll_avg = 0.5 * Np * N_avg * particles_.W * collision_rate * dt;
        std::poisson_distribution<size_type>     poisson(N_coll_avg);
        const auto                               N_coll = poisson(gen_);
        std::uniform_int_distribution<size_type> rand_index(0, Np - 1);
        for (size_type collision{}; collision < N_coll; ++collision) {
            const auto i   = rand_index(gen_);
            auto       tmp = rand_index(gen_);
            while (tmp == i) {
                tmp = rand_index(gen_);
            }
            const auto j = tmp;
            // Вычисляем относительную скорость
            const auto dux      = particles_.ux[i] - particles_.ux[j];
            const auto duy      = particles_.uy[i] - particles_.uy[j];
            const auto duz      = particles_.uz[i] - particles_.uz[j];
            const auto gij_2    = dux * dux + duy * duy + duz * duz;
            const auto gij      = std::sqrt(gij_2);
            const auto half_gij = 0.5 * gij;
            const auto sigma_g  = particles_.sigma(gij) * gij;
            if (sigma_g > sigma_g_max_[c]) sigma_g_max_[c] = sigma_g;
            if (u01(gen_) < sigma_g / sigma_g_max_[c]) {  // Принимаем столкновения с вероятностью
                /// Считаем скорость центра масс
                const auto ucm_x = 0.5 * (particles_.ux[i] + particles_.ux[j]);
                const auto ucm_y = 0.5 * (particles_.uy[i] + particles_.uy[j]);
                const auto ucm_z = 0.5 * (particles_.uz[i] + particles_.uz[j]);
                /// Выбираем случайное направление
                auto       nx   = normal(gen_);
                auto       ny   = normal(gen_);
                auto       nz   = normal(gen_);
                const auto norm = std::sqrt(nx * nx + ny * ny + nz * nz);
                nx /= norm;
                ny /= norm;
                nz /= norm;
                /// Обновляем скорости
                particles_.ux[i] = ucm_x + half_gij * nx;
                particles_.uy[i] = ucm_y + half_gij * ny;
                particles_.uz[i] = ucm_z + half_gij * nz;

                particles_.ux[j] = ucm_x - half_gij * nx;
                particles_.uy[j] = ucm_y - half_gij * ny;
                particles_.uz[j] = ucm_z - half_gij * nz;
            }
        }
    }

    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Collide particles " << duration_cast<std::chrono::milliseconds>(stop - start) << " ms\n";
}

void mf::Domain::saveXYZ(std::string file_name) const {
    std::ofstream fout(file_name);
    if (!fout) {
        throw std::runtime_error("writeXYZ: Cannot open file " + file_name);
    }

    const auto N = particles_.getNParticles();

    fout << N << "\n";

    const auto        now = std::chrono::system_clock::now();
    const std::time_t t   = std::chrono::system_clock::to_time_t(now);

    fout << "file generated by DSMC program ver." << KINETICA_VERSION_STRING << " "
         << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "\n";

    for (size_type i{}; i < N; ++i) {
        fout << particles_.m << " " << particles_.x[i] << " " << particles_.y[i] << " " << particles_.z[i] << " "
             << particles_.ux[i] << " " << particles_.uy[i] << " " << particles_.uz[i] << '\n';
    }

    fout.close();
}

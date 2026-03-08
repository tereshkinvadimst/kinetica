#include "kinetica/Domain/domain.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <ranges>

#include "kinetica/Boundaries/boundaries.hh"
#include "kinetica/DSMC/utils.hh"
#include "kinetica/Random/random.hh"
#include "kinetica/constants.hh"
#include "kinetica/init/generate_positions.hh"
#include "kinetica/init/generate_velocities.hh"
#include "version.hh"

mf::Domain::Domain(Box domain_box, value_type m, value_type W, value_type molecule_size, value_type scale_factor)
    : gen_{std::random_device{}()}
    , domain_box_(domain_box)
    , particles_{}
    , n_cells_x_{}
    , n_cells_y_{}
    , n_cells_z_{}
    , sigma_g_max_{}
    , cells_{}
    , flow_properties_{}
    , cell_list_{}
    , stats_{}
    , is_diffuse_walls_{false, false, false, false, false, false}
    , diffuse_wall_temperature_{}
    , lambda_{}
    , time_step_{}
    , cell_size_{}
    , scale_factor_{scale_factor}
    , max_velocity_{}
    , xprofiler_(&domain_box, &flow_properties_) {
    particles_.m             = m;
    particles_.W             = W;
    particles_.molecule_size = molecule_size;
    stats_.initStat("t_coll [ms]");
    stats_.initStat("t_clist [ms]");
    stats_.initStat("t_mp [ms]");
    stats_.initStat("t_b [ms]");
    stats_.initStat("Navg");
}

void mf::Domain::generateParticles(value_type n_density, value_type T) {
    const auto velocity_scale_factor = computeVelocityScaleFactor(T, particles_.m);
    const auto Np                    = computeNParticles(n_density, domain_box_.volume(), particles_.W, gen_);
    std::cout << "Np = " << Np << '\n';
    particles_ = Particles(Np, particles_.m, particles_.W, particles_.molecule_size);
    lambda_    = meanFreePath(n_density, particles_.molecule_size);
    // Генерация положений
    generateUniformPositions(particles_.x, domain_box_.x0, domain_box_.x0 + domain_box_.Lx, gen_);
    generateUniformPositions(particles_.y, domain_box_.y0, domain_box_.y0 + domain_box_.Ly, gen_);
    generateUniformPositions(particles_.z, domain_box_.z0, domain_box_.z0 + domain_box_.Lz, gen_);
    // Генерация скоростей
    generateMaxwellVelocity(particles_.ux, velocity_scale_factor, gen_);
    generateMaxwellVelocity(particles_.uy, velocity_scale_factor, gen_);
    generateMaxwellVelocity(particles_.uz, velocity_scale_factor, gen_);

    const auto max_velocity_x = std::max_element(particles_.ux.begin(), particles_.ux.end());
    const auto max_velocity_y = std::max_element(particles_.uy.begin(), particles_.uy.end());
    const auto max_velocity_z = std::max_element(particles_.uz.begin(), particles_.uz.end());
    max_velocity_             = std::sqrt((*max_velocity_x) * (*max_velocity_x) + (*max_velocity_y) * (*max_velocity_y) +
                              (*max_velocity_z) * (*max_velocity_z));
    time_step_                = computeTimeStep(lambda_, max_velocity_, scale_factor_);
    sigma_g_max_              = particles_.sigma(2. * max_velocity_) * max_velocity_;
}

void mf::Domain::generateMesh() {
    cell_size_ = computeCellSize(lambda_, scale_factor_);
    // Вычисляем число ячеек
    n_cells_x_ = static_cast<size_type>(domain_box_.Lx / cell_size_);
    n_cells_y_ = static_cast<size_type>(domain_box_.Ly / cell_size_);
    n_cells_z_ = static_cast<size_type>(domain_box_.Lz / cell_size_);
    cells_     = std::vector<Box>(n_cells_x_ * n_cells_y_ * n_cells_z_);
    std::cout << "N cells x, y, z " << n_cells_x_ << '\t' << n_cells_y_ << '\t' << n_cells_z_ << '\n';
    flow_properties_ = FlowProperties(cells_.size());
    // Заполняем ячейки
    for (size_type k = 0; k < n_cells_z_; ++k) {
        for (size_type j = 0; j < n_cells_y_; ++j) {
            for (size_type i = 0; i < n_cells_x_; ++i) {
                const size_type  id = i + n_cells_x_ * (j + n_cells_y_ * k);
                const value_type x  = domain_box_.x0 + i * cell_size_;
                const value_type y  = domain_box_.y0 + j * cell_size_;
                const value_type z  = domain_box_.z0 + k * cell_size_;
                cells_[id].x0       = x;
                cells_[id].y0       = y;
                cells_[id].z0       = z;
                cells_[id].Lx       = cell_size_;
                cells_[id].Ly       = cell_size_;
                cells_[id].Lz       = cell_size_;
            }
        }
    }
}

void mf::Domain::makeCellList() {
    cell_list_ = CellList(cells_.size(), particles_.getNParticles());
    updateCellList();
}

void mf::Domain::updateCellList() {
    const auto start = std::chrono::high_resolution_clock::now();

    // Определяем, в каких ячейках находятся частицы
    for (size_type p{}; p < particles_.getNParticles(); ++p) {
        if (!particles_.isAlive(p)) continue;
        particles_.cell_id[p] = cellIndex(particles_.x[p], particles_.y[p], particles_.z[p]);
    }
    // Конструируем список ячеек
    cell_list_.build(particles_.cell_id);

    const auto stop = std::chrono::high_resolution_clock::now();
    stats_.addStat("t_clist [ms]", std::chrono::duration<value_type, std::milli>(stop - start).count());
}

auto mf::Domain::cellIndex(double x, double y, double z) -> size_type const {
    size_type i =
        std::min(static_cast<size_type>(std::floor((x - domain_box_.x0) / cell_size_)), static_cast<size_type>(n_cells_x_ - 1));
    size_type j =
        std::min(static_cast<size_type>(std::floor((y - domain_box_.y0) / cell_size_)), static_cast<size_type>(n_cells_y_ - 1));
    size_type k =
        std::min(static_cast<size_type>(std::floor((z - domain_box_.z0) / cell_size_)), static_cast<size_type>(n_cells_z_ - 1));
    i = std::max<size_type>(0, i);
    j = std::max<size_type>(0, j);
    k = std::max<size_type>(0, k);
    return i + n_cells_x_ * (j + n_cells_y_ * k);
}

void mf::Domain::applyPeriodicBoundaries(bool px, bool py, bool pz) { applyPeriodic(particles_, domain_box_, px, py, pz); }

void mf::Domain::moveParticles() {
    const auto start = std::chrono::high_resolution_clock::now();

    for (size_type p = 0; p < particles_.getNParticles(); ++p) {
        if (!particles_.isAlive(p)) continue;

        value_type x         = particles_.x[p];
        value_type y         = particles_.y[p];
        value_type z         = particles_.z[p];

        value_type ux        = particles_.ux[p];
        value_type uy        = particles_.uy[p];
        value_type uz        = particles_.uz[p];

        value_type new_x     = x + ux * time_step_;
        value_type new_y     = y + uy * time_step_;
        value_type new_z     = z + uz * time_step_;

        bool       scattered = false;
        value_type t_hit     = time_step_;

        char       axis      = 0;
        int        sign      = 0;
        value_type Tw        = 0;

        // ---------- X walls ----------
        if (ux != 0) {
            const value_type xw0 = domain_box_.x0;
            const value_type xw1 = domain_box_.x0 + domain_box_.Lx;

            if (new_x < xw0 && is_diffuse_walls_[0]) {
                const value_type t = (xw0 - x) / ux;
                if (t >= 0 && t <= time_step_) {
                    t_hit     = t;
                    axis      = 'x';
                    sign      = +1;
                    Tw        = diffuse_wall_temperature_[0];
                    scattered = true;
                }
            } else if (new_x > xw1 && is_diffuse_walls_[1]) {
                const value_type t = (xw1 - x) / ux;
                if (t >= 0 && t <= time_step_) {
                    t_hit     = t;
                    axis      = 'x';
                    sign      = -1;
                    Tw        = diffuse_wall_temperature_[1];
                    scattered = true;
                }
            }
        }

        // ---------- Y walls ----------
        if (!scattered && uy != 0) {
            const value_type yw0 = domain_box_.y0;
            const value_type yw1 = domain_box_.y0 + domain_box_.Ly;

            if (new_y < yw0 && is_diffuse_walls_[2]) {
                const value_type t = (yw0 - y) / uy;
                if (t >= 0 && t <= time_step_) {
                    t_hit     = t;
                    axis      = 'y';
                    sign      = +1;
                    Tw        = diffuse_wall_temperature_[2];
                    scattered = true;
                }
            } else if (new_y > yw1 && is_diffuse_walls_[3]) {
                const value_type t = (yw1 - y) / uy;
                if (t >= 0 && t <= time_step_) {
                    t_hit     = t;
                    axis      = 'y';
                    sign      = -1;
                    Tw        = diffuse_wall_temperature_[3];
                    scattered = true;
                }
            }
        }

        // ---------- Z walls ----------
        if (!scattered && uz != 0) {
            const value_type zw0 = domain_box_.z0;
            const value_type zw1 = domain_box_.z0 + domain_box_.Lz;

            if (new_z < zw0 && is_diffuse_walls_[4]) {
                const value_type t = (zw0 - z) / uz;
                if (t >= 0 && t <= time_step_) {
                    t_hit     = t;
                    axis      = 'z';
                    sign      = +1;
                    Tw        = diffuse_wall_temperature_[4];
                    scattered = true;
                }
            } else if (new_z > zw1 && is_diffuse_walls_[5]) {
                const value_type t = (zw1 - z) / uz;
                if (t >= 0 && t <= time_step_) {
                    t_hit     = t;
                    axis      = 'z';
                    sign      = -1;
                    Tw        = diffuse_wall_temperature_[5];
                    scattered = true;
                }
            }
        }

        if (!scattered) {
            particles_.x[p] = new_x;
            particles_.y[p] = new_y;
            particles_.z[p] = new_z;

        } else {
            // точка удара
            x += ux * t_hit;
            y += uy * t_hit;
            z += uz * t_hit;

            std::tie(ux, uy, uz) = scatterDiffuse(axis, sign, Tw, particles_.m, gen_);

            const value_type dt2 = time_step_ - t_hit;

            x += ux * dt2;
            y += uy * dt2;
            z += uz * dt2;

            particles_.x[p]  = x;
            particles_.y[p]  = y;
            particles_.z[p]  = z;

            particles_.ux[p] = ux;
            particles_.uy[p] = uy;
            particles_.uz[p] = uz;
        }
    }

    const auto stop = std::chrono::high_resolution_clock::now();
    stats_.addStat("t_mp [ms]", std::chrono::duration<value_type, std::milli>(stop - start).count());
}

void mf::Domain::collideParticles() {
    const auto                                 start = std::chrono::high_resolution_clock::now();

    std::uniform_real_distribution<value_type> u01(0., 1.);
    std::normal_distribution<double>           normal(0., 1.);
    value_type                                 Navg = {};
    for (size_type c{}; c < cells_.size(); ++c) {
        auto       particles_in_cell = cell_list_.getParticlesInCell(c);
        const auto Np                = particles_in_cell.size();
        Navg += Np;
        if (Np < 2) continue;  // пропускаем ячейки с <2 частицами
        const auto collision_rate = sigma_g_max_ / cells_[c].volume();
        // Вычисляем число сталкивающихся частиц
        const auto                               N_avg      = Np;
        const auto                               N_coll_avg = 0.5 * Np * N_avg * particles_.W * collision_rate * time_step_;
        std::poisson_distribution<size_type>     poisson(N_coll_avg);
        const auto                               N_coll = poisson(gen_);
        std::uniform_int_distribution<size_type> rand_index(0, Np - 1);
        for (size_type collision{}; collision < N_coll; ++collision) {
            const auto i   = particles_in_cell[rand_index(gen_)];
            auto       tmp = particles_in_cell[rand_index(gen_)];
            while (tmp == i) {
                tmp = particles_in_cell[rand_index(gen_)];
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
            if (sigma_g > sigma_g_max_) sigma_g_max_ = sigma_g;
            if (u01(gen_) < sigma_g / sigma_g_max_) {  // Принимаем столкновения с вероятностью
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
    stats_.addStat("t_coll [ms]", std::chrono::duration<value_type, std::milli>(stop - start).count());
    stats_.addStat("Navg", Navg / (n_cells_x_ * n_cells_y_ * n_cells_z_));
}

void mf::Domain::saveXYZ(std::string file_name) const {
    std::ofstream fout(file_name);
    if (!fout) {
        throw std::runtime_error("writeXYZ: Cannot open file " + file_name);
    }

    const auto N = std::count(particles_.is_alive.begin(), particles_.is_alive.end(), 1);

    fout << N << "\n";

    const auto        now = std::chrono::system_clock::now();
    const std::time_t t   = std::chrono::system_clock::to_time_t(now);

    fout << "file generated by DSMC program ver." << KINETICA_VERSION_STRING << " "
         << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "\n";

    for (size_type i{}; i < particles_.getNParticles(); ++i) {
        if (!particles_.isAlive(i)) continue;
        fout << particles_.m << " " << particles_.x[i] << " " << particles_.y[i] << " " << particles_.z[i] << " "
             << particles_.ux[i] << " " << particles_.uy[i] << " " << particles_.uz[i] << '\n';
    }

    fout.close();
}

void mf::Domain::computeFlowProperties() {
    const value_type mass_over_kB = particles_.m / K_B;

    for (size_type c = 0; c < cells_.size(); ++c) {
        auto              particles_ids = cell_list_.getParticlesInCell(c);

        const std::size_t N             = particles_ids.size();

        if (N == 0) {
            flow_properties_.n_density[c] = 0.0;

            flow_properties_.ux[c]        = 0.0;
            flow_properties_.uy[c]        = 0.0;
            flow_properties_.uz[c]        = 0.0;

            flow_properties_.u2x[c]       = 0.0;
            flow_properties_.u2y[c]       = 0.0;
            flow_properties_.u2z[c]       = 0.0;

            flow_properties_.Ttrx[c]      = 0.0;
            flow_properties_.Ttry[c]      = 0.0;
            flow_properties_.Ttrz[c]      = 0.0;
            flow_properties_.Ttr[c]       = 0.0;

            continue;
        }

        value_type ux = 0.0, uy = 0.0, uz = 0.0;
        value_type u2x = 0.0, u2y = 0.0, u2z = 0.0;

        for (auto i : particles_ids) {
            const auto vx = particles_.ux[i];
            const auto vy = particles_.uy[i];
            const auto vz = particles_.uz[i];

            ux += vx;
            uy += vy;
            uz += vz;

            u2x += vx * vx;
            u2y += vy * vy;
            u2z += vz * vz;
        }

        const value_type invN         = 1.0 / static_cast<value_type>(N);

        const value_type ux_m         = ux * invN;
        const value_type uy_m         = uy * invN;
        const value_type uz_m         = uz * invN;

        const value_type u2x_m        = u2x * invN;
        const value_type u2y_m        = u2y * invN;
        const value_type u2z_m        = u2z * invN;

        flow_properties_.n_density[c] = N * particles_.W / cells_[c].volume();

        flow_properties_.ux[c]        = ux_m;
        flow_properties_.uy[c]        = uy_m;
        flow_properties_.uz[c]        = uz_m;

        flow_properties_.u2x[c]       = u2x_m;
        flow_properties_.u2y[c]       = u2y_m;
        flow_properties_.u2z[c]       = u2z_m;

        // ВАЖНО: дисперсия
        const value_type cx2     = u2x_m - ux_m * ux_m;
        const value_type cy2     = u2y_m - uy_m * uy_m;
        const value_type cz2     = u2z_m - uz_m * uz_m;

        flow_properties_.Ttrx[c] = mass_over_kB * cx2;
        flow_properties_.Ttry[c] = mass_over_kB * cy2;
        flow_properties_.Ttrz[c] = mass_over_kB * cz2;

        flow_properties_.Ttr[c]  = (flow_properties_.Ttrx[c] + flow_properties_.Ttry[c] + flow_properties_.Ttrz[c]) / 3.0;
    }
}

void mf::Domain::printStatsHeader() { stats_.printHeader(); }

void mf::Domain::printStats(value_type time) { stats_.printStats(time); }

void mf::Domain::setDiffuseWall(size_type side, value_type Tw) {
    is_diffuse_walls_[side]         = true;
    diffuse_wall_temperature_[side] = Tw;
}

void mf::Domain::writeVTU(std::string file_name) const {
    std::ofstream out(file_name);
    if (!out) throw std::runtime_error("writeVTU: cannot open file " + file_name);

    const std::size_t nx      = n_cells_x_;
    const std::size_t ny      = n_cells_y_;
    const std::size_t nz      = n_cells_z_;

    const std::size_t nPoints = (nx + 1) * (ny + 1) * (nz + 1);

    const std::size_t nCells  = cells_.size();

    auto              pointId = [&](std::size_t i, std::size_t j, std::size_t k) { return i + (nx + 1) * (j + (ny + 1) * k); };

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << nPoints << "\" NumberOfCells=\"" << nCells << "\">\n";

    // ------------------------------------------------------------
    // Points
    // ------------------------------------------------------------
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (std::size_t k = 0; k <= nz; ++k) {
        double z = domain_box_.z0 + k * cell_size_;
        if (z > domain_box_.z0 + domain_box_.Lz) z = domain_box_.z0 + domain_box_.Lz;

        for (std::size_t j = 0; j <= ny; ++j) {
            double y = domain_box_.y0 + j * cell_size_;
            if (y > domain_box_.y0 + domain_box_.Ly) y = domain_box_.y0 + domain_box_.Ly;

            for (std::size_t i = 0; i <= nx; ++i) {
                double x = domain_box_.x0 + i * cell_size_;
                if (x > domain_box_.x0 + domain_box_.Lx) x = domain_box_.x0 + domain_box_.Lx;

                out << x << " " << y << " " << z << "\n";
            }
        }
    }

    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // ------------------------------------------------------------
    // Cells
    // ------------------------------------------------------------
    out << "      <Cells>\n";

    // connectivity
    out << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t p0 = pointId(i, j, k);
                const std::size_t p1 = pointId(i + 1, j, k);
                const std::size_t p2 = pointId(i + 1, j + 1, k);
                const std::size_t p3 = pointId(i, j + 1, k);
                const std::size_t p4 = pointId(i, j, k + 1);
                const std::size_t p5 = pointId(i + 1, j, k + 1);
                const std::size_t p6 = pointId(i + 1, j + 1, k + 1);
                const std::size_t p7 = pointId(i, j + 1, k + 1);

                out << p0 << " " << p1 << " " << p2 << " " << p3 << " " << p4 << " " << p5 << " " << p6 << " " << p7 << "\n";
            }
        }
    }

    out << "        </DataArray>\n";

    // offsets
    out << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";

    std::size_t offset = 0;
    for (std::size_t c = 0; c < nCells; ++c) {
        offset += 8;
        out << offset << "\n";
    }

    out << "        </DataArray>\n";

    // types (VTK_HEXAHEDRON = 12)
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < nCells; ++c) out << 12 << "\n";
    out << "        </DataArray>\n";

    out << "      </Cells>\n";

    // ------------------------------------------------------------
    // CellData
    // ------------------------------------------------------------
    out << "      <CellData>\n";

    auto writeScalar = [&](const char* name, const std::vector<double>& v) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (double x : v) out << x << "\n";
        out << "        </DataArray>\n";
    };

    auto writeVector =
        [&](const char* name, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz) {
            out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";

            for (std::size_t i = 0; i < vx.size(); ++i) out << vx[i] << " " << vy[i] << " " << vz[i] << "\n";

            out << "        </DataArray>\n";
        };

    writeScalar("n_density", flow_properties_.n_density);

    writeVector("u", flow_properties_.ux, flow_properties_.uy, flow_properties_.uz);

    writeVector("u2", flow_properties_.u2x, flow_properties_.u2y, flow_properties_.u2z);

    writeScalar("Ttrx", flow_properties_.Ttrx);
    writeScalar("Ttry", flow_properties_.Ttry);
    writeScalar("Ttrz", flow_properties_.Ttrz);
    writeScalar("Ttr", flow_properties_.Ttr);

    out << "      </CellData>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

auto mf::Domain::getTimeStep() const noexcept -> value_type { return time_step_; }

void mf::Domain::writeXProfile(std::string file_name) { xprofiler_(file_name, n_cells_x_, n_cells_y_, n_cells_z_, cell_size_); }

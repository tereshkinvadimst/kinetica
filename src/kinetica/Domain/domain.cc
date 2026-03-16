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
#include "kinetica/DSMC/collider.hh"
#include "kinetica/DSMC/mover.hh"
#include "kinetica/DSMC/utils.hh"
#include "kinetica/Random/random.hh"
#include "kinetica/constants.hh"
#include "kinetica/init/generate_positions.hh"
#include "kinetica/init/generate_velocities.hh"
#include "version.hh"

mf::Domain::Domain(
    Box domain_box, value_type m, value_type W, value_type molecule_size, value_type scale_factor, value_type time_scale_factor)
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
    , lambda_{std::numeric_limits<value_type>::max()}
    , time_step_{std::numeric_limits<value_type>::max()}
    , cell_size_{}
    , scale_factor_{scale_factor}
    , max_velocity_{}
    , xprofiler_(&domain_box_, &flow_properties_)
    , walls_{}
    , time_scale_factor_{time_scale_factor} {
    particles_.m             = m;
    particles_.W             = W;
    particles_.molecule_size = molecule_size;
    stats_.initStat("t_coll [ms]");
    stats_.initStat("t_clist [ms]");
    stats_.initStat("t_mp [ms]");
    stats_.initStat("t_b [ms]");
    stats_.initStat("Navg");
}

void mf::Domain::generateParticles(value_type n_density, value_type T, value_type x_min, value_type x_max) {
    const auto       velocity_scale_factor = computeVelocityScaleFactor(T, particles_.m);
    const value_type x0                    = std::max(x_min, domain_box_.x0);
    const value_type x_end                 = std::min(x_max, domain_box_.x0 + domain_box_.Lx);

    const value_type Lx                    = x_end - x0;
    const Box        local_box{
               .x0 = x0,
               .y0 = domain_box_.y0,
               .z0 = domain_box_.z0,
               .Lx = Lx,
               .Ly = domain_box_.Ly,
               .Lz = domain_box_.Lz,
    };
    const auto Np = computeNParticles(n_density, local_box.volume(), particles_.W, gen_);

    std::cout << "Added Np = " << Np << '\n';

    const size_type old_size = particles_.getNParticles();
    particles_.addParticles(Np);

    const value_type new_lambda = meanFreePath(n_density, particles_.molecule_size);

    if (new_lambda < lambda_) {
        lambda_ = new_lambda;
    }

    // Генерация положений
    generateUniformPositions(std::span{particles_.x.data() + old_size, Np}, x0, x0 + Lx, gen_);
    generateUniformPositions(
        std::span{particles_.y.data() + old_size, Np}, domain_box_.y0, domain_box_.y0 + domain_box_.Ly, gen_);
    generateUniformPositions(
        std::span{particles_.z.data() + old_size, Np}, domain_box_.z0, domain_box_.z0 + domain_box_.Lz, gen_);
    // Генерация скоростей
    generateMaxwellVelocity(std::span{particles_.ux.data() + old_size, Np}, velocity_scale_factor, gen_);
    generateMaxwellVelocity(std::span{particles_.uy.data() + old_size, Np}, velocity_scale_factor, gen_);
    generateMaxwellVelocity(std::span{particles_.uz.data() + old_size, Np}, velocity_scale_factor, gen_);

    computeMaxVelocity();
    const value_type new_time_step   = computeTimeStep(scale_factor_ * lambda_, max_velocity_, time_scale_factor_);
    const value_type new_sigma_g_max = particles_.sigma(2. * max_velocity_) * max_velocity_;
    if (new_time_step < time_step_) {
        time_step_ = new_time_step;
    }
    if (new_sigma_g_max > sigma_g_max_) {
        sigma_g_max_ = new_sigma_g_max;
    }
}

void mf::Domain::generateMesh() {
    cell_size_ = computeCellSize(lambda_, scale_factor_);
    std::cout << "cell size is " << cell_size_ << '\n';
    // Вычисляем число ячеек
    n_cells_x_ = static_cast<size_type>(std::ceil(domain_box_.Lx / cell_size_));
    n_cells_y_ = static_cast<size_type>(std::ceil(domain_box_.Ly / cell_size_));
    n_cells_z_ = static_cast<size_type>(std::ceil(domain_box_.Lz / cell_size_));
    std::cout << "N cells x, y, z " << n_cells_x_ << '\t' << n_cells_y_ << '\t' << n_cells_z_ << '\n';
    cells_           = std::vector<Box>(n_cells_x_ * n_cells_y_ * n_cells_z_);
    flow_properties_ = FlowProperties(cells_.size());
    // Заполняем ячейки
    for (size_type k = 0; k < n_cells_z_; ++k) {
        const value_type z = domain_box_.z0 + k * cell_size_;
        for (size_type j = 0; j < n_cells_y_; ++j) {
            const value_type y = domain_box_.y0 + j * cell_size_;
            for (size_type i = 0; i < n_cells_x_; ++i) {
                const size_type  id = i + n_cells_x_ * (j + n_cells_y_ * k);
                const value_type x  = domain_box_.x0 + i * cell_size_;
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
    const value_type rx = (x - domain_box_.x0) / cell_size_;
    const value_type ry = (y - domain_box_.y0) / cell_size_;
    const value_type rz = (z - domain_box_.z0) / cell_size_;

    size_type        i  = std::clamp(static_cast<size_type>(rx), size_type{}, n_cells_x_ - 1);
    size_type        j  = std::clamp(static_cast<size_type>(ry), size_type{}, n_cells_y_ - 1);
    size_type        k  = std::clamp(static_cast<size_type>(rz), size_type{}, n_cells_z_ - 1);

    i                   = std::max<size_type>(0, i);
    j                   = std::max<size_type>(0, j);
    k                   = std::max<size_type>(0, k);
    return i + n_cells_x_ * (j + n_cells_y_ * k);
}

void mf::Domain::applyPeriodicBoundaries(bool px, bool py, bool pz) { applyPeriodic(particles_, domain_box_, px, py, pz); }

void mf::Domain::moveParticles() {
    const auto start = std::chrono::high_resolution_clock::now();
    // Вычисляем новый шаг по времени
    computeMaxVelocity();
    // Вычисляем длину свободного пробега
    const value_type new_time_step   = computeTimeStep(scale_factor_ * lambda_, max_velocity_, time_scale_factor_);
    const value_type new_sigma_g_max = particles_.sigma(2. * max_velocity_) * max_velocity_;
    if (new_time_step < time_step_) {
        time_step_ = new_time_step;
    }
    if (new_sigma_g_max > sigma_g_max_) {
        sigma_g_max_ = new_sigma_g_max;
    }
    // Перемещаем частицы
    for (size_type c{}; c < cells_.size(); ++c) {
        auto particles_in_cell = cell_list_.getParticlesInCell(c);
        // Если стенок нет!
        if (walls_.empty()) {
            mover(particles_, particles_in_cell, time_step_);
            continue;
        }
        bool has_walls   = false;
        auto time_remain = time_step_;
        for (const auto id : particles_in_cell) {   // Цикл по всем частицам в ячейке
            for (auto& wall : walls_) {             // Цикл по всем стенкам
                if (wall->intersects(cells_[c])) {  // Проверяем, что стенка пересекает ячейку
                    wall->collide(particles_, id, time_remain);
                }
            }
            // Перемещаем частицу за оставшееся время
            particles_.x[id] += particles_.ux[id] * time_remain;
            particles_.y[id] += particles_.uy[id] * time_remain;
            particles_.z[id] += particles_.uz[id] * time_remain;
        }

        // Если стенок в ячейках нет, то просто перемещаем частицы
        if (!has_walls) mover(particles_, particles_in_cell, time_step_);
    }

    // Перемещаем стенки
    for (auto& wall : walls_) {
        wall->move(time_step_);
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
        auto particles_in_cell = cell_list_.getParticlesInCell(c);
        Navg += particles_in_cell.size();
        sigma_g_max_ = collider(particles_, particles_in_cell, cells_[c].volume(), time_step_, sigma_g_max_, gen_);
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
            flow_properties_.n_particles[c] = 0.0;

            flow_properties_.n_density[c]   = 0.0;

            flow_properties_.ux[c]          = 0.0;
            flow_properties_.uy[c]          = 0.0;
            flow_properties_.uz[c]          = 0.0;

            flow_properties_.u2x[c]         = 0.0;
            flow_properties_.u2y[c]         = 0.0;
            flow_properties_.u2z[c]         = 0.0;

            flow_properties_.Ttrx[c]        = 0.0;
            flow_properties_.Ttry[c]        = 0.0;
            flow_properties_.Ttrz[c]        = 0.0;
            flow_properties_.Ttr[c]         = 0.0;

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

        const value_type invN           = 1.0 / static_cast<value_type>(N);

        const value_type ux_m           = ux * invN;
        const value_type uy_m           = uy * invN;
        const value_type uz_m           = uz * invN;

        const value_type u2x_m          = u2x * invN;
        const value_type u2y_m          = u2y * invN;
        const value_type u2z_m          = u2z * invN;

        flow_properties_.n_particles[c] = N;

        flow_properties_.n_density[c]   = N * particles_.W / cells_[c].volume();

        flow_properties_.ux[c]          = ux_m;
        flow_properties_.uy[c]          = uy_m;
        flow_properties_.uz[c]          = uz_m;

        flow_properties_.u2x[c]         = u2x_m;
        flow_properties_.u2y[c]         = u2y_m;
        flow_properties_.u2z[c]         = u2z_m;

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

void mf::Domain::addWall(std::shared_ptr<Wall> wall) { walls_.emplace_back(wall); }

void mf::Domain::computeMaxVelocity() {
    // Максимальный модуль скорости по каждой координате
    double max_ux = *std::max_element(
        particles_.ux.begin(), particles_.ux.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    double max_uy = *std::max_element(
        particles_.uy.begin(), particles_.uy.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    double max_uz = *std::max_element(
        particles_.uz.begin(), particles_.uz.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });

    // Максимальная скорость частицы в системе
    const auto particle_max_velocity = std::sqrt(max_ux * max_ux + max_uy * max_uy + max_uz * max_uz);

    value_type wall_max_velocity     = {};

    for (const auto& wall : walls_) {
        const auto& v     = wall->getVelocity();  // v — std::array<double,3> или Eigen::Vector3d
        double      vnorm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (vnorm > wall_max_velocity) wall_max_velocity = vnorm;
    }

    // --- 3. Общая максимальная скорость в системе ---
    max_velocity_ = std::max(particle_max_velocity, wall_max_velocity);
}
#include "kinetica/DSMC/wall.hh"

#include <algorithm>
#include <cmath>
#include <execution>
#include <limits>

#include "kinetica/constants.hh"

mf::Wall::Wall(Vector3 position, Vector3 normal, Vector2 size, Vector3 velocity, [[maybe_unused]] value_type Tw)
    : plane_(normal.normalized(), position), center_{position}, axis1_{}, axis2_{}, velocity_{velocity}, size_{size} {
    compute_plane_axes();
}

void mf::Wall::compute_plane_axes() {
    axis1_ = plane_.normal().unitOrthogonal();
    axis2_ = plane_.normal().cross(axis1_).normalized();
}

auto mf::Wall::timeTo(const Particles& particles, size_type i) const -> value_type {
    const Vector3    pos{particles.x[i], particles.y[i], particles.z[i]};
    const Vector3    vel{particles.ux[i], particles.uy[i], particles.uz[i]};

    const Vector3    v_rel = vel - velocity_;

    const value_type denom = plane_.normal().dot(v_rel);

    if (denom >= 0) return std::numeric_limits<value_type>::infinity();

    const value_type t = plane_.signedDistance(pos) / denom;

    if (t < 0) return std::numeric_limits<value_type>::infinity();

    const Vector3    pos_coll = pos + v_rel * t;

    const Vector3    d        = pos_coll - center_;

    const value_type u        = d.dot(axis1_);
    const value_type v        = d.dot(axis2_);

    if (std::abs(u) > 0.5 * size_.x() || std::abs(v) > 0.5 * size_.y()) return std::numeric_limits<value_type>::infinity();

    return t;
}

void mf::Wall::collide(Particles& particles, size_type id, value_type& time_remain) const {
    specularReflection(particles, id, time_remain);
}

void mf::Wall::specularReflection(Particles& particles, size_type id, value_type& time_remain) const {
    if (!particles.isAlive(id)) return;
    // Вычисляем время до столкновения
    const auto t_wall  = timeTo(particles, id);
    const auto dt_wall = time_remain - t_wall;
    // если стенка слишком далеко, то просто пропускаем
    if (t_wall > time_remain) {
        return;
    }
    // В другом случае перемещаем частицу к стенке
    particles.x[id] += particles.ux[id] * t_wall;
    particles.y[id] += particles.uy[id] * t_wall;
    particles.z[id] += particles.uz[id] * t_wall;

    // Вычисляем скорость после отражения
    const Vector3    v{particles.ux[id], particles.uy[id], particles.uz[id]};

    const Vector3    v_rel = v - velocity_;

    const value_type dot   = v_rel.dot(plane_.normal());

    const Vector3    v_new = v_rel - value_type{2} * dot * plane_.normal() + velocity_;
    // Записываем скорость после отражения
    particles.ux[id] = v_new.x();
    particles.uy[id] = v_new.y();
    particles.uz[id] = v_new.z();

    // Вычисляем оставшееся время после столкновения
    const auto dt_remain = time_remain - t_wall;
}

void mf::Wall::setVelocity(Vector3 velocity) { velocity_ = velocity; }

void mf::Wall::move(value_type dt) {
    center_ += velocity_ * dt;
    plane_ = Plane(plane_.normal(), center_);
}

auto mf::Wall::intersects(const Box& cell) const noexcept -> bool {
    const auto&   n        = plane_.normal();

    const Vector3 verts[8] = {
        {cell.x0, cell.y0, cell.z0},
        {cell.x0 + cell.Lx, cell.y0, cell.z0},
        {cell.x0, cell.y0 + cell.Ly, cell.z0},
        {cell.x0, cell.y0, cell.z0 + cell.Lz},
        {cell.x0 + cell.Lx, cell.y0 + cell.Ly, cell.z0},
        {cell.x0 + cell.Lx, cell.y0, cell.z0 + cell.Lz},
        {cell.x0, cell.y0 + cell.Ly, cell.z0 + cell.Lz},
        {cell.x0 + cell.Lx, cell.y0 + cell.Ly, cell.z0 + cell.Lz},
    };

    value_type min_d = std::numeric_limits<value_type>::max();
    value_type max_d = -std::numeric_limits<value_type>::max();

    for (const auto& p : verts) {
        const auto d = n.dot(p - center_);
        min_d        = std::min(min_d, d);
        max_d        = std::max(max_d, d);
    }

    if (min_d > 0 || max_d < 0) return false;

    const Vector3    cell_center(cell.cx(), cell.cy(), cell.cz());

    const Vector3    proj    = cell_center - n * n.dot(cell_center - center_);

    const Vector3    d       = proj - center_;

    const value_type u_coord = d.dot(axis1_);
    const value_type v_coord = d.dot(axis2_);

    return std::abs(u_coord) <= value_type{0.5} * size_.x() && std::abs(v_coord) <= value_type{0.5} * size_.y();
}

auto mf::Wall::getCenterPosition() const noexcept -> Vector3 { return center_; }

auto mf::Wall::getVelocity() const noexcept -> Vector3 { return velocity_; }

mf::DiffuseWall::DiffuseWall(Vector3 position, Vector3 normal, Vector2 size, Vector3 velocity, value_type Tw)
    : Wall(position, normal, size, velocity), Tw_(Tw) {}

void mf::DiffuseWall::collide(Particles& particles, size_type id, value_type& time_remain, random& gen) const {
    const auto t_wall = timeTo(particles, id);

    if (t_wall > time_remain) {
        return;
    }

    // движение до стенки
    particles.x[id] += particles.ux[id] * t_wall;
    particles.y[id] += particles.uy[id] * t_wall;
    particles.z[id] += particles.uz[id] * t_wall;

    time_remain = time_remain - t_wall;

    diffuseReflection(particles, id, gen);
}

void mf::DiffuseWall::diffuseReflection(Particles& particles, size_type id, random& gen) const {
    static std::normal_distribution<value_type>       normal_dist(0.0, 1.0);
    static std::uniform_real_distribution<value_type> u01(0.0, 1.0);

    const Vector3&                                    n   = plane_.normal();
    const Vector3&                                    t1  = axis1_;
    const Vector3&                                    t2  = axis2_;

    const value_type                                  RTw = K_B * Tw_ / particles.m;

    // тангенциальные компоненты
    const value_type vt1 = std::normal_distribution<double>(0., std::sqrt(RTw))(gen);
    const value_type vt2 = std::normal_distribution<double>(0., std::sqrt(RTw))(gen);

    // нормальная компонента (Rayleigh)
    const value_type vn = Rayleigh(std::sqrt(RTw))(gen);

    // относительная скорость
    const value_type vx = vn * n.x() + vt1 * t1.x() + vt2 * t2.x();

    const value_type vy = vn * n.y() + vt1 * t1.y() + vt2 * t2.y();

    const value_type vz = vn * n.z() + vt1 * t1.z() + vt2 * t2.z();

    // перевод в лабораторную систему
    particles.ux[id] = vx + velocity_.x();
    particles.uy[id] = vy + velocity_.y();
    particles.uz[id] = vz + velocity_.z();
}

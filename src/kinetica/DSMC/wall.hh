#ifndef MF_KINETICA_DSMC_WALL_H
#define MF_KINETICA_DSMC_WALL_H
#pragma once

#include <Eigen/Dense>
#include <span>

#include "kinetica/Particles/particles.hh"
#include "kinetica/Random/random.hh"

namespace mf {

class Wall {
   public:
    using value_type = Particles::value_type;
    using size_type  = Particles::size_type;
    using Vector3    = Eigen::Vector3<value_type>;
    using Vector2    = Eigen::Vector2<value_type>;
    using Plane      = Eigen::Hyperplane<value_type, 3>;

    Wall(Vector3 position, Vector3 normal, Vector2 size, Vector3 velocity = Vector3::Zero(), value_type Tw = {});
    virtual ~Wall() = default;

    void collide(Particles& particles, std::span<const size_type> particles_id, value_type dt) const;

    void setVelocity(Vector3 velocity);

    void move(value_type dt);

   protected:
    void compute_plane_axes();

    auto timeTo(const Particles& particles, size_type i) const -> value_type;

    void specularReflection(Particles& particles, std::span<const size_type> particles_id, value_type dt) const;

   protected:
    Plane   plane_;
    Vector3 center_;
    Vector3 axis1_;
    Vector3 axis2_;
    Vector3 velocity_;
    Vector2 size_;
};

class DiffuseWall final : public Wall {
   public:
    DiffuseWall(Vector3 position, Vector3 normal, Vector2 size, Vector3 velocity, value_type Tw);

    void collide(Particles& particles, std::span<const size_type> particles_id, value_type dt, random& gen) const;

   private:
    void diffuseReflection(Particles& particles, size_type id, random& gen) const;

   private:
    value_type Tw_;
};

}  // namespace mf

#endif  // MF_KINETICA_DSMC_WALL_H
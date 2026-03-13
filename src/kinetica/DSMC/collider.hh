#ifndef MF_KINETICA_DSMC_COLLIDER_H
#define MF_KINETICA_DSMC_COLLIDER_H
#pragma once
#include <span>

#include "kinetica/Particles/particles.hh"
#include "kinetica/Random/random.hh"

namespace mf {

[[nodiscard]] auto collider(Particles&                            particles,
                            std::span<const Particles::size_type> particles_id,
                            Particles::value_type                 volume,
                            Particles::value_type                 dt,
                            Particles::value_type                 sigma_g_max,
                            random&                               gen) -> Particles::value_type;

}

#endif  // MF_KINETICA_DSMC_COLLIDER_H
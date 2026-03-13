#ifndef MF_KINETICA_DSMC_MOVER_H
#define MF_KINETICA_DSMC_MOVER_H
#pragma once

#include <span>

#include "kinetica/Particles/particles.hh"

namespace mf {

void mover(Particles& particles, std::span<const Particles::size_type> particles_id, Particles::value_type dt) noexcept;


}  // namespace mf

#endif  // MF_KINETICA_DSMC_MOVER_H
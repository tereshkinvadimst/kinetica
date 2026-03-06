#ifndef MF_KINETICA_RANDOM_H
#define MF_KINETICA_RANDOM_H
#pragma once
#include <random>

#include "kinetica/Random/generators.hh"
#include "kinetica/Random/xoshiro256.hh"

namespace mf {

using random = xoshiro256;

}

#endif  // MF_KINETICA_RANDOM_H
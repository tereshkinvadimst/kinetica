#include "kinetica/Random/xoshiro256.hh"

mf::xoshiro256::xoshiro256(result_type seed) { stateFromSeed(seed); }

auto mf::xoshiro256::operator()() -> result_type {
    const result_type result = rotateLeft(state_[1] * 5, 7) * 9;

    const result_type t      = state_[1] << 17;

    state_[2] ^= state_[0];
    state_[3] ^= state_[1];
    state_[1] ^= state_[2];
    state_[0] ^= state_[3];

    state_[2] ^= t;
    state_[3] = rotateLeft(state_[3], 45);

    return result;
}

auto mf::xoshiro256::rotateLeft(const result_type x, int k) noexcept -> result_type { return (x << k) | (x >> (64 - k)); }

void mf::xoshiro256::stateFromSeed(result_type seed) {
    splitmix64 sm(seed);
    state_[0] = sm.next();
    state_[1] = sm.next();
    state_[2] = sm.next();
    state_[3] = sm.next();
}

mf::xoshiro256::splitmix64::splitmix64(result_type seed) : state_(seed) {}

auto mf::xoshiro256::splitmix64::next() noexcept -> result_type {
    result_type z = (state_ += 0x9e3779b97f4a7c15ULL);
    z             = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z             = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
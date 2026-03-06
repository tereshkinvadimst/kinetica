#ifndef MF_KINETICA_XOSHIRO256_H
#define MF_KINETICA_XOSHIRO256_H
#pragma once

#include <array>
#include <cstdint>
#include <limits>

namespace mf {

class xoshiro256 final {
   public:
    using result_type = std::uint64_t;

    explicit xoshiro256(result_type seed);
    // -----------------------------------------------
    // Запрет на копирование
    xoshiro256(const xoshiro256&)            = delete;
    xoshiro256& operator=(const xoshiro256&) = delete;
    // -----------------------------------------------

    /// Генерация случайного числа
    auto                         operator()() -> result_type;

    static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }
    static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }

   private:
    // Битовый циклических сдвиг влево
    static auto rotateLeft(const result_type x, int k) noexcept -> result_type;

    void        stateFromSeed(result_type seed);
    class splitmix64 {
       public:
        explicit splitmix64(result_type seed);

        auto next() noexcept -> result_type;

       private:
        result_type state_;
    };

   private:
    std::array<result_type, 4> state_;
};

}  // namespace mf

#endif  // MF_KINETICA_XOSHIRO256_H
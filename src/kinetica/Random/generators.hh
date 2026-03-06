#ifndef MF_KINETICA_GENERATORS_H
#define MF_KINETICA_GENERATORS_H
#include <cmath>
#include <random>

namespace mf {

class Rayleigh {
   public:
    explicit Rayleigh(double sigma) : sigma_(sigma), dist_(0., 1.) {}

    template <class RNG>
    double operator()(RNG& rng) {
        double U = dist_(rng);  // равномерное [0,1)
        return sigma_ * std::sqrt(-2.0 * std::log(U));
    }

   private:
    double                                 sigma_;
    std::uniform_real_distribution<double> dist_;
};
}  // namespace mf

#endif  // MF_KINETICA_GENERATORS_H
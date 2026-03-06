#ifndef MF_KINETICA_FLOW_PROPERTIES_H
#define MF_KINETICA_FLOW_PROPERTIES_H
#pragma once
#include <vector>

namespace mf {

struct FlowProperties final {
    using value_type = double;
    using size_type  = std::size_t;

    FlowProperties() = default;
    explicit FlowProperties(size_type N_cells);

    std::vector<value_type> n_density;
    std::vector<value_type> ux;
    std::vector<value_type> uy;
    std::vector<value_type> uz;
    std::vector<value_type> u2x;
    std::vector<value_type> u2y;
    std::vector<value_type> u2z;
    std::vector<value_type> Ttrx;
    std::vector<value_type> Ttry;
    std::vector<value_type> Ttrz;
    std::vector<value_type> Ttr;
};

}  // namespace mf

#endif  // MF_KINETICA_FLOW_PROPERTIES_H
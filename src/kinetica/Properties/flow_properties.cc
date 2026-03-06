#include "kinetica/Properties/flow_properties.hh"

mf::FlowProperties::FlowProperties(size_type N_cells)
    : n_density(N_cells)
    , ux(N_cells)
    , uy(N_cells)
    , uz(N_cells)
    , u2x(N_cells)
    , u2y(N_cells)
    , u2z(N_cells)
    , Ttrx(N_cells)
    , Ttry(N_cells)
    , Ttrz(N_cells)
    , Ttr(N_cells) {}

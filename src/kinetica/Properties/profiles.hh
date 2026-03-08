#ifndef MF_KINETICA_PROPERTIES_PROFILES_H
#define MF_KINETICA_PROPERTIES_PROFILES_H
#pragma once

#include <string>

#include "kinetica/Box/box.hh"
#include "kinetica/CellList/cell_list.hh"
#include "kinetica/Properties/flow_properties.hh"

namespace mf {

class XProfiler final {
   public:
    using value_type = FlowProperties::value_type;
    using size_type  = FlowProperties::size_type;

    XProfiler(Box* domain_box, FlowProperties* flow_properties);

    void operator()(std::string file_name, size_type n_cells_x, size_type n_cells_y, size_type n_cells_z, value_type cell_size);

   private:
    Box*            domain_box_;
    FlowProperties* flow_properties_;
};

}  // namespace mf

#endif  // MF_KINETICA_PROPERTIES_PROFILES_H
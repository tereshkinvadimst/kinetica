#ifndef MF_KINETICA_DOMAIN_H
#define MF_KINETICA_DOMAIN_H
#pragma once
#include <string>

#include "kinetica/Boundaries/boundaries.hh"
#include "kinetica/Box/box.hh"
#include "kinetica/CellList/cell_list.hh"
#include "kinetica/Particles/particles.hh"
#include "kinetica/Properties/flow_properties.hh"
#include "kinetica/Properties/stats.hh"
#include "kinetica/Random/random.hh"


namespace mf {

class Domain {
   public:
    using size_type  = std::size_t;
    using value_type = double;

    explicit Domain(Box domain_box, value_type m, value_type W, value_type molecule_size);

    void generateParticles(value_type n_density, value_type T);
    void generateMesh(value_type hx, value_type hy, value_type hz);
    void makeCellList();
    void updateCellList();
    void applyPeriodicBoundaries(bool px, bool py, bool pz);
    void moveParticles(value_type dt);
    void collideParticles(value_type dt);
    void saveXYZ(std::string file_name) const;
    void computeFlowProperties();
    void printStatsHeader();
    void printStats(value_type time);
    void setDiffuseWall(size_type side, value_type Tw);
    void writeVTU(std::string file_name) const;

   private:
    auto cellIndex(double x, double y, double z) -> size_type const;

   private:
    random                    gen_;
    Box                       domain_box_;
    Particles                 particles_;
    size_type                 n_cells_x_;
    size_type                 n_cells_y_;
    size_type                 n_cells_z_;
    value_type                hx_;
    value_type                hy_;
    value_type                hz_;
    value_type                sigma_g_max_;
    std::vector<Box>          cells_;
    FlowProperties            flow_properties_;
    CellList                  cell_list_;
    Stats                     stats_;
    std::array<bool, 6>       is_diffuse_walls_;
    std::array<value_type, 6> diffuse_wall_temperature_;
};

}  // namespace mf

#endif  // MF_KINETICA_DOMAIN_H
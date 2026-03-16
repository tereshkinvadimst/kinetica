#ifndef MF_KINETICA_DOMAIN_H
#define MF_KINETICA_DOMAIN_H
#pragma once
#include <string>

#include "kinetica/Boundaries/boundaries.hh"
#include "kinetica/Box/box.hh"
#include "kinetica/CellList/cell_list.hh"
#include "kinetica/DSMC/wall.hh"
#include "kinetica/Particles/particles.hh"
#include "kinetica/Properties/flow_properties.hh"
#include "kinetica/Properties/profiles.hh"
#include "kinetica/Properties/stats.hh"
#include "kinetica/Random/random.hh"

namespace mf {

class Domain {
   public:
    using size_type  = std::size_t;
    using value_type = double;

    explicit Domain(Box        domain_box,
                    value_type m,
                    value_type W,
                    value_type molecule_size,
                    value_type scale_factor,
                    value_type time_scale_factor);

    void generateParticles(value_type n_density, value_type T, value_type x_min, value_type x_max);
    void generateMesh();
    void makeCellList();
    void updateCellList();
    void applyPeriodicBoundaries(bool px, bool py, bool pz);
    void moveParticles();
    void collideParticles();
    void saveXYZ(std::string file_name) const;
    void computeFlowProperties();
    void printStatsHeader();
    void printStats(value_type time);
    void addWall(std::shared_ptr<Wall> wall);
    void writeVTU(std::string file_name) const;
    auto getTimeStep() const noexcept -> value_type;
    void writeXProfile(std::string file_name);

   private:
    auto cellIndex(double x, double y, double z) -> size_type const;

    void computeMaxVelocity();

   private:
    random                             gen_;
    Box                                domain_box_;
    Particles                          particles_;
    size_type                          n_cells_x_;
    size_type                          n_cells_y_;
    size_type                          n_cells_z_;
    value_type                         sigma_g_max_;
    std::vector<Box>                   cells_;
    FlowProperties                     flow_properties_;
    CellList                           cell_list_;
    Stats                              stats_;
    value_type                         lambda_;
    value_type                         time_step_;
    value_type                         cell_size_;
    value_type                         scale_factor_;
    value_type                         max_velocity_;
    XProfiler                          xprofiler_;
    std::vector<std::shared_ptr<Wall>> walls_;
    value_type                         time_scale_factor_;
};

}  // namespace mf

#endif  // MF_KINETICA_DOMAIN_H
#ifndef MF_KINETICA_DOMAIN_H
#define MF_KINETICA_DOMAIN_H
#pragma once
#include <string>

#include "lab/Box/box.hh"
#include "lab/CellList/cell_list.hh"
#include "lab/Particles/particles.hh"
#include "lab/Random/xoshiro256.hh"

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

   private:
    auto cellIndex(double x, double y, double z) -> size_type const;

   private:
    xoshiro256          gen_;
    Box                 domain_box_;
    Particles           particles_;
    size_type           n_cells_x_;
    size_type           n_cells_y_;
    size_type           n_cells_z_;
    value_type          hx_;
    value_type          hy_;
    value_type          hz_;
    value_type          avg_sigma_g_max_;
    std::vector<Box>    cells_;
    std::vector<double> sigma_g_max_;
    CellList            cell_list_;
};

}  // namespace mf

#endif  // MF_KINETICA_DOMAIN_H
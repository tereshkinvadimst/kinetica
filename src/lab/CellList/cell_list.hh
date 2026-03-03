#ifndef MF_KINETICA_CELL_LIST_H
#define MF_KINETICA_CELL_LIST_H
#pragma once
#include <span>
#include <vector>

namespace mf {

class CellList final {
   public:
    using size_type = std::size_t;

    CellList()      = default;

    CellList(size_type n_cells, size_type n_particles);

    void build(const std::vector<size_type>& particle_cell_id);

    auto getParticlesInCell(size_type cell_id) const -> std::span<const size_type>;

   private:
    size_type              n_cells_;
    std::vector<size_type> cell_begin_;
    std::vector<size_type> particles_ids;
};

}  // namespace mf

#endif  // MF_KINETICA_CELL_LIST_H
#include "kinetica/CellList/cell_list.hh"

#include <algorithm>

mf::CellList::CellList(size_type n_cells, size_type n_particles)
    : n_cells_{n_cells}, cell_begin_(n_cells + 1), work_(n_cells + 1), particles_ids(n_particles) {}

void mf::CellList::build(const std::vector<size_type>& particle_cell_id) {
    cell_begin_.assign(cell_begin_.size(), size_type{});

    for (const auto cell_id : particle_cell_id) {
        ++cell_begin_[cell_id + 1];
    }

    for (size_type i = 1; i <= n_cells_; ++i) cell_begin_[i] += cell_begin_[i - 1];

    std::copy(cell_begin_.begin(), cell_begin_.begin() + n_cells_, work_.begin());

    for (size_type p{}; p < particle_cell_id.size(); ++p) {
        auto c                    = particle_cell_id[p];
        particles_ids[work_[c]++] = p;
    }
}

auto mf::CellList::getParticlesInCell(size_type cell_id) const -> std::span<const size_type> {
    return {particles_ids.data() + cell_begin_[cell_id], particles_ids.data() + cell_begin_[cell_id + 1]};
}
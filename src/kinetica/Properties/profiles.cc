#include "kinetica/Properties/profiles.hh"

#include <fstream>

#include "kinetica/Properties/flow_properties.hh"

mf::XProfiler::XProfiler(Box* domain_box, FlowProperties* flow_properties)
    : domain_box_{domain_box}, flow_properties_{flow_properties} {}

void mf::XProfiler::operator()(
    std::string file_name, size_type n_cells_x, size_type n_cells_y, size_type n_cells_z, value_type cell_size) {
    std::ofstream file(file_name);

    file << "X[m]\tNd[m^-3]\tux[m/s]\tuy[m/s]\tuz[m/s]\tu2x[(m/s)^2]\tu2y[(m/s)^2]\tu2z[(m/"
            "s)^2]\tTtrx[K]\tTtry[K]\tTtrz[K]\tTtr[K]\n";

    const size_type n_cells_yz = n_cells_y * n_cells_z;
    for (size_type i{}; i < n_cells_x; ++i) {
        auto n_density = value_type{};
        auto ux        = value_type{};
        auto uy        = value_type{};
        auto uz        = value_type{};
        auto u2x       = value_type{};
        auto u2y       = value_type{};
        auto u2z       = value_type{};
        auto Ttrx      = value_type{};
        auto Ttry      = value_type{};
        auto Ttrz      = value_type{};
        auto Ttr       = value_type{};

        for (size_type j{}; j < n_cells_y; ++j) {
            for (size_type k{}; k < n_cells_z; ++k) {
                const size_type c = i + n_cells_x * (j + n_cells_y * k);
                n_density += flow_properties_->n_density[c] / n_cells_yz;
                ux += flow_properties_->ux[c] / n_cells_yz;
                uy += flow_properties_->uy[c] / n_cells_yz;
                uz += flow_properties_->uz[c] / n_cells_yz;
                u2x += flow_properties_->u2x[c] / n_cells_yz;
                u2y += flow_properties_->u2y[c] / n_cells_yz;
                u2z += flow_properties_->u2z[c] / n_cells_yz;
                Ttrx += flow_properties_->Ttrx[c] / n_cells_yz;
                Ttry += flow_properties_->Ttry[c] / n_cells_yz;
                Ttrz += flow_properties_->Ttrz[c] / n_cells_yz;
                Ttr += flow_properties_->Ttr[c] / n_cells_yz;
            }
        }

        const value_type x = domain_box_->x0 + i * cell_size;

        file << x << '\t' << n_density << '\t' << ux << '\t' << uy << '\t' << uz << '\t' << u2x << '\t' << u2y << '\t' << u2z
             << '\t' << Ttrx << '\t' << Ttry << '\t' << Ttrz << '\t' << Ttr << '\n';
    }

    file.close();
}

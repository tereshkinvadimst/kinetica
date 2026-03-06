
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>

#include "kinetica/Box/box.hh"
#include "kinetica/Domain/domain.hh"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_kinetica, m) {
    // --- Box ---
    nb::class_<mf::Box>(m, "Box")
        .def(nb::init<double, double, double, double, double, double>(), "x0"_a, "y0"_a, "z0"_a, "Lx"_a, "Ly"_a, "Lz"_a)

        // координаты центра
        .def("cx", &mf::Box::cx)
        .def("cy", &mf::Box::cy)
        .def("cz", &mf::Box::cz)

        // проверка попадания точки
        .def("contains", &mf::Box::contains, "x"_a, "y"_a, "z"_a)

        // объём коробки
        .def("volume", &mf::Box::volume)

        // публичные поля как свойства
        .def_prop_rw(
            "x0", [](const mf::Box& b) { return b.x0; }, [](mf::Box& b, double v) { b.x0 = v; })
        .def_prop_rw(
            "y0", [](const mf::Box& b) { return b.y0; }, [](mf::Box& b, double v) { b.y0 = v; })
        .def_prop_rw(
            "z0", [](const mf::Box& b) { return b.z0; }, [](mf::Box& b, double v) { b.z0 = v; })
        .def_prop_rw(
            "Lx", [](const mf::Box& b) { return b.Lx; }, [](mf::Box& b, double v) { b.Lx = v; })
        .def_prop_rw(
            "Ly", [](const mf::Box& b) { return b.Ly; }, [](mf::Box& b, double v) { b.Ly = v; })
        .def_prop_rw("Lz", [](const mf::Box& b) { return b.Lz; }, [](mf::Box& b, double v) { b.Lz = v; });

    // --- Domain ---
    nb::class_<mf::Domain>(m, "Domain")
        .def(nb::init<mf::Box, double, double, double>(), "domain_box"_a, "mass"_a, "W"_a, "molecule_size"_a)

        // Методы для генерации частиц и сетки
        .def("generateParticles", &mf::Domain::generateParticles, "n_density"_a, "T"_a)
        .def("generateMesh", &mf::Domain::generateMesh, "hx"_a, "hy"_a, "hz"_a)
        .def("makeCellList", &mf::Domain::makeCellList)
        .def("updateCellList", &mf::Domain::updateCellList)
        .def("applyPeriodicBoundaries", &mf::Domain::applyPeriodicBoundaries, "px"_a, "py"_a, "pz"_a)
        .def("moveParticles", &mf::Domain::moveParticles, "dt"_a)
        .def("collideParticles", &mf::Domain::collideParticles, "dt"_a)
        .def("saveXYZ", &mf::Domain::saveXYZ, "file_name"_a)
        .def("computeFlowProperties", &mf::Domain::computeFlowProperties)
        .def("printStatsHeader", &mf::Domain::printStatsHeader)
        .def("printStats", &mf::Domain::printStats, "time"_a)
        .def("setDiffuseWall", &mf::Domain::setDiffuseWall, "side"_a, "Tw"_a)
        .def("writeVTU", &mf::Domain::writeVTU, "file_name"_a);
}

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>

#include "lab/Random/xoshiro256.hh"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(kinetics_lab, m) {
    nb::class_<mf::xoshiro256>(m, "Xoshiro256")

        .def(nb::init<mf::xoshiro256::result_type>(), "seed"_a)

        // вызов как функции: rng()
        .def("__call__", &mf::xoshiro256::operator())

        // аналоги стандартного интерфейса C++ RNG
        .def_static("min", &mf::xoshiro256::min)
        .def_static("max", &mf::xoshiro256::max);
}
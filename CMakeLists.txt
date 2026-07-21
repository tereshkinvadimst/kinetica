cmake_minimum_required(VERSION 3.30 FATAL_ERROR)

project(kinetica LANGUAGES C CXX Fortran)

add_executable(
  kinetica
  src/kinds.f90
  src/constants.f90
  src/utils.f90
  src/random.f90
  src/species.f90
  src/particles.f90
  src/domain.f90
  src/collider.f90
  src/boundary.f90
  src/piston.f90
  src/mover.f90
  src/global.f90
  src/init.f90
  src/xyz.f90
  src/kinetica.f90)

set_target_properties(
  kinetica
  PROPERTIES Fortran_STANDARD 2023
  Fortran_STANDARD_REQUIRED on
  Fortran_EXTENSIONS off)

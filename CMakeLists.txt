# Минимальная версия cmake проекта
cmake_minimum_required(VERSION 3.30 FATAL_ERROR)

# Настройки проекта
set(PYTHON_VERSION 3.13.1)

message(STATUS "CMAKE_CXX_COMPILER   = ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_TOOLCHAIN_FILE = ${CMAKE_TOOLCHAIN_FILE}")
message(STATUS "VCPKG_TARGET_TRIPLET = ${VCPKG_TARGET_TRIPLET}")
message(STATUS "VCPKG_HOST_TRIPLET   = ${VCPKG_HOST_TRIPLET}")

project(
    _kinetica
    VERSION 0.1.5
    LANGUAGES CXX
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/version.hh @ONLY
)

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    execute_process(
        COMMAND
            ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_BINARY_DIR}/compile_commands.json
            ${CMAKE_SOURCE_DIR}/compile_commands.json
    )
endif()

find_package(
    Python ${PYTHON_VERSION} REQUIRED
    COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule
)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_ROOT
)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING "Choose the type of build." FORCE
    )
    set_property(
        CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                        "RelWithDebInfo"
    )
endif()

find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
    # Name of the extension
    _kinetica
    # Target the stable ABI for Python 3.12+, which reduces the number of binary
    # wheels that must be built. This does nothing on older Python versions
    STABLE_ABI
    # Source code goes here
    src/kinetica/config.hh
    src/kinetica/Random/generators.hh
    src/kinetica/Random/xoshiro256.hh
    src/kinetica/Random/xoshiro256.cc
    src/kinetica/Random/random.hh
    src/kinetica/CellList/cell_list.hh
    src/kinetica/CellList/cell_list.cc
    src/kinetica/Particles/particles.hh
    src/kinetica/Particles/particles.cc
    src/kinetica/Box/box.hh
    src/kinetica/Boundaries/boundaries.hh
    src/kinetica/Boundaries/boundaries.cc
    src/kinetica/Properties/flow_properties.hh
    src/kinetica/Properties/flow_properties.cc
    src/kinetica/Properties/stats.hh
    src/kinetica/Properties/stats.cc
    src/kinetica/DSMC/mover.hh
    src/kinetica/DSMC/mover.cc
    src/kinetica/DSMC/collider.hh
    src/kinetica/DSMC/collider.cc
    src/kinetica/DSMC/wall.hh
    src/kinetica/DSMC/wall.cc
    src/kinetica/DSMC/utils.hh
    src/kinetica/Properties/profiles.hh
    src/kinetica/Properties/profiles.cc
    src/kinetica/Domain/domain.hh
    src/kinetica/Domain/domain.cc
    src/kinetica/core.cc
)

set_target_properties(
    _kinetica
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED on
               CXX_EXTENSIONS off
)

target_include_directories(
    _kinetica
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
           $<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
           ${CMAKE_CURRENT_BINARY_DIR}/generated
)

install(TARGETS _kinetica LIBRARY DESTINATION kinetica)

nanobind_add_stub(
    kinetica_stub
    MODULE
    _kinetica
    OUTPUT
    ${CMAKE_SOURCE_DIR}/src/kinetica/__init__.pyi
    PYTHON_PATH
    $<TARGET_FILE_DIR:_kinetica>
    DEPENDS
    _kinetica
)

find_package(Eigen3 CONFIG REQUIRED NO_MODULE)

target_link_libraries(_kinetica PUBLIC Eigen3::Eigen)

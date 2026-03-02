# Минимальная версия cmake проекта
cmake_minimum_required(VERSION 3.30 FATAL_ERROR)

# Настройки проекта
set(target_name kinetics-lab)
set(PYTHON_VERSION 3.13)

message(STATUS "CMAKE_CXX_COMPILER   = ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_TOOLCHAIN_FILE = ${CMAKE_TOOLCHAIN_FILE}")
message(STATUS "VCPKG_TARGET_TRIPLET = ${VCPKG_TARGET_TRIPLET}")
message(STATUS "VCPKG_HOST_TRIPLET   = ${VCPKG_HOST_TRIPLET}")

project(
    ${target_name}
    VERSION 0.0.1
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

set(lib_sources src/lab/core.cc)

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

find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
    # Name of the extension
    ${target_name}
    # Target the stable ABI for Python 3.12+, which reduces the number of binary
    # wheels that must be built. This does nothing on older Python versions
    STABLE_ABI
    # Source code goes here
    ${lib_sources}
)

install(TARGETS ${target_name} LIBRARY DESTINATION kinetica)

# nanobind_add_stub( mdca_stub MODULE ${target_name} OUTPUT
# ${CMAKE_SOURCE_DIR}/src/kinetica/__init__.pyi PYTHON_PATH
# $<TARGET_FILE_DIR:${target_name}> DEPENDS ${target_name} )

# Генерируем файл с макросами для экспорта функций
include(GenerateExportHeader)
generate_export_header(
    ${target_name} EXPORT_FILE_NAME
    ${CMAKE_CURRENT_BINARY_DIR}/generated/lib_export.hh
)

set_target_properties(
    ${target_name}
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED on
               CXX_EXTENSIONS off
)

target_include_directories(
    ${target_name}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
           $<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
           ${CMAKE_CURRENT_BINARY_DIR}/generated
)

find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(${target_name} PUBLIC Eigen3::Eigen)

# GridFire for CMake projects that fetch it (FetchContent) or add it as a subdirectory.
#
# GridFire is a Meson project; this file builds and installs it with Meson into this build tree, at configure time,
# and then defines the same imported targets that find_package(GridFire) does:
#
#   GridFire::gridfire          C++ library (C++23)
#   GridFire::gridfire_extern   C API: #include "gridfire/extern/gridfire_extern.h"
#   GridFire::gridfire_fortran  Fortran module `gridfire_mod`, compiled with your Fortran compiler (only if Fortran
#                               is enabled in your project before GridFire is made available)
#
# Example:
#
#   include(FetchContent)
#   FetchContent_Declare(GridFire
#       GIT_REPOSITORY https://github.com/4D-STAR/GridFire.git
#       GIT_TAG        v1.2.0
#       GIT_SHALLOW    TRUE
#       FIND_PACKAGE_ARGS 1.2)   # optional: use an installed GridFire (e.g. the pip wheel) when there is one
#   FetchContent_MakeAvailable(GridFire)
#   target_link_libraries(my_code PRIVATE GridFire::gridfire_extern)
#
# Needs Meson (>= 1.5) and Ninja on PATH (pip install meson ninja), C and C++23 compilers, and network access the
# first time (Meson downloads GridFire's dependencies). The first configure builds GridFire and takes a few minutes;
# later configures only rebuild what changed.
#
# Options (cache variables):
#   GRIDFIRE_MESON_OPTIONS          extra arguments for `meson setup`, e.g. "-Dopenmp_support=true"
#   GRIDFIRE_BUILDTYPE              Meson build type (default: release)
#   GRIDFIRE_USE_PROJECT_COMPILERS  build GridFire with this project's C and C++ compilers (default: ON)
#   GRIDFIRE_INSTALL                install GridFire (libraries, headers, ...) with this project (default: OFF)
cmake_minimum_required(VERSION 3.21)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/meson.build" _gf_project_line REGEX "^project\\('GridFire'")
string(REGEX MATCH "version: *'v?([0-9]+\\.[0-9]+\\.[0-9]+)" _gf_match "${_gf_project_line}")
project(GridFire VERSION ${CMAKE_MATCH_1} LANGUAGES NONE)

set(GRIDFIRE_MESON_OPTIONS "" CACHE STRING "Extra arguments for meson setup when building GridFire")
set(GRIDFIRE_BUILDTYPE "release" CACHE STRING "Meson build type for GridFire")
option(GRIDFIRE_USE_PROJECT_COMPILERS "Build GridFire with this project's C and C++ compilers" ON)
option(GRIDFIRE_INSTALL "Install GridFire together with this project" OFF)

find_program(GRIDFIRE_MESON_EXECUTABLE meson)
find_program(GRIDFIRE_NINJA_EXECUTABLE NAMES ninja ninja-build)
if(NOT GRIDFIRE_MESON_EXECUTABLE OR NOT GRIDFIRE_NINJA_EXECUTABLE)
    message(FATAL_ERROR
        "GridFire is built with Meson and Ninja, and they were not found (meson: ${GRIDFIRE_MESON_EXECUTABLE}, "
        "ninja: ${GRIDFIRE_NINJA_EXECUTABLE}). Install them (pip install meson ninja) or set "
        "GRIDFIRE_MESON_EXECUTABLE / GRIDFIRE_NINJA_EXECUTABLE. Alternatively install the gridfire Python package "
        "and use find_package(GridFire) with CMAKE_PREFIX_PATH set to the output of `gf-config --cmake-prefix`.")
endif()
execute_process(COMMAND "${GRIDFIRE_MESON_EXECUTABLE}" --version OUTPUT_VARIABLE _gf_meson_version
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if(_gf_meson_version VERSION_LESS 1.5)
    message(FATAL_ERROR "GridFire needs Meson 1.5 or newer; ${GRIDFIRE_MESON_EXECUTABLE} is ${_gf_meson_version}")
endif()

set(_gf_build_dir "${CMAKE_CURRENT_BINARY_DIR}/meson-build")
set(_gf_install_dir "${CMAKE_CURRENT_BINARY_DIR}/install")

# Only the libraries: no tests, tools, examples or Python bindings. b_ndebug=if-release because some Meson and
# Apple libc++ combinations reject the assertion macro Meson adds when b_ndebug=false.
set(_gf_setup_args
    --prefix=${_gf_install_dir} --libdir=lib --buildtype=${GRIDFIRE_BUILDTYPE} -Db_ndebug=if-release
    -Dbuild_tests=false -Dbuild_tools=false -Dbuild_examples=false -Dbuild_benchmarks=false -Dbuild_python=false
    -Dbuild_c_api=true -Dbuild_fortran=false -Dpkg_config=true)
separate_arguments(_gf_extra_args UNIX_COMMAND "${GRIDFIRE_MESON_OPTIONS}")
list(APPEND _gf_setup_args ${_gf_extra_args})

set(_gf_env "")
if(GRIDFIRE_USE_PROJECT_COMPILERS)
    if(CMAKE_C_COMPILER)
        list(APPEND _gf_env "CC=${CMAKE_C_COMPILER}")
    endif()
    if(CMAKE_CXX_COMPILER)
        list(APPEND _gf_env "CXX=${CMAKE_CXX_COMPILER}")
    endif()
endif()
list(APPEND _gf_env "NINJA=${GRIDFIRE_NINJA_EXECUTABLE}")

function(_gf_run what)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E env ${_gf_env} ${ARGN}
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" RESULT_VARIABLE _rc)
    if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "GridFire: ${what} failed (${_rc}). The Meson build directory is ${_gf_build_dir}.")
    endif()
endfunction()

# Set up again from scratch whenever the arguments or compilers change; otherwise reuse the build directory, so that
# reconfiguring the parent project only rebuilds what changed.
set(_gf_stamp "${_gf_build_dir}/cmake-setup-args.txt")
string(REPLACE ";" "\n" _gf_stamp_content "${_gf_env};${_gf_setup_args}")
set(_gf_previous "")
if(EXISTS "${_gf_stamp}")
    file(READ "${_gf_stamp}" _gf_previous)
endif()
if(NOT EXISTS "${_gf_build_dir}/build.ninja" OR NOT _gf_previous STREQUAL _gf_stamp_content)
    message(STATUS "GridFire: configuring with Meson in ${_gf_build_dir}")
    if(EXISTS "${_gf_build_dir}/build.ninja")
        set(_gf_wipe --wipe)
    else()
        set(_gf_wipe "")
    endif()
    _gf_run("meson setup" "${GRIDFIRE_MESON_EXECUTABLE}" setup ${_gf_wipe} ${_gf_setup_args}
            "${_gf_build_dir}" "${CMAKE_CURRENT_SOURCE_DIR}")
    file(WRITE "${_gf_stamp}" "${_gf_stamp_content}")
endif()
message(STATUS "GridFire: building and installing into ${_gf_install_dir} (the first build takes a few minutes)")
_gf_run("meson compile" "${GRIDFIRE_MESON_EXECUTABLE}" compile -C "${_gf_build_dir}")
_gf_run("meson install" "${GRIDFIRE_MESON_EXECUTABLE}" install -C "${_gf_build_dir}" --only-changed --quiet)

# The installed GridFireConfig.cmake defines the targets; global, so the project that fetched GridFire sees them
set(_gf_config "${_gf_install_dir}/lib/cmake/GridFire/GridFireConfig.cmake")
set(GridFire_IMPORT_GLOBAL TRUE)
set(GridFire_FIND_COMPONENTS "")
include("${_gf_config}")
if(NOT GridFire_FOUND)
    message(FATAL_ERROR "GridFire was built but its targets could not be defined: ${GridFire_NOT_FOUND_MESSAGE}")
endif()
set(GridFire_DIR "${_gf_install_dir}/lib/cmake/GridFire" CACHE PATH "GridFire package configuration" FORCE)
# Variables set here stay in this directory's scope; cache them so the project that fetched GridFire sees them too
set(GridFire_VERSION "${GridFire_VERSION}" CACHE INTERNAL "GridFire version")
set(GridFire_FORTRAN_SOURCE "${GridFire_FORTRAN_SOURCE}" CACHE INTERNAL "GridFire Fortran module source")

# Later find_package(GridFire) calls in the parent project resolve to this build (CMake >= 3.24)
if(CMAKE_FIND_PACKAGE_REDIRECTS_DIR)
    file(WRITE "${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/gridfire-extra.cmake" "include(\"${_gf_config}\")\n")
    file(WRITE "${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/gridfire-config-version.cmake"
         "include(\"${_gf_install_dir}/lib/cmake/GridFire/GridFireConfigVersion.cmake\")\n")
endif()

if(GRIDFIRE_INSTALL)
    install(DIRECTORY "${_gf_install_dir}/" DESTINATION "." USE_SOURCE_PERMISSIONS)
endif()

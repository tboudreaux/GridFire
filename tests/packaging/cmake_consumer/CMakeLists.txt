# Builds the packaging consumers (../consumer.{cpp,c,f90}) against GridFire from CMake, in one of two ways:
#
#   find_package (default): an installed GridFire, e.g. the pip wheel:
#       cmake -S . -B build -DCMAKE_PREFIX_PATH="$(gf-config --cmake-prefix)"
#   FetchContent (-DGRIDFIRE_CONSUMER_FETCH=ON): GridFire fetched from git and built with Meson. To build a local
#   checkout instead of the tagged release, add -DFETCHCONTENT_SOURCE_DIR_GRIDFIRE=/path/to/GridFire.
#
# Then: cmake --build build && ctest --test-dir build
cmake_minimum_required(VERSION 3.24)
project(gridfire_cmake_consumer LANGUAGES C CXX Fortran)

option(GRIDFIRE_CONSUMER_FETCH "Fetch GridFire with FetchContent instead of find_package" OFF)
set(GRIDFIRE_CONSUMER_GIT_TAG "v1.2.0" CACHE STRING "GridFire release to fetch")

if(GRIDFIRE_CONSUMER_FETCH)
    include(FetchContent)
    FetchContent_Declare(GridFire
        GIT_REPOSITORY https://github.com/4D-STAR/GridFire.git
        GIT_TAG        ${GRIDFIRE_CONSUMER_GIT_TAG}
        GIT_SHALLOW    TRUE
        # Use an installed GridFire when one is found; fetch and build it otherwise
        FIND_PACKAGE_ARGS 1.2 COMPONENTS CXX C Fortran)
    FetchContent_MakeAvailable(GridFire)
    # Later find_package calls (e.g. in subdirectories) resolve to what was just fetched
    find_package(GridFire 1.2 REQUIRED COMPONENTS CXX C Fortran)
else()
    find_package(GridFire 1.2 REQUIRED COMPONENTS CXX C Fortran)
endif()
message(STATUS "GridFire ${GridFire_VERSION}; Fortran module source: ${GridFire_FORTRAN_SOURCE}")

set(src "${CMAKE_CURRENT_SOURCE_DIR}/..")

add_executable(consumer_cpp "${src}/consumer.cpp")
target_link_libraries(consumer_cpp PRIVATE GridFire::gridfire)

add_executable(consumer_c "${src}/consumer.c")
target_link_libraries(consumer_c PRIVATE GridFire::gridfire_extern)

add_executable(consumer_f "${src}/consumer.f90")
target_link_libraries(consumer_f PRIVATE GridFire::gridfire_fortran)

enable_testing()
foreach(t IN ITEMS cpp c f)
    add_test(NAME consumer_${t} COMMAND consumer_${t})
endforeach()

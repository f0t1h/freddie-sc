cmake_minimum_required(VERSION 3.20)

# scFreddie superbuild: build the three component tools, each via its own CMake,
# into its own build/ dir (the binary paths the Snakemake workflows expect):
#     freddie-core/build/freddie-core
#     sclr-putagene/build/putagene
#     extern/tksm2/build/tksm2
#
#   cmake -S . -B build [-DFREDDIE_ENABLE_GUROBI=ON]
#   cmake --build build -j
#
# htslib/zlib are discovered from the active conda env; activate it first
# (micromamba activate scfreddie) or pass -DCMAKE_PREFIX_PATH=/path/to/env.
project(scfreddie NONE)

include(ExternalProject)

option(FREDDIE_ENABLE_GUROBI "Build freddie-core's Gurobi backend" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Prefer the active conda env for htslib/zlib discovery.
set(_prefix_path ${CMAKE_PREFIX_PATH})
if(DEFINED ENV{CONDA_PREFIX})
  list(PREPEND _prefix_path "$ENV{CONDA_PREFIX}")
endif()

set(_common_cmake_args
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  "-DCMAKE_PREFIX_PATH=${_prefix_path}")

# Each component is an independent ExternalProject built in-place under its own
# build/ dir. BUILD_ALWAYS re-runs the component build so source edits are picked
# up; no install step (the Snakefiles run the binaries straight from build/).
#
# No custom BUILD_COMMAND: the default ($(MAKE)) shares the top-level make
# jobserver, so `cmake --build build -jN` caps total compile jobs at N across all
# three components. A custom `cmake --build --parallel` would instead give each
# component all cores at once (3x oversubscription -> OOM in tight containers).
function(scfreddie_component name src)
  ExternalProject_Add(${name}
    SOURCE_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/${src}
    BINARY_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/${src}/build
    CMAKE_ARGS   ${_common_cmake_args} ${ARGN}
    BUILD_ALWAYS TRUE
    INSTALL_COMMAND ""
    USES_TERMINAL_CONFIGURE TRUE
    USES_TERMINAL_BUILD TRUE)
endfunction()

scfreddie_component(freddie-core freddie-core
  -DFREDDIE_ENABLE_GUROBI=${FREDDIE_ENABLE_GUROBI})
scfreddie_component(putagene sclr-putagene)
scfreddie_component(tksm2 extern/tksm2)

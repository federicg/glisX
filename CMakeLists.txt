# ==============================================================================
#  Project Setup
# ==============================================================================
cmake_minimum_required(VERSION 3.17)
project(glisX LANGUAGES C CXX Fortran)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# ==============================================================================
# Stage 1: Dependencies
# ==============================================================================
include(ExternalProject)

set(DEPS_PREFIX ${CMAKE_BINARY_DIR}/_deps)

# ----------------------------
# Build LIS
# ----------------------------
ExternalProject_Add(
  lis_external
  URL https://www.ssisc.org/lis/dl/lis-2.1.10.zip
  PREFIX ${DEPS_PREFIX}/lis
  BUILD_IN_SOURCE TRUE

  CONFIGURE_COMMAND ./configure --enable-mpi --enable-omp --prefix=${DEPS_PREFIX}/lis/install
  BUILD_COMMAND $(MAKE)
  INSTALL_COMMAND $(MAKE) install

  UPDATE_COMMAND ""
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# ----------------------------
# Build OCTAVE
# ----------------------------
ExternalProject_Add(
  octave_external
  URL https://ftp.gnu.org/gnu/octave/octave-6.2.0.tar.gz
  PREFIX ${DEPS_PREFIX}/octave
  BUILD_IN_SOURCE TRUE

  CONFIGURE_COMMAND ./configure --enable-mpi --enable-omp --prefix=${DEPS_PREFIX}/octave/install
  BUILD_COMMAND $(MAKE)
  INSTALL_COMMAND $(MAKE) install

  UPDATE_COMMAND ""
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# ----------------------------
# Build OCTAVE_FILE_IO
# ----------------------------
ExternalProject_Add(
  octave_file_io_external
  GIT_REPOSITORY https://github.com/carlodefalco/octave_file_io.git
  PREFIX ${DEPS_PREFIX}/octave_file_io
  BUILD_IN_SOURCE TRUE

  CONFIGURE_COMMAND ./autogen.sh && ./configure --with-octave-home=${DEPS_PREFIX}/octave/install/bin CPPFLAGS=-I${DEPS_PREFIX}/octave/install/include/octave-6.2.0
  --prefix=${DEPS_PREFIX}/octave_file_io/install
  BUILD_COMMAND $(MAKE)
  INSTALL_COMMAND $(MAKE) install

  UPDATE_COMMAND ""
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
add_dependencies(octave_file_io_external octave_external)

# ----------------------------
# Build BIMPP
# ----------------------------
ExternalProject_Add(
  bimpp_external
  GIT_REPOSITORY https://github.com/carlodefalco/bimpp.git
  GIT_TAG landslides
  PREFIX ${DEPS_PREFIX}/bimpp
  BUILD_IN_SOURCE TRUE

  CONFIGURE_COMMAND
    ./autogen.sh &&
    ./configure
    --prefix=${DEPS_PREFIX}/bimpp/install
    --disable-shared
    --enable-static
    FC=mpif90
    CC=mpicc
    CXX=mpicxx
    CFLAGS=-fopenmp
    CXXFLAGS=-fopenmp
    FFLAGS=-fopenmp
    LDFLAGS=-L${DEPS_PREFIX}/octave/install/lib\ -L${DEPS_PREFIX}/octave_file_io/install/lib\ -L${DEPS_PREFIX}/lis/install/lib
    CPPFLAGS=-I${DEPS_PREFIX}/octave/install/include\ -I${DEPS_PREFIX}/octave_file_io/install/include\ -I${DEPS_PREFIX}/lis/install/include\ -DOMPI_SKIP_MPICXX\ -DBIM_TIMING
    --with-lis-home=${DEPS_PREFIX}/lis/install
    --with-mumps-extra-libs=-ldmumps\ -lmumps_common\ -lmetis\ -lpord\ -lesmumps\ -lscotch\ -lscotcherr\ -lscalapack\ -lopenblas\ -lgfortran\ -lpthread\ -lgomp
    --with-octave-home=${DEPS_PREFIX}/octave/install/bin
    --with-octave_file_io-home=${DEPS_PREFIX}/octave_file_io/install

  BUILD_COMMAND $(MAKE) -C src
  INSTALL_COMMAND $(MAKE) -C src install
  UPDATE_COMMAND ""

  DEPENDS lis_external octave_external octave_file_io_external
)

# ----------------------------
# Stage 1 target: build all deps
# ----------------------------
add_custom_target(deps
  DEPENDS lis_external octave_external octave_file_io_external bimpp_external
)

# ==============================================================================
# Import BIM++ as a CMake library target (available after 'make deps')
# ==============================================================================
add_library(bimpp INTERFACE IMPORTED GLOBAL)
add_dependencies(bimpp bimpp_external)

# Tell CMake these dirs will exist at build time
file(MAKE_DIRECTORY ${DEPS_PREFIX}/bimpp/install/include)
file(MAKE_DIRECTORY ${DEPS_PREFIX}/octave_file_io/install/include)
file(MAKE_DIRECTORY ${DEPS_PREFIX}/octave/install/include/octave-6.2.0)
file(MAKE_DIRECTORY ${DEPS_PREFIX}/lis/install/include)

set_target_properties(bimpp PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES
    "${DEPS_PREFIX}/bimpp/install/include;\
${DEPS_PREFIX}/octave_file_io/install/include;\
${DEPS_PREFIX}/octave/install/include/octave-6.2.0;\
${DEPS_PREFIX}/lis/install/include"
  INTERFACE_LINK_LIBRARIES
    "-Wl,--start-group;\
${DEPS_PREFIX}/bimpp/install/lib/libbimsystemsolver.a;\
${DEPS_PREFIX}/bimpp/install/lib/libbimmumps.a;\
${DEPS_PREFIX}/bimpp/install/lib/libbimlis.a;\
${DEPS_PREFIX}/bimpp/install/lib/libbimlinalg.a;\
${DEPS_PREFIX}/bimpp/install/lib/libbim.a;\
${DEPS_PREFIX}/bimpp/install/lib/libbimp4est.a;\
-Wl,--end-group;\
${DEPS_PREFIX}/octave_file_io/install/lib/liboctave_file_io.so;\
${DEPS_PREFIX}/octave/install/lib/octave/6.2.0/liboctinterp.so;\
${DEPS_PREFIX}/octave/install/lib/octave/6.2.0/liboctave.so;\
-L${DEPS_PREFIX}/lis/install/lib;\
-ldmumps;-lmumps_common;-lmetis;-lpord;\
-lesmumps;-lscotch;-lscotcherr;\
-lscalapack;-lopenblas;-lgfortran;-lpthread;-lgomp;\
-lp4est;-lsc"
)

# ==============================================================================
# Fetch nlohmann/json (single header)
# ==============================================================================
include(FetchContent)

FetchContent_Declare(
  nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp
  DOWNLOAD_NO_EXTRACT TRUE
  DOWNLOAD_DIR ${DEPS_PREFIX}/nlohmann
)
FetchContent_MakeAvailable(nlohmann_json)

# ==============================================================================
# Stage 2: glisX executables
# ==============================================================================

# Compile-time test-case selectors (see source/*_macro_input.h). Leave unset to
# use the default baked in the header, or pass e.g. -DTWO_PHASE_TEST=3.
set(SINGLE_PHASE_TEST "" CACHE STRING "Single-phase test case (see single_phase_macro_input.h)")
set(TWO_PHASE_TEST    "" CACHE STRING "Two-phase test case (see two_phase_macro_input.h)")

# Single-phase depth-averaged granular flow solver.
add_executable(glisX-single-phase
  source/main_TG2IMEX-RKC_Strang_balanced.cpp
  source/Taylor_Galerkin_IMEX-RKC_Strang_balanced.cpp
)

# Two-phase (solid + fluid) depth-averaged debris flow solver.
add_executable(glisX-two-phase
  source/main_TG2IMEX-RKC_Strang_balanced_two_phase.cpp
  source/Taylor_Galerkin_IMEX-RKC_Strang_balanced_two_phase.cpp
)

foreach(tgt glisX-single-phase glisX-two-phase)
  add_dependencies(${tgt} deps)
  target_include_directories(${tgt} PRIVATE ${DEPS_PREFIX}/nlohmann)
  target_link_libraries(${tgt} PRIVATE bimpp MPI::MPI_CXX MPI::MPI_C OpenMP::OpenMP_CXX)
endforeach()

if(NOT SINGLE_PHASE_TEST STREQUAL "")
  target_compile_definitions(glisX-single-phase PRIVATE SINGLE_PHASE_TEST=${SINGLE_PHASE_TEST})
endif()
if(NOT TWO_PHASE_TEST STREQUAL "")
  target_compile_definitions(glisX-two-phase PRIVATE TWO_PHASE_TEST=${TWO_PHASE_TEST})
endif()

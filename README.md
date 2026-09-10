# glisX

A finite element simulation code on adaptive quadtrees for depth-averaged
**landslide / debris flow** dynamics, built on top of
[BIM++](https://github.com/carlodefalco/bimpp) (`landslides` branch).

Two solvers are provided, both using a second-order optimally-stable IMEX
(pseudo-) staggered Galerkin discretization with Strang splitting and
well-balancing on adaptive quadtree meshes:

| Executable | Source | Model |
|---|---|---|
| `glisX-single-phase` | `main_TG2IMEX-RKC_Strang_balanced.cpp` | Single-phase depth-averaged granular / mud flow (Herschel–Bulkley type rheology, bed friction, turbulence). |
| `glisX-two-phase`    | `main_TG2IMEX-RKC_Strang_balanced_two_phase.cpp` | Two-phase (solid + interstitial fluid) depth-averaged debris flow with excess pore-water pressure, erosion and consolidation. |

Both solvers share the same code lineage as the
[`lava-flow`](https://github.com/federicg/lava-flow) code; `glisX` is the
landslide branch of that work, developed during the author's PhD.

## Dependencies

All dependencies are downloaded and built automatically by CMake:

| Library | Version / Source | Purpose |
|---|---|---|
| [LIS](https://www.ssisc.org/lis/) | 2.1.10 | Iterative linear solvers |
| [GNU Octave](https://www.gnu.org/software/octave/) | 6.2.0 | Scripting and post-processing |
| [octave_file_io](https://github.com/carlodefalco/octave_file_io) | git HEAD | Octave-based file I/O bridge |
| [BIM++](https://github.com/carlodefalco/bimpp) | `landslides` branch | Core finite-element/volume library |
| [nlohmann/json](https://github.com/nlohmann/json) | 3.11.3 | JSON input parsing (single header) |

## System Requirements

The following must be available on your system before building:

- CMake ≥ 3.17 (loaded: cmake/3.30.5)
- MPI compilers: `mpicc`, `mpicxx`, `mpif90` (loaded: openmpi/4.1.7)
- OpenMP-capable compilers — GCC recommended (loaded: gcc/12.2.0)
- Fortran compiler: `gfortran` (included with GCC)
- Standard build tools: `make`, `autoconf`, `automake`, `libtool`
- MUMPS and its dependencies (loaded: mumps):
  - `dmumps`, `mumps_common`, `pord`
  - `metis/5.1.0`, `parmetis/4.0.3`
  - `scotch/7.0.4`, `scotcherr`, `esmumps`
  - `scalapack` (loaded: netlib-scalapack/2.2.0)
  - `openblas/0.3.28`
- p4est/2.8 (adaptive mesh refinement)
- Lua/5.4.6 (if used for configuration/scripting)
- `ncurses/6.5`, `readline/8.2`, `pcre/8.45` (terminal/scripting dependencies)

## Building

```bash
# 1. Clone the repository
git clone <repo-url>
cd glisX

# 2. Configure
cmake -S . -B build

# 3. Build dependencies (first time only, this will take a while)
cd build
make deps -j$(nproc)

# 4. Build the glisX executables
make -j$(nproc)
```

After the first `make deps`, subsequent builds only require `make -j$(nproc)`.

## Test case data

The example input files expect two rasters in `inputs/`:

- `inputs/dem.octbin.gz` — bedrock elevation, Octave variable `dem`
- `inputs/mask_in.octbin.gz` — initial release indicator, Octave variable `mask_in`

An idealised inclined-plane test case (101 × 101 raster, 20° slope, circular
release) is shipped in `inputs/`. To regenerate it or produce your own:

```bash
cd scripts/octave
octave --no-gui make_testcase.m
```

For a real scenario, replace these files with your own DEM and release mask
(same Octave variable names, matching the `number raster columns/rows` and
`raster resolution` fields in the JSON input).

## Running

The simulations are configured via JSON input files at the repository root.
Run from the `build` directory:

```bash
mkdir -p results

# Single-phase granular / mud flow
mpirun -np NPROCS glisX-single-phase ../glisX_input-single-phase.json

# Two-phase debris flow
mpirun -np NPROCS glisX-two-phase ../glisX_input-two-phase.json
```

Results are written as `results/swe_*_<step>_<rank>.octbin.gz`.

## Post-processing

From the `build` directory, convert the octbin output to ParaView `.vtu`:

```bash
# single-phase run
../convertTovtk.sh NPROCS . single

# two-phase run
../convertTovtk.sh NPROCS . two
```

This produces `glisX_<step>.vtu` files.

## Citation

If you use this code in your research, please cite the associated PhD thesis /
paper:

```bibtex
% TODO: replace with the final reference
@phdthesis{gatti_phd_landslide,
  title  = {Numerical modelling of landslide and debris flow dynamics with
            adaptive staggered Galerkin schemes},
  author = {Gatti, Federico},
  school = {Politecnico di Milano},
  year   = {2024}
}
```

See also the companion lava-flow paper:

```bibtex
@article{gatti2025second,
  title={Second-order Optimally Stable IMEX (pseudo-) staggered Galerkin discretization: application to lava flow modeling},
  author={Gatti, Federico and Orlando, Giuseppe},
  journal={arXiv preprint arXiv:2509.09460},
  year={2025}
}
```

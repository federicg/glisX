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

## Reproducing the published benchmarks

The initial condition and geometry of **both** solvers are selected at
**compile time** through a macro (same idea as `lava_macro_input.h` in
lava-flow); every physical parameter stays in the JSON. Set the macro — by
editing the header, or with `cmake -S . -B build -D<MACRO>=<n>` — rebuild the
target, and run it with the matching JSON.

### Single-phase — `SINGLE_PHASE_TEST` (`source/single_phase_macro_input.h`)

> F. Gatti, C. de Falco, S. Perotto, L. Formaggia,
> *A scalable well-balanced numerical scheme for the simulation of fast
> landslides with efficient time stepping*, Appl. Math. Comput. 468 (2024) 128525.

| `SINGLE_PHASE_TEST` | Paper | Description | JSON |
|:--:|:--|:--|:--|
| `0` | Sec. 4.3 | Generic run — `Z`, release read from rasters (Bindo-Cortenova real case, `H0 = 38 m` where mask == 1) | `glisX_input-single-phase.json`, `glisX_input-bindo-cortenova.json` |
| `1` | Sec. 4.1.1 | Viscous dam break, frictionless (Eq. 30) | `glisX_input-single-phase-dambreak.json` |
| `2` | Sec. 4.1.1 | Smooth solution / convergence, `H(x,0) = Z(x)` (Eq. 31–32) | `glisX_input-single-phase-smooth.json` |
| `3` | Sec. 4.1.2 | Well-balancing, smooth topography (Eq. 33) | `glisX_input-single-phase-wb-smooth.json` |
| `4` | Sec. 4.1.2 | Well-balancing, discontinuous topography (Eq. 34) | `glisX_input-single-phase-wb-discontinuous.json` |
| `5` | Sec. 4.2.1 | Radial dam break, flat frictionless bed, Newtonian (Eq. 36) | `glisX_input-single-phase-radial-newtonian.json` |
| `6` | Sec. 4.2.2 | Radial dam break, Bingham rheology (Eq. 36) | `glisX_input-single-phase-radial-bingham.json` |
| `7` | Sec. 4.2.3 | Granular slide on a ~22° inclined plane (Eq. 37) | `glisX_input-single-phase-inclined-plane.json` |

```bash
cmake -S . -B build -DSINGLE_PHASE_TEST=3
cmake --build build --target glisX-single-phase -j$(nproc)
mpirun -np 4 build/glisX-single-phase glisX_input-single-phase-wb-smooth.json
```

### Two-phase — `TWO_PHASE_TEST` (`source/two_phase_macro_input.h`)

All physical parameters stay in the JSON file, so reproducing a
benchmark from

> F. Gatti, C. de Falco, S. Perotto, L. Formaggia, M. Pastor,
> *A scalable well-balanced numerical scheme for the modeling of two-phase
> shallow granular landslide consolidation*, J. Comput. Phys. 501 (2024) 112798.

means: set the macro, rebuild `glisX-two-phase`, and run it with the matching
JSON.

| `TWO_PHASE_TEST` | Paper | Description | JSON |
|:--:|:--|:--|:--|
| `0` | Sec. 4.2 | Generic run — `Z` and initial `h` read from the octbin rasters in the JSON | `glisX_input-two-phase.json`, `glisX_input-sham-tseng.json` |
| `1` | Sec. 4.1.1 | Well-balancing, smooth topography `Z = 5 e^{-2/5 (x-5)^2}` | `glisX_input-two-phase-wb-smooth.json` |
| `2` | Sec. 4.1.1 | Well-balancing, discontinuous topography (`Z = 4` on `4 ≤ x ≤ 8`) | `glisX_input-two-phase-wb-discontinuous.json` |
| `3` | Sec. 4.1.2 | Loss of hyperbolicity, null excess pwp (Eq. 40) | `glisX_input-two-phase-hyperbolicity-null.json` |
| `4` | Sec. 4.1.2 | Loss of hyperbolicity, non-null excess pwp (Eq. 40 + linear pwp) | `glisX_input-two-phase-hyperbolicity-pwp.json` |
| `5` | Sec. 4.1.3 | Dam break, Riemann problem, infinite permeability (Eq. 42) | `glisX_input-two-phase-dambreak-riemann.json` |
| `6` | Sec. 4.1.3 | Wet-dry dam break, 2D square release | `glisX_input-two-phase-dambreak-wetdry-2d.json` |
| `7` | Sec. 4.1.3 | Wet-dry dam break, 1D release | `glisX_input-two-phase-dambreak-wetdry-1d.json` |
| `8` | Sec. 4.1.4 | Efficiency test, discontinuous topography (Eq. 43) | `glisX_input-two-phase-efficiency.json` |

```bash
# example: well-balancing test with the discontinuous bed
cmake -S . -B build -DTWO_PHASE_TEST=2
cmake --build build --target glisX-two-phase -j$(nproc)
mpirun -np 4 build/glisX-two-phase glisX_input-two-phase-wb-discontinuous.json
```

Notes (both solvers):
- The analytic tests **ignore** the `dem file` / `mask file` entries (they may
  be omitted from those JSONs). Only `*_TEST 0` reads rasters.
- The domain is the box `[0, L] × [0, H]` with `L = res·(Nx−1)`,
  `H = res·(Ny−1)` taken from the JSON raster fields; adjust those to match the
  paper domain.
- The JSON files shipped here reproduce each paper's **geometry and initial
  data**; the remaining physical coefficients (friction, viscosity, `C_v`,
  `E_m`, `V_T`, thresholds, densities) are set to the paper values as far as
  they are stated — double-check them against the article before a production
  run.

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

If you use this code in your research, please cite the relevant paper(s).

**Single-phase solver** (`glisX-single-phase`):

```bibtex
@article{gatti2024singlephase,
  title   = {A scalable well-balanced numerical scheme for the simulation of
             fast landslides with efficient time stepping},
  author  = {Gatti, Federico and de Falco, Carlo and Perotto, Simona and
             Formaggia, Luca},
  journal = {Applied Mathematics and Computation},
  volume  = {468},
  pages   = {128525},
  year    = {2024},
  doi     = {10.1016/j.amc.2023.128525}
}
```

**Two-phase solver** (`glisX-two-phase`):

```bibtex
@article{gatti2024twophase,
  title   = {A scalable well-balanced numerical scheme for the modeling of
             two-phase shallow granular landslide consolidation},
  author  = {Gatti, Federico and de Falco, Carlo and Perotto, Simona and
             Formaggia, Luca and Pastor, Manuel},
  journal = {Journal of Computational Physics},
  volume  = {501},
  pages   = {112798},
  year    = {2024},
  doi     = {10.1016/j.jcp.2024.112798}
}
```

Companion lava-flow paper (shared scheme lineage):

```bibtex
@article{gatti2025lava,
  title   = {Second-order Optimally Stable IMEX (pseudo-) staggered Galerkin
             discretization: application to lava flow modeling},
  author  = {Gatti, Federico and Orlando, Giuseppe},
  journal = {arXiv preprint arXiv:2509.09460},
  year    = {2025}
}
```

#ifndef SINGLE_PHASE_MACRO_INPUT_H
#define SINGLE_PHASE_MACRO_INPUT_H

// ============================================================================
//  Compile-time selection of the single-phase test case.
//
//  As for the two-phase solver, all runtime parameters (density, viscosity,
//  yield stress, friction, turbulence, thresholds, final time, mesh sizes,
//  ...) stay in the JSON input file.  This macro only selects the geometry
//  and the initial condition (topography Z, material height H, mass fluxes
//  U_x, U_y), so that the benchmarks of
//
//    F. Gatti, C. de Falco, S. Perotto, L. Formaggia,
//    "A scalable well-balanced numerical scheme for the simulation of fast
//     landslides with efficient time stepping",
//    Appl. Math. Comput. 468 (2024) 128525.  doi:10.1016/j.amc.2023.128525
//
//  can be reproduced by only changing SINGLE_PHASE_TEST and rebuilding.
//
//  SINGLE_PHASE_TEST:
//    0  Generic run: topography and the initial release are read from the
//       "dem file"/"mask file" octbin rasters named in the JSON.  Used for the
//       Bindo-Cortenova real case study (Sec. 4.3): H0 = 38 m where mask == 1.
//       -> glisX_input-single-phase.json  /  glisX_input-bindo-cortenova.json
//
//    1  Sec. 4.1.1  Viscous dam break, frictionless (Eq. 30)
//                    H = 1,  U_x = 1 for x <= L/2 else 0.5.
//                    -> glisX_input-single-phase-dambreak.json
//    2  Sec. 4.1.1  Smooth solution / convergence (Eq. 31-32)
//                    Z(x) = 1 + 1/10 exp(-50/L^2 (x-L/2)^2),  H(x,0) = Z(x).
//                    -> glisX_input-single-phase-smooth.json
//    3  Sec. 4.1.2  Well-balancing, smooth topography (Eq. 33)
//                    Z(x) = 5 exp(-2/5 (x-5)^2),  H = 10 - Z,  at rest.
//                    -> glisX_input-single-phase-wb-smooth.json
//    4  Sec. 4.1.2  Well-balancing, discontinuous topography (Eq. 34)
//                    Z(x) = 4 for 4 <= x <= 8 else 0,  H = 10 - Z,  at rest.
//                    -> glisX_input-single-phase-wb-discontinuous.json
//    5  Sec. 4.2.1  Example 1: radial dam break, flat frictionless bed,
//                    Newtonian rheology (Eq. 36)   H = 2 for r <= 1/2 else 1.
//                    -> glisX_input-single-phase-radial-newtonian.json
//    6  Sec. 4.2.2  Example 2: same release as 5, Bingham rheology.
//                    -> glisX_input-single-phase-radial-bingham.json
//    7  Sec. 4.2.3  Example 3: granular slide on a ~22 deg inclined plane,
//                    initial release of Eq. (37).
//                    -> glisX_input-single-phase-inclined-plane.json
//
//  Notes:
//   * Tests 1-7 are analytic and ignore the "dem file"/"mask file" entries.
//     Only SINGLE_PHASE_TEST 0 reads rasters.
//   * The domain is [0, L] x [0, H] with L = res*(Nx-1), H = res*(Ny-1) from
//     the JSON raster fields.
//   * All these tests start from rest (U_x = U_y = 0) except test 1.
//   * Test 7 assumes the inclined plane descends along +x,
//     Z(x) = (L - x) tan(22 deg); flip the sign if your set-up is mirrored.
// ============================================================================

// Edit the default here, or override from CMake with -DSINGLE_PHASE_TEST=<n>
// (cmake -S . -B build -DSINGLE_PHASE_TEST=3).
#ifndef SINGLE_PHASE_TEST
#define SINGLE_PHASE_TEST 0
#endif

#endif // SINGLE_PHASE_MACRO_INPUT_H

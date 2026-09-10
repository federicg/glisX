#ifndef TWO_PHASE_MACRO_INPUT_H
#define TWO_PHASE_MACRO_INPUT_H

// ============================================================================
//  Compile-time selection of the two-phase test case.
//
//  All the runtime parameters (densities, friction, thresholds, C_v, E_m,
//  final time, mesh sizes, ...) stay in the JSON input file.  This macro only
//  selects the *geometry and the initial condition* (topography Z, porosity
//  n_w, total material height h, and the phase mass fluxes U_w, U_s), so that
//  the benchmarks of
//
//    F. Gatti, C. de Falco, S. Perotto, L. Formaggia, M. Pastor,
//    "A scalable well-balanced numerical scheme for the modeling of two-phase
//     shallow granular landslide consolidation",
//    J. Comput. Phys. 501 (2024) 112798.  doi:10.1016/j.jcp.2024.112798
//
//  can be reproduced by only changing TWO_PHASE_TEST and rebuilding.
//
//  TWO_PHASE_TEST:
//    0  Generic run: topography and initial material height are read from the
//       "dem file"/"mask file" octbin rasters named in the JSON.  Use this for
//       real scenarios (e.g. the Sham Tseng San Tsuen debris flow, Sec. 4.2)
//       via  glisX_input-two-phase.json  /  glisX_input-sham-tseng.json.
//
//    1  Sec. 4.1.1  Well-balancing, smooth topography
//                    Z(x) = 5 exp(-2/5 (x-5)^2),  h = 10 - Z,  at rest.
//                    -> glisX_input-two-phase-wb-smooth.json
//    2  Sec. 4.1.1  Well-balancing, discontinuous topography
//                    Z(x) = 4 for 4 <= x <= 8 else 0,  h = 10 - Z,  at rest.
//                    -> glisX_input-two-phase-wb-discontinuous.json
//    3  Sec. 4.1.2  Loss of hyperbolicity, null excess pwp (Eq. 40)
//                    -> glisX_input-two-phase-hyperbolicity-null.json
//    4  Sec. 4.1.2  Loss of hyperbolicity, non-null excess pwp (Eq. 40 + linear pwp)
//                    -> glisX_input-two-phase-hyperbolicity-pwp.json
//    5  Sec. 4.1.3  Dam break over a flat plane, Riemann problem,
//                    infinite permeability (Eq. 42)
//                    -> glisX_input-two-phase-dambreak-riemann.json
//    6  Sec. 4.1.3  Wet-dry dam break, 2D square release
//                    h = 10 on |x-L/2| <= L/10 and |y-H/2| <= H/10 else 0
//                    -> glisX_input-two-phase-dambreak-wetdry-2d.json
//    7  Sec. 4.1.3  Wet-dry dam break, 1D release  (h = 10 for x <= 10 else 0)
//                    -> glisX_input-two-phase-dambreak-wetdry-1d.json
//    8  Sec. 4.1.4  Efficiency test, discontinuous topography (Eq. 43),
//                    1D release  (h = 10 for x <= 10 else 0)
//                    -> glisX_input-two-phase-efficiency.json
//
//  Notes:
//   * Tests 1,2,6,7,8 keep the constant porosity  n_w = (rho_s - rho)/(rho_s - rho_w)
//     from the JSON densities; tests 3,4,5 impose a piecewise porosity consistent
//     with the reference q(x,0).
//   * The phase initial data U_w, U_s are the depth-averaged *mass fluxes*
//     (U = h v), matching the conserved variables of the paper.
//   * Whether the excess pwp is actually initialised/solved is still governed by
//     the JSON flag "do you want the pore water pressure?"; tests 4, 6 and 8 need
//     it set to true.
// ============================================================================

// Edit the default here, or override from CMake with -DTWO_PHASE_TEST=<n>
// (cmake -S . -B build -DTWO_PHASE_TEST=3).
#ifndef TWO_PHASE_TEST
#define TWO_PHASE_TEST 0
#endif

#endif // TWO_PHASE_MACRO_INPUT_H

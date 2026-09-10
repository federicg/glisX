/*
  Copyright (C) 2020 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>
 
#include <bim_distributed_vector.h>
#include <bim_timing.h>
#include <mumps_class.h> 
#include <tmesh.h>
#include <quad_operators.h>
#include <fstream>

#include "json.hpp"
#include "Taylor_Galerkin_IMEX-RKC_Strang_balanced.h"
#include "single_phase_macro_input.h"

using json = nlohmann::json;

// mpirun -np 4 main_TG2IMEXRKC glisX_input_tg2RKC_Bindo.json >out_TG2RKC.txt
// mpirun -np 1 main_TG2IMEXRKC $PWD inputs/dem_riemann.octbin.gz inputs/mask_in_vladi.octbin.gz 
// mpirun -np 1 main_TG2IMEXRKC $PWD inputs/dem_acheron.octbin.gz inputs/mask_in_acheron.octbin.gz

static constexpr char VARNAME_1[255] = "dem"; 
static constexpr char VARNAME_2[255] = "mask_in";  


static std::vector<double> dem;
static std::vector<double> basin_mask;

double h_min, L, H, res;
int NUM_REFINEMENTS, Nx, Ny;


// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
hanging_refinement (tmesh::quadrant_iterator quadrant)
{
  
  double x_minus, x_plus, y_minus, y_plus;
  x_minus = quadrant->p (0, 0);
  x_plus  = quadrant->p (0, 1);
  y_minus = quadrant->p (1, 0);
  y_plus  = quadrant->p (1, 2);
  
  double x_center, y_center;
  x_center = quadrant->centroid (0);
  y_center = quadrant->centroid (1);
  
  const auto marker = x_center<3./4.*L && x_center>L/4. && y_center<3./4.*H && y_center>H/4.;
//  const auto marker = (x_minus+x_plus)/2<L/2;//(x_minus+x_plus)/2<L/2 && (y_minus+y_plus)/2>H/2 ? 1 : 0;
  return marker; 
}


static int
raster_2_vector(const double& i_x,
                const double& i_y)
{
  static int ii;
  ii = i_y + Ny*i_x;
  
  return(ii);
}

static std::array<int,3>
global_coord_2_raster(const double& x,
                      const double& y)
{
  static double i_x;
  static double i_y;
  static int ii;
  
  
  // nearest neighbor
  i_x =   std::round(x / res);
  i_y = - std::round(y / res) + (Ny-1);

  ii = raster_2_vector(i_x,i_y);
  
  return(std::array<int,3>{{ ii,int(i_x),int(i_y) }});
}


static double
raster_value(const double& x,
             const double& y,
             std::vector<double> DD)
{
  // bilinear interp.
  const double ix = x/res;
  const double iy = -y/res + (Ny-1);

  const std::array<double,2> ix_q = {std::floor(ix), std::ceil (ix)}; 
  const std::array<double,2> iy_q = {std::floor(iy), std::ceil (iy)};

  const auto Dx_adi = ix_q[1]-ix_q[0];
  const auto Dy_adi = iy_q[1]-iy_q[0];

  std::array<double,4> gamma = {0,0,0,0};
  for (int i=0; i<2; i++)
  {
    for (int j=0; j<2; j++)
    {
      gamma[i+j*2] = DD[raster_2_vector(ix_q[(i+1)%2],iy_q[(j+1)%2])];

      gamma[i+j*2] *= Dx_adi!= 0 ? std::abs(ix_q[i]-ix)/Dx_adi : .5;
      gamma[i+j*2] *= Dy_adi!= 0 ? std::abs(iy_q[j]-iy)/Dy_adi : .5;
    }
  }
  
  return(gamma[0]+gamma[1]+gamma[2]+gamma[3]);
}

// ---------------------------------------------------------------------------
//  Initial condition and geometry, selected at compile time through
//  SINGLE_PHASE_TEST (see single_phase_macro_input.h).  dem_fun returns the
//  topography Z; h0_fun the material height H; Ux0_fun / Uy0_fun the mass
//  fluxes U = H v.
// ---------------------------------------------------------------------------

double dem_fun (const double& xx, const double& yy)
{
#if   SINGLE_PHASE_TEST == 2                    // Sec. 4.1.1 - smooth bump, Eq. (31)
  return 1. + (1./10.) * std::exp (-50./(L*L) * (xx - L/2.) * (xx - L/2.));
#elif SINGLE_PHASE_TEST == 3                    // Sec. 4.1.2 - smooth topography, Eq. (33)
  return 5. * std::exp (-2./5. * (xx - 5.) * (xx - 5.));
#elif SINGLE_PHASE_TEST == 4                    // Sec. 4.1.2 - discontinuous topography, Eq. (34)
  return (xx >= 4. && xx <= 8.) ? 4. : 0.;
#elif SINGLE_PHASE_TEST == 7                    // Sec. 4.2.3 - ~22 deg inclined plane
  return (L - xx) * std::tan (22. * M_PI / 180.);
#elif SINGLE_PHASE_TEST >= 1 && SINGLE_PHASE_TEST <= 6   // flat bottom
  return 0.;
#else                                           // SINGLE_PHASE_TEST == 0 - raster
  return raster_value (xx, yy, dem);
#endif
}



using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = distributed_vector; //distributed_vector; //std::vector<double>;         // Typedef for local q_0 vector // distributed_vector

double h0_fun (const double& xx, const double& yy)
{
#if   SINGLE_PHASE_TEST == 1                    // Sec. 4.1.1 - viscous dam break, H = 1
  return 1.;
#elif SINGLE_PHASE_TEST == 2                    // Sec. 4.1.1 - smooth solution, H(x,0) = Z(x), Eq. (32)
  return dem_fun (xx, yy);
#elif SINGLE_PHASE_TEST == 3 || SINGLE_PHASE_TEST == 4   // Sec. 4.1.2 - lake at rest, free surface at 10 m, Eq. (35)
  return 10. - dem_fun (xx, yy);
#elif SINGLE_PHASE_TEST == 5 || SINGLE_PHASE_TEST == 6   // Sec. 4.2.1 / 4.2.2 - radial dam break, Eq. (36)
  return (std::sqrt ((xx - L/2.) * (xx - L/2.) + (yy - L/2.) * (yy - L/2.)) <= 0.5) ? 2. : 1.;
#elif SINGLE_PHASE_TEST == 7                    // Sec. 4.2.3 - granular slide, Eq. (37)
  {
    const double amp = 0.2 + 0.01 * std::sin (10. * (yy - L/2.) / L * M_PI / L);
    const bool in_V = ((xx - L/2.) * (xx - L/2.) / (L*L)
                     + (yy - L/2.) * (yy - L/2.) / (L*L)) <= amp * amp;
    return in_V ? std::max (0., std::min (500. * xx / L - 200., 30.)) : 0.;
  }
#else                                           // SINGLE_PHASE_TEST == 0 - raster release, H0 = 38 m
  return (basin_mask[global_coord_2_raster (xx, yy)[0]] == 1 ? 38. : 0.);
#endif
}
double Ux0_fun (double xx, double yy)
{
#if SINGLE_PHASE_TEST == 1                      // Sec. 4.1.1 - Eq. (30)
  return (xx <= L/2.) ? 1. : 0.5;
#else
  return 0.;
#endif
}
double Uy0_fun (double xx, double yy)
{
  return 0.;
}


// Assemble vector from mesh.
// FIXME  the following two functions are copied over from
// "quad_operators.cpp" as they were not exported in an header,
// should find better way to avoid code duplication
void
assemble_vector (tmesh::quadrant_iterator& quadrant,
                 const std::array<double, 4>& locrhs,
                 Q1& rhs,
                 const ordering& ord = default_ord)
{
  
  std::vector<unsigned int> rows;
  rows.reserve (2);
  int i, r;
  
  for (i = 0; i < 4; ++i)
  {
    rows.clear ();
    
    if (! quadrant->is_hanging (i))
      rows.push_back (quadrant->gt (i));
    else
    {
      rows.push_back (quadrant->gparent (0, i));
      rows.push_back (quadrant->gparent (1, i));
    }
    
    for (r = 0; r < rows.size (); ++r)
    rhs[ord (rows[r])] +=
    locrhs[i] / rows.size ();
  }
}



template <class T>
void
quadrant_marker_list (tmesh::quadrant_iterator& q, 
                      const T& only_h, const T& only_Ux, const T& only_Uy, const double& dt, std::set<int>& output)
{

  std::array<double,4> h_current = {0,0,0,0};
  std::array<double,2> vel       = {0,0};
  for (int ii = 0; ii < 4; ++ii)
    {
      if (! q->is_hanging (ii)){
        const auto & h_candidate = only_h[q->gt (ii)];
        h_current[ii] = h_candidate;

        vel[0] += h_candidate>h_min ? only_Ux[q->gt (ii)]/h_candidate : 0.;
        vel[1] += h_candidate>h_min ? only_Uy[q->gt (ii)]/h_candidate : 0.;
      }
      else
      {
        const auto & h_candidate = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
        h_current[ii] = h_candidate;

        vel[0] += h_candidate>h_min ? .5*(only_Ux[q->gparent (0,ii)] + only_Ux[q->gparent (1,ii)])/h_candidate : 0.;

        vel[1] += h_candidate>h_min ? .5*(only_Uy[q->gparent (0,ii)] + only_Uy[q->gparent (1,ii)])/h_candidate : 0.;
      }
    }
  vel[0] /= 4.;
  vel[1] /= 4.;

  const bool basin_check = ((h_current[0]+h_current[1]+h_current[2]+h_current[3])> h_min && 
                            (h_current[0]*h_current[1]*h_current[2]*h_current[3])<=h_min) ? true : false;


  if (basin_check)
  {
    std::vector<std::tuple<tmesh::quadrant_iterator, int> > quadrant_list;
    output.insert(q->get_global_quad_idx ());
    for (auto quadrant_nei  = q->begin_neighbor_sweep();
              quadrant_nei != q->end_neighbor_sweep (); ++quadrant_nei)
    {
      quadrant_list.push_back(std::tuple<tmesh::quadrant_iterator, int>{quadrant_nei, quadrant_nei->get_global_quad_idx ()});
    }

    const double xx_ini = q->centroid (0);
    const double yy_ini = q->centroid (1);

    const double xx_fin = xx_ini + vel[0]*dt;
    const double yy_fin = yy_ini + vel[1]*dt;

    for (auto & qq : quadrant_list)
    {
      auto & quadrant = std::get<0>(qq);

      //std::cout << quadrant->p(0,0) << " " << std::get<1>(qq) << std::endl;

      const double x1 = quadrant->p(0,0);
      const double x2 = quadrant->p(0,1);
      const double y1 = quadrant->p(1,0);
      const double y2 = quadrant->p(1,2); 

      //std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;

      const bool cond1 = ( yy_fin - y1)>=0;
      const bool cond2 = (-yy_fin + y2)>=0;
      const bool cond3 = ( xx_fin - x1)>=0;
      const bool cond4 = (-xx_fin + x2)>=0;


      if (cond1 && cond2 && cond3 && cond4) // the point is internal to the current quadrant
      {
        output.insert(std::get<1>(qq));
      }

    }
  }

}






// Re-Define tic and toc to add an MPI_Barrier
#define TIC()  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }


int
main (int argc, char **argv)
{

  // parse input file,
  MPI_Init (&argc, &argv);
  std::ifstream input_file(argv[1]);
  json input_data = json::parse(input_file);

                 res                                      = input_data["raster resolution"];
                 Nx                                       = input_data["number raster columns"];
                 Ny                                       = input_data["number raster rows"];
                 NUM_REFINEMENTS                          = input_data["initial level of refinement"];
  const double & REDCDT                                   = input_data["CFL condition"];
  const double & T                                        = input_data["final time in seconds"];
  const double & SPACE_ADAPTDT                            = input_data["space adaptation procedure interval in seconds"];
  const double & SAVEDT                                   = input_data["saving interval in seconds"];
  const double & DELTAT                                   = input_data["maximum time step allowed"];
  const double & mesh_size_dry                            = input_data["desired resolution in meters of the mesh size in dry regions"];
  const double & mesh_size_wet                            = input_data["minimum resolution in meters of the mesh size in wet regions"];
  const double & mesh_size_interface                      = input_data["desired resolution in meters of the mesh size in wet-dry interface regions"];
  const bool   & is_time_adaptivity                       = input_data["do you want the time step predictor?"];
  const bool   & is_initial_refinement                    = input_data["do you want to refine the mesh initially?"];
  const bool   & is_space_adaptivity                      = input_data["do you want the space adaptation with interface tracking?"];
  const bool   & is_non_reflBC                            = input_data["do you want non reflecting BC?"];
  const bool   & is_bed_friction                          = input_data["do you want the bed friction?"];
  const bool   & is_stress_tensor                         = input_data["do you want the stress tensor?"];
  const bool   & is_max_time_step_from_CFL                = input_data["do you want the maximum time step given by CFL condition for the transport term?"];
                 h_min                                    = input_data["minimum material height threshold"];
  const double & grav                                     = input_data["gravitational field"];
  const double & density                                  = input_data["material density"];
  const double & turbulence_coeff                         = input_data["turbulence coefficient"];
  const double & surface_pressure                         = input_data["surface atmospheric pressure"];
        double   bed_friction_angle_rad                   = input_data["bed friction angle in degrees"];
                 bed_friction_angle_rad                   *= M_PI/180;
  const double & fluid_viscosity                          = input_data["fluid dynamic viscosity"];
  const double & yield_shear_stress                       = input_data["yield shear stress"];
  const double & tolerance_space_adapt                    = input_data["tolerance space adaptation"];

  const std::string & SAVE_DIR    = input_data["home saving directory, i.e., where we can find the directory results"];
#if SINGLE_PHASE_TEST == 0
  const std::string & DEM_DIR     = input_data["dem file, complete path"];
  const std::string & MASK_DIR    = input_data["mask file, complete path"];
#endif

  L = res*(Nx-1);
  H = res*(Ny-1);


  // Connectivity of local element
  constexpr p4est_topidx_t simple_conn_num_vertices = 4;
  constexpr p4est_topidx_t simple_conn_num_trees = 1;
  const double simple_conn_p[simple_conn_num_vertices*2] =
  {0,  0,
  0,  H,
  L,  0,
  L,  H};

  const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {  1,    3,    4,    2,    1 };



  // Management of solutions ordering
  ordering ordh  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 0> (gt); };
  ordering ordUx = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 1> (gt); };
  ordering ordUy = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 2> (gt); };

  
  
  // Initialize MPI
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);


  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  
  
  TIC ();
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  
  // for (int ii=0; ii<1; ii++)
  // {
  //   tmsh.set_refine_marker (hanging_refinement);
  //   tmsh.refine (recursive, 1);
  // }
  TOC ("Uniform refinement");
  
  
  
  // ln_nodes sono i dof non gli hanging node!! (sono esclusi dal calcolo)
  tmesh::idx_t gn_nodes    = tmsh.num_global_nodes (); // Return total number of nodes owned by all process
  tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes (); // Return number of nodes owned by local process
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();  // Return number of quadrants owned by local process across all trees
  tmesh::idx_t gn_elements = tmsh.num_global_quadrants (); // Return number of quadrants owned by all processes across all trees
  

  /// Allocate initial data container
  Q1 sol  (ln_nodes * 3);
  Q1 incr (ln_nodes * 3);
  sol.get_owned_data  ().assign (sol.get_owned_data  ().size (), 0.0);
  incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);

  
  Q1 mass (ln_nodes * 3);
  bim2a_mass_vector (tmsh, mass, ordh);
  bim2a_mass_vector (tmsh, mass, ordUx);
  bim2a_mass_vector (tmsh, mass, ordUy);
  mass.assemble ();
  
  Q0 sol_onehalf (ln_elements * 3);
  sol_onehalf.get_owned_data  ().assign (sol_onehalf.get_owned_data  ().size (), 0.0);

  Q0 Z_onehalf (ln_elements);
  Z_onehalf.get_owned_data  ().assign (Z_onehalf.get_owned_data  ().size (), 0.0);

  std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 3);
  
  Q1 Z (ln_nodes);
  Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);


  Q1 Newton_it (ln_nodes);
  Newton_it.get_owned_data ().assign (Newton_it.get_owned_data ().size (), 0.0);

  // Q1 mask_fin (ln_nodes);
  // mask_fin.get_owned_data ().assign (mask_fin.get_owned_data ().size (), 0.0);


  std::string str = ""; 
  char filename[255]="", arr[255]="";

#if SINGLE_PHASE_TEST == 0
  TIC();
  str = std::string(DEM_DIR);
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);


  octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
  octave_value v;

  octave_io_open (filename, m_in, &m_out);
  octave_load (VARNAME_1, v);
  Matrix M = v.matrix_value ();
  dem.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), dem.begin ());

  str = std::string(MASK_DIR);
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);

  octave_io_open (filename, m_in, &m_out);
  octave_load (VARNAME_2, v);
  M = v.matrix_value ();
  basin_mask.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), basin_mask.begin ());
  TOC("Load data matrix");
#else
  (void) VARNAME_1; (void) VARNAME_2;
#endif



  // Initialize 
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
  {
    double xx_c=quadrant->centroid(0);
    double yy_c=quadrant->centroid(1); 
    

    for (int ii = 0; ii < 4; ++ii)
    {
      if (! quadrant->is_hanging (ii)){
        double xx=quadrant->p(0,ii);
        double yy=quadrant->p(1,ii); 
        
        sol [ordh     (quadrant->gt (ii))] = h0_fun  (xx, yy);
        sol [ordUx    (quadrant->gt (ii))] = Ux0_fun (xx, yy);
        sol [ordUy    (quadrant->gt (ii))] = Uy0_fun (xx, yy);
        
        Z           [quadrant->gt (ii)] = dem_fun(xx,yy); 
	      Newton_it   [quadrant->gt (ii)] = 0.;
      }
      
      else
      {
        // touch parent nodes to set up distributed vector structure
        sol [ordh   (quadrant->gparent(0,ii))] += 0.;
        sol [ordh   (quadrant->gparent(1,ii))] += 0.;
        sol [ordUx  (quadrant->gparent(0,ii))] += 0.;
        sol [ordUx  (quadrant->gparent(1,ii))] += 0.;
        sol [ordUy  (quadrant->gparent(0,ii))] += 0.;
        sol [ordUy  (quadrant->gparent(1,ii))] += 0.;
        
        Z   [quadrant->gparent(0,ii)] += 0.;
        Z   [quadrant->gparent(1,ii)] += 0.;

	      Newton_it   [quadrant->gparent(0,ii)] += 0.;
        Newton_it   [quadrant->gparent(1,ii)] += 0.;
      }
    }
  }


  // bim2a_solution_with_ghosts in quad_operators.cpp
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);
  
  bim2a_solution_with_ghosts (tmsh, Z, replace_op);

  bim2a_solution_with_ghosts (tmsh, Newton_it, replace_op);

  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy);

  

  if (is_initial_refinement)
  {
    
    TIC();
    Q1 only_h (ln_nodes);
    bim2a_solution_with_ghosts (tmsh, only_h);
    for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
    {
      only_h(idx) = sol(ordh(idx));
    }
    only_h.assemble (replace_op);
    TOC("get separated sol.");

  
  
   
    TIC();
    auto dh = bim2c_quadtree_pde_recovered_gradient (tmsh, only_h);
    //q2_vec h_star = bim2c_quadtree_pde_recovered_solution (tmsh, only_h, dh);
    TOC ("gradient and hstar");


    TIC();
    // auto estimator = [& h_star, & only_h] (tmesh::quadrant_iterator q)
    // {
    //   return estimator_sol (q, h_star, only_h);
    // };
    auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
    {
      return estimator_grad(q, dh, only_h);
    };
    auto estimator_flux = [& only_h] (tmesh::quadrant_iterator q)
    {

      std::array<double,4> h_mesh = {0,0,0,0};
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! q->is_hanging (ii)){
          h_mesh[ii] = only_h[q->gt (ii)];
        }
        else
        {
          h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
        }
      }

      const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 

      return (basin_check); 
    };


    auto dry_function = [& only_h] (tmesh::quadrant_iterator q)
    {

      std::array<double,4> h_mesh = {0,0,0,0};
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! q->is_hanging (ii)){
          h_mesh[ii] = only_h[q->gt (ii)];
        }
        else
        {
          h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
        }
      }

      const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])==0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 

      return (basin_check); 
    };


    tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, tolerance_space_adapt, 6, 0, 0);
    //tmsh.set_metrics_marker (estimator, 1e-5, 4, 3, 1);
    tmsh.metrics_refine (1e7);  // RAFFINAMENTO (arg is max element)

    // tmsh.set_coarsen_marker (coarsen_function);
    // tmsh.set_refine_marker  (refine_function);
    // tmsh.coarsen (recursive, 1, 0);
    // tmsh.refine  (recursive, 1);
    TOC ("refine");
  
    // Ottengo i parametri della mesh corrente
    TIC();
    gn_nodes    = tmsh.num_global_nodes ();
    ln_nodes    = tmsh.num_owned_nodes ();
    ln_elements = tmsh.num_local_quadrants ();
    gn_elements = tmsh.num_global_quadrants ();
    TOC ("Obtaining new parameters");
  
  
    // Interpolo sol sulla nuova mesh
    TIC();
  
    Q1 incr_ (ln_nodes * 3);
    incr_.get_owned_data ().assign (incr_.get_owned_data ().size(), 0.0);
    incr_.assemble ();
  
    Q1 mass_ (ln_nodes * 3);
    bim2a_mass_vector (tmsh, mass_, ordh );
    bim2a_mass_vector (tmsh, mass_, ordUx);
    bim2a_mass_vector (tmsh, mass_, ordUy);
    mass_.assemble ();
  
    Q0 sol_onehalf_ (ln_elements * 3);
    sol_onehalf_.get_owned_data ().assign (sol_onehalf_.get_owned_data ().size(), 0.0);

    Q0 Z_onehalf_ (ln_elements);
    Z_onehalf_.get_owned_data ().assign (Z_onehalf_.get_owned_data ().size(), 0.0);


    std::vector<std::array<double,4>> incr_anti_diff_ (ln_elements * 3);


    Q1 sol_ (ln_nodes * 3);
    Q1 Z_ (ln_nodes);
    // Q1 mask_fin_ (ln_nodes);
    Q1 Newton_it_ (ln_nodes);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {
      double xx_c=quadrant->centroid(0);
      double yy_c=quadrant->centroid(1); 

      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii)){
          double xx=quadrant->p(0,ii);
          double yy=quadrant->p(1,ii);
          
          sol_ [ordh     (quadrant->gt (ii))] = h0_fun  (xx, yy);
          sol_ [ordUx    (quadrant->gt (ii))] = Ux0_fun (xx, yy);
          sol_ [ordUy    (quadrant->gt (ii))] = Uy0_fun (xx, yy);

          Z_[quadrant->gt (ii)] = dem_fun(xx,yy);

	        Newton_it_[quadrant->gt (ii)] = 0.;
        }
        
        else
        {
          sol_ [ordh   (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordh   (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordUx  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUx  (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordUy  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUy  (quadrant->gparent(1,ii))] += 0.;

          Z_[quadrant->gparent(0,ii)] += 0.;
          Z_[quadrant->gparent(1,ii)] += 0.;
	  
	        Newton_it_[quadrant->gparent(0,ii)] += 0.;
          Newton_it_[quadrant->gparent(1,ii)] += 0.;
        }
      }
    }

    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordh,  false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUx, false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUy);
    
    bim2a_solution_with_ghosts (tmsh, Z_, replace_op);

    bim2a_solution_with_ghosts (tmsh, Newton_it_, replace_op);
    
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordh,  false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUx, false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUy);


    sol                 = sol_;
    incr                = incr_;
    incr_anti_diff      = incr_anti_diff_;
    mass                = mass_;
    sol_onehalf         = sol_onehalf_;
    Z                   = Z_;
    Z_onehalf           = Z_onehalf_;
    Newton_it 		      = Newton_it_;
  
    TOC ("compute initial condition");
  }
  
  Q1 sol_dyn                 = sol;
  Q1 sold_dyn                = sol;
  Q1 soldd_dyn               = sol;
  Q1 sold_rkc_dyn            = sol;
  Q1 soldd_rkc_dyn           = sol;
  Q1 sol_ini_rkc_dyn         = sol;
  Q1 incr_dyn                = incr;
  Q1 incr_initial_source_dyn = incr;
  Q1 incr_source_dyn         = incr;
  Q1 stress_initial_step_dyn = incr;
  Q1 stress_step_dyn         = incr;
  Q1 P_plus_dyn              = incr;
  Q1 P_minus_dyn             = incr;  
  Q1 spec_radius_nodal_dyn   = incr;
  Q1 mass_dyn                = mass;
  Q1 Z_dyn                   = Z;
  Q1 Newton_it_dyn 	         = Newton_it;
  Q0 sol_onehalf_dyn         = sol_onehalf;
  Q0 Z_onehalf_dyn           = Z_onehalf;

  std::vector<std::array<double,4>> incr_anti_diff_dyn = incr_anti_diff;

  
  TG2_scheme stp(sol_dyn, 
                 sold_dyn, 
                 soldd_dyn, 
                 sold_rkc_dyn,
                 soldd_rkc_dyn,
                 sol_ini_rkc_dyn,
                 incr_dyn,
                 incr_initial_source_dyn,
                 incr_source_dyn, 
                 incr_anti_diff_dyn,
                 stress_initial_step_dyn,
                 stress_step_dyn,
                 P_plus_dyn, 
                 P_minus_dyn, 
                 spec_radius_nodal_dyn,
                 sol_onehalf_dyn, 
                 mass_dyn,
                 ordh, ordUx, ordUy, 
                 Z_dyn,
                 Z_onehalf_dyn,
		             Newton_it_dyn, 
                 DELTAT, h_min, is_non_reflBC, is_bed_friction, is_stress_tensor, grav,
                 density, turbulence_coeff, surface_pressure, bed_friction_angle_rad, fluid_viscosity, yield_shear_stress);
  
  
  // Save initial conditions

  str = std::string(SAVE_DIR) + "/results/swe_h_%4.4d"; 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordh);

  str = std::string(SAVE_DIR) + "/results/swe_Ux_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordUx); 

  str = std::string(SAVE_DIR) + "/results/swe_Uy_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordUy);
  
  str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0); 
  tmsh.octbin_export (filename, Z_dyn);

  str = std::string(SAVE_DIR) + "/results/swe_Newton_it_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, Newton_it_dyn);
  


  std::vector<double> full_time_vector;
  full_time_vector.reserve (static_cast<int> (T/DELTAT));
  std::vector<double> save_time_vector;
  save_time_vector.reserve (static_cast<int> (T/SAVEDT));
  std::vector<double> RKC_steps;
  RKC_steps.reserve (static_cast<int> (T/DELTAT));
  
  // Time loop
  double time      = 0.0;
  double time_old  = 0.0;
  double time_oldd = 0.0;

  
  stp.set_dt (DELTAT);
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
  {
    stp.compute_dt(quadrant);
  }
  double max_dt = REDCDT * stp.dt;
  MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&max_dt), 1, MPI_DOUBLE, MPI_MIN, tmsh.comm);
  stp.set_dt(max_dt);
  stp.set_old_dt(0.);
  time_old  -= stp.dt;
  time_oldd -= 2*stp.dt;
  stp.set_times(time, time_old, time_oldd);
  
  
  if(rank==0) {
    full_time_vector.push_back (0.0);
    save_time_vector.push_back (0.0);
  }
  int count = 0;
  
  double savecount = 0.0, space_adapt_count = 0.0;
  
  if (rank == 0)
  {
    std::cout << "start loop" << std::endl;
  }

  int counter_savings = 0, tot_number_savings = std::round(T/SAVEDT);
  
  TIC();
  while (counter_savings != tot_number_savings)
  {
    
    
    // Reset increment, and limiter terms
    //TIC();
    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);

    P_plus_dyn.get_owned_data ().assign (P_plus_dyn.get_owned_data ().size (), 0.0);
    P_plus_dyn.assemble (replace_op);

    P_minus_dyn.get_owned_data ().assign (P_minus_dyn.get_owned_data ().size (), 0.0);
    P_minus_dyn.assemble (replace_op);

    stress_initial_step_dyn.get_owned_data ().assign (stress_initial_step_dyn.get_owned_data ().size (), 0.0);
    stress_initial_step_dyn.assemble (replace_op);

    spec_radius_nodal_dyn.get_owned_data ().assign (spec_radius_nodal_dyn.get_owned_data ().size (), 0.0);
    spec_radius_nodal_dyn.assemble (replace_op);
    //TOC("Reset");
    //TIC();
    
  
    // compute time step, 
    stp.Fr = 0.;
    stp.set_dt (DELTAT);
    if (is_max_time_step_from_CFL)
    {
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
      {
        stp.compute_dt(quadrant);
      }
    }
    max_dt = REDCDT * stp.dt;

    stp.g_coeff = 1./(1.-stp.Fr*stp.Fr);

    stp.set_dt(max_dt); // deltat max
    MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&stp.dt), 1, MPI_DOUBLE, MPI_MIN, tmsh.comm);
 
    // Print current time
    if(rank==0)
    {
      std::cout << "MAXIMUM TIME STEP = " << stp.dt << std::endl;
    }   

    // time adaptivity
    if (is_time_adaptivity)
    { 
    
      stp.nu_htot = 0.;
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
        quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
      {
        stp.compute_dt_adaptive(quadrant);
      }
      MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&stp.nu_htot), 1, MPI_DOUBLE, MPI_SUM, tmsh.comm);

      const double local_estimator_time_tolerance = 1e-5;//5e-3*(stp.time-stp.timed)*std::sqrt(stp.time-stp.timed)/std::sqrt(stp.nu_htot);

      const double candidate_dt = local_estimator_time_tolerance/std::sqrt(stp.nu_htot)*(stp.time-stp.timed);
      stp.set_dt( (stp.nu_htot>0 && candidate_dt<stp.dt) ? candidate_dt : stp.dt );
      //stp.set_dt( (stp.nu_htot>0 && candidate_dt<stp.dt) ? std::max(candidate_dt, stp.dt/2.) : stp.dt ); // we set a minimum dt (std::max(,)), stp.dt/2.
    }

    

    if (stp.dt == 0 && rank == 0)
    {
      std::cout << "dt has gone to zero, sorry, STOP!" << std::endl;
      exit( -1. );
    }

    // check save with given frequency
    //std::cout << (savecount) << " " << stp.dt << " " << (savecount+stp.dt)/SAVEDT << " " << (stp.dt - std::fmod(savecount+stp.dt,SAVEDT) - SAVEDT*(std::floor(savecount+stp.dt/SAVEDT)-1)) << std::endl;
    // stp.set_dt((savecount+stp.dt)/SAVEDT>1 ? (stp.dt - std::fmod(savecount+stp.dt,SAVEDT) - SAVEDT*(std::floor(savecount+stp.dt/SAVEDT)-1))-SAVEDT : stp.dt);
    //stp.set_dt((savecount+stp.dt)/SAVEDT>1 ? (stp.dt - std::fmod(savecount+stp.dt,SAVEDT) - SAVEDT*(std::floor(savecount+stp.dt/SAVEDT)-1)) : stp.dt);
    stp.set_dt((savecount+stp.dt)/SAVEDT>1 ? SAVEDT-savecount : stp.dt);
    //stp.set_dt((time+stp.dt)>T ? T-(time+stp.dt) : stp.dt);



    time_oldd = time_old;
    time_old = time;
    time += stp.dt; 
    savecount += stp.dt;
    space_adapt_count += stp.dt;

    // credo questi qui non servano
    // MPI_Bcast (static_cast<void*> (&time),              1, MPI_DOUBLE, 0, tmsh.comm);
    // MPI_Bcast (static_cast<void*> (&savecount),         1, MPI_DOUBLE, 0, tmsh.comm);
    // MPI_Bcast (static_cast<void*> (&space_adapt_count), 1, MPI_DOUBLE, 0, tmsh.comm);
    // MPI_Barrier (tmsh.comm); // tmsh.comm = MPI_COMM_WORLD
    
    // Print current time
    if(rank==0) 
    {
      std::cout << "TIME = " << time << ", dt = " << stp.dt << std::endl;
      full_time_vector.push_back (time);
    }
    
    //std::cout << (time<T) << " " << time << " " << T << " " << time-T << std::endl;
    
    
    
    // first step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.first_step(quadrant);
    }
    


    // 
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.compute_nodal_anti_diffusive_fluxes(quadrant);
    }
    incr_dyn.assemble ();
    P_plus_dyn.assemble ();
    P_minus_dyn.assemble ();




    stp.set_times(time, time_old, time_oldd);
    soldd_dyn = sold_dyn;
    sold_dyn  = sol_dyn;
    

    // low order solution
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
    {
      sol_dyn.get_owned_data ()[kk] += (stp.dt + stp.dt_old)*.5*incr_dyn.get_owned_data ()[kk]/mass_dyn.get_owned_data ()[kk];
    }
    sol_dyn.assemble(replace_op); 



    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii) && sol_dyn [ordh    (quadrant->gt (ii))]<0){
          sol_dyn [ordh    (quadrant->gt (ii))] = 0.; //h_min; //0.;
        }
      }
    }
    //sol_dyn.assemble (replace_op);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordh,  false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUx, false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUy);

    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);


    
    // second order correction
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.second_step(quadrant);
    }
    incr_dyn.assemble ();
    //TOC("Compute step");

    
    

    //TIC();
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
    {
      sol_dyn.get_owned_data ()[kk] += (stp.dt + stp.dt_old)*.5*incr_dyn.get_owned_data ()[kk] / mass_dyn.get_owned_data ()[kk];
    }

    
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii) && sol_dyn [ordh (quadrant->gt (ii))] < 0){
          sol_dyn [ordh    (quadrant->gt (ii))] = 0.; //h_min; //0.;
        }
      }
    }
    sol_dyn.assemble (replace_op);
    //TOC("Apply increment");


    // Verwer IMEX-RKC
    sol_ini_rkc_dyn = sol_dyn; // copy
    soldd_rkc_dyn   = sol_dyn; // copy
    sold_rkc_dyn    = sol_dyn; // copy


    for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.compute_stress_slope(quadrant, true);
    }
    stress_initial_step_dyn.assemble();


    // for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=3)
    // {
    //   std::cout << stress_initial_step_dyn.get_owned_data ()[kk] << " " << stress_initial_step_dyn.get_owned_data ()[kk+1] << " " << stress_initial_step_dyn.get_owned_data ()[kk+2] << std::endl;
    //   //std::cout << sol_dyn.get_owned_data ()[kk] << " " << sol_dyn.get_owned_data ()[kk+1] << " " << sol_dyn.get_owned_data ()[kk+2] << std::endl;//" " << sol_dyn.get_owned_data ()[kk+1] << std::endl;
    // }
    // exit(1);




    for (int kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=3)
    {
      stp.loop_step(kk, true);
    }
    incr_initial_source_dyn.assemble (replace_op);



    // compute the spec_radius_nodal
    for (auto kk = 0; kk < spec_radius_nodal_dyn.get_owned_data ().size (); ++kk)
    {
      spec_radius_nodal_dyn.get_owned_data ()[kk] /= mass_dyn.get_owned_data ()[kk]; 
    }
    spec_radius_nodal_dyn.assemble(replace_op);


 
    double spec_radius = 0.;
    for (auto kk = 0; kk < spec_radius_nodal_dyn.get_owned_data ().size (); ++kk)
    {
      spec_radius = std::max(spec_radius, spec_radius_nodal_dyn.get_owned_data ()[kk]);
    }
    MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&spec_radius), 1, MPI_DOUBLE, MPI_MAX, tmsh.comm);


    double s = 1. + std::round(std::sqrt(1 + stp.dt*spec_radius/.653));//5;//std::round(std::max(std::sqrt(stp.dt*spec_radius/.653), 2.));


    if(rank==0)
    {
      RKC_steps.push_back(s);
    }


    if (rank==0) std::cout << "number of steps and spectral radius, " << s << " " << spec_radius << std::endl;

    // compute here the coefficients!, it is to prepare the following loop
    stp.prepare_IMEXRKC_coefficients(s);

    

    //double s = 1 + std::round(std::sqrt(1 + stp.dt*spec_radius/.653)); // # of stages minimum is 2!!! otherwise errors inside for the recursion!
    for (int jj = 1; jj <= s; jj++)
    {

      if (rank==0) std::cout << "current IMEX-RKC step, " << jj << std::endl;

      for (int kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=3)
      {
        stp.rkc(jj, s, kk);
      }
      sol_dyn.assemble(replace_op);
      


      soldd_rkc_dyn = sold_rkc_dyn; // copy
      sold_rkc_dyn  = sol_dyn;      // copy


      stress_step_dyn.get_owned_data ().assign (stress_step_dyn.get_owned_data ().size (), 0.0);
      stress_step_dyn.assemble (replace_op);


      for (auto quadrant = tmsh.begin_quadrant_sweep ();
        quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
      {
        stp.compute_stress_slope(quadrant, false);
      }
      stress_step_dyn.assemble();

      
      // // Possible internal update of the spectral radius //
      // spec_radius_nodal_dyn.get_owned_data ().assign (spec_radius_nodal_dyn.get_owned_data ().size (), 0.0);
      // spec_radius_nodal_dyn.assemble (replace_op);
      // for (auto kk = 0; kk < spec_radius_nodal_dyn.get_owned_data ().size (); ++kk)
      // {
      //   spec_radius_nodal_dyn.get_owned_data ()[kk] /= mass_dyn.get_owned_data ()[kk]; 
      // }
      // spec_radius_nodal_dyn.assemble(replace_op); 

      // spec_radius = 0.;
      // for (auto kk = 0; kk < spec_radius_nodal_dyn.get_owned_data ().size (); ++kk)
      // {
      //   spec_radius = std::max(spec_radius, spec_radius_nodal_dyn.get_owned_data ()[kk]);
      // }
      // MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&spec_radius), 1, MPI_DOUBLE, MPI_MAX, tmsh.comm);

      // int number_remaining_stages = std::round(std::max(std::sqrt(stp.dt*(1. - c_fun(j,s))*spec_radius/.653), 2.)); 

      // if (number_remaining_stages!=(s-jj))
      // {
      //   // o calcolare nuove condizioni iniziali 
      //   // o interpolare con polinomio di Hermite le soluzioni che già abbiamo ma fare comunque update di s!!
      // }
      // //----------------------------//
      

      for (int kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=3)
      {
        stp.loop_step(kk, false);
      }
      incr_source_dyn.assemble (replace_op);
        
    }
    Newton_it_dyn.assemble (replace_op);

    //return 0;

    

    stp.set_old_dt(stp.dt);
    stp.set_old_dt(0.);

    // first, Strang half step!

    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);

    P_plus_dyn.get_owned_data ().assign (P_plus_dyn.get_owned_data ().size (), 0.0);
    P_plus_dyn.assemble (replace_op);

    P_minus_dyn.get_owned_data ().assign (P_minus_dyn.get_owned_data ().size (), 0.0);
    P_minus_dyn.assemble (replace_op);

    // first step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
     quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.first_step(quadrant);
    }



    // 
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
     quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.compute_nodal_anti_diffusive_fluxes(quadrant);
    }
    incr_dyn.assemble ();
    P_plus_dyn.assemble ();
    P_minus_dyn.assemble ();


    // low order solution
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
    {
      sol_dyn.get_owned_data ()[kk  ] += (stp.dt + stp.dt_old)*.5*incr_dyn.get_owned_data ()[kk  ]/mass_dyn.get_owned_data ()[kk  ];
    }
    //sol_dyn.assemble(replace_op);
    


    for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep ();
      ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii) && sol_dyn [ordh    (quadrant->gt (ii))]<0)
        {
          sol_dyn [ordh    (quadrant->gt (ii))] = 0.; //h_min; //0.;
        }
      }
    }
    //sol_dyn.assemble (replace_op);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordh,  false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUx, false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUy);


    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);


    // second order correction
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
     quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.second_step(quadrant);
    }
    incr_dyn.assemble ();
    //TOC("Compute step");


    //TIC();
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
    {
      sol_dyn.get_owned_data ()[kk] += (stp.dt + stp.dt_old)*.5*incr_dyn.get_owned_data ()[kk] / mass_dyn.get_owned_data ()[kk];
    }


    for (auto quadrant = tmsh.begin_quadrant_sweep ();
     quadrant != tmsh.end_quadrant_sweep ();
     ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii) && sol_dyn [ordh (quadrant->gt (ii))] < 0)
        {
          sol_dyn [ordh    (quadrant->gt (ii))] = 0.; //h_min; //0.;
        }
      }
    }
    sol_dyn.assemble (replace_op);
    //TOC("Apply increment");


    // Save solution
    if ((savecount-SAVEDT) >= -std::numeric_limits<double>::epsilon()*SAVEDT) 
    {
      //TIC();
      if (rank == 0)
        std::cout << "savecount = " << savecount << std::endl;
      count++;
      save_time_vector.push_back (time);

     
      str = std::string(SAVE_DIR) + "/results/swe_h_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,   count);
      tmsh.octbin_export (filename, sol_dyn, ordh);
      
      str = std::string(SAVE_DIR) + "/results/swe_Ux_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_dyn, ordUx);
      
      str = std::string(SAVE_DIR) + "/results/swe_Uy_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_dyn, ordUy);
      
      str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, Z_dyn);

      str = std::string(SAVE_DIR) + "/results/swe_Newton_it_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, Newton_it_dyn);
      savecount = 0.0;
      //TOC("Exporting solution");

      counter_savings++;

    }



    
    if (is_space_adaptivity &&  ((space_adapt_count-SPACE_ADAPTDT) >= -std::numeric_limits<double>::epsilon()*SPACE_ADAPTDT))
    {

      //TIC();
      Q1 only_h (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_h);
      for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
      {
        only_h(idx) = sol_dyn(ordh(idx));
      }
      only_h.assemble (replace_op);

      Q1 only_Ux (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_Ux);
      for (auto idx = only_Ux.get_range_start (); idx != only_Ux.get_range_end (); ++idx)
      {
        only_Ux(idx) = sol_dyn(ordUx(idx));
      }
      only_Ux.assemble (replace_op);

      Q1 only_Uy (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_Uy);
      for (auto idx = only_Uy.get_range_start (); idx != only_Uy.get_range_end (); ++idx)
      {
        only_Uy(idx) = sol_dyn(ordUy(idx));
      }
      only_Uy.assemble (replace_op);
      //TOC("get separated sol.");

      //TIC();
      std::set<int> global_index_quad;
      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep ();
           ++q)
      {
        quadrant_marker_list(q, only_h,  only_Ux, only_Uy, stp.dt, global_index_quad);
      }      
      //TOC("front track.");  
      
      
      //TIC();
      auto dh = bim2c_quadtree_pde_recovered_gradient (tmsh, only_h);
      //q2_vec h_star = bim2c_quadtree_pde_recovered_solution (tmsh, only_h, dh);
      //TOC ("gradient and hstar");
      
      //TIC();
      // auto estimator = [& h_star, & only_h] (tmesh::quadrant_iterator q)
      // {
      //   return estimator_sol (q, h_star, only_h);
      // };
      auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
      {
        return estimator_grad(q, dh, only_h);
      };

      // auto estimator_flux = [& h_star, & only_h] (tmesh::quadrant_iterator q)
      // {
      //   std::array<double,4> h_mesh = {0,0,0,0};
      //   for (int ii = 0; ii < 4; ++ii)
      //   {
      //     if (! q->is_hanging (ii)){
      //       h_mesh[ii] = only_h[q->gt (ii)];
      //     }
      //     else
      //     {
      //       h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
      //     } 
      //   }

      //   const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 
      //   //const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>h_min && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])<=h_min) ? 1 : 0; 

      //   return (basin_check); 
      // };


      auto estimator_flux = [& global_index_quad] (tmesh::quadrant_iterator q)
      {
        if ( global_index_quad.find(q->get_global_quad_idx ()) != global_index_quad.end() )
        {
          return 1;
        }
        return 0;
      };


      auto dry_function = [& only_h] (tmesh::quadrant_iterator q)
      {

        std::array<double,4> h_mesh = {0,0,0,0};
        for (int ii = 0; ii < 4; ++ii)
        {
          if (! q->is_hanging (ii)){
            h_mesh[ii] = only_h[q->gt (ii)];
          }
          else
          {
            h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
          }
        }

        const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])<h_min && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])<h_min) ? 1 : 0; 
        //const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])==0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 

        return (basin_check); 
      };


      tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, tolerance_space_adapt, 6, 0, 0);
      //tmsh.set_metrics_marker (estimator, 1e-5, 4, 3, 1); 
      tmsh.metrics_refine (1e7);  // RAFFINAMENTO (arg is max element)

      // tmsh.set_coarsen_marker (coarsen_function);
      // tmsh.set_refine_marker  (refine_function);
      // tmsh.coarsen (recursive, 1, 0);
      // tmsh.refine  (recursive, 1);
      //TOC ("refine");

      // Ottengo i parametri della mesh corrente
      //TIC();
      gn_nodes    = tmsh.num_global_nodes ();
      ln_nodes    = tmsh.num_owned_nodes ();
      ln_elements = tmsh.num_local_quadrants ();
      gn_elements = tmsh.num_global_quadrants ();
      //TOC ("Obtaining new parameters");
      
      
      // Interpolo sol sulla nuova mesh
      //TIC();
      Q1 sol (ln_nodes * 3);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);
      interpolate_vector (tmsh, sol_dyn, sol, ordh);
      interpolate_vector (tmsh, sol_dyn, sol, ordUx);
      interpolate_vector (tmsh, sol_dyn, sol, ordUy);
      sol.assemble (replace_op);
      
      Q1 sold (ln_nodes * 3);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUy);
      interpolate_vector (tmsh, sold_dyn, sold, ordh);
      interpolate_vector (tmsh, sold_dyn, sold, ordUx);
      interpolate_vector (tmsh, sold_dyn, sold, ordUy);
      sold.assemble (replace_op);
      
      
      Q1 soldd (ln_nodes * 3);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUy);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordh );
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUy);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUx);
      soldd.assemble (replace_op);
      
      
      Q1 incr (ln_nodes * 3);
      incr.get_owned_data ().assign (incr.get_owned_data ().size(), 0.0);
      incr.assemble ();

      std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 3);

      
      Q1 mass (ln_nodes * 3);
      bim2a_mass_vector (tmsh, mass, ordh );
      bim2a_mass_vector (tmsh, mass, ordUx);
      bim2a_mass_vector (tmsh, mass, ordUy);
      mass.assemble ();
      
      Q0 sol_onehalf (ln_elements * 3);
      sol_onehalf.get_owned_data ().assign (sol_onehalf.get_owned_data ().size(), 0.0);


      Q0 Z_onehalf (ln_elements);
      Z_onehalf.get_owned_data ().assign (Z_onehalf.get_owned_data ().size(), 0.0);


      Q1 Z (ln_nodes);
      Q1 Newton_it (ln_nodes);
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
      {
        double xx_c=quadrant->centroid(0);
        double yy_c=quadrant->centroid(1); 
        
        for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            double xx=quadrant->p(0,ii);
            double yy=quadrant->p(1,ii);
            Z           [quadrant->gt (ii)] = dem_fun(xx,yy); 
	          Newton_it   [quadrant->gt (ii)] = 0.;
          }
           
          else
          {
            Z[quadrant->gparent(0,ii)] += 0.;
            Z[quadrant->gparent(1,ii)] += 0.;

	          Newton_it[quadrant->gparent(0,ii)] += 0.;
            Newton_it[quadrant->gparent(1,ii)] += 0.;
          }
        }
      }
      
      bim2a_solution_with_ghosts (tmsh, Newton_it, replace_op);
      
      bim2a_solution_with_ghosts (tmsh, Z, replace_op);
      
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy);
      
       
      sol_dyn                 = sol;
      sold_dyn                = sold;
      soldd_dyn               = soldd;
      sold_rkc_dyn            = soldd;
      soldd_rkc_dyn           = soldd;
      sol_ini_rkc_dyn         = soldd;
      incr_dyn                = incr;
      incr_source_dyn         = incr;
      incr_initial_source_dyn = incr;
      incr_anti_diff_dyn      = incr_anti_diff;
      stress_initial_step_dyn = incr;
      stress_step_dyn         = incr;
      P_plus_dyn              = incr;
      P_minus_dyn             = incr;
      spec_radius_nodal_dyn   = incr;
      mass_dyn                = mass;
      sol_onehalf_dyn         = sol_onehalf;
      Z_onehalf_dyn           = Z_onehalf;
      Z_dyn                   = Z;
      Newton_it_dyn 	        = Newton_it;	


      space_adapt_count = 0.0;

      
      //TOC ("Interpolation");
      
    }
      
    
    
    
  }
  
  
  if (rank == 0)
  {
    str = std::string(SAVE_DIR) + "/results/timesteps.octbin"; 
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);

    ColumnVector vtmp (save_time_vector.size ());
    std::copy (save_time_vector.begin (), save_time_vector.end (), vtmp.fortran_vec ());
    octave_io_mode m;
    octave_io_open (filename, gz_write_mode, &m);
    octave_save ("save_time", vtmp);
    vtmp.resize (full_time_vector.size ());
    std::copy (full_time_vector.begin (), full_time_vector.end (), vtmp.fortran_vec ());
    octave_save ("full_time", vtmp);
    vtmp.resize (RKC_steps.size ());
    std::copy (RKC_steps.begin (), RKC_steps.end (), vtmp.fortran_vec ());
    octave_save ("RKC_steps", vtmp);
    octave_io_close ();
  }
  
  TOC ("loop completed");
  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  MPI_Finalize ();
  return 0;
  
}







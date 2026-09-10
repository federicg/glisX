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
#include "Taylor_Galerkin_IMEX-RKC_Strang_balanced_two_phase.h"


using json = nlohmann::json;

// mpirun -np 4 main_TG2IMEXRKC2PHASE glisX_input_two-phase_wet-wet.json >out_two-phase.txt

// sqrt((Uxw+Uxs)/((hw+hs)*((hw+hs)>1e-2))*(Uxw+Uxs)/((hw+hs)*((hw+hs)>1e-2)) + (Uyw+Uys)/((hw+hs)*((hw+hs)>1e-2))*(Uyw+Uys)/((hw+hs)*((hw+hs)>1e-2)))
// (Uxw+Uxs)/((hw+hs)*((hw+hs)>1e-2))
// sqrt((Uxw+Uxs)*(Uxw+Uxs) + (Uyw+Uys)*(Uyw+Uys))


static constexpr char VARNAME_1[255] = "dem"; 
static constexpr char VARNAME_2[255] = "mask_in"; 


static std::vector<double> dem;
static std::vector<double> h_initial_cond;

double h_min, L, H, res, density, density_w, density_s;
int NUM_REFINEMENTS, Nx, Ny, number_FD_points;



// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }


static int
raster_2_vector(const double& i_x,
                const double& i_y)
{
  static int ii;
  ii = i_y + Ny*i_x;
  
  return(ii);
}


size_t
dof_ordering_pressure (p4est_gloidx_t gt, size_t n)
{ return number_FD_points*gt+n; }


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

using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = distributed_vector; //distributed_vector; //std::vector<double>;         // Typedef for local q_0 vector // distributed_vector


double dem_fun (const double& xx, const double& yy)
{ 
  //double aa = (xx>5 && xx<50 ? -xx+100 : xx<=5 ? 95. : 50.);
  //aa += (xx>60 && xx<100 && yy>80 && yy<120 ? 30 : 0.);
  //return( aa ); 
  //return(0);
  //return(5.*std::exp(-2./5*(xx-5.)*(xx-5.)));
  //return(xx>4 && xx<8 ? 4. : 0.);
  //return(-xx+L);
  return(raster_value(xx,yy,dem));
}

double poro_0_fun (const double& xx, const double& yy)
{
  //return(xx>.5 && xx<1.5 ? .4 : .5);
  //return(xx<L/2. ? .3 : .6); 
  //std::cout << (density_s - density)/(density_s - density_w) << std::endl;
  return((density_s - density)/(density_s - density_w));
} 


double h0_fun (const double& xx, const double& yy) 
{ 
  //return(xx<10. ? 10. : 0.);
  //return(xx>4.5 && xx<5.5 ? 1. : .5);
  //return(1.);
  return(std::sqrt( (xx-L/2.)*(xx-L/2.) + (yy-H/2.)*(yy-H/2.) )<=L/10 ? 10 : 0. );
  return(std::abs(xx-L/2.)<=L/10. && std::abs(yy-H/2.)<=H/10. ? 10. : 0.  );
  //return(xx<10. ? 1.e-8 : 0.);
  //return(yy>L/2. ? 10. : 0.);
  //return (xx<=L/2. ? 100. : 50.);
  //return(10.-dem_fun(xx,yy));
  //return(std::abs(xx-L/2.)<=L/10. && std::abs(yy-H/2.)<=H/10. ? 10. : 0.  );
  //return(std::sqrt( (xx-L/2.)*(xx-L/2.) + (yy-H/2.)*(yy-H/2.) )<=L/10 ? 10 : 0. );
  return(raster_value(xx,yy,h_initial_cond));
}

double Ux0_w_fun (const double& xx, const double& yy) 
{ 
  //return(.2*poro_0_fun(xx,yy));
  //return(xx<L/2. ? .3*poro_0_fun(xx,yy)*h0_fun(xx,yy) : -.1*poro_0_fun(xx,yy)*h0_fun(xx,yy));
  //return(-.3);
  return 0.; 
}
double Uy0_w_fun (const double& xx, const double& yy) 
{ 
  return 0.; 
}


double Ux0_s_fun (const double& xx, const double& yy) 
{ 
  //return(-.3*(1.-poro_0_fun(xx,yy)));
  //return(xx<L/2. ? -1.4*(1.-poro_0_fun(xx,yy))*h0_fun(xx,yy) : -.9*(1.-poro_0_fun(xx,yy))*h0_fun(xx,yy));
  //return(.2);
  return 0.; 
}
double Uy0_s_fun (const double& xx, const double& yy) 
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
  const bool   & is_pore_water_pressure                   = input_data["do you want the pore water pressure?"];
  const bool   & is_max_time_step_from_CFL                = input_data["do you want the maximum time step given by CFL condition for the transport term?"];
  const double & h_min_p                                  = input_data["minimum material height threshold for the pressure equation"];               
                 h_min                                    = input_data["minimum material height threshold"];
  const double & grav                                     = input_data["gravitational field"];
                 density                                  = input_data["material density"];
                 density_s                                = input_data["solid density"];
                 density_w                                = input_data["fluid density"];
  const double & turbulence_coeff                         = input_data["turbulence coefficient"];
  const double & surface_pressure                         = input_data["surface atmospheric pressure"];
        double   bed_friction_angle_rad                   = input_data["bed friction angle in degrees"];
                 bed_friction_angle_rad                  *= M_PI/180;
  const double & erosion_coefficient                      = input_data["erosion coefficient"];
  const double & terminal_velocity                        = input_data["terminal velocity"];
        double   odometric_coeff                          = input_data["odometric coefficient"];
        double   consolidation_coefficient                = input_data["consolidation coefficient"];
  const double & m_coeff                                  = input_data["m coefficient"];
  const double & Young_modulus                            = input_data["Young modulus"];
  const double & Poissons_ratio                           = input_data["Poisson's ratio"];
  const double & alfa_coeff                               = input_data["alpha coefficient"];
                 number_FD_points                         = input_data["number of points in the FD mesh"];
  const double & thickness_basal_layer                    = input_data["thickness basal layer"];
  const double & thr_erodible_layer                       = input_data["threshold of the orography height in meters under which we have no erodible layer"];
  const double & tolerance_space_adapt                    = input_data["tolerance space adaptation"];

  const std::string & SAVE_DIR    = input_data["home saving directory, i.e., where we can find the directory results"];
  const std::string & DEM_DIR     = input_data["dem file, complete path"]; 
  const std::string & MASK_DIR    = input_data["mask file, complete path"];

  L = res*(Nx-1);
  H = res*(Ny-1);

  //const auto K_v = Young_modulus/(3*(1.-2.*Poissons_ratio));
  const auto number_FD_elements = number_FD_points-1;


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



  // Management of solutions ordering for mass and momentum
  ordering ordhw  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<6, 0> (gt); };
  ordering ordhs  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<6, 1> (gt); };
  ordering ordUxw = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<6, 2> (gt); };
  ordering ordUyw = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<6, 3> (gt); };
  ordering ordUxs = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<6, 4> (gt); };
  ordering ordUys = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<6, 5> (gt); };


  // Management of solutions ordering for excess pwp
  ordering ordBottom  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, 0); };
  ordering ordSurface = [&number_FD_elements] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, number_FD_elements); };

  
  
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
  TOC ("Uniform refinement");
  
  
  
  // ln_nodes sono i dof non gli hanging node!! (sono esclusi dal calcolo)
  tmesh::idx_t gn_nodes    = tmsh.num_global_nodes (); // Return total number of nodes owned by all process
  tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes (); // Return number of nodes owned by local process
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();  // Return number of quadrants owned by local process across all trees
  tmesh::idx_t gn_elements = tmsh.num_global_quadrants (); // Return number of quadrants owned by all processes across all trees
  

  /// Allocate initial data container
  Q1 sol  (ln_nodes * 6);
  sol.get_owned_data  ().assign (sol.get_owned_data  ().size (), 0.0);


  Q1 incr (ln_nodes * 6);
  incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
  incr.assemble();

  
  Q1 mass (ln_nodes * 6);
  bim2a_mass_vector (tmsh, mass, ordhw);
  bim2a_mass_vector (tmsh, mass, ordhs);
  bim2a_mass_vector (tmsh, mass, ordUxw);
  bim2a_mass_vector (tmsh, mass, ordUyw);
  bim2a_mass_vector (tmsh, mass, ordUxs);
  bim2a_mass_vector (tmsh, mass, ordUys);
  mass.assemble ();
  
  Q0 sol_onehalf (ln_elements * 6);
  sol_onehalf.get_owned_data ().assign (sol_onehalf.get_owned_data ().size (), 0.0);
  sol_onehalf.assemble();

  Q0 Z_onehalf (ln_elements);
  Z_onehalf.get_owned_data ().assign (Z_onehalf.get_owned_data ().size (), 0.0);
  Z_onehalf.assemble();

  std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 6);

  std::vector<std::array<double,4>> incr_anti_diff_pressure (ln_elements * number_FD_points);
  
  Q1 Z (ln_nodes);
  Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);

  Q1 excess_pore_water_pressure(ln_nodes * number_FD_points);
  excess_pore_water_pressure.get_owned_data  ().assign (excess_pore_water_pressure.get_owned_data  ().size (), 0.0);

  Q0 excess_pore_water_pressure_onehalf(ln_elements * number_FD_points);
  excess_pore_water_pressure_onehalf.get_owned_data  ().assign (excess_pore_water_pressure_onehalf.get_owned_data  ().size (), 0.0);

  std::string str = ""; 
  char filename[255]="", arr[255]="";

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
  h_initial_cond.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), h_initial_cond.begin ());
  TOC("Load data matrix");
  

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

        
        const double initial_porosity_coeff = poro_0_fun(xx,yy); 

        sol [ordhw    (quadrant->gt (ii))] = h0_fun    (xx, yy)*initial_porosity_coeff;
        sol [ordhs    (quadrant->gt (ii))] = h0_fun    (xx, yy)*(1.-initial_porosity_coeff);
        sol [ordUxw   (quadrant->gt (ii))] = Ux0_w_fun (xx, yy);
        sol [ordUyw   (quadrant->gt (ii))] = Uy0_w_fun (xx, yy);
        sol [ordUxs   (quadrant->gt (ii))] = Ux0_s_fun (xx, yy);
        sol [ordUys   (quadrant->gt (ii))] = Uy0_s_fun (xx, yy);
        
        Z   [          quadrant->gt (ii) ] = dem_fun (xx, yy);

        for (int kk=ordBottom(quadrant->gt (ii)); kk<=ordSurface(quadrant->gt (ii)); kk++)
        {
          const auto current_z_coord = h0_fun(xx, yy)/number_FD_elements*(kk%number_FD_points);
          excess_pore_water_pressure[kk] = is_pore_water_pressure ? (1.-initial_porosity_coeff)*(density_s - density_w)*grav*(h0_fun(xx, yy)-current_z_coord) : 0.; //current_z_coord<=thickness_basal_layer ? (1.-initial_porosity_coeff)*(density_s - density_w)*grav*h0_fun(xx, yy)*1. /*0.65*/ : 0.;
        }
      }
      
      else
      {
        // touch parent nodes to set up distributed vector structure
        sol [ordhw   (quadrant->gparent(0,ii))] += 0.;
        sol [ordhw   (quadrant->gparent(1,ii))] += 0.;
        sol [ordhs   (quadrant->gparent(0,ii))] += 0.;
        sol [ordhs   (quadrant->gparent(1,ii))] += 0.;
        sol [ordUxw  (quadrant->gparent(0,ii))] += 0.;
        sol [ordUxw  (quadrant->gparent(1,ii))] += 0.;
        sol [ordUyw  (quadrant->gparent(0,ii))] += 0.;
        sol [ordUyw  (quadrant->gparent(1,ii))] += 0.;
        sol [ordUxs  (quadrant->gparent(0,ii))] += 0.;
        sol [ordUxs  (quadrant->gparent(1,ii))] += 0.;
        sol [ordUys  (quadrant->gparent(0,ii))] += 0.;
        sol [ordUys  (quadrant->gparent(1,ii))] += 0.;
        
        Z   [         quadrant->gparent(0,ii) ] += 0.;
        Z   [         quadrant->gparent(1,ii) ] += 0.;

        for (int jj=0; jj<=1; jj++)
          for (int kk=ordBottom(quadrant->gparent(jj,ii)); kk<=ordSurface(quadrant->gparent(jj,ii)); kk++)
          {
            excess_pore_water_pressure[kk] += 0.;
          }
      }
    }
  }

  // bim2a_solution_with_ghosts in quad_operators.cpp
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordhw,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordhs,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUxw, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUyw, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUxs, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUys);

  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordhw,  false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordhs,  false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUxw, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUyw, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUxs, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUys);
  
  bim2a_solution_with_ghosts (tmsh, Z, replace_op);

  bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure, replace_op, ordBottom, false);
  for (int kk=1; kk<number_FD_elements; kk++)
  {
    ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
    bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure, replace_op, ordInternal, false);
  } 
  bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure, replace_op, ordSurface);



  if (is_initial_refinement)
  {
    
    TIC();
    Q1 only_h (ln_nodes);
    bim2a_solution_with_ghosts (tmsh, only_h);
    for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
    {
      only_h(idx) = sol(ordhw(idx))+sol(ordhs(idx));
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


    tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, 1e-5, 4, 0, 0);
    //tmsh.set_metrics_marker (estimator, 1e-5, 4, 3, 1);
    tmsh.metrics_refine (1e4);  // RAFFINAMENTO (arg is max element)

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
  
    Q1 incr_ (ln_nodes * 6);
    incr_.get_owned_data ().assign (incr_.get_owned_data ().size(), 0.0);
    incr_.assemble ();
  
    Q1 mass_ (ln_nodes * 6);
    bim2a_mass_vector (tmsh, mass_, ordhw );
    bim2a_mass_vector (tmsh, mass_, ordhs );
    bim2a_mass_vector (tmsh, mass_, ordUxw);
    bim2a_mass_vector (tmsh, mass_, ordUyw);
    bim2a_mass_vector (tmsh, mass_, ordUxs);
    bim2a_mass_vector (tmsh, mass_, ordUys);
    mass_.assemble ();
  
    Q0 sol_onehalf_ (ln_elements * 6);
    sol_onehalf_.get_owned_data ().assign (sol_onehalf_.get_owned_data ().size(), 0.0);
    sol_onehalf_.assemble();

    Q0 Z_onehalf_ (ln_elements);
    Z_onehalf_.get_owned_data ().assign (Z_onehalf_.get_owned_data ().size(), 0.0);
    Z_onehalf_.assemble();

    std::vector<std::array<double,4>> incr_anti_diff_ (ln_elements * 6);

    std::vector<std::array<double,4>> incr_anti_diff_pressure_ (ln_elements * number_FD_points);

    Q1 excess_pore_water_pressure_(ln_nodes * number_FD_points);
    excess_pore_water_pressure_.get_owned_data  ().assign (excess_pore_water_pressure_.get_owned_data  ().size (), 0.0);

    Q0 excess_pore_water_pressure_onehalf_(ln_elements * number_FD_points);
    excess_pore_water_pressure_onehalf_.get_owned_data  ().assign (excess_pore_water_pressure_onehalf_.get_owned_data  ().size (), 0.0);


    Q1 sol_ (ln_nodes * 6);
    Q1 Z_ (ln_nodes);
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


          const double initial_porosity_coeff = poro_0_fun (xx,yy); //(density_s - density)/(density_s - density_w);
          
          sol_ [ordhw    (quadrant->gt (ii))] = h0_fun  (xx, yy)*initial_porosity_coeff;
          sol_ [ordhs    (quadrant->gt (ii))] = h0_fun  (xx, yy)*(1.-initial_porosity_coeff);
          sol_ [ordUxw   (quadrant->gt (ii))] = Ux0_w_fun (xx, yy);
          sol_ [ordUyw   (quadrant->gt (ii))] = Uy0_w_fun (xx, yy);
          sol_ [ordUxs   (quadrant->gt (ii))] = Ux0_s_fun (xx, yy);
          sol_ [ordUys   (quadrant->gt (ii))] = Uy0_s_fun (xx, yy);
          

          Z_[quadrant->gt (ii)] = dem_fun (xx, yy); 


          for (int kk=ordBottom(quadrant->gt (ii)); kk<=ordSurface(quadrant->gt (ii)); kk++)
          {
            const auto current_z_coord = h0_fun(xx, yy)/number_FD_elements*(kk%number_FD_points);
            excess_pore_water_pressure_[kk] = is_pore_water_pressure ? (1.-initial_porosity_coeff)*(density_s - density_w)*grav*(h0_fun(xx, yy)-current_z_coord) : 0.; //is_pore_water_pressure ? (current_z_coord<=thickness_basal_layer ? (1.-initial_porosity_coeff)*(density_s - density_w)*grav*h0_fun(xx, yy)*0. /*0.65*/ : 0.) : 0.;
          }
        }
        
        else
        {
          sol_ [ordhw   (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordhw   (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordhs   (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordhs   (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordUxw  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUxw  (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordUyw  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUyw  (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordUxs  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUxs  (quadrant->gparent(1,ii))] += 0.;
          sol_ [ordUys  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUys  (quadrant->gparent(1,ii))] += 0.;

          Z_[quadrant->gparent(0,ii)] += 0.;
          Z_[quadrant->gparent(1,ii)] += 0.;


          for (int jj=0; jj<=1; jj++)
            for (int kk=ordBottom(quadrant->gparent(jj,ii)); kk<=ordSurface(quadrant->gparent(jj,ii)); kk++)
            {
              excess_pore_water_pressure_[kk] += 0.;
            }
        }
      }
    }

    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordhw,  false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordhs,  false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUxw, false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUyw, false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUxs, false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUys);

    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordhw,  false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordhs,  false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUxw, false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUyw, false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUxs, false);
    bim2a_solution_with_ghosts (tmsh, incr_, replace_op, ordUys);

    bim2a_solution_with_ghosts (tmsh, Z_, replace_op);

    bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_, replace_op, ordBottom, false);
    for (int kk=1; kk<number_FD_elements; kk++)
    {
      ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
      bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_, replace_op, ordInternal, false);
    } 
    bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_, replace_op, ordSurface);


    sol                 = sol_;
    incr                = incr_;
    incr_anti_diff      = incr_anti_diff_;
    mass                = mass_;
    sol_onehalf         = sol_onehalf_;
    Z                   = Z_;
    Z_onehalf           = Z_onehalf_;

    excess_pore_water_pressure         = excess_pore_water_pressure_;
    excess_pore_water_pressure_onehalf = excess_pore_water_pressure_onehalf_;

    incr_anti_diff_pressure = incr_anti_diff_pressure_;
  
    TOC ("compute initial condition");
  }
  
  Q1 sol_dyn                 = sol;
  Q1 sold_dyn                = sol;
  Q1 soldd_dyn               = sol;
  Q1 incr_dyn                = incr;
  Q1 P_plus_dyn              = incr;
  Q1 P_minus_dyn             = incr;
  Q1 mass_dyn                = mass;
  Q1 Z_dyn                   = Z;
  Q0 sol_onehalf_dyn         = sol_onehalf;
  Q0 Z_onehalf_dyn           = Z_onehalf;

  std::vector<std::array<double,4>> incr_anti_diff_dyn = incr_anti_diff;
  std::vector<std::array<double,4>> incr_anti_diff_pressure_dyn = incr_anti_diff_pressure;

  Q1 P_plus_pressure_dyn              = excess_pore_water_pressure;
  Q1 P_minus_pressure_dyn             = excess_pore_water_pressure;

  Q1 excess_pore_water_pressure_dyn = excess_pore_water_pressure;
  Q0 excess_pore_water_pressure_onehalf_dyn = excess_pore_water_pressure_onehalf;
  Q1 excess_pore_water_pressure_incr_dyn    = excess_pore_water_pressure;

  
  TG2_scheme stp(sol_dyn, 
                 sold_dyn, 
                 soldd_dyn, 
                 incr_dyn,
                 incr_anti_diff_dyn,
                 incr_anti_diff_pressure_dyn,
                 P_plus_dyn, 
                 P_minus_dyn, 
                 P_plus_pressure_dyn, 
                 P_minus_pressure_dyn, 
                 sol_onehalf_dyn, 
                 mass_dyn,
                 excess_pore_water_pressure_dyn,
                 excess_pore_water_pressure_onehalf_dyn,
                 excess_pore_water_pressure_incr_dyn,
                 ordhw, 
                 ordhs, 
                 ordUxw, 
                 ordUyw, 
                 ordUxs, 
                 ordUys, 
                 ordBottom,
                 ordSurface,
                 Z_dyn,
                 Z_onehalf_dyn,
                 DELTAT, 
                 h_min, 
                 is_non_reflBC, 
                 is_bed_friction, 
                 is_pore_water_pressure,
                 grav,
                 density_w, 
                 density_s, 
                 turbulence_coeff, 
                 bed_friction_angle_rad, 
                 erosion_coefficient, 
                 m_coeff, 
                 terminal_velocity,
                 odometric_coeff,
                 consolidation_coefficient,
                 thr_erodible_layer,
                 number_FD_points);

  stp.set_r_coeff();
  stp.resize_vectors();


  // Save initial conditions
  str = std::string(SAVE_DIR) + "/results/swe_hw_%4.4d"; 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordhw);

  str = std::string(SAVE_DIR) + "/results/swe_hs_%4.4d"; 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordhs);

  str = std::string(SAVE_DIR) + "/results/swe_Uxw_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordUxw); 

  str = std::string(SAVE_DIR) + "/results/swe_Uyw_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordUyw);

  str = std::string(SAVE_DIR) + "/results/swe_Uxs_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordUxs); 

  str = std::string(SAVE_DIR) + "/results/swe_Uys_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol_dyn, ordUys);
  
  str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, Z_dyn);

  for (int kk=0; kk<number_FD_points; kk++)
  {
    str = std::string(SAVE_DIR) + "/results/swe_p_" + std::to_string(kk) + "_%4.4d";
    strcpy(arr, str.c_str());
    sprintf(filename, arr,  0);
    ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
    tmsh.octbin_export (filename, excess_pore_water_pressure_dyn, ordInternal);
  } 

  std::vector<double> full_time_vector;
  full_time_vector.reserve (static_cast<int> (T/DELTAT));
  std::vector<double> save_time_vector;
  save_time_vector.reserve (static_cast<int> (T/SAVEDT));

  
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

  int restart_loop = 0;
  double max_dt_post = max_dt;

  TIC();
  while (counter_savings != tot_number_savings)
  {
    // Reset increment, and limiter terms
    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);

    P_plus_dyn.get_owned_data ().assign (P_plus_dyn.get_owned_data ().size (), 0.0);
    P_plus_dyn.assemble (replace_op);

    P_minus_dyn.get_owned_data ().assign (P_minus_dyn.get_owned_data ().size (), 0.0);
    P_minus_dyn.assemble (replace_op);

    P_plus_pressure_dyn.get_owned_data ().assign (P_plus_pressure_dyn.get_owned_data ().size (), 0.0);
    P_plus_pressure_dyn.assemble (replace_op);

    P_minus_pressure_dyn.get_owned_data ().assign (P_minus_pressure_dyn.get_owned_data ().size (), 0.0);
    P_minus_pressure_dyn.assemble (replace_op);

    excess_pore_water_pressure_incr_dyn.get_owned_data ().assign (excess_pore_water_pressure_incr_dyn.get_owned_data ().size (), 0.0);
    excess_pore_water_pressure_incr_dyn.assemble (replace_op);


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

      const double local_estimator_time_tolerance = 1e-5;

      const double candidate_dt = local_estimator_time_tolerance/std::sqrt(stp.nu_htot)*(stp.time-stp.timed);
      stp.set_dt( (stp.nu_htot>0 && candidate_dt<stp.dt) ? candidate_dt : stp.dt );
    }

    

    if (stp.dt == 0 && rank == 0)
    {
      std::cout << "dt has gone to zero, sorry, STOP!" << std::endl;
      exit( -1. );
    }

    // check save with given frequency
    stp.set_dt((savecount+stp.dt)/SAVEDT>1 ? SAVEDT-savecount : stp.dt);
    stp.set_tau();



    time_oldd = time_old;
    time_old = time;
    time += stp.dt; 
    savecount += stp.dt;
    space_adapt_count += stp.dt;

    
    // Print current time
    if(rank==0) 
    {
      std::cout << "TIME = " << time << ", dt = " << stp.dt << std::endl;
      full_time_vector.push_back (time);
    }

    
    TIC();
    // first step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.first_step_consolidation(quadrant);
      stp.first_step(quadrant);

      //restart_loop = stp.error>stp.tolerance ? 1 : restart_loop;
    }
    //MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&restart_loop), 1, MPI_INT, MPI_SUM, tmsh.comm);
/*
    if (restart_loop>0)
    {
      max_dt_post = max_dt*.5;
      restart_loop = 0;
      continue;
    }
    else
    {
      max_dt_post = DELTAT*REDCDT;
    }*/

    bim2a_solution_with_ghosts_center (tmsh, sol_onehalf_dyn, replace_op, ordhw,  false);
    bim2a_solution_with_ghosts_center (tmsh, sol_onehalf_dyn, replace_op, ordhs,  false);
    bim2a_solution_with_ghosts_center (tmsh, sol_onehalf_dyn, replace_op, ordUxw, false);
    bim2a_solution_with_ghosts_center (tmsh, sol_onehalf_dyn, replace_op, ordUyw, false);
    bim2a_solution_with_ghosts_center (tmsh, sol_onehalf_dyn, replace_op, ordUxs, false);
    bim2a_solution_with_ghosts_center (tmsh, sol_onehalf_dyn, replace_op, ordUys);

    bim2a_solution_with_ghosts_center (tmsh, Z_onehalf_dyn, replace_op);


    bim2a_solution_with_ghosts_center (tmsh, excess_pore_water_pressure_onehalf_dyn, replace_op, ordBottom, false);
    for (int kk=1; kk<number_FD_elements; kk++)
    {
      ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
      bim2a_solution_with_ghosts_center (tmsh, excess_pore_water_pressure_onehalf_dyn, replace_op, ordInternal, false);
    } 
    bim2a_solution_with_ghosts_center (tmsh, excess_pore_water_pressure_onehalf_dyn, replace_op, ordSurface);
    TOC("first step");


    //
    TIC();
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
    

    // low order solution, 
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
    {
      sol_dyn.get_owned_data ()[kk] += stp.dt*incr_dyn.get_owned_data ()[kk]/mass_dyn.get_owned_data ()[kk];
    }
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordhw,  false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordhs,  false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUxw, false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUyw, false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUxs, false);
    bim2a_solution_with_ghosts (tmsh, sol_dyn, replace_op, ordUys);


    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);
    TOC("compute_nodal_anti_diffusive_fluxes");


    // second order correction
    TIC();
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.second_step(quadrant);
    }
    incr_dyn.assemble ();
    TOC("Compute step");


    // get the updated solution
    TIC();
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=6)
    {
      stp.solve_non_lin_h(kk);
      stp.solve_non_lin_U(kk);
    }


    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii) && sol_dyn [ordhw    (quadrant->gt (ii))]<0){
          sol_dyn [ordhw    (quadrant->gt (ii))] = 0.; //h_min; //0.;
        }
        if (! quadrant->is_hanging (ii) && sol_dyn [ordhs    (quadrant->gt (ii))]<0){
          sol_dyn [ordhs    (quadrant->gt (ii))] = 0.; //h_min; //0.;
        }
      }
    }
    sol_dyn.assemble(replace_op);

    // solve the consolidation problem here,
    incr_dyn.get_owned_data ().assign (incr_dyn.get_owned_data ().size (), 0.0);
    incr_dyn.assemble (replace_op);


    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.terminate_second_step(quadrant);
    }
    TOC("Newton");


    TIC();
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.solve_second_step_cons_equation(quadrant); 
    }
    excess_pore_water_pressure_incr_dyn.assemble();
    P_plus_pressure_dyn.assemble (); 
    P_minus_pressure_dyn.assemble ();
    TOC("excess_pore_water_pressure_incr_dyn");

    

    TIC();
    for (auto kk = 0; kk < excess_pore_water_pressure_incr_dyn.get_owned_data ().size (); kk+=number_FD_points)
    {
      const int kk_single = kk/number_FD_points;
      const int kk_ = kk_single*6; 

      const double hdofold = sold_dyn.get_owned_data ()[kk_] + sold_dyn.get_owned_data ()[kk_+1];
      const double hdof_c  = sol_dyn .get_owned_data ()[kk_] + sol_dyn .get_owned_data ()[kk_+1];

      const double delta_h_old = hdofold/(number_FD_points-1);

      //consolidation_coefficient = hdofold>h_min ? odometric_coeff/((sold_dyn.get_owned_data ()[kk_]>h_min && hdofold>h_min && sold_dyn.get_owned_data ()[kk_+1]>h_min) ? sold_dyn.get_owned_data ()[kk_+1]/hdofold/std::pow(sold_dyn.get_owned_data ()[kk_]/hdofold, m_coeff)/terminal_velocity*(density_s-density_w)*grav : 0.) : 1.e-10;

      if (stp.dt<=.5*delta_h_old*delta_h_old/consolidation_coefficient*stp.cfl_dp || consolidation_coefficient==0) // explicit case
      {
        //std::cout << "explicit" << std::endl;
        // In case, Neumann BC 
        //std::cout << excess_pore_water_pressure_incr_dyn.get_owned_data ()[kk] << std::endl;
        excess_pore_water_pressure_dyn.get_owned_data ()[kk] = Z_dyn.get_owned_data ()[kk_single]<thr_erodible_layer ? stp.dt*excess_pore_water_pressure_incr_dyn.get_owned_data ()[kk]/mass_dyn.get_owned_data ()[kk_] : excess_pore_water_pressure_dyn.get_owned_data ()[kk];
        for (int kkk=kk+1; kkk<kk+number_FD_points-1; kkk++) // eliminate the boundaries 
        { 
          excess_pore_water_pressure_dyn.get_owned_data ()[kkk] = hdof_c>h_min ? stp.dt*excess_pore_water_pressure_incr_dyn.get_owned_data ()[kkk]/mass_dyn.get_owned_data ()[kk_] : 0.;
        }
      } 
      else // implicit case   
      {
        //std::cout << "implicit" << std::endl;
        //return 0;
        const double delta_h = hdof_c/(number_FD_points-1);

        //consolidation_coefficient = hdof_c>h_min ? odometric_coeff/((sol_dyn.get_owned_data ()[kk_]>h_min && hdof_c>h_min && sol_dyn.get_owned_data ()[kk_+1]>h_min) ? sol_dyn.get_owned_data ()[kk_+1]/hdof_c/std::pow(sol_dyn.get_owned_data ()[kk_]/hdof_c, m_coeff)/terminal_velocity*(density_s-density_w)*grav : 0.) : 1.e-10;

        const double mu_coeff = delta_h>h_min ? consolidation_coefficient*stp.dt/delta_h/delta_h : 0.;

        // ND and DD depending on the thr_erodible_layer,
        Z_dyn.get_owned_data ()[kk_single]<thr_erodible_layer ? stp.thomas_algorithm(mu_coeff, kk_, kk, kk+number_FD_points-2, kk) : stp.thomas_algorithm(mu_coeff, kk_, kk+1, kk+number_FD_points-2, kk);
      }
      
    }
    bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_dyn, replace_op, ordBottom, false);
    for (int kk=1; kk<number_FD_elements; kk++)
    {
      ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
      bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_dyn, replace_op, ordInternal, false);
    } 
    bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_dyn, replace_op, ordSurface);
    TOC ("excess_pore_water_pressure_dyn.assemble(replace_op)");



    TIC();
    excess_pore_water_pressure_incr_dyn.get_owned_data ().assign (excess_pore_water_pressure_incr_dyn.get_owned_data ().size (), 0.0);
    excess_pore_water_pressure_incr_dyn.assemble (replace_op);

    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      for (int kkk=1; kkk<number_FD_points-1; kkk++) // eliminate the boundaries 
      {
        stp.second_step_pressure(quadrant, kkk);
      }
    }
    excess_pore_water_pressure_incr_dyn.assemble ();

    for (auto kk = 0; kk < excess_pore_water_pressure_incr_dyn.get_owned_data ().size (); kk++)
    {
      const int kk_ = kk/number_FD_points;
      excess_pore_water_pressure_dyn.get_owned_data ()[kk] += stp.dt*excess_pore_water_pressure_incr_dyn.get_owned_data ()[kk]/mass_dyn.get_owned_data ()[kk_*6];
    }

    for (auto kk = 0; kk < excess_pore_water_pressure_incr_dyn.get_owned_data ().size (); kk+=number_FD_points)
    {
      double dp_mean = 0.;
      stp.numerical_integration_pressure(kk, dp_mean, excess_pore_water_pressure_dyn, 1.);

      const int kk_ = kk/number_FD_points;
      const double hw_current_node = sol_dyn.get_owned_data ()[kk_*6];
      const double hs_current_node = sol_dyn.get_owned_data ()[kk_*6+1];

      const double h_current_node = hw_current_node + hs_current_node;
      //const double gamma_coeff = std::min(dp_mean + density_w*grav*h_current_node, 0.)*stp.sf + std::max(dp_mean - (hw_current_node>h_min ? density_w*grav*h_current_node*(1.+stp.r_coeff)*.5*hs_current_node/hw_current_node : 0.), 0.);//*(2. - stp.sf);

      for (int kkk=kk; kkk<kk+number_FD_points; kkk++) 
      { 
        //const double zeta_greek_current = (kkk%number_FD_points)/double(number_FD_elements);
        //const double func_distr = 6.*zeta_greek_current*(1. - zeta_greek_current);
        //excess_pore_water_pressure_dyn.get_owned_data ()[kkk] -= gamma_coeff*func_distr;
        auto & dp_p_c = excess_pore_water_pressure_dyn.get_owned_data ()[kkk];
        dp_p_c -= std::min(dp_p_c + density_w*grav*h_current_node, 0.)*stp.sf + std::max(dp_p_c - (hw_current_node>h_min ? density_w*grav*h_current_node*(1.+stp.r_coeff)*.5*hs_current_node/hw_current_node : 0.), 0.);
        dp_p_c = h_current_node>h_min_p ? dp_p_c : 0.;
      }

      //dp_mean = 0;
      //stp.numerical_integration_pressure(kk, dp_mean, excess_pore_water_pressure_dyn, 1.);
      //if (gamma_coeff<0 && std::min(dp_mean + density_w*grav*h_current_node, 0.)<0)
      //{
      //  std::cout << gamma_coeff << " " << std::min(dp_mean + density_w*grav*h_current_node, 0.) << std::endl;
      //}

    }
    excess_pore_water_pressure_dyn.assemble(replace_op);
    TOC("limiting step");
    
  
    TIC();
    // now that you have the dp, update the mass flux
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.terminate_second_step(quadrant);
    }
    incr_dyn.assemble ();
    TOC ("terminate_second_step");

    /*
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk++)
    {
      if (std::isnan(sol_dyn.get_owned_data ()[kk]))
      {
        std::cout << sol_dyn.get_owned_data ()[kk] << " " << "stoppppp" << std::endl;
        exit(1);
      }
    }
    */


    TIC();
    for (auto kk = 0; kk < incr_dyn.get_owned_data ().size (); kk+=6)
    {
      sol_dyn.get_owned_data ()[kk+2] += stp.dt*incr_dyn.get_owned_data ()[kk+2]/mass_dyn.get_owned_data ()[kk+2];
      sol_dyn.get_owned_data ()[kk+3] += stp.dt*incr_dyn.get_owned_data ()[kk+3]/mass_dyn.get_owned_data ()[kk+3];
      sol_dyn.get_owned_data ()[kk+4] += stp.dt*incr_dyn.get_owned_data ()[kk+4]/mass_dyn.get_owned_data ()[kk+4];
      sol_dyn.get_owned_data ()[kk+5] += stp.dt*incr_dyn.get_owned_data ()[kk+5]/mass_dyn.get_owned_data ()[kk+5];

      stp.stabilization_term(kk);
    }
    sol_dyn.assemble(replace_op);
    TOC("stabilization_term");



    // Save solution
    if ((savecount-SAVEDT) >= -std::numeric_limits<double>::epsilon()*SAVEDT) 
    {
      //TIC();
      if (rank == 0)
        std::cout << "savecount = " << savecount << std::endl;
      count++;
      save_time_vector.push_back (time); 

     
      str = std::string(SAVE_DIR) + "/results/swe_hw_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,   count);
      tmsh.octbin_export (filename, sol_dyn, ordhw);

      str = std::string(SAVE_DIR) + "/results/swe_hs_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,   count);
      tmsh.octbin_export (filename, sol_dyn, ordhs);
      
      str = std::string(SAVE_DIR) + "/results/swe_Uxw_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_dyn, ordUxw);
      
      str = std::string(SAVE_DIR) + "/results/swe_Uyw_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_dyn, ordUyw);

      str = std::string(SAVE_DIR) + "/results/swe_Uxs_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_dyn, ordUxs);
      
      str = std::string(SAVE_DIR) + "/results/swe_Uys_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_dyn, ordUys);
      
      str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, Z_dyn);

      for (int kk=0; kk<number_FD_points; kk++)
      {
        str = std::string(SAVE_DIR) + "/results/swe_p_" + std::to_string(kk) + "_%4.4d";
        strcpy(arr, str.c_str());
        sprintf(filename, arr,  count);
        ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
        tmsh.octbin_export (filename, excess_pore_water_pressure_dyn, ordInternal);
      } 

      savecount = 0.0;
      //TOC("Exporting solution");

      counter_savings++;

    }


    // sistemare adattazione spaziale!!
    if (is_space_adaptivity && ((space_adapt_count-SPACE_ADAPTDT) >= -std::numeric_limits<double>::epsilon()*SPACE_ADAPTDT))
    //(is_space_adaptivity && counter_savings%8==0)//  ((space_adapt_count-SPACE_ADAPTDT) >= -std::numeric_limits<double>::epsilon()*SPACE_ADAPTDT))
    {

      //TIC();
      Q1 only_h (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_h);
      for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
      {
        only_h(idx) = sol_dyn(ordhw(idx))+sol_dyn(ordhs(idx));
      }
      only_h.assemble (replace_op);

      Q1 only_Ux (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_Ux);
      for (auto idx = only_Ux.get_range_start (); idx != only_Ux.get_range_end (); ++idx)
      {
        only_Ux(idx) = sol_dyn(ordUxw(idx)) + sol_dyn(ordUxs(idx));
      }
      only_Ux.assemble (replace_op);

      Q1 only_Uy (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_Uy);
      for (auto idx = only_Uy.get_range_start (); idx != only_Uy.get_range_end (); ++idx)
      {
        only_Uy(idx) = sol_dyn(ordUyw(idx)) + sol_dyn(ordUys(idx));
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
      auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
      {
        return estimator_grad(q, dh, only_h);
      };


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


      tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, 1e-5, 4, 0, 0);
      tmsh.metrics_refine (1e6);  // RAFFINAMENTO (arg is max element)


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
      Q1 sol (ln_nodes * 6);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordhw,  false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordhs,  false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUxw, false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUyw, false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUxs, false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUys);
      interpolate_vector (tmsh, sol_dyn, sol, ordhw );
      interpolate_vector (tmsh, sol_dyn, sol, ordhs );
      interpolate_vector (tmsh, sol_dyn, sol, ordUxw);
      interpolate_vector (tmsh, sol_dyn, sol, ordUyw);
      interpolate_vector (tmsh, sol_dyn, sol, ordUxs);
      interpolate_vector (tmsh, sol_dyn, sol, ordUys);
      sol.assemble (replace_op);
      
      Q1 sold (ln_nodes * 6);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordhw,  false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordhs,  false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUxw, false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUyw, false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUxs, false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUys);
      interpolate_vector (tmsh, sold_dyn, sold, ordhw );
      interpolate_vector (tmsh, sold_dyn, sold, ordhs );
      interpolate_vector (tmsh, sold_dyn, sold, ordUxw);
      interpolate_vector (tmsh, sold_dyn, sold, ordUyw);
      interpolate_vector (tmsh, sold_dyn, sold, ordUxs);
      interpolate_vector (tmsh, sold_dyn, sold, ordUys);
      sold.assemble (replace_op);
      
      
      Q1 soldd (ln_nodes * 6);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordhw,  false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordhs,  false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUxw, false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUyw, false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUxs, false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUys);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordhw );
      interpolate_vector (tmsh, soldd_dyn, soldd, ordhs );
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUxw);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUyw);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUxs);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUys);
      soldd.assemble (replace_op);


      Q1 excess_pore_water_pressure (ln_nodes * number_FD_points);
      for (int kk=0; kk<number_FD_elements; kk++)
      {
        ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
        bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure, replace_op, ordInternal,  false);
      } 
      bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure, replace_op, ordSurface);
      for (int kk=0; kk<number_FD_points; kk++)
      {
        ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
        interpolate_vector (tmsh, excess_pore_water_pressure_dyn, excess_pore_water_pressure, ordInternal );
      } 
      excess_pore_water_pressure.assemble(replace_op);


      Q1 excess_pore_water_pressure_incr (ln_nodes * number_FD_points);
      excess_pore_water_pressure_incr.get_owned_data ().assign(excess_pore_water_pressure_incr.get_owned_data ().size(), 0.0);
      excess_pore_water_pressure_incr.assemble();
      for (int kk=0; kk<number_FD_elements; kk++)
      {
        ordering ordInternal = [&kk] (tmesh::idx_t gt) -> size_t { return dof_ordering_pressure(gt, kk); };
        bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_incr, replace_op, ordInternal,  false);
      } 
      bim2a_solution_with_ghosts (tmsh, excess_pore_water_pressure_incr, replace_op, ordSurface);


      Q0 excess_pore_water_pressure_onehalf (ln_elements * number_FD_points);
      excess_pore_water_pressure_onehalf.get_owned_data ().assign (excess_pore_water_pressure_onehalf.get_owned_data ().size(), 0.0);
      excess_pore_water_pressure_onehalf.assemble ();
      
      
      Q1 incr (ln_nodes * 6);
      incr.get_owned_data ().assign (incr.get_owned_data ().size(), 0.0);
      incr.assemble ();

      std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 6);

      std::vector<std::array<double,4>> incr_anti_diff_pressure (ln_elements * number_FD_points);

      
      Q1 mass (ln_nodes * 6);
      bim2a_mass_vector (tmsh, mass, ordhw );
      bim2a_mass_vector (tmsh, mass, ordhs );
      bim2a_mass_vector (tmsh, mass, ordUxw);
      bim2a_mass_vector (tmsh, mass, ordUyw);
      bim2a_mass_vector (tmsh, mass, ordUxs);
      bim2a_mass_vector (tmsh, mass, ordUys);
      mass.assemble ();
      
      Q0 sol_onehalf (ln_elements * 6);
      sol_onehalf.get_owned_data ().assign (sol_onehalf.get_owned_data ().size(), 0.0);
      sol_onehalf.assemble();


      Q0 Z_onehalf (ln_elements);
      Z_onehalf.get_owned_data ().assign (Z_onehalf.get_owned_data ().size(), 0.0);
      Z_onehalf.assemble();


      Q1 Z (ln_nodes);
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
            Z[quadrant->gt (ii)] = dem_fun (xx, yy); 
          }
           
          else
          {
            Z[quadrant->gparent(0,ii)] += 0.;
            Z[quadrant->gparent(1,ii)] += 0.;
          }
        }
      }
      bim2a_solution_with_ghosts (tmsh, Z, replace_op);
      
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordhw,  false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordhs,  false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUxw, false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUyw, false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUxs, false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUys);

      
       
      sol_dyn                 = sol;
      sold_dyn                = sold;
      soldd_dyn               = soldd;
      incr_dyn                = incr;
      incr_anti_diff_dyn      = incr_anti_diff;
      P_plus_dyn              = incr;
      P_minus_dyn             = incr;
      mass_dyn                = mass;
      sol_onehalf_dyn         = sol_onehalf;
      Z_onehalf_dyn           = Z_onehalf;
      Z_dyn                   = Z;


      incr_anti_diff_pressure_dyn      = incr_anti_diff_pressure;

      excess_pore_water_pressure_dyn = excess_pore_water_pressure;
      excess_pore_water_pressure_incr_dyn = excess_pore_water_pressure_incr;
      excess_pore_water_pressure_onehalf_dyn = excess_pore_water_pressure_onehalf;

      P_plus_pressure_dyn  = excess_pore_water_pressure_incr;
      P_minus_pressure_dyn = excess_pore_water_pressure_incr;

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
    octave_io_close ();
  }
  
  TOC ("loop completed");
  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  MPI_Finalize ();
  return 0;
  
}







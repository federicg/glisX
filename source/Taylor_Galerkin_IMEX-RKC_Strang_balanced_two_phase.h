#ifndef TAYLOR_GALERKIN_H
#define TAYLOR_GALERKIN_H

#include <numeric> 
#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>

 


class TG2_scheme   
{
  using Q1  = q1_vec<distributed_vector>;
  using Q0  = distributed_vector;
  
public:
  
  TG2_scheme(Q1& sol,
             Q1& sold,
             Q1& soldd,
             Q1& incr,
             std::vector<std::array<double,4>>& incr_anti_diff,
             std::vector<std::array<double,4>>& incr_anti_diff_pressure,
             Q1& P_plus,
             Q1& P_minus,
             Q1& P_plus_pressure,
             Q1& P_minus_pressure,
             Q0& sol_onehalf,
             Q1& mass,
             Q1& excess_pore_water_pressure,
             Q0& excess_pore_water_pressure_onehalf,
             Q1& excess_pore_water_pressure_incr,
             const ordering& ohw,
             const ordering& ohs,
             const ordering& oUxw,
             const ordering& oUyw,
             const ordering& oUxs,
             const ordering& oUys,
             const ordering& oBottom,
             const ordering& oSurface,
             Q1& Z,
             Q0& Z_onehalf,
             const double& DELTAT,
             const double& h_min,
             const bool& is_non_reflBC,
             const bool& is_bed_friction,
             const bool& is_pore_water_pressure,
             const double& grav,
             const double& density_w,
             const double& density_s,
             const double& turbulence_coeff,
             const double& bed_friction_angle_rad,
             const double& erosion_coefficient,
             const double& m_coeff,
             const double& terminal_velocity,
                   double& odometric_coeff,
                   double& consolidation_coefficient,
             const double& thr_erodible_layer,
             const int& number_FD_points);
  
  TG2_scheme() = delete;
  
  ~TG2_scheme() = default;
  

  std::array<double,2>
  max_eigen (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean);
  
  void
  compute_dt (tmesh::quadrant_iterator quadrant);
  
  void
  compute_dt_adaptive (tmesh::quadrant_iterator quadrant);

  void
  solve_non_lin_h(const int& kk);

  void
  solve_non_lin_U(const int& kk);

  void
  Newton_mass_balance(double& hw_c, double& hs_c, const double& Uxw_c, const double& Uyw_c, const double& Uxs_c, const double& Uys_c, const double& tau_);

  void
  Newton_momentum_balance(const double& hw_c, const double& hs_c, double& Uxw_c, double& Uyw_c, double& Uxs_c, double& Uys_c, const double& bed_excess_pore_water_pressure, const double& tau_, const double& dt_);
  
  void
  first_step (tmesh::quadrant_iterator quadrant);

  void
  first_step_consolidation (tmesh::quadrant_iterator quadrant);

  void
  compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant);
  
  void
  second_step (tmesh::quadrant_iterator quadrant);

  void
  second_step_pressure (tmesh::quadrant_iterator quadrant, const int& kkk);

  void
  terminate_second_step (tmesh::quadrant_iterator quadrant);

  void
  solve_second_step_cons_equation (tmesh::quadrant_iterator quadrant);
  
  void
  flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& mass_node, double& phi_cell_Q);

  void
  set_dt (const double dt_);

  void
  set_tau ();

  void
  set_old_dt (const double dt_);
  
  void
  set_times(const double& time, const double& time_old, const double& time_oldd);

  void
  communication_part(tmesh::quadrant_iterator quadrant, std::vector<MPI_Request>& reqs, int& count_req, int& count_send, int& shift_send, int& shift_rec);

  void
  communication_part(tmesh::quadrant_iterator quadrant, int& shift_send, int& shift_rec);

  void
  loop_step (const int& kk, const bool& isInitial);

  void
  stabilization_term(const int& kk);

  template <class T>
  void
  numerical_integration_pressure(const int& kk, double& dp_mean, T& vv, const double& is_dof_or_hanging)
  {
    for (int kkk=kk; kkk<kk+number_FD_points-1; kkk++) 
    {
      dp_mean += .5*(vv.get_owned_data ()[kkk] + vv.get_owned_data ()[kkk+1]);
    }
    dp_mean /= double(number_FD_points-1);
    dp_mean *= is_dof_or_hanging;

    /*
    for (int kkk=kk; kkk<kk+number_FD_points-1; kkk+=2) 
    {
      //std::cout << vv.get_owned_data ()[kkk] << " " << vv.get_owned_data ()[kkk+1] << " " << vv.get_owned_data ()[kkk+2] << std::endl;
      dp_mean += (1./3.)*(vv.get_owned_data ()[kkk] + 4.*vv.get_owned_data ()[kkk+1] + vv.get_owned_data ()[kkk+2]);
    }
    dp_mean /= double(number_FD_points-1);
    dp_mean *= is_dof_or_hanging;*/
  }

  template <class T>
  void
  numerical_integration_pressure_2(const int& kk, double& dp_mean, T& vv, const double& is_dof_or_hanging)
  {
    for (int kkk=kk; kkk<kk+number_FD_points-1; kkk++) 
    {
      dp_mean += .5*(vv[kkk] + vv[kkk+1]);
    }
    dp_mean /= double(number_FD_points-1);
    dp_mean *= is_dof_or_hanging;

    /*
    for (int kkk=kk; kkk<kk+number_FD_points-1; kkk+=2) 
    {
      dp_mean += (1./3.)*(vv[kkk] + 4.*vv[kkk+1] + vv[kkk+2]);
    }
    dp_mean /= double(number_FD_points-1);
    dp_mean *= is_dof_or_hanging;*/
  }


  double
  linear_interpolation(const int& kk, const int& index_quadrant_c, const double& h_tot_old, const double& hdof_c, const double& Ux_tot_c, const double & Uy_tot_c);

  std::array<double,2>
  linear_interpolation(const int& kk, const double& h_tot_old, const double& hdof_c, const double& Ux_tot_c, const double & Uy_tot_c);

  double
  linear_interpolation(const int& kkk_ini, const double& delta_h_old, const double& Z_moved);

  double
  linear_interpolation_second(const double& Z_node_c, const double& delta_h_old, const double& Pxi, const double& Z_2, const int& kkk_ini, const double& hdofold_c, const double& vxs_cell, const double& vys_cell);
  
  double
  get_dt ();
  
  
  double dt, dt_old, tau, tau_c, tau_cc, tau_ccc;
  
  double Dx, Dy, area;
  
  
    ///  quadrant vertex (dofs) coordinates
    ///  The assumed numbering for quadrant nodes is
    ///  the following :
    ///    ^
    ///   yI
    ///   2------------------3
    ///   |                  |
    ///   |                  |
    ///   |                  |
    ///   |                  |
    ///   0------------------1 -->x
  
  std::array<double, 4> xn = {0, 0, 0, 0};
  std::array<double, 4> yn = {0, 0, 0, 0};
  
  
  // local dofs for state vector components
  std::array<double, 4> etawdof = {0, 0, 0, 0};
  std::array<double, 4> etasdof = {0, 0, 0, 0};
  std::array<double, 4> hdofold = {0, 0, 0, 0};
  std::array<double, 4> hwdofold = {0, 0, 0, 0};
  std::array<double, 4> hsdofold = {0, 0, 0, 0};
  std::array<double, 4> dp_kk_dof    = {0, 0, 0, 0};
  std::array<double, 4> hdof    = {0, 0, 0, 0};
  std::array<double, 4> hwdof   = {0, 0, 0, 0};
  std::array<double, 4> hsdof   = {0, 0, 0, 0};
  std::array<double, 4> Uxwdof  = {0, 0, 0, 0};
  std::array<double, 4> Uywdof  = {0, 0, 0, 0};
  std::array<double, 4> Uxsdof  = {0, 0, 0, 0};
  std::array<double, 4> Uysdof  = {0, 0, 0, 0};
  std::array<double, 4> hwdof_nei   = {0, 0, 0, 0};
  std::array<double, 4> hsdof_nei   = {0, 0, 0, 0};
  std::array<double, 4> Uxwdof_nei  = {0, 0, 0, 0};
  std::array<double, 4> Uywdof_nei  = {0, 0, 0, 0};
  std::array<double, 4> Uxsdof_nei  = {0, 0, 0, 0};
  std::array<double, 4> Uysdof_nei  = {0, 0, 0, 0};
  std::array<double, 4> Z_node  = {0, 0, 0, 0};
  std::array<double, 4> n_node  = {0, 0, 0, 0};
  std::array<double, 4> ns_node = {0, 0, 0, 0};
  std::array<double, 4> solid_vel_x = {0, 0, 0, 0};
  std::array<double, 4> solid_vel_y = {0, 0, 0, 0};
  std::array<double, 4> liquid_vel_x = {0, 0, 0, 0};
  std::array<double, 4> liquid_vel_y = {0, 0, 0, 0};
  std::array<double, 4> P_plus_hw_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_hw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_hs_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_hs_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uxw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uxw_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uyw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uyw_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uxs_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uxs_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uys_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uys_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_pressure_dp_kk_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_pressure_dp_kk_dof = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_hw_node    = {0, 0, 0, 0}, fluxy_hw_node    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_hs_node    = {0, 0, 0, 0}, fluxy_hs_node    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uxw_node   = {0, 0, 0, 0}, fluxy_Uxw_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uyw_node   = {0, 0, 0, 0}, fluxy_Uyw_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uxs_node   = {0, 0, 0, 0}, fluxy_Uxs_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uys_node   = {0, 0, 0, 0}, fluxy_Uys_node   = {0, 0, 0, 0};
  

  
  
  // flux functions
  double
  hw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  hw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  hs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  hs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean);
  
  double
  Uxw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean);

  double
  Uxs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean);
  
  double
  Uxs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uys_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uys_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean);
  



  // slope source terms
  double 
  src_slope_formula (const double& h, const double& S);

  double 
  src_slope_formula (const double& h, const double& S_x, const double& S_y, const int& kk);

  
  // source terms
  double
  hw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  hs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure);
  
  double
  Uys_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure);

  double
  Uxs_src_formula_2 (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure);
  
  double
  Uys_src_formula_2 (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure);

  double
  Uxw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  void
  thomas_algorithm(const double& mu_coeff, const int& kk_, const int& kk_ini, const int& kk_fin, const int& kk_ini_bottom);

  double
  signum (const double& x);

  void
  resize_vectors();

  void
  set_r_coeff();
  
  double time, timed, timedd;
  double nu_htot = 0.;
  double Fr = 0.;
  double g_coeff = 0.;

  double cfl_dp = .9;
  double nthr = 0.01;

  double sf = 1.2;

  double r_coeff = 0.;

  double error;
  const double tolerance = 1.e-8; 

  Q1& sol;
  Q1& sold;
  Q1& soldd;
  Q1& incr;
  std::vector<std::array<double,4>>& incr_anti_diff;
  std::vector<std::array<double,4>>& incr_anti_diff_pressure;
  Q1& P_plus;
  Q1& P_minus;
  Q1& P_plus_pressure;
  Q1& P_minus_pressure;
  Q0& sol_onehalf;
  Q1& Z;
  Q0& Z_onehalf;
  Q1& mass;
  Q1& excess_pore_water_pressure;
  Q0& excess_pore_water_pressure_onehalf;
  Q1& excess_pore_water_pressure_incr;
  
private:

  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging, der_coeffs_x, der_coeffs_y, der_coeffs_x_s, der_coeffs_y_s, D_U;
  std::array<double, 2> grad_cell_vsx, grad_cell_vsy, grad_cell_dp, grad_cell_Z, grad_cell_hw, grad_cell_hs, grad_cell_Uxw, grad_cell_Uyw, grad_cell_Uxs, grad_cell_Uys, grad_cell_dp_mean;
  
  const ordering& ordhw;
  const ordering& ordhs;
  const ordering& ordUxw;
  const ordering& ordUyw;
  const ordering& ordUxs;
  const ordering& ordUys;
  const ordering& ordBottom;
  const ordering& ordSurface;
  const double& DELTAT;
  const double& epsilon;
  const bool& is_non_reflBC;
  const double& grav;
  const double& erosion_coefficient;
  const double& density_w;
  const double& density_s;
  const double& turbulence_coeff;
  const bool& is_bed_friction;
  const bool& is_pore_water_pressure;
  const double& bed_friction_angle_rad;
  const double& m_coeff;
  const double& terminal_velocity;
  double& odometric_coeff;
  double& consolidation_coefficient;
  const double& thr_erodible_layer;
  const int& number_FD_points;


  int index_quadrant, index_quadrant_local;

  double hw_cell_average = 0., hs_cell_average = 0., erosion_contribution = 0.;

  std::array<double,4> dp_mean_vec = {0., 0., 0., 0.};

  std::vector<std::array<double,4>> contr_x;
  std::vector<std::array<double,4>> contr_y;

  std::vector<double> alfa_vec;
  std::vector<double> beta_vec;
  std::vector<double>    y_vec;

  // 5 arrays of storage as in Verwer's paper IMEX-RKCs,
  //std::vector<double> b_vect, mu_tilde_vect, gamma_tilde_vect, v_vect, mu_vect;

  const double tolerance_sign = 1.e-2; 

  const int Nmax = 1e3;
  int count;

  double dp_old = 0;
  std::array<double,2> grad_dp_old = {0,0};
  std::array<double,4> pressure = {0,0,0,0};

};




#endif

#include "Taylor_Galerkin_IMEX-RKC_Strang_balanced_two_phase.h"
#include <algorithm>
#include <cassert>

TG2_scheme::TG2_scheme(Q1& sol,
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
             const int& number_FD_points)
: sol(sol), 
sold(sold), 
soldd(soldd), 
incr(incr), 
incr_anti_diff(incr_anti_diff), 
incr_anti_diff_pressure(incr_anti_diff_pressure),
P_plus(P_plus), 
P_minus(P_minus), 
P_plus_pressure(P_plus_pressure), 
P_minus_pressure(P_minus_pressure), 
sol_onehalf(sol_onehalf), 
mass(mass), 
excess_pore_water_pressure(excess_pore_water_pressure),
excess_pore_water_pressure_onehalf(excess_pore_water_pressure_onehalf),
excess_pore_water_pressure_incr(excess_pore_water_pressure_incr),
ordhw(ohw), 
ordhs(ohs), 
ordUxw(oUxw), 
ordUyw(oUyw), 
ordUxs(oUxs), 
ordUys(oUys), 
ordBottom(oBottom),
ordSurface(oSurface),
Z(Z), 
Z_onehalf(Z_onehalf), 
DELTAT(DELTAT), 
epsilon(h_min), 
is_non_reflBC(is_non_reflBC), 
is_bed_friction(is_bed_friction), 
is_pore_water_pressure(is_pore_water_pressure),
grav(grav),
density_w(density_w), 
density_s(density_s), 
turbulence_coeff(turbulence_coeff), 
bed_friction_angle_rad(bed_friction_angle_rad), 
erosion_coefficient(erosion_coefficient), 
m_coeff(m_coeff), 
terminal_velocity(terminal_velocity),
odometric_coeff(odometric_coeff),
thr_erodible_layer(thr_erodible_layer),
consolidation_coefficient(consolidation_coefficient),
number_FD_points(number_FD_points)
{ }
 

std::array<double,2>
TG2_scheme::max_eigen (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean)
{
  double h = hw+hs,

         vel_x  = h >epsilon ? (Uxw+Uxs)/h : 0., 
         vel_y  = h >epsilon ? (Uyw+Uys)/h : 0., 
         velw_x = hw>epsilon ? Uxw/hw : 0., 
         velw_y = hw>epsilon ? Uyw/hw : 0., 
         vels_x = hs>epsilon ? Uxs/hs : 0., 
         vels_y = hs>epsilon ? Uys/hs : 0., 

         n = h>epsilon ? hw/h : 0.,

         celerity = std::sqrt(grav*h),
         beta_coeff = std::sqrt(.5*n*(1.-r_coeff)),
         beta_coeff_square = beta_coeff*beta_coeff,
         celerity_square = celerity*celerity,

         b_coeff = 4.*beta_coeff_square*dp_mean + 2.*celerity_square*density_w*(1.+beta_coeff_square),
         c_coeff = 16.*dp_mean*dp_mean*beta_coeff_square*beta_coeff_square + 16.*celerity_square*beta_coeff_square*beta_coeff_square*density_w*dp_mean + 4.*celerity_square*celerity_square*density_w*density_w*(1.-beta_coeff_square)*(1.-beta_coeff_square);



  const auto c_coeff_1 = .5*std::sqrt((b_coeff+std::sqrt(c_coeff))/density_w);

  const auto lambda_x = h>epsilon ? std::max(std::abs(velw_x), std::abs(vels_x))+c_coeff_1 : 0.;
  const auto lambda_y = h>epsilon ? std::max(std::abs(velw_y), std::abs(vels_y))+c_coeff_1 : 0.;

  //const auto lambda_x = h>epsilon ? std::max(std::abs(velw_x), std::abs(vels_x))+std::sqrt(grav*h) : 0.;
  //const auto lambda_y = h>epsilon ? std::max(std::abs(velw_y), std::abs(vels_y))+std::sqrt(grav*h) : 0.;

  //std::cout << lambda_x << " " << (h>epsilon ? std::max(std::abs(velw_x), std::abs(vels_x))+std::sqrt(grav*h) : 0.) << std::endl;

  return(std::array<double,2>{{lambda_x, lambda_y}});
}
 

void
TG2_scheme::compute_dt (tmesh::quadrant_iterator quadrant)
{
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii); 
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];

  dp_mean_vec = {0., 0., 0., 0.};
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hwdof [ii] = sol [ordhw  (quadrant->gt (ii) )];
      hsdof [ii] = sol [ordhs  (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      numerical_integration_pressure_2(ordBottom(quadrant->gt (ii)), dp_mean_vec[ii], excess_pore_water_pressure, 1.);

    }
    else
    {
      hwdof [ii] = .5 * (sol [ordhw  (quadrant->gparent (0, ii) )] +
                         sol [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof [ii] = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
      Uxwdof[ii] = .5 * (sol [ordUxw (quadrant->gparent (0, ii) )] +
                         sol [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof[ii] = .5 * (sol [ordUyw (quadrant->gparent (0, ii) )] +
                         sol [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);

      for (int jj=0; jj<=1; jj++)
        numerical_integration_pressure_2(ordBottom(quadrant->gparent(jj,ii)), dp_mean_vec[ii], excess_pore_water_pressure, .5);

    }

    const auto h_c  = hwdof [ii] + hsdof [ii];
    const auto Ux_c = Uxwdof[ii] + Uxsdof[ii];
    const auto Uy_c = Uywdof[ii] + Uysdof[ii];

    const double vel_abs = h_c>epsilon ? std::sqrt(Ux_c*Ux_c + Uy_c*Uy_c)/h_c : 0.; 
    Fr = std::max(h_c>epsilon ? vel_abs/std::sqrt(grav*h_c) : 0., Fr);
  }
  
  for (int ii = 0; ii < 4; ++ii){

    const auto lambdas = max_eigen (hwdof[ii], hsdof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii], dp_mean_vec[ii]);
    const auto & vel_rusanov_cell_x = lambdas[0];
    const auto & vel_rusanov_cell_y = lambdas[1];

    const auto h_c = hwdof[ii]+hsdof[ii];

    const auto dtoptx = h_c>epsilon ? Dx/vel_rusanov_cell_x : DELTAT;
    const auto dtopty = h_c>epsilon ? Dy/vel_rusanov_cell_y : DELTAT;
    const auto dtopt = dtoptx > dtopty ? dtopty : dtoptx;

    if (dt > dtopt) set_dt (dtopt);
  }
}

void
TG2_scheme::compute_dt_adaptive (tmesh::quadrant_iterator quadrant)
{
  
  double Nu_hmean_cell = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hwdof, hwdof_old, hwdof_oldd, hsdof, hsdof_old, hsdof_oldd, hdof, hdof_old, hdof_oldd, dh_t, h1, h2, h3, a_coeff, b_coeff;
    if (! quadrant->is_hanging (ii) )
    {
      hwdof      = sol   [ordhw  (quadrant->gt (ii) )];
      hwdof_old  = sold  [ordhw  (quadrant->gt (ii) )];
      hwdof_oldd = soldd [ordhw  (quadrant->gt (ii) )];

      hsdof      = sol   [ordhs  (quadrant->gt (ii) )];
      hsdof_old  = sold  [ordhs  (quadrant->gt (ii) )];
      hsdof_oldd = soldd [ordhs  (quadrant->gt (ii) )];
    }
    else
    {
      hwdof      = .5 * (sol   [ordhw  (quadrant->gparent (0, ii) )] +
                         sol   [ordhw  (quadrant->gparent (1, ii) )]);
      hwdof_old  = .5 * (sold  [ordhw  (quadrant->gparent (0, ii) )] +
                         sold  [ordhw  (quadrant->gparent (1, ii) )]);
      hwdof_oldd = .5 * (soldd [ordhw  (quadrant->gparent (0, ii) )] +
                         soldd [ordhw  (quadrant->gparent (1, ii) )]);

      hsdof      = .5 * (sol   [ordhs  (quadrant->gparent (0, ii) )] +
                         sol   [ordhs  (quadrant->gparent (1, ii) )]);
      hsdof_old  = .5 * (sold  [ordhs  (quadrant->gparent (0, ii) )] +
                         sold  [ordhs  (quadrant->gparent (1, ii) )]);
      hsdof_oldd = .5 * (soldd [ordhs  (quadrant->gparent (0, ii) )] +
                         soldd [ordhs  (quadrant->gparent (1, ii) )]);
    }
    hdof = hwdof+hsdof;
    hdof_old = hwdof_old+hsdof_old;
    hdof_oldd = hwdof_oldd+hsdof_oldd;

    dh_t = (hdof - hdof_old)/(time - timed);
    
    h1   = hdof_oldd/((timedd - timed )*(timedd - time ));
    h2   = hdof_old /((timed  - timedd)*(timed  - time ));
    h3   = hdof     /((time   - timedd)*(time   - timed));
    
    a_coeff = h1+h2+h3;
    b_coeff = - (h1*(time+timed) + h2*(time+timedd) + h3*(timed+timedd));
    
    Nu_hmean_cell += (1./3.*a_coeff*a_coeff*(time*time+time*timed+timed*timed) + a_coeff*(b_coeff-dh_t)*(time+timed) + (b_coeff-dh_t)*(b_coeff-dh_t));
    
  }
  Nu_hmean_cell /= 4.;
  nu_htot += Nu_hmean_cell*(time-timed)*(time-timed); // the dimension is length^2, this is eta^2
  
}

void
TG2_scheme::first_step_consolidation (tmesh::quadrant_iterator quadrant)
{

  index_quadrant = quadrant->get_global_quad_idx (); 
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }


  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;

  double solid_vel_x_average = 0., solid_vel_y_average = 0., liquid_vel_x_average = 0., liquid_vel_y_average = 0.;

  dp_mean_vec = {0., 0., 0., 0.};
  hw_cell_average = 0., hs_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hwdof [ii] = sol [ordhw  (quadrant->gt (ii) )];
      hsdof [ii] = sol [ordhs  (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      Z_node[ii]  = Z [quadrant->gt (ii)];

      numerical_integration_pressure_2(ordBottom(quadrant->gt (ii)), dp_mean_vec[ii], excess_pore_water_pressure, 1.);
    }
    else
    {
      hwdof [ii] = .5 * (sol [ordhw  (quadrant->gparent (0, ii) )] +
                         sol [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof [ii] = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
      Uxwdof[ii] = .5 * (sol [ordUxw (quadrant->gparent (0, ii) )] +
                         sol [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof[ii] = .5 * (sol [ordUyw (quadrant->gparent (0, ii) )] +
                         sol [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);

      Z_node[ii]  = .5 * (Z [quadrant->gparent(0,ii)] +
                          Z [quadrant->gparent(1,ii)]);

      for (int jj=0; jj<=1; jj++)
        numerical_integration_pressure_2(ordBottom(quadrant->gparent(jj,ii)), dp_mean_vec[ii], excess_pore_water_pressure, .5);

    }

    const auto & hwdof_c  = hwdof  [ii];
    const auto & hsdof_c  = hsdof  [ii];
    const auto & Uxwdof_c = Uxwdof [ii];
    const auto & Uywdof_c = Uywdof [ii];
    const auto & Uxsdof_c = Uxsdof [ii];
    const auto & Uysdof_c = Uysdof [ii];

    solid_vel_x[ii] = hsdof_c>epsilon ? Uxsdof_c/hsdof_c : 0.;
    solid_vel_y[ii] = hsdof_c>epsilon ? Uysdof_c/hsdof_c : 0.;

    liquid_vel_x[ii] = hwdof_c>epsilon ? Uxwdof_c/hwdof_c : 0.;
    liquid_vel_y[ii] = hwdof_c>epsilon ? Uywdof_c/hwdof_c : 0.;

    solid_vel_x_average += solid_vel_x[ii];
    solid_vel_y_average += solid_vel_y[ii];

    liquid_vel_x_average += liquid_vel_x[ii];
    liquid_vel_y_average += liquid_vel_y[ii];

    hw_cell_average += hwdof_c;
    hs_cell_average += hsdof_c;
    
    fluxx_hw_node [ii] = hw_flux_formula_x (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_hw_node [ii] = hw_flux_formula_y (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_hs_node [ii] = hs_flux_formula_x (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_hs_node [ii] = hs_flux_formula_y (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
  }
  hw_cell_average  /= 4.;
  hs_cell_average  /= 4.;

  solid_vel_x_average /= 4.;
  solid_vel_y_average /= 4.;

  liquid_vel_x_average /= 4.;
  liquid_vel_y_average /= 4.;


  const auto div_solid_x = .5*((solid_vel_x[1]-solid_vel_x[0]) + (solid_vel_x[3]-solid_vel_x[2]));
  const auto div_solid_y = .5*((solid_vel_y[2]-solid_vel_y[0]) + (solid_vel_y[3]-solid_vel_y[1]));
  const auto div_solid_vel = Dy*div_solid_x + Dx*div_solid_y;

  const auto div_liquid_x = .5*((liquid_vel_x[1]-liquid_vel_x[0]) + (liquid_vel_x[3]-liquid_vel_x[2]));
  const auto div_liquid_y = .5*((liquid_vel_y[2]-liquid_vel_y[0]) + (liquid_vel_y[3]-liquid_vel_y[1]));
  const auto div_liquid_vel = Dy*div_liquid_x + Dx*div_liquid_y;



  const auto grad_hw_x = .5*((hwdof[1] - hwdof[0]) + (hwdof[3] - hwdof[2]))/Dx;
  const auto grad_hw_y = .5*((hwdof[2] - hwdof[0]) + (hwdof[3] - hwdof[1]))/Dy;
  

  const auto div_Fhw_x = .5*((fluxx_hw_node[1]-fluxx_hw_node[0]) + (fluxx_hw_node[3]-fluxx_hw_node[2]));
  const auto div_Fhw_y = .5*((fluxy_hw_node[2]-fluxy_hw_node[0]) + (fluxy_hw_node[3]-fluxy_hw_node[1]));
  const auto div_Fhw_cell = Dy*div_Fhw_x + Dx*div_Fhw_y;

  const auto div_Fhs_x = .5*((fluxx_hs_node[1]-fluxx_hs_node[0]) + (fluxx_hs_node[3]-fluxx_hs_node[2]));
  const auto div_Fhs_y = .5*((fluxy_hs_node[2]-fluxy_hs_node[0]) + (fluxy_hs_node[3]-fluxy_hs_node[1]));
  const auto div_Fhs_cell = Dy*div_Fhs_x + Dx*div_Fhs_y;


  auto & hw_c = sol_onehalf[ordhw   (index_quadrant)];
  auto & hs_c = sol_onehalf[ordhs   (index_quadrant)];

  const auto v_hw = hw_cell_average - tau*div_Fhw_cell/area;
  const auto v_hs = hs_cell_average - tau*div_Fhs_cell/area;

  hw_c  = v_hw>0. ? v_hw : 0.;
  hs_c  = v_hs>0. ? v_hs : 0.;

  // add source term,
  const auto Uxw_c = (Uxwdof[0] + Uxwdof[1] + Uxwdof[2] + Uxwdof[3])*.25;
  const auto Uyw_c = (Uywdof[0] + Uywdof[1] + Uywdof[2] + Uywdof[3])*.25;
  const auto Uxs_c = (Uxsdof[0] + Uxsdof[1] + Uxsdof[2] + Uxsdof[3])*.25;
  const auto Uys_c = (Uysdof[0] + Uysdof[1] + Uysdof[2] + Uysdof[3])*.25;

  Newton_mass_balance(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, tau);

  // now solve the consolidation equation, first step,

  // compute the previous half step position, should be inside the quadrant because 1/2<1/sqrt(2)
  //std::vector<double> increment_array(number_FD_points);
  //increment_array.assign(number_FD_points, 0.);

  const double h_tot_c = hw_c + hs_c;
  const double delta_h_current = h_tot_c/(number_FD_points-1);
  const double n_s = h_tot_c>epsilon ? hs_c/h_tot_c : 0.;
  const double n_w = h_tot_c>epsilon ? hw_c/h_tot_c : 0.;
  const double density_prime = n_s*(density_s - density_w);

  //const double C_d = (hw_c>epsilon && h_tot_c>epsilon && hs_c>epsilon) ? n_w/terminal_velocity/std::pow(n_w, m_coeff)*density_prime*grav : 0.;
  //odometric_coeff = hw_c*hw_c>epsilon ? consolidation_coefficient*C_d/(n_w*n_w) : 0.;

  const double C_1_star = h_tot_c>epsilon ? odometric_coeff/h_tot_c : 0.;
  const double C_2_star = hs_c>epsilon && h_tot_c*hs_c>epsilon ? C_1_star*hw_c/hs_c : 0.;

  for (int kk=ordBottom(index_quadrant); kk<=ordSurface(index_quadrant); kk++) 
  {
    // Vertical ALE-material derivative contribution,

    const double z_current_mesh_node = (ordBottom(index_quadrant)!=0 ? kk%ordBottom(index_quadrant) : kk)*delta_h_current;
    //const double z_current_mesh_node = kk*delta_h_current;
    const double constant_ratio = h_tot_c>epsilon ? z_current_mesh_node/h_tot_c : 0.;
    const double vel_z_current_node = h_tot_c>epsilon ? erosion_contribution*(1. - constant_ratio) : 0.;

    double z_moved = z_current_mesh_node - tau*vel_z_current_node;
    z_moved = z_moved<0 ? 0. : z_moved; // here we are imposing BC

    //excess_pore_water_pressure_onehalf[kk] = 0.;

    
    pressure = {0,0,0,0};
    for (int ii = 0; ii < 4; ++ii)
    {

      if (! quadrant->is_hanging (ii) )
      {
        const double h_tot_old = sol [ordhw  (quadrant->gt (ii) )] + sol [ordhs  (quadrant->gt (ii) )];

        const double delta_h_old = h_tot_old/(number_FD_points-1);
        const double Z_moved = h_tot_c>epsilon ? z_moved*(h_tot_old/h_tot_c) : 0.;

        pressure[ii] = linear_interpolation(ordBottom(quadrant->gt (ii)), delta_h_old, Z_moved); 
      }
      else
      {
        for (int jj=0; jj<=1; jj++)
        {
          const double h_tot_old = sol [ordhw  (quadrant->gparent (jj, ii) )] + sol [ordhs  (quadrant->gparent (jj, ii) )];

          const double delta_h_old = h_tot_old/(number_FD_points-1);
          const double Z_moved = h_tot_c>epsilon ? z_moved*(h_tot_old/h_tot_c) : 0.;

          // the .5 is for the hanging node
          pressure[ii] += .5*linear_interpolation(ordBottom(quadrant->gparent (jj,ii)), delta_h_old, Z_moved); 
        }

      }
    }

    // mean over the cell
    excess_pore_water_pressure_onehalf[kk] = (pressure[0] + pressure[1] + pressure[2] + pressure[3])*.25;

    // add the horizontal transport contribution,
    std::array<double,2> grad_pres = {.5 * ( (pressure[3] - pressure[2]) + (pressure[1] - pressure[0]) )/Dx, .5 * ( (pressure[2] - pressure[0]) + (pressure[3] - pressure[1]) )/Dy};
    excess_pore_water_pressure_onehalf[kk] -= tau* (solid_vel_x_average*grad_pres[0] + solid_vel_y_average*grad_pres[1]);

    // solid_vel_x
    //std::array<double,2> grad_pres = {.5 * ( (pressure[3]*solid_vel_x[3] - pressure[2]*solid_vel_x[2]) + (pressure[1]*solid_vel_x[1] - pressure[0]*solid_vel_x[0]) )/Dx, .5 * ( (pressure[2]*solid_vel_y[2] - pressure[0]*solid_vel_y[0]) + (pressure[3]*solid_vel_y[3] - pressure[1]*solid_vel_y[1]) )/Dy};
    //excess_pore_water_pressure_onehalf[kk] -= tau* (grad_pres[0] + grad_pres[1]);

    // add the source term,
    const double C_1 = density_prime*grav*(1.-constant_ratio) -  C_1_star;
    const double C_2 = density_prime*grav*(1.-constant_ratio) +  C_2_star;
    excess_pore_water_pressure_onehalf[kk] += tau*(density_prime*grav*vel_z_current_node - C_1*hw_c*div_liquid_vel - C_2*hs_c*div_solid_vel - C_1*((liquid_vel_x_average - solid_vel_x_average)*grad_hw_x + (liquid_vel_y_average - solid_vel_y_average)*grad_hw_y)); //tau* ( density_prime*grav*vel_z_current_node - C_1*hw_c*div_liquid_vel - C_2*hs_c*div_solid_vel - C_1*((liquid_vel_x_average - solid_vel_x_average)*grad_hw_x + (liquid_vel_y_average - solid_vel_y_average)*grad_hw_y) );

    excess_pore_water_pressure_onehalf[kk] *= is_pore_water_pressure;
    
  }

}






void
TG2_scheme::first_step (tmesh::quadrant_iterator quadrant)
{
  index_quadrant = quadrant->get_global_quad_idx ();

  double dp_mean_cell = 0.;
  numerical_integration_pressure_2(ordBottom(index_quadrant), dp_mean_cell, excess_pore_water_pressure_onehalf, 1.);


  double Uxw_cell_average = 0., Uyw_cell_average = 0., Uxs_cell_average = 0., Uys_cell_average = 0., dp_mean_average = 0., h_cell_average = 0., nw_cell_average = 0., ns_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    const auto & hwdof_c   = hwdof  [ii];
    const auto & hsdof_c   = hsdof  [ii];
    const auto & Uxwdof_c  = Uxwdof [ii];
    const auto & Uywdof_c  = Uywdof [ii];
    const auto & Uxsdof_c  = Uxsdof [ii];
    const auto & Uysdof_c  = Uysdof [ii];
    const auto & dp_mean_c = dp_mean_vec[ii];

    Uxw_cell_average += Uxwdof_c;
    Uyw_cell_average += Uywdof_c;
    Uxs_cell_average += Uxsdof_c;
    Uys_cell_average += Uysdof_c;
    dp_mean_average  += dp_mean_c;

    fluxx_Uxw_node[ii] = Uxw_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c, dp_mean_cell);
    fluxy_Uxw_node[ii] = Uxw_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uyw_node[ii] = Uyw_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uyw_node[ii] = Uyw_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c, dp_mean_cell);
    fluxx_Uxs_node[ii] = Uxs_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c, dp_mean_cell);
    fluxy_Uxs_node[ii] = Uxs_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uys_node[ii] = Uys_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uys_node[ii] = Uys_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c, dp_mean_cell);
  }
  Uxw_cell_average /= 4.;
  Uyw_cell_average /= 4.;
  Uxs_cell_average /= 4.;
  Uys_cell_average /= 4.;
  dp_mean_average  /= 4.;


  h_cell_average  = hw_cell_average + hs_cell_average;
  nw_cell_average = h_cell_average>epsilon ? hw_cell_average/h_cell_average : 0.;
  ns_cell_average = h_cell_average>epsilon ? hs_cell_average/h_cell_average : 0.;
  
  
  const auto div_FUxw_x = .5*((fluxx_Uxw_node[1]-fluxx_Uxw_node[0]) + (fluxx_Uxw_node[3]-fluxx_Uxw_node[2]));
  const auto div_FUxw_y = .5*((fluxy_Uxw_node[2]-fluxy_Uxw_node[0]) + (fluxy_Uxw_node[3]-fluxy_Uxw_node[1]));
  const auto div_FUxw_cell = Dy*div_FUxw_x + Dx*div_FUxw_y;
  
  const auto div_FUyw_x = .5*((fluxx_Uyw_node[1]-fluxx_Uyw_node[0]) + (fluxx_Uyw_node[3]-fluxx_Uyw_node[2]));
  const auto div_FUyw_y = .5*((fluxy_Uyw_node[2]-fluxy_Uyw_node[0]) + (fluxy_Uyw_node[3]-fluxy_Uyw_node[1]));
  const auto div_FUyw_cell = Dy*div_FUyw_x + Dx*div_FUyw_y;

  const auto div_FUxs_x = .5*((fluxx_Uxs_node[1]-fluxx_Uxs_node[0]) + (fluxx_Uxs_node[3]-fluxx_Uxs_node[2]));
  const auto div_FUxs_y = .5*((fluxy_Uxs_node[2]-fluxy_Uxs_node[0]) + (fluxy_Uxs_node[3]-fluxy_Uxs_node[1]));
  const auto div_FUxs_cell = Dy*div_FUxs_x + Dx*div_FUxs_y;
  
  const auto div_FUys_x = .5*((fluxx_Uys_node[1]-fluxx_Uys_node[0]) + (fluxx_Uys_node[3]-fluxx_Uys_node[2]));
  const auto div_FUys_y = .5*((fluxy_Uys_node[2]-fluxy_Uys_node[0]) + (fluxy_Uys_node[3]-fluxy_Uys_node[1]));
  const auto div_FUys_cell = Dy*div_FUys_x + Dx*div_FUys_y;


  const auto slope_x_c = ((Z_node[1] - Z_node[0]) + (Z_node[3] - Z_node[2]))/Dx/2.;
  const auto slope_y_c = ((Z_node[2] - Z_node[0]) + (Z_node[3] - Z_node[1]))/Dy/2.;

  const auto grad_hw_x = ((hwdof[1] - hwdof[0]) + (hwdof[3] - hwdof[2]))/Dx/2.;
  const auto grad_hw_y = ((hwdof[2] - hwdof[0]) + (hwdof[3] - hwdof[1]))/Dy/2.;

  const auto grad_hs_x = ((hsdof[1] - hsdof[0]) + (hsdof[3] - hsdof[2]))/Dx/2.;
  const auto grad_hs_y = ((hsdof[2] - hsdof[0]) + (hsdof[3] - hsdof[1]))/Dy/2.;


  Z_onehalf[index_quadrant] = (Z_node[0]+Z_node[1]+Z_node[2]+Z_node[3])*.25;

  auto & hw_c  = sol_onehalf[ordhw   (index_quadrant)];
  auto & hs_c  = sol_onehalf[ordhs   (index_quadrant)];
  auto & Uxw_c = sol_onehalf[ordUxw  (index_quadrant)];
  auto & Uyw_c = sol_onehalf[ordUyw  (index_quadrant)];
  auto & Uxs_c = sol_onehalf[ordUxs  (index_quadrant)];
  auto & Uys_c = sol_onehalf[ordUys  (index_quadrant)];

  auto & bed_excess_pore_water_pressure = excess_pore_water_pressure_onehalf[ordBottom (index_quadrant)];

  Uxw_c = Uxw_cell_average - tau * (hw_cell_average>epsilon ? div_FUxw_cell/area - src_slope_formula (hw_cell_average, slope_x_c) +         (grav*hw_cell_average+nw_cell_average*dp_mean_average/density_w)*grad_hs_x - ns_cell_average*dp_mean_average/density_w*grad_hw_x : 0.);
  Uyw_c = Uyw_cell_average - tau * (hw_cell_average>epsilon ? div_FUyw_cell/area - src_slope_formula (hw_cell_average, slope_y_c) +         (grav*hw_cell_average+nw_cell_average*dp_mean_average/density_w)*grad_hs_y - ns_cell_average*dp_mean_average/density_w*grad_hw_y : 0.);
  Uxs_c = Uxs_cell_average - tau * (hs_cell_average>epsilon ? div_FUxs_cell/area - src_slope_formula (hs_cell_average, slope_x_c) + (r_coeff*grav*hs_cell_average+ns_cell_average*dp_mean_average/density_s)*grad_hw_x - nw_cell_average*dp_mean_average/density_s*grad_hs_x : 0.);
  Uys_c = Uys_cell_average - tau * (hs_cell_average>epsilon ? div_FUys_cell/area - src_slope_formula (hs_cell_average, slope_y_c) + (r_coeff*grav*hs_cell_average+ns_cell_average*dp_mean_average/density_s)*grad_hw_y - nw_cell_average*dp_mean_average/density_s*grad_hs_y : 0.);

  /*
  Uxw_c += tau_cc*Uxw_src_formula(hw_cell_average, hs_cell_average, Uxw_cell_average, Uyw_cell_average, Uxs_cell_average, Uys_cell_average);
  Uyw_c += tau_cc*Uyw_src_formula(hw_cell_average, hs_cell_average, Uxw_cell_average, Uyw_cell_average, Uxs_cell_average, Uys_cell_average);
  Uxs_c += tau_cc*Uxs_src_formula(hw_cell_average, hs_cell_average, Uxw_cell_average, Uyw_cell_average, Uxs_cell_average, Uys_cell_average, bed_excess_pore_water_pressure);
  Uys_c += tau_cc*Uys_src_formula(hw_cell_average, hs_cell_average, Uxw_cell_average, Uyw_cell_average, Uxs_cell_average, Uys_cell_average, bed_excess_pore_water_pressure);

  Newton_momentum_balance(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure, tau_c);
  */

  Newton_momentum_balance(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure, tau, tau);

}


void
TG2_scheme::Newton_mass_balance(double& hw_c, double& hs_c, const double& Uxw_c, const double& Uyw_c, const double& Uxs_c, const double& Uys_c, const double& tau_)
{
  // solve non-lin,
  const auto v_hw = hw_c;
  const auto v_hs = hs_c;

  const auto U_tot_x = Uxs_c + Uxw_c;
  const auto U_tot_y = Uys_c + Uyw_c;

  const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);
  erosion_contribution = erosion_coefficient*abs_mass_flux;


  count = -1;
  error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  {
    const auto h_c = hw_c+hs_c;
    const auto nw_c = h_c>epsilon ? hw_c/h_c : 0.;
    const auto ns_c = h_c>epsilon ? hs_c/h_c : 0.;
 

    // Newton Jacobian matrix, 
    // diagonalizzazione dei termini extra-diag delle matrici g e f
    const double big_A = 1.-(h_c>epsilon ? tau_*ns_c/h_c*erosion_contribution : 0.);
    const double big_B =     h_c>epsilon ? tau_*nw_c/h_c*erosion_contribution : 0.;
    const double big_C =     h_c>epsilon ? tau_*ns_c/h_c*erosion_contribution : 0.;
    const double big_D = 1.-(h_c>epsilon ? tau_*nw_c/h_c*erosion_contribution : 0.);

    // rhs terms of the Newton matrix,
    const auto f_hw = v_hw + tau_*hw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - hw_c;
    const auto f_hs = v_hs + tau_*hs_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - hs_c;

    const double big_det = big_A*big_D-big_B*big_C;

    const auto delta_hw = (f_hw*big_D-f_hs*big_B)/big_det;
    const auto delta_hs = (f_hs*big_A-f_hw*big_C)/big_det;
    
    error = std::sqrt(delta_hw*delta_hw + delta_hs*delta_hs);

    hw_c += delta_hw;
    hs_c += delta_hs;
  }

  if (error>tolerance)
  {
    std::cout << "No convergence erosion!! " << error << std::endl;
  }
}



void
TG2_scheme::Newton_momentum_balance(const double& hw_c, const double& hs_c, double& Uxw_c, double& Uyw_c, double& Uxs_c, double& Uys_c, const double& bed_excess_pore_water_pressure, const double& tau_, const double& dt_)
{

  const double h_c = hw_c+hs_c;
  
  const double n_c  = h_c>epsilon ? hw_c/h_c : 0.;
  const double ns_c = h_c>epsilon ? hs_c/h_c : 0.;
  
  const double common_coeff = hw_c>epsilon ? dt_/terminal_velocity*(density_s - density_w)*grav/std::pow(n_c, m_coeff) : 0.;

  const double a_coeff = common_coeff*ns_c;
  const double b_coeff = common_coeff*n_c;


  const auto v_Ux_w = Uxw_c;
  const auto v_Uy_w = Uyw_c;
  const auto v_Ux_s = Uxs_c;
  const auto v_Uy_s = Uys_c;

  /*
  const double density_prime = ns_c*(density_s - density_w);
  const double density = ns_c*density_s + n_c*density_w;

  const double delta_coeff = (density_prime*grav*h_c - bed_excess_pore_water_pressure)*std::tan(bed_friction_angle_rad);
  const double gamma_coeff = (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/h_c/h_c : 0.;
*/

  const auto A_11 = 1. + a_coeff/density_w;
  const auto A_22 = 1. + b_coeff/density_s;
  const auto A_21 =    - a_coeff/density_s;
  const auto A_12 =    - b_coeff/density_w;


  const auto det_A = A_11*A_22 - A_12*A_21;

  const auto Ainv_11 = A_22/det_A;
  const auto Ainv_22 = A_11/det_A;
  const auto Ainv_21 =-A_21/det_A;
  const auto Ainv_12 =-A_12/det_A;
  
  
  // solve non-linearities
  double Uxw_new_c, Uyw_new_c, Uxs_new_c, Uys_new_c;
  count = -1;
  error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  { 
    const auto U_tot_x = Uxs_c + Uxw_c;
    const auto U_tot_y = Uys_c + Uyw_c;


    Uxw_new_c = Ainv_11*v_Ux_w + Ainv_12*(v_Ux_s + tau_*Uxs_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure));
    Uyw_new_c = Ainv_11*v_Uy_w + Ainv_12*(v_Uy_s + tau_*Uys_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure));

    Uxs_new_c = Ainv_21*v_Ux_w + Ainv_22*(v_Ux_s + tau_*Uxs_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure));
    Uys_new_c = Ainv_21*v_Uy_w + Ainv_22*(v_Uy_s + tau_*Uys_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure));

    const auto delta_Uwx = Uxw_new_c - Uxw_c;
    const auto delta_Uwy = Uyw_new_c - Uyw_c;
    const auto delta_Usx = Uxs_new_c - Uxs_c;
    const auto delta_Usy = Uys_new_c - Uys_c;

    error = std::sqrt(delta_Uwx*delta_Uwx + delta_Uwy*delta_Uwy + delta_Usx*delta_Usx + delta_Usy*delta_Usy);
 
    Uxw_c = Uxw_new_c;
    Uyw_c = Uyw_new_c;
    Uxs_c = Uxs_new_c;
    Uys_c = Uys_new_c;
  }

  if (error>tolerance)
  {
    std::cout << "No convergence momentum!! " << error << std::endl;
  }

/*
    const double delta_x = 1./(1.+U_tot_x*U_tot_x*(M_PI*.5)*(M_PI*.5)); //std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign;
    const double delta_y = 1./(1.+U_tot_y*U_tot_y*(M_PI*.5)*(M_PI*.5)); //std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign;

    const auto J_11_x = A_11;
    const auto J_12_x = A_12;
    const auto J_22_x = A_22 - tau_*(- delta_coeff*delta_x - 2.*gamma_coeff*std::abs(U_tot_x))/density_s;
    const auto J_21_x = A_21 - tau_*(- delta_coeff*delta_x - 2.*gamma_coeff*std::abs(U_tot_x))/density_s;

    const auto J_11_y = A_11;
    const auto J_12_y = A_12;
    const auto J_22_y = A_22 - tau_*(- delta_coeff*delta_y - 2.*gamma_coeff*std::abs(U_tot_y))/density_s;
    const auto J_21_y = A_21 - tau_*(- delta_coeff*delta_y - 2.*gamma_coeff*std::abs(U_tot_y))/density_s;


    const auto rhs_1_x = 
    const auto rhs_2_x = 

    const auto rhs_1_y = 
    const auto rhs_2_y = 

    const auto f_1 = tau_*Uxs_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure);
*/


/*

    const auto numerator_x = tau_/(1.+a_coeff)*v_Ux_w*(-delta_coeff*delta_x - 2.*gamma_coeff*std::abs(U_tot_x)) - a_coeff/(1.+a_coeff)*v_Ux_w - v_Ux_s - tau_*Uxs_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure);
    const auto denominator_x = 1.+b_coeff - (b_coeff*a_coeff)/(1.+a_coeff) + tau_*(1.+b_coeff/(1.+a_coeff))*(-delta_coeff*delta_x-2.*gamma_coeff*std::abs(U_tot_x);
    
    const auto numerator_y = tau_/(1.+a_coeff)*v_Uy_w*(-delta_coeff*delta_y - 2.*gamma_coeff*std::abs(U_tot_y)) - a_coeff/(1.+a_coeff)*v_Uy_w - v_Uy_s - tau_*Uys_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure);
    const auto denominator_y = 1.+b_coeff - (b_coeff*a_coeff)/(1.+a_coeff) + tau_*(1.+b_coeff/(1.+a_coeff))*(-delta_coeff*delta_y-2.*gamma_coeff*std::abs(U_tot_y);

    
    const double delta_Usx = numerator_x/denominator_x;
    const double delta_Usy = numerator_y/denominator_y;
 
    Uxs_c += delta_Usx;
    Uys_c += delta_Usy;


  }
*/

/*
  const double h_c = hw_c+hs_c;

  const double n_c  = h_c>epsilon ? hw_c/h_c : 0.;
  const double ns_c = h_c>epsilon ? hs_c/h_c : 0.;

  const double density = ns_c*density_s + n_c*density_w;
  const double density_prime = ns_c*(density_s-density_w);

  const double C_d = (hw_c>epsilon && h_c>epsilon && hs_c>epsilon) ? n_c/terminal_velocity/std::pow(n_c, m_coeff)*density_prime*grav : 0.;

  double big_Ax = 1. + (hw_c>epsilon ? tau_*C_d*h_c/hw_c/density_w : 0.);
  double big_Bx = hs_c>epsilon ? -tau_*C_d*h_c/hs_c/density_w : 0.;
  double big_Cx = hw_c>epsilon ? -tau_*C_d*h_c/hw_c/density_s : 0.;
  double big_Dx = 1. + (hs_c>epsilon ? tau_*C_d*h_c/hs_c/density_s : 0.);

  const double rhs_1x = Uxw_c;
  const double rhs_2x = Uxs_c;

  double big_detx = big_Ax*big_Dx-big_Bx*big_Cx; 

  double big_Ay = 1. + (hw_c>epsilon ? tau_*C_d*h_c/hw_c/density_w : 0.);
  double big_By = hs_c>epsilon ? -tau_*C_d*h_c/hs_c/density_w : 0.;
  double big_Cy = hw_c>epsilon ? -tau_*C_d*h_c/hw_c/density_s : 0.;
  double big_Dy = 1. + (hs_c>epsilon ? tau_*C_d*h_c/hs_c/density_s : 0.);

  const double rhs_1y = Uyw_c;
  const double rhs_2y = Uys_c;

  double big_dety = big_Ay*big_Dy-big_By*big_Cy;

  double big_dety = big_Ay*big_Dy-big_By*big_Cy;

  Uxw_c = (rhs_1x*big_Dx-rhs_2x*big_Bx)/big_detx;
  Uyw_c = (rhs_1y*big_Dy-rhs_2y*big_By)/big_dety;
  Uxs_c = (rhs_2x*big_Ax-rhs_1x*big_Cx)/big_detx;
  Uys_c = (rhs_2y*big_Ay-rhs_1y*big_Cy)/big_dety;

  const auto v_Ux_w = Uxw_c;
  const auto v_Uy_w = Uyw_c;
  const auto v_Ux_s = Uxs_c;
  const auto v_Uy_s = Uys_c;


  // solve non-linearities
  count = -1;
  error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  {
    const auto U_tot_x = Uxs_c + Uxw_c;
    const auto U_tot_y = Uys_c + Uyw_c;

    const double delta_x = 1./(1.+U_tot_x*U_tot_x*(M_PI*.5)*(M_PI*.5)); //std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign;
    const double delta_y = 1./(1.+U_tot_y*U_tot_y*(M_PI*.5)*(M_PI*.5)); //std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign;

    const double fric_x = (h_c*h_c)>epsilon ? (density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x)) : 0.;
    const double fric_y = (h_c*h_c)>epsilon ? (density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y)) : 0.;


    big_Dx = 1. + tau_* (is_bed_friction && hs_c>epsilon ? ((density_prime*grav*h_c-bed_excess_pore_water_pressure)*std::tan(bed_friction_angle_rad)*delta_x + fric_x) : 0.) /density_s;
    big_Dy = 1. + tau_* (is_bed_friction && hs_c>epsilon ? ((density_prime*grav*h_c-bed_excess_pore_water_pressure)*std::tan(bed_friction_angle_rad)*delta_y + fric_y) : 0.) /density_s;


    const double f_Usx = - Uxs_c + v_Ux_s + tau_*Uxs_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure);
    const double f_Usy = - Uys_c + v_Uy_s + tau_*Uys_src_formula_2(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure);
    
    const double delta_Usx = f_Usx/big_Dx;
    const double delta_Usy = f_Usy/big_Dy;

    error = std::sqrt(delta_Usx*delta_Usx + delta_Usy*delta_Usy);

    Uxs_c += delta_Usx; 
    Uys_c += delta_Usy;

  }

  if (error>tolerance)
  {
    std::cout << "No convergence momentum!! " << error << std::endl;
  }
*/
}




void
TG2_scheme::solve_non_lin_h(const int& kk)
{
  auto & hw_c  = sol.get_owned_data ()[kk  ];
  auto & hs_c  = sol.get_owned_data ()[kk+1];

  const auto & hw_cc  = sold.get_owned_data ()[kk  ];
  const auto & hs_cc  = sold.get_owned_data ()[kk+1];
  const auto & Uxw_cc = sold.get_owned_data ()[kk+2];
  const auto & Uyw_cc = sold.get_owned_data ()[kk+3];
  const auto & Uxs_cc = sold.get_owned_data ()[kk+4];
  const auto & Uys_cc = sold.get_owned_data ()[kk+5];

  // solve non-linearities like the first step of the TG2 method to get the complete low order solution,
  // terminated this part add the nodal correction to get the hyperbolicity,
  hw_c += dt*incr.get_owned_data ()[kk  ]/mass.get_owned_data ()[kk  ] + tau*hw_src_formula(hw_cc, hs_cc, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc);
  hs_c += dt*incr.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] + tau*hs_src_formula(hw_cc, hs_cc, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc);

  hw_c *= (hw_c>0);
  hs_c *= (hs_c>0);

  Newton_mass_balance(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc, tau);

}

void
TG2_scheme::solve_non_lin_U(const int& kk)
{
  const auto & hw_c  = sol.get_owned_data ()[kk  ];
  const auto & hs_c  = sol.get_owned_data ()[kk+1];
        auto & Uxw_c = sol.get_owned_data ()[kk+2];
        auto & Uyw_c = sol.get_owned_data ()[kk+3];
        auto & Uxs_c = sol.get_owned_data ()[kk+4];
        auto & Uys_c = sol.get_owned_data ()[kk+5];


  /*
  //const auto & hw_cc  = sold.get_owned_data ()[kk  ];
  //const auto & hs_cc  = sold.get_owned_data ()[kk+1];
  const auto Uxw_cc = Uxw_c; //sold.get_owned_data ()[kk+2];
  const auto Uyw_cc = Uyw_c; //sold.get_owned_data ()[kk+3];
  const auto Uxs_cc = Uxs_c; //sold.get_owned_data ()[kk+4];
  const auto Uys_cc = Uys_c; //sold.get_owned_data ()[kk+5];
*/
 
  const auto Uxw_cc = sold.get_owned_data ()[kk+2];
  const auto Uyw_cc = sold.get_owned_data ()[kk+3];
  const auto Uxs_cc = sold.get_owned_data ()[kk+4];
  const auto Uys_cc = sold.get_owned_data ()[kk+5];
 
 
  const auto & bed_excess_pore_water_pressure = excess_pore_water_pressure.get_owned_data ()[(kk/6)*number_FD_points];


  Uxw_c += dt*incr.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2];// + tau*Uxw_src_formula_2(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc);
  Uyw_c += dt*incr.get_owned_data ()[kk+3]/mass.get_owned_data ()[kk+3];// + tau*Uyw_src_formula_2(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc);
  Uxs_c += dt*incr.get_owned_data ()[kk+4]/mass.get_owned_data ()[kk+4] + tau*Uxs_src_formula_2(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc, bed_excess_pore_water_pressure);
  Uys_c += dt*incr.get_owned_data ()[kk+5]/mass.get_owned_data ()[kk+5] + tau*Uys_src_formula_2(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc, bed_excess_pore_water_pressure);

  Newton_momentum_balance(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure, tau, dt);

/*
  Uxw_c += dt*incr.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + tau_c*Uxw_src_formula(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc);
  Uyw_c += dt*incr.get_owned_data ()[kk+3]/mass.get_owned_data ()[kk+3] + tau_c*Uyw_src_formula(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc);
  Uxs_c += dt*incr.get_owned_data ()[kk+4]/mass.get_owned_data ()[kk+4] + tau_c*Uxs_src_formula(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc, bed_excess_pore_water_pressure);
  Uys_c += dt*incr.get_owned_data ()[kk+5]/mass.get_owned_data ()[kk+5] + tau_c*Uys_src_formula(hw_c, hs_c, Uxw_cc, Uyw_cc, Uxs_cc, Uys_cc, bed_excess_pore_water_pressure);

  Newton_momentum_balance(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c, bed_excess_pore_water_pressure, tau_c);
*/


}



void
TG2_scheme::solve_second_step_cons_equation (tmesh::quadrant_iterator quadrant)
{
  // look at tmesh.h 
  index_quadrant = quadrant->get_global_quad_idx (); 
  index_quadrant_local = quadrant->get_forest_quad_idx ();

  std::array<int,4> bimpp_to_rev_ord = {0, 1, 3, 2};
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;


  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hwdof [ii] = sol  [ordhw  (quadrant->gt (ii) )];
      hsdof [ii] = sol  [ordhs  (quadrant->gt (ii) )];

      hwdofold [ii] = sold [ordhw  (quadrant->gt (ii) )];
      hsdofold [ii] = sold [ordhs  (quadrant->gt (ii) )];
      Uxwdof   [ii] = sold [ordUxw (quadrant->gt (ii) )];
      Uywdof   [ii] = sold [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof   [ii] = sold [ordUxs (quadrant->gt (ii) )];
      Uysdof   [ii] = sold [ordUys (quadrant->gt (ii) )];

      isdof_or_hanging[ii] = 1.;
      
    }
    else
    {
      hwdof [ii] = .5 * (sol  [ordhw  (quadrant->gparent (0, ii) )] +
                         sol  [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof [ii] = .5 * (sol  [ordhs  (quadrant->gparent (0, ii) )] +
                         sol  [ordhs  (quadrant->gparent (1, ii) )]);

      hwdofold [ii] = .5 * (sold [ordhw  (quadrant->gparent (0, ii) )] +
                            sold [ordhw  (quadrant->gparent (1, ii) )]);
      hsdofold [ii] = .5 * (sold [ordhs  (quadrant->gparent (0, ii) )] +
                            sold [ordhs  (quadrant->gparent (1, ii) )]);
      Uxwdof   [ii] = .5 * (sold [ordUxw (quadrant->gparent (0, ii) )] +
                            sold [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof   [ii] = .5 * (sold [ordUyw (quadrant->gparent (0, ii) )] +
                            sold [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof   [ii] = .5 * (sold [ordUxs (quadrant->gparent (0, ii) )] +
                            sold [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof   [ii] = .5 * (sold [ordUys (quadrant->gparent (0, ii) )] +
                            sold [ordUys (quadrant->gparent (1, ii) )]);

      isdof_or_hanging[ii] = .5;
    }

    hdof   [ii] = hwdof   [ii] + hsdof   [ii];
    hdofold[ii] = hwdofold[ii] + hsdofold[ii];
  } 

  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
    -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
   
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
    +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};

  const double & hw_cell    = sol_onehalf[ordhw    (index_quadrant)];
  const double & Uxw_cell   = sol_onehalf[ordUxw   (index_quadrant)];
  const double & Uyw_cell   = sol_onehalf[ordUyw   (index_quadrant)];

  const double & hs_cell    = sol_onehalf[ordhs    (index_quadrant)];
  const double & Uxs_cell   = sol_onehalf[ordUxs   (index_quadrant)];
  const double & Uys_cell   = sol_onehalf[ordUys   (index_quadrant)];


  const auto h_cell = hw_cell + hs_cell;

  const double delta_h_cell = h_cell/(number_FD_points-1);

  const auto vxs_cell = hs_cell>epsilon ? Uxs_cell/hs_cell : 0.;
  const auto vys_cell = hs_cell>epsilon ? Uys_cell/hs_cell : 0.;

  const auto nw_cell = h_cell>epsilon ? hw_cell/h_cell : 0.; 
  const auto ns_cell = h_cell>epsilon ? hs_cell/h_cell : 0.; 
  const double density_prime = ns_cell*(density_s-density_w);

  //const double C_d = (hw_cell>epsilon && h_cell>epsilon && hs_cell>epsilon) ? nw_cell/terminal_velocity/std::pow(nw_cell, m_coeff)*density_prime*grav : 0.;
  //odometric_coeff = hw_cell*hw_cell>epsilon ? consolidation_coefficient*C_d/(nw_cell*nw_cell) : 0.;

  auto C = [&] (const double& z)
  {
    return(density_prime*grav*(1.- (h_cell>epsilon ? z/h_cell : 0.)));
  };

  auto C_1 = [&] (const double& z)
  {
    return(C(z) - (h_cell>epsilon ? odometric_coeff/h_cell : 0.));
  };

  auto C_2 = [&] (const double& z)
  {
    return(C(z) + ((hs_cell>epsilon && h_cell*hs_cell>epsilon) ? odometric_coeff/h_cell*hw_cell/hs_cell : 0.));
  };


  contr_x.assign(number_FD_points, std::array<double,4>{{0,0,0,0}});
  contr_y.assign(number_FD_points, std::array<double,4>{{0,0,0,0}});


  // PC for the horizontal transport term,
  bool is_boundary_edge = true;
  for (int iEdge = 0; iEdge < 4; ++iEdge){

    is_boundary_edge = true;

    const auto i_1 = bimpp_to_rev_ord[iEdge];
    const auto i_2 = bimpp_to_rev_ord[(iEdge+1)%4];

    const auto edge_length = std::sqrt(std::pow((xn[i_1]-xn[i_2]),2.) + std::pow((yn[i_1]-yn[i_2]),2.));
    const std::array<double,2> outward_normal_edge = {(-yn[i_1]+yn[i_2])/edge_length, ( xn[i_1]-xn[i_2])/edge_length};  

    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {
      std::array<double,4> Xn, Yn;

      for (int ii = 0; ii < 4; ++ii) {
        Xn[ii] = quadrant_nei->p(0, ii);
        Yn[ii] = quadrant_nei->p(1, ii);
      }

      const auto & index_quadrant_nei = quadrant_nei->get_global_quad_idx ();


      for (int jEdge = 0; jEdge < 4; ++jEdge) { // cycle neigh edges 

        const auto j_1 = bimpp_to_rev_ord[jEdge];
        const auto j_2 = bimpp_to_rev_ord[(jEdge+1)%4];

        const auto edge_length_nei = std::sqrt(std::pow((Xn[j_1]-Xn[j_2]),2.) + std::pow((Yn[j_1]-Yn[j_2]),2.));
        const std::array<double,2> outward_normal_edge_nei = {(-Yn[j_1]+Yn[j_2])/edge_length_nei, ( Xn[j_1]-Xn[j_2])/edge_length_nei};  
        const bool check_orthogonality = std::inner_product(outward_normal_edge_nei.begin(), outward_normal_edge_nei.end(), outward_normal_edge.begin(), 0.) == -1;


        if ( (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) ||
               (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
              ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) ||
               (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant!=index_quadrant_nei )
        {
          is_boundary_edge = false;

          // scrivere qui la somma dei contributi per i termini non-cons.!
          double h_cell_nei, hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei;

          hw_cell_nei  = sol_onehalf[ordhw    (index_quadrant_nei)];
          hs_cell_nei  = sol_onehalf[ordhs    (index_quadrant_nei)];
          Uxw_cell_nei = sol_onehalf[ordUxw   (index_quadrant_nei)];
          Uyw_cell_nei = sol_onehalf[ordUyw   (index_quadrant_nei)];
          Uxs_cell_nei = sol_onehalf[ordUxs   (index_quadrant_nei)];
          Uys_cell_nei = sol_onehalf[ordUys   (index_quadrant_nei)];

          h_cell_nei = hw_cell_nei + hs_cell_nei;

          const auto vxs_cell_nei = hs_cell_nei>epsilon ? Uxs_cell_nei/hs_cell_nei : 0.;
          const auto vys_cell_nei = hs_cell_nei>epsilon ? Uys_cell_nei/hs_cell_nei : 0.;

          for (int kk=0; kk<number_FD_points; kk++)
          {
            const double z = kk*delta_h_cell;

            const auto dp_cell_i_1     = linear_interpolation(kk, index_quadrant,     h_cell,     hdof[i_1], Uxwdof[i_1]+Uxsdof[i_1], Uywdof[i_1]+Uysdof[i_1]);
            const auto dp_cell_i_2     = linear_interpolation(kk, index_quadrant,     h_cell,     hdof[i_2], Uxwdof[i_2]+Uxsdof[i_2], Uywdof[i_2]+Uysdof[i_2]);

            const auto dp_cell_nei_i_1 = linear_interpolation(kk, index_quadrant_nei, h_cell_nei, hdof[i_1], Uxwdof[i_1]+Uxsdof[i_1], Uywdof[i_1]+Uysdof[i_1]);
            const auto dp_cell_nei_i_2 = linear_interpolation(kk, index_quadrant_nei, h_cell_nei, hdof[i_2], Uxwdof[i_2]+Uxsdof[i_2], Uywdof[i_2]+Uysdof[i_2]);

            //contr_x[kk][i_1] += .5*signum(outward_normal_edge[0])*(vxs_cell*(dp_cell_nei_i_1 - dp_cell_i_1) + C_2(z)*hs_cell*(vxs_cell_nei - vxs_cell) + C_1(z)*(Uxw_cell_nei - Uxw_cell) - C_1(z)*vxs_cell*(hw_cell_nei-hw_cell) )*isdof_or_hanging[i_1];
            //contr_x[kk][i_2] += .5*signum(outward_normal_edge[0])*(vxs_cell*(dp_cell_nei_i_2 - dp_cell_i_2) + C_2(z)*hs_cell*(vxs_cell_nei - vxs_cell) + C_1(z)*(Uxw_cell_nei - Uxw_cell) - C_1(z)*vxs_cell*(hw_cell_nei-hw_cell) )*isdof_or_hanging[i_2];

            //contr_y[kk][i_1] += .5*signum(outward_normal_edge[1])*(vys_cell*(dp_cell_nei_i_1 - dp_cell_i_1) + C_2(z)*hs_cell*(vys_cell_nei - vys_cell) + C_1(z)*(Uyw_cell_nei - Uyw_cell) - C_1(z)*vys_cell*(hw_cell_nei-hw_cell) )*isdof_or_hanging[i_1];
            //contr_y[kk][i_2] += .5*signum(outward_normal_edge[1])*(vys_cell*(dp_cell_nei_i_2 - dp_cell_i_2) + C_2(z)*hs_cell*(vys_cell_nei - vys_cell) + C_1(z)*(Uyw_cell_nei - Uyw_cell) - C_1(z)*vys_cell*(hw_cell_nei-hw_cell) )*isdof_or_hanging[i_2];

            //if (dp_cell_i_1!=excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)+kk] || dp_cell_i_2!=excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)+kk])
            //std::cout << dp_cell_i_1 << " " << dp_cell_i_2 << " " << excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)+kk] << " " << dp_cell_i_1-excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)+kk] << std::endl;

            // .5 salta fuori dall'integrazione per trapezi tra 0 e 1 in coordinata \xi (è il valore in LHS da metter qui sotto!) 
            contr_x[kk][i_1] += .5*signum(outward_normal_edge[0])*(h_cell>epsilon ? -dp_cell_i_1*(vxs_cell_nei - vxs_cell) + C_2(z)*hs_cell*(vxs_cell_nei - vxs_cell) + C_1(z)*(Uxw_cell_nei - Uxw_cell) - C_1(z)*vxs_cell*(hw_cell_nei-hw_cell) : 0. )*isdof_or_hanging[i_1];
            contr_x[kk][i_2] += .5*signum(outward_normal_edge[0])*(h_cell>epsilon ? -dp_cell_i_2*(vxs_cell_nei - vxs_cell) + C_2(z)*hs_cell*(vxs_cell_nei - vxs_cell) + C_1(z)*(Uxw_cell_nei - Uxw_cell) - C_1(z)*vxs_cell*(hw_cell_nei-hw_cell) : 0. )*isdof_or_hanging[i_2];

            contr_y[kk][i_1] += .5*signum(outward_normal_edge[1])*(h_cell>epsilon ? -dp_cell_i_1*(vys_cell_nei - vys_cell) + C_2(z)*hs_cell*(vys_cell_nei - vys_cell) + C_1(z)*(Uyw_cell_nei - Uyw_cell) - C_1(z)*vys_cell*(hw_cell_nei-hw_cell) : 0. )*isdof_or_hanging[i_1];
            contr_y[kk][i_2] += .5*signum(outward_normal_edge[1])*(h_cell>epsilon ? -dp_cell_i_2*(vys_cell_nei - vys_cell) + C_2(z)*hs_cell*(vys_cell_nei - vys_cell) + C_1(z)*(Uyw_cell_nei - Uyw_cell) - C_1(z)*vys_cell*(hw_cell_nei-hw_cell) : 0. )*isdof_or_hanging[i_2];

            //contr_x[kk][i_1] += .5*signum(outward_normal_edge[0])*(-dp_cell_i_1*(vxs_cell_nei - vxs_cell) - C_1(z)*vxs_cell*(hw_cell_nei - hw_cell)*0 - C_2(z)*vxs_cell*(hs_cell_nei - hs_cell) )*isdof_or_hanging[i_1];
            //contr_x[kk][i_2] += .5*signum(outward_normal_edge[0])*(-dp_cell_i_2*(vxs_cell_nei - vxs_cell) - C_1(z)*vxs_cell*(hw_cell_nei - hw_cell)*0 - C_2(z)*vxs_cell*(hs_cell_nei - hs_cell) )*isdof_or_hanging[i_2];

            //contr_y[kk][i_1] += .5*signum(outward_normal_edge[1])*(-dp_cell_i_1*(vys_cell_nei - vys_cell) - C_1(z)*vys_cell*(hw_cell_nei - hw_cell)*0 - C_2(z)*vys_cell*(hs_cell_nei - hs_cell) )*isdof_or_hanging[i_1];
            //contr_y[kk][i_2] += .5*signum(outward_normal_edge[1])*(-dp_cell_i_2*(vys_cell_nei - vys_cell) - C_1(z)*vys_cell*(hw_cell_nei - hw_cell)*0 - C_2(z)*vys_cell*(hs_cell_nei - hs_cell) )*isdof_or_hanging[i_2];

            //std::cout << C_1(z) << " " << C_2(z) << std::endl;
          }

        }
      }
    }
      
    if (is_boundary_edge) // set boundary conditions
    { 

      auto hs_cell_nei  = hs_cell;
      auto Uxs_cell_nei = Uxs_cell;
      auto Uys_cell_nei = Uys_cell;

      Uxs_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[0];
      Uys_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[1];

      const auto vxs_cell_nei = hs_cell_nei>epsilon ? Uxs_cell_nei/hs_cell_nei : 0.;
      const auto vys_cell_nei = hs_cell_nei>epsilon ? Uys_cell_nei/hs_cell_nei : 0.;

      const std::array<double,2> speed     = {std::abs(vxs_cell    ), std::abs(vys_cell    )};
      const std::array<double,2> speed_nei = {std::abs(vxs_cell_nei), std::abs(vys_cell_nei)};

      const auto smax = std::max(speed[0]*outward_normal_edge[0]+speed[1]*outward_normal_edge[1], speed_nei[0]*outward_normal_edge[0]+speed_nei[1]*outward_normal_edge[1]); 


      // .5 is the base function evaluated in the middle, mid-point integration
      if (! quadrant->is_hanging (i_1))
      {
        for (int kk=ordBottom(quadrant->gt (i_1)); kk<ordSurface(quadrant->gt (i_1)); kk++) // excluded the top because it is not a dof for the DBC
        {
          int kkk = kk%number_FD_points;

          const auto dp_cell_i_1 = linear_interpolation(kkk, index_quadrant, h_cell, hdof[i_1], Uxwdof[i_1]+Uxsdof[i_1], Uywdof[i_1]+Uysdof[i_1]);
          const auto dp_cell_nei_i_1 = dp_cell_i_1;
          const auto flux_int_dp_kk  = .5*((vxs_cell*dp_cell_i_1+vxs_cell_nei*dp_cell_nei_i_1)*outward_normal_edge[0] + (vys_cell*dp_cell_i_1+vys_cell_nei*dp_cell_nei_i_1)*outward_normal_edge[1]) - .5*smax*(dp_cell_nei_i_1 - dp_cell_i_1 );
          
          excess_pore_water_pressure_incr[kk] += -edge_length*flux_int_dp_kk *.5;

          //std::cout << flux_int_dp_kk << std::endl;
        }
      }
      else
      {
        for (int jj=0; jj<=1; jj++)
        {
          for (int kk=ordBottom(quadrant->gparent(jj,i_1)); kk<ordSurface(quadrant->gparent(jj,i_1)); kk++) // excluded the top because it is not a dof for the DBC
          {
            int kkk = kk%number_FD_points;

            const auto dp_cell_i_1 = linear_interpolation(kkk, index_quadrant, h_cell, hdof[i_1], Uxwdof[i_1]+Uxsdof[i_1], Uywdof[i_1]+Uysdof[i_1]);
            const auto dp_cell_nei_i_1 = dp_cell_i_1; //linear_interpolation(kkk, index_quadrant_nei, h_cell_nei, hdof[i_1], Uxwdof[i_1]+Uxsdof[i_1], Uywdof[i_1]+Uysdof[i_1]);
            const auto flux_int_dp_kk  = .5*((vxs_cell*dp_cell_i_1+vxs_cell_nei*dp_cell_nei_i_1)*outward_normal_edge[0] + (vys_cell*dp_cell_i_1+vys_cell_nei*dp_cell_nei_i_1)*outward_normal_edge[1]) - .5*smax*(dp_cell_nei_i_1 - dp_cell_i_1 );
          
            excess_pore_water_pressure_incr[kk] += -edge_length*flux_int_dp_kk *.5*.5;
          }
        }
      }

      if (! quadrant->is_hanging (i_2))
      {
        for (int kk=ordBottom(quadrant->gt (i_2)); kk<ordSurface(quadrant->gt (i_2)); kk++) // excluded the top because it is not a dof for the DBC
        {
          int kkk = kk%number_FD_points;

          const auto dp_cell_i_2 = linear_interpolation(kkk, index_quadrant, h_cell, hdof[i_2], Uxwdof[i_2]+Uxsdof[i_2], Uywdof[i_2]+Uysdof[i_2]);
          const auto dp_cell_nei_i_2 = dp_cell_i_2; //linear_interpolation(kkk, index_quadrant_nei, h_cell_nei, hdof[i_2], Uxwdof[i_2]+Uxsdof[i_2], Uywdof[i_2]+Uysdof[i_2]);
          const auto flux_int_dp_kk  = .5*((vxs_cell*dp_cell_i_2+vxs_cell_nei*dp_cell_nei_i_2)*outward_normal_edge[0] + (vys_cell*dp_cell_i_2+vys_cell_nei*dp_cell_nei_i_2)*outward_normal_edge[1]) - .5*smax*(dp_cell_nei_i_2 - dp_cell_i_2 );
          
          excess_pore_water_pressure_incr[kk] += -edge_length*flux_int_dp_kk *.5;
        }
      }
      else
      {
        for (int jj=0; jj<=1; jj++)
        {
          for (int kk=ordBottom(quadrant->gparent(jj,i_2)); kk<ordSurface(quadrant->gparent(jj,i_2)); kk++) // excluded the top because it is not a dof for the DBC
          {
            int kkk = kk%number_FD_points;

            const auto dp_cell_i_2 = linear_interpolation(kkk, index_quadrant, h_cell, hdof[i_2], Uxwdof[i_2]+Uxsdof[i_2], Uywdof[i_2]+Uysdof[i_2]);
            const auto dp_cell_nei_i_2 = dp_cell_i_2; //linear_interpolation(kkk, index_quadrant_nei, h_cell_nei, hdof[i_2], Uxwdof[i_2]+Uxsdof[i_2], Uywdof[i_2]+Uysdof[i_2]);
            const auto flux_int_dp_kk  = .5*((vxs_cell*dp_cell_i_2+vxs_cell_nei*dp_cell_nei_i_2)*outward_normal_edge[0] + (vys_cell*dp_cell_i_2+vys_cell_nei*dp_cell_nei_i_2)*outward_normal_edge[1]) - .5*smax*(dp_cell_nei_i_2 - dp_cell_i_2 );
          
            excess_pore_water_pressure_incr[kk] += -edge_length*flux_int_dp_kk *.5*.5;
          }
        }
      }
    }

  }

  const double Ux_cell = Uxw_cell+Uxs_cell;
  const double Uy_cell = Uyw_cell+Uys_cell;
  const auto erosion_contribution_cell = erosion_coefficient*std::sqrt(Ux_cell*Ux_cell + Uy_cell*Uy_cell);

  for (int ii=0; ii<4; ii++)
  {
    solid_vel_x[ii] = hsdof[ii]>epsilon ? Uxsdof[ii]/hsdof[ii] : 0.;
    solid_vel_y[ii] = hsdof[ii]>epsilon ? Uysdof[ii]/hsdof[ii] : 0.;
  }
  const double vs_med_x = std::abs(solid_vel_x[0] + solid_vel_x[1] + solid_vel_x[2] + solid_vel_x[3])*.25;
  const double vs_med_y = std::abs(solid_vel_y[0] + solid_vel_y[1] + solid_vel_y[2] + solid_vel_y[3])*.25;
  const double hs_med   = (hsdofold[0] + hsdofold[1] + hsdofold[2] + hsdofold[3])*.25;


  grad_cell_hw    = {.5 * ( (hwdofold[3]-hwdofold[2]) + (hwdofold[1]-hwdofold[0]) ), .5 * ( (hwdofold[2]-hwdofold[0]) + (hwdofold[3]-hwdofold[1]) )};
  grad_cell_Uxw   = {.5 * ( (Uxwdof[3]-Uxwdof[2]) + (Uxwdof[1]-Uxwdof[0]) ), .5 * ( (Uxwdof[2]-Uxwdof[0]) + (Uxwdof[3]-Uxwdof[1]) )};
  grad_cell_Uyw   = {.5 * ( (Uywdof[3]-Uywdof[2]) + (Uywdof[1]-Uywdof[0]) ), .5 * ( (Uywdof[2]-Uywdof[0]) + (Uywdof[3]-Uywdof[1]) )};

  grad_cell_vsx = {.5 * ( (solid_vel_x[3]-solid_vel_x[2]) + (solid_vel_x[1]-solid_vel_x[0]) ), .5 * ( (solid_vel_x[2]-solid_vel_x[0]) + (solid_vel_x[3]-solid_vel_x[1]) )};
  grad_cell_vsy = {.5 * ( (solid_vel_y[3]-solid_vel_y[2]) + (solid_vel_y[1]-solid_vel_y[0]) ), .5 * ( (solid_vel_y[2]-solid_vel_y[0]) + (solid_vel_y[3]-solid_vel_y[1]) )};

  //const auto diff_term_Uys_x = grad_cell_Uys [0]*vel_rusanov_cell_y*.5;
  //const auto diff_term_Uys_y = grad_cell_Uys [1]*vel_rusanov_cell_x*.5;




  for (int ii = 0; ii < 4; ++ii){

    double dp_;
    int kkk;

    //consolidation_coefficient = hdof[ii]>epsilon ? odometric_coeff/((hwdof[ii]>epsilon && hdof[ii]>epsilon && hsdof[ii]>epsilon) ? hsdof[ii]/hdof[ii]/std::pow(hwdof[ii]/hdof[ii], m_coeff)/terminal_velocity*(density_s-density_w)*grav : 0.) : 1.e-10;

    if (! quadrant->is_hanging (ii)){
      for (int kk=ordBottom(quadrant->gt (ii)) + (Z [quadrant->gt (ii)]<thr_erodible_layer ? 0 : 1); kk<ordSurface(quadrant->gt (ii)); kk++) // excluded the top because it is not a dof for the DBC
      {
        kkk = kk%number_FD_points;

        const auto dp_cell_ii = linear_interpolation(kkk, index_quadrant, h_cell, hdof[ii], Uxwdof[ii]+Uxsdof[ii], Uywdof[ii]+Uysdof[ii]);

        const double delta_h_c = hdofold[ii]/(number_FD_points-1);
        const double z = kkk*delta_h_c;
        const double zz = kkk*delta_h_cell;
        const double vel_z_current_node = h_cell>epsilon ? erosion_contribution_cell*(1. - zz/h_cell) : 0.;
        const double Z_1 = z - tau*vel_z_current_node;
        const double Z_2 = Z_1 - dt*vel_z_current_node;
        
        const auto Pxi = linear_interpolation(ordBottom(quadrant->gt (ii)), delta_h_c, Z_2);

        pressure = {0,0,0,0};
        for (int iii=0; iii<4; iii++)
        {
          const double delta_h_old_c = hdofold[iii]/(number_FD_points-1);
          const double Z_moved_c = hdofold[ii]>epsilon ? Z_2*(hdofold[iii]/hdofold[ii]) : 0.;

          if (! quadrant->is_hanging (iii))
          {
            pressure[iii] = linear_interpolation(ordBottom(quadrant->gt (iii)), delta_h_old_c, Z_moved_c);
          }
          else
          {
            for (int jjj=0; jjj<=1; jjj++)
            {
              pressure[iii] += .5*linear_interpolation(ordBottom(quadrant->gparent (jjj,iii)), delta_h_old_c, Z_moved_c);
            }
          }
          
        }

        dp_old = linear_interpolation_second(Z [quadrant->gt (ii)], delta_h_c, Pxi, Z_2, ordBottom(quadrant->gt (ii)), hdofold[ii], vxs_cell, vys_cell);

        grad_dp_old = {.5 * ( (pressure[3]-pressure[2]) + (pressure[1]-pressure[0]) ), .5 * ( (pressure[2]-pressure[0]) + (pressure[3]-pressure[1]) )};

        const auto dp_diffx_kk = grad_dp_old[0]*vs_med_x*.5;
        const auto dp_diffy_kk = grad_dp_old[1]*vs_med_y*.5;
        //std::cout << dp_cell_ii << " " << excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)+kkk] << std::endl;

        const auto dp_al_kk = der_coeffs_x[ii]*(vxs_cell*dp_cell_ii - dp_diffx_kk) + der_coeffs_y[ii]*(vys_cell*dp_cell_ii - dp_diffy_kk);
        const auto dp_kk = der_coeffs_x[ii]*dp_diffx_kk + der_coeffs_y[ii]*dp_diffy_kk;

        incr_anti_diff_pressure[ordBottom(index_quadrant_local) + kkk][ii] = dp_kk;

        P_plus_pressure  [kk] += std::max(0., dp_kk );
        P_minus_pressure [kk] += std::min(0., dp_kk );

        //if (grad_dp_old[0]!=0)
        //{
        //  std::cout << grad_dp_old[0] << " " << diff_term_x << " " << grad_dp_old[1] << std::endl;
        //}
        //if (grad_cell_hw[0]!=0) //(grad_dp_old[0]*vs_med_x!=0)
        //std::cout << vs_med_x << " " << grad_dp_old[0] << " " << (grad_cell_Uxw[0] - grad_cell_hw[0]*vs_med_x)*C_1(z) + hs_med*grad_cell_vsx[0]*C_2(z) << std::endl;

        dp_ = -.5*(Dy*contr_x[kkk][ii] + Dx*contr_y[kkk][ii]) + dp_al_kk + .25*area*(Pxi/dt + C(zz)*erosion_contribution_cell + dp_old); 
        excess_pore_water_pressure_incr[kk] += is_pore_water_pressure*dp_; //.25*area*(Pxi/dt + dp_old) + excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)+kkk]; //is_pore_water_pressure*dp_;
      } 
    } 
    else {
      for (int jj=0; jj<=1; jj++)
      {
        for (int kk=ordBottom(quadrant->gparent (jj,ii)) + (Z [quadrant->gparent(jj,ii)]<thr_erodible_layer ? 0 : 1); kk<ordSurface(quadrant->gparent (jj,ii)); kk++) // excluded the top because it is not a dof for the DBC
        {

          kkk = kk%number_FD_points;

          const auto dp_cell_ii = linear_interpolation(kkk, index_quadrant, h_cell, hdof[ii], Uxwdof[ii]+Uxsdof[ii], Uywdof[ii]+Uysdof[ii]);

          const double delta_h_c = hdofold[ii]/(number_FD_points-1);
          const double z = kkk*delta_h_c;
          const double zz = kkk*delta_h_cell;
          const double vel_z_current_node = h_cell>epsilon ? erosion_contribution_cell*(1. - zz/h_cell) : 0.;
          const double Z_1 = z - tau*vel_z_current_node;
          const double Z_2 = Z_1 - dt*vel_z_current_node;
          const auto Pxi = linear_interpolation(ordBottom(quadrant->gt (ii)), delta_h_c, Z_2);


          pressure = {0,0,0,0};
          for (int iii=0; iii<4; iii++)
          {
            const double delta_h_old_c = hdofold[iii]/(number_FD_points-1);
            const double Z_moved_c = hdofold[ii]>epsilon ? Z_2*(hdofold[iii]/hdofold[ii]) : 0.;
          
            if (! quadrant->is_hanging (iii))
            {
              pressure[iii] = linear_interpolation(ordBottom(quadrant->gt (iii)), delta_h_old_c, Z_moved_c);
            }
            else
            {
              for (int jjj=0; jjj<=1; jjj++)
              {
                pressure[iii] += .5*linear_interpolation(ordBottom(quadrant->gparent (jjj,iii)), delta_h_old_c, Z_moved_c);
              }
            }
          }

          dp_old = linear_interpolation_second(Z [quadrant->gparent(jj,ii)], delta_h_c, Pxi, Z_2, ordBottom(quadrant->gparent (jj,ii)), sold [ordhw  (quadrant->gparent (jj, ii) )] + sold [ordhs  (quadrant->gparent (jj, ii) )], vxs_cell, vys_cell);

          grad_dp_old = {.5 * ( (pressure[3]-pressure[2]) + (pressure[1]-pressure[0]) ), .5 * ( (pressure[2]-pressure[0]) + (pressure[3]-pressure[1]) )};

          const auto dp_diffx_kk = grad_dp_old[0]*vs_med_x*.5;
          const auto dp_diffy_kk = grad_dp_old[1]*vs_med_y*.5;

          const auto dp_al_kk = der_coeffs_x[ii]*(vxs_cell*dp_cell_ii - dp_diffx_kk) + der_coeffs_y[ii]*(vys_cell*dp_cell_ii - dp_diffy_kk);
          const auto dp_kk = der_coeffs_x[ii]*dp_diffx_kk + der_coeffs_y[ii]*dp_diffy_kk;

          incr_anti_diff_pressure[ordBottom(index_quadrant_local) + kkk][ii] = dp_kk;

          P_plus_pressure  [kk] += std::max(0., dp_kk );
          P_minus_pressure [kk] += std::min(0., dp_kk );

          dp_ = -.5*(Dy*contr_x[kkk][ii] + Dx*contr_y[kkk][ii]) + dp_al_kk + .25*area*(Pxi/dt + C(zz)*erosion_contribution_cell + dp_old)*.5;
          excess_pore_water_pressure_incr[kk] += is_pore_water_pressure*dp_;
        }
      }
      
    }
    

  }

}




double
TG2_scheme::linear_interpolation_second(const double& Z_node_c, const double& delta_h_old, const double& P_xi, const double& Z_moved, const int& kkk_ini, const double& hdofold_c, const double& vxs_cell, const double& vys_cell)
{

  double Z_plus  = Z_moved+delta_h_old;
  Z_plus = Z_plus>hdofold_c ? hdofold_c : Z_plus; // BC at the top
  const auto P_xiplus = linear_interpolation(kkk_ini, delta_h_old, Z_plus);

  double Z_minus = Z_moved-delta_h_old;
  //Z_minus = Z_minus<0 ? 0. : Z_minus; // Dirichlet BC at the bottom
  //auto P_ximinus = linear_interpolation(kkk_ini, delta_h_old, Z_minus);

  auto P_ximinus = Z_minus<0 ? (Z_node_c<thr_erodible_layer ? P_xiplus : linear_interpolation(kkk_ini, delta_h_old, 0.)) : linear_interpolation(kkk_ini, delta_h_old, Z_minus);

  const double cfl_number_mu = dt/std::max(Dx,Dy)*std::sqrt(vxs_cell*vxs_cell + vys_cell*vys_cell);

  const double second_der_dp = (delta_h_old>epsilon && dt<=.5*delta_h_old*delta_h_old/consolidation_coefficient*(1.- cfl_number_mu*cfl_number_mu )) ? consolidation_coefficient*(P_xiplus - 2.*P_xi + P_ximinus)/delta_h_old/delta_h_old : 0.;

  return(second_der_dp);
}


double
TG2_scheme::linear_interpolation(const int& kkk_ini, const double& delta_h_old, const double& Z_moved)
{
  // 1d linear interpolation, to change the space accuracy change this
  double i_z = delta_h_old>epsilon ? Z_moved/delta_h_old : 0.;
  i_z = i_z>(number_FD_points-1) ? (number_FD_points-1) : i_z;
  i_z = i_z<0 ? 0 : i_z; // here we are imposing BC

  const int iz_minus = std::floor(i_z);
  const int iz_plus  = std::ceil (i_z);

  const double delta_iz = (iz_plus-iz_minus);

  const double w1 = delta_iz>0 ? (iz_plus-i_z) /(iz_plus-iz_minus) : 0.;
  const double w2 = delta_iz>0 ? (i_z-iz_minus)/(iz_plus-iz_minus) : 1.;

  //std::cout << kkk_ini << " " << iz_minus << " " << iz_plus << " " << i_z << " " << Z_moved << std::endl;

  return(w1*excess_pore_water_pressure[kkk_ini+iz_minus] + w2*excess_pore_water_pressure[kkk_ini+iz_plus]);
}




double
TG2_scheme::linear_interpolation(const int& kk, const int& index_quadrant_c, const double& h_tot_old, const double& hdof_c, const double& Ux_tot_c, const double & Uy_tot_c)
{
  // 
  const double delta_h_node_c = hdof_c/(number_FD_points-1);
  const double z_current_mesh_node = kk*delta_h_node_c;
  const double constant_ratio = hdof_c>epsilon ? z_current_mesh_node/hdof_c : 0.;
  const auto erosion_contribution_c = erosion_coefficient*std::sqrt(Ux_tot_c*Ux_tot_c + Uy_tot_c*Uy_tot_c);
  const double vel_z_current_node = hdof_c>epsilon ? erosion_contribution_c*(1. - constant_ratio) : 0.;

  double z_moved = z_current_mesh_node - tau*vel_z_current_node;
  z_moved = z_moved<0 ? 0. : z_moved; // here we are imposing BC


  const double delta_h_old = h_tot_old/(number_FD_points-1);
  const double Z_moved = hdof_c>epsilon ? z_moved*(h_tot_old/hdof_c) : 0.;


  // 1d linear interpolation, to change the space accuracy change this, this part below is similar to the function linear_interpolation with three inputs defined above, the difference is just in the excess_pore_water_pressure_onehalf the other one has excess_pore_water_pressure
  double i_z = delta_h_old>epsilon ? Z_moved/delta_h_old : 0.;
  i_z = i_z>(number_FD_points-1) ? (number_FD_points-1) : i_z;

  const int iz_minus = std::floor(i_z);
  const int iz_plus  = std::ceil (i_z);

  const double delta_iz = (iz_plus-iz_minus);

  const double w1 = delta_iz>0 ? (iz_plus-i_z) /(iz_plus-iz_minus) : 0.;
  const double w2 = delta_iz>0 ? (i_z-iz_minus)/(iz_plus-iz_minus) : 1.;

  return(w1*excess_pore_water_pressure_onehalf[ordBottom(index_quadrant_c)+iz_minus] + w2*excess_pore_water_pressure_onehalf[ordBottom(index_quadrant_c)+iz_plus]);
}



void 
TG2_scheme::thomas_algorithm(const double& mu_coeff, const int& kk_, const int& kk_ini, const int& kk_fin, const int& kk_ini_bottom) 
{

  const double b_coeff = 1.+2.*mu_coeff;
  const double a_coeff = -mu_coeff;
  const double c_coeff = -mu_coeff;

  const auto kkk_ini = kk_ini%number_FD_points;
  const auto kkk_fin = kk_fin%number_FD_points;

  const double hdof_c  = sol.get_owned_data ()[kk_] + sol.get_owned_data ()[kk_+1];

  //std::cout << mass.get_owned_data().size() << " " << kk_ << " " << mass.get_owned_data() [kk_] << std::endl;

  alfa_vec[kkk_ini] = b_coeff;
  y_vec[kkk_ini] = dt*excess_pore_water_pressure_incr.get_owned_data() [kk_ini]/mass.get_owned_data() [kk_] + kkk_ini*mu_coeff*excess_pore_water_pressure.get_owned_data() [kk_ini_bottom];

  beta_vec[kkk_ini+1] = c_coeff/alfa_vec[kkk_ini];
  alfa_vec[kkk_ini+1] = b_coeff - beta_vec[kkk_ini+1]*c_coeff*(2-kkk_ini);
  y_vec[kkk_ini+1] = dt*excess_pore_water_pressure_incr.get_owned_data() [kk_ini+1]/mass.get_owned_data() [kk_] - beta_vec[kkk_ini+1]*y_vec[kkk_ini];

  for (int i=kk_ini+2; i<=kk_fin; i++) 
  {
    const auto i_ = i%number_FD_points;

    beta_vec[i_] = c_coeff/alfa_vec[i_-1];
    alfa_vec[i_] = b_coeff - beta_vec[i_]*c_coeff;

    y_vec[i_] = dt*excess_pore_water_pressure_incr.get_owned_data() [i]/mass.get_owned_data() [kk_] - beta_vec[i_]*y_vec[i_-1];
  }

  // This is the reverse sweep, 
  excess_pore_water_pressure.get_owned_data() [kk_fin] = y_vec[kkk_fin]/alfa_vec[kkk_fin];                                                                                                                                       
  for (int i=kk_fin-1; i>kk_ini; i--) 
  {
    const auto i_ = i%number_FD_points;

    excess_pore_water_pressure.get_owned_data() [i] = hdof_c>epsilon ? (y_vec[i_] - c_coeff*excess_pore_water_pressure.get_owned_data() [i+1])/alfa_vec[i_] : 0.;
    //if (excess_pore_water_pressure.get_owned_data() [i]!=0)
    //std::cout << excess_pore_water_pressure.get_owned_data() [i] << std::endl;
  }
  excess_pore_water_pressure.get_owned_data() [kk_ini] = hdof_c>epsilon ? (y_vec[kkk_ini] - (2-kkk_ini)*c_coeff*excess_pore_water_pressure.get_owned_data() [kk_ini+1])/alfa_vec[kkk_ini] : 0.;

/*
#include <iostream>
#include <vector>

using namespace std;

int main()
{
  std::vector<double> result(10), beta_vec(10), alfa_vec(10), y_vec(10);
  
  const double b_coeff = 2;
  const double a_coeff = -1;
  const double c_coeff = -1;

  alfa_vec[0] = b_coeff;
  y_vec[0] = 1.;

  beta_vec[1] = c_coeff/b_coeff;
  alfa_vec[1] = b_coeff - beta_vec[1]*c_coeff*2; // - beta_vec[kkk_ini+1]*c_coeff*2*(1-kkk_ini)
  y_vec[1] = 1. - beta_vec[1]*y_vec[0];

  for (int i=2; i<=9; i++) 
  {
    beta_vec[i] = c_coeff/alfa_vec[i-1];
    alfa_vec[i] = b_coeff - beta_vec[i]*c_coeff;

    y_vec[i] = 1. - beta_vec[i]*y_vec[i-1];
  }

  // This is the reverse sweep, 
  result [9] = y_vec[9]/alfa_vec[9];                                                                                                                                       
  for (int i=8; i>0; i--) 
  {
    result [i] = (y_vec[i] - c_coeff*result [i+1])/alfa_vec[i];
  }
  //result [0] = (y_vec[0] - c_coeff*result [1])/alfa_vec[0];
  result [0] = (y_vec[0] - 2*c_coeff*result [1])/alfa_vec[0];
  
  for (const auto & it : result)
  {
      std::cout << it << std::endl;
  }


    return 0;
}
*/


}




void
TG2_scheme::compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant)
{

  // look at tmesh.h 
  index_quadrant = quadrant->get_global_quad_idx (); 
  index_quadrant_local = quadrant->get_forest_quad_idx ();

  std::array<int,4> bimpp_to_rev_ord = {0, 1, 3, 2};
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;

  dp_mean_vec = {0., 0., 0., 0.};
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hwdof [ii] = sol [ordhw  (quadrant->gt (ii) )];
      hsdof [ii] = sol [ordhs  (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      Z_node[ii]  = Z [quadrant->gt (ii)];
      
      isdof_or_hanging[ii] = 1.;

      numerical_integration_pressure_2(ordBottom(quadrant->gt (ii)), dp_mean_vec[ii], excess_pore_water_pressure, 1.);
    }
    else
    {
      hwdof [ii] = .5 * (sol [ordhw  (quadrant->gparent (0, ii) )] +
                         sol [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof [ii] = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
      Uxwdof[ii] = .5 * (sol [ordUxw (quadrant->gparent (0, ii) )] +
                         sol [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof[ii] = .5 * (sol [ordUyw (quadrant->gparent (0, ii) )] +
                         sol [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);

      Z_node[ii]  = .5 * (Z [quadrant->gparent(0,ii)] +
                          Z [quadrant->gparent(1,ii)]);

      isdof_or_hanging[ii] = .5;

      for (int jj=0; jj<=1; jj++)
        numerical_integration_pressure_2(ordBottom(quadrant->gparent(jj,ii)), dp_mean_vec[ii], excess_pore_water_pressure, .5);

    }

    hdof   [ii] = hwdof[ii] + hsdof[ii];
  }

  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
    -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
   
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
    +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};

  double vel_rusanov_cell_x = 0., vel_rusanov_cell_y = 0.;
  for (int ii = 0; ii < 4; ++ii){
    const auto lambdas = max_eigen (hwdof[ii], hsdof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii], dp_mean_vec[ii]);
    vel_rusanov_cell_x += lambdas[0];
    vel_rusanov_cell_y += lambdas[1];
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.;


  grad_cell_Z     = {.5 * ( (Z_node  [3] - Z_node  [2]) + (Z_node  [1] - Z_node  [0]) ), .5 * ( (Z_node  [2] - Z_node  [0]) + (Z_node  [3] - Z_node  [1]) )};
 
  grad_cell_hw    = {.5 * ( (hwdof   [3] - hwdof   [2]) + (hwdof   [1] - hwdof   [0]) ), .5 * ( (hwdof   [2] - hwdof   [0]) + (hwdof   [3] - hwdof   [1]) )};
  grad_cell_Uxw   = {.5 * ( (Uxwdof  [3] - Uxwdof  [2]) + (Uxwdof  [1] - Uxwdof  [0]) ), .5 * ( (Uxwdof  [2] - Uxwdof  [0]) + (Uxwdof  [3] - Uxwdof  [1]) )};
  grad_cell_Uyw   = {.5 * ( (Uywdof  [3] - Uywdof  [2]) + (Uywdof  [1] - Uywdof  [0]) ), .5 * ( (Uywdof  [2] - Uywdof  [0]) + (Uywdof  [3] - Uywdof  [1]) )};

  grad_cell_hs    = {.5 * ( (hsdof   [3] - hsdof   [2]) + (hsdof   [1] - hsdof   [0]) ), .5 * ( (hsdof   [2] - hsdof   [0]) + (hsdof   [3] - hsdof   [1]) )};
  grad_cell_Uxs   = {.5 * ( (Uxsdof  [3] - Uxsdof  [2]) + (Uxsdof  [1] - Uxsdof  [0]) ), .5 * ( (Uxsdof  [2] - Uxsdof  [0]) + (Uxsdof  [3] - Uxsdof  [1]) )};
  grad_cell_Uys   = {.5 * ( (Uysdof  [3] - Uysdof  [2]) + (Uysdof  [1] - Uysdof  [0]) ), .5 * ( (Uysdof  [2] - Uysdof  [0]) + (Uysdof  [3] - Uysdof  [1]) )};


  const double & hw_cell    = sol_onehalf[ordhw    (index_quadrant)];
  const double & Uxw_cell   = sol_onehalf[ordUxw   (index_quadrant)];
  const double & Uyw_cell   = sol_onehalf[ordUyw   (index_quadrant)];

  const double & hs_cell    = sol_onehalf[ordhs    (index_quadrant)];
  const double & Uxs_cell   = sol_onehalf[ordUxs   (index_quadrant)];
  const double & Uys_cell   = sol_onehalf[ordUys   (index_quadrant)];

  const double & Z_cell     = Z_onehalf[index_quadrant];

  double dp_mean_cell = 0.;
  numerical_integration_pressure_2(ordBottom(index_quadrant), dp_mean_cell, excess_pore_water_pressure_onehalf, 1.);

  const auto slope_x_cell = 0.;//grad_cell_Z[0];
  const auto slope_y_cell = 0.;//grad_cell_Z[1];

  const auto contr_slope_xs = .5*src_slope_formula (hs_cell, slope_x_cell)*Dy;
  const auto contr_slope_ys = .5*src_slope_formula (hs_cell, slope_y_cell)*Dx;

  const auto contr_slope_xw = .5*src_slope_formula (hw_cell, slope_x_cell)*Dy;
  const auto contr_slope_yw = .5*src_slope_formula (hw_cell, slope_y_cell)*Dx;


  const auto hw_cell_  = .25*(hwdof [0]+hwdof [1]+hwdof [2]+hwdof [3]);
  const auto hs_cell_  = .25*(hsdof [0]+hsdof [1]+hsdof [2]+hsdof [3]);
  //const auto Uxw_cell_ = .25*(Uxwdof[0]+Uxwdof[1]+Uxwdof[2]+Uxwdof[3]);
  //const auto Uxs_cell_ = .25*(Uxsdof[0]+Uxsdof[1]+Uxsdof[2]+Uxsdof[3]);
  //const auto Uyw_cell_ = .25*(Uywdof[0]+Uywdof[1]+Uywdof[2]+Uywdof[3]);
  //const auto Uys_cell_ = .25*(Uysdof[0]+Uysdof[1]+Uysdof[2]+Uysdof[3]);
  

  //double dp_mean_cell_ = .25*(dp_mean_vec[0]+dp_mean_vec[1]+dp_mean_vec[2]+dp_mean_vec[3]);


  const auto h_cell_  = hw_cell_+hs_cell_;
  const auto n_cell_  = h_cell_>epsilon ? hw_cell_/h_cell_ : 0.; 
  const auto ns_cell_ = h_cell_>epsilon ? hs_cell_/h_cell_ : 0.; 


  const auto h_cell  = hw_cell+hs_cell;
  const auto n_cell  = h_cell>epsilon ? hw_cell/h_cell : 0.; 
  const auto ns_cell = h_cell>epsilon ? hs_cell/h_cell : 0.; 


  std::array<double,4> contr_x_w  = {0., 0., 0., 0.}, contr_y_w  = {0., 0., 0., 0.}, contr_x_s  = {0., 0., 0., 0.}, contr_y_s  = {0., 0., 0., 0.};
  //std::array<double,4> contr_x_w_ = {0., 0., 0., 0.}, contr_y_w_ = {0., 0., 0., 0.}, contr_x_s_ = {0., 0., 0., 0.}, contr_y_s_ = {0., 0., 0., 0.};


  // boundary conditions, just for the transport term!
  bool is_boundary_edge = true;
  for (int iEdge = 0; iEdge < 4; ++iEdge){

    is_boundary_edge = true;

    const auto i_1 = bimpp_to_rev_ord[iEdge];
    const auto i_2 = bimpp_to_rev_ord[(iEdge+1)%4];

    const auto edge_length = std::sqrt(std::pow((xn[i_1]-xn[i_2]),2.) + std::pow((yn[i_1]-yn[i_2]),2.));
    const std::array<double,2> outward_normal_edge = {(-yn[i_1]+yn[i_2])/edge_length, ( xn[i_1]-xn[i_2])/edge_length};  

    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {
      std::array<double,4> Xn, Yn;

      for (int ii = 0; ii < 4; ++ii) {
        Xn[ii] = quadrant_nei->p(0, ii);
        Yn[ii] = quadrant_nei->p(1, ii);

        if (! quadrant_nei->is_hanging (ii)){
          hwdof_nei [ii]  = sol [ordhw   (quadrant_nei->gt (ii) )];
          hsdof_nei [ii]  = sol [ordhs   (quadrant_nei->gt (ii) )];
          Uxsdof_nei[ii]  = sol [ordUxs  (quadrant_nei->gt (ii) )];
          Uysdof_nei[ii]  = sol [ordUys  (quadrant_nei->gt (ii) )];
          Uxwdof_nei[ii]  = sol [ordUxw  (quadrant_nei->gt (ii) )];
          Uywdof_nei[ii]  = sol [ordUyw  (quadrant_nei->gt (ii) )];
        } else {
          hwdof_nei[ii]   = .5 * (sol [ordhw   (quadrant_nei->gparent (0, ii) )] +
                                  sol [ordhw   (quadrant_nei->gparent (1, ii) )]);
          hsdof_nei[ii]   = .5 * (sol [ordhs   (quadrant_nei->gparent (0, ii) )] +
                                  sol [ordhs   (quadrant_nei->gparent (1, ii) )]);
          Uxsdof_nei[ii]  = .5 * (sol [ordUxs  (quadrant_nei->gparent (0, ii) )] +
                                  sol [ordUxs  (quadrant_nei->gparent (1, ii) )]);
          Uysdof_nei[ii]  = .5 * (sol [ordUys  (quadrant_nei->gparent (0, ii) )] +
                                  sol [ordUys  (quadrant_nei->gparent (1, ii) )]);
          Uxwdof_nei[ii]  = .5 * (sol [ordUxw  (quadrant_nei->gparent (0, ii) )] +
                                  sol [ordUxw  (quadrant_nei->gparent (1, ii) )]);
          Uywdof_nei[ii]  = .5 * (sol [ordUyw  (quadrant_nei->gparent (0, ii) )] +
                                  sol [ordUyw  (quadrant_nei->gparent (1, ii) )]);
        }

      }

      const auto & index_quadrant_nei = quadrant_nei->get_global_quad_idx ();


      for (int jEdge = 0; jEdge < 4; ++jEdge) { // cycle neigh edges 

        const auto j_1 = bimpp_to_rev_ord[jEdge];
        const auto j_2 = bimpp_to_rev_ord[(jEdge+1)%4];

        const auto edge_length_nei = std::sqrt(std::pow((Xn[j_1]-Xn[j_2]),2.) + std::pow((Yn[j_1]-Yn[j_2]),2.));
        const std::array<double,2> outward_normal_edge_nei = {(-Yn[j_1]+Yn[j_2])/edge_length_nei, ( Xn[j_1]-Xn[j_2])/edge_length_nei};  
        const bool check_orthogonality = std::inner_product(outward_normal_edge_nei.begin(), outward_normal_edge_nei.end(), outward_normal_edge.begin(), 0.) == -1;


        if ( (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) ||
               (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
              ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) ||
               (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant!=index_quadrant_nei )
        {
          is_boundary_edge = false;

          // scrivere qui la somma dei contributi per i termini non-cons.!
          double Z_cell_nei, hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei, hw_cell_nei_, hs_cell_nei_, Uxw_cell_nei_, Uyw_cell_nei_, Uxs_cell_nei_, Uys_cell_nei_;

          Z_cell_nei     = Z_onehalf  [          index_quadrant_nei ];
          hw_cell_nei    = sol_onehalf[ordhw    (index_quadrant_nei)];
          hs_cell_nei    = sol_onehalf[ordhs    (index_quadrant_nei)];
          Uxw_cell_nei   = sol_onehalf[ordUxw   (index_quadrant_nei)];
          Uyw_cell_nei   = sol_onehalf[ordUyw   (index_quadrant_nei)];
          Uxs_cell_nei   = sol_onehalf[ordUxs   (index_quadrant_nei)];
          Uys_cell_nei   = sol_onehalf[ordUys   (index_quadrant_nei)];


          hw_cell_nei_  = (hwdof_nei [0]+hwdof_nei [1]+hwdof_nei [2]+hwdof_nei [3])*.25;
          hs_cell_nei_  = (hsdof_nei [0]+hsdof_nei [1]+hsdof_nei [2]+hsdof_nei [3])*.25;
          Uxs_cell_nei_ = (Uxsdof_nei[0]+Uxsdof_nei[1]+Uxsdof_nei[2]+Uxsdof_nei[3])*.25;
          Uys_cell_nei_ = (Uysdof_nei[0]+Uysdof_nei[1]+Uysdof_nei[2]+Uysdof_nei[3])*.25;
          Uxw_cell_nei_ = (Uxwdof_nei[0]+Uxwdof_nei[1]+Uxwdof_nei[2]+Uxwdof_nei[3])*.25;
          Uyw_cell_nei_ = (Uywdof_nei[0]+Uywdof_nei[1]+Uywdof_nei[2]+Uywdof_nei[3])*.25;


          //std::cout << (Z_cell_nei - Z_cell) << std::endl;

          // .5 salta fuori dall'integrazione per trapezi tra 0 e 1 in coordinata \xi (è il valore in LHS da metter qui sotto!) 
          contr_x_w[i_1] += hw_cell>epsilon ? .5*signum(outward_normal_edge[0])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         (grav*hw_cell+n_cell*dp_mean_cell/density_w)*(hs_cell_nei - hs_cell) + n_cell*dp_mean_cell/density_w*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_1] : 0.;
          contr_x_w[i_2] += hw_cell>epsilon ? .5*signum(outward_normal_edge[0])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         (grav*hw_cell+n_cell*dp_mean_cell/density_w)*(hs_cell_nei - hs_cell) + n_cell*dp_mean_cell/density_w*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_2] : 0.;

          contr_y_w[i_1] += hw_cell>epsilon ? .5*signum(outward_normal_edge[1])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         (grav*hw_cell+n_cell*dp_mean_cell/density_w)*(hs_cell_nei - hs_cell) + n_cell*dp_mean_cell/density_w*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_1] : 0.;
          contr_y_w[i_2] += hw_cell>epsilon ? .5*signum(outward_normal_edge[1])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         (grav*hw_cell+n_cell*dp_mean_cell/density_w)*(hs_cell_nei - hs_cell) + n_cell*dp_mean_cell/density_w*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_2] : 0.;


          contr_x_s[i_1] += hs_cell>epsilon ? .5*signum(outward_normal_edge[0])*(grav*hs_cell*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell-n_cell*dp_mean_cell/density_s)*(hw_cell_nei - hw_cell) - n_cell*dp_mean_cell/density_s*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_1] : 0.;
          contr_x_s[i_2] += hs_cell>epsilon ? .5*signum(outward_normal_edge[0])*(grav*hs_cell*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell-n_cell*dp_mean_cell/density_s)*(hw_cell_nei - hw_cell) - n_cell*dp_mean_cell/density_s*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_2] : 0.;

          contr_y_s[i_1] += hs_cell>epsilon ? .5*signum(outward_normal_edge[1])*(grav*hs_cell*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell-n_cell*dp_mean_cell/density_s)*(hw_cell_nei - hw_cell) - n_cell*dp_mean_cell/density_s*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_1] : 0.;
          contr_y_s[i_2] += hs_cell>epsilon ? .5*signum(outward_normal_edge[1])*(grav*hs_cell*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell-n_cell*dp_mean_cell/density_s)*(hw_cell_nei - hw_cell) - n_cell*dp_mean_cell/density_s*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_2] : 0.;

/*
          // .5 salta fuori dall'integrazione per trapezi tra 0 e 1 in coordinata \xi (è il valore in LHS da metter qui sotto!) 
          contr_x_w_[i_1] += .5*signum(outward_normal_edge[0])*(grav*hw_cell_*(Z_cell_nei - Z_cell) +         (grav*hw_cell_+n_cell_*dp_mean_cell_/density_w)*(hs_cell_nei_ - hs_cell_) + n_cell_*dp_mean_cell_/density_w*(hw_cell_nei_ - hw_cell_))*isdof_or_hanging[i_1];
          contr_x_w_[i_2] += .5*signum(outward_normal_edge[0])*(grav*hw_cell_*(Z_cell_nei - Z_cell) +         (grav*hw_cell_+n_cell_*dp_mean_cell_/density_w)*(hs_cell_nei_ - hs_cell_) + n_cell_*dp_mean_cell_/density_w*(hw_cell_nei_ - hw_cell_))*isdof_or_hanging[i_2];

          contr_y_w_[i_1] += .5*signum(outward_normal_edge[1])*(grav*hw_cell_*(Z_cell_nei - Z_cell) +         (grav*hw_cell_+n_cell_*dp_mean_cell_/density_w)*(hs_cell_nei_ - hs_cell_) + n_cell_*dp_mean_cell_/density_w*(hw_cell_nei_ - hw_cell_))*isdof_or_hanging[i_1];
          contr_y_w_[i_2] += .5*signum(outward_normal_edge[1])*(grav*hw_cell_*(Z_cell_nei - Z_cell) +         (grav*hw_cell_+n_cell_*dp_mean_cell_/density_w)*(hs_cell_nei_ - hs_cell_) + n_cell_*dp_mean_cell_/density_w*(hw_cell_nei_ - hw_cell_))*isdof_or_hanging[i_2];


          contr_x_s_[i_1] += .5*signum(outward_normal_edge[0])*(grav*hs_cell_*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell_-n_cell_*dp_mean_cell_/density_s)*(hw_cell_nei_ - hw_cell_) - n_cell_*dp_mean_cell_/density_s*(hs_cell_nei_ - hs_cell_))*isdof_or_hanging[i_1];
          contr_x_s_[i_2] += .5*signum(outward_normal_edge[0])*(grav*hs_cell_*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell_-n_cell_*dp_mean_cell_/density_s)*(hw_cell_nei_ - hw_cell_) - n_cell_*dp_mean_cell_/density_s*(hs_cell_nei_ - hs_cell_))*isdof_or_hanging[i_2];

          contr_y_s_[i_1] += .5*signum(outward_normal_edge[1])*(grav*hs_cell_*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell_-n_cell_*dp_mean_cell_/density_s)*(hw_cell_nei_ - hw_cell_) - n_cell_*dp_mean_cell_/density_s*(hs_cell_nei_ - hs_cell_))*isdof_or_hanging[i_1];
          contr_y_s_[i_2] += .5*signum(outward_normal_edge[1])*(grav*hs_cell_*(Z_cell_nei - Z_cell) + (r_coeff*grav*hs_cell_-n_cell_*dp_mean_cell_/density_s)*(hw_cell_nei_ - hw_cell_) - n_cell_*dp_mean_cell_/density_s*(hs_cell_nei_ - hs_cell_))*isdof_or_hanging[i_2];
*/

          //break; // this just goes outside the jEdge cycle 
        }
      }

    }

    if (is_boundary_edge) // set boundary conditions
    { 

      auto hw_cell_nei  = hw_cell;
      auto Uxw_cell_nei = Uxw_cell;
      auto Uyw_cell_nei = Uyw_cell;

      auto hs_cell_nei  = hs_cell;
      auto Uxs_cell_nei = Uxs_cell;
      auto Uys_cell_nei = Uys_cell;

      double dp_mean_cell_c = 0; //dp_mean_cell;
      auto dp_mean_cell_nei_c = dp_mean_cell_c;

      Uxw_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxw_cell + outward_normal_edge[1]*Uyw_cell)*outward_normal_edge[0];
      Uyw_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxw_cell + outward_normal_edge[1]*Uyw_cell)*outward_normal_edge[1];

      Uxs_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[0];
      Uys_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[1];

      const auto speed     = max_eigen(hw_cell,     hs_cell,     Uxw_cell,     Uyw_cell,     Uxs_cell,     Uys_cell,     dp_mean_cell_c);
      const auto speed_nei = max_eigen(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei, dp_mean_cell_nei_c);

      const auto smax = std::max(speed[0]*outward_normal_edge[0]+speed[1]*outward_normal_edge[1], speed_nei[0]*outward_normal_edge[0]+speed_nei[1]*outward_normal_edge[1]); 

      const auto flux_int_hw  = .5*((hw_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+hw_flux_formula_x (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[0] + (hw_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+hw_flux_formula_y (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[1]) - .5*smax*(hw_cell_nei -hw_cell );
      const auto flux_int_Uxw = .5*((Uxw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, dp_mean_cell_c)+Uxw_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei, dp_mean_cell_nei_c))*outward_normal_edge[0] + (Uxw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+Uxw_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[1]) - .5*smax*(Uxw_cell_nei-Uxw_cell);
      const auto flux_int_Uyw = .5*((Uyw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+Uyw_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[0] + (Uyw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, dp_mean_cell_c)+Uyw_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei, dp_mean_cell_nei_c))*outward_normal_edge[1]) - .5*smax*(Uyw_cell_nei-Uyw_cell);

      const auto flux_int_hs  = .5*((hs_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+hs_flux_formula_x (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[0] + (hs_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+hs_flux_formula_y (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[1]) - .5*smax*(hs_cell_nei -hs_cell );
      const auto flux_int_Uxs = .5*((Uxs_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, dp_mean_cell_c)+Uxs_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei, dp_mean_cell_nei_c))*outward_normal_edge[0] + (Uxs_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+Uxs_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[1]) - .5*smax*(Uxs_cell_nei-Uxs_cell);
      const auto flux_int_Uys = .5*((Uys_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell                )+Uys_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei                    ))*outward_normal_edge[0] + (Uys_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, dp_mean_cell_c)+Uys_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei, dp_mean_cell_nei_c))*outward_normal_edge[1]) - .5*smax*(Uys_cell_nei-Uys_cell);


      auto outward_normal_edge_abs = outward_normal_edge;
      outward_normal_edge_abs[0] = std::abs(outward_normal_edge[0]);
      outward_normal_edge_abs[1] = std::abs(outward_normal_edge[1]);


      // .5 is the base function evaluated in the middle, mid-point intergration
      if (! quadrant->is_hanging (i_1))
      {
        incr[ordhw   (quadrant->gt (i_1))] += -edge_length*flux_int_hw *.5;
        incr[ordUxw  (quadrant->gt (i_1))] += -edge_length*flux_int_Uxw*.5 + contr_slope_xw*outward_normal_edge_abs[0];
        incr[ordUyw  (quadrant->gt (i_1))] += -edge_length*flux_int_Uyw*.5 + contr_slope_yw*outward_normal_edge_abs[1];

        incr[ordhs   (quadrant->gt (i_1))] += -edge_length*flux_int_hs *.5;
        incr[ordUxs  (quadrant->gt (i_1))] += -edge_length*flux_int_Uxs*.5 + contr_slope_xs*outward_normal_edge_abs[0];
        incr[ordUys  (quadrant->gt (i_1))] += -edge_length*flux_int_Uys*.5 + contr_slope_ys*outward_normal_edge_abs[1];
      }
      else
      {
        incr [ordhw  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_hw *.5*.5;
        incr [ordhw  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_hw *.5*.5;
      
        incr [ordUxw (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge_abs[0];
        incr [ordUxw (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge_abs[0];
      
        incr [ordUyw (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge_abs[1];
        incr [ordUyw (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge_abs[1];


        incr [ordhs  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_hs *.5*.5;
        incr [ordhs  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_hs *.5*.5;
      
        incr [ordUxs (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge_abs[0];
        incr [ordUxs (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge_abs[0];
      
        incr [ordUys (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge_abs[1];
        incr [ordUys (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge_abs[1];
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr[ordhw   (quadrant->gt (i_2))] += -edge_length*flux_int_hw *.5;
        incr[ordUxw  (quadrant->gt (i_2))] += -edge_length*flux_int_Uxw*.5 + contr_slope_xw*outward_normal_edge_abs[0];
        incr[ordUyw  (quadrant->gt (i_2))] += -edge_length*flux_int_Uyw*.5 + contr_slope_yw*outward_normal_edge_abs[1];

        incr[ordhs   (quadrant->gt (i_2))] += -edge_length*flux_int_hs *.5;
        incr[ordUxs  (quadrant->gt (i_2))] += -edge_length*flux_int_Uxs*.5 + contr_slope_xs*outward_normal_edge_abs[0];
        incr[ordUys  (quadrant->gt (i_2))] += -edge_length*flux_int_Uys*.5 + contr_slope_ys*outward_normal_edge_abs[1];
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr [ordhw  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_hw *.5*.5;
        incr [ordhw  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_hw *.5*.5;
      
        incr [ordUxw (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge_abs[0];
        incr [ordUxw (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge_abs[0];
      
        incr [ordUyw (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge_abs[1];
        incr [ordUyw (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge_abs[1];


        incr [ordhs  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_hs *.5*.5;
        incr [ordhs  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_hs *.5*.5;
      
        incr [ordUxs (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge_abs[0];
        incr [ordUxs (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge_abs[0];
      
        incr [ordUys (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge_abs[1];
        incr [ordUys (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge_abs[1];
      }
    }

  }

  
  
  const auto diff_term_hw_x  = h_cell_>epsilon ? (grad_cell_Z[0]*n_cell_ *g_coeff+grad_cell_hw[0])*vel_rusanov_cell_y*.5 : grad_cell_hw[0]*vel_rusanov_cell_y*.5;
  const auto diff_term_hw_y  = h_cell_>epsilon ? (grad_cell_Z[1]*n_cell_ *g_coeff+grad_cell_hw[1])*vel_rusanov_cell_x*.5 : grad_cell_hw[1]*vel_rusanov_cell_x*.5;

  const auto diff_term_hs_x  = h_cell_>epsilon ? (grad_cell_Z[0]*ns_cell_*g_coeff+grad_cell_hs[0])*vel_rusanov_cell_y*.5 : grad_cell_hs[0]*vel_rusanov_cell_y*.5;
  const auto diff_term_hs_y  = h_cell_>epsilon ? (grad_cell_Z[1]*ns_cell_*g_coeff+grad_cell_hs[1])*vel_rusanov_cell_x*.5 : grad_cell_hs[1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uxw_x = grad_cell_Uxw [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uxw_y = grad_cell_Uxw [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uyw_x = grad_cell_Uyw [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uyw_y = grad_cell_Uyw [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uxs_x = grad_cell_Uxs [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uxs_y = grad_cell_Uxs [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uys_x = grad_cell_Uys [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uys_y = grad_cell_Uys [1]*vel_rusanov_cell_x*.5;

/*
  //
  const auto F_star_hw_x_  = hw_flux_formula_x (hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_hw_x;
  const auto F_star_hw_y_  = hw_flux_formula_y (hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_hw_y;

  const auto F_star_hs_x_  = hs_flux_formula_x (hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_hs_x;
  const auto F_star_hs_y_  = hs_flux_formula_y (hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_hs_y;

  const auto F_star_Uxw_x_ = Uxw_flux_formula_x(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_, 0.) - diff_term_Uxw_x;
  const auto F_star_Uxw_y_ = Uxw_flux_formula_y(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_Uxw_y;

  const auto F_star_Uyw_x_ = Uyw_flux_formula_x(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_Uyw_x;
  const auto F_star_Uyw_y_ = Uyw_flux_formula_y(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_, 0.) - diff_term_Uyw_y;

  const auto F_star_Uxs_x_ = Uxs_flux_formula_x(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_, 0.) - diff_term_Uxs_x;
  const auto F_star_Uxs_y_ = Uxs_flux_formula_y(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_Uxs_y;

  const auto F_star_Uys_x_ = Uys_flux_formula_x(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_    ) - diff_term_Uys_x;
  const auto F_star_Uys_y_ = Uys_flux_formula_y(hw_cell_, hs_cell_, Uxw_cell_, Uyw_cell_, Uxs_cell_, Uys_cell_, 0.) - diff_term_Uys_y;
*/

  //
  const auto F_star_hw_x  = hw_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );
  const auto F_star_hw_y  = hw_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );

  const auto F_star_hs_x  = hs_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );
  const auto F_star_hs_y  = hs_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );

  const auto F_star_Uxw_x = Uxw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, 0.);
  const auto F_star_Uxw_y = Uxw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );

  const auto F_star_Uyw_x = Uyw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );
  const auto F_star_Uyw_y = Uyw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, 0.);

  const auto F_star_Uxs_x = Uxs_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, 0.);
  const auto F_star_Uxs_y = Uxs_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );

  const auto F_star_Uys_x = Uys_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell    );
  const auto F_star_Uys_y = Uys_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, 0.);



  for (int ii = 0; ii < 4; ++ii){
/*
    // high order fluxes
    const auto hw_h    = der_coeffs_x[ii]*F_star_hw_x +der_coeffs_y[ii]*F_star_hw_y;
    const auto Uxw_h   = der_coeffs_x[ii]*F_star_Uxw_x+der_coeffs_y[ii]*F_star_Uxw_y - .5*Dy*contr_x_w[ii]; // + .25*area*src_slope_formula(hw_cell, slope_x_cell);
    const auto Uyw_h   = der_coeffs_x[ii]*F_star_Uyw_x+der_coeffs_y[ii]*F_star_Uyw_y - .5*Dx*contr_y_w[ii];

    const auto hs_h    = der_coeffs_x[ii]*F_star_hs_x +der_coeffs_y[ii]*F_star_hs_y;
    const auto Uxs_h   = der_coeffs_x[ii]*F_star_Uxs_x+der_coeffs_y[ii]*F_star_Uxs_y - .5*Dy*contr_x_s[ii]; // + .25*area*src_slope_formula(hs_cell, slope_x_cell);
    const auto Uys_h   = der_coeffs_x[ii]*F_star_Uys_x+der_coeffs_y[ii]*F_star_Uys_y - .5*Dx*contr_y_s[ii];

    // low order fluxes
    const auto hw_    = der_coeffs_x[ii]*F_star_hw_x_ +der_coeffs_y[ii]*F_star_hw_y_;
    const auto Uxw_   = der_coeffs_x[ii]*F_star_Uxw_x_+der_coeffs_y[ii]*F_star_Uxw_y_ - .5*Dy*contr_x_w_[ii];
    const auto Uyw_   = der_coeffs_x[ii]*F_star_Uyw_x_+der_coeffs_y[ii]*F_star_Uyw_y_ - .5*Dx*contr_y_w_[ii];

    const auto hs_    = der_coeffs_x[ii]*F_star_hs_x_ +der_coeffs_y[ii]*F_star_hs_y_;
    const auto Uxs_   = der_coeffs_x[ii]*F_star_Uxs_x_+der_coeffs_y[ii]*F_star_Uxs_y_ - .5*Dy*contr_x_s_[ii];
    const auto Uys_   = der_coeffs_x[ii]*F_star_Uys_x_+der_coeffs_y[ii]*F_star_Uys_y_ - .5*Dx*contr_y_s_[ii];

    // anti-diffusion part
    const auto hw_al  = hw_h - hw_; 
    const auto Uxw_al = Uxw_h - Uxw_; 
    const auto Uyw_al = Uyw_h - Uyw_;

    const auto hs_al  = hs_h - hs_; 
    const auto Uxs_al = Uxs_h - Uxs_; 
    const auto Uys_al = Uys_h - Uys_; 
*/

    //if (contr_x_w[ii]!=contr_x_s[ii])
    //std::cout << contr_x_w[ii] << " " << contr_x_s[ii] << std::endl;

    const auto hw_    = der_coeffs_x[ii]*(F_star_hw_x -diff_term_hw_x )+der_coeffs_y[ii]*(F_star_hw_y -diff_term_hw_y);
    const auto Uxw_   = der_coeffs_x[ii]*(F_star_Uxw_x-diff_term_Uxw_x)+der_coeffs_y[ii]*(F_star_Uxw_y-diff_term_Uxw_y) - .5*Dy*contr_x_w[ii];
    const auto Uyw_   = der_coeffs_x[ii]*(F_star_Uyw_x-diff_term_Uyw_x)+der_coeffs_y[ii]*(F_star_Uyw_y-diff_term_Uyw_y) - .5*Dx*contr_y_w[ii]; 

    const auto hs_    = der_coeffs_x[ii]*(F_star_hs_x -diff_term_hs_x )+der_coeffs_y[ii]*(F_star_hs_y -diff_term_hs_y);
    const auto Uxs_   = der_coeffs_x[ii]*(F_star_Uxs_x-diff_term_Uxs_x)+der_coeffs_y[ii]*(F_star_Uxs_y-diff_term_Uxs_y) - .5*Dy*contr_x_s[ii];
    const auto Uys_   = der_coeffs_x[ii]*(F_star_Uys_x-diff_term_Uys_x)+der_coeffs_y[ii]*(F_star_Uys_y-diff_term_Uys_y) - .5*Dx*contr_y_s[ii];


    const auto hw_al  = der_coeffs_x[ii]*diff_term_hw_x  + der_coeffs_y[ii]*diff_term_hw_y; 
    const auto Uxw_al = der_coeffs_x[ii]*diff_term_Uxw_x + der_coeffs_y[ii]*diff_term_Uxw_y;
    const auto Uyw_al = der_coeffs_x[ii]*diff_term_Uyw_x + der_coeffs_y[ii]*diff_term_Uyw_y;

    const auto hs_al  = der_coeffs_x[ii]*diff_term_hs_x  + der_coeffs_y[ii]*diff_term_hs_y; 
    const auto Uxs_al = der_coeffs_x[ii]*diff_term_Uxs_x + der_coeffs_y[ii]*diff_term_Uxs_y;
    const auto Uys_al = der_coeffs_x[ii]*diff_term_Uys_x + der_coeffs_y[ii]*diff_term_Uys_y;


    incr_anti_diff[ordhw (index_quadrant_local)][ii] = hw_al;
    incr_anti_diff[ordUxw(index_quadrant_local)][ii] = Uxw_al;
    incr_anti_diff[ordUyw(index_quadrant_local)][ii] = Uyw_al;

    incr_anti_diff[ordhs (index_quadrant_local)][ii] = hs_al;
    incr_anti_diff[ordUxs(index_quadrant_local)][ii] = Uxs_al;
    incr_anti_diff[ordUys(index_quadrant_local)][ii] = Uys_al;


    if (! quadrant->is_hanging (ii)){

      incr [ordhw  (quadrant->gt (ii))] += hw_;
      incr [ordUxw (quadrant->gt (ii))] += Uxw_;
      incr [ordUyw (quadrant->gt (ii))] += Uyw_;

      P_plus [ordhw  (quadrant->gt (ii))] += std::max(0., hw_al );
      P_plus [ordUxw (quadrant->gt (ii))] += std::max(0., Uxw_al);
      P_plus [ordUyw (quadrant->gt (ii))] += std::max(0., Uyw_al);

      P_minus [ordhw  (quadrant->gt (ii))] += std::min(0., hw_al );
      P_minus [ordUxw (quadrant->gt (ii))] += std::min(0., Uxw_al);
      P_minus [ordUyw (quadrant->gt (ii))] += std::min(0., Uyw_al);


      incr [ordhs  (quadrant->gt (ii))] += hs_;
      incr [ordUxs (quadrant->gt (ii))] += Uxs_;
      incr [ordUys (quadrant->gt (ii))] += Uys_;

      P_plus [ordhs  (quadrant->gt (ii))] += std::max(0., hs_al );
      P_plus [ordUxs (quadrant->gt (ii))] += std::max(0., Uxs_al);
      P_plus [ordUys (quadrant->gt (ii))] += std::max(0., Uys_al);

      P_minus [ordhs  (quadrant->gt (ii))] += std::min(0., hs_al );
      P_minus [ordUxs (quadrant->gt (ii))] += std::min(0., Uxs_al);
      P_minus [ordUys (quadrant->gt (ii))] += std::min(0., Uys_al);
      
      
    } else {

      // w part
      incr [ordhw  (quadrant->gparent(0,ii))] += hw_;
      incr [ordhw  (quadrant->gparent(1,ii))] += hw_;
      
      incr [ordUxw (quadrant->gparent(0,ii))] += Uxw_;
      incr [ordUxw (quadrant->gparent(1,ii))] += Uxw_;
      
      incr [ordUyw (quadrant->gparent(0,ii))] += Uyw_;
      incr [ordUyw (quadrant->gparent(1,ii))] += Uyw_;



      P_plus [ordhw  (quadrant->gparent(0,ii))] += std::max(0., hw_al);
      P_plus [ordhw  (quadrant->gparent(1,ii))] += std::max(0., hw_al);
      
      P_plus [ordUxw (quadrant->gparent(0,ii))] += std::max(0., Uxw_al);
      P_plus [ordUxw (quadrant->gparent(1,ii))] += std::max(0., Uxw_al);
      
      P_plus [ordUyw (quadrant->gparent(0,ii))] += std::max(0., Uyw_al);
      P_plus [ordUyw (quadrant->gparent(1,ii))] += std::max(0., Uyw_al);



      P_minus [ordhw  (quadrant->gparent(0,ii))] += std::min(0., hw_al);
      P_minus [ordhw  (quadrant->gparent(1,ii))] += std::min(0., hw_al);
      
      P_minus [ordUxw (quadrant->gparent(0,ii))] += std::min(0., Uxw_al);
      P_minus [ordUxw (quadrant->gparent(1,ii))] += std::min(0., Uxw_al);
      
      P_minus [ordUyw (quadrant->gparent(0,ii))] += std::min(0., Uyw_al);
      P_minus [ordUyw (quadrant->gparent(1,ii))] += std::min(0., Uyw_al);



      // s part
      incr [ordhs  (quadrant->gparent(0,ii))] += hs_;
      incr [ordhs  (quadrant->gparent(1,ii))] += hs_;
      
      incr [ordUxs (quadrant->gparent(0,ii))] += Uxs_;
      incr [ordUxs (quadrant->gparent(1,ii))] += Uxs_;
      
      incr [ordUys (quadrant->gparent(0,ii))] += Uys_;
      incr [ordUys (quadrant->gparent(1,ii))] += Uys_;



      P_plus [ordhs  (quadrant->gparent(0,ii))] += std::max(0., hs_al);
      P_plus [ordhs  (quadrant->gparent(1,ii))] += std::max(0., hs_al);
      
      P_plus [ordUxs (quadrant->gparent(0,ii))] += std::max(0., Uxs_al);
      P_plus [ordUxs (quadrant->gparent(1,ii))] += std::max(0., Uxs_al);
      
      P_plus [ordUys (quadrant->gparent(0,ii))] += std::max(0., Uys_al);
      P_plus [ordUys (quadrant->gparent(1,ii))] += std::max(0., Uys_al);



      P_minus [ordhs  (quadrant->gparent(0,ii))] += std::min(0., hs_al);
      P_minus [ordhs  (quadrant->gparent(1,ii))] += std::min(0., hs_al);
      
      P_minus [ordUxs (quadrant->gparent(0,ii))] += std::min(0., Uxs_al);
      P_minus [ordUxs (quadrant->gparent(1,ii))] += std::min(0., Uxs_al);
      
      P_minus [ordUys (quadrant->gparent(0,ii))] += std::min(0., Uys_al);
      P_minus [ordUys (quadrant->gparent(1,ii))] += std::min(0., Uys_al);

    }
  }

}

void
TG2_scheme::terminate_second_step (tmesh::quadrant_iterator quadrant)
{
  index_quadrant = quadrant->get_global_quad_idx (); 

  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;

  dp_mean_vec = {0., 0., 0., 0.};
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      numerical_integration_pressure_2(ordBottom(quadrant->gt (ii)), dp_mean_vec[ii], excess_pore_water_pressure, 1.);
    }
    else
    {
      for (int jj=0; jj<=1; jj++)
        numerical_integration_pressure_2(ordBottom(quadrant->gparent(jj,ii)), dp_mean_vec[ii], excess_pore_water_pressure, .5);
    }
  }

  grad_cell_dp_mean = {.5 * ( (dp_mean_vec[3] - dp_mean_vec[2]) + (dp_mean_vec[1] - dp_mean_vec[0]) )/Dx, .5 * ( (dp_mean_vec[2] - dp_mean_vec[0]) + (dp_mean_vec[3] - dp_mean_vec[1]) )/Dy};

  const double & hw_cell = sol_onehalf[ordhw (index_quadrant)];
  const double & hs_cell = sol_onehalf[ordhs (index_quadrant)];

  const auto contr_x_w = hw_cell>epsilon ? -.5*hw_cell*grad_cell_dp_mean[0]/density_w*.25*area : 0.;
  const auto contr_y_w = hw_cell>epsilon ? -.5*hw_cell*grad_cell_dp_mean[1]/density_w*.25*area : 0.;
  const auto contr_x_s = hs_cell>epsilon ?  .5*hw_cell*grad_cell_dp_mean[0]/density_s*.25*area : 0.;
  const auto contr_y_s = hs_cell>epsilon ?  .5*hw_cell*grad_cell_dp_mean[1]/density_s*.25*area : 0.;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){

      incr [ordUxw (quadrant->gt (ii))] += contr_x_w;
      incr [ordUyw (quadrant->gt (ii))] += contr_y_w;
      incr [ordUxs (quadrant->gt (ii))] += contr_x_s;
      incr [ordUys (quadrant->gt (ii))] += contr_y_s;
    }
    else
    {
      // w
      incr [ordUxw (quadrant->gparent(0,ii))] += contr_x_w*.5;
      incr [ordUxw (quadrant->gparent(1,ii))] += contr_x_w*.5;
      
      incr [ordUyw (quadrant->gparent(0,ii))] += contr_y_w*.5;
      incr [ordUyw (quadrant->gparent(1,ii))] += contr_y_w*.5;

      // s part
      incr [ordUxs (quadrant->gparent(0,ii))] += contr_x_s*.5;
      incr [ordUxs (quadrant->gparent(1,ii))] += contr_x_s*.5;
      
      incr [ordUys (quadrant->gparent(0,ii))] += contr_y_s*.5;
      incr [ordUys (quadrant->gparent(1,ii))] += contr_y_s*.5;
    }
  }

}




void
TG2_scheme::second_step_pressure (tmesh::quadrant_iterator quadrant, const int& kkk)
{

  index_quadrant = quadrant->get_forest_quad_idx (); 

  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }

  
  for (int ii = 0; ii < 4; ++ii){ 

    double dp_kk_dof_c, P_plus_pressure_dp_kk_c, P_minus_pressure_dp_kk_c;

    if (! quadrant->is_hanging (ii)){
      hsdof [ii] = sol [ordhs  (quadrant->gt (ii) )];
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )];

      dp_kk_dof_c      = excess_pore_water_pressure [ordBottom    (quadrant->gt (ii)) + kkk];

      P_plus_pressure_dp_kk_c   = P_plus_pressure [ordBottom   (quadrant->gt (ii)) + kkk];
      P_minus_pressure_dp_kk_c  = P_minus_pressure[ordBottom   (quadrant->gt (ii)) + kkk];

    } else {
      hsdof [ii] = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);

      dp_kk_dof_c      = .5 * (excess_pore_water_pressure [ordBottom    (quadrant->gparent (0,ii)) + kkk] + 
                               excess_pore_water_pressure [ordBottom    (quadrant->gparent (1,ii)) + kkk]);

      P_plus_pressure_dp_kk_c   = .5 * (P_plus_pressure  [ordBottom (quadrant->gparent(0,ii)) + kkk] +
                                        P_plus_pressure  [ordBottom (quadrant->gparent(1,ii)) + kkk]);
      P_minus_pressure_dp_kk_c  = .5 * (P_minus_pressure [ordBottom (quadrant->gparent(0,ii)) + kkk] +
                                        P_minus_pressure [ordBottom (quadrant->gparent(1,ii)) + kkk]);
      
    }
    dp_kk_dof [ii] = dp_kk_dof_c;

    P_plus_pressure_dp_kk_dof  [ii] = P_plus_pressure_dp_kk_c;
    P_minus_pressure_dp_kk_dof [ii] = P_minus_pressure_dp_kk_c;
  }
  
  
  // compute local extrema
  const auto dp_kk_min_cell  = *std::min_element(dp_kk_dof.begin(), dp_kk_dof.end());
  const auto dp_kk_max_cell  = *std::max_element(dp_kk_dof.begin(), dp_kk_dof.end());

  bool is_node_in_element = false;

  std::array<double,4> dp_kk_min  = {dp_kk_min_cell, dp_kk_min_cell, dp_kk_min_cell, dp_kk_min_cell }, dp_kk_max  = {dp_kk_max_cell, dp_kk_max_cell, dp_kk_max_cell, dp_kk_max_cell };

  for (int ii = 0; ii < 4; ++ii){
    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {

      for (int jj = 0; jj < 4; ++jj) {
        if (xn[ii] == quadrant_nei->p(0, jj) && yn[ii] == quadrant_nei->p(1, jj))
        {
          is_node_in_element = true;
          break;
        }
      }

      if (is_node_in_element)
      {
        for (int jj = 0; jj < 4; ++jj) { 

          double dp_kk_current_cell;

          if (! quadrant_nei->is_hanging (jj)){
            dp_kk_current_cell  = excess_pore_water_pressure [ordBottom(quadrant_nei->gt(jj)) + kkk];
          } else {
            dp_kk_current_cell  = .5 * (excess_pore_water_pressure [ordBottom(quadrant_nei->gparent(0,jj)) + kkk] +
                                        excess_pore_water_pressure [ordBottom(quadrant_nei->gparent(1,jj)) + kkk]);
          }

          dp_kk_min [ii] = std::min(dp_kk_min [ii], dp_kk_current_cell);
          dp_kk_max [ii] = std::max(dp_kk_max [ii], dp_kk_current_cell);

        }

        is_node_in_element = false;
      }


    }
  }

  //std::cout << dp_kk_min[0] << " " << dp_kk_min[1] << " " << dp_kk_min[2] << " " << dp_kk_min[3] << " " << dp_kk_dof[0] << " " << dp_kk_dof[1] << " " << dp_kk_dof[2] << " " << dp_kk_dof[3] << " " << dp_kk_max[0] << " " << dp_kk_max[1] << " " << dp_kk_max[2] << " " << dp_kk_max[3] << std::endl;



  
  // compute flux correction
  double phi_cell_dp_kk = 1.; 
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_dp_kk  = incr_anti_diff_pressure[ordBottom(index_quadrant) + kkk][ii];

    //if (flux_on_the_node_dp_kk!=0)
    //std::cout << flux_on_the_node_dp_kk << std::endl;

    const auto vel_rusanov_cell_x = hsdof[ii]>epsilon ? Uxsdof[ii]/hsdof[ii] : 0.; 
    const auto vel_rusanov_cell_y = hsdof[ii]>epsilon ? Uysdof[ii]/hsdof[ii] : 0.; 
    const auto vel_square_rusanov_cell = std::abs(vel_rusanov_cell_x * vel_rusanov_cell_y);

    flux_limiter(dp_kk_min [ii], dp_kk_max [ii], dp_kk_dof[ii], P_plus_pressure_dp_kk_dof [ii], P_minus_pressure_dp_kk_dof  [ii], flux_on_the_node_dp_kk,  vel_square_rusanov_cell, phi_cell_dp_kk );
  }

  //if (phi_cell_dp_kk!=1 && phi_cell_dp_kk!=0)
  //std::cout << phi_cell_dp_kk << std::endl;


  for (int ii = 0; ii < 4; ++ii){

    const auto flux_on_the_node_dp_kk  = incr_anti_diff_pressure[ordBottom(index_quadrant) + kkk][ii]*phi_cell_dp_kk;

    if (! quadrant->is_hanging (ii)){

      excess_pore_water_pressure_incr [ordBottom(quadrant->gt(ii)) + kkk] += flux_on_the_node_dp_kk;

    } else {
      for (int jj=0; jj<=1; jj++)
      {
        excess_pore_water_pressure_incr [ordBottom(quadrant->gparent(jj,ii)) + kkk] += flux_on_the_node_dp_kk;
      }

    }

  }


}




void
TG2_scheme::second_step (tmesh::quadrant_iterator quadrant)
{
  index_quadrant = quadrant->get_global_quad_idx (); 
  index_quadrant_local = quadrant->get_forest_quad_idx (); 

  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;

  dp_mean_vec = {0., 0., 0., 0.};
  for (int ii = 0; ii < 4; ++ii){ 

    double hwdof_c, Uxwdof_c, Uywdof_c, P_plus_hw_c, P_minus_hw_c, P_plus_Uxw_c, P_minus_Uxw_c, P_plus_Uyw_c, P_minus_Uyw_c;
    double hsdof_c, Uxsdof_c, Uysdof_c, P_plus_hs_c, P_minus_hs_c, P_plus_Uxs_c, P_minus_Uxs_c, P_plus_Uys_c, P_minus_Uys_c;

    if (! quadrant->is_hanging (ii)){
      hwdof_c      = sol [ordhw    (quadrant->gt (ii))];
      Uxwdof_c     = sol [ordUxw   (quadrant->gt (ii))];
      Uywdof_c     = sol [ordUyw   (quadrant->gt (ii))];

      hsdof_c      = sol [ordhs    (quadrant->gt (ii))];
      Uxsdof_c     = sol [ordUxs   (quadrant->gt (ii))];
      Uysdof_c     = sol [ordUys   (quadrant->gt (ii))];

      Z_node[ii]    = Z [quadrant->gt (ii)];

      P_plus_hw_c   = P_plus [ordhw   (quadrant->gt (ii))];
      P_minus_hw_c  = P_minus[ordhw   (quadrant->gt (ii))];

      P_plus_Uxw_c   = P_plus [ordUxw   (quadrant->gt (ii))];
      P_minus_Uxw_c  = P_minus[ordUxw   (quadrant->gt (ii))];

      P_plus_Uyw_c   = P_plus [ordUyw   (quadrant->gt (ii))];
      P_minus_Uyw_c  = P_minus[ordUyw   (quadrant->gt (ii))];


      P_plus_hs_c   = P_plus [ordhs   (quadrant->gt (ii))];
      P_minus_hs_c  = P_minus[ordhs   (quadrant->gt (ii))];

      P_plus_Uxs_c   = P_plus [ordUxs   (quadrant->gt (ii))];
      P_minus_Uxs_c  = P_minus[ordUxs   (quadrant->gt (ii))];

      P_plus_Uys_c   = P_plus [ordUys   (quadrant->gt (ii))];
      P_minus_Uys_c  = P_minus[ordUys   (quadrant->gt (ii))];

      isdof_or_hanging[ii] = 1.;
      
      numerical_integration_pressure_2(ordBottom(quadrant->gt (ii)), dp_mean_vec[ii], excess_pore_water_pressure, 1.);

    } else {
      hwdof_c   = .5 * (sol [ordhw  (quadrant->gparent(0,ii))] +
                        sol [ordhw  (quadrant->gparent(1,ii))]);
      Uxwdof_c  = .5 * (sol [ordUxw (quadrant->gparent(0,ii))] +
                        sol [ordUxw (quadrant->gparent(1,ii))]);
      Uywdof_c  = .5 * (sol [ordUyw (quadrant->gparent(0,ii))] +
                        sol [ordUyw (quadrant->gparent(1,ii))]);


      hsdof_c   = .5 * (sol [ordhs  (quadrant->gparent(0,ii))] +
                        sol [ordhs  (quadrant->gparent(1,ii))]);
      Uxsdof_c  = .5 * (sol [ordUxs (quadrant->gparent(0,ii))] +
                        sol [ordUxs (quadrant->gparent(1,ii))]);
      Uysdof_c  = .5 * (sol [ordUys (quadrant->gparent(0,ii))] +
                        sol [ordUys (quadrant->gparent(1,ii))]);

      Z_node[ii] = .5 * (Z [quadrant->gparent(0,ii)] +
                         Z [quadrant->gparent(1,ii)]);

      P_plus_hw_c   = .5 * (P_plus [ordhw (quadrant->gparent(0,ii))] +
                            P_plus [ordhw (quadrant->gparent(1,ii))]);
      P_minus_hw_c  = .5 * (P_minus [ordhw (quadrant->gparent(0,ii))] +
                            P_minus [ordhw (quadrant->gparent(1,ii))]);

      P_plus_Uxw_c   = .5 * (P_plus [ordUxw (quadrant->gparent(0,ii))] +
                             P_plus [ordUxw (quadrant->gparent(1,ii))]);
      P_minus_Uxw_c  = .5 * (P_minus [ordUxw (quadrant->gparent(0,ii))] +
                             P_minus [ordUxw (quadrant->gparent(1,ii))]);

      P_plus_Uyw_c   = .5 * (P_plus [ordUyw (quadrant->gparent(0,ii))] +
                             P_plus [ordUyw (quadrant->gparent(1,ii))]);
      P_minus_Uyw_c  = .5 * (P_minus [ordUyw (quadrant->gparent(0,ii))] +
                             P_minus [ordUyw (quadrant->gparent(1,ii))]);



      P_plus_hs_c   = .5 * (P_plus [ordhs (quadrant->gparent(0,ii))] +
                            P_plus [ordhs (quadrant->gparent(1,ii))]);
      P_minus_hs_c  = .5 * (P_minus [ordhs (quadrant->gparent(0,ii))] +
                            P_minus [ordhs (quadrant->gparent(1,ii))]);

      P_plus_Uxs_c   = .5 * (P_plus [ordUxs (quadrant->gparent(0,ii))] +
                             P_plus [ordUxs (quadrant->gparent(1,ii))]);
      P_minus_Uxs_c  = .5 * (P_minus [ordUxs (quadrant->gparent(0,ii))] +
                             P_minus [ordUxs (quadrant->gparent(1,ii))]);

      P_plus_Uys_c   = .5 * (P_plus [ordUys (quadrant->gparent(0,ii))] +
                             P_plus [ordUys (quadrant->gparent(1,ii))]);
      P_minus_Uys_c  = .5 * (P_minus [ordUys (quadrant->gparent(0,ii))] +
                             P_minus [ordUys (quadrant->gparent(1,ii))]);

      isdof_or_hanging[ii] = .5;

      for (int jj=0; jj<=1; jj++)
        numerical_integration_pressure_2(ordBottom(quadrant->gparent (jj,ii)), dp_mean_vec[ii], excess_pore_water_pressure, .5);
      
    }

    double hdof_c  = hwdof_c+hsdof_c;
    double ndof_c  = hdof_c>epsilon ? hwdof_c/hdof_c : 0.;
    double nsdof_c = hdof_c>epsilon ? hsdof_c/hdof_c : 0.;

    hdof       [ii] = hdof_c;
    hwdof      [ii] = hwdof_c;
    hsdof      [ii] = hsdof_c;
    etawdof    [ii] = hdof_c>epsilon ? hwdof_c+Z_node[ii]*ndof_c *g_coeff : hwdof_c;
    etasdof    [ii] = hdof_c>epsilon ? hsdof_c+Z_node[ii]*nsdof_c*g_coeff : hsdof_c;
    Uxwdof     [ii] = Uxwdof_c;
    Uywdof     [ii] = Uywdof_c;
    Uxsdof     [ii] = Uxsdof_c;
    Uysdof     [ii] = Uysdof_c;



    P_plus_hw_dof  [ii] = P_plus_hw_c;
    P_minus_hw_dof [ii] = P_minus_hw_c;

    P_plus_hs_dof  [ii] = P_plus_hs_c;
    P_minus_hs_dof [ii] = P_minus_hs_c;

    P_plus_Uxw_dof [ii] = P_plus_Uxw_c;
    P_minus_Uxw_dof[ii] = P_minus_Uxw_c;

    P_plus_Uyw_dof [ii] = P_plus_Uyw_c;
    P_minus_Uyw_dof[ii] = P_minus_Uyw_c;

    P_plus_Uxs_dof [ii] = P_plus_Uxs_c;
    P_minus_Uxs_dof[ii] = P_minus_Uxs_c;

    P_plus_Uys_dof [ii] = P_plus_Uys_c;
    P_minus_Uys_dof[ii] = P_minus_Uys_c;

  }
  
  
  // compute local extrema
  const auto hw_min_cell  = *std::min_element(etawdof.begin(), etawdof.end());
  const auto hw_max_cell  = *std::max_element(etawdof.begin(), etawdof.end());

  const auto hs_min_cell  = *std::min_element(etasdof.begin(), etasdof.end());
  const auto hs_max_cell  = *std::max_element(etasdof.begin(), etasdof.end());

  const auto Uxw_min_cell = *std::min_element(Uxwdof.begin(),  Uxwdof.end());
  const auto Uxw_max_cell = *std::max_element(Uxwdof.begin(),  Uxwdof.end());

  const auto Uyw_min_cell = *std::min_element(Uywdof.begin(),  Uywdof.end());
  const auto Uyw_max_cell = *std::max_element(Uywdof.begin(),  Uywdof.end());

  const auto Uxs_min_cell = *std::min_element(Uxsdof.begin(),  Uxsdof.end());
  const auto Uxs_max_cell = *std::max_element(Uxsdof.begin(),  Uxsdof.end());

  const auto Uys_min_cell = *std::min_element(Uysdof.begin(),  Uysdof.end());
  const auto Uys_max_cell = *std::max_element(Uysdof.begin(),  Uysdof.end());   

  bool is_node_in_element = false;

  std::array<double,4> hw_min  = {hw_min_cell, hw_min_cell, hw_min_cell, hw_min_cell }, hw_max  = {hw_max_cell, hw_max_cell, hw_max_cell, hw_max_cell },
                       hs_min  = {hs_min_cell, hs_min_cell, hs_min_cell, hs_min_cell }, hs_max  = {hs_max_cell, hs_max_cell, hs_max_cell, hs_max_cell },
                       Uxw_min = {Uxw_min_cell,Uxw_min_cell,Uxw_min_cell,Uxw_min_cell}, Uxw_max = {Uxw_max_cell,Uxw_max_cell,Uxw_max_cell,Uxw_max_cell},
                       Uyw_min = {Uyw_min_cell,Uyw_min_cell,Uyw_min_cell,Uyw_min_cell}, Uyw_max = {Uyw_max_cell,Uyw_max_cell,Uyw_max_cell,Uyw_max_cell},
                       Uxs_min = {Uxs_min_cell,Uxs_min_cell,Uxs_min_cell,Uxs_min_cell}, Uxs_max = {Uxs_max_cell,Uxs_max_cell,Uxs_max_cell,Uxs_max_cell},
                       Uys_min = {Uys_min_cell,Uys_min_cell,Uys_min_cell,Uys_min_cell}, Uys_max = {Uys_max_cell,Uys_max_cell,Uys_max_cell,Uys_max_cell};

  for (int ii = 0; ii < 4; ++ii){
    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {

      for (int jj = 0; jj < 4; ++jj) {
        if (xn[ii] == quadrant_nei->p(0, jj) && yn[ii] == quadrant_nei->p(1, jj))
        {
          is_node_in_element = true;
          break;
        }
      }

      if (is_node_in_element)
      {
        for (int jj = 0; jj < 4; ++jj) { 

          double n_current_cell, ns_current_cell, Z_current_cell, h_current_cell, hw_current_cell, hs_current_cell, Uxw_current_cell, Uyw_current_cell, Uxs_current_cell, Uys_current_cell;

          if (! quadrant_nei->is_hanging (jj)){

            hw_current_cell  = sol [ordhw  (quadrant_nei->gt (jj))];
            hs_current_cell  = sol [ordhs  (quadrant_nei->gt (jj))];
            Uxw_current_cell = sol [ordUxw (quadrant_nei->gt (jj))];
            Uyw_current_cell = sol [ordUyw (quadrant_nei->gt (jj))];
            Uxs_current_cell = sol [ordUxs (quadrant_nei->gt (jj))];
            Uys_current_cell = sol [ordUys (quadrant_nei->gt (jj))];
            Z_current_cell   = Z   [        quadrant_nei->gt (jj) ];
          } else {

            hw_current_cell  = .5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))] +
                                     sol [ordhw  (quadrant_nei->gparent(1,jj))]);
            hs_current_cell  = .5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] +
                                     sol [ordhs  (quadrant_nei->gparent(1,jj))]);
            Uxw_current_cell = .5 * (sol [ordUxw (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUxw (quadrant_nei->gparent(1,jj))]);
            Uyw_current_cell = .5 * (sol [ordUyw (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUyw (quadrant_nei->gparent(1,jj))]);
            Uxs_current_cell = .5 * (sol [ordUxs (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUxs (quadrant_nei->gparent(1,jj))]);
            Uys_current_cell = .5 * (sol [ordUys (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUys (quadrant_nei->gparent(1,jj))]);
            Z_current_cell   = .5 * (Z   [        quadrant_nei->gparent(0,jj) ] +
                                     Z   [        quadrant_nei->gparent(1,jj) ]);
          }

          h_current_cell  = hw_current_cell + hs_current_cell;
          n_current_cell  = h_current_cell>epsilon ? hw_current_cell/h_current_cell : 0.; 
          ns_current_cell = h_current_cell>epsilon ? hs_current_cell/h_current_cell : 0.;

          hw_current_cell = hw_current_cell + Z_current_cell*n_current_cell *g_coeff;
          hs_current_cell = hs_current_cell + Z_current_cell*ns_current_cell*g_coeff;

          hw_min [ii] = std::min(hw_min [ii], hw_current_cell);
          hw_max [ii] = std::max(hw_max [ii], hw_current_cell);

          hs_min [ii] = std::min(hs_min [ii], hs_current_cell);
          hs_max [ii] = std::max(hs_max [ii], hs_current_cell);

          Uxw_min[ii] = std::min(Uxw_min[ii], Uxw_current_cell);
          Uxw_max[ii] = std::max(Uxw_max[ii], Uxw_current_cell);

          Uyw_min[ii] = std::min(Uyw_min[ii], Uyw_current_cell);
          Uyw_max[ii] = std::max(Uyw_max[ii], Uyw_current_cell);

          Uxs_min[ii] = std::min(Uxs_min[ii], Uxs_current_cell);
          Uxs_max[ii] = std::max(Uxs_max[ii], Uxs_current_cell);

          Uys_min[ii] = std::min(Uys_min[ii], Uys_current_cell);
          Uys_max[ii] = std::max(Uys_max[ii], Uys_current_cell);

        }

        is_node_in_element = false;
      }


    }
  }



  
  // compute flux correction
  double phi_cell_hw = 1., phi_cell_hs = 1., phi_cell_Uxw = 1., phi_cell_Uyw = 1., phi_cell_Uxs = 1., phi_cell_Uys = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_hw  = incr_anti_diff[ordhw (index_quadrant_local)][ii];
    const auto & flux_on_the_node_hs  = incr_anti_diff[ordhs (index_quadrant_local)][ii];
    const auto & flux_on_the_node_Uxw = incr_anti_diff[ordUxw(index_quadrant_local)][ii];
    const auto & flux_on_the_node_Uyw = incr_anti_diff[ordUyw(index_quadrant_local)][ii];
    const auto & flux_on_the_node_Uxs = incr_anti_diff[ordUxs(index_quadrant_local)][ii];
    const auto & flux_on_the_node_Uys = incr_anti_diff[ordUys(index_quadrant_local)][ii];


    const auto speed = max_eigen(hwdof[ii], hsdof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii], dp_mean_vec[ii]);
    const auto vel_rusanov_cell_x = speed[0]; //hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = speed[1]; //hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_square_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(hw_min [ii], hw_max [ii], etawdof[ii], P_plus_hw_dof [ii], P_minus_hw_dof  [ii], flux_on_the_node_hw,  vel_square_rusanov_cell, phi_cell_hw );
    flux_limiter(hs_min [ii], hs_max [ii], etasdof[ii], P_plus_hs_dof [ii], P_minus_hs_dof  [ii], flux_on_the_node_hs,  vel_square_rusanov_cell, phi_cell_hs );
    flux_limiter(Uxw_min[ii], Uxw_max[ii], Uxwdof [ii], P_plus_Uxw_dof[ii], P_minus_Uxw_dof [ii], flux_on_the_node_Uxw, vel_square_rusanov_cell, phi_cell_Uxw);
    flux_limiter(Uyw_min[ii], Uyw_max[ii], Uywdof [ii], P_plus_Uyw_dof[ii], P_minus_Uyw_dof [ii], flux_on_the_node_Uyw, vel_square_rusanov_cell, phi_cell_Uyw);
    flux_limiter(Uxs_min[ii], Uxs_max[ii], Uxsdof [ii], P_plus_Uxs_dof[ii], P_minus_Uxs_dof [ii], flux_on_the_node_Uxs, vel_square_rusanov_cell, phi_cell_Uxs);
    flux_limiter(Uys_min[ii], Uys_max[ii], Uysdof [ii], P_plus_Uys_dof[ii], P_minus_Uys_dof [ii], flux_on_the_node_Uys, vel_square_rusanov_cell, phi_cell_Uys);
  }

  //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;

  //phi_cell_hw = 0., phi_cell_hs = 0., phi_cell_Uxw = 0., phi_cell_Uyw = 0., phi_cell_Uxs = 0., phi_cell_Uys = 0.; 
  //if (phi_cell_hw!=1) std::cout << phi_cell_hw << std::endl;
  //phi_cell_hw = 1., phi_cell_hs = 1.;
  //phi_cell_Uxw = 1., phi_cell_Uyw = 1., phi_cell_Uxs = 1., phi_cell_Uys = 1.;
  //phi_cell_hw = 0.5, phi_cell_hs = 0.5, phi_cell_Uxw = .5, phi_cell_Uyw = 1., phi_cell_Uxs = .5, phi_cell_Uys = 1.; 
  //std::cout << "we put limiter equal to zero" << std::endl;


  //if (phi_cell_hw!=1 && phi_cell_hw!=0)
  //std::cout << phi_cell_hw << " " << phi_cell_Uxw << std::endl;

  const double & hw_cell    = sol_onehalf[ordhw    (index_quadrant)];
  const double & Uxw_cell   = sol_onehalf[ordUxw   (index_quadrant)];
  const double & Uyw_cell   = sol_onehalf[ordUyw   (index_quadrant)];

  const double & hs_cell    = sol_onehalf[ordhs    (index_quadrant)];
  const double & Uxs_cell   = sol_onehalf[ordUxs   (index_quadrant)];
  const double & Uys_cell   = sol_onehalf[ordUys   (index_quadrant)];

  const double & bed_excess_pore_water_pressure = excess_pore_water_pressure_onehalf[ordBottom(index_quadrant)];

  for (int ii = 0; ii < 4; ++ii){

    //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;

    const auto flux_on_the_node_hw  = incr_anti_diff[ordhw (index_quadrant_local)][ii]*phi_cell_hw;
    const auto flux_on_the_node_hs  = incr_anti_diff[ordhs (index_quadrant_local)][ii]*phi_cell_hs;
    const auto flux_on_the_node_Uxw = incr_anti_diff[ordUxw(index_quadrant_local)][ii]*phi_cell_Uxw + .25*area*isdof_or_hanging[ii]*tau_ccc*Uxw_src_formula(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell);
    const auto flux_on_the_node_Uyw = incr_anti_diff[ordUyw(index_quadrant_local)][ii]*phi_cell_Uyw + .25*area*isdof_or_hanging[ii]*tau_ccc*Uyw_src_formula(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell);
    const auto flux_on_the_node_Uxs = incr_anti_diff[ordUxs(index_quadrant_local)][ii]*phi_cell_Uxs + .25*area*isdof_or_hanging[ii]*tau_ccc*Uxs_src_formula(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, bed_excess_pore_water_pressure);
    const auto flux_on_the_node_Uys = incr_anti_diff[ordUys(index_quadrant_local)][ii]*phi_cell_Uys + .25*area*isdof_or_hanging[ii]*tau_ccc*Uys_src_formula(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell, bed_excess_pore_water_pressure);


    if (! quadrant->is_hanging (ii)){

      incr [ordhw  (quadrant->gt (ii))] += flux_on_the_node_hw;
      incr [ordhs  (quadrant->gt (ii))] += flux_on_the_node_hs;
      incr [ordUxw (quadrant->gt (ii))] += flux_on_the_node_Uxw;
      incr [ordUyw (quadrant->gt (ii))] += flux_on_the_node_Uyw;
      incr [ordUxs (quadrant->gt (ii))] += flux_on_the_node_Uxs;
      incr [ordUys (quadrant->gt (ii))] += flux_on_the_node_Uys;

    } else {

      incr [ordhw  (quadrant->gparent(0,ii))] += flux_on_the_node_hw;
      incr [ordhw  (quadrant->gparent(1,ii))] += flux_on_the_node_hw;

      incr [ordhs  (quadrant->gparent(0,ii))] += flux_on_the_node_hs;
      incr [ordhs  (quadrant->gparent(1,ii))] += flux_on_the_node_hs;
      
      incr [ordUxw (quadrant->gparent(0,ii))] += flux_on_the_node_Uxw;
      incr [ordUxw (quadrant->gparent(1,ii))] += flux_on_the_node_Uxw;
      
      incr [ordUyw (quadrant->gparent(0,ii))] += flux_on_the_node_Uyw;
      incr [ordUyw (quadrant->gparent(1,ii))] += flux_on_the_node_Uyw;

      incr [ordUxs (quadrant->gparent(0,ii))] += flux_on_the_node_Uxs;
      incr [ordUxs (quadrant->gparent(1,ii))] += flux_on_the_node_Uxs;
      
      incr [ordUys (quadrant->gparent(0,ii))] += flux_on_the_node_Uys;
      incr [ordUys (quadrant->gparent(1,ii))] += flux_on_the_node_Uys;

    }

  }


}


void
TG2_scheme::flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& vel_square_rusanov_cell, double& phi_cell_Q)
{
  const auto Q_plus  = (Q_max-Q_dof)*dt*vel_square_rusanov_cell;
  const auto Q_minus = (Q_min-Q_dof)*dt*vel_square_rusanov_cell;

  const auto R_plus  = P_plus_Q ==0 ? 1. : std::min(1., Q_plus /P_plus_Q );
  const auto R_minus = P_minus_Q==0 ? 1. : std::min(1., Q_minus/P_minus_Q);

  //const auto R_plus  = (P_plus_Q ==0 || Q_plus ==0) ? 1 : std::min(1., Q_plus /P_plus_Q );
  //const auto R_minus = (P_minus_Q==0 || Q_minus==0) ? 1 : std::min(1., Q_minus/P_minus_Q);

  phi_cell_Q = std::min(phi_cell_Q, flux_on_the_node>=0 ? R_plus : R_minus);

  //std::cout << Q_min << " " << Q_max << " " << Q_dof << std::endl;

}


void
TG2_scheme::stabilization_term (const int& kk)
{
  // apply the corrector step now, it's like adding an extra interphase drag
  const auto & hw_c  = sol.get_owned_data ()[kk  ];
  const auto & hs_c  = sol.get_owned_data ()[kk+1];
  const auto & Uxw_c = sol.get_owned_data ()[kk+2];
  const auto & Uyw_c = sol.get_owned_data ()[kk+3];
  const auto & Uxs_c = sol.get_owned_data ()[kk+4];
  const auto & Uys_c = sol.get_owned_data ()[kk+5];


  double dp_mean = 0.;
  numerical_integration_pressure((kk/6)*number_FD_points, dp_mean, excess_pore_water_pressure, 1.);

  const double h_c = hw_c+hs_c;

  const double vel_w_x = hw_c>epsilon ? Uxw_c/hw_c : 0.;
  const double vel_s_x = hs_c>epsilon ? Uxs_c/hs_c : 0.;
  const double vel_w_y = hw_c>epsilon ? Uyw_c/hw_c : 0.;
  const double vel_s_y = hs_c>epsilon ? Uys_c/hs_c : 0.;

  const auto n  = h_c>epsilon ? hw_c/h_c : 0.;
  const auto ns = h_c>epsilon ? hs_c/h_c : 0.;

  const double kinematic_speed_wave = std::sqrt(grav*h_c);
  const double beta_coeff  = std::sqrt(.5*n*(1.-r_coeff));

  const auto density = ns*density_s + n*density_w;

  const auto abs_delta_vel_x = std::abs(vel_w_x - vel_s_x);
  const auto abs_delta_vel_y = std::abs(vel_w_y - vel_s_y);


  const double beta_coeff_square = beta_coeff*beta_coeff;
  const double celerity_square = kinematic_speed_wave*kinematic_speed_wave;

  const double b_coeff = 4.*beta_coeff_square*dp_mean + 2.*celerity_square*density_w*(1.+beta_coeff_square);
  const double c_coeff = b_coeff*b_coeff - 16*celerity_square*beta_coeff_square*density_w*density_w*(celerity_square + dp_mean/density_w); //16.*dp_mean*dp_mean*beta_coeff_square*beta_coeff_square + 16.*celerity_square*beta_coeff_square*beta_coeff_square*density_w*dp_mean + 4.*celerity_square*celerity_square*density_w*density_w*(1.-beta_coeff_square)*(1.-beta_coeff_square);


  //const auto common_coeff = b_coeff/density_w;
  //const auto second_common_coeff = std::sqrt(common_coeff*common_coeff - 16.*celerity_square*beta_coeff_square*(celerity_square+dp_mean/density_w));
  const auto common_coeff = std::sqrt(c_coeff);
  const auto c_1 = .5*std::sqrt(b_coeff + common_coeff)/std::sqrt(density_w);
  const auto c_2 = .5*std::sqrt(b_coeff - common_coeff)/std::sqrt(density_w);

  //if (c_coeff<0 || b_coeff<common_coeff)
  //{
  //  std::cout << c_coeff << " " << b_coeff - common_coeff << " " << celerity_square + dp_mean/density_w << " " << c_2 << std::endl;
  //  exit(1);
  //}

  const auto x_coeff_1 = 2.*c_1;
  const auto x_coeff_2 = 2.*c_2;

  const auto hyp_diff_x  = abs_delta_vel_x - x_coeff_2;
  const auto hyp_diff_y  = abs_delta_vel_y - x_coeff_2;
  //const auto hyp_diff_x_ = abs_delta_vel_x - x_coeff_1;
  //const auto hyp_diff_y_ = abs_delta_vel_y - x_coeff_1;


  //const auto hyp_diff_x  = abs_delta_vel_x - 2.*kinematic_speed_wave*beta_coeff;
  //const auto hyp_diff_y  = abs_delta_vel_y - 2.*kinematic_speed_wave*beta_coeff;
  //const auto hyp_diff_x_ = abs_delta_vel_x - 2.*kinematic_speed_wave;
  //const auto hyp_diff_y_ = abs_delta_vel_y - 2.*kinematic_speed_wave;


  //const double cx_sgn = std::max( h_c>epsilon && hyp_diff_x_<0 && hw_c>epsilon ? n*ns*hyp_diff_x/(dt*(n*r_coeff+ns))*.5/kinematic_speed_wave/beta_coeff : 0., 0.);
  //const double cy_sgn = std::max( h_c>epsilon && hyp_diff_y_<0 && hw_c>epsilon ? n*ns*hyp_diff_y/(dt*(n*r_coeff+ns))*.5/kinematic_speed_wave/beta_coeff : 0., 0.);




  const double cx_sgn = std::max( (h_c>epsilon && hw_c>epsilon && hs_c>epsilon) ? n*ns*hyp_diff_x/(dt*(n*r_coeff+ns))/x_coeff_2 : 0., 0.);
  const double cy_sgn = std::max( (h_c>epsilon && hw_c>epsilon && hs_c>epsilon) ? n*ns*hyp_diff_y/(dt*(n*r_coeff+ns))/x_coeff_2 : 0., 0.);

  //const double cx_sgn_2 = std::min( (h_c>epsilon && hw_c>epsilon && hs_c>epsilon && hyp_diff_x>0) ? n*ns*hyp_diff_x_/(dt*(n*r_coeff+ns))/x_coeff_1 : 0., 0.);
  //const double cy_sgn_2 = std::min( (h_c>epsilon && hw_c>epsilon && hs_c>epsilon && hyp_diff_y>0) ? n*ns*hyp_diff_y_/(dt*(n*r_coeff+ns))/x_coeff_1 : 0., 0.);

  //cx_sgn = std::min(std::abs(cx_sgn), std::abs(cx_sgn_2));
  //cy_sgn = std::min(std::abs(cy_sgn), std::abs(cy_sgn_2));


  //const double cx_sgn = (h_c>epsilon && hw_c>epsilon && hs_c>epsilon) ? 1e3 : 0.;
  //const double cy_sgn = (h_c>epsilon && hw_c>epsilon && hs_c>epsilon) ? 1e3 : 0.;

  //if (cx_sgn!=0)
  //{
  //  std::cout << cx_sgn << std::endl;
  //}


  //const double cx_sgn = std::max( h_c>epsilon && hyp_diff_x>0 && hyp_diff_x_<0 ? 1.e7 : 0., 0.);
  //const double cy_sgn = std::max( h_c>epsilon && hyp_diff_y>0 && hyp_diff_y_<0 ? 1.e7 : 0., 0.);


  const double big_Ax = 1. + (hw_c>epsilon ? dt*cx_sgn*h_c/hw_c : 0.);
  const double big_Bx = hs_c>epsilon ? -dt*cx_sgn*h_c/hs_c : 0.;
  const double big_Cx = hw_c>epsilon ? -dt*cx_sgn*h_c*r_coeff/hw_c : 0.;
  const double big_Dx = 1. + (hs_c>epsilon ? dt*cx_sgn*h_c*r_coeff/hs_c : 0.);

  const double rhs_1x = Uxw_c;
  const double rhs_2x = Uxs_c;

  const double big_detx = big_Ax*big_Dx-big_Bx*big_Cx; 

  const double big_Ay = 1. + (hw_c>epsilon ? dt*cy_sgn*h_c/hw_c : 0.);
  const double big_By = hs_c>epsilon ? -dt*cy_sgn*h_c/hs_c : 0.;
  const double big_Cy = hw_c>epsilon ? -dt*cy_sgn*h_c*r_coeff/hw_c : 0.;
  const double big_Dy = 1. + (hs_c>epsilon ? dt*cy_sgn*h_c*r_coeff/hs_c : 0.);

  const double rhs_1y = Uyw_c;
  const double rhs_2y = Uys_c;

  const double big_dety = big_Ay*big_Dy-big_By*big_Cy;

  sol.get_owned_data ()[kk+2] = (rhs_1x*big_Dx-rhs_2x*big_Bx)/big_detx;
  sol.get_owned_data ()[kk+3] = (rhs_1y*big_Dy-rhs_2y*big_By)/big_dety;
  sol.get_owned_data ()[kk+4] = (rhs_2x*big_Ax-rhs_1x*big_Cx)/big_detx;
  sol.get_owned_data ()[kk+5] = (rhs_2y*big_Ay-rhs_1y*big_Cy)/big_dety;

}




void
TG2_scheme::set_dt (const double dt_)
{ dt = dt_; }

void
TG2_scheme::set_tau ()
{ 
  tau = dt*.5; 
  tau_c = dt*(1. - std::sqrt(2.)*.5);
  tau_cc = tau*(std::sqrt(2.)-1.);
  tau_ccc = (std::sqrt(2.)-1.);
}

void
TG2_scheme::set_old_dt (const double dt_)
{ dt_old = dt_; }

void
TG2_scheme::set_times(const double& t, const double& td, const double& tdd)
{ time = t; timed = td; timedd = tdd; }

double
TG2_scheme::get_dt ()
{ return dt; }



// flux functions
double
TG2_scheme::hw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys) 
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (hw>epsilon ? Uxw : 0.); 
}

double
TG2_scheme::hw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (hw>epsilon ? Uyw : 0.); 
}

double
TG2_scheme::hs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  return (hs>epsilon ? Uxs : 0.); 
}

double
TG2_scheme::hs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  return (hs>epsilon ? Uys : 0.); 
}

double
TG2_scheme::Uxw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean)
{ 
  const auto vel_x = hw>epsilon ? Uxw/hw : 0.;
  return (hw>epsilon ? Uxw*vel_x + grav*hw*hw/2. + hw/density_w*dp_mean : 0.);
}
 
double
TG2_scheme::Uxw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hw>epsilon ? Uyw*Uxw/hw : 0.); }

double
TG2_scheme::Uyw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hw>epsilon ? Uyw*Uxw/hw : 0.); }

double
TG2_scheme::Uyw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean)
{ 
  const auto vel_y = hw>epsilon ? Uyw/hw : 0.;
  return (hw>epsilon ? Uyw*vel_y + grav*hw*hw/2. + hw/density_w*dp_mean : 0.); 
}


double
TG2_scheme::Uxs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean)
{ 
  const auto vel_x = hs>epsilon ? Uxs/hs : 0.;
  return (hs>epsilon ? Uxs*vel_x + grav*hs*hs/2. + grav*(1-r_coeff)*hs*hw/2. - hw/density_s*dp_mean : 0.);  
}
 
double
TG2_scheme::Uxs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hs>epsilon ? Uys*Uxs/hs : 0.); }

double
TG2_scheme::Uys_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hs>epsilon ? Uys*Uxs/hs : 0.); }

double
TG2_scheme::Uys_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& dp_mean)
{ 
  const auto vel_y = hs>epsilon ? Uys/hs : 0.;
  return (hs>epsilon ? Uys*vel_y + grav*hs*hs/2. + grav*(1.-r_coeff)*hs*hw/2. - hw/density_s*dp_mean : 0.); 
}


// source terms
double
TG2_scheme::hw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto h = hw+hs;
  const auto Ux = Uxw + Uxs;
  const auto Uy = Uyw + Uys;
  const auto nw = h>epsilon ? hw/h : 0.;
  return (nw*erosion_coefficient*std::sqrt(Ux*Ux + Uy*Uy)); 
}

double
TG2_scheme::hs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto h = hw+hs;
  const auto Ux = Uxw + Uxs;
  const auto Uy = Uyw + Uys;
  const auto ns = h>epsilon ? hs/h : 0.;
  return (ns*erosion_coefficient*std::sqrt(Ux*Ux + Uy*Uy)); 
}

double
TG2_scheme::Uxw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;

  const double vel_w_x = hw>epsilon ? Uxw/hw : 0.;
  const double vel_s_x = hs>epsilon ? Uxs/hs : 0.;

  const double C_d = (hw>epsilon && h>epsilon && hs>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(density_s/density_w-1.)*grav : 0.;
  const double R_x = C_d*(vel_w_x-vel_s_x);


  return ( -h*R_x );
}

double
TG2_scheme::Uyw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;

  const double vel_w_y = hw>epsilon ? Uyw/hw : 0.;
  const double vel_s_y = hs>epsilon ? Uys/hs : 0.;

  const double C_d = (hw>epsilon && h>epsilon && hs>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(density_s/density_w-1.)*grav : 0.;
  const double R_y = C_d*(vel_w_y-vel_s_y);

  return ( -h*R_y );
}

double
TG2_scheme::Uxs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;
  const double bed_pressure = grav*hs*(1.-r_coeff) - bed_excess_pore_water_pressure/density_s; 
  const auto Ux = Uxw + Uxs;
  const auto Uy = Uyw + Uys;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::abs( vel_x );


  const double vel_w_x = hw>epsilon ? Uxw/hw : 0.;
  const double vel_s_x = hs>epsilon ? Uxs/hs : 0.;

  const double density = ns + n*r_coeff;

  const double C_d = (hw>epsilon && h>epsilon && hs>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(1.-r_coeff)*grav : 0.;
  const double R_x = C_d*(vel_w_x-vel_s_x);

  //std::cout << h << " " << dhdx << " " << dZdx << " " << grav*h*(dZdx+dhdx) << std::endl;

  const double vel_x_sign = 2./M_PI*std::atan(M_PI*.5*Ux); //std::abs(Ux)>tolerance_sign ? Ux/std::abs(Ux) : Ux/tolerance_sign;
  //const double vel_x_sign = abs_vel>tolerance_sign ? vel_x/abs_vel : 0.;
  //const double vel_x_sign = abs_vel>tolerance_sign ? vel_x/abs_vel : vel_x/tolerance_sign;

  //const double bed_fric_contr = is_bed_friction ? vel_x_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = (h*h)>epsilon && is_bed_friction && hs>epsilon ? density*Ux*grav*std::abs(Ux)/turbulence_coeff/h/h : 0.; //is_bed_friction ? vel_x*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction && hs>epsilon ? vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  return ( - bed_fric_contr_one - bed_fric_contr_two + h*R_x );
}

double
TG2_scheme::Uys_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;
  const double bed_pressure = grav*hs*(1.-r_coeff) - bed_excess_pore_water_pressure/density_s;
  const auto Ux = Uxw + Uxs;
  const auto Uy = Uyw + Uys;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::abs( vel_y );



  const double vel_w_y = hw>epsilon ? Uyw/hw : 0.;
  const double vel_s_y = hs>epsilon ? Uys/hs : 0.;

  const double density = ns + n*r_coeff;

  const double C_d = (hw>epsilon && h>epsilon && hs>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(1.-r_coeff)*grav : 0.;

  const double R_y = C_d*(vel_w_y-vel_s_y);

  const double vel_y_sign = 2./M_PI*std::atan(M_PI*.5*Uy); //std::abs(Uy)>tolerance_sign ? Uy/std::abs(Uy) : Uy/tolerance_sign;
  //const double vel_y_sign = abs_vel>tolerance_sign ? vel_y/abs_vel : vel_y/tolerance_sign;
  //const double vel_y_sign = (vel_y > ) ? 1.0 : (vel_y < 0) ? -1.0 : 0.0;

  //const double bed_fric_contr = is_bed_friction ? vel_y_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;
  const double bed_fric_contr_one = (h*h)>epsilon && is_bed_friction && hs>epsilon ? density*Uy*grav*std::abs(Uy)/turbulence_coeff/h/h : 0.;
  const double bed_fric_contr_two = is_bed_friction && hs>epsilon ? vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;


  return ( - bed_fric_contr_one - bed_fric_contr_two + h*R_y );
}


double
TG2_scheme::Uxs_src_formula_2 (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;
  const double bed_pressure = grav*hs*(1.-r_coeff) - bed_excess_pore_water_pressure/density_s; 
  const auto Ux = Uxw + Uxs;

  const double density = ns + n*r_coeff;

  const double vel_x_sign = 2./M_PI*std::atan(M_PI*.5*Ux); //std::abs(Ux)>tolerance_sign ? Ux/std::abs(Ux) : Ux/tolerance_sign;

  const double bed_fric_contr_one = h*h>epsilon && is_bed_friction && hs>epsilon ? density*Ux*grav*std::abs(Ux)/turbulence_coeff/h/h : 0.; //is_bed_friction ? vel_x*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction && hs>epsilon ? vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  return ( - bed_fric_contr_one - bed_fric_contr_two );
}

double
TG2_scheme::Uys_src_formula_2 (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys, const double& bed_excess_pore_water_pressure)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;
  const double bed_pressure = grav*hs*(1.-r_coeff) - bed_excess_pore_water_pressure/density_s;
  const auto Uy = Uyw + Uys;


  const double density = ns + n*r_coeff;


  const double vel_y_sign = 2./M_PI*std::atan(M_PI*.5*Uy); //std::abs(Uy)>tolerance_sign ? Uy/std::abs(Uy) : Uy/tolerance_sign;

  const double bed_fric_contr_one = (h*h)>epsilon && is_bed_friction && hs>epsilon ? density*Uy*grav*std::abs(Uy)/turbulence_coeff/h/h : 0.; 
  const double bed_fric_contr_two = is_bed_friction && hs>epsilon ? vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  return ( - bed_fric_contr_one - bed_fric_contr_two );
}


double
TG2_scheme::src_slope_formula (const double& h, const double& S_x, const double& S_y, const int& kk)
{
  // nodal slope term for the second step in the main

  double S = 0.;

  if (kk%3==1) // Ux
  {
    S = S_x;
  }
  else if (kk%3==2) // Uy
  {
    S = S_y;
  }

  //std::cout << S << " " << kk << " " << kk%3 << std::endl;

  return (-grav*S*h);

}


double 
TG2_scheme::src_slope_formula (const double& h, const double& S)
{
  // cell-wise source term to build the incr vector
  return (-grav*S*h);
}


double
TG2_scheme::signum (const double& x)
{ return ((x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0); }


void
TG2_scheme::resize_vectors ()
{
  contr_x.resize(number_FD_points);
  contr_y.resize(number_FD_points);

  alfa_vec.resize(number_FD_points-1);
  beta_vec.resize(number_FD_points-1);
     y_vec.resize(number_FD_points-1);
}

void
TG2_scheme::set_r_coeff()
{
  r_coeff = density_w/density_s;
}









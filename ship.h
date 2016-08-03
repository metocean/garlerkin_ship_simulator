/*  The Galerkin Ship Simulator
 *  Copyright (C) 2014 MetOcean Solutions Limited
 *
 *  This file is part of the Galerkin Ship Simulator.        
 *
 *  The Galerkin Ship Simulator is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The Galerkin Ship Simulator is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the Galerkin Ship Simulator.  If not, see <http://www.gnu.org/licenses/>.
*/

typedef gpointer     (* LUFactoriseFunc)              (gsl_matrix * A);
typedef void         (* LUFactorisedSolveFunc)        (gpointer p,
						       gsl_vector * rhs);
typedef void         (* LUDestroyFunc)                (gpointer p);

typedef void         (* SchemeFunc)                   (Simulation * sim,
						       gdouble t,
						       gboolean prediction);
typedef gsl_vector * (* RHSFunc)                      (Simulation * sim, 
						       gdouble t,
						       gsl_vector * rhs);
typedef BoundarySubProblem *  (* BuildSubproblemFunc) (Spline2D * sp1,
						       Spline2D * sp2,
						       SelfInfluenceFunc self_influence_coeffs);

typedef struct {
  gdouble u, v;
  gboolean is_wet;
  Point p;
} CornerPoint;

typedef CornerPoint * (* CornerPointNewFunc)          (Spline2D * sp, Hull * h, 
						       gdouble u, gdouble v,
						       HeightCurve hz, double t,
						       gpointer data);
typedef void         (* PanelIntegrationFunc)         (SPPanel * spp,
						       Simulation * sim,
						       gpointer data,
						       HeightCurve hz,
						       gdouble t,
						       gpointer hz_data);

typedef void         (* LocalGaussIntegrationFunc)    (SPPanel * spp,
						       Simulation * sim,
						       gpointer data,
						       gint m, gint n,
						       HeightCurve hz,
						       gdouble t,
						       gpointer hz_data);

typedef gboolean     (* ToleranceIntegrationFunc)     (CornerPoint * p00,
						       CornerPoint * p10,
						       CornerPoint * p11,
						       CornerPoint * p01,
						       gpointer data);
typedef void         (* LocalIntegrationFunc)         (Spline2D * sp,
						       Simulation * sim,
						       gpointer data,
						       gdouble u, gdouble v,
						       gdouble weight,
						       HeightCurve hz,
						       gdouble t,
						       gpointer hz_data);

typedef struct {
  gdouble a, b, c, d, e;
} Poly4;

typedef struct {
  double xdir, ship_heading;
  gboolean rotated;
  int nwe, ndir;
  double * we, * dir;
  double * amp, * phase;
} Spectrum;

typedef struct {
  int nwe, ndir;
  double * we, * dir;
  double * re[6], * im[6];

  interp2d_spline * interp_re[6], * interp_im[6];
  gsl_interp_accel * xa, * ya;
} DiffractionCoeffs;

struct _Simulation {
  Hull * hull;
  FreeSurface * fs;
  Bathymetry * bathy;
  gint N, M;
  gdouble h, t, w, A;

  // Background translation velocity
  Vector U;
  // Gravity
  gdouble g;
  // Density
  gdouble rho;
  // Time parameters for the simulation
  Time time;
  // Length of boat/object
  gdouble lpp;
  // Sub boundary problems
  GPtrArray * sub_problems;
  // Neumann
  BoundaryProblem * neumann_problem;
  BoundaryProblem * dirichlet_problem;
  BoundaryProblem * mixed_problem;
  BoundaryProblem * nospeed_problem;
  // LU problems
  gpointer neumann_lu;
  gpointer dirichlet_lu;
  gpointer mixed_lu;
  gpointer nospeed_lu;
  // LU methods
  LUFactoriseFunc lu_factorise;
  LUFactorisedSolveFunc lu_factorised_solve;
  LUDestroyFunc lu_destroy;
  // Parameters of the wave forcing
  WaveParams wp;
  // Boundary problem method
  BuildSubproblemFunc build_boundary_subproblem;
  SelfInfluenceFunc self_influence_function;
  // Numerical beaches
  XYZFunc numerical_beaches;
  // Free-surface elevation advection scheme for first order potential
  SchemeFunc fs_elevation_update;
  SchemeFunc fs_elevation_rhs_store;
  gsl_vector * fseu_rhs1, * fseu_rhs2, * fseu_rhs3, * fseu_rhs4; // Storage for 4th order abm
  // Free-surface potential advection scheme for first order potential
  SchemeFunc fs_potential_update;
  SchemeFunc fs_potential_rhs_store;
  gdouble beta;
  gsl_vector * fspu_rhs1, * fspu_rhs2, * fspu_rhs3, * fspu_rhs4; // Storage for 4th order abm

  // Free-surface elevation advection scheme for first order potential
  SchemeFunc fs_elevation2_update;
  SchemeFunc fs_elevation2_rhs_store;
  gsl_vector * fseu2_rhs1, * fseu2_rhs2, * fseu2_rhs3, * fseu2_rhs4; // Storage for 4th order abm
  // Free-surface potential advection scheme for first order potential
  SchemeFunc fs_potential2_update;
  SchemeFunc fs_potential2_rhs_store;
  gsl_vector * fspu2_rhs1, * fspu2_rhs2, * fspu2_rhs3, * fspu2_rhs4; // Storage for 4th order abm

  // Forces
  //ForcesHistory * fh;
  // External forces (function pointer)
  GSList * forces;
  GSList * forces_2;
  // Inversed mass matrix
  gpointer mass_lu;
  // Netcdf forcing
  NetCDFForcing * ncdf;
  // Verbose
  gpointer verbose;
  gboolean continuity;
  // Moore masters
  GPtrArray * MM;
  gdouble force_coeff;
  // MooringLines
  GPtrArray * MooringLines;
  Poly4 mooring_forces[6][6];
  // Fenders
  GPtrArray * fenders;
  // Mode of motion for forced radiation problem
  gint d;
  gdouble forcing_scaling_factor;
  Spline2D * tmp;
  // Spectral forcing
  Spectrum * spec;
  DiffractionCoeffs * idc;
};





Spline2D * periodic_sphere              (gint M, gint N, gint k,
					 gint ninner, gint nouter,
					 gdouble r);
void       initialise_motion            (Simulation * sim);
void       solve_equation_of_motion     (Hull * hull,
					 Time * time,
					 Simulation * sim);
void       calculate_added_mass_matrix  (Simulation * sim, Point xg);
typedef    void      (* ForceFunc)      (Simulation * sim,
					 Forces * f, gdouble t,
					 gdouble u[6], gdouble x[6],
					 gboolean prediction);
typedef    void      (* ForceFunc1)     (Simulation * sim,
					 Forces * f, gdouble t,
					 Motion m,
					 gboolean prediction);
void       find_equilibrium_position_z  (Hull * hull,
					 Simulation * sim);
void       adjust_linear_equilibrium_position_xgx (Hull * hull,
						   Simulation * sim);


Simulation * simulation_new ();
void simulation_destroy (Simulation * sim);
void simulation_build_problems_nospeed_implicit (Simulation * sim);
void simulation_set_to_zero (Simulation * sim);
GSList * simulation_all_patches_list (Simulation * sim);
void simulation_build_galerkin_fit_matrixes (Simulation * sim);
void simulation_build_problems (Simulation * sim);

void solve_netcdf_boundary_problem_for_disturbance_flow_kim_1 (Simulation * sim);
void solve_boundary_problem_for_disturbance_flow_1 (Simulation * sim);
void solve_boundary_problem_for_disturbance_flow_nospeed_kim (Simulation * sim);
void solve_boundary_problem_for_basis_flow (Simulation * sim);
void solve_boundary_problem_for_disturbance_flow_kim (Simulation * sim);

void build_free_surface (Simulation * sim,
			 gint k, gint ninner, gint nouter);


gdouble numerical_beaches_radiation (gdouble x, gdouble y, gdouble z,
				     gpointer data);
gdouble trial_numerical_beaches_implicit (SPPanel * spp, gint m, gint n, gpointer data);

void semi_implicit_no_speed_potential_update (Simulation * sim, gdouble t, gboolean prediction);
void leapfrog_no_speed_elevation_update (Simulation * sim, gdouble t, gboolean prediction);
void abm4_no_speed_elevation_update (Simulation * sim, gdouble t, gboolean prediction);
void abm4_no_speed_elevation_store (Simulation * sim, gdouble t, gboolean prediction);
void abm4_no_speed_potential_update (Simulation * sim, gdouble t, gboolean prediction);
void abm4_no_speed_potential_store (Simulation * sim, gdouble t, gboolean prediction);

gdouble flat_sea (gdouble x, gdouble y, gdouble t, gpointer data);

void spline2d_filter_variable (Spline2D * sp, gint var);
void periodic_fs_filter_variable (Spline2D * sp, gint var);

gdouble calculate_mass_from_position_at_rest (Hull * hull,
					      Simulation * sim);
gdouble calculate_mass_from_position_at_rest_linear (Hull * hull,
						     Simulation * sim);
void calculate_hydrostatic_restoring_coeffs (Hull * hull,
					     Simulation * sim,
					     HeightCurve hz, gdouble t,
					     gpointer hz_data);

void solve_equation_of_motion_1 (Hull * hull,
				 Time * time,
				 Simulation * sim,
				 GSList * forces);

Forces forces_evaluation (Hull * hull,
			  gdouble t,
			  Time * time,
			  Motion m,
			  Simulation * sim,
			  GSList * forces,
			  gboolean prediction);
Forces * wet_hull_pressure_force_integration (Hull * hull,
					      Simulation * sim,
					      HeightCurve hz, gdouble t,
					      gpointer hz_data);

Point find_center_of_buoyancy (Hull * hull,
			       Simulation * sim,
			       HeightCurve hz, gdouble t,
			       gpointer hz_data);

void adjust_equilibrium_position_xgx (Hull * hull, Simulation * sim);

void add_netcdf_fk_fh_forces (Simulation * sim, Forces * f, gdouble t,
			      gdouble u[6], gdouble x[6]);
void add_hydrostatic_restoring_forces (Simulation * sim, Forces * f, gdouble t, gdouble u[6], gdouble x[6]);
void add_full_netcdf_fk_fh_forces (Simulation * sim,
				   Forces * f, gdouble t,
				   gdouble u[6], gdouble x[6],
				   gboolean prediction);
void add_roll_damping_force (Simulation * sim, Forces * f, gdouble t,
			     gdouble u[6], gdouble x[6]);
void add_full_netcdf_fk_fh_forces_implicit (Simulation * sim,
					    Forces * f, gdouble t,
					    gdouble u[6], gdouble x[6]);
void add_gravity_force (Simulation * sim, Forces * f, gdouble t,
			gdouble u[6], gdouble x[6]);
void add_freqdomain_damping_general_1 (Simulation * sim,
				       Forces * f,
				       gdouble t,
				       Motion m,
				       gboolean prediction);
void add_pdstrip_netcdf_fk_fh_forces_1 (Simulation * sim,
					Forces * f,
					gdouble t,
					Motion m,
					gboolean prediction);
void add_hydrostatic_restoring_force_1 (Simulation * sim,
					Forces * f,
					gdouble t,
					Motion m,
					gboolean prediction);
void add_roll_damping_force_1  (Simulation * sim,
				Forces * f,
				gdouble t,
				Motion m,
				gboolean prediction);
void add_full_netcdf_fk_fh_forces_1 (Simulation * sim,
				     Forces * f,
				     gdouble t,
				     Motion m,
				     gboolean prediction);
void add_diffraction_forces_from_spectra_1 (Simulation * sim,
					    Forces * f,
					    gdouble t,
					    Motion m,
					    gboolean prediction);
void add_pdstrip_initial_damping_1 (Simulation * sim,
				    Forces * f,
				    gdouble t,
				    Motion m,
				    gboolean prediction);




Spectrum * spectrum_new (int nwe, int ndir);
void spectrum_rotate_to_ship_heading (Spectrum * spec, double heading);
void spectrum_diwasp_generate_time_series (Spectrum * spec);
void spectrum_array_print (Spectrum * spec, char * out_file_name);
Spectrum * spectrum_diwasp_from_file (char * file_name, gdouble scaling);


void calculate_nemoh_radiation_matrix (Simulation * sim, gchar * f_a_name, gchar * f_b_name);
void calculate_pdstrip_radiation_matrix (Simulation * sim, FILE * fam);

DiffractionCoeffs * diffraction_coeffs_new (int nwe, int ndir);
DiffractionCoeffs * diffraction_coeffs_from_pdstrip_file (gchar * file_name);
DiffractionCoeffs * diffraction_coeffs_from_nemoh_file (gchar * file_name);
DiffractionCoeffs * diffraction_coeffs_interpolated (DiffractionCoeffs * dc, double * we,
						     double * dir, int nwe, int ndir);
void diffraction_coeffs_print (DiffractionCoeffs * dc, gchar * out_file_name, int k);
void diffraction_coeffs_destroy (DiffractionCoeffs * dc);


GSList * sppanel_point_location (SPPanel * spp,
				 gdouble x, gdouble z,
				 Hull * h);

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

typedef struct {
  Vector v, omega;
  Point p;
} KinematicTorseur;

typedef struct {
  
} DynamicTorseur;

typedef struct {
  gdouble fv, fx;
  Vector Fm, Fl, Fh;
  Vector Mm, Ml, Mh;
  // Memory, local, hydrostatic
  gdouble forces_m[6], forces_l[6], forces_h[6];
  // Froude-Krylov, Non-linear restoring
  gdouble forces_fk[6], forces_r[6];
  // Entrainment forces
  gdouble forces_e[6];
  // External forces
  gdouble forces_ext[6];
  gdouble v[6];
  gdouble u[6];
  gdouble x[6];
  gdouble phi1[6], phi2[6], a[6];
  Vector gradphi1[6];

  gdouble x1[6], x2[6];
  gdouble u1[6], u2[6];
  gdouble v1[6], v2[6];
  gdouble forces_1[6], forces_2[2];
  gdouble forces_h1[6], forces_h2[6];
  gdouble forces_fk1[6], forces_fk2[6];
  gdouble forces_ext1[6], forces_ext2[6];
  gdouble forces_v1[6], forces_v2[6];
  // RHS of fs potential and elevation
  //gsl_vector * fse_rhs, * fsp_rhs;
} Forces;

Forces   forces_sum                           (Forces f1, Forces f2);
void     forces_set_to_zero                   (Forces * f);
void     forces_print                         (Forces * f, FILE * fp);
void     forces_total_print                   (Forces * f, FILE * fp, Time * t);
void     forces_add_contribution_forces       (Forces * f, Vector * d,
					       gdouble wm, gdouble wl, gdouble wh,
					       gdouble wfk, gdouble wr);
void     forces_add_contribution_moments      (Forces * f, Vector * d,
					       gdouble wm, gdouble wl, gdouble wh,
					       gdouble wfk, gdouble wr);

void     test_rk                              ();

void     test_abm                             ();

void     add_force_to_history                 (ForcesHistory * fh,
					       Forces f);
void     forces_history_destroy               (ForcesHistory * fh);

Motion   motion_ABM4_time_integration         (gdouble t, gdouble dt,
					       Motion m0, ForcesHistory * fh,
					       Forces f1);
void     update_rotation_matrix               (Motion * m,
					       gdouble x[6],
					       gdouble u[6]);

void     solve_equation_of_motion_pure_RK4    (Hull * hull,
					       Time * time,
					       Simulation * sim,
					       Forces * f1,
					       GSList * forces);
gdouble  hull_motion_time_local_bc_huang      (SPPanel * spp,
					       gint m, gint n,
					       gpointer data);

/* gdouble  hull_motion_time_local_bc_linear_kim (SPPanel * spp, */
/* 					       gint m, gint n, */
/* 					       gpointer data); */
Forces * solve_equation_of_motion_RK4_ABM4    (Hull * hull,
					       Time * time,
					       Simulation * sim,
					       Forces * f1,
					       GSList * forces);
/*****************************************************/
/*                 CLEAN FORMULATION                 */
/*****************************************************/

void      motion_update_rotation_matrix_1     (Motion * m,
					       gdouble x1[6]);
void      motion_update_rotation_matrix_2     (Motion * m,
					       gdouble x1[6],
					       gdouble x2[6]);

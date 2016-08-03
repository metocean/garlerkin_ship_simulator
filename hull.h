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
  GSList * f;
  
} ForcesHistory;

typedef struct {
  gdouble x0;
  gdouble v0;

  gdouble x[6];
  gdouble v[6];
  gdouble u[6];

  /* Clean */
  Matrix3 Rbi1, Rbi2;
  Matrix3 Rib1, Rib2;
  gdouble x1[6], x2[6];
  gdouble u1[6], u2[6];

  Point equilibrium;
  Vector t;
  Matrix3 euler_m;
  Matrix3 euler_r;
} Motion;

void     motion_update_euler_matrix       (Motion * m);
Point    motion_transform_point           (Point p, Motion * m);
Vector   motion_transform_vector          (Vector v,
					   Motion * m);

typedef struct {
  GSList * patches;
  GSList * wet_patches;
  GSList * hr_patches;
  /* Transformation t; */
  ForcesHistory * fh;
  Point xg; // Center of mass
  gdouble mg; // mass
  gdouble Ig[3][3]; // Matrix of interia
  //gdouble B[6][6]; // Change of coordinate system matrix
  gdouble M[6][6]; // Mass matrix
  gdouble A[6][6]; // Added mass matrix
  gdouble D[6][6]; // Damping coefficients matrix
  gdouble R[6][6]; // Restoring coefficients matrix
  Motion m;
  Vector U;

  /* Clean */
  gpointer mass_lu1; // Mass matrix with added mass
  gpointer mass_lu2; // Mass matrix without added mass
  /* ForcesHistory * fh1; */
  /* ForcesHistory * fh2; */
  gdouble M1[6][6], M2[6][6];
  gdouble diff[6][2];
  gdouble we;

  // For pdstrip
  gdouble * RET;
  int nt;
  GList * u_history;

  // Hull characterisitic dimensions
  double lbp, beam, max_draft, raise, vessel_height, deck_height;
} Hull;

GPtrArray * sort_panels                      (GArray * all_panels);

Hull   *  hull_new                           ();
void      hull_destroy                       (Hull * h);
void      hull_read                          (Hull * h, FILE * fp,
					      gint M, gint N,
					      gboolean flip_u,
					      gboolean flip_v,
					      gboolean centripetal_reparam,
					      gboolean swap);
void      hull_read_old                      (Hull * h, FILE * fp,
					      gint M, gint N, gboolean flip,
					      gboolean centripetal_reparam);
void      hull_print                         (Hull * h, FILE * fout);
void      hull_print_gnuplot                 (Hull * h, FILE * fout,
					      gint var, gdouble t);
void      hull_motion_print                  (Hull * h, FILE * fp);
DCurve *  hull_intersect_with_free_surface   (Hull * h, HeightCurve hz,
					      gdouble t, gpointer data,
					      gint N);
gdouble   hull_normal_velocity               (SPPanel * spp,
					      gint m, gint n,
					      gpointer data);
Vector    hull_normal_time_derivative        (Hull * hull,
					       gdouble u, gdouble v);
Point     spline2d_hull_eval_point           (Spline2D * sp,
					      Hull * h,
					      gdouble u, gdouble v);
Point     spline2d_hull_eval_point_gauss_point (Spline2D * sp,
						GaussPoints * gp,
						Hull * h,
						gint m, gint n);
Point     hull_transformed_point              (Hull * h, Point p);
Vector    hull_transformed_vector             (Hull * h, Vector v);
Point     hull_transformed_eval_point         (Hull * h, Spline2D * sp,
					       gdouble u, gdouble v);
Vector    hull_velocity_at_point              (Hull * h, Point p);
Vector    hull_velocity_at_point_1            (Hull * h, Point p);
Spline2D * spline2d_wet_part                  (Spline2D * sp,
					       Hull * hull,
					       HeightCurve hz,
					       gdouble t,
					       gpointer data);
void      hull_print_wet                      (Hull * h,
					       FILE * fout);
void hull_generate_wet_hull (Hull * hull, HeightCurve hz,
			     gdouble t, gpointer data);
/*****************************************************/
/*                 CLEAN FORMULATION                 */
/*****************************************************/

Vector    hull_normal_1                       (Hull * hull,
					       Vector * n0);
Vector    hull_normal_2                       (Hull * hull,
					       Vector * n0);
void      hull_initialise_motion              (Hull * h);
void      hull_print_mayavi                   (Hull * h,
					       FILE * fout);

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <string.h>
#include <complex.h>
#include <sys/resource.h>

#include <getopt.h>
#include <unistd.h>
#define _GNU_SOURCE
#include <sched.h>

/* Glib */
#include <glib.h>

/* GSL library for spline fitting */
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>

/* Some compiling options */
#define CUDA FALSE
#define OPENCL FALSE
#define PLASMA TRUE

#define DEBUG FALSE
#define LINEAR TRUE
#define NOFLUX FALSE

#if CUDA
/* Magma GPU dense system solver library */
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <magma.h>
#endif

#if OPENCL
#include <clmagma.h>
#endif

/* Multi-thread dense system solver library */
#include <plasma.h>

/* SuperLU sparse system direct solver library */
#include <slu_ddefs.h>

/* NetCDF */
#include <netcdf.h>

/* 2D interpolation https://github.com/diazona/interp2d/ */
#include <interp2d.h>
#include <interp2d_spline.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))
#define max(x1,x2) ((x1) > (x2) ? (x1):(x2))

gint     sign                              (gdouble x);

typedef struct {
  gdouble x, y, z, xi;
} Point;

typedef struct {
  gdouble x, y, z;
} Vector;

gdouble  vector_norm                       (Vector v);
Vector   vector_normalise                  (Vector v);
gdouble  vector_scalar_product             (Vector * v0, Vector * v1);
gdouble  vector_normalized_scalar_product  (Vector * v0, Vector * v1);
Vector   vector_vector_product             (Vector * v0, Vector * v1);
Vector   vector_sum                        (Vector v0, Vector v1);
Vector   vector_times_constant             (Vector v, gdouble c);
void     vector_set_to_zero                (Vector * v);

typedef struct {
  gdouble x[6];
} Vector6;

void     vector6_set_to_zero               (Vector6 * v);
Vector6  vector6_sum                       (Vector6 v0, Vector6 v1);
Vector6  vector6_times_constant            (Vector6 v, gdouble c);
void     vector6_add                       (Vector6 * v0, Vector6 * v1);
void     vector_multiply_by_constant       (Vector6 * v, gdouble c);
typedef struct {
  gdouble a[3][3];
} Matrix3;

Vector   vector_rotate                     (Matrix3 * m, Vector v0);
void     matrix3_inverse                   (Matrix3 * m);

typedef struct {
  Point xg;
  Vector t;
  Vector euler;
  Matrix3 euler_m;
} Transformation;

void     transformation_euler_matrix       (Transformation * t);
Point    transform_point                   (Point p, Point * xg,
					    Vector * t, Matrix3 * euler_m);
Vector   transform_vector                  (Vector v, Matrix3 * euler_m);

void     point_print                       (Point p, FILE * fp,
					    Transformation * t);
gdouble  point_distance                    (Point p1, Point p2);
gdouble  point_distance_squared            (Point p1, Point p2);
typedef gdouble (* HeightCurve)          (gdouble x, gdouble y, gdouble t,
					  gpointer data);

typedef struct {
  gdouble start, end;
  gdouble t, dt;
  gint iend, itime;
} Time;

typedef struct {
  Point p[4];
  gdouble var [12];
  gint i, j;
  gdouble h[2];
  gboolean border;
} Panel;

Panel *  panel_new                ();
void     panel_destroy            (Panel * p);
gboolean panel_read               (Panel * p, FILE * fp);
Panel *  panel_new_read           (FILE * fp);
Vector   panel_first_order_normal (Panel * p);
gdouble  panel_eval               (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_y             (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_x             (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dx            (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dx_x          (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dy            (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dy_y          (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dxdx          (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dxx           (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dxx_x         (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dyy           (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dyy_y         (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_dxy           (Panel * p, gdouble i,
				   gdouble j, gint var);
#define PANEL_VAL(p,i)     (p->var[i])

gdouble  panel_eval_int           (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_int_y         (Panel * p, gdouble i,
				   gdouble j, gint var);
gdouble  panel_eval_int_x         (Panel * p, gdouble i,
				   gdouble j, gint var);

typedef struct {
  gint ncoeffs;
  gsl_vector * bi, * ci;
  gsl_matrix * bdi, * cov;
  gsl_bspline_workspace * bw;
  gsl_bspline_deriv_workspace * bdw;
} BSpline;

/* Contains a discrete representation of the boudary curve */
typedef struct {
  GArray * p;
} DCurve;

/* Contains a spline representation of the boundary curve */
typedef struct {
  BSpline * bsx, * bsy, * bsz;
} BCurve;

DCurve *      dcurve_new                       (gint N);
void          dcurve_destroy                   (DCurve * c);
void          dcurve_print                     (DCurve * dc, FILE * fp);
DCurve *      resample_point_list_centripetal  (GSList * l, gint m);
DCurve *      resample_point_list_hybrid       (GSList * l, gint m);

BSpline *     bspline_new                      ();
void          bspline_destroy                  (BSpline * bs);
gdouble       bspline_eval                     (BSpline * bs, gdouble xi);
gdouble       bspline_eval_first_derivative    (BSpline * bs, gdouble xi);

Point         bcurve_eval                      (gdouble u, BCurve * bc);
BCurve *      bcurve_fit_point_list            (GSList * l);

typedef struct {
  GArray * matrix;
  GArray * column;
  GArray * index;

  // For superlu (lu decomposition storage)
  gint * perm_r, * perm_c;
  SuperMatrix L, U;
  SuperLUStat_t stat;
  gint info;
  gboolean factorised;
} CCSProblem;

typedef struct {
  GArray * ui, * vj;
  GArray * Pi, * Ni, * Ji;
  GArray * ustart, * vstart;
  GPtrArray * Bu, * Bv;
  GArray * fsdata;
} GrevillePoints;

GrevillePoints * greville_points_new           (gint M, gint N);
void             greville_points_destroy       (GrevillePoints * gp);

typedef struct {
  gdouble v[45];
} SplineCoeffs;

typedef struct _Spline2D Spline2D;
typedef struct _SPPanel SPPanel;
typedef struct _GaussPoints GaussPoints;

typedef Vector          (* NeumannFunc)                    (Spline2D * sp,
							    gdouble u, gdouble v,
							    gpointer data);
typedef gdouble         (* DirichletFunc)                  (Spline2D * sp,
							    gdouble u, gdouble v,
							    gpointer data);
typedef gdouble         (* Spline2DFunc)                   (Spline2D * sp,
							    gdouble u, gdouble v,
							    gpointer data);
typedef gdouble         (* GrevilleFunc)                   (Spline2D * sp,
							    gint m, gint n,
							    gpointer data);
typedef gdouble         (* XYZFunc)                        (gdouble x,
							    gdouble y,
							    gdouble z,
							    gpointer data);
typedef gdouble         (* SPPanelFunc)                     (SPPanel * spp,
							    gdouble u, gdouble v,
							    gpointer data);

typedef gdouble         (* GaussFunc)                      (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
typedef CCSProblem *    (* BuildFitMatrixFunc)             (Spline2D * sp);
typedef gsl_vector *    (* BuildFitRhsFunc)                (Spline2D * sp,
						            GaussFunc Func,
							    gpointer data,
							    Spline2DFunc bc_func,
							    gpointer bc_data,
							    gsl_vector * rhs);
typedef void            (* CopyFitSolutionFunc)            (Spline2D * sp,
							    gsl_vector * lhs,
							    gint var);
typedef void            (* AddFitSolutionFunc)             (Spline2D * sp,
							    gsl_vector * lhs,
							    gint var,
							    gint var0);
typedef void            (* ReinitPanelsFunc)               (Spline2D * sp);
typedef gint            (* SizeFunc)                       (Spline2D * sp);
typedef gsl_vector *    (* RhsAllocFunc)                   (Spline2D * sp,
							    gsl_vector * rhs);

struct _Spline2D {
  // Spline order for geometry, for variables
  gint k;
  // number of planel in each directions
  gint M, N;
  // number of spline coefficients in each directions
  gint NXU, NXV, NUT; // NUT = total in u direction for free-surface
  gint NU, NV;
  // Is it a periodic spline?
  gboolean periodic;
  gboolean noflux;
  // Splines
  gsl_bspline_workspace * wx_u;
  gsl_bspline_workspace * w_u, * w_v;
  // Derivatives
  gsl_bspline_deriv_workspace * wxd_u;
  gsl_bspline_deriv_workspace * wd_u, * wd_v;
  // Splines coeffs
  GPtrArray * coeffs;
  GPtrArray * panels;
  // Need to add the rotation/translation matrix
  gint ninner, nouter;
  //Transformation * t;
  gint istart;
  // Fitting matrix
  CCSProblem * fit;
  CCSProblem * fit_noflux;
  CCSProblem * fit_greville;
  // Fitting methods
  BuildFitMatrixFunc  build_fit_matrix;
  BuildFitRhsFunc     build_fit_rhs;
  BuildFitMatrixFunc  build_fit_noflux_matrix;
  BuildFitRhsFunc     build_fit_noflux_rhs;
  CopyFitSolutionFunc copy_fit_solution;
  AddFitSolutionFunc  add_fit_solution;
  ReinitPanelsFunc    reinit_panels;
  SizeFunc            size;
  RhsAllocFunc        rhs_vector;
  // Greville Points for collocation method
  GrevillePoints * gr;
  // 
  Spline2D * next;
  gint fs_index, fs_index_x;
  gboolean fully_submerged;
  Spline2D * hull_patch;
  gsl_vector * rhs;
  CCSProblem * fit_test;
  CCSProblem * rhs_test;
  CCSProblem * fit_tmp;
};

typedef struct {
  GArray * x, * X;
  gsl_matrix * x2X, * X2x;
  size_t istart, jstart;
  GArray * w, * n, * om;
  gdouble ue, ve;
} S2P;

typedef struct {
#if LINEAR
  //Vector gradPhi;
  gdouble dzdzphi;
  //#else
  Vector gradphi, gradPhi, gradPhi2, gradphi0, gradzeta0, gradzeta;
  Vector graddzPhi, graddzphi0, graddzphi;
  gdouble zeta, zeta0, dtzeta0, dtPhi, dtphi0, dtphi;
  gdouble dtdzPhi, dtdzphi, dtdzphi0;
  /* gdouble phidt, Phidt, phi0dt, zetadt, zeta0dt; */
#endif
} FSData;

typedef struct {
  gint ix1, ix2;
  gint iy1, iy2;
  float w11, w12, w21, w22;
} Interpol;

struct _GaussPoints {
  GArray * ui, * vj;
  GArray * Pi, * Ni, * wJij, * wij;
  GArray * c1, * c2, * c3, * c4, * c5, * c6;
  GPtrArray * Bu, * Bv;
  GPtrArray * Bux;
  gint istart, jstart;
  gint istart_x;
  GArray * fsdata;
  GPtrArray * interpol;
};

GaussPoints * gausspoints_new                     (gint size);
void          gausspoints_destroy                 (GaussPoints * gp);
Point         spline2d_eval_gauss_point_point     (Spline2D * s,
						   GaussPoints * gp,
						   gint m, gint n);
gdouble       spline2d_eval_gauss_point           (Spline2D * s,
						   GaussPoints * gp,
						   gint m, gint n,
						   gint var);

gdouble       spline2d_eval_greville_point        (Spline2D * s,
						   GrevillePoints * gr,
						   gint m, gint n,
						   gint var);
void          spline2d_store_greville_points_data (Spline2D * sp);

Spline2D *    spline2d_symmetrical_y              (Spline2D * sp,
						   gdouble y);
void          spline2d_translate                  (Spline2D * sp,
						   gdouble x,
						   gdouble y,
						   gdouble z);
void          spline2d_rescale                    (Spline2D * sp,
						   gdouble sx,
						   gdouble sy,
						   gdouble sz);
void          spline2d_rotate                     (Spline2D * sp,
						   Point xc,
						   gdouble angle);

typedef struct {
  // Source influence coefficients
  gsl_matrix * psi;
  // Dipole influence coefficients
  gsl_matrix * phi;
} InfluenceCoeffs;

/* typedef struct { */
/*   Vector gradphi, gradPhi, gradPhi2, gradphi0, gradeta; */
/* } Gradients; */

struct _SPPanel {
  Spline2D * sp;
  gint k;
  gdouble u0, u1, v0, v1;
  gdouble ue, ve;
  Point pe;
  GaussPoints * outer;
};

SPPanel *       sppanel_new                           (gint k);
void            sppanel_destroy                       (SPPanel * spp);
GaussPoints *   sppanel_store_gauss_legendre_points   (SPPanel * spp, gint ng);
void            sppanel_store_gauss_legendre_data     (GaussPoints * gp,
						       SPPanel * spp);

gdouble      coeff                                       (Spline2D * s, gint i,
							  gint j, gint var);
void         coeff_assign                                (Spline2D * s, gint i,
							  gint j, gint var,
							  gdouble coeff);
void         coeff_set_var_to_zero                        (Spline2D * s, gint var);
void         coeff_set_var_to_one                         (Spline2D * s, gint var);
void         coeff_set_var_to_constant                    (Spline2D * s, gint var, gdouble val);
gsl_vector * spline2d_rhs_vector                          (Spline2D * sp,
							   gsl_vector * rhs);
Spline2D *   spline2d_new                                 (gint N, gint M, gint k,
							   gint ninner, gint nouter);
void         spline2d_destroy                             (Spline2D * s);
gdouble      spline2d_eval                                (Spline2D * s, gdouble u,
							   gdouble v, gint var);
Point        spline2d_eval_point                          (Spline2D * s, gdouble u,
							   gdouble v);
gdouble      spline2d_eval_spline                         (Spline2D * s, gint i, gint j,
							   gdouble u, gdouble v);
gdouble      spline2d_derivative_eval                     (Spline2D * s, gdouble u,
							   gdouble v, gint m,
							   gint n, gint var);
gdouble      spline2d_derivative_eval_gauss_point         (Spline2D * s,
							   GaussPoints * gp,
							   gint m, gint n,
							   gint dm, gint dn,
							   gint var);
gdouble      spline2d_derivative_eval_greville_point      (Spline2D * s,
							   GrevillePoints * gr,
							   gint m, gint n,
							   gint dm, gint dn,
							   gint var);
void         spline2d_copy_var                            (Spline2D * sp,
							   gint var1, gint var2);
void         spline2d_list_copy_var                       (GSList * list,
							   gint var1, gint var2);
void         spline2d_fit_geometry                        (Spline2D * s);
gdouble      spline2d_integration                         (Spline2D * s, gint var);
gdouble      spline2d_gauss_integration                   (Spline2D * s,
							   GaussFunc func,
							   gpointer data);
gdouble      spline2d_adaptive_gauss_integration          (Spline2D * sp,
							   SPPanelFunc func,
							   gpointer data,
							   gint level);
void         spline2d_init_panels                         (Spline2D * sp);
void         spline2d_reinit_panels_physical_quantities   (Spline2D * sp);
void         sppanel_print                                (SPPanel * spp, FILE * fp);
void         spline2d_print_panels                        (Spline2D * sp, FILE * fp);
void         spline2d_print_transformed_panels            (Spline2D * sp, FILE * fp,
							   Point * xg, Vector * t,
							   Matrix3 * euler_m);
void         spline2d_print_panels_gnuplot                (Spline2D * sp, FILE * fp,
							   gint var);
void         spline2d_print_normals                       (Spline2D * sp, FILE * fp);
void         spline2d_print_transformed_normals           (Spline2D * sp, FILE * fp,
							   Point * xg, Vector * t,
							   Matrix3 * euler_m);
GSList *     spline2d_intersect_with_free_surface         (Spline2D * patch,
							   HeightCurve hz,
							   gdouble t, gpointer data);
/* Spline2D *   spline2d_wet_part                            (Spline2D * sp, HeightCurve hz, */
/* 							   gdouble t, gpointer data); */
/* gboolean     sppanel_is_wet                               (Spline2D * sp, */
/* 							   gdouble u0, gdouble u1, */
/* 							   gdouble v0, gdouble v1, */
/* 							   HeightCurve hz, */
/* 							   gdouble t, gpointer data, */
/* 							   gboolean * intersect); */
gboolean     sppanel_intersects_curve                     (SPPanel * spp, HeightCurve hz,
							   gdouble t, gpointer data);
double       spline2d_jacobian                            (Spline2D * sp, gdouble u,
							   gdouble v);
Vector       spline2d_normal                              (Spline2D * sp, gdouble u,
							   gdouble v);
Vector spline2d_dimensional_normal                        (Spline2D * sp, gdouble u,
							   gdouble v);



InfluenceCoeffs *  influencecoeffs_new              (gint k);
void               influencecoeffs_destroy          (InfluenceCoeffs * ic);    
InfluenceCoeffs *  influencecoeffs_add_and_destroy  (InfluenceCoeffs * ic1,
						     InfluenceCoeffs * ic2);
void               spline2d_clear_self_influence_coefficients (Spline2D * sp);
SPPanel *          spline_to_poly_matrix                   (Spline2D * s, gdouble ue,
							    gdouble ve);

typedef struct {
  gdouble u0, u1;
  gdouble v0, v1;
  gdouble ue, ve;

  gpointer children[4];
  gint level;
  Point pe;
  gint k;
  gdouble d;
  gint ng;

  GArray * Pi, * Ni, * wJij;
  GPtrArray * Bu, * Bv;
  GPtrArray * Bux;
} GPCell;

GPCell *           gpcell_tree_new                          (SPPanel * spp);

void               gpcell_tree_destroy                      (GPCell * gpc);

void               spline_near_field_influence_coeff_recursive (SPPanel * spp,
								GPCell * gpc,
								Point p,
								gdouble * psi,
								gdouble * phi,
								gboolean edge);
InfluenceCoeffs *  sppanel_self_influence_coeff            (SPPanel * spp,
							    gint im, gint in);
InfluenceCoeffs *  spline_near_field_influence_coeff       (SPPanel * spp, Point p);
InfluenceCoeffs *  spline_far_field_influence_coeff        (SPPanel * spp,
							    Point p);
InfluenceCoeffs *  sppanel_nearfar_influence_coefficients  (SPPanel * spp, Spline2D * sp,
							    gdouble um, gdouble vn, Point p,
							    gboolean belongs_to_patch);
InfluenceCoeffs *  sppanel_spline_nearfar_influence_coefficients  (SPPanel * spp,
								   Spline2D * sp,
								   gdouble um,
								   gdouble vn, Point p,
								   gboolean belongs_to_patch);


typedef void       (* SelfInfluenceFunc)                    (SPPanel * spp,
							      gdouble up, gdouble vp,
							      Point p,
							      gdouble * psi,
							      gdouble * phi);
void               lachat_watson_self_influence_coefficients_qag (SPPanel * spp,
								  gdouble up, gdouble vp,
								  Point p,
								  gdouble * psi,
								  gdouble * phi);
void               lachat_watson_self_influence_coefficients (SPPanel * spp,
							      gdouble up, gdouble vp,
							      Point p,
							      gdouble * psi,
							      gdouble * phi);
void               rong_self_influence_coefficients          (SPPanel * spp,
							      gdouble up, gdouble vp,
							      Point p,
							      gdouble * psi,
							      gdouble * phi);
void               wamit_self_influence_coefficients         (SPPanel * spp,
							      gdouble up, gdouble vp,
							      Point p,
							      gdouble * psi,
							      gdouble * phi);
void               centered_wamit_self_influence_coefficients (SPPanel * spp,
							       gdouble up, gdouble vp,
							       Point p,
							       gdouble * psi,
							       gdouble * phi);
void               wamit_self_influence_coefficients_cauchy   (SPPanel * spp,
							       gdouble up, gdouble vp,
							       Point p,
							       gdouble * psi,
							       gdouble * phi);

GArray *           sppanel_change_polynomial_origin        (SPPanel * spp, gdouble uc,
							    gdouble vc);
S2P *              sppanel_spline_to_poly_matrix           (SPPanel * spp, Spline2D * sp,
				                            gdouble ue, gdouble ve);
gdouble            spp_characteristic_length               (SPPanel * spp, Spline2D * sp);

gdouble            calculate_added_mass                    (GSList * list);
void               print_potential                         (GSList * list);
void               print_potential_gauss                   (GSList * list);
void               print_potential_sphere                  (GSList * list);
void               print_potential_spheroid                (GSList * list);
void               print_potential_cylinder                (GSList * list);
void               print_cylinder_profile                  (GSList * list);


gdouble            zero_potential                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
Vector             zero_velocity                           (Spline2D * spp,
							    gdouble u, gdouble v,
							    gpointer data);
Vector             m_terms_dirichlet_bc                    (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            m1_term_dirichlet_bc                    (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            m2_term_dirichlet_bc                    (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            m3_term_dirichlet_bc                    (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            mode_forcing_1                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            mode_forcing_2                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            mode_forcing_3                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            mode_forcing_4                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            mode_forcing_5                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            mode_forcing_6                          (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
Vector             uniform_velocity                        (Spline2D * sp,
							    gdouble u, gdouble v,
							    gpointer data);
gdouble            uniform_normal_velocity                 (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            zero_normal_velocity                    (SPPanel * spp,
							    gint m, gint n,
							    gpointer data);
gdouble            neumann_forcing_term                    (SPPanel * spp,
							    Point p,
							    NeumannFunc Vn,
							    gpointer data);
gdouble            dirichlet_forcing_term                  (SPPanel * spp,
							    gdouble u, gdouble v,
							    Point p, gint var, NeumannFunc Vn, gpointer data);
gdouble            old_normal_velocity_forcing_term        (SPPanel * spp, Point p,
							    Vector * Vn);
void               print_velocity_on_surface               (GSList * list,
							    NeumannFunc Vn);
Vector             potential_gradient_on_surface           (Spline2D * sp,
							    gdouble u, gdouble v,
							    gint var);
Vector             potential_gradient_on_surface_gauss_point (Spline2D * sp,
							     GaussPoints * gp,
							     gint m, gint n,
							     gint var);
Vector           potential_gradient2d_on_surface_gauss_point (Spline2D * sp,
							      GaussPoints * gp,
							      gint m, gint n,
							      gint var);
Vector             gradient_on_surface_2d                  (Spline2D * sp,
							    gdouble u, gdouble v,
							    gint var);
void              apply_neumann_conditions                 (Spline2D * sp,
							    NeumannFunc func,
							    gpointer data,
							    gint var);
void              spline_set_var_to_constant               (Spline2D * sp,
							    gint var,
							    gdouble val);

void lachat_watson_self_influence_coefficients2 (SPPanel * spp,
						 GaussPoints * gp,
						 /* GArray * up, */
						 /* GArray * vp, */
						 /* GArray * pi, */
						 gdouble * psi,
						 gdouble * phi);

void spline_near_field_influence_coeff_recursive2 (SPPanel * spp,
						   GPCell * gpc,
						   GaussPoints * gp,
						   /* GArray * p, */
						   gdouble * psi,
						   gdouble * phi);

typedef struct {
  gdouble Pm, Pl, Ph, Pfk, Pr, phi2, phi, Pm2, Pl2, Pfk2;
  gdouble P1, Ph1, Pfk1, phi1;
} Pressure;

CCSProblem * spline2d_build_galerkin_fit_matrix   (Spline2D * sp);
gsl_vector * build_galerkin_rhs_gauss             (Spline2D * sp,
						   GaussFunc func,
						   gpointer data,
						   Spline2DFunc bc_func,
						   gpointer bc_data,
						   gsl_vector * rhs);
void         spline2d_copy_problem_solution       (Spline2D * sp,
						   gsl_vector * lhs,
						   gint var);

void         spline2d_add_problem_solution        (Spline2D * sp,
						   gsl_vector * lhs,
						   gint var,
						   gint var0);

typedef struct _Simulation Simulation;

void lachat_watson_ye_self_influence_coefficients (SPPanel * spp,
						   gdouble up, gdouble vp,
						   Point p,
						   gdouble * psi,
						   gdouble * phi);
void qin_self_influence_coefficients              (SPPanel * spp,
						   gdouble up, gdouble vp,
						   Point p,
						   gdouble * psi,
						   gdouble * phi);
Point newton_raphson (SPPanel * spp, Point p);

gdouble zero_spline2d_func                       (Spline2D * sp,
						  gdouble u, gdouble v,
						  gpointer data);
void print_potential_square (GSList * list);

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
  /* Matrix */
  GPtrArray * A, * A2;
  /* Indices */
  GPtrArray * Ai, * Ai2;
  GArray * rhs, * lhs, * rhs2, * lhs2;
  /* Ghost cells */
  GArray * gxi0, * gxiN;
  GArray * geta0, * getaN;
  gint istart, jstart;
} LinearProblem;

/* LinearProblem methods */
LinearProblem *  linear_problem_new                ();
void             linear_problem_destroy            (LinearProblem * lp);
void             linear_problem_init_size          (LinearProblem * lp,
						    gint N, gint M);
void             linear_problem_add_stencil        (LinearProblem * lp,
						    gint i, gint j,
						    gdouble coeff);
gboolean         linear_problem_solve_using_SOR    (LinearProblem * lp,
						    gdouble tolerance,
						    gboolean verbose);
gboolean         linear_problem_SOR_iteration      (LinearProblem * lp);
void             build_potential_sparse_problem    (GSList * list,
						    LinearProblem * lp,
						    gint SIZE);
gboolean         linear_problem_solve_using_gauss  (LinearProblem * lp);
void             copy_solution_to_patches          (GSList * list,
						    LinearProblem * lp);
gint             spline_numbering                  (GSList * list);
gint             problem_size                      (GSList * list);
void             compute_self_influence_coefficients (Spline2D * sp);

typedef struct {
  gsl_matrix * A;
  gsl_vector * rhs;
  gint istart, jstart;
} BoundaryProblem;

typedef struct {
  BoundaryProblem * neumann;
  BoundaryProblem * dirichlet;
} BoundarySubProblem;

BoundaryProblem *   boundary_problem_new                           (gint N, gint M);
void                 boundary_problem_destroy                      (BoundaryProblem * bp);
BoundarySubProblem * boundary_subproblem_new                       (gint N, gint M);
void                 boundary_subproblem_destroy                   (BoundarySubProblem * b);


BoundarySubProblem * build_boundary_subproblem_galerkin            (Spline2D * sp1,
								    Spline2D * sp2,
								    SelfInfluenceFunc self_influence_coeffs);
BoundarySubProblem * build_boundary_subproblem_collocation         (Spline2D * sp1,
								    Spline2D * sp2,
								    SelfInfluenceFunc self_influence_coeffs);
void              boundary_problem_assemble_neumann                (BoundaryProblem * bp,
								    GPtrArray * sub_problems);
void              boundary_problem_assemble_neumann_rhs            (BoundaryProblem * bp,
								    GPtrArray * sub_problems);
void              boundary_problem_assemble_dirichlet              (BoundaryProblem * bp,
								    GPtrArray * sub_problems);
void              boundary_problem_assemble_dirichlet_rhs          (BoundaryProblem * bp,
								    GPtrArray * sub_problems);
void              boundary_problem_direct_LU_solve                 (BoundaryProblem * bp);
void              boundary_problem_copy_solution_to_patches        (GSList * list,
								    BoundaryProblem * bp,
								    gint var);

void              boundary_subproblem_build_rhs_neumann            (Spline2D * sp1,
								    Spline2D * sp2,
								    BoundarySubProblem * bsp,
								    gint var);
void              boundary_subproblem_build_rhs_dirichlet          (Spline2D * sp1,
								    Spline2D * sp2,
								    BoundarySubProblem * bsp,
								    gint var);
//CCSProblem *      spline2d_build_galerkin_fit_matrix                (Spline2D * sp);
/* void              spline2d_copy_problem_solution                    (Spline2D * sp, */
/* 								     gsl_vector * lhs, */
/* 								     gint var); */
CCSProblem *      spline2d_build_galerkin_fit_matrix_no_metric      (Spline2D * sp);
void              spline2d_build_freesurface_galerkin_fit_matrix    (Spline2D * sp);
void              spline2d_build_freesurface_noflux_galerkin_fit_matrix (Spline2D * sp);
CCSProblem * periodic_fs_build_galerkin_fit_noflux_test_matrix/* _neumann */ (Spline2D * splines);
gsl_vector * periodic_fs_build_galerkin_noflux_test_rhs_gauss (Spline2D * splines, 
							       GaussFunc func, 
							       gpointer data, 
							       Spline2DFunc bc_func, 
							       gpointer bc_data,
							       gsl_vector * rhs);
Spline2D *        rectangular_grid                                  (gint M, gint N);
Spline2D *        parametric_grid                                   (gint M, gint N,
								     GaussFunc func_x,
								     GaussFunc func_y,
								     gpointer data);
Spline2D *        spline2d_parametric_patch                         (gint M, gint N,
								     GaussFunc func_x,
								     GaussFunc func_y,
								     GaussFunc func_z,
								     gpointer data,
								     gint k, gint ninner,
								     gint nouter);
Spline2D *        spline2d_parametric_periodic_patch                (gint M, gint N,
								     GaussFunc func_x,
								     GaussFunc func_y,
								     GaussFunc func_z,
								     gpointer data);
Spline2D *        parametric_grid2                                  (gint M, gint N,
								     GaussFunc func_x,
								     GaussFunc func_y,
								     gpointer data);
Spline2D *        parametric_grid3                                  (gint M, gint N,
								     GrevilleFunc func_x,
								     GrevilleFunc func_y,
								     gpointer data);
gboolean          gsl_problem_solve_using_SOR                       (gsl_matrix * Aij,
								     gsl_vector * rhsi,
								     gsl_vector * lhsi,
								     gdouble tolerance,
								     gboolean verbose);

void             spline2d_fit_greville                              (Spline2D * sp,
								     GrevilleFunc func,
								     gpointer data,
								     gint var);
void             spline2d_fit_greville_periodic                     (Spline2D * sp,
								     GrevilleFunc func,
								     gpointer data,
								     gint var);

void             spline2d_fit_galerkin                              (Spline2D * sp,
								     GaussFunc func,
								     gpointer data,
								     gint var);
void             boundary_problem_direct_sparse_solve               (gdouble ** A,
								     gdouble * RHS,
								     gint size);

void             apply_dirichlet_conditions                         (Spline2D * sp,
								     GaussFunc func,
								     gpointer data,
								     gint var);
/* gsl_vector *     build_galerkin_rhs_gauss                           (Spline2D * sp, */
/* 								     GaussFunc func, */
/* 								     gpointer data); */

gsl_vector *     build_galerkin_rhs_gauss_no_metric                 (Spline2D * sp,
								     GaussFunc func,
								     gpointer data,
								     Spline2DFunc bc_func,
								     gpointer bc_data);

#if CUDA
/*** CUDA ***/
void              cuda_init                                      (CUdevice * dev,
								  CUcontext * context);
void              cuda_finalize                                  (CUcontext * context);

typedef struct {
  double * A, * rhs;
  magma_int_t * ipiv;
} CudaLUProblem;

CudaLUProblem *   cuda_lu_problem_new                            (gint M, gint N);
void              cuda_lu_problem_destroy                        (CudaLUProblem * cp);
CudaLUProblem *   cuda_lu_factorise                              (gsl_matrix * A);
void              cuda_solve_factorised_lu                       (CudaLUProblem * cp,
								  gsl_vector * rhs);
void              cuda_lu_solve                                  (gsl_matrix * A,
								  gsl_vector * rhs);

typedef struct {
  double * du, * dvt, * drhs, * dw; 
  double * u, * vt, * s;
  gsl_vector * w;
  gint M, N;
} CudaSVDProblem;

CudaSVDProblem *  cuda_svd_problem_new                           (gint M, gint N);
void              cuda_svd_problem_destroy                       (CudaSVDProblem * cp);
CudaSVDProblem *  cuda_svd_factorise                             (gsl_matrix * A);
void              cuda_solve_factorised_svd                      (CudaSVDProblem * cp,
								  gsl_vector * rhs);


#endif
#if OPENCL
/*** CUDA ***/
void              opencl_init                                      ();
void              opencl_finalize                                  ();

typedef struct {
  double * A, * rhs;
  magma_int_t * ipiv;
  magma_queue_t  queue;
  magma_device_t device;
  int num;
} OpenCLLUProblem;

OpenCLLUProblem * opencl_lu_problem_new                   (gint M,
							   gint N);
void              opencl_lu_problem_destroy               (OpenCLLUProblem * cp);
OpenCLLUProblem * opencl_lu_factorise                     (gsl_matrix * A);
void              opencl_solve_factorised_lu              (OpenCLLUProblem * cp,
							   gsl_vector * rhs);
void             opencl_lu_solve                          (gsl_matrix * A,
							   gsl_vector * rhs);

typedef struct {
  double * du, * dvt, * drhs, * dw; 
  double * u, * vt, * s;
  gsl_vector * w;
  gint M, N;
  magma_queue_t  queue;
  magma_device_t device;
  int num;
} OpenCLSVDProblem;
#endif
void              cuda_boundary_problem_copy_solution_to_patches (GSList * list,
								  BoundaryProblem * bp,
								  gint var);

/**** SUPERLU ****/


CCSProblem *     ccs_problem_new                                 ();
void             ccs_problem_destroy                             (CCSProblem * ccs);
void             ccs_problem_lu_solve                            (CCSProblem * ccs,
								  gsl_vector * rhs);
void             spline2d_fit_greville_border                    (Spline2D * sp,
								  Spline2DFunc func,
								  gpointer data,
								  gint var);

/***** PLASMA *****/

typedef struct {
  double * A, * rhs;
  double * L;
  int * ipiv;
} PlasmaLUProblem;

PlasmaLUProblem * plasma_lu_problem_new                             (gint M, gint N);
void              plasma_lu_problem_destroy                         (PlasmaLUProblem * pp);
PlasmaLUProblem * plasma_lu_factorise                               (gsl_matrix * A);
void              plasma_solve_factorised_lu                        (PlasmaLUProblem * pp,
								     gsl_vector * rhs);
void              plasma_lu_solve                                   (gsl_matrix * A,
								     gsl_vector * rhs);
void              plasma_boundary_problem_copy_solution_to_patches  (GSList * list,
								     BoundaryProblem * bp,
								     gint var);

/***** Periodic splines *****/


Spline2D *     periodic_fs_new                                 (gint M, gint N,
								gint k,
								gint ninner,
								gint nouter);
gint           periodic_fs_numbering                           (Spline2D * fs);
void           periodic_fs_init_panels                         (Spline2D * sp);
void           periodic_fs_reinit_panels_physical_quantities   (Spline2D * sp);
void           periodic_fs_copy_problem_solution               (Spline2D * sp,
								gsl_vector * lhs,
								gint var);
void           periodic_fs_add_problem_solution                (Spline2D * splines,
								gsl_vector * lhs,
								gint var,
								gint var0);
CCSProblem *   periodic_fs_build_galerkin_fit_matrix           (Spline2D * sp);
gsl_vector *   periodic_fs_build_galerkin_rhs_gauss            (Spline2D * sp,
								GaussFunc func,
								gpointer data,
								Spline2DFunc bc_func,
								gpointer bc_data,
								gsl_vector * rhs);

gsl_vector *   periodic_fs_build_galerkin_rhs_gauss_no_metric  (Spline2D * sp,
								GaussFunc func,
								gpointer data,
								Spline2DFunc bc_func,
								gpointer bc_data);
CCSProblem *   periodic_fs_build_galerkin_fit_matrix_no_metric (Spline2D * sp);
void           periodic_fs_copy_problem_solution_no_metric     (Spline2D * splines,
								gsl_vector * lhs,
								gint var);
void           periodic_fs_elliptic_smoothing                  (Spline2D * splines);
Spline2D *     spline2d_regrid                                 (Spline2D * sp,
								gint M, gint N);
void           calculate_dzz_direct                            (GSList * patches);
gsl_vector *   build_galerkin_rhs_gauss_no_metric_bc_dirichlet (Spline2D * sp,
								GaussFunc func,
								gpointer data,
								Spline2DFunc bc_func,
								gpointer bc_data);
CCSProblem *   spline2d_build_galerkin_fit_matrix_no_metric_bc_dirichlet 
                                                               (Spline2D * sp);
CCSProblem * spline2d_build_galerkin_fit_matrix_bc           (Spline2D * sp);

gsl_vector * build_galerkin_rhs_gauss_bc                     (Spline2D * sp,
							      GaussFunc func,
							      gpointer data,
							      Spline2DFunc bc_func,
							      gpointer bc_data,
							      gsl_vector * rhs);

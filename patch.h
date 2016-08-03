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
  GPtrArray * rows;
  Transformation * t;
} Patch;

void     patch_destroy                        (Patch * p);
void     patch_add_row                        (Patch * p);
Patch *  patch_new                            ();
void     patch_add_panel                      (Patch * p, Panel panel);
void     patch_tag_borders                    (Patch * p);
void     patch_print                          (Patch * p, FILE * fp);
void     patch_center_print                   (Patch * p, FILE * fp);
gdouble  patch_eval                           (Patch * p, gdouble i,
					       gdouble j, gint var);
gdouble  patch_eval_dx                        (Patch * p, gdouble i,
					       gdouble j, gint var);
gdouble  patch_eval_dy                        (Patch * p, gdouble i,
					       gdouble j, gint var);
gdouble  patch_eval_dxx                       (Patch * p, gdouble i,
					       gdouble j, gint var);
gdouble  patch_eval_dyy                       (Patch * p, gdouble i,
					       gdouble j, gint var);
gdouble  patch_eval_dxy                       (Patch * p, gdouble i,
					       gdouble j, gint var);
gdouble  patch_eval_int                       (Patch * p, gdouble i,
					       gdouble j, gint var);
Point    patch_eval_point                     (Patch * patch, gdouble i,
					       gdouble j);
Point    patch_eval_transformed_point         (Patch * patch, gdouble i,
					       gdouble j);

typedef void (* PanelFunction) (Panel * p, gpointer data);
Panel *  patch_panel_get                      (Patch * patch, gint i,
					       gint j);
void     patch_forall_panels                  (Patch * patch, PanelFunction func,
					       gpointer data);

GSList * sort_list_by_angle                   (GSList * l);
GSList * patch_intersect_with_free_surface    (Patch * patch);

Vector   patch_normal                         (Patch * p, gdouble i,
					       gdouble j);
Vector   patch_unit_normal                    (Patch * p, gdouble i,
					       gdouble j);
gdouble  patch_local_metric                   (Patch * p, gdouble i,
					       gdouble j);
void     patch_normal_print                   (Patch * p, FILE * fp);
Vector   patch_adaptive_panel_integral        (Patch * patch, HeightCurve hz,
					       gint var);
gdouble  patch_adaptive_panel_scalar_integral (Patch * patch, HeightCurve hz,
					       gint var);
void     patch_compute_metric                 (Patch * patch);
gdouble  patch_fast_panel_scalar_integral     (Patch * patch, HeightCurve hz,
					       int var);
void     patch_compute_volume_times_metric    (Patch * patch, HeightCurve hz);
gdouble  patch_fast_panel_scalar_integral2    (Patch * patch, HeightCurve hz,
					       int var);
void     patch_fit_panels                     (Patch * patch);
void     patch_prefit_panels                  (Patch * patch);


//************************************************************//
Spline2D *   spline2d_fit_geometry2_old  (Patch * patch,
					  gint M, gint N,
					  gint order,
					  gint inner, gint outer,
					  gboolean flip);
Spline2D *   spline2d_fit_geometry2      (Patch * patch,
					  gint M, gint N,
					  gint order,
					  gint inner, gint outer,
					  gboolean flip_u,
					  gboolean flip_v,
					  gboolean centripetal_reparam,
					  gboolean swap);
Spline2D *  spline2d_interpolate_geometry (Patch * patch, gint M, gint N,
					   gint order, gint inner, gint outer,
					   gboolean flip_u, gboolean flip_v);
Spline2D *   spline2d_fit_geometry3      (FILE * finput,
					  gint M, gint N,
					  gint order, gint inner,
					  gint outer,
					  gint NX, gint NY);
Spline2D *   spline2d_fit_geometry4      (FILE * finput,
					  gint M, gint N,
					  gint order,
					  gint inner,
					  gint outer,
					  gint NX, gint NY);
Spline2D *   spline2d_fit_geometry5      (FILE * finput,
					  gint M, gint N,
					  gint order,
					  gint inner,
					  gint outer,
					  gint NX, gint NY);
Spline2D *   spline2d_fit_geometry6      (GArray * input,
					  gint M, gint N,
					  gint order,
					  gint inner,
					  gint outer,
					  gint NX, gint NY);
Spline2D * spline2d_fit_geometry8 (FILE * finput, gint M, gint N, gint order, gint inner, gint outer, gint NX, gint NY);

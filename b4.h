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

gdouble  b4_patch_eval                           (Patch * p, gdouble i,
						  gdouble j, gint var);
gdouble  b4_patch_eval_dx                        (Patch * p, gdouble i,
						  gdouble j, gint var);
gdouble  b4_patch_eval_dy                        (Patch * p, gdouble i,
						  gdouble j, gint var);
void     b4_patch_fit_panels                     (Patch * patch);

Point    b4_patch_eval_point                     (Patch * patch,
						  gdouble i, gdouble j);
Point    b4_patch_eval_transformed_point         (Patch * patch,
						  gdouble i, gdouble j);

GSList * b4_patch_intersect_with_free_surface    (Patch * patch,
						  HeightCurve hz);

Vector   b4_patch_normal                         (Patch * p, gdouble i,
						  gdouble j);
Vector   b4_patch_unit_normal                    (Patch * p, gdouble i,
						  gdouble j);
gdouble  b4_patch_local_metric                   (Patch * p, gdouble i,
						  gdouble j);
void     b4_patch_normal_print                   (Patch * p, FILE * fp);
void     b4_patch_store_derivatives              (Patch * patch);

Vector   b4_patch_adaptive_panel_integral        (Patch * patch, HeightCurve hz,
						  gint var);
gdouble  b4_patch_adaptive_panel_scalar_integral (Patch * patch, HeightCurve hz,
						  gint var);
void     b4_patch_compute_metric                 (Patch * patch);
gdouble  b4_patch_fast_panel_scalar_integral     (Patch * patch, HeightCurve hz,
						  gint var);
Vector   b4_patch_adaptive_panel_flux_integral   (Patch * patch, HeightCurve hz,
						  gint var);
void     b4_patch_compute_volume_times_metric    (Patch * patch, HeightCurve hz);

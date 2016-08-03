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
  GSList * patches;
  Boundaries * b;
  HeightCurve hz;
} Surface;

typedef struct _WaveParams WaveParams; 

typedef gdouble    (* ScalarWaveFunc)                          (WaveParams * wp,
								Point p, gdouble t);
typedef Vector     (* VectorWaveFunc)                          (WaveParams * wp,
								Point p, gdouble t);
typedef gdouble    (* XYTWaveFunc)                             (gdouble x, gdouble y,
								gdouble t,
								gpointer data
								/* WaveParams * wp */);

struct _WaveParams {
  gdouble g, A, w, h, k;
  gdouble cosb, sinb;
  gdouble r1, r2, Cs, Cw;
  gdouble r1_inner, r2_inner;
  
  ScalarWaveFunc wave_potential;
  ScalarWaveFunc wave_potential_dt;
  VectorWaveFunc wave_potential_gradient;
  ScalarWaveFunc wave_potential_dz_dt;
  VectorWaveFunc wave_potential_z_derivative_gradient;
  XYTWaveFunc wave_elevation;
  VectorWaveFunc wave_elevation_gradient;
  ScalarWaveFunc wave_elevation_time_derivative;
  VectorWaveFunc wave_normal_time_derivative;

  ScalarWaveFunc wave_potential2;
  ScalarWaveFunc wave_potential2_dt;
  VectorWaveFunc wave_potential2_gradient;
  ScalarWaveFunc wave_potential2_dz_dt;
  VectorWaveFunc wave_potential2_z_derivative_gradient;
  XYTWaveFunc wave_elevation2;
  VectorWaveFunc wave_elevation2_gradient;
  ScalarWaveFunc wave_elevation2_time_derivative;
  VectorWaveFunc wave_normal2_time_derivative;
};

Surface     * surface_new                                       ();
void          surface_destroy                                   (Surface * s);
gdouble       surface_eval                                      (Surface * s, gdouble x,
								 gdouble y);
void          surface_generate_grid                             (Surface * s, gboolean flip);
void          surface_print_grid                                (Surface * s, FILE * fp);
void          spline2d_surface_print_grid                       (Surface * s, FILE * fp);
void          spline2d_discretize_free_surface                  (Spline2D * sp,
								 WaveParams * wp,
								 gdouble t);
void          spline2d_discretize_bathymetry                    (Spline2D * sp,
								 WaveParams * wp,
								 gdouble t);
typedef struct {
  Surface * s; 
} FreeSurface;

FreeSurface * freesurface_new                                   ();
void          freesurface_init                                  (FreeSurface * f,
								 WaveParams * wp);
void          freesurface_destroy                               (FreeSurface * f);

typedef struct {
  Surface * s;
} Bathymetry;

Bathymetry *  bathymetry_new                                    ();
void          bathymetry_init                                   (Bathymetry * b,
								 DCurve * dcb,
								 WaveParams * wp);
void          bathymetry_destroy                                (Bathymetry * b);


gdouble       finite_depth_wave_potential                       (WaveParams * wp,
								 Point p, gdouble t);
gdouble       finite_depth_wave_potential_dt                    (WaveParams * wp,
								 Point p, gdouble t);
Vector        finite_depth_wave_potential_gradient              (WaveParams * wp,
								 Point p, gdouble t);
gdouble       finite_depth_wave_potential_dz_dt                 (WaveParams * wp,
								 Point p, gdouble t);
Vector        finite_depth_wave_potential_z_derivative_gradient (WaveParams * wp,
								 Point p, gdouble t);
gdouble       finite_depth_wave_elevation                       (gdouble x, gdouble y,
								 gdouble t,
								 gpointer data);
Vector        finite_depth_wave_elevation_gradient              (WaveParams * wp,
								 Point p, gdouble t);
gdouble       finite_depth_wave_elevation_time_derivative       (WaveParams * wp,
								 Point p, gdouble t);
Vector        finite_depth_wave_normal_time_derivative          (WaveParams * wp,
								 Point p, gdouble t);

gdouble       infinite_depth_wave_potential                     (WaveParams * wp,
								 Point p, gdouble t);
gdouble       infinite_depth_wave_potential_dt                  (WaveParams * wp,
								 Point p, gdouble t);
Vector        infinite_depth_wave_potential_gradient            (WaveParams * wp,
								 Point p, gdouble t);
gdouble       infinite_depth_wave_potential_dz_dt               (WaveParams * wp,
								 Point p, gdouble t);
Vector        infinite_depth_wave_potential_z_derivative_gradient (WaveParams * wp,
								   Point p, gdouble t);

//gdouble       infinite_depth_wave_elevation                     (WaveParams * wp,
//								 Point p, gdouble t);
gdouble       infinite_depth_wave_elevation                     (gdouble x,
								 gdouble y,
								 gdouble t,
								 gpointer data);
Vector        infinite_depth_wave_elevation_gradient            (WaveParams * wp,
								 Point p, gdouble t);
gdouble       infinite_depth_wave_elevation_time_derivative     (WaveParams * wp,
								 Point p, gdouble t);
Vector        infinite_depth_wave_normal_time_derivative        (WaveParams * wp,
								 Point p, gdouble t);


gdouble       solve_dispersion_relation                         (WaveParams * wp);
void          print_free_surface                                (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void          print_free_surface_tmp                                (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void         print_free_surface_tecplot                         (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void         print_free_surface_potential                       (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void         print_hull_potential                               (GSList * list,
  WaveParams * wp,
								 gdouble t);
void         free_surface_galerkin_fit                          (Spline2D * sp,
								 Spline2DFunc func,
								 gpointer data,
								 gint var);
void         elliptic_smoothing_periodic                        (Spline2D * sp);
gdouble      zero_scalar_wave_func                              (WaveParams * wp,
								 Point p, gdouble t);
Vector       zero_vector_wave_func                              (WaveParams * wp,
								 Point p, gdouble t);
gdouble      zero_wave_elevation                                (gdouble x, gdouble y,
								 gdouble t,
								 gpointer data
								 /* WaveParams * wp */);
void         print_free_surface_hr                              (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void         print_waterline                                    (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void         print_free_surface_tecplot_hr                      (GSList * list,
								 WaveParams * wp,
								 gdouble t);
void         print_free_surface_mayavi                          (GSList * list,
								 WaveParams * wp,
								 gdouble t);


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
  float * val;
  gint nx, ny;
  gint varid;
  gint it;
  float t;
} NetCDFData;

typedef struct {
  gchar * file_name;
  gint ncid;
  gint it;
  gdouble t1, t2, time;

  gint xid, yid, tid;
  size_t nx, ny, nt;
  gint xvarid, yvarid, tvarid;
  float * x, * y, * t;
  float tmax;
  gint varp1, varp2;

  gdouble x0, y0;
  gdouble calpha, salpha;

  NetCDFData * phi_1, * elevation_1, * u_1, * v_1;
  NetCDFData * phi_2, * elevation_2, * u_2, * v_2;
  
} NetCDFForcing;

typedef struct {
  gsl_vector * rhs_phi;
  gsl_vector * rhs_elevation;
  gsl_vector * rhs_u;
  gsl_vector * rhs_v;
} DoubleRhs;

NetCDFForcing *  netcdf_forcing_new                            (gchar * file_name,
								gdouble xs,
								gdouble ys,
								gdouble angles);
void             netcdf_forcing_destroy                        (NetCDFForcing * ncf);
void             netcdf_build_spatial_interpolation_operators  (Spline2D * splines,
								NetCDFForcing * ncdf);

gint             readSurffnc                                   (char * surf_file,
								gdouble ctime);
void             nc_test                                       ();
/* void             add_netcdf_fk_forces                          (Simulation * sim, */
/* 								Forces * f, gdouble t, */
/* 								gdouble u[6], gdouble x[6]); */
NetCDFData *     netcdf_load_array                             (NetCDFForcing * ncf,
								NetCDFData * ncd,
								gchar * field_name,
								gint it);
DoubleRhs *      netcdf_build_double_galerkin_rhs              (Spline2D * sp,
								NetCDFForcing * ncdf,
								Hull * hull);
DoubleRhs *      netcdf_build_linear_double_galerkin_rhs       (Spline2D * sp,
								NetCDFForcing * ncdf,
								Hull * hull,
								gdouble coeff);

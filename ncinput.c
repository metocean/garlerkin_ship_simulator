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

#include "structures.h"
#include "patch.h"
#include "linearproblem.h"
/* #include "boundaries.h" */
/* #include "surface.h" */
#include "hull.h"
#include "motion.h"
#include "ncinput.h"

int nc_get_vara_time (int ncid, int tvarid, size_t it, size_t i1, gdouble t2)
{
  int stat;
  float tdum;

  nc_type nctype;

  stat = nc_inq_vartype (ncid, tvarid, &nctype);

  if ( stat < 0) {
    return stat;
  }

  if (nctype == NC_FLOAT) {
    stat = nc_get_vara_float (ncid, tvarid, &it, &i1, &tdum);
    t2 = tdum;
  }
  else {
    stat = nc_get_vara_double (ncid, tvarid, &it, &i1, &t2);
  }

  return stat;
}

/**
   ncid : Id of netcdf file
   it : Index of time
   arr : Array of data?
   ni : first dimension of arr ?
   nj : second dimension of arr ?
 **/
gint read_ncarray (gint ncid, gint it, float * arr, gint ni, gint nj, float * fac, gint * inti, gint * intj, gint ivr, gint nx, gint ny)
{
  size_t strt[3], icnt[3];
  gint i, j, k, stat;
  float bufin[nx][ny];
  
  strt[0] = 0;
  strt[1] = 0;
  strt[2] = it;
  icnt[0] = nx;
  icnt[1] = ny;
  icnt[2] = 1;

  if (ivr < 0)
    return;

  // Reads the field starting at 0,0,it over a length nx,ny
  // This is the all the data
  stat = nc_get_vara_float (ncid, ivr, strt, icnt, &bufin[0][0]);

  // Applies the interpolation to the data
  // and stores the result
  for ( i = 0; i < ni; i++) {
    for ( j = 0; j < nj; j++) {
      *(arr + i*nj + j) = 0.; // This is were the results is stored
      float fac0 = 0.;
      for ( k = 0; k < 4; k++) {
  	if ( *(fac + (i*nj + j)*4 + k) > 0. &&
	     bufin[(*(inti + (i*nj + j)*4 + k))][(*(intj + (i*nj + j)*4 + k))] > -99. &&
	     bufin[(*(inti + (i*nj + j)*4 + k))][(*(intj + (i*nj + j)*4 + k))] < 1.e10 ) {
  	  *(arr + i*nj + j) = *(arr + i*nj + j) + *(fac + (i*nj + j)*4 + k)*bufin[(*(inti + (i*nj + j)*4 + k))][(*(intj + (i*nj + j)*4 + k))]; // sum of data*weights
	  fac0 += *(fac + (i*nj + j)*4 + k); // sum of weights
	}
      }
      if ( fac0 > 0.)
	*(arr + i*nj + j) = *(arr + i*nj + j)/fac0; // sum of data*weights/sum of weights
    }
  }


}

/* subroutine build_interp(ncid,nx,ny, */
/*      &factor,inti,intj,xx,yy,mask,ni,nj,mskd) */
/**
 ncid : id of netcdf file
 nx : dimension of x data in netcdf file
 ny : dimension of y data in netcdf file
 factor : interpolation coefficients (result)
 inti : index of x data to use for interpolation
 intj : index of y data to use for interpolation
 xx : x coordinates for interpolated values
 yy : y coordinates for interpolated values
 mask : not applicable
 ni : length of xx ?
 nj : length of yy ?
 mskd : not applicable
 **/
gint build_interp (gint ncid, gint nx, gint ny, gint ni, gint nj, float * xx, float * yy, float * factor, gint * inti, gint * intj)
{
  gint stat;
  gint xvarid, yvarid;

  float xax[nx], yax[ny];

  /* float factor[ni][nj][4]; */

  gint i, j;
  /* gint inti[ni][nj][4],  intj[ni][nj][4]; */
  /* gdouble xx[ni][nj], yy[ni][nj]; */

  gdouble x1, y1;

  // returns id of var given its name
  stat = nc_inq_varid (ncid, "longitude", &xvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid (ncid, "x", &xvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid (ncid, "lon", &xvarid);

  stat = nc_inq_varid (ncid, "latitude", &yvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid (ncid, "y", &yvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid (ncid, "lat", &yvarid);

  // Read all x and y at once and store in xax, yax
  stat = nc_get_var_float (ncid, xvarid, xax);
  stat = nc_get_var_float (ncid, yvarid, yax);

  for ( i = 0; i < ni; i++) {
    for ( j = 0; j < nj; j++) {
      // Find position of x
      gint ix = 0;
      while (xax[ix] < *(xx + i*nj + j) ) {
  	ix++;
  	if (ix > nx-1)
  	  break;
      }

      gint ix1 = MAX(0, ix-1);
      gint ix2 = MIN(ix, nx-1);
      /* inti[i][j][0] = ix1; */
      /* inti[i][j][1] = ix2; */
      /* inti[i][j][2] = ix1; */
      /* inti[i][j][3] = ix2; */
      *(inti + (i*nj + j)*4 + 0) = ix1;
      *(inti + (i*nj + j)*4 + 1) = ix2;
      *(inti + (i*nj + j)*4 + 2) = ix1;
      *(inti + (i*nj + j)*4 + 3) = ix2;

      gdouble dxint=xax[ix2]-xax[ix1];
      // Find lever
      if (dxint > 0)
  	x1 = (/* xx[i][j] */*(xx + i*nj + j)-xax[ix1])/dxint;
      else
  	x1 = 1.;
      
      // Same for y
      gint iy = 0;
      while (yax[iy] < /* yy[i][j] */*(yy + i*nj + j)) {
  	iy++;
  	if ( iy > (ny-1))
  	  break;
      }
      gint iy1 = MAX(0, iy-1);
      gint iy2 = MIN(iy-1, ny-1);
      /* intj[i][j][0] = iy1; */
      /* intj[i][j][1] = iy1; */
      /* intj[i][j][2] = iy2; */
      /* intj[i][j][3] = iy2; */
      *(intj + (i*nj + j)*4 + 0) = iy1;
      *(intj + (i*nj + j)*4 + 1) = iy1;
      *(intj + (i*nj + j)*4 + 2) = iy2;
      *(intj + (i*nj + j)*4 + 3) = iy2;

      gdouble dyint = yax[iy2] - yax[iy1];
      if (dyint > 0.)
  	y1 = (/* yy[i][j] */*(yy + i*nj + j)-yax[iy1])/dyint;
      else
  	y1 = 1.;

      // Interpolation operators
      /* factor [i][j][0] = (1.-x1)*(1.-y1); */
      /* factor [i][j][1] = x1*(1.-y1); */
      /* factor [i][j][2] = (1.-x1)*y1; */
      /* factor [i][j][3] = x1*y1; */

      *(factor + (i*nj + j)*4 + 0) = (1.-x1)*(1.-y1);
      *(factor + (i*nj + j)*4 + 1) = x1*(1.-y1);
      *(factor + (i*nj + j)*4 + 2) = (1.-x1)*y1;
      *(factor + (i*nj + j)*4 + 3) = x1*y1;
      // This must be the end product	
    }
  }
}







gint readSurffnc (char * surf_file, gdouble ctime)
{
  gint stat;
  gint ncid;

  gint iswitch = 1;

  gint im, jm, nvar;
  gdouble outarray[im][jm][nvar];

  int tvarid, xdimid, ydimid, tdimid;
  size_t nxdim, nydim, ntdim;

  gdouble t1, t2;

  // Open netcdf file and initial setup
  //	iswitch=1:   Surface file
  // iswitch=0:	 Hotfile

  /* if ( iswitch == 1 || iswitch == 0 ) { */
  /*   if ( iswitch == 1 ) { */
  /*     stat = nc_open (surf_file, 0, &ncid); */
  /*     if ( stat != NC_NOERR ) { */
  /* 	fprintf (stderr, "  *** Cannot open file: %s *** \n", surf_file); */
  /* 	g_assert_not_reached (); */
  /*     } */
	
  /*   } */
  /*   else if ( iswitch == 0 ) { */
      
  /*   } */

  /* } */

  stat = nc_open (surf_file, 0, &ncid);
  if ( stat != NC_NOERR ) {
    fprintf (stderr, "  *** Cannot open file: %s *** \n", surf_file);
    g_assert_not_reached ();
  }

  // Alloc
  float bufarray1[im][jm][nvar];
  float bufarray2[im][jm];
  float newarray[im][jm];
  float efacint[im][jm][4], ufacint[im][jm][4], vfacint[im][jm][4];
  int einti[im][jm][4], eintj[im][jm][4];
  int uinti[im][jm][4], uintj[im][jm][4];
  int vinti[im][jm][4], vintj[im][jm][4];
  int it;

  // Get id of x variable
  stat = nc_inq_dimid (ncid, "longitude", &xdimid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (ncid, "lon", &xdimid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (ncid, "x", &xdimid);

  fprintf (stderr, "x: %i \n", xdimid);

  // Get id of y variable
  stat = nc_inq_dimid (ncid, "latitude", &ydimid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (ncid, "lat", &ydimid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (ncid, "y", &ydimid);

  fprintf (stderr, "y: %i \n", ydimid);

  // Get id of time variable
  stat = nc_inq_dimid (ncid, "time", &tdimid);

  fprintf (stderr, "t: %i \n", tdimid);

  /* gint nxdim, nydim, ntdim; */
  // Get size of x, y, time dimensions
  stat = nc_inq_dimlen (ncid, xdimid, &nxdim);
  stat = nc_inq_dimlen (ncid, ydimid, &nydim);
  stat = nc_inq_dimlen (ncid, tdimid, &ntdim);

  fprintf (stderr, "nx ny nt: %i %i %i \n", (gint) nxdim, (gint) nydim, (gint) ntdim);

  // time variable id and read first value
  stat = nc_inq_varid (ncid, "time", &tvarid);

  fprintf (stderr, "tvar: %i \n", tvarid);

  gint x0id;
  nc_type x0type;
  size_t lenp;
  //  stat = nc_inq_att (ncid, x0id, "x0", &x0type, &lenp);

  // stat = nc_inq_varid (ncid, "x0", NC_GLOBAL);

  
  // Get informations on frame of reference of data.
  float x0, y0, angle;
  stat = nc_get_att_float (ncid, NC_GLOBAL, "x0", &x0);
  stat = nc_get_att_float (ncid, NC_GLOBAL, "y0", &y0);
  stat = nc_get_att_float (ncid, NC_GLOBAL, "angle", &angle);
  fprintf (stderr, "x0 %f y0: %f Angle: %f \n", x0, y0, angle);


  g_assert ( stat == NC_NOERR);

  it = 0;
  t1 = 0.;
  // Get time ?
  while ( it < ntdim && t1 < ctime ) {
    it++;
    stat = nc_get_vara_time (ncid, tvarid, it, 1, t1);
  }
  if ( t1 >= ctime && it > 1) {
    it--;
    stat = nc_get_vara_time (ncid, tvarid, it, 1, t1);
  }
  t2 = t1;

  // l 545: Get variable ids
  gint phiid, zid;
  stat = nc_inq_varid (ncid, "phi", &phiid);

  fprintf (stderr, "phi: %i \n", phiid);

  //stat = nc_inq_varid (ncid, "elevation", &zid);

  float m2[im][jm];
  // Build interpolation array / operator
  
  
  

  return 0;
}

void netcdf_build_spatial_interpolation_operators (Spline2D * splines,  NetCDFForcing * ncdf)
						   /* gint ncid, gint nx, gint ny) */
{ 
  // The interpolation points are the gauss-points
  gint i, j, k, ii, m, n;
  gint ng = splines->nouter;

  gint nx = ncdf->nx;
  gint ny = ncdf->ny;

  float xax[nx], yax[ny];

  for ( i = 0; i < nx; i++)
    xax[i] = *(ncdf->x + i);

  for ( i = 0; i < ny; i++)
    yax[i] = *(ncdf->y + i);

  gdouble x1, y1;

  Spline2D * sp = splines;
  while (sp) {
    for ( ii = 0; ii < sp->panels->len; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      GaussPoints * gp = spp->outer;

      if (gp->interpol->len == 0) {
      	for ( m = 0; m < ng*ng; m++) {
      	  Interpol * new = g_malloc (sizeof(Interpol));
      	  g_ptr_array_add (gp->interpol, new);
      	}
      }

      /* if (gp->interpol->len != 0) { */
      /* 	for ( m = 0; m < gp->interpol->len; m++) */
      /* 	  g_free (g_ptr_array_index (gp->interpol, m)); */
      /* 	g_ptr_array_free (gp->interpol, TRUE); */
      /* 	gp->interpol = g_ptr_array_new (); */
      /* } */


      // The GaussPoints are the interpolation points
      for ( n = 0; n < ng; n++ ) {
  	for ( m = 0; m < ng; m++ ) {
	  Interpol * interpol = g_ptr_array_index (gp->interpol, m + n*ng);
  	  Point ptmp = g_array_index (gp->Pi, Point, m + n*ng);

  	  // Get p in the inertial frame of reference
	
  	  // Here we should apply a transformation to get the
  	  // position of the point in the coordinate system
  	  // used by the NetCDF file
	
  	  // We can build the linear-interpolation operator
  	  // for point p

	  Point p;
	  p.x = ptmp.x /* + 2600. */;
	  p.y = ptmp.y /* + 2600. */;
	  p.z = 0.;
	
  	  // Find position of x
  	  gint ix = 0;
  	  while (xax[ix] < p.x) {
  	    ix++;
  	    if (ix > nx-1)
  	      break;
  	  }

  	  gint ix1 = MAX(0, ix-1);
  	  gint ix2 = MIN(ix, nx-1);

	  //fprintf (stdout, "ix1: %i ix2: %i \n", ix1, ix2);

  	  gdouble dxint=xax[ix2]-xax[ix1];
  	  // Find lever
  	  if (dxint > 0)
  	    x1 = (p.x-xax[ix1])/dxint;
  	  else
  	    x1 = 1.;

	  

  	  // Same for y
  	  gint iy = 0;
  	  while (yax[iy] < p.y) {
  	    iy++;
  	    if ( iy > (ny-1))
  	      break;
  	  }
  	  gint iy1 = MAX(0, iy-1);
  	  gint iy2 = MIN(iy-1, ny-1);

  	  gdouble dyint = yax[iy2] - yax[iy1];
  	  if (dyint > 0.)
  	    y1 = (p.y-yax[iy1])/dyint;
  	  else
  	    y1 = 1.;

	  /* Interpol * interpol = g_malloc (sizeof(Interpol)); */
  	  interpol->ix1 = ix1;
  	  interpol->ix2 = ix2;
  	  interpol->iy1 = iy1;
  	  interpol->iy2 = iy2;
  	  interpol->w11 = (1.-x1)*(1.-y1);
  	  interpol->w21 = x1*(1.-y1);
  	  interpol->w12 = (1.-x1)*y1;
  	  interpol->w22 = x1*y1;
	  /* g_ptr_array_add (gp->interpol, interpol); */
  	}
      }

    }
    sp = sp->next;
  }
  
}

NetCDFForcing * netcdf_forcing_new (gchar * file_name, gdouble xs, gdouble ys, gdouble angles)
{
  NetCDFForcing * new = g_malloc (sizeof (NetCDFForcing));
  new->it = 0;
  new->t1 = new->t2 = new->time = 0.;
  new->varp1 = -1;
  new->varp2 = -1;
  new->file_name = file_name;

  // Opens the netcdf file
  gint stat = nc_open (file_name, 0, &new->ncid);
  if ( stat != NC_NOERR ) {
    fprintf (stderr, "  *** Cannot open file: %s *** \n", file_name);
    g_assert_not_reached ();
  }

  // Get id of x variable
  stat = nc_inq_dimid (new->ncid, "longitude", &new->xid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "lon", &new->xid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "x", &new->xid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "X", &new->xid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "Easting", &new->xid);
  g_assert (stat == NC_NOERR);
  

  fprintf (stderr, "x: %i \n", new->xid);

  // Get id of y variable
  stat = nc_inq_dimid (new->ncid, "latitude", &new->yid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "lat", &new->yid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "y", &new->yid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "Y", &new->yid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "Northing", &new->yid);

  g_assert (stat == NC_NOERR);

  fprintf (stderr, "y: %i \n", new->yid);

  // Get id of time variable
  stat = nc_inq_dimid (new->ncid, "time", &new->tid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "Time", &new->tid);
  if (stat != NC_NOERR)
    stat = nc_inq_dimid (new->ncid, "ModelTime", &new->tid);  
  fprintf (stderr, "t: %i \n", new->tid);

  g_assert (stat == NC_NOERR);

  // Get size of x, y, time dimensions
  stat = nc_inq_dimlen (new->ncid, new->xid, &new->nx);
  stat = nc_inq_dimlen (new->ncid, new->yid, &new->ny);
  stat = nc_inq_dimlen (new->ncid, new->tid, &new->nt);

  fprintf (stderr, "nx ny nt: %i %i %i \n", (gint) new->nx, (gint) new->ny, (gint) new->nt);

  // Find position of time index
  gint tvarid;

  stat = nc_inq_varid (new->ncid, "longitude", &new->xvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "x", &new->xvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "lon", &new->xvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "X", &new->xvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "Easting", &new->xvarid);
  g_assert (stat == NC_NOERR);

  stat = nc_inq_varid(new->ncid, "latitude", &new->yvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "y", &new->yvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "lat", &new->yvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "Y", &new->yvarid);
  if (stat != NC_NOERR)
    stat = nc_inq_varid(new->ncid, "Northing", &new->yvarid);

  g_assert (stat == NC_NOERR);
  
  stat = nc_inq_varid (new->ncid, "time", &new->tvarid);
  stat = nc_inq_varid (new->ncid, "Time", &new->tvarid);
  stat = nc_inq_varid (new->ncid, "ModelTime", &new->tvarid);
  
  fprintf (stderr, "tvar: %i \n", new->tvarid);
  
  // Read all the values of x, y and t
  new->x = g_malloc (new->nx*sizeof(float));
  stat = nc_get_var_float (new->ncid, new->xvarid, new->x);
  new->y = g_malloc (new->ny*sizeof(float));
  stat = nc_get_var_float (new->ncid, new->yvarid, new->y);
  new->t = g_malloc (new->nt*sizeof(float));
  stat = nc_get_var_float (new->ncid, new->tvarid, new->t);
  new->tmax = *(new->t + new->nt - 1) + 1.;

  /* gint i; */
  /* for ( i = 0; i < new->ny; i++) { */
  /*   fprintf (stderr, "%i %f \n", i, *(new->y +i)); */
  /* } */

  // Get informations on frame of reference of data.
  float x0, y0, angle0;
  stat = nc_get_att_float (new->ncid, NC_GLOBAL, "x0", &x0);
  stat = nc_get_att_float (new->ncid, NC_GLOBAL, "y0", &y0);
  stat = nc_get_att_float (new->ncid, NC_GLOBAL, "angle", &angle0);

  angle0 *= M_PI/180.;
  
  /* gdouble xs = 1689540.7934; */
  /* gdouble ys  = 5676464.0678; */
  /* gdouble angles = 21.128*M_PI/180.; */
  

  new->x0 = (xs-x0)*cos(angle0) + (ys-y0)*sin(angle0);
  new->y0 = -(xs-x0)*sin(angle0) + (ys-y0)*cos(angle0);


  fprintf (stderr, "X0: %f Y0: %f  prop: %f %f \n", new->x0, new->y0, new->x0/597, new->y0/745);

  /* // Position of the center of the domain in Dave's coordinate system */
  /* new->x0 = (1689540.7934 - x0)*cos(angle/180*M_PI) - (5676464.0678 - y0)*sin(angle/180*M_PI); */
  /* new->y0 = (1689540.7934 - x0)*sin(angle/180*M_PI) + (5676464.0678 - y0)*cos(angle/180*M_PI); */

  /* new->x0 = ; */
  /* new->y0 = ; */

  fprintf (stderr, "x0 %f %f %f \n", x0, 1689540.7934, new->x0);
  fprintf (stderr, "y0 %f %f %f \n", y0, 5676464.0678, new->y0);
  fprintf (stderr, "angle %f %f %f \n", angle0, 21.128, angle0 - 21.128);

  // Angle of Dave's domain with respect to NZTM 2000 East axis
  // minus angle of the pier with respect to NZTM 2000 East axis
  //angle = (21.128-angle);
  

  angle0 = angles-angle0;

  fprintf (stderr, "ANGLE: %f %f \n", angle0, angle0*180./M_PI);

  new->calpha = cos(angle0);
  new->salpha = sin(angle0);

  fprintf (stderr, "x0 %f y0: %f Angle: %f \n", x0, y0, angle0);

  new->phi_1 = NULL;
  new->phi_2 = NULL;
  new->elevation_1 = NULL;
  new->elevation_2 = NULL;
  new->u_1 = NULL;
  new->u_2 = NULL;
  new->v_1 = NULL;
  new->v_2 = NULL;

  return new;
}

void netcdf_forcing_destroy (NetCDFForcing * ncf)
{
  if (ncf->phi_1)
    g_assert_not_reached ();
  // Needs to be finished

  g_free (ncf->x);
  g_free (ncf->y);
  g_free (ncf->t);
  g_free (ncf);
}

gdouble netcdf_fit_rhs_gauss (SPPanel * spp, gint m, gint n, gpointer data)
{
  NetCDFData * ncd = (NetCDFData *) data;
  gint ng = spp->sp->nouter;
  
  Interpol * ip = g_ptr_array_index (spp->outer->interpol, m + n*ng);
  gint ny = ncd->ny;
  gint nx = ncd->nx;
  float * val = ncd->val;

  /* fprintf (stdout, "ix1: %i ix2: %i iy1: %i iy2: %i \n", */
  /* 	   ip->ix1, ip->ix2, ip->iy1, ip->iy2); */
  /* fprintf (stdout, "%f %f %f %f \n", */
  /* 	   ip->w11, ip->w21, ip->w12, ip->w22); */
  /* fprintf (stdout, "%f %f %f %f\n", *(val + ip->ix1*ny + ip->iy1), *(val + ip->ix1*ny + ip->iy2), */
  /* 	   *(val + ip->ix2*ny + ip->iy2), *(val + ip->ix2*ny + ip->iy1)); */

  /* gint i, j; */
  /* for ( i = 0; i < ncdf->nx; i++) { */
  /*   for ( j = 0; j < ncdf->ny; j++) { */
  /*     *(ncd->val + i*n->ny + j) = *(ncdf->x + i); */
  /*   } */
  /* } */

  /* g_assert_not_reached (); */

  return ip->w11* *(val + ip->ix1 + ip->iy1*nx)
    + ip->w12* *(val + ip->ix1 + ip->iy2*nx)
    + ip->w22* *(val + ip->ix2 + ip->iy2*nx)
    + ip->w21* *(val + ip->ix2 + ip->iy1*nx);

  /* return ip->w11* *(val + ip->ix1*ny + ip->iy1) */
  /*   + ip->w12* *(val + ip->ix1*ny + ip->iy2) */
  /*   + ip->w22* *(val + ip->ix2*ny + ip->iy2) */
  /*   + ip->w21* *(val + ip->ix2*ny + ip->iy1); */
}

/* gint sort_data (gconstpointer a, gconstpointer b) */
/* { */
/*   if ( *(float *) a < *(float *) b) */
/*     return -1; */
/* } */

gsl_vector * netcdf_build_galerkin_rhs (Spline2D * sp, NetCDFForcing * ncdf, NetCDFData * ncd)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble RHS[size];
  gdouble x1, y1;

  gint ny = ncd->ny;
  gint nx = ncd->nx;
  float * val = ncd->val;
  
  float xax[nx], yax[ny];

  /* GTree * xtree = g_tree_new (sort_data); */
  /* for ( i = 0; i < nx; i++) */
  /*   g_tree_insert (xtree, , (ncdf->x + i)); */

  /* GTree * ytree = g_tree_new (sort_data); */
  /* for ( i = 0; i < ny; i++) { */
  /*   gint * value = g_malloc (sizeof(gint)); */
  /*   *value = i; */
  /*   g_tree_insert (ytree, (ncdf->y + i), value); */
  /* } */
  
  for ( i = 0; i < nx; i++)
    xax[i] = *(ncdf->x + i);

  for ( i = 0; i < ny; i++)
    yax[i] = *(ncdf->y + i);

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_matrix * Bu, * Bv;
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    ustart = gp->istart;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);


	Point ptmp = g_array_index (gp->Pi, Point, m + n*ng);

	// Get p in the inertial frame of reference
	
	// Here we should apply a transformation to get the
	// position of the point in the coordinate system
	// used by the NetCDF file
	
	// We can build the linear-interpolation operator
	// for point p

	Point p;
	p.x = 300 + ptmp.x*cos(M_PI/4) - ptmp.y*sin(M_PI/4) /* + 2600. */;
	p.y = 300 + ptmp.y*cos(M_PI/4) + ptmp.x*sin(M_PI/4) /* + 2600. */;
	p.z = 0.;

	/* p.x = ncdf->x0 + ncdf->calpha*ptmp.x - ncdf->salpha*ptmp.y; */
	/* p.y = ncdf->y0 + ncdf->salpha*ptmp.x + ncdf->calpha*ptmp.y; */
	/* p.z = 0.; */
	
	// Find position of x
	gint ix = 0;
	while (xax[ix] < p.x) {
	  ix++;
	  if (ix > nx-1)
	    break;
	}

	gint ix1 = MAX(0, ix-1);
	gint ix2 = MIN(ix, nx-1);

	//fprintf (stdout, "ix1: %i ix2: %i \n", ix1, ix2);

	gdouble dxint=xax[ix2]-xax[ix1];
	// Find lever
	if (dxint > 0)
	  x1 = (p.x-xax[ix1])/dxint;
	else
	  x1 = 1.;

	  

	// Same for y
	gint iy = 0;
	while (yax[iy] < p.y) {
	  iy++;
	  if ( iy > (ny-1))
	    break;
	}
	gint iy1 = MAX(0, iy-1);
	gint iy2 = MIN(iy-1, ny-1);

	gdouble dyint = yax[iy2] - yax[iy1];
	if (dyint > 0.)
	  y1 = (p.y-yax[iy1])/dyint;
	else
	  y1 = 1.;

	gdouble fmn = (1.-x1)*((1.-y1)* *(val + ix1 + iy1*nx)
			       + y1* *(val + ix1 + iy2*nx))
	  + x1*(y1* *(val + ix2 + iy2*nx)
		+ (1.-y1)* *(val + ix2 + iy1*nx));


	// Loop over the splines whose support is included in the panel
	for ( j = vstart; j < vstart + sp->k; j++) {
	  gdouble wmnj = wmn*gsl_matrix_get (Bv, j-vstart, 0);
	  gint indexj = j*NU;
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmnij = wmnj*gsl_matrix_get (Bu, i-ustart, 0);	 
	    RHS[indexj + i] += wmnij*fmn;
	  }
	}


      }
    }
  }

  // Copy lhs and rhs to gsl structures
  /* FILE * fp = fopen ("rhs-galerkin.tmp","w"); */
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs, i, RHS[i]);
    /* fprintf (fp, "%i %f \n", i, RHS[i]); */
  }
  /* fclose (fp); */

  return rhs;
}

gsl_vector * netcdf_build_galerkin_rhs_and_integrate_forces (Spline2D * sp,
							     NetCDFForcing * ncdf,
							     gdouble * forces_fk)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble RHS[size];
  gdouble x1, y1;

  gint ny = ncdf->ny;
  gint nx = ncdf->nx;

  gdouble t1 = ncdf->phi_1->t;
  gdouble t2 = ncdf->phi_2->t;
  gdouble dt = t2-t1;
  gdouble t = ncdf->time;

  /* fprintf (stderr, "Time %f \n", t); */
  
  /* gdouble forces_fk[6]; */
  /* for ( i = 0; i < 6; i++) */
  /*   forces_fk[6] = 0.; */

  float xax[nx], yax[ny];

  /* GTree * xtree = g_tree_new (sort_data); */
  /* for ( i = 0; i < nx; i++) */
  /*   g_tree_insert (xtree, , (ncdf->x + i)); */

  /* GTree * ytree = g_tree_new (sort_data); */
  /* for ( i = 0; i < ny; i++) { */
  /*   gint * value = g_malloc (sizeof(gint)); */
  /*   *value = i; */
  /*   g_tree_insert (ytree, (ncdf->y + i), value); */
  /* } */
  
  for ( i = 0; i < nx; i++)
    xax[i] = *(ncdf->x + i);

  for ( i = 0; i < ny; i++)
    yax[i] = *(ncdf->y + i);

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_matrix * Bu, * Bv;
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    ustart = gp->istart;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
	Point ptmp = g_array_index (gp->Pi, Point, m + n*ng);

	// Get p in the inertial frame of reference
	
	// Here we should apply a transformation to get the
	// position of the point in the coordinate system
	// used by the NetCDF file
	
	// We can build the linear-interpolation operator
	// for point p

	Point p;
	/* p.x = 300 + ptmp.x*cos(M_PI/4) - ptmp.y*sin(M_PI/4) /\* + 2600. *\/; */
	/* p.y = 300 + ptmp.y*cos(M_PI/4) + ptmp.x*sin(M_PI/4) /\* + 2600. *\/; */
	/* p.z = 0.; */

	p.x = 300 + ptmp.x;
	p.y = 300 + ptmp.y;
	p.z = 0.;

	/* p.x = ncdf->x0 + ncdf->calpha*ptmp.x - ncdf->salpha*ptmp.y; */
	/* p.y = ncdf->y0 + ncdf->salpha*ptmp.x + ncdf->calpha*ptmp.y; */
	/* p.z = 0.; */
	
	// Find position of x
	gint ix = 0;
	while (xax[ix] < p.x) {
	  ix++;
	  if (ix > nx-1)
	    break;
	}

	gint ix1 = MAX(0, ix-1);
	gint ix2 = MIN(ix, nx-1);

	//fprintf (stdout, "ix1: %i ix2: %i \n", ix1, ix2);

	gdouble dxint=xax[ix2]-xax[ix1];
	// Find lever
	if (dxint > 0)
	  x1 = (p.x-xax[ix1])/dxint;
	else
	  x1 = 1.;

	  

	// Same for y
	gint iy = 0;
	while (yax[iy] < p.y) {
	  iy++;
	  if ( iy > (ny-1))
	    break;
	}
	gint iy1 = MAX(0, iy-1);
	gint iy2 = MIN(iy-1, ny-1);

	gdouble dyint = yax[iy2] - yax[iy1];
	if (dyint > 0.)
	  y1 = (p.y-yax[iy1])/dyint;
	else
	  y1 = 1.;

	gdouble w11 = (1.-x1)*(1.-y1);
	gdouble w21 = x1*(1.-y1);
  	gdouble w12 = (1.-x1)*y1;
  	gdouble w22 = x1*y1;

	gdouble phi1 = w11* *(ncdf->phi_1->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->phi_1->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->phi_1->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->phi_1->val + ix2 + iy1*nx);

	gdouble phi2 = w11* *(ncdf->phi_2->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->phi_2->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->phi_2->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->phi_2->val + ix2 + iy1*nx);

	/* gdouble elevation1 = w11* *(ncdf->elevation_1->val + ix1 + iy1*nx) */
	/*   + w12 * *(ncdf->elevation_1->val + ix1 + iy2*nx) */
	/*   + w21 * *(ncdf->elevation_1->val + ix2 + iy2*nx) */
	/*   + w22 * *(ncdf->elevation_1->val + ix2 + iy1*nx); */

	/* gdouble elevation2 = w11* *(ncdf->elevation_2->val + ix1 + iy1*nx) */
	/*   + w12 * *(ncdf->elevation_2->val + ix1 + iy2*nx) */
	/*   + w21 * *(ncdf->elevation_2->val + ix2 + iy2*nx) */
	/*   + w22 * *(ncdf->elevation_2->val + ix2 + iy1*nx); */

	/* gdouble fmn = (1.-x1)*((1.-y1)* *(val + ix1 + iy1*nx) */
	/* 		       + y1* *(val + ix1 + iy2*nx)) */
	/*   + x1*(y1* *(val + ix2 + iy2*nx) */
	/* 	+ (1.-y1)* *(val + ix2 + iy1*nx)); */

	gdouble rho /* = sim->rho */;

	gdouble pressure = -rho*wmn*(phi2-phi1)/dt;

	*(forces_fk + 0) += pressure*N.x;
	*(forces_fk + 1) += pressure*N.y;
	*(forces_fk + 2) += pressure*N.z;

	Vector x;
	Point xg /* = sim->hull->xg */;
	x.x = ptmp.x-xg.x; x.y = ptmp.y-xg.y; x.z = ptmp.z-xg.z;
	x = vector_vector_product (&x, &N);
	
	*(forces_fk + 3) += pressure*x.x;
	*(forces_fk + 4) += pressure*x.y;
	*(forces_fk + 5) += pressure*x.z;

	/* gdouble fmn = ((t2-t)*elevation1 + (t-t1)*elevation2)/dt; */

	gdouble fmn = ((t2-t)*phi1 + (t-t1)*phi2)/dt;
	
	//gdouble fmn = (phi2-phi1)/dt;

	// Loop over the splines whose support is included in the panel
	for ( j = vstart; j < vstart + sp->k; j++) {
	  gdouble wmnj = wmn*gsl_matrix_get (Bv, j-vstart, 0);
	  gint indexj = j*NU;
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmnij = wmnj*gsl_matrix_get (Bu, i-ustart, 0);	 
	    RHS[indexj + i] += wmnij*fmn;
	  }
	}


      }
    }
  }

  // Copy lhs and rhs to gsl structures
  /* FILE * fp = fopen ("rhs-galerkin.tmp","w"); */
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs, i, RHS[i]);
    /* fprintf (fp, "%i %f \n", i, RHS[i]); */
  }
  /* fclose (fp); */

  return rhs;
}

gint find_index_dicho (float * y_data, gdouble y, gint ny)
{
  gint y1 = 0, y2 = ny-1, ym = (y2+y1)/2;

  if (y <= *(y_data))
    return 0;

  if (y >= *(y_data+ny-1))
    return ny-1;

  while (y2-y1 >= 2) {
    if ( y < *(y_data + ym))
      y2 = ym;
    else
      y1 = ym;

    ym = (y1+y2)/2;
  }

  return y2;
}

DoubleRhs * netcdf_build_double_galerkin_rhs (Spline2D * sp,
					      NetCDFForcing * ncdf,
					      Hull * hull)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  //gdouble RHS_PHI[size], RHS_ELEVATION[size];
  gdouble RHS_ELEVATION[size];
  gdouble x1, y1;

  gint ny = ncdf->ny;
  gint nx = ncdf->nx;

  gdouble t1 = ncdf->elevation_1->t;
  gdouble t2 = ncdf->elevation_2->t;
  gdouble dt = t2-t1;
  gdouble one_over_dt = 1./dt;
  gdouble t = ncdf->time;

  gdouble ramp = t < 200. ? t/200.:1.;

  gdouble a1 = (t2-t)*one_over_dt*ramp, a2 = (t-t1)*one_over_dt*ramp;

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    /* RHS_PHI[i] =  */RHS_ELEVATION[i] = 0.;
  }
  
  gsl_matrix * Bu, * Bv;
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    ustart = gp->istart;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
	Point pg = g_array_index (gp->Pi, Point, m + n*ng);

	// Get p in the inertial frame of reference
	Point p_inertial = hull_transformed_point (hull, pg);
	//p_inertial = pg;
	
	// Here we should apply a transformation to get the
	// position of the point in the coordinate system
	// used by the NetCDF file
	
	// We can build the linear-interpolation operator
	// for point p

	Point p;
	p.x = ncdf->x0 + ncdf->calpha*p_inertial.x - ncdf->salpha*p_inertial.y;
	p.y = ncdf->y0 + ncdf->salpha*p_inertial.x + ncdf->calpha*p_inertial.y;
	p.z = 0.;

	gint ix = find_index_dicho (ncdf->x, p.x, nx);

	gint ix1 = MAX(0, ix-1);
	gint ix2 = MIN(ix, nx-1);

	//gdouble dxint=xax[ix2]-xax[ix1];
	gdouble dxint=*(ncdf->x + ix2)-*(ncdf->x + ix1);
	// Find lever
	if (dxint > 0)
	  x1 = (p.x-*(ncdf->x + ix1))/dxint;
	else
	  x1 = 1.;

	  

	// Same for y
	gint iy = find_index_dicho (ncdf->y, p.y, ny);

	gint iy1 = MAX(0, iy-1);
	gint iy2 = MIN(iy-1, ny-1);

	gdouble dyint = *(ncdf->y + iy2) - *(ncdf->y + iy1);
	if (dyint > 0.)
	  y1 = (p.y-*(ncdf->y + iy1))/dyint;
	else
	  y1 = 1.;

	gdouble w11 = (1.-x1)*(1.-y1);
	gdouble w21 = x1*(1.-y1);
  	gdouble w12 = (1.-x1)*y1;
  	gdouble w22 = x1*y1;

	/* gdouble phi1 = w11* *(ncdf->phi_1->val + ix1 + iy1*nx) */
	/*   + w12 * *(ncdf->phi_1->val + ix1 + iy2*nx) */
	/*   + w21 * *(ncdf->phi_1->val + ix2 + iy2*nx) */
	/*   + w22 * *(ncdf->phi_1->val + ix2 + iy1*nx); */

	/* gdouble phi2 = w11* *(ncdf->phi_2->val + ix1 + iy1*nx) */
	/*   + w12 * *(ncdf->phi_2->val + ix1 + iy2*nx) */
	/*   + w21 * *(ncdf->phi_2->val + ix2 + iy2*nx) */
	/*   + w22 * *(ncdf->phi_2->val + ix2 + iy1*nx); */

	gdouble elevation1 = w11* *(ncdf->elevation_1->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->elevation_1->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->elevation_1->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->elevation_1->val + ix2 + iy1*nx);

	gdouble elevation2 = w11* *(ncdf->elevation_2->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->elevation_2->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->elevation_2->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->elevation_2->val + ix2 + iy1*nx);


	/* gdouble fmn_elevation = ((t2-t)*phi1 + (t-t1)*phi2)/dt; */
	gdouble fmn_elevation = a1*elevation1 + a2*elevation2;
	//gdouble fmn_phi = (phi2-phi1)*one_over_dt;

	// Loop over the splines whose support is included in the panel
	for ( j = vstart; j < vstart + sp->k; j++) {
	  gdouble wmnj = wmn*gsl_matrix_get (Bv, j-vstart, 0);
	  gint indexj = j*NU;
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmnij = wmnj*gsl_matrix_get (Bu, i-ustart, 0);	 
	    // RHS_PHI[indexj + i] += wmnij*fmn_phi;
	    RHS_ELEVATION[indexj + i] += wmnij*fmn_elevation;
	  }
	}


      }
    }
  }

  // Copy lhs and rhs to gsl structures
  /* FILE * fp = fopen ("rhs-galerkin.tmp","w"); */
  DoubleRhs * rhs = g_malloc (sizeof(DoubleRhs));
  //rhs->rhs_phi = gsl_vector_alloc (size);
  rhs->rhs_elevation = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++) {
    //gsl_vector_set (rhs->rhs_phi, i, RHS_PHI[i]);
    gsl_vector_set (rhs->rhs_elevation, i, RHS_ELEVATION[i]);
  }

  return rhs;
}

DoubleRhs * netcdf_build_linear_double_galerkin_rhs (Spline2D * sp,
						     NetCDFForcing * ncdf,
						     Hull * hull,
						     gdouble coeff)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  //gdouble RHS_PHI[size], RHS_ELEVATION[size];
  gdouble RHS_ELEVATION[size], RHS_U[size], RHS_V[size];
  gdouble x1, y1;

  gint ny = ncdf->ny;
  gint nx = ncdf->nx;

  gdouble t1 = ncdf->elevation_1->t;
  gdouble t2 = ncdf->elevation_2->t;
  gdouble dt = t2-t1;
  gdouble one_over_dt = 1./dt;
  gdouble t = ncdf->time;

  gdouble ramp = t < 200. ? t/200.:1.;

  gdouble a1 = (t2-t)*one_over_dt*ramp, a2 = (t-t1)*one_over_dt*ramp;

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    /* RHS_PHI[i] =  */RHS_ELEVATION[i] = RHS_U[i] = RHS_V[i] = 0.;
  }
  
  /* FILE * ftest = fopen ("bk.tmp","w"); */

  gsl_matrix * Bu, * Bv;
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    ustart = gp->istart;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
	Point pg = g_array_index (gp->Pi, Point, m + n*ng);
	
	// Here we should apply a transformation to get the
	// position of the point in the coordinate system
	// used by the NetCDF file
	
	// We can build the linear-interpolation operator
	// for point p

	Point p;
	p.x = ncdf->x0 + ncdf->calpha*pg.x - ncdf->salpha*pg.y;
	p.y = ncdf->y0 + ncdf->salpha*pg.x + ncdf->calpha*pg.y;
	p.z = 0.;

	/* fprintf (ftest, "%f %f %f\n", p.x, p.y, p.z); */

	gint ix = find_index_dicho (ncdf->x, p.x, nx);

	gint ix1 = MAX(0, ix-1);
	gint ix2 = MIN(ix, nx-1);

	//gdouble dxint=xax[ix2]-xax[ix1];
	gdouble dxint=*(ncdf->x + ix2)-*(ncdf->x + ix1);
	// Find lever
	if (dxint > 0)
	  x1 = (p.x-*(ncdf->x + ix1))/dxint;
	else
	  x1 = 1.;

	  

	// Same for y
	gint iy = find_index_dicho (ncdf->y, p.y, ny);

	gint iy1 = MAX(0, iy-1);
	gint iy2 = MIN(iy-1, ny-1);

	gdouble dyint = *(ncdf->y + iy2) - *(ncdf->y + iy1);
	if (dyint > 0.)
	  y1 = (p.y-*(ncdf->y + iy1))/dyint;
	else
	  y1 = 1.;

	gdouble w11 = (1.-x1)*(1.-y1);
	gdouble w21 = x1*(1.-y1);
  	gdouble w12 = (1.-x1)*y1;
  	gdouble w22 = x1*y1;

	/* gdouble phi1 = w11* *(ncdf->phi_1->val + ix1 + iy1*nx) */
	/*   + w12 * *(ncdf->phi_1->val + ix1 + iy2*nx) */
	/*   + w21 * *(ncdf->phi_1->val + ix2 + iy2*nx) */
	/*   + w22 * *(ncdf->phi_1->val + ix2 + iy1*nx); */

	/* gdouble phi2 = w11* *(ncdf->phi_2->val + ix1 + iy1*nx) */
	/*   + w12 * *(ncdf->phi_2->val + ix1 + iy2*nx) */
	/*   + w21 * *(ncdf->phi_2->val + ix2 + iy2*nx) */
	/*   + w22 * *(ncdf->phi_2->val + ix2 + iy1*nx); */

	gdouble elevation1 = w11* *(ncdf->elevation_1->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->elevation_1->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->elevation_1->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->elevation_1->val + ix2 + iy1*nx);

	gdouble elevation2 = w11* *(ncdf->elevation_2->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->elevation_2->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->elevation_2->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->elevation_2->val + ix2 + iy1*nx);

	gdouble u1 = w11* *(ncdf->u_1->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->u_1->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->u_1->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->u_1->val + ix2 + iy1*nx);

	gdouble u2 = w11* *(ncdf->u_2->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->u_2->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->u_2->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->u_2->val + ix2 + iy1*nx);

	gdouble v1 = w11* *(ncdf->v_1->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->v_1->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->v_1->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->v_1->val + ix2 + iy1*nx);

	gdouble v2 = w11* *(ncdf->v_2->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->v_2->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->v_2->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->v_2->val + ix2 + iy1*nx);

	/* gdouble fmn_elevation = ((t2-t)*phi1 + (t-t1)*phi2)/dt; */
	gdouble fmn_elevation = a1*elevation1 + a2*elevation2;
	gdouble fmn_u = a1*u1 + a2*u2;
	gdouble fmn_v = a1*v1 + a2*v2;
	//gdouble fmn_phi = (phi2-phi1)*one_over_dt;

	// Loop over the splines whose support is included in the panel
	for ( j = vstart; j < vstart + sp->k; j++) {
	  gdouble wmnj = wmn*gsl_matrix_get (Bv, j-vstart, 0);
	  gint indexj = j*NU;
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmnij = wmnj*gsl_matrix_get (Bu, i-ustart, 0);	 
	    // RHS_PHI[indexj + i] += wmnij*fmn_phi;
	    RHS_ELEVATION[indexj + i] += coeff*wmnij*fmn_elevation;
	    RHS_U[indexj + i] += coeff*wmnij*fmn_u;
	    RHS_V[indexj + i] += coeff*wmnij*fmn_v;
	  }
	}


      }
    }
  }

  /* fclose (ftest); */
  /* g_assert_not_reached (); */

  // Copy lhs and rhs to gsl structures
  /* FILE * fp = fopen ("rhs-galerkin.tmp","w"); */
  DoubleRhs * rhs = g_malloc (sizeof(DoubleRhs));
  //rhs->rhs_phi = gsl_vector_alloc (size);
  rhs->rhs_elevation = gsl_vector_alloc (size);
  rhs->rhs_u = gsl_vector_alloc (size);
  rhs->rhs_v = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++) {
    //gsl_vector_set (rhs->rhs_phi, i, RHS_PHI[i]);
    gsl_vector_set (rhs->rhs_elevation, i, RHS_ELEVATION[i]);
    gsl_vector_set (rhs->rhs_u, i, RHS_U[i]);
    gsl_vector_set (rhs->rhs_v, i, RHS_V[i]);
  }

  return rhs;
}

void netcdf_integrate_fk_forces (Spline2D * sp,
				 NetCDFForcing * ncdf,
				 gdouble * forces_fk,
				 Point xg,
				 gdouble rho)
{
  gint i, j, m, n, ii;
  gdouble x1, y1;

  gint ny = ncdf->ny;
  gint nx = ncdf->nx;

  gdouble t1 = ncdf->elevation_1->t;
  gdouble t2 = ncdf->elevation_2->t;
  gdouble dt = t2-t1;
  gdouble t = ncdf->time;

  float xax[nx], yax[ny];
  
  for ( i = 0; i < nx; i++)
    xax[i] = *(ncdf->x + i);

  for ( i = 0; i < ny; i++)
    yax[i] = *(ncdf->y + i);

  g_assert (sp != NULL);
  
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    for ( m = 0; m < ng; m++) {
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
	Point ptmp = g_array_index (gp->Pi, Point, m + n*ng);

	// Get p in the inertial frame of reference
	
	// Here we should apply a transformation to get the
	// position of the point in the coordinate system
	// used by the NetCDF file
	
	// We can build the linear-interpolation operator
	// for point p

	Point p;
	/* p.x = 300 + ptmp.x*cos(M_PI/4) - ptmp.y*sin(M_PI/4) /\* + 2600. *\/; */
	/* p.y = 300 + ptmp.y*cos(M_PI/4) + ptmp.x*sin(M_PI/4) /\* + 2600. *\/; */
	/* p.z = 0.; */

	p.x = 300 + ptmp.x;
	p.y = 300 + ptmp.y;
	p.z = 0.;

	/* p.x = ncdf->x0 + ncdf->calpha*ptmp.x - ncdf->salpha*ptmp.y; */
	/* p.y = ncdf->y0 + ncdf->salpha*ptmp.x + ncdf->calpha*ptmp.y; */
	/* p.z = 0.; */
	
	// Find position of x
	gint ix = 0;
	while (xax[ix] < p.x) {
	  ix++;
	  if (ix > nx-1)
	    break;
	}

	gint ix1 = MAX(0, ix-1);
	gint ix2 = MIN(ix, nx-1);

	//fprintf (stdout, "ix1: %i ix2: %i \n", ix1, ix2);

	gdouble dxint=xax[ix2]-xax[ix1];
	// Find lever
	if (dxint > 0)
	  x1 = (p.x-xax[ix1])/dxint;
	else
	  x1 = 1.;

	  
	// Same for y
	gint iy = 0;
	while (yax[iy] < p.y) {
	  iy++;
	  if ( iy > (ny-1))
	    break;
	}
	gint iy1 = MAX(0, iy-1);
	gint iy2 = MIN(iy-1, ny-1);

	gdouble dyint = yax[iy2] - yax[iy1];
	if (dyint > 0.)
	  y1 = (p.y-yax[iy1])/dyint;
	else
	  y1 = 1.;

	gdouble w11 = (1.-x1)*(1.-y1);
	gdouble w21 = x1*(1.-y1);
  	gdouble w12 = (1.-x1)*y1;
  	gdouble w22 = x1*y1;

	gdouble phi1 = w11* *(ncdf->phi_1->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->phi_1->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->phi_1->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->phi_1->val + ix2 + iy1*nx);

	gdouble phi2 = w11* *(ncdf->phi_2->val + ix1 + iy1*nx)
	  + w12 * *(ncdf->phi_2->val + ix1 + iy2*nx)
	  + w21 * *(ncdf->phi_2->val + ix2 + iy2*nx)
	  + w22 * *(ncdf->phi_2->val + ix2 + iy1*nx);

	gdouble pressure = rho*wmn*(phi2-phi1)/dt;

	*(forces_fk + 0) += pressure*N.x;
	*(forces_fk + 1) += pressure*N.y;
	*(forces_fk + 2) += pressure*N.z;

	Vector x;
	x.x = ptmp.x-xg.x; x.y = ptmp.y-xg.y; x.z = ptmp.z-xg.z;
	x = vector_vector_product (&x, &N);
	
	*(forces_fk + 3) += pressure*x.x;
	*(forces_fk + 4) += pressure*x.y;
	*(forces_fk + 5) += pressure*x.z;
      }
    }
  }
}

void netcdf_data_destroy (NetCDFData * ncd)
{
  g_free (ncd->val);
  g_free (ncd);
}

NetCDFData * netcdf_load_array (NetCDFForcing * ncf,
				NetCDFData * ncd,
				gchar * field_name,
				gint it)
{
  gint stat;
  size_t strt[3], icnt[3];

  if (ncd == NULL) {
    ncd = g_malloc (sizeof(NetCDFData));
    ncd->nx = ncf->nx;
    ncd->ny = ncf->ny;
    ncd->val = g_malloc (ncd->nx*ncd->ny*sizeof(float));

      // Find id of variable
    stat = nc_inq_varid (ncf->ncid, field_name, &ncd->varid);
    //fprintf (stderr, "%s : %i \n", field_name, ncd->varid);
    g_assert (stat == NC_NOERR);
  }

  ncd->it = it;
  ncd->t = *(ncf->t + it%ncf->nt -1) + ( ((int)it)/((int)ncf->nt))*ncf->nt;

  // Load the array of data
  strt[0] = it % ncf->nt;
  strt[1] = 0;
  strt[2] = 0;
  icnt[0] = 1;
  icnt[1] = ncd->ny;
  icnt[2] = ncd->nx;

  //fprintf (stdout, "%f %i %i\n", ncd->t, ncd->it, strt[0]);

  stat = nc_get_vara_float (ncf->ncid, ncd->varid, strt, icnt, ncd->val);

  //fprintf (stderr, "Time: %f it:%i \n", ncd->t, ncd->it);
  /* gint i, j; */
  /* FILE * fa = fopen("array.tmp","w"); */
  /* for ( i = 0; i < ncd->nx; i++) { */
  /*   for ( j = 0; j < ncd->ny; j++) { */
  /*     fprintf (/\* stdout *\/fa, "%f %f %f \n", *(ncf->x + i), *(ncf->y + j), *(ncd->val + i +j*ncf->nx)); */
  /*   } */
  /*   fprintf (fa, "\n"); */
  /* } */
  /* fclose (fa); */
  /* g_assert_not_reached (); */

  return ncd;
}

void load_netcdf_data (/* Simulation * sim */GSList * all_splines, NetCDFForcing * ncf, gchar * field_name, gdouble time)
{
  /* GSList * all_splines = NULL; */

  Spline2D * sp;

  gint stat;
  int ncid = ncf->ncid;

  ncf->it = 0;
  float t1 = 0., t2 = 0.;
  
  while ( ncf->it < ncf->nt && ncf->t1 < time ) {
    ncf->it++;
    t1 = t2;
    t2 = *(ncf->t + ncf->it);
  }
  ncf->it--;

  if ( t1 != ncf->t1 && t2 != ncf->t2) {
    // Need to load more data
    if (t1 == ncf->t2) {
      // One array of data to load
      // First copy data for t1 to t2

      // Load new data
      NetCDFData * ncd2 = netcdf_load_array (ncf, ncd2, "potential", ncf->it+1);
      
      // Fit of data
      GSList * splines = all_splines;
      while (splines) {
      	sp = splines->data;
	
	gsl_vector * gsl_rhs = sp->build_fit_rhs (sp,
						  netcdf_fit_rhs_gauss,
						  ncd2, NULL, NULL, gsl_rhs);
	ccs_problem_lu_solve (sp->fit, gsl_rhs);
	sp->copy_fit_solution (sp, gsl_rhs, ncf->varp2);
	gsl_vector_free (gsl_rhs);

	splines = splines->next;
      }
      
      netcdf_data_destroy (ncd2);
    }
    else {
      // Two arrays of data to load
      NetCDFData * ncd1 = netcdf_load_array (ncf, ncd1, "potential", ncf->it);

      gsl_vector * gsl_rhs = sp->build_fit_rhs (sp,
						netcdf_fit_rhs_gauss,
						ncd1, NULL, NULL,
						gsl_rhs);
      ccs_problem_lu_solve (sp->fit, gsl_rhs);
      sp->copy_fit_solution (sp, gsl_rhs, ncf->varp1);
      gsl_vector_free (gsl_rhs);

      /* netcdf_data_destroy (ncd1); */
      ncd1 = netcdf_load_array (ncf, ncd1, "potential", ncf->it+1);

      gsl_rhs = sp->build_fit_rhs (sp,
				   netcdf_fit_rhs_gauss,
				   ncd1, NULL, NULL,
				   gsl_rhs);
      ccs_problem_lu_solve (sp->fit, gsl_rhs);
      sp->copy_fit_solution (sp, gsl_rhs, ncf->varp2);
      gsl_vector_free (gsl_rhs);

      netcdf_data_destroy (ncd1);
    }
  }
  ncf->t1 = t1;
  ncf->t2 = t2;

  /* while ( ncf->it < ncf->nt && ncf->t1 < time ) { */
  /*   ncf->it++; */
  /*   stat = nc_get_vara_time (ncid, ncf->tvarid, ncf->it, 1, ncf->t1); */
  /* } */
  /* ncf->t2 = ncf->t1; */
  /* if ( ncf->t1 >= ncf->time && ncf->it > 1) { */
  /*   ncf->it--; */
  /*   stat = nc_get_vara_time (ncid, ncf->tvarid, ncf->it, 1, ncf->t1); */
  /* } */

  // Time interpolation will be between t1 and t2
  // i.e. it and it+1
  /* gint eid, pid; */
  /* stat = nc_inq_varid (ncid, "potential", &pid); */
  /* stat = nc_inq_varid (ncid, "elevation", &eid); */

  // Should be done somewhere else
  /* netcdf_build_spatial_interpolation_operators (sp, ncid, ncf->nx, ncf->ny); */

  

  // The fit has to be done at some point

  // We need something that makes sure that we are not between t1 and t2 anymore
  // and that reading is really required

}

static gdouble grid_x (SPPanel * spp, gint m, gint n, gpointer data)
{
  gdouble u = g_array_index (spp->outer->ui, gdouble, m);
  gdouble v = g_array_index (spp->outer->vj, gdouble, n);

  return /* 500 */400*(u-0.5);
}

static gdouble grid_y (SPPanel * spp, gint m, gint n, gpointer data)
{
  gdouble u = g_array_index (spp->outer->ui, gdouble, m);
  gdouble v = g_array_index (spp->outer->vj, gdouble, n);

  return 400*(v-0.5);
}



/* static void add_fh_fk_forces (Simulation * sim, Forces * f, gdouble t */
/* 			      , gdouble u[6], gdouble x[6]) */
/* { */
/*   gdouble ff[12]; */
/*   gint i; */
/*   for ( i = 0; i < 12; i++) */
/*     ff[i] = 0.; */

/*   wet_hull_integration (sim->hull, sim, &ff[0], */
/* 			local_fh_fk_forces_contribution, */
/* 			sppanel_fh_fk_forces_integral_gauss, */
/* 			pressure_force_integration_tolerance, */
/* 			sim->wp.wave_elevation, t, &sim->wp); */

/*   gdouble ramp = (1-exp(-t/40.)); */

/*   f->forces_h[0] += ff[0]; */
/*   f->forces_h[1] += ff[1]; */
/*   f->forces_h[2] += ff[2]; */

/*   f->forces_fk[0] += ff[3]*ramp; */
/*   f->forces_fk[1] += ff[4]*ramp; */
/*   f->forces_fk[2] += ff[5]*ramp; */


/*   f->forces_h[3] += ff[6]; */
/*   f->forces_h[4] += ff[7]; */
/*   f->forces_h[5] += ff[8]; */

/*   f->forces_fk[3] += ff[9]*ramp; */
/*   f->forces_fk[4] += ff[10]*ramp; */
/*   f->forces_fk[5] += ff[11]*ramp; */
/* } */

/* void add_netcdf_fk_forces (Simulation * sim, Forces * f, gdouble t, */
/* 			   gdouble u[6], gdouble x[6]) */
/* { */
/*   Hull * hull /\* = sim->hull *\/; */
/*   NetCDFForcing * ncdf /\* = sim->ncdf *\/; */
/*   Point xg /\* = sim->hull->xg *\/; */
/*   gdouble rho /\*= sim->rho*\/; */
  
/*   g_assert (ncdf != NULL); */
  
/*   ncdf->time = t; */
    
/*   // Find time position */
/*   gint it = 1; */
/*   gdouble t1 = 0.; */
/*   // Get time ? */
/*   while ( it < ncdf->nt && *(ncdf->t+it-1) < t ) */
/*     it++; */
    
/*   if ( *(ncdf->t+it-1) >= t && it > 1) */
/*     it--; */

/*   t1 = *(ncdf->t + it -1); */
/*   //fprintf (stderr, "t1: %f it: %i \n", t1, it); */

/*   /\* g_test_timer_start (); *\/ */
/*   // READ DATA IF REQUIRED */
/*   if (ncdf->phi_1 == NULL) { // No data has ever been loaded */
/*     ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it); */
/*     //ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it); */
      
/*     ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */
/*     //ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1); */
/*   } */
/*   else if (ncdf->phi_1->it != it) { */
/*     if (ncdf->phi_2->it == it) { // Can keep half of the data and move them */
/*       NetCDFData * tmp = ncdf->phi_1; */
/*       ncdf->phi_1 = ncdf->phi_2; */
/*       ncdf->phi_2 = tmp; */
/*       ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */

/*       /\* tmp = ncdf->elevation_1; *\/ */
/*       /\* ncdf->elevation_1 = ncdf->elevation_2; *\/ */
/*       /\* ncdf->elevation_2 = tmp; *\/ */
/*       /\* ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "u", it+1); *\/ */
/*     } */
/*     else { // Need to read two slices */
/*       ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it); */
/*       //ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it); */
      
/*       ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */
/*       //ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1); */
/*     } */
/*   } */
/*   g_assert (ncdf->phi_2->it == it+1); */
/*   g_assert (ncdf->phi_2 != NULL); */
/*   /\* fprintf (stderr, "LOAD DATA %f \n", g_test_timer_elapsed()); *\/ */


/*   GSList * patches = hull->patches; */

/*   gdouble forces[6]; */
/*   gint i; */
/*   for ( i = 0; i < 6; i++) */
/*     forces[i] = 0.; */
  
/*   for ( i = 0; i < 6; i++) */
/*     f->forces_fk[i] = 0.; */

/*   while ( patches ) { */
/*     Spline2D * sp = patches->data; */
/*     netcdf_integrate_fk_forces (sp, ncdf, f->forces_fk, xg, rho); */
/*     patches = patches->next; */
/*   } */

/*   // What we really want is the integral of dt phi on the surface and the value of the elevation */


/*   //gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs (grid, ncdf, ncdf->phi_1); */


/*   /\* gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs_and_integrate_forces (grid, ncdf); *\/ */
  

/*   /\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\/ */

/*   /\* /\\* g_test_timer_start (); *\\/ *\/ */
/*   /\* /\\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\\/ *\/ */
/*   /\* /\\* fprintf (stderr, "ccs_problem_lu_solve %f \n", g_test_timer_elapsed()); *\\/ *\/ */
      
/*   /\* g_test_timer_start (); *\/ */
/*   /\* grid->copy_fit_solution (grid, gsl_rhs, 12); *\/ */
/*   /\* fprintf (stderr, "copy_fit_solution %f \n", g_test_timer_elapsed()); *\/ */

/*   /\* gsl_vector_free (gsl_rhs); *\/ */

/* } */

/* void add_netcdf_fk_forces (Simulation * sim, Forces * f, gdouble t, */
/* 			   gdouble u[6], gdouble x[6]) */
/* { */
/*   Hull * hull /\* = sim->hull *\/; */
/*   NetCDFForcing * ncdf /\* = sim->ncdf *\/; */
/*   Point xg /\* = sim->hull->xg *\/; */
/*   gdouble rho /\*= sim->rho*\/; */
  
/*   g_assert (ncdf != NULL); */
  
/*   ncdf->time = t; */
    
/*   // Find time position */
/*   gint it = 1; */
/*   gdouble t1 = 0.; */
/*   // Get time ? */
/*   while ( it < ncdf->nt && *(ncdf->t+it-1) < t ) */
/*     it++; */
    
/*   if ( *(ncdf->t+it-1) >= t && it > 1) */
/*     it--; */

/*   t1 = *(ncdf->t + it -1); */
/*   //fprintf (stderr, "t1: %f it: %i \n", t1, it); */

/*   /\* g_test_timer_start (); *\/ */
/*   // READ DATA IF REQUIRED */
/*   if (ncdf->phi_1 == NULL) { // No data has ever been loaded */
/*     ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it); */
/*     //ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it); */
      
/*     ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */
/*     //ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1); */
/*   } */
/*   else if (ncdf->phi_1->it != it) { */
/*     if (ncdf->phi_2->it == it) { // Can keep half of the data and move them */
/*       NetCDFData * tmp = ncdf->phi_1; */
/*       ncdf->phi_1 = ncdf->phi_2; */
/*       ncdf->phi_2 = tmp; */
/*       ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */

/*       /\* tmp = ncdf->elevation_1; *\/ */
/*       /\* ncdf->elevation_1 = ncdf->elevation_2; *\/ */
/*       /\* ncdf->elevation_2 = tmp; *\/ */
/*       /\* ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "u", it+1); *\/ */
/*     } */
/*     else { // Need to read two slices */
/*       ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it); */
/*       //ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it); */
      
/*       ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */
/*       //ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1); */
/*     } */
/*   } */
/*   g_assert (ncdf->phi_2->it == it+1); */
/*   g_assert (ncdf->phi_2 != NULL); */
/*   /\* fprintf (stderr, "LOAD DATA %f \n", g_test_timer_elapsed()); *\/ */

/*   // Put the interpolated dt phi and elevation on the hull */
/*   GSList * patches = hull->patches; */
/*   while ( patches ) { */
/*     Spline2D * sp = patches->data; */

/*     DoubleRhs * rhs = netcdf_build_double_galerkin_rhs (sp, ncdf); */

/*     if (sp->fit == NULL) */
/*       sp->fit = sp->build_fit_matrix (sp); */

/*     ccs_problem_lu_solve (sp->fit, rhs->rhs_phi); */
/*     sp->copy_fit_solution (sp, rhs->rhs_phi, 30); */
    
/*     ccs_problem_lu_solve (sp->fit, rhs->rhs_elevation); */
/*     sp->copy_fit_solution (sp, rhs->rhs_elevation, 31); */

/*     gsl_vector_free (rhs->rhs_phi); */
/*     gsl_vector_free (rhs->rhs_elevation); */
/*     g_free (rhs); */

/*     patches = patches->next; */
/*   } */

/*   // Integration Froude-Krylov forces */



/*   /\* gdouble forces[6]; *\/ */
/*   /\* gint i; *\/ */
/*   /\* for ( i = 0; i < 6; i++) *\/ */
/*   /\*   forces[i] = 0.; *\/ */
  
/*   /\* for ( i = 0; i < 6; i++) *\/ */
/*   /\*   f->forces_fk[i] = 0.; *\/ */

/*   /\* while ( patches ) { *\/ */
/*   /\*   Spline2D * sp = patches->data; *\/ */
/*   /\*   netcdf_integrate_fk_forces (sp, ncdf, f->forces_fk, xg, rho); *\/ */
/*   /\*   patches = patches->next; *\/ */
/*   /\* } *\/ */

/*   // What we really want is the integral of dt phi on the surface and the value of the elevation */


/*   //gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs (grid, ncdf, ncdf->phi_1); */


/*   /\* gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs_and_integrate_forces (grid, ncdf); *\/ */
  

/*   /\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\/ */

/*   /\* /\\* g_test_timer_start (); *\\/ *\/ */
/*   /\* /\\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\\/ *\/ */
/*   /\* /\\* fprintf (stderr, "ccs_problem_lu_solve %f \n", g_test_timer_elapsed()); *\\/ *\/ */
      
/*   /\* g_test_timer_start (); *\/ */
/*   /\* grid->copy_fit_solution (grid, gsl_rhs, 12); *\/ */
/*   /\* fprintf (stderr, "copy_fit_solution %f \n", g_test_timer_elapsed()); *\/ */

/*   /\* gsl_vector_free (gsl_rhs); *\/ */

/* } */

/* void nc_test () */
/* { */
/*   gdouble t = 9.5; */

/*   //readSurffnc ("phi_test.nc", t); */

/*   g_test_timer_start (); */
/*   /\* NetCDFForcing * ncdf =  netcdf_forcing_new ("phi_2.nc"); *\/ */
/*   NetCDFForcing * ncdf =  netcdf_forcing_new ("phi_real_test.nc"); */
/*   fprintf (stderr, "netcdf_forcing_new %f \n", g_test_timer_elapsed()); */

/*   FILE * fin = fopen("1704deck-flat-mesh-12.2.GDF","r"); */
/*   Hull * hull = hull_new (); */
/*   hull_read (hull, fin, 25, 15, TRUE, TRUE, FALSE); */
/*   Spline2D * grid = hull->patches->data; */
/*   fclose (fin); */
/*   hull->patches = g_slist_append (hull->patches, */
/*   				  spline2d_symmetrical_y (hull->patches->data, 0)); */
/*   spline2d_translate  (hull->patches->data, -86., -16.-2., 0.); */
/*   spline2d_translate  (hull->patches->next->data , -86., -16.-2., 0.); */

/*   /\* hull->patches = g_slist_append (hull->patches, *\/ */
/*   /\* 				  parametric_grid (60, 60, grid_x, grid_y, NULL)); *\/ */

/*   /\* if (grid->fit == NULL) *\/ */
/*   /\*   grid->fit = grid->build_fit_matrix (grid); *\/ */

/*   GSList * patches = hull->patches; */
/*   while (patches) { */
/*     Spline2D * sp = patches->data; */
/*     if (sp->fit == NULL) */
/*       sp->fit = sp->build_fit_matrix (sp); */
/*     patches = patches->next; */
/*   } */

/*   /\* g_test_timer_start (); *\/ */
/*   /\* netcdf_build_spatial_interpolation_operators (grid, ncdf); *\/ */
/*   /\* fprintf (stderr, "netcdf_build_spatial_interpolation_operators %f \n", g_test_timer_elapsed()); *\/ */


/*   for ( t = 1.1; t < /\* 11.52 *\/1.2; t += 0.1) { */

/*     ncdf->time = t; */
    
/*     // Find time position */
/*     gint it = 1; */
/*     gdouble t1 = 0.; */
/*     // Get time ? */
/*     while ( it < ncdf->nt && *(ncdf->t+it-1) < t ) */
/*       it++; */
    
/*     if ( *(ncdf->t+it-1) >= t && it > 1) */
/*       it--; */

/*     t1 = *(ncdf->t + it -1); */
/*     //fprintf (stderr, "t1: %f it: %i \n", t1, it); */

/*     g_test_timer_start (); */
/*     // READ DATA IF REQUIRED */
/*     if (ncdf->phi_1 == NULL) { // No data has ever been loaded */
/*       ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it); */
/*       //ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it); */
      
/*       ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */
/*       //ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1); */
/*     } */
/*     else if (ncdf->phi_1->it != it) { */
/*       if (ncdf->phi_2->it == it) { // Can keep half of the data and move them */
/* 	NetCDFData * tmp = ncdf->phi_1; */
/* 	ncdf->phi_1 = ncdf->phi_2; */
/* 	ncdf->phi_2 = tmp; */
/* 	ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */

/* 	/\* tmp = ncdf->elevation_1; *\/ */
/* 	/\* ncdf->elevation_1 = ncdf->elevation_2; *\/ */
/* 	/\* ncdf->elevation_2 = tmp; *\/ */
/* 	/\* ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "u", it+1); *\/ */
/*       } */
/*       else { // Need to read two slices */
/* 	ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it); */
/* 	//ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it); */

/* 	ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1); */
/* 	//ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1); */
/*       } */
/*     } */
/*     g_assert (ncdf->phi_2->it == it+1); */
/*     g_assert (ncdf->phi_2 != NULL); */
/*     //fprintf (stderr, "LOAD DATA %f \n", g_test_timer_elapsed()); */


/*     // What we really want is the integral of dt phi on the surface and the value of the elevation */

/*     /\* g_test_timer_start (); *\/ */
/*     /\* NetCDFData * phi = NULL; *\/ */
/*     /\* phi = netcdf_load_array (ncdf, phi, "u", 80); *\/ */
/*     /\* fprintf (stderr, "netcdf_load_array %f \n", g_test_timer_elapsed()); *\/ */

/*     /\* Time of data *\/ */
/*     /\* fprintf (stderr, "time: %f \n", *(ncdf->t+80-1)); *\/ */

/*     /\* g_test_timer_start (); *\/ */
/*     /\* gsl_vector * gsl_rhs = grid->build_fit_rhs (grid, netcdf_fit_rhs_gauss, phi, NULL, NULL); *\/ */
/*     /\* fprintf (stderr, "grid->build_fit_rhs %f \n", g_test_timer_elapsed()); *\/ */

/*     g_test_timer_start (); */

/*     patches = hull->patches; */
/*     while (patches) { */
/*       Spline2D * sp = patches->data; */
    
/*       DoubleRhs * rhs = netcdf_build_double_galerkin_rhs (sp, ncdf); */
      
/*       if (sp->fit == NULL) */
/* 	sp->fit = sp->build_fit_matrix (sp); */
      
/*       ccs_problem_lu_solve (sp->fit, rhs->rhs_phi); */
/*       sp->copy_fit_solution (sp, rhs->rhs_phi, 30); */
      
/*       ccs_problem_lu_solve (sp->fit, rhs->rhs_elevation); */
/*       sp->copy_fit_solution (sp, rhs->rhs_elevation, 31); */
      
/*       gsl_vector_free (rhs->rhs_phi); */
/*       gsl_vector_free (rhs->rhs_elevation); */
/*       g_free (rhs); */
      
/*       patches = patches->next; */
/*     } */


/*     /\* DoubleRhs * rhs = netcdf_build_double_galerkin_rhs (grid, ncdf); *\/ */

/*     /\* if (grid->fit == NULL) *\/ */
/*     /\*   grid->fit = grid->build_fit_matrix (grid); *\/ */

/*     /\* ccs_problem_lu_solve (grid->fit, rhs->rhs_phi); *\/ */
/*     /\* grid->copy_fit_solution (grid, rhs->rhs_phi, 30); *\/ */
    
/*     /\* ccs_problem_lu_solve (grid->fit, rhs->rhs_elevation); *\/ */
/*     /\* grid->copy_fit_solution (grid, rhs->rhs_elevation, 31); *\/ */

/*     /\* gsl_vector_free (rhs->rhs_phi); *\/ */
/*     /\* gsl_vector_free (rhs->rhs_elevation); *\/ */
/*     /\* g_free (rhs); *\/ */



/*     /\* //gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs (grid, ncdf, ncdf->phi_1); *\/ */
/*     /\* gdouble forces_fk[6]; *\/ */
/*     /\* gint i; *\/ */
/*     /\* for ( i = 0; i < 6; i++) *\/ */
/*     /\*   forces_fk[6] = 0.; *\/ */

/*     /\* gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs_and_integrate_forces (grid, ncdf, &forces_fk[0]); *\/ */
/*     /\* fprintf (stderr, "netcdf_build_galerkin_rhs %f \n", g_test_timer_elapsed()); *\/ */

/*     /\* g_test_timer_start (); *\/ */
/*     /\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\/ */
/*     /\* fprintf (stderr, "ccs_problem_lu_solve %f \n", g_test_timer_elapsed()); *\/ */

/*     /\* /\\* g_test_timer_start (); *\\/ *\/ */
/*     /\* /\\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\\/ *\/ */
/*     /\* /\\* fprintf (stderr, "ccs_problem_lu_solve %f \n", g_test_timer_elapsed()); *\\/ *\/ */
      
/*     /\* g_test_timer_start (); *\/ */
/*     /\* grid->copy_fit_solution (grid, gsl_rhs, 12); *\/ */
/*     /\* fprintf (stderr, "copy_fit_solution %f \n", g_test_timer_elapsed()); *\/ */

/*     /\* gsl_vector_free (gsl_rhs); *\/ */



/*   } */

/*   FILE * fp = fopen ("grid.tmp","w"); */

/*   /\* spline2d_print_panels (grid, fp); *\/ */
/*   patches = hull->patches; */
/*   while (patches) { */
/*     Spline2D * sp = patches->data; */
    
/*     spline2d_print_panels (sp, fp); */

/*     patches = patches->next; */
/*   } */

/*   fclose (fp); */
  
/*   fp = fopen ("nc.tmp","w"); */

/*   patches = hull->patches; */
/*   while (patches) { */
/*     Spline2D * sp = patches->data; */
    
/*     spline2d_print_panels_gnuplot (sp, fp, 30); */

/*     patches = patches->next; */
/*   } */

/*   // spline2d_print_panels_gnuplot (grid, fp, 30); */
/*   fclose (fp); */

/*   patches = hull->patches; */
/*   while (patches) { */
/*     Spline2D * sp = patches->data; */
    
/*     spline2d_destroy (sp); */

/*     patches = patches->next; */
/*   } */

/*   /\* spline2d_destroy (grid); *\/ */

/*   netcdf_forcing_destroy (ncdf); */

/*   g_assert_not_reached (); */
/* } */

void nc_test ()
{
  gdouble t = 9.5;

  //readSurffnc ("phi_test.nc", t);

  g_test_timer_start ();
  /* NetCDFForcing * ncdf =  netcdf_forcing_new ("phi_2.nc"); */
  NetCDFForcing * ncdf =  netcdf_forcing_new ("phi_real_test.nc", 1689540.7934, 5676464.0678, 21.128*M_PI/180.);
  fprintf (stderr, "netcdf_forcing_new %f \n", g_test_timer_elapsed());

  FILE * fin = fopen("1704deck-flat-mesh-12.2.GDF","r");
  Hull * hull = hull_new ();
  hull_read (hull, fin, 25, 15, TRUE, FALSE, TRUE, FALSE);
  Spline2D * grid = hull->patches->data;
  fclose (fin);
  hull->patches = g_slist_append (hull->patches,
  				  spline2d_symmetrical_y (hull->patches->data, 0));
  spline2d_translate  (hull->patches->data, -86., -16.-2., 0.);
  spline2d_translate  (hull->patches->next->data , -86., -16.-2., 0.);

  /* hull->patches = g_slist_append (hull->patches, */
  /* 				  parametric_grid (60, 60, grid_x, grid_y, NULL)); */

  /* if (grid->fit == NULL) */
  /*   grid->fit = grid->build_fit_matrix (grid); */

  GSList * patches = hull->patches;
  while (patches) {
    Spline2D * sp = patches->data;
    if (sp->fit == NULL)
      sp->fit = sp->build_fit_matrix (sp);
    patches = patches->next;
  }

  /* g_test_timer_start (); */
  /* netcdf_build_spatial_interpolation_operators (grid, ncdf); */
  /* fprintf (stderr, "netcdf_build_spatial_interpolation_operators %f \n", g_test_timer_elapsed()); */


  for ( t = 1.1; t < /* 11.52 */1.2; t += 0.1) {

    ncdf->time = t;
    
    // Find time position
    gint it = 1;
    gdouble t1 = 0.;
    // Get time ?
    while ( it < ncdf->nt && *(ncdf->t+it-1) < t )
      it++;
    
    if ( *(ncdf->t+it-1) >= t && it > 1)
      it--;

    t1 = *(ncdf->t + it -1);
    //fprintf (stderr, "t1: %f it: %i \n", t1, it);

    g_test_timer_start ();
    // READ DATA IF REQUIRED
    if (ncdf->phi_1 == NULL) { // No data has ever been loaded
      ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it);
      //ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it);
      
      ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1);
      //ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1);
    }
    else if (ncdf->phi_1->it != it) {
      if (ncdf->phi_2->it == it) { // Can keep half of the data and move them
	NetCDFData * tmp = ncdf->phi_1;
	ncdf->phi_1 = ncdf->phi_2;
	ncdf->phi_2 = tmp;
	ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1);

	/* tmp = ncdf->elevation_1; */
	/* ncdf->elevation_1 = ncdf->elevation_2; */
	/* ncdf->elevation_2 = tmp; */
	/* ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "u", it+1); */
      }
      else { // Need to read two slices
	ncdf->phi_1 = netcdf_load_array (ncdf, ncdf->phi_1, "u", it);
	//ncdf->elevation_1 = netcdf_load_array (ncdf, ncdf->elevation_1, "elevation", it);

	ncdf->phi_2 = netcdf_load_array (ncdf, ncdf->phi_2, "u", it+1);
	//ncdf->elevation_2 = netcdf_load_array (ncdf, ncdf->elevation_2, "elevation", it+1);
      }
    }
    g_assert (ncdf->phi_2->it == it+1);
    g_assert (ncdf->phi_2 != NULL);
    //fprintf (stderr, "LOAD DATA %f \n", g_test_timer_elapsed());


    // What we really want is the integral of dt phi on the surface and the value of the elevation

    /* g_test_timer_start (); */
    /* NetCDFData * phi = NULL; */
    /* phi = netcdf_load_array (ncdf, phi, "u", 80); */
    /* fprintf (stderr, "netcdf_load_array %f \n", g_test_timer_elapsed()); */

    /* Time of data */
    /* fprintf (stderr, "time: %f \n", *(ncdf->t+80-1)); */

    /* g_test_timer_start (); */
    /* gsl_vector * gsl_rhs = grid->build_fit_rhs (grid, netcdf_fit_rhs_gauss, phi, NULL, NULL); */
    /* fprintf (stderr, "grid->build_fit_rhs %f \n", g_test_timer_elapsed()); */

    g_test_timer_start ();

    patches = hull->patches;
    while (patches) {
      Spline2D * sp = patches->data;
    
      DoubleRhs * rhs = netcdf_build_double_galerkin_rhs (sp, ncdf, hull);
      
      if (sp->fit == NULL)
	sp->fit = sp->build_fit_matrix (sp);
      
      ccs_problem_lu_solve (sp->fit, rhs->rhs_phi);
      sp->copy_fit_solution (sp, rhs->rhs_phi, 30);
      
      ccs_problem_lu_solve (sp->fit, rhs->rhs_elevation);
      sp->copy_fit_solution (sp, rhs->rhs_elevation, 31);
      
      gsl_vector_free (rhs->rhs_phi);
      gsl_vector_free (rhs->rhs_elevation);
      g_free (rhs);
      
      patches = patches->next;
    }


    /* DoubleRhs * rhs = netcdf_build_double_galerkin_rhs (grid, ncdf); */

    /* if (grid->fit == NULL) */
    /*   grid->fit = grid->build_fit_matrix (grid); */

    /* ccs_problem_lu_solve (grid->fit, rhs->rhs_phi); */
    /* grid->copy_fit_solution (grid, rhs->rhs_phi, 30); */
    
    /* ccs_problem_lu_solve (grid->fit, rhs->rhs_elevation); */
    /* grid->copy_fit_solution (grid, rhs->rhs_elevation, 31); */

    /* gsl_vector_free (rhs->rhs_phi); */
    /* gsl_vector_free (rhs->rhs_elevation); */
    /* g_free (rhs); */



    /* //gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs (grid, ncdf, ncdf->phi_1); */
    /* gdouble forces_fk[6]; */
    /* gint i; */
    /* for ( i = 0; i < 6; i++) */
    /*   forces_fk[6] = 0.; */

    /* gsl_vector * gsl_rhs = netcdf_build_galerkin_rhs_and_integrate_forces (grid, ncdf, &forces_fk[0]); */
    /* fprintf (stderr, "netcdf_build_galerkin_rhs %f \n", g_test_timer_elapsed()); */

    /* g_test_timer_start (); */
    /* ccs_problem_lu_solve (grid->fit, gsl_rhs); */
    /* fprintf (stderr, "ccs_problem_lu_solve %f \n", g_test_timer_elapsed()); */

    /* /\* g_test_timer_start (); *\/ */
    /* /\* ccs_problem_lu_solve (grid->fit, gsl_rhs); *\/ */
    /* /\* fprintf (stderr, "ccs_problem_lu_solve %f \n", g_test_timer_elapsed()); *\/ */
      
    /* g_test_timer_start (); */
    /* grid->copy_fit_solution (grid, gsl_rhs, 12); */
    /* fprintf (stderr, "copy_fit_solution %f \n", g_test_timer_elapsed()); */

    /* gsl_vector_free (gsl_rhs); */



  }

  FILE * fp = fopen ("grid.tmp","w");

  /* spline2d_print_panels (grid, fp); */
  patches = hull->patches;
  while (patches) {
    Spline2D * sp = patches->data;
    
    spline2d_print_panels (sp, fp);

    patches = patches->next;
  }

  fclose (fp);
  
  fp = fopen ("nc.tmp","w");

  patches = hull->patches;
  while (patches) {
    Spline2D * sp = patches->data;
    
    spline2d_print_panels_gnuplot (sp, fp, 30);

    patches = patches->next;
  }

  // spline2d_print_panels_gnuplot (grid, fp, 30);
  fclose (fp);

  patches = hull->patches;
  while (patches) {
    Spline2D * sp = patches->data;
    
    spline2d_destroy (sp);

    patches = patches->next;
  }

  /* spline2d_destroy (grid); */

  netcdf_forcing_destroy (ncdf);

  g_assert_not_reached ();
}


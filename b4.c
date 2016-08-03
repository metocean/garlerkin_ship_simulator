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


/* Higher-order spline */
static gdouble b4 (gdouble x, gdouble h)
{
  gdouble x0 = -5.*h/2.;
  gdouble x1 = -3.*h/2.;
  gdouble x2 = -h/2.;
  gdouble x3 = h/2.;
  gdouble x4 = 3.*h/2.;
  gdouble x5 = 5.*h/2.;

  if ( x <= -5.*h/2.)
    return 0.;
  if ( x > 5.*h/2.)
    return 0.;

  if ( -5.*h/2 < x && x <= -3.*h/2.)
    return 1./(24.*pow(h,4.))*pow(x-x0, 4.);

  if ( -3.*h/2 < x && x <= -h/2.)
    return 1./(24.*pow(h,4.)) * ((x2-x)*pow(x-x0,3.)
				 + pow(x-x0,2.)*(x-x1)*(x3-x)
				 + (x-x0)*pow(x-x1,2.)*(x4-x)
				 + pow(x-x1,3.)*(x5-x));

  if ( -h/2 < x && x <= h/2.)
    return 1./(24.*pow(h,4.)) * ( pow(x-x0,2.)*pow(x3-x,2.)
				  + (x-x0)*(x-x1)*(x3-x)*(x4-x)
				  + pow(x-x1,2.)*(x3-x)*(x5-x)
				  + (x-x0)*(x-x2)*pow(x4-x,2.)
				  + (x-x1)*(x-x2)*(x4-x)*(x5-x)
				  + pow(x-x2,2.)*pow(x5-x,2.) );

  if ( h/2 < x && x <= 3.*h/2.)
    return 1./(24.*pow(h,4.))*( (x-x0)*pow(x4-x,3.)
				+ (x-x1)*pow(x4-x, 2.)*(x5-x)
				+ (x-x2)*(x4-x)*pow(x5-x,2.)
				+ (x-x3)*pow(x5-x,3.) );

  if ( 3.*h/2 < x && x <= 5.*h/2.)
    return 1./(24.*pow(h,4.))*pow(x5-x, 4.);
  return 0.;
}

/**
 * Returns the contribution of a B4 panel
 **/
static gdouble b4_panel_eval (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b4 (j-p->j, 1.)*b4 (i-p->i, 1.);
}

/* Unidirectional bspline interpolation for border panels */
static  gdouble b4_panel_eval_y (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (i-p->i) < 0.5)
    return PANEL_VAL(p,var)*b4 (j-p->j, 1.);
  else
    return 0.;
}

static  gdouble b4_panel_eval_x (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (j-p->j) < 0.5)
    return PANEL_VAL(p,var)*b4 (i-p->i, 1.);
  else
    return 0.;
}

/* Returns TRUE if the point i, j belong to a border cell of patch */
static gboolean patch_border (Patch * patch, gdouble i, gdouble j)
{
  GArray * row = g_ptr_array_index (patch->rows, 0);

  if ( i < 1.5 || j < 1.5 || j > patch->rows->len-2.5 || i > row->len-2.5 )
    return TRUE;
  return FALSE;
}

static double patch_eval_y (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval_y (&panel, i, j, var);
    }
  }

  return sum;
}

static double patch_eval_x (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval_x (&panel, i, j, var);
    }
  }

  return sum;
}

/* High order 2D polynomial interpolation for the corners of the patch */
static gdouble interpolation (gdouble v[5][5], gdouble x, gdouble y)
{
  gdouble a0 = -1./6.*v[0][0] + 0.5*v[1][0] - 0.5*v[2][0] + 1./6.*v[3][0];
  gdouble a1 = -1./6.*v[0][0] + 0.5*v[0][1] - 0.5*v[0][2] + 1./6.*v[0][3];
  gdouble a2 = -0.5*v[0][0] + v[1][0]     -0.5*v[2][0] + 0.5*v[0][1] - v[1][1] + 0.5*v[2][1];
  gdouble a3 = -0.5*v[0][0] + 0.5*v[1][0] + v[0][1] - v[1][1] -0.5*v[0][2] + 0.5*v[1][2];
  gdouble a4 = v[0][0] - 2.5*v[1][0] + 2.0*v[2][0] - 0.5*v[3][0];
  gdouble a5 = v[0][0] -2.5*v[0][1] +2.*v[0][2] -0.5*v[0][3];
  gdouble a6 = 2.*v[0][0] - 2.5*v[1][0] + 0.5*v[2][0] -2.5*v[0][1] + 3.0*v[1][1] - 0.5*v[2][1] + 0.5*v[0][2] - 0.5*v[1][2];
  gdouble a7 = -11./6.*v[0][0] + 3.*v[1][0] -1.5*v[2][0] + 1./3.*v[3][0];
  gdouble a8 = -11./6.*v[0][0] + 3.*v[0][1] -1.5*v[0][2] +1./3.*v[0][3];
  gdouble a9 = v[0][0];

  return a0*x*x*x + a1*y*y*y + a2*x*x*y + a3*x*y*y + a4*x*x + a5*y*y + a6*x*y + a7*x + a8*y + a9;
}

/* High order 1D polynomial interpolation for the corners of the patch */
static gdouble interpolation_1D (gdouble v[5], gdouble x)
{
  gdouble a0 = v[0]/24. - v[1]/6. + 0.25*v[2] - v[3]/6. + v[4]/24.;
  gdouble a1 = -v[0]/2.4 + 1.5*v[1] - 2.*v[2] + 7./6.*v[3] - 0.25*v[4];
  gdouble a2 = 35./24.*v[0] - 13./3.*v[1] + 19./4*v[2] - 7./3.*v[3] + 11./24.*v[4];
  gdouble a3 = -25./12.*v[0] + 4.*v[1] - 3.*v[2] + 4./3.*v[3] - 0.25*v[4];
  gdouble a4 = v[0];

  return a0*pow(x,4.) + a1*pow(x,3.) + a2*pow(x,2.) + a3*x + a4;
}

/**
 * Evaluates the variable var at a given location (i,j) in the Patch p
 **/
gdouble b4_patch_eval (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;
  gint N = row->len;
  gint M = p->rows->len;

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-8;
  if ( j == -0.5 || j == 0.5)
    j += 1e-8;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-8;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-8;

  /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+2), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-2), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+2));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-2));

    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [5][5];
      gdouble dx = ((N-1) - i);
      gdouble dy = ((M-1) - j);
      
      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, N-1-l, M-1-k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(N-1-l),  rint(M-1-k), var);
	}
      
      return interpolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v [5][5];
      gdouble dx = ((N-1) - i);
      gdouble dy = j;

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, N-1-l, k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(N-1-l),  rint(k), var);
	}

      return interpolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v [5][5];
      gdouble dx = i;
      gdouble dy = j;

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, l, k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(l),  rint(k), var);
	}
      
      return interpolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v [5][5];
      gdouble dx = i;
      gdouble dy = ((M-1) - j);

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, l, M-1-k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(l),  rint(M-1-k), var);
	}

      return interpolation (v, dx, dy);
    }
    else {
      gdouble dx;
      gdouble v[5];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = ((N-1) - i);
	v[0] = patch_eval_y (p, N-1, j, var);
	v[1] = patch_eval_y (p, N-2, j, var);
	v[2] = patch_eval_y (p, N-3, j, var);
	v[3] = patch_eval_y (p, N-4, j, var);
	v[4] = patch_eval_y (p, N-5, j, var);
      }
      
      if (pl == NULL) {
	dx = i;
	v[0] = patch_eval_y (p, 0, j, var);
	v[1] = patch_eval_y (p, 1, j, var);
	v[2] = patch_eval_y (p, 2, j, var);
	v[3] = patch_eval_y (p, 3, j, var);
	v[4] = patch_eval_y (p, 4, j, var);
      }
      
      if (pt == NULL) {
	dx = (M-1 - j);
	v[0] = patch_eval_x (p, i, M-1, var);
	v[1] = patch_eval_x (p, i, M-2, var);
	v[2] = patch_eval_x (p, i, M-3, var);
	v[3] = patch_eval_x (p, i, M-4, var);
	v[4] = patch_eval_x (p, i, M-5, var);
      }

      if (pb == NULL) {
	dx = (j);
	v[0] = patch_eval_x (p, i, 0, var);
	v[1] = patch_eval_x (p, i, 1, var);
	v[2] = patch_eval_x (p, i, 2, var);
	v[3] = patch_eval_x (p, i, 3, var);
	v[4] = patch_eval_x (p, i, 4, var);
      }
      return interpolation_1D (v, dx);
    }
  }
  
  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval (&panel, i, j, var);
    }
  }

  return sum;
}

/* First derivative of b4 element */
static gdouble b4_x (gdouble x, gdouble h)
{
  gdouble x0 = -5.*h/2.;
  gdouble x1 = -3.*h/2.;
  gdouble x2 = -h/2.;
  gdouble x3 = h/2.;
  gdouble x4 = 3.*h/2.;
  gdouble x5 = 5.*h/2.;

  if ( x <= -5.*h/2.)
    return 0.;
  if ( x > 5.*h/2.)
    return 0.;

  if ( -5.*h/2 < x && x <= -3.*h/2.)
    return 1./(6.*pow(h,4.))*pow(x-x0, 3.);

  if ( -3.*h/2 < x && x <= -h/2.)
    return 1./(24.*pow(h,4.)) * (- pow(x-x0,3.) + 3.*(x2-x)*pow(x-x0,2.)
				 + 2.*(x-x0)*(x-x1)*(x3-x) + pow(x-x0,2.)*(x3-x) - pow(x-x0,2.)*(x-x1)
				 + pow(x-x1,2.)*(x4-x) + 2.*(x-x0)*(x-x1)*(x4-x) - (x-x0)*pow(x-x1,2.)
				 + 3.*pow(x-x1,2.)*(x5-x) - pow(x-x1,3.));

  if ( -h/2 < x && x <= h/2.)
     return 1./(24.*pow(h,4.)) * ( 2.*(x-x0)*pow(x3-x,2.) - 2.*(x3-x)*pow(x-x0,2.)
				   + (x-x1)*(x3-x)*(x4-x) + (x-x0)*(x3-x)*(x4-x)
				   - (x-x0)*(x-x1)*(x4-x) - (x-x0)*(x-x1)*(x3-x)
				   + 2*(x-x1)*(x3-x)*(x5-x) - pow(x-x1,2.)*(x5-x) - pow(x-x1,2.)*(x3-x)
				   + (x-x2)*pow(x4-x,2.) + (x-x0)*pow(x4-x,2.) - 2.*(x-x0)*(x-x2)*(x4-x)
				   +  (x-x2)*(x4-x)*(x5-x) + (x-x1)*(x4-x)*(x5-x)
				   - (x-x1)*(x-x2)*(x5-x) - (x-x1)*(x-x2)*(x4-x)
				   + 2.*(x-x2)*pow(x5-x,2.) - 2.*pow(x-x2,2.)*(x5-x) );

  if ( h/2 < x && x <= 3.*h/2.)
    return 1./(24.*pow(h,4.))*( pow(x4-x,3.) - 3.*(x-x0)*pow(x4-x,2.)
				+ pow(x4-x, 2.)*(x5-x) - 2.*(x-x1)*(x4-x)*(x5-x)
				- (x-x1)*pow(x4-x, 2.) + (x4-x)*pow(x5-x,2.)
				- (x-x2)*pow(x5-x,2.) - 2.*(x-x2)*(x4-x)*(x5-x)
				+ pow(x5-x,3.) - 3.*(x-x3)*pow(x5-x,2.)  );

  if ( 3.*h/2 < x && x <= 5.*h/2.)
    return -1./(6.*pow(h,4.))*pow(x5-x, 3.);
  return 0.;
}

static gdouble b4_panel_eval_dx (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b4 (j-p->j, 1.)*b4_x (i-p->i, 1.);
}

/* Unidirectional bspline interpolation for border panels */
static  gdouble b4_panel_eval_dx_x (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (j-p->j) < 0.5)
    return PANEL_VAL(p,var)*b4_x (i-p->i, 1.);
  else
    return 0.;
}

static gdouble interpolation_dx (gdouble v[5][5], gdouble x, gdouble y)
{
  gdouble a0 = -1./6.*v[0][0] + 0.5*v[1][0] - 0.5*v[2][0] + 1./6.*v[3][0];
  gdouble a2 = -0.5*v[0][0] + v[1][0]     -0.5*v[2][0] + 0.5*v[0][1] - v[1][1] + 0.5*v[2][1];
  gdouble a3 = -0.5*v[0][0] + 0.5*v[1][0] + v[0][1] - v[1][1] -0.5*v[0][2] + 0.5*v[1][2];
  gdouble a4 = v[0][0] - 2.5*v[1][0] + 2.0*v[2][0] - 0.5*v[3][0];
  gdouble a6 = 2.*v[0][0] - 2.5*v[1][0] + 0.5*v[2][0] -2.5*v[0][1] + 3.0*v[1][1]
    - 0.5*v[2][1] + 0.5*v[0][2] - 0.5*v[1][2];
  gdouble a7 = -11./6.*v[0][0] + 3.*v[1][0] -1.5*v[2][0] + 1./3.*v[3][0];

  return 3.*a0*x*x + 2.*a2*x*y + a3*y*y + 2.*a4*x + a6*y + a7;
}

static gdouble interpolation_1D_dx (gdouble v[5], gdouble x)
{
  gdouble a0 = v[0]/24. - v[1]/6. + 0.25*v[2] - v[3]/6. + v[4]/24.;
  gdouble a1 = -v[0]/2.4 + 1.5*v[1] - 2.*v[2] + 7./6.*v[3] - 0.25*v[4];
  gdouble a2 = 35./24.*v[0] - 13./3.*v[1] + 19./4*v[2] - 7./3.*v[3] + 11./24.*v[4];
  gdouble a3 = -25./12.*v[0] + 4.*v[1] - 3.*v[2] + 4./3.*v[3] - 0.25*v[4];

  return 4.*a0*pow(x,3.) + 3.*a1*pow(x,2.) + 2.*a2*x + a3;
}

static double patch_eval_dx_x (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval_dx_x (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Evaluates the x derivative of variable var at a given location (i,j) in the Patch p
 **/
gdouble b4_patch_eval_dx (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;
  gint N = row->len;
  gint M = p->rows->len;

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-8;
  if ( j == -0.5 || j == 0.5)
    j += 1e-8;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-8;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-8;

  /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+2), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-2), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+2));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-2));

    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [5][5];
      gdouble dx = ((N-1) - i);
      gdouble dy = ((M-1) - j);
      
      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, N-1-l, M-1-k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(N-1-l),  rint(M-1-k), var);
	}
      
      return -interpolation_dx (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v [5][5];
      gdouble dx = ((N-1) - i);
      gdouble dy = j;

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, N-1-l, k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(N-1-l),  rint(k), var);
	}

      return -interpolation_dx (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v [5][5];
      gdouble dx = i;
      gdouble dy = j;

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, l, k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(l),  rint(k), var);
	}
      
      return interpolation_dx (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v [5][5];
      gdouble dx = i;
      gdouble dy = ((M-1) - j);

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, l, M-1-k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(l),  rint(M-1-k), var);
	}

      return interpolation_dx (v, dx, dy);
    }
    else {
      gdouble dx;
      gdouble v[5];
      gdouble a0, a1, a2, a3, a4;
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = ((N-1) - i);
	v[0] = patch_eval_y (p, N-1, j, var);
	v[1] = patch_eval_y (p, N-2, j, var);
	v[2] = patch_eval_y (p, N-3, j, var);
	v[3] = patch_eval_y (p, N-4, j, var);
	v[4] = patch_eval_y (p, N-5, j, var);
	return -interpolation_1D_dx (v, dx);
      }
      
      if (pl == NULL) {
	dx = i;
	v[0] = patch_eval_y (p, 0, j, var);
	v[1] = patch_eval_y (p, 1, j, var);
	v[2] = patch_eval_y (p, 2, j, var);
	v[3] = patch_eval_y (p, 3, j, var);
	v[4] = patch_eval_y (p, 4, j, var);
	return interpolation_1D_dx (v, dx);
      }
      
      if (pt == NULL) {
	dx = (M-1 - j);
	v[0] = patch_eval_dx_x (p, i, M-1, var);
	v[1] = patch_eval_dx_x (p, i, M-2, var);
	v[2] = patch_eval_dx_x (p, i, M-3, var);
	v[3] = patch_eval_dx_x (p, i, M-4, var);
	v[4] = patch_eval_dx_x (p, i, M-5, var);
	return interpolation_1D (v, dx);
      }

      if (pb == NULL) {
	dx = (j);
	v[0] = patch_eval_dx_x (p, i, 0, var);
	v[1] = patch_eval_dx_x (p, i, 1, var);
	v[2] = patch_eval_dx_x (p, i, 2, var);
	v[3] = patch_eval_dx_x (p, i, 3, var);
	v[4] = patch_eval_dx_x (p, i, 4, var);
	return interpolation_1D (v, dx);
      }

      return a0*pow(dx,4.) + a1*pow(dx,3.) + a2*pow(dx,2.) + a3*dx + a4;
    }
  }
  
  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval_dx (&panel, i, j, var);
    }
  }

  return sum;
}

/* First derivative in the y direction */
static gdouble b4_panel_eval_dy (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b4_x (j-p->j, 1.)*b4 (i-p->i, 1.);
}

/* Unidirectional bspline interpolation for border panels */
static  gdouble b4_panel_eval_dy_y (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (i-p->i) < 0.5)
    return PANEL_VAL(p,var)*b4_x (j-p->j, 1.);
  else
    return 0.;
}

static gdouble interpolation_dy (gdouble v[5][5], gdouble x, gdouble y)
{
  gdouble a1 = -1./6.*v[0][0] + 0.5*v[0][1] - 0.5*v[0][2] + 1./6.*v[0][3];
  gdouble a2 = -0.5*v[0][0] + v[1][0]     -0.5*v[2][0] + 0.5*v[0][1] - v[1][1] + 0.5*v[2][1];
  gdouble a3 = -0.5*v[0][0] + 0.5*v[1][0] + v[0][1] - v[1][1] -0.5*v[0][2] + 0.5*v[1][2];
  gdouble a5 = v[0][0] -2.5*v[0][1] +2.*v[0][2] -0.5*v[0][3];
  gdouble a6 = 2.*v[0][0] - 2.5*v[1][0] + 0.5*v[2][0]
    -2.5*v[0][1] + 3.0*v[1][1] - 0.5*v[2][1] + 0.5*v[0][2] - 0.5*v[1][2];
  gdouble a8 = -11./6.*v[0][0] + 3.*v[0][1] -1.5*v[0][2] +1./3.*v[0][3];

  return 3.*a1*y*y + a2*x*x + 2.*a3*x*y + 2.*a5*y + a6*x + a8;
}

static double patch_eval_dy_y (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval_dy_y (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Evaluates the x derivative of variable var at a given location (i,j) in the Patch p
 **/
gdouble b4_patch_eval_dy (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;
  gint N = row->len;
  gint M = p->rows->len;

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-8;
  if ( j == -0.5 || j == 0.5)
    j += 1e-8;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-8;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-8;

  /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+2), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-2), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+2));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-2));

    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [5][5];
      gdouble dx = ((N-1) - i);
      gdouble dy = ((M-1) - j);
      
      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, N-1-l, M-1-k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(N-1-l),  rint(M-1-k), var);
	}
      
      return -interpolation_dy (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v [5][5];
      gdouble dx = ((N-1) - i);
      gdouble dy = j;

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, N-1-l, k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(N-1-l),  rint(k), var);
	}

      return interpolation_dy (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v [5][5];
      gdouble dx = i;
      gdouble dy = j;

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, l, k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(l),  rint(k), var);
	}
      
      return interpolation_dy (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v [5][5];
      gdouble dx = i;
      gdouble dy = ((M-1) - j);

      for (l = 0; l < 5 ; l++)
	for (k = 0; k < 5; k++) {
	  if ( k < 2 || l < 2 )
	    v[l][k] = PANEL_VAL(patch_panel_get (p, l, M-1-k), var);
	  else if ( k <= l )
	    v[l][k] = b4_patch_eval (p, rint(l),  rint(M-1-k), var);
	}

      return -interpolation_dy (v, dx, dy);
    }
    else {
      gdouble dx;
      gdouble v[5];
      gdouble a0, a1, a2, a3, a4;
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = ((N-1) - i);
	v[0] = patch_eval_dy_y (p, N-1, j, var);
	v[1] = patch_eval_dy_y (p, N-2, j, var);
	v[2] = patch_eval_dy_y (p, N-3, j, var);
	v[3] = patch_eval_dy_y (p, N-4, j, var);
	v[4] = patch_eval_dy_y (p, N-5, j, var);
	return interpolation_1D (v, dx);
      }
      
      if (pl == NULL) {
	dx = i;
	v[0] = patch_eval_dy_y (p, 0, j, var);
	v[1] = patch_eval_dy_y (p, 1, j, var);
	v[2] = patch_eval_dy_y (p, 2, j, var);
	v[3] = patch_eval_dy_y (p, 3, j, var);
	v[4] = patch_eval_dy_y (p, 4, j, var);
	return interpolation_1D (v, dx);
      }
      
      if (pt == NULL) {
	dx = (M-1 - j);
	v[0] = patch_eval_x (p, i, M-1, var);
	v[1] = patch_eval_x (p, i, M-2, var);
	v[2] = patch_eval_x (p, i, M-3, var);
	v[3] = patch_eval_x (p, i, M-4, var);
	v[4] = patch_eval_x (p, i, M-5, var);
	return -interpolation_1D_dx (v, dx);
      }

      if (pb == NULL) {
	dx = (j);
	v[0] = patch_eval_x (p, i, 0, var);
	v[1] = patch_eval_x (p, i, 1, var);
	v[2] = patch_eval_x (p, i, 2, var);
	v[3] = patch_eval_x (p, i, 3, var);
	v[4] = patch_eval_x (p, i, 4, var);
	a0 = v[0]/24. - v[1]/6. + 0.25*v[2] - v[3]/6. + v[4]/24.;
	a1 = -v[0]/2.4 + 1.5*v[1] - 2.*v[2] + 7./6.*v[3] - 0.25*v[4];
	a2 = 35./24.*v[0] - 13./3.*v[1] + 19./4*v[2] - 7./3.*v[3] + 11./24.*v[4];
	a3 = -25./12.*v[0] + 4.*v[1] - 3.*v[2] + 4./3.*v[3] - 0.25*v[4];
	a4 = v[0];
	return interpolation_1D_dx (v, dx);
      }

      return a0*pow(dx,4.) + a1*pow(dx,3.) + a2*pow(dx,2.) + a3*dx + a4;
    }
  }
  
  for ( k = MAX (0, j-3); k <= MIN (p->rows->len-1, j+3) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-3); l <= MIN (row->len-1, i+3); l++) {
      panel = g_array_index (row, Panel, l);
      sum += b4_panel_eval_dy (&panel, i, j, var);
    }
  }

  return sum;
}

static gdouble b4_2D_bspline (gdouble i, gdouble j, gdouble x, gdouble y)
{
  return b4 (x-i, 1.)*b4 (y-j, 1.);
}

/** Initial surface fitting */
/**
 * Finds the best patch fit using the corner values of the panels.
 * The center values of the panels are determined that way.
 **/
void b4_patch_fit_panels (Patch * patch)
{
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  gint N = row->len + 1;
  gint M = patch->rows->len + 1;
  Point p[N][M]; /* Points to fit */
  gdouble px[N][M], py[N][M], pz[N][M];

  patch_prefit_panels (patch);

  Point p0;
  p0.x = p0.y = p0.z = 0.;

  for ( i = 0; i < N; i++)
    for ( j = 0; j < M; j++)
      p[i][j] = p0;
  
  /* Copy data to p */
  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      p[i][j] = panel.p[2];
    }
  }

  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    i = row->len-1;
    Panel panel = g_array_index (row, Panel, i);
    p[N-1][j] = panel.p[1];
  }

  j = patch->rows->len-1;
  row = g_ptr_array_index (patch->rows, j);
  for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      p[i][M-1] = panel.p[3];
  }

  i = row->len-1;
  Panel panel = g_array_index (row, Panel, i);
  p[N-1][M-1] = panel.p[0];

  
  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *w, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;
  gint n = (N-4)*(M-4);

  cx = gsl_vector_alloc (n);
  cy = gsl_vector_alloc (n);
  cz = gsl_vector_alloc (n);
  cov_x = gsl_matrix_alloc (n, n);
  cov_y = gsl_matrix_alloc (n, n);
  cov_z = gsl_matrix_alloc (n, n);

  X = gsl_matrix_alloc (n, n);
  x = gsl_vector_alloc (n);
  y = gsl_vector_alloc (n);
  z = gsl_vector_alloc (n);

  /* Initial guess */
  for ( i = 2; i < N-2; i++)
    for ( j = 2; j < M-2; j++) {
      gsl_vector_set (cx, (i-2)+(j-2)*(N-4), (p[i][j].x+p[i+1][j].x+p[i][j+1].x+p[i+1][j+1].x)/4.);
      gsl_vector_set (cy, (i-2)+(j-2)*(N-4), (p[i][j].y+p[i+1][j].y+p[i][j+1].y+p[i+1][j+1].y)/4.);
      gsl_vector_set (cz, (i-2)+(j-2)*(N-4), (p[i][j].z+p[i+1][j].z+p[i][j+1].z+p[i+1][j+1].z)/4.);
    }

  gint k, l;
  for ( i = 2; i < N-2; i++)
    for ( j = 2; j < M-2; j++) {
      for ( k = 2; k < N-2; k++)
  	for ( l = 2; l < M-2; l++) {
  	  gsl_matrix_set (X, (k-2)+(l-2)*(N-4), (i-2)+(j-2)*(N-4),
			  b4_2D_bspline ((gdouble) i, (gdouble) j, (gdouble) k, (gdouble) l));
  	}
    }
  
  for ( k = 2; k < N-2; k++)
    for ( l = 2; l < M-2; l++) {
      gdouble xi = p[k][l].x;
      gdouble yi = p[k][l].y;
      gdouble zi = p[k][l].z;

      for ( i = 0; i < N; i++) {
      	xi -= p[i][0].x*b4_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	xi -= p[i][M-1].x*b4_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
  	xi -= p[i][1].x*b4_2D_bspline ((gdouble) i, (gdouble) 1, (gdouble) k, (gdouble) l);
      	xi -= p[i][M-2].x*b4_2D_bspline ((gdouble) i, (gdouble) M-2, (gdouble) k, (gdouble) l);

  	yi -= p[i][0].y*b4_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	yi -= p[i][M-1].y*b4_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
  	yi -= p[i][1].y*b4_2D_bspline ((gdouble) i, (gdouble) 1, (gdouble) k, (gdouble) l);
      	yi -= p[i][M-2].y*b4_2D_bspline ((gdouble) i, (gdouble) M-2, (gdouble) k, (gdouble) l);

  	zi -= p[i][0].z*b4_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	zi -= p[i][M-1].z*b4_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
  	zi -= p[i][1].z*b4_2D_bspline ((gdouble) i, (gdouble) 1, (gdouble) k, (gdouble) l);
      	zi -= p[i][M-2].z*b4_2D_bspline ((gdouble) i, (gdouble) M-2, (gdouble) k, (gdouble) l);
      }

      for ( j = 0; j < M; j++) {
      	xi -= p[0][j].x*b4_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	xi -= p[N-1][j].x*b4_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
  	xi -= p[1][j].x*b4_2D_bspline ((gdouble) 1, (gdouble) j, (gdouble) k, (gdouble) l);
      	xi -= p[N-2][j].x*b4_2D_bspline ((gdouble) N-2, (gdouble) j, (gdouble) k, (gdouble) l);

  	yi -= p[0][j].y*b4_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	yi -= p[N-1][j].y*b4_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
  	yi -= p[1][j].y*b4_2D_bspline ((gdouble) 1, (gdouble) j, (gdouble) k, (gdouble) l);
      	yi -= p[N-2][j].y*b4_2D_bspline ((gdouble) N-2, (gdouble) j, (gdouble) k, (gdouble) l);

  	zi -= p[0][j].z*b4_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	zi -= p[N-1][j].z*b4_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
  	zi -= p[1][j].z*b4_2D_bspline ((gdouble) 1, (gdouble) j, (gdouble) k, (gdouble) l);
      	zi -= p[N-2][j].z*b4_2D_bspline ((gdouble) N-2, (gdouble) j, (gdouble) k, (gdouble) l);
      }

      xi += p[0][0].x*b4_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      xi += p[0][M-1].x*b4_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      xi += p[N-1][M-1].x*b4_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      xi += p[N-1][0].x*b4_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      xi += p[1][1].x*b4_2D_bspline ((gdouble) 1, (gdouble) 1, (gdouble) k, (gdouble) l);
      xi += p[1][M-2].x*b4_2D_bspline ((gdouble) 1, (gdouble) M-2, (gdouble) k, (gdouble) l);
      xi += p[N-2][M-2].x*b4_2D_bspline ((gdouble) N-2, (gdouble) M-2, (gdouble) k, (gdouble) l);
      xi += p[N-2][1].x*b4_2D_bspline ((gdouble) N-2, (gdouble) 1, (gdouble) k, (gdouble) l);

      yi += p[0][0].y*b4_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      yi += p[0][M-1].y*b4_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      yi += p[N-1][M-1].y*b4_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      yi += p[N-1][0].y*b4_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      yi += p[1][1].y*b4_2D_bspline ((gdouble) 1, (gdouble) 1, (gdouble) k, (gdouble) l);
      yi += p[1][M-2].y*b4_2D_bspline ((gdouble) 1, (gdouble) M-2, (gdouble) k, (gdouble) l);
      yi += p[N-2][M-2].y*b4_2D_bspline ((gdouble) N-2, (gdouble) M-2, (gdouble) k, (gdouble) l);
      yi += p[N-2][1].y*b4_2D_bspline ((gdouble) N-2, (gdouble) 1, (gdouble) k, (gdouble) l);

      zi += p[0][0].z*b4_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      zi += p[0][M-1].z*b4_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      zi += p[N-1][M-1].z*b4_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      zi += p[N-1][0].z*b4_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      zi += p[1][1].z*b4_2D_bspline ((gdouble) 1, (gdouble) 1, (gdouble) k, (gdouble) l);
      zi += p[1][M-2].z*b4_2D_bspline ((gdouble) 1, (gdouble) M-2, (gdouble) k, (gdouble) l);
      zi += p[N-2][M-2].z*b4_2D_bspline ((gdouble) N-2, (gdouble) M-2, (gdouble) k, (gdouble) l);
      zi += p[N-2][1].z*b4_2D_bspline ((gdouble) N-2, (gdouble) 1, (gdouble) k, (gdouble) l);

      gsl_vector_set (x, (k-2)+(l-2)*(N-4), xi);
      gsl_vector_set (y, (k-2)+(l-2)*(N-4), yi);
      gsl_vector_set (z, (k-2)+(l-2)*(N-4), zi);
    }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, n);

  /* Multidimensional fitting usind singular values decomposition method */
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, 1e-6, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, 1e-6, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, 1e-6, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b4 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 2; i < N-2; i++)
    for ( j = 2; j < M-2; j++) {
      p[i][j].x = gsl_vector_get(cx, i-2+(j-2)*(N-4));
      p[i][j].y = gsl_vector_get(cy, i-2+(j-2)*(N-4));
      p[i][j].z = gsl_vector_get(cz, i-2+(j-2)*(N-4));
    }
  
  gsl_matrix_free (X);
  gsl_vector_free (x);
  gsl_vector_free (y);
  gsl_vector_free (z);
  gsl_vector_free (cx);
  gsl_vector_free (cy);
  gsl_vector_free (cz);
  gsl_matrix_free (cov_x);
  gsl_matrix_free (cov_y);
  gsl_matrix_free (cov_z);	 

  /* Interpolate the spline to the center points of the panels */
  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);

      if ( j > 1 && j < patch->rows->len-2 && i > 1 && row->len-2 > i ) {
	panel.var[0] = 0.;
	panel.var[1] = 0.;
	panel.var[2] = 0.;
	for ( k = 0; k < N; k++)
	  for ( l = 0; l < M; l++) {
	    panel.var[0] += p[k][l].x
	      *b4_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
	    panel.var[1] += p[k][l].y
	      *b4_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
	    panel.var[2] += p[k][l].z
	      *b4_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
	  }
      }
      panel.p[0] = p[i+1][j+1];
      panel.p[1] = p[i+1][j];
      panel.p[2] = p[i][j];
      panel.p[3] = p[i][j+1];
      g_array_index (row, Panel, i) = panel;
    }
  }

  FILE * fp = fopen("fit.tmp","w");
  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      fprintf(fp, "%g %g %g \n", panel.var[0], panel.var[1], panel.var[2]);
    }
  }
  fclose (fp);
}


Point b4_patch_eval_point (Patch * patch, gdouble i, gdouble j)
{
  Point p;

  p.x = b4_patch_eval (patch, i, j, 0);
  p.y = b4_patch_eval (patch, i, j, 1);
  p.z = b4_patch_eval (patch, i, j, 2);

  return p;
}

Point b4_patch_eval_transformed_point (Patch * patch, gdouble i, gdouble j)
{
  Point p;

  g_assert (patch->t != NULL);

  p.x = b4_patch_eval (patch, i, j, 0);
  p.y = b4_patch_eval (patch, i, j, 1);
  p.z = b4_patch_eval (patch, i, j, 2);

  return transform_point (p, &patch->t->xg, &patch->t->t, &patch->t->euler_m);
}

/**
 * Return a vector containing the normal at the surface at point (i,j)
 * The normal is expressed in the (x,y,z) coordinate system.
 * The normal is oriented towards the inside of the hull.
 * The returned normal is non-normalised.
 **/
Vector b4_patch_normal (Patch * p, gdouble i, gdouble j)
{
  gdouble xi = b4_patch_eval_dx (p,i,j,0);
  gdouble xj = b4_patch_eval_dy (p,i,j,0);
  gdouble yi = b4_patch_eval_dx (p,i,j,1);
  gdouble yj = b4_patch_eval_dy (p,i,j,1);
  gdouble zi = b4_patch_eval_dx (p, i, j, 2);
  gdouble zj = b4_patch_eval_dy (p, i, j, 2);

  gdouble norm = 0., J = 0.;
  Vector n, vn = panel_first_order_normal (patch_panel_get (p, i, j));

  /* return panel_first_order_normal (patch_panel_get (p, i, j)); */

  if (fabs(vn.x) > fabs(vn.y) && fabs(vn.x) > fabs(vn.z)) {
    J = yi*zj-yj*zi;
    n.x = 1.;
    n.y = -1./J*(xi*zj-xj*zi);
    n.z = -1./J*(-xi*yj+xj*yi);
  }
  else if (fabs(vn.y) > fabs(vn.x) && fabs(vn.y) > fabs(vn.z)) {
    J = zi*xj-zj*xi;
    n.x = -1./J*(-yi*zj+yj*zi);
    n.y = 1.;
    n.z = -1./J*(yi*xj-yj*xi);
  }
  else {
    J = xi*yj-xj*yi;
    n.x = -1./J*(zi*yj-zj*yi);
    n.y = -1./J*(-zi*xj+zj*xi);
    n.z = 1.;
  }

  g_assert (J != 0.);

  /* If the normal to the spline panel does not point in the same direction as the normal  */
  /* to the flat panel, then we flip it over                                               */
  if ( vector_scalar_product (&vn, &n) < 0 ) {
    n.x *= -1.;
    n.y *= -1.;
    n.z *= -1.;
  }

  return n;
}

static Vector dx (Patch * p, gdouble i, gdouble j)
{
  gdouble xi = b4_patch_eval_dx (p,i,j,0);
  gdouble xj = b4_patch_eval_dy (p,i,j,0);
  gdouble yi = b4_patch_eval_dx (p,i,j,1);
  gdouble yj = b4_patch_eval_dy (p,i,j,1);
  gdouble zi = b4_patch_eval_dx (p, i, j, 2);
  gdouble zj = b4_patch_eval_dy (p, i, j, 2);

  Vector dx;

  double Jxy = xi*yj-xj*yi;
  double Jyz = yi*zj-yj*zi;
  double Jzx = zi*xj-zj*xi;

  fprintf(stderr, "%f %f | %f %f | %f %f \n", yj/Jxy, -zi/Jzx, zj/Jyz, -xi/Jxy, xj/Jzx, -yi/Jyz);

  return dx;
}

Vector b4_patch_unit_normal (Patch * p, gdouble i, gdouble j)
{
  return vector_normalise (b4_patch_normal (p, i, j));
}

gdouble b4_patch_local_metric (Patch * p, gdouble i, gdouble j)
{
  gdouble xi = b4_patch_eval_dx (p,i,j,0);
  gdouble xj = b4_patch_eval_dy (p,i,j,0);
  gdouble yi = b4_patch_eval_dx (p,i,j,1);
  gdouble yj = b4_patch_eval_dy (p,i,j,1);
  gdouble zi = b4_patch_eval_dx (p,i,j,2);
  gdouble zj = b4_patch_eval_dy (p,i,j,2);
  Vector g;

  g.x = yi*zj-yj*zi;
  g.y = zi*xj-zj*xi;
  g.z = xi*yj-xj*yi;

  return vector_norm (g);
}

void b4_patch_normal_print (Patch * p, FILE * fp)
{
  gdouble i, j;
  GArray * row = g_ptr_array_index (p->rows, 0);

  for ( j = -0.5 ; j <= p->rows->len-0.5 ; j += 0.5) {
    for ( i = -0.5; i <= row->len-0.5; i += 0.5) {
      Vector n = b4_patch_unit_normal (p, i, j);
      Point p0 = b4_patch_eval_point (p, i, j);
      fprintf(fp, "%g %g %g \n %g %g %g \n\n\n", p0.x, p0.y, p0.z, p0.x+n.x, p0.y+n.y, p0.z+n.z);
    }
  }
}

void b4_patch_store_derivatives (Patch * patch)
{
  gint i, j;
  GArray * row = g_ptr_array_index (patch->rows, 0);
  GArray * row0, * row2;
  
  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 1; i < row->len-1; i++) {
      Panel p = g_array_index (row, Panel, i);

      p.var[3] = b4_patch_eval_dx (patch, p.i, p.j, 0);
      p.var[4] = b4_patch_eval_dy (patch, p.i, p.j, 0);
      p.var[5] = b4_patch_eval_dx (patch, p.i, p.j, 1);
      p.var[6] = b4_patch_eval_dy (patch, p.i, p.j, 1);
      p.var[7] = b4_patch_eval_dx (patch, p.i, p.j, 2);
      p.var[8] = b4_patch_eval_dy (patch, p.i, p.j, 2);
    }
  }
}


/**********************************************************************/

static gboolean is_under_freesurface (Patch * patch, HeightCurve hz,
				      gdouble il, gdouble jl,
				      gdouble ih, gdouble jh)
{
  /* if (hz != NULL) { */
  /*   gdouble dz[4]; */
  /*   dz[0] = b4_patch_eval (patch, il, jl, 2) */
  /*     - hz (b4_patch_eval (patch, il, jl, 0), b4_patch_eval (patch, il, jl, 1)); */
  /*   dz[1] = b4_patch_eval (patch, ih, jl, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jl, 0), b4_patch_eval (patch, ih, jl, 1)); */
  /*   dz[2] = b4_patch_eval (patch, ih, jh, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jh, 0), b4_patch_eval (patch, ih, jh, 1)); */
  /*   dz[3] = patch_eval (patch, il, jh, 2) */
  /*     - hz (b4_patch_eval (patch, il, jh, 0), b4_patch_eval (patch, il, jh, 1)); */
    
  /*   if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.)) */
  /*     return FALSE; */
  /* } */
  return TRUE;
}

static GSList * intersection (Patch * patch, HeightCurve hz, gdouble il,
			      gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 1e-2;
  GSList * l = NULL;

  /* Check for intersection */
  gdouble dz[4];
  Point c[4];
  gint i;
  c[0] = b4_patch_eval_transformed_point (patch, il, jl);
  c[1] = b4_patch_eval_transformed_point (patch, ih, jl);
  c[2] = b4_patch_eval_transformed_point (patch, ih, jh);
  c[3] = b4_patch_eval_transformed_point (patch, il, jh);

  /* for ( i = 0; i < 4; i++) */
  /*   dz[i] = c[i].z - hz (c[i].x, c[i].y); */

  if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.) ||
       (dz[0] < 0. && dz[1] < 0. && dz[2] < 0. && dz[3] < 0.) )
    return l;

  /* We refine */
  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  if (fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance) {
    Point * p = g_malloc (sizeof(Point));
    *p = b4_patch_eval_transformed_point (patch, ihalf, jhalf);
    l = g_slist_append (l, p);
    return l;
  }

  /* l = g_slist_concat (l, intersection (patch, hz, il, ihalf, jl, jhalf)); */
  /* l = g_slist_concat (l, intersection (patch, hz, ihalf, ih, jl, jhalf)); */
  /* l = g_slist_concat (l, intersection (patch, hz, il, ihalf, jhalf, jh)); */
  /* l = g_slist_concat (l, intersection (patch, hz, ihalf, ih, jhalf, jh)); */
  return NULL;
}

GSList * b4_patch_intersect_with_free_surface (Patch * patch, HeightCurve hz)
{
  GSList * l = NULL;
  gint i, j;

  /* for ( j = 0; j < patch->rows->len; j++) { */
  /*   GArray * row = g_ptr_array_index (patch->rows, j); */
  /*   for ( i = 0; i < row->len; i++) { */
  /*     Panel panel = g_array_index (row, Panel, i); */
  /*     l = g_slist_concat (l, intersection (patch, hz, panel.i - 0.5, panel.i + 0.5, panel.j - 0.5, panel.j + 0.5)); */
  /*   } */
  /* } */
  return l;
}

/**
 * Returns the value of the integral of a given variable times the local unit normal
 * var over the bspline surface of the patch.
 **/
static Vector adaptive_panel_integral (Patch * patch, HeightCurve hz, gint var,
				       gdouble il, gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 0.0625;
  /* gdouble border_tolerance = 0.001; */
  Vector v;

  /* Initialize to zero */
  v.x = v.y = v.z = 0.;

  if (!is_under_freesurface (patch, hz, il, jl, ih, jh))
    return v;

  /* if (hz != NULL) { */
  /*   /\* Check if we are under the free-surface *\/ */
  /*   gdouble dz[4]; */
  /*   dz[0] = b4_patch_eval (patch, il, jl, 2) */
  /*     - hz (b4_patch_eval (patch, il, jl, 0), b4_patch_eval (patch, il, jl, 1)); */
  /*   dz[1] = b4_patch_eval (patch, ih, jl, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jl, 0), b4_patch_eval (patch, ih, jl, 1)); */
  /*   dz[2] = b4_patch_eval (patch, ih, jh, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jh, 0), b4_patch_eval (patch, ih, jh, 1)); */
  /*   dz[3] = patch_eval (patch, il, jh, 2) */
  /*     - hz (b4_patch_eval (patch, il, jh, 0), b4_patch_eval (patch, il, jh, 1)); */
    
  /*   if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.)) */
  /*     return v; */
  /* } */

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Unit normal at the center of the small piece of panel */
    Vector n = b4_patch_normal (patch, (il+ih)/2., (jl+jh)/2.);
    
    /* gdouble J = vector_norm (n); */
    /* n = vector_normalise (v); */

    gdouble val = patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)*fabs(ih-il)*fabs(jh-jl);

    v.x = val*n.x;
    v.y = val*n.y;
    v.z = val*n.z;
    
    return v;
  }

  /* We refine */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, il, ihalf, jl, jhalf)); */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, ihalf, ih, jl, jhalf)); */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, il, ihalf, jhalf, jh)); */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, ihalf, ih, jhalf, jh)); */

  return v;
}

Vector b4_patch_adaptive_panel_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  Vector sum;
  sum.x = sum.y = sum.z = 0.;

  /* for ( j = 0; j < patch->rows->len ; j++) { */
  /*   GArray * row = g_ptr_array_index (patch->rows, j); */
  /*   for ( i = 0; i < row->len; i++) { */
  /*     Panel p = g_array_index (row, Panel, i); */
  /*     sum = vector_sum (sum, adaptive_panel_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5)); */
  /*   } */
  /* } */
  
  if (patch->t != NULL)
    return transform_vector (sum, patch->t);
  else
    return sum;
}

/**
 * Returns the value of the integral of a given variable times the local unit normal
 * var over the bspline surface of the patch.
 **/
static gdouble adaptive_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var,
					       gdouble il, gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 1./* 0.625 */;
  /* gdouble border_tolerance = 0.001; */
  gdouble v = 0.;

  if (!is_under_freesurface (patch, hz, il, jl, ih, jh))
    return v;
  
  /* if (hz != NULL) { */
  /*   /\* Check if we are under the free-surface *\/ */
  /*   gdouble dz[4]; */
  /*   dz[0] = b4_patch_eval (patch, il, jl, 2) */
  /*     - hz (b4_patch_eval (patch, il, jl, 0), b4_patch_eval (patch, il, jl, 1)); */
  /*   dz[1] = b4_patch_eval (patch, ih, jl, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jl, 0), b4_patch_eval (patch, ih, jl, 1)); */
  /*   dz[2] = b4_patch_eval (patch, ih, jh, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jh, 0), b4_patch_eval (patch, ih, jh, 1)); */
  /*   dz[3] = b4_patch_eval (patch, il, jh, 2) */
  /*     - hz (b4_patch_eval (patch, il, jh, 0), b4_patch_eval (patch, il, jh, 1)); */
    
  /*   if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.)) */
  /*     return v; */
  /* } */

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Surface integral of var over the small piece of panel */
    return /* patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)* */
      b4_patch_local_metric (patch, (il+ih)/2., (jl+jh)/2.)*fabs(ih-il)*fabs(jh-jl);
  }

  /* We refine */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jl, jhalf); */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jl, jhalf); */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jhalf, jh); */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jhalf, jh); */

  return v;
}

gdouble b4_patch_adaptive_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  gdouble sum = 0.;

  /* for ( j = 0; j < patch->rows->len ; j++) { */
  /*   GArray * row = g_ptr_array_index (patch->rows, j); */
  /*   for ( i = 0; i < row->len; i++) { */
  /*     Panel p = g_array_index (row, Panel, i); */
  /*     sum += adaptive_panel_scalar_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5); */
  /*   } */
  /* } */
  return sum;
}

static void panel_compute_metric (Panel * panel, gpointer data)
{
  Patch * patch = (Patch *) data;
  gdouble xi = b4_patch_eval_dx (patch, panel->i, panel->j, 0);
  gdouble xj = b4_patch_eval_dy (patch, panel->i, panel->j, 0);
  gdouble yi = b4_patch_eval_dx (patch, panel->i, panel->j, 1);
  gdouble yj = b4_patch_eval_dy (patch, panel->i, panel->j, 1);
  gdouble zi = b4_patch_eval_dx (patch, panel->i, panel->j, 2);
  gdouble zj = b4_patch_eval_dy (patch, panel->i, panel->j, 2);
  
  Vector g;
  g.x = yi*zj-yj*zi;
  g.y = zi*xj-zj*xi;
  g.z = xi*yj-xj*yi;

  PANEL_VAL (panel, 9) = vector_norm (g);
}

void b4_patch_compute_metric (Patch * patch)
{
  patch_forall_panels (patch, panel_compute_metric, patch);
}

/**
 * Return the integral of a scalar over a patch.
 * This assume that the function patch_compute_metric has been run 
 * beforehand.
 **/
static gdouble fast_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var,
					   gdouble il, gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 1./* 0.0625 */;
  /* gdouble border_tolerance = 0.001; */
  gdouble v = 0.;

  if (!is_under_freesurface (patch, hz, il, jl, ih, jh))
    return v;
  /* if (hz != NULL) { */
  /*   /\* Check if we are under the free-surface *\/ */
  /*   gdouble dz[4]; */
  /*   dz[0] = patch_eval (patch, il, jl, 2) */
  /*     - hz (patch_eval (patch, il, jl, 0), patch_eval (patch, il, jl, 1)); */
  /*   dz[1] = patch_eval (patch, ih, jl, 2) */
  /*     - hz (patch_eval (patch, ih, jl, 0), patch_eval (patch, ih, jl, 1)); */
  /*   dz[2] = patch_eval (patch, ih, jh, 2) */
  /*     - hz (patch_eval (patch, ih, jh, 0), patch_eval (patch, ih, jh, 1)); */
  /*   dz[3] = patch_eval (patch, il, jh, 2) */
  /*     - hz (patch_eval (patch, il, jh, 0), patch_eval (patch, il, jh, 1)); */
    
  /*   if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.)) */
  /*     return v; */
  /* } */

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Surface integral of var over the small piece of panel */
    /* var=3 == metric */

    return b4_patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)
      *b4_patch_eval (patch, (il+ih)/2., (jl+jh)/2., 3)*fabs(ih-il)*fabs(jh-jl);
  }

  /* /\* We refine *\/ */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jl, jhalf); */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jl, jhalf); */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jhalf, jh); */
  /* v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jhalf, jh); */

  return v;
}

gdouble b4_patch_fast_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  gdouble sum = 0.;

  /* for ( j = 0; j < patch->rows->len ; j++) { */
  /*   GArray * row = g_ptr_array_index (patch->rows, j); */
  /*   for ( i = 0; i < row->len; i++) { */
  /*     Panel p = g_array_index (row, Panel, i); */
  /*     sum += fast_panel_scalar_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5); */
  /*   } */
  /* } */
  return sum;
}

static Vector b4_patch_local_gradient (Patch * patch, gdouble i,
				       gdouble j, gint var)
{
  Vector grad;

  b4_patch_eval_dx (patch, i, j, var);

  /* grad.x = patch_eval (); */

  return grad;
}

/**
 * Returns the value of the integral of a given variable times the local unit normal
 * var over the bspline surface of the patch.
 **/
static Vector adaptive_panel_flux_integral (Patch * patch, HeightCurve hz, gint var,
					    gdouble il, gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 0.0625;
  /* gdouble border_tolerance = 0.001; */
  Vector v;

  /* Initialize to zero */
  v.x = v.y = v.z = 0.;

  if (!is_under_freesurface (patch, hz, il, jl, ih, jh))
    return v;

  /* if (hz != NULL) { */
  /*   /\* Check if we are under the free-surface *\/ */
  /*   gdouble dz[4]; */
  /*   dz[0] = b4_patch_eval (patch, il, jl, 2) */
  /*     - hz (b4_patch_eval (patch, il, jl, 0), b4_patch_eval (patch, il, jl, 1)); */
  /*   dz[1] = b4_patch_eval (patch, ih, jl, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jl, 0), b4_patch_eval (patch, ih, jl, 1)); */
  /*   dz[2] = b4_patch_eval (patch, ih, jh, 2) */
  /*     - hz (b4_patch_eval (patch, ih, jh, 0), b4_patch_eval (patch, ih, jh, 1)); */
  /*   dz[3] = patch_eval (patch, il, jh, 2) */
  /*     - hz (b4_patch_eval (patch, il, jh, 0), b4_patch_eval (patch, il, jh, 1)); */
    
  /*   if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.)) */
  /*     return v; */
  /* } */

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Unit normal at the center of the small piece of panel */
    Vector n = b4_patch_unit_normal (patch, (il+ih)/2., (jl+jh)/2.);
    Vector grad;
  
    
    /* gdouble J = vector_norm (n); */
    /* n = vector_normalise (v); */

    gdouble val = patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)*fabs(ih-il)*fabs(jh-jl);

    v.x = val*n.x;
    v.y = val*n.y;
    v.z = val*n.z;
    
    return v;
  }

  /* /\* We refine *\/ */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, il, ihalf, jl, jhalf)); */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, ihalf, ih, jl, jhalf)); */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, il, ihalf, jhalf, jh)); */
  /* v = vector_sum (v, adaptive_panel_integral (patch, hz, var, ihalf, ih, jhalf, jh)); */

  return v;
}

Vector b4_patch_adaptive_panel_flux_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  Vector sum;
  sum.x = sum.y = sum.z = 0.;

  /* for ( j = 0; j < patch->rows->len ; j++) { */
  /*   GArray * row = g_ptr_array_index (patch->rows, j); */
  /*   for ( i = 0; i < row->len; i++) { */
  /*     Panel p = g_array_index (row, Panel, i); */
  /*     sum = vector_sum (sum, adaptive_panel_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5)); */
  /*   } */
  /* } */
  
  if (patch->t != NULL)
    return transform_vector (sum, patch->t);
  else
    return sum;
}

/* static gdouble b4_single_source_integral (Patch * patch, HeightCurve hz, */
/* 					  gint i, gint j, gint var) */
/* { */
/*   gdouble sum = 0.; */
/*   gint k, l; */
/*   GArray * row; */

/*   for ( l = MAX (0, j-2); l < MIN (j+3, patch->rows->len); l++) { */
/*     row = g_array_ptr_array (patch->rows, l); */
/*     for ( k = MAX (0, i-2); k < MIN (i+3, row->len); k++) { */
/*       Panel p = g_array_index (row, Panel, k); */
/*       /\* NOT  REALLY TRUE *\/ */
/*       sum += b4_panel_eval (p, i, j, var); */
/*     } */
/*   } */
/*   return sum; */
/* } */



typedef struct {
  Patch * patch;
  HeightCurve hz;
} VTMParams;

static void panel_compute_volume_times_metric (Panel * panel, gpointer data)
{
  VTMParams * par = (VTMParams *) data;

  /* if (par->hz == NULL) */
  /*   PANEL_VAL (panel, 9) = PANEL_VAL (panel, 9); */
  /* else */
  /*   PANEL_VAL (panel, 10) = adaptive_panel_scalar_integral (par->patch, par->hz, 3, */
  /* 							    panel->i-0.5, panel->i+0.5, */
  /* 							    panel->j-0.5, panel->j+0.5); */
}

void b4_patch_compute_volume_times_metric (Patch * patch, HeightCurve hz)
{
  /* VTMParams par; */
  /* par.patch = patch; */
  /* par.hz = hz; */

  /* patch_forall_panels (patch, panel_compute_volume_times_metric, &par); */
}

static gdouble b4_int (gdouble x, gdouble h)
{
  if ( x <= -5.*h/2.)
    return 0.;
  if ( x > 5.*h/2.)
    return 0.;

  if ( -2.5*h < x && x <= -1.5*h)
    return 1./120.*h;

  if ( -3.*h/2 < x && x <= -h/2.)
    return 13.*h/60.;

  if ( -h/2 < x && x <= h/2.)
    return 11./20.*h;

  if ( h/2 < x && x <= 3.*h/2.)
    return 13./60.*h;

  if ( 1.5*h < x && x <= 2.5*h)
    return 1./120.*h;
  return 0.;
}


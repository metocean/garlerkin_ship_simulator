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
#include "boundaries.h"

void patch_destroy (Patch * p)
{
  /* g_ptr_array_set_free_func (p->rows, ); */
  g_ptr_array_free (p->rows, TRUE);
  g_free (p);
}

void patch_add_row (Patch * p)
{
  GArray * new = g_array_new (FALSE, FALSE, sizeof(Panel));

  g_ptr_array_add (p->rows, new);
}

Patch * patch_new ()
{
  Patch * new = g_malloc (sizeof(Patch));
  new->rows = g_ptr_array_new ();
  new->t = NULL;
  return new;
}

void patch_add_panel (Patch * p, Panel panel)
{
  GArray * row = g_ptr_array_index (p->rows, p->rows->len-1);
  panel.i = row->len;
  panel.j = p->rows->len-1;
  g_array_append_val (row, panel);
}

void patch_tag_borders (Patch * p)
{
  gint i, j;
  GArray * row;

  i = 0;
  row = g_ptr_array_index (p->rows, i);
  for ( j = 0; j < row->len; j++) {
    g_array_index (row, Panel, j).border = TRUE;
  }
    
  i = p->rows->len - 1;
  row = g_ptr_array_index (p->rows, i);
  for ( j = 0; j < row->len; j++) {
    g_array_index (row, Panel, j).border = TRUE;
  }

  for ( i = 0; i < p->rows->len ; i++) {
    row = g_ptr_array_index (p->rows, i);
    g_array_index (row, Panel, 0).border = TRUE;
    g_array_index (row, Panel, row->len-1).border = TRUE;
  }
}

void patch_print (Patch * patch, FILE * fp)
{
  gint i, j, k;

  for ( i = 0; i < patch->rows->len ; i++) {
    GArray * row = g_ptr_array_index (patch->rows, i);
    for ( j = 0; j < row->len; j++) {
      Panel p = g_array_index (row, Panel, j);
      for ( k = 0; k < 4; k++)
	point_print (p.p[k], fp, patch->t);
      point_print (p.p[0], fp, patch->t);
      fprintf(fp,"\n\n");
    }
  }
}

void patch_center_print (Patch * p, FILE * fp)
{
  gint i, j;

  for ( i = 0; i < p->rows->len ; i++) {
    GArray * row = g_ptr_array_index (p->rows, i);
    for ( j = 0; j < row->len; j++) {
      Panel p = g_array_index (row, Panel, j);
      fprintf(fp, "%g %g %g %i %i\n", p.var[0],  p.var[1], p.var[2], p.i, p.j);
    }
  }
}

/* Returns the panel that contains the point i, j */
Panel * patch_panel_get (Patch * patch, gint i, gint j)
{
  GArray * row = g_ptr_array_index (patch->rows, 0);

  if ( i < 0 || j < 0 || j > patch->rows->len-1 || i > row->len-1 )
    return NULL;
  
  row = g_ptr_array_index (patch->rows, j);
  return &g_array_index (row, Panel, i);
}

/* Returns TRUE if the point i, j belong to a border cell of patch */
static gboolean patch_border (Patch * patch, gdouble i, gdouble j)
{
  GArray * row = g_ptr_array_index (patch->rows, 0);

  if ( i < 0.5 || j < 0.5 || j > patch->rows->len-1.5 || i > row->len-1.5 )
    return TRUE;
  return FALSE;
}

/* Returns a second order extrapolation of the corner point located at */
/* (dx, dy) form the center of the stencil                             */
static gdouble corner_extrapolation (gdouble v[3][3], gdouble dx, gdouble dy)
{
  return v[1][1] + dx*(v[0][1]-v[2][1])/2. + dy*(v[1][0]-v[1][2])/2.
	+ dx*dx/2.*(v[0][1]+v[2][1]-2.*v[1][1]) + dy*dy/2.*(v[1][0]+v[1][2]-2.*v[1][1])
	+ dx*dy/4.*(v[0][0]+v[2][2]-v[0][2]-v[2][0]);
}

static double patch_eval_y (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_y (&panel, i, j, var);
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

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_x (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Evaluates the variable var at a given location (i,j) in the Patch p
 **/
gdouble patch_eval (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;
  gint N = row->len;
  gint M = p->rows->len;

  /* if ( i < -0.5 || j < -0.5 || i > row->len-0.5 || j > p->rows->len-0.5){ */
  /*   fprintf(stderr,"*** Warning: Trying to evaluate a point that does not belong to this patch ***\n"); */
  /* } */

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
    Panel * pr = patch_panel_get (p, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-1));
    
    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [3][3];
      gdouble dx = (i - (N-2));
      gdouble dy = (j - (M-2));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = PANEL_VAL(p0, var);
	  else
	    v[l][k] = patch_eval (p, rint(i-l), rint(j-k), var);
	}

      return corner_extrapolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v [3][3];
      gdouble dx = (i - (N-2));
      gdouble dy = -(j - 1);

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = PANEL_VAL(p0, var);
	  else
	    v[l][k] = patch_eval (p, rint(i-l),  rint(j+k), var);
	}
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v [3][3];
      gdouble dx = -(i - 1);
      gdouble dy = -(j - 1);

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = PANEL_VAL(p0, var);
	  else
	    v[l][k] = patch_eval (p, rint(i+l),  rint(j+k), var);
	}
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v [3][3];
      gdouble dx = -(i - 1);
      gdouble dy = (j - (M-2));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = PANEL_VAL(p0, var);
	  else
	    v[l][k] = patch_eval (p, rint(i+l),  rint(j-k), var);
	}
      
      return corner_extrapolation (v, dx, dy);
    }
    else {
      gdouble dx;
      gdouble v[3];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = (i - rint(i-1));
	v[0] = patch_eval_y (p, rint(i), j, var);
	v[1] = patch_eval (p, rint(i-1), j, var);
	v[2] = patch_eval (p, rint(i-2), j, var);
      }
      
      if (pl == NULL) {
	dx = (rint(i+1)-i);
	v[0] = patch_eval_y (p, rint(i), j, var);
	v[1] = patch_eval (p, rint(i+1), j, var);
	v[2] = patch_eval (p, rint(i+2), j, var);
      }
      
      if (pt == NULL) {
	dx = (j - rint(j-1));
	v[0] = patch_eval_x (p, i, rint(j), var);
	v[1] = patch_eval (p, i, rint(j-1), var);
	v[2] = patch_eval (p, i, rint(j-2), var);
      }

      if (pb == NULL) {
	dx = (rint(j+1)-j);
	v[0] = patch_eval_x (p, i, rint(j), var);
	v[1] = patch_eval (p, i, rint(j+1), var);
	v[2] = patch_eval (p, i, rint(j+2), var);
      }

      return v[1] + dx*(v[0] - v[2])/2.	+ dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval (&panel, i, j, var);
    }
  }

  return sum;
}

static double patch_eval_dx_x (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dx_x (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Returns the first derivative along the x direction in i,j
 **/
gdouble patch_eval_dx (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;

   if ( i < -0.5 || j < -0.5 || i > row->len-0.5 || j > p->rows->len-0.5){
    fprintf(stderr,"*** Warning: Trying to evaluate a point that does not belong to this patch ***\n");
  }

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-11;
  if ( j == -0.5 || j == 0.5)
    j += 1e-11;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-11;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-11;

   /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-1));
    
    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dx (p, rint(i-l), rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dx (p, rint(i-l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;

      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dx (p, rint(i+l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dx (p, rint(i+l),  rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else {
      gdouble dx, dxx = 1.;
      gdouble v[3];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = (i - rint(i-1));
	dxx = 1.;
	v[0] = patch_eval_y (p, rint(i), j, var);
	v[1] = patch_eval (p, rint(i-1), j, var);
	v[2] = patch_eval (p, rint(i-2), j, var);
      }
      
      if (pl == NULL) {
	dx = (rint(i+1)-i);
	dxx = -1.;
	v[0] = patch_eval_y (p, rint(i), j, var);
	v[1] = patch_eval (p, rint(i+1), j, var);
	v[2] = patch_eval (p, rint(i+2), j, var);
      }
      
      if (pt == NULL) {
	dx = (j - rint(j-1));
	dxx = 0.;
	v[0] = patch_eval_dx_x (p, i, rint(j), var);
	v[1] = patch_eval_dx (p, i, rint(j-1), var);
	v[2] = patch_eval_dx (p, i, rint(j-2), var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }

      if (pb == NULL) {
	dx = (rint(j+1)-j);
	dxx = 0.;
	v[0] = patch_eval_dx_x (p, i, rint(j), var);
	v[1] = patch_eval_dx (p, i, rint(j+1), var);
	v[2] = patch_eval_dx (p, i, rint(j+2), var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }

      return dxx*(v[0]-v[2])/2. + dx*dxx*(v[0] + v[2]- 2.*v[1]);
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dx (&panel, i, j, var);
    }
  }

  return sum;
}

static double patch_eval_dy_y (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dy_y (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Returns the first derivative along the y direction in i,j
 **/
gdouble patch_eval_dy (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;

   if ( i < -0.5 || j < -0.5 || i > row->len-0.5 || j > p->rows->len-0.5){
    fprintf(stderr,"*** Warning: Trying to evaluate a point that does not belong to this patch ***\n");
  }

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-11;
  if ( j == -0.5 || j == 0.5)
    j += 1e-11;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-11;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-11;

   /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-1));
    
    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dy (p, rint(i-l), rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dy (p, rint(i-l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;

      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dy (p, rint(i+l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dy (p, rint(i+l),  rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else {
      gdouble dx, dxx = 1.;
      gdouble v[3];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = (i - rint(i-1));
	v[0] = patch_eval_dy_y (p, rint(i), j, var);
	v[1] = patch_eval_dy (p, rint(i-1), j, var);
	v[2] = patch_eval_dy (p, rint(i-2), j, var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }
      
      if (pl == NULL) {
	dx = (rint(i+1)-i);
	v[0] = patch_eval_dy_y (p, rint(i), j, var);
	v[1] = patch_eval_dy (p, rint(i+1), j, var);
	v[2] = patch_eval_dy (p, rint(i+2), j, var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }
      
      if (pt == NULL) {
	dx = (j - rint(j-1));
	dxx = 1.;
	v[0] = patch_eval_x (p, i, rint(j), var);
	v[1] = patch_eval (p, i, rint(j-1), var);
	v[2] = patch_eval (p, i, rint(j-2), var);
      }

      if (pb == NULL) {
	dx = (rint(j+1)-j);
	dxx = -1;
	v[0] = patch_eval_x (p, i, rint(j), var);
	v[1] = patch_eval (p, i, rint(j+1), var);
	v[2] = patch_eval (p, i, rint(j+2), var);
      }

      return dxx*(v[0]-v[2])/2. + dx*dxx*(v[0] + v[2]- 2.*v[1]);
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dy (&panel, i, j, var);
    }
  }

  return sum;
}

static double patch_eval_dxx_x (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dxx_x (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Returns the second derivative along the x direction in i,j
 **/
gdouble patch_eval_dxx (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;

   if ( i < -0.5 || j < -0.5 || i > row->len-0.5 || j > p->rows->len-0.5){
    fprintf(stderr,"*** Warning: Trying to evaluate a point that does not belong to this patch ***\n");
  }

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-11;
  if ( j == -0.5 || j == 0.5)
    j += 1e-11;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-11;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-11;

   /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-1));
    
    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxx (p, rint(i-l), rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxx (p, rint(i-l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;

      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxx (p, rint(i+l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxx (p, rint(i+l),  rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else {
      gdouble dx, dxx = 1.;
      gdouble v[3];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = (i - rint(i-1));
	dxx = 1.;
	v[0] = patch_eval_y (p, rint(i), j, var);
	v[1] = patch_eval (p, rint(i-1), j, var);
	v[2] = patch_eval (p, rint(i-2), j, var);
      }
      
      if (pl == NULL) {
	dx = (rint(i+1)-i);
	dxx = -1.;
	v[0] = patch_eval_y (p, rint(i), j, var);
	v[1] = patch_eval (p, rint(i+1), j, var);
	v[2] = patch_eval (p, rint(i+2), j, var);
      }
      
      if (pt == NULL) {
	dx = (j - rint(j-1));
	dxx = 0.;
	v[0] = patch_eval_dxx_x (p, i, rint(j), var);
	v[1] = patch_eval_dxx (p, i, rint(j-1), var);
	v[2] = patch_eval_dxx (p, i, rint(j-2), var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }

      if (pb == NULL) {
	dx = (rint(j+1)-j);
	dxx = 0.;
	v[0] = patch_eval_dxx_x (p, i, rint(j), var);
	v[1] = patch_eval_dxx (p, i, rint(j+1), var);
	v[2] = patch_eval_dxx (p, i, rint(j+2), var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }

      return dxx*dxx*(v[0] + v[2]- 2.*v[1]);
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dxx (&panel, i, j, var);
    }
  }

  return sum;
}

static double patch_eval_dyy_y (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dyy_y (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Returns the first derivative along the y direction in i,j
 **/
gdouble patch_eval_dyy (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;

   if ( i < -0.5 || j < -0.5 || i > row->len-0.5 || j > p->rows->len-0.5){
    fprintf(stderr,"*** Warning: Trying to evaluate a point that does not belong to this patch ***\n");
  }

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-11;
  if ( j == -0.5 || j == 0.5)
    j += 1e-11;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-11;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-11;

   /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-1));
    
    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dyy (p, rint(i-l), rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dyy (p, rint(i-l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;

      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dyy (p, rint(i+l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dyy (p, rint(i+l),  rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else {
      gdouble dx, dxx = 1.;
      gdouble v[3];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = (i - rint(i-1));
	v[0] = patch_eval_dyy_y (p, rint(i), j, var);
	v[1] = patch_eval_dyy (p, rint(i-1), j, var);
	v[2] = patch_eval_dyy (p, rint(i-2), j, var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }
      
      if (pl == NULL) {
	dx = (rint(i+1)-i);
	v[0] = patch_eval_dyy_y (p, rint(i), j, var);
	v[1] = patch_eval_dyy (p, rint(i+1), j, var);
	v[2] = patch_eval_dyy (p, rint(i+2), j, var);
	return  v[1] + dx*(v[0] - v[2])/2. + dx*dx/2.*((v[0] + v[2]- 2.*v[1]));
      }
      
      if (pt == NULL) {
	dx = (j - rint(j-1));
	dxx = 1.;
	v[0] = patch_eval_x (p, i, rint(j), var);
	v[1] = patch_eval (p, i, rint(j-1), var);
	v[2] = patch_eval (p, i, rint(j-2), var);
      }

      if (pb == NULL) {
	dx = (rint(j+1)-j);
	dxx = -1;
	v[0] = patch_eval_x (p, i, rint(j), var);
	v[1] = patch_eval (p, i, rint(j+1), var);
	v[2] = patch_eval (p, i, rint(j+2), var);
      }

      return dxx*dxx*(v[0] + v[2]- 2.*v[1]);
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dyy (&panel, i, j, var);
    }
  }

  return sum;
}

/**
 * Returns the second derivative along the x direction in i,j
 **/
gdouble patch_eval_dxy (Patch * p, gdouble i, gdouble j, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (p->rows, 0);
  Panel panel;

   if ( i < -0.5 || j < -0.5 || i > row->len-0.5 || j > p->rows->len-0.5){
    fprintf(stderr,"*** Warning: Trying to evaluate a point that does not belong to this patch ***\n");
  }

  /* Avoids roun-off problems with rint near borders*/
  if ( i == -0.5 || i == 0.5)
    i += 1e-11;
  if ( j == -0.5 || j == 0.5)
    j += 1e-11;
  if ( i == row->len-0.5 || i == row->len-1.5)
    i -= 1e-11;
  if ( j == p->rows->len-0.5 || j == p->rows->len-1.5)
    j -= 1e-11;

   /* Border/corners special treatment */
  if (patch_border (p, i, j)) {
    Panel * p0 = patch_panel_get (p, rint(i), rint(j));
    g_assert (p0 != NULL);
    Panel * pr = patch_panel_get (p, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (p, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (p, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (p, rint(i), rint(j-1));
    
    /* FIRST order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) {
      gdouble v [3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxy (p, rint(i-l), rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pr == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = (i - rint(i-1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxy (p, rint(i-l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;

      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pb == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = -(j - rint(j+1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxy (p, rint(i+l),  rint(j+k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else if (pl == NULL && pt == NULL)  {
      gdouble v[3][3];
      gdouble dx = -(i - rint(i+1));
      gdouble dy = (j - rint(j-1));

      for (l = 0; l < 3 ; l++)
	for (k = 0; k < 3; k++) {
	  if ( k == 0 && l == 0)
	    v[l][k] = 0.;
	  else
	    v[l][k] = patch_eval_dxy (p, rint(i+l),  rint(j-k), var);
	}
      v[0][0] = (v[1][0] + v[0][1] + v[1][1])/3.;
      
      return corner_extrapolation (v, dx, dy);
    }
    else {
      gdouble dx;
      gdouble v[3];
      
      /* Non corner cells */
      if (pr == NULL) {
	dx = (i - rint(i-1));
	v[1] = patch_eval_dxy (p, rint(i-1), j, var);
	v[2] = patch_eval_dxy (p, rint(i-2), j, var);
      }
      
      if (pl == NULL) {
	dx = (rint(i+1)-i);
	v[1] = patch_eval_dxy (p, rint(i+1), j, var);
	v[2] = patch_eval_dxy (p, rint(i+2), j, var);
      }
      
      if (pt == NULL) {
	dx = (j - rint(j-1));
	v[1] = patch_eval_dxy (p, i, rint(j-1), var);
	v[2] = patch_eval_dxy (p, i, rint(j-2), var);
      }

      if (pb == NULL) {
	dx = (rint(j+1)-j);
	v[1] = patch_eval_dxy (p, i, rint(j+1), var);
	v[2] = patch_eval_dxy (p, i, rint(j+2), var);
      }

      return  v[1] + dx*(v[1]-v[2]);
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_dxy (&panel, i, j, var);
    }
  }

  return sum;
}


static gdouble patch_eval_int_y (Patch * p, gint i, gint j, gint var)
{
  gint k;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (p->rows, k);
    /* for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) { */
      panel = g_array_index (row, Panel, i);
      sum += panel_eval_y (&panel, i, j, var);
    /* } */
  }

  return sum;
}

static gdouble patch_eval_int_x (Patch * p, gint i, gint j, gint var)
{
  gint l;
  gdouble sum = 0.;
  GArray * row;
  Panel panel;

  /* for ( k = MAX (0, j-2); k <= MIN (p->rows->len-1, j+2) ;k++) { */
  row = g_ptr_array_index (p->rows, j);
  for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
    panel = g_array_index (row, Panel, l);
    sum += panel_eval_int_x (&panel, i, j, var);
  }
  /* } */

  return sum;
}

/**
 * Evaluates the integral of var over the panel p in a spline sense
 **/
gdouble patch_panel_int (Patch * patch, Panel * p0, gint var)
{
  gint k, l;
  gdouble sum = 0.;
  GArray * row = g_ptr_array_index (patch->rows, 0);
  Panel panel;

  gint i = p0->i;
  gint j = p0->j;

  /* Border/corners special treatment */
  if (patch_border (patch, i, j)) {
    Panel * pr = patch_panel_get (patch, rint(i+1), rint(j));
    Panel * pl = patch_panel_get (patch, rint(i-1), rint(j));
    Panel * pt = patch_panel_get (patch, rint(i), rint(j+1));
    Panel * pb = patch_panel_get (patch, rint(i), rint(j-1));
    
    /* Second order border and corner cells extrapolation */
    if (pr == NULL && pt == NULL) 
      return PANEL_VAL (p0, var);
    else if (pr == NULL && pb == NULL)
      return PANEL_VAL (p0, var);
    else if (pl == NULL && pb == NULL)
      return PANEL_VAL (p0, var);
    else if (pl == NULL && pt == NULL)
      return PANEL_VAL (p0, var);
    else {
      
      /* Non corner cells */
      if (pr == NULL)
	return patch_eval_y (patch, i, j, var);
      
      if (pl == NULL)
	return patch_eval_y (patch, i, j, var);
       
      if (pt == NULL)
	return patch_eval_x (patch, i, j, var);

      if (pb == NULL)
	return patch_eval_x (patch, i, j, var);

      g_assert_not_reached ();
    }
  }
  
  for ( k = MAX (0, j-2); k <= MIN (patch->rows->len-1, j+2) ;k++) {
    row = g_ptr_array_index (patch->rows, k);
    for ( l = MAX (0, i-2); l <= MIN (row->len-1, i+2); l++) {
      panel = g_array_index (row, Panel, l);
      sum += panel_eval_int (&panel, i, j, var);
    }
  }

  return sum;
}

Point patch_eval_point (Patch * patch, gdouble i, gdouble j)
{
  Point p;

  p.x = patch_eval (patch, i, j, 0);
  p.y = patch_eval (patch, i, j, 1);
  p.z = patch_eval (patch, i, j, 2);

  return p;
}

Point patch_eval_transformed_point (Patch * patch, gdouble i, gdouble j)
{
  Point p;
  Transformation * t = patch->t;

  g_assert (patch->t != NULL);

  p.x = patch_eval (patch, i, j, 0);
  p.y = patch_eval (patch, i, j, 1);
  p.z = patch_eval (patch, i, j, 2);

  return transform_point (p, &t->xg, &t->t, &t->euler_m);
}

void patch_forall_panels (Patch * patch, PanelFunction func, gpointer data)
{
  gint i, j;

  for ( j = 0; j < patch->rows->len; j++) {
    GArray * row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      func (&panel, data);
      g_array_index (row, Panel, i) = panel;
    }
  }
}

static gdouble free_surface (gdouble x, gdouble y)
{
  return cos(x/5.)/* 0.5 */;
}

static GSList * intersection (Patch * patch, gdouble il, gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 1e-2;
  GSList * l = NULL;

  /* Check for intersection */
  gdouble dz[4];
  Point c[4];
  gint i;
  c[0] = patch_eval_transformed_point (patch, il, jl);
  c[1] = patch_eval_transformed_point (patch, ih, jl);
  c[2] = patch_eval_transformed_point (patch, ih, jh);
  c[3] = patch_eval_transformed_point (patch, il, jh);

  for ( i = 0; i < 4; i++)
    dz[i] = c[i].z - free_surface (c[i].x, c[i].y);

  if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.) ||
       (dz[0] < 0. && dz[1] < 0. && dz[2] < 0. && dz[3] < 0.) )
    return l;

  /* We refine */
  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  if (fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance) {
    Point * p = g_malloc (sizeof(Point));
    *p = patch_eval_transformed_point (patch, ihalf, jhalf);
    l = g_slist_append (l, p);
    return l;
  }

  l = g_slist_concat (l, intersection (patch, il, ihalf, jl, jhalf));
  l = g_slist_concat (l, intersection (patch, ihalf, ih, jl, jhalf));
  l = g_slist_concat (l, intersection (patch, il, ihalf, jhalf, jh));
  l = g_slist_concat (l, intersection (patch, ihalf, ih, jhalf, jh));

  return l; // ??? Used to work without it
}

static gint list_angle_sort (Point * p1, Point * p2)
{
  gdouble angle1 = p1->y/fabs(p1->y)*acos(p1->x/sqrt(p1->x*p1->x+p1->y*p1->y));
  gdouble angle2 = p2->y/fabs(p2->y)*acos(p2->x/sqrt(p2->x*p2->x+p2->y*p2->y));

  return (angle1 > angle2 ? -1 : 1);
}

GSList * sort_list_by_angle (GSList * l)
{
  return g_slist_sort (l, (GCompareFunc) list_angle_sort);
}

GSList * patch_intersect_with_free_surface (Patch * patch)
{
  GSList * l = NULL;
  gint i, j;

  for ( j = 0; j < patch->rows->len; j++) {
    GArray * row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      l = g_slist_concat (l, intersection (patch, panel.i - 0.5, panel.i + 0.5, panel.j - 0.5, panel.j + 0.5));
    }
  }

  /* l = sort_list_by_angle (l); */

  return l;
}

/**
 * Return a vector containing the normal at the surface at point (i,j)
 * The normal is expressed in the (x,y,z) coordinate system.
 * The normal is oriented towards the inside of the hull.
 * The returned normal is non-normalised.
 **/
Vector patch_normal (Patch * p, gdouble i, gdouble j)
{
  gdouble xi = patch_eval_dx (p,i,j,0);
  gdouble xj = patch_eval_dy (p,i,j,0);
  gdouble yi = patch_eval_dx (p,i,j,1);
  gdouble yj = patch_eval_dy (p,i,j,1);
  gdouble zi = patch_eval_dx (p, i, j, 2);
  gdouble zj = patch_eval_dy (p, i, j, 2);

  gdouble  J = 0.;
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

Vector patch_unit_normal (Patch * p, gdouble i, gdouble j)
{
  return vector_normalise (patch_normal (p, i, j));
}

gdouble patch_local_metric (Patch * p, gdouble i, gdouble j)
{
  gdouble xi = patch_eval_dx (p,i,j,0);
  gdouble xj = patch_eval_dy (p,i,j,0);
  gdouble yi = patch_eval_dx (p,i,j,1);
  gdouble yj = patch_eval_dy (p,i,j,1);
  gdouble zi = patch_eval_dx (p,i,j,2);
  gdouble zj = patch_eval_dy (p,i,j,2);
  Vector g;

  g.x = yi*zj-yj*zi;
  g.y = zi*xj-zj*xi;
  g.z = xi*yj-xj*yi;

  return vector_norm (g);
}

void patch_normal_print (Patch * p, FILE * fp)
{
  gdouble i, j;
  GArray * row = g_ptr_array_index (p->rows, 0);
  /* fprintf(stderr,"SIZE: %i %i\n",p->rows->len, row->len); */

  /* /\* for ( j = /\\* -0.5 *\\/1.5; j <= 1.6/\\* p->rows->len-0.5 *\\/ ; j += 0.5) { *\/ */
  /* for ( j = -0.5 ; j <= p->rows->len-0.5 ; j += 0.5) { */
  /*   /\* fprintf(fp,"#Row: %f\n", j); *\/ */
  /*   for ( i = -0.5; i <= row->len-0.5; i += 0.5) { */
  /*   /\* for ( i = -0.2; i <= -0.1; i += 0.5) { *\/ */
  /*     /\* fprintf(fp,"#Line: %f\n", i); *\/ */
  /*     /\* fprintf(fp, "%g %g %g\n", patch_eval (p, i,j,0), patch_eval (p, i,j,1), patch_eval (p, i,j,2)); *\/ */
  /*     /\* fprintf(fp, "%g %g %g\n", j, patch_eval (p, i,j,1), patch_eval_dx (p, i,j,1)); *\/ */
  /*     fprintf(fp, "%g %g %g %g\n", i/\* patch_eval (p, i,j,0) *\/, j/\* patch_eval (p, i,j,1) *\/, patch_eval (p, i,j,2), patch_eval_dxy (p, i,j,1)); */
  /*   } */
  /* } */

  for ( j = -0.5 ; j <= p->rows->len-0.5 ; j += 0.5) {
    for ( i = -0.5; i <= row->len-0.5; i += 0.5) {
      Vector n = patch_unit_normal (p, i, j);
      /* n = vector_normalise (panel_first_order_normal (patch_panel_get (p, i, j))); */
      Point p0 = patch_eval_point (p, i, j);

      fprintf(fp, "%g %g %g \n %g %g %g \n\n\n", p0.x, p0.y, p0.z, p0.x+n.x, p0.y+n.y, p0.z+n.z);
    }
  }
}

/* Unidirectional cubic spline */
static gdouble b2 (gdouble x, gdouble h)
{
  if (x <= -3.*h/2.)
    return 0.;
  if (x > 3.*h/2.)
    return 0.;

  if (-3.*h/2. < x && x <= -h/2.)
    return 1./(2.*h*h) * pow(x+3.*h/2., 2.);

  if (-h/2. < x && h/2. >= x)
    return 1./(h*h) * (-x*x + 3.*h*h/4.);

  if (3.*h/2. >= x && x > h/2.)
    return 1./(2.*h*h) * pow(-x+3.*h/2., 2.);
  g_assert_not_reached ();
}

static gdouble generic_2D_bspline (gdouble i, gdouble j, gdouble x, gdouble y)
{
  return b2 (x-i, 1.)*b2 (y-j, 1.);
}

/**
 * Finds the best patch fit using the corner values of the panels.
 * The center values of the panels are determined that way.
 **/
void patch_fit_panels (Patch * patch)
{
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  gint N = row->len + 1;
  gint M = patch->rows->len + 1;
  Point p[N][M]; /* Points to fit */

  Point p0;
  p0.x = p0.y = p0.z = 0.;

  for ( i = 0; i < N; i++)
    for ( j = 0; j < M; j++)
      p[i][j] = p0;
  
  /* Copy data to p */
  for ( j = 0; j <  patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len ; i++) {
      Panel panel = g_array_index (row, Panel, i);
      p[i][j] = panel.p[2];
    }
  }

  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    i = row->len-1;
    Panel panel = g_array_index (row, Panel, i);
    p[row->len][j] = panel.p[1];
  }

  j = patch->rows->len-1;
  row = g_ptr_array_index (patch->rows, j);
  for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      p[i][patch->rows->len] = panel.p[3];
  }

  i = row->len-1;
  Panel panel = g_array_index (row, Panel, i);
  p[row->len][patch->rows->len] = panel.p[0];
  
  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;
  gint n = (N-2)*(M-2);

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
  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++) {
      gsl_vector_set (cx, (i-1)+(j-1)*(N-2), (p[i][j].x+p[i+1][j].x+p[i][j+1].x+p[i+1][j+1].x)/4.);
      gsl_vector_set (cy, (i-1)+(j-1)*(N-2), (p[i][j].y+p[i+1][j].y+p[i][j+1].y+p[i+1][j+1].y)/4.);
      gsl_vector_set (cz, (i-1)+(j-1)*(N-2), (p[i][j].z+p[i+1][j].z+p[i][j+1].z+p[i+1][j+1].z)/4.);
    }

  gint k, l;
  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++) {
      for ( k = 1; k < N-1; k++)
  	for ( l = 1; l < M-1; l++) {
  	  gsl_matrix_set (X, (k-1)+(l-1)*(N-2), (i-1)+(j-1)*(N-2), generic_2D_bspline ((gdouble) i, (gdouble) j, (gdouble) k, (gdouble) l));
  	}
    }
  
  for ( k = 1; k < N-1; k++)
    for ( l = 1; l < M-1; l++) {
      gdouble xi = p[k][l].x;
      gdouble yi = p[k][l].y;
      gdouble zi = p[k][l].z;

      for ( i = 0; i < N; i++) {
      	xi -= p[i][0].x*generic_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	xi -= p[i][M-1].x*generic_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
  	yi -= p[i][0].y*generic_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	yi -= p[i][M-1].y*generic_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
  	zi -= p[i][0].z*generic_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	zi -= p[i][M-1].z*generic_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
      }

      for ( j = 0; j < M; j++) {
      	xi -= p[0][j].x*generic_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	xi -= p[N-1][j].x*generic_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
  	yi -= p[0][j].y*generic_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	yi -= p[N-1][j].y*generic_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
  	zi -= p[0][j].z*generic_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	zi -= p[N-1][j].z*generic_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
      }

      xi += p[0][0].x*generic_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      xi += p[0][M-1].x*generic_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      xi += p[N-1][M-1].x*generic_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      xi += p[N-1][0].x*generic_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      yi += p[0][0].y*generic_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      yi += p[0][M-1].y*generic_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      yi += p[N-1][M-1].y*generic_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      yi += p[N-1][0].y*generic_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      zi += p[0][0].z*generic_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      zi += p[0][M-1].z*generic_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      zi += p[N-1][M-1].z*generic_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      zi += p[N-1][0].z*generic_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);

      gsl_vector_set (x, (k-1)+(l-1)*(N-2), xi);
      gsl_vector_set (y, (k-1)+(l-1)*(N-2), yi);
      gsl_vector_set (z, (k-1)+(l-1)*(N-2), zi);
    }


  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, n);

  /* Multidimensional fitting usind singular values decomposition method */
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, 1e-6, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, 1e-6, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, 1e-6, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++) {
      p[i][j].x = gsl_vector_get(cx, i-1+(j-1)*(N-2));
      p[i][j].y = gsl_vector_get(cy, i-1+(j-1)*(N-2));
      p[i][j].z = gsl_vector_get(cz, i-1+(j-1)*(N-2));
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

  FILE * fp = fopen("fit.tmp","w");
  gdouble xx, yy;
  for ( yy = 0.5; yy <= M-1.5; yy+=0.125) {
    for ( xx = 0.5; xx <= N-1.5; xx+=0.125) {
      
      gdouble px = 0., py = 0., pz = 0.;
      for ( k = 0; k < N; k++)
      	for ( l = 0; l < M; l++) {
      	  px += p[k][l].x*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) xx, (gdouble) yy);
      	  py += p[k][l].y*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) xx, (gdouble) yy);
      	  pz += p[k][l].z*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) xx, (gdouble) yy);
      	}
     
      fprintf(fp, "%g %g %g \n", px, py, pz);
    }
  }
  
  
  /* fclose (fp); */
  
  fp = fopen("fit2.tmp","w");
  /* Interpolate the spline to the center points of the panels */
  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      panel.var[0] = 0.;
      panel.var[1] = 0.;
      panel.var[2] = 0.;
      for ( k = 0; k < N; k++)
      	for ( l = 0; l < M; l++) {
      	  panel.var[0] += p[k][l].x*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
      	  panel.var[1] += p[k][l].y*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
      	  panel.var[2] += p[k][l].z*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
      	}
      panel.p[0] = p[i+1][j+1];
      panel.p[1] = p[i+1][j];
      panel.p[2] = p[i][j];
      panel.p[3] = p[i][j+1];
      g_array_index (row, Panel, i) = panel;
      fprintf(fp, "%g %g %g \n", panel.var[0], panel.var[1], panel.var[2]);
    }
  }
  fclose (fp);
}

void patch_prefit_panels (Patch * patch)
{
  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint i, j;
  gint N = row->len + 1;
  gint M = patch->rows->len + 1;
  Point p[N][M]; /* Points to fit */

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
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;
  gint n = (N-2)*(M-2);

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

  gint k, l;
  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++) {
      for ( k = 1; k < N-1; k++)
	for ( l = 1; l < M-1; l++) {
	  gsl_matrix_set (X, (k-1)+(l-1)*(N-2), (i-1)+(j-1)*(N-2), generic_2D_bspline ((gdouble) i, (gdouble) j, (gdouble) k, (gdouble) l));
	}
    }
  
  for ( k = 1; k < N-1; k++)
    for ( l = 1; l < M-1; l++) {
      gdouble xi = p[k][l].x;
      gdouble yi = p[k][l].y;
      gdouble zi = p[k][l].z;

      for ( i = 0; i < N; i++) {
      	xi -= p[i][0].x*generic_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	xi -= p[i][M-1].x*generic_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
	yi -= p[i][0].y*generic_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	yi -= p[i][M-1].y*generic_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
	zi -= p[i][0].z*generic_2D_bspline ((gdouble) i, (gdouble) 0, (gdouble) k, (gdouble) l);
      	zi -= p[i][M-1].z*generic_2D_bspline ((gdouble) i, (gdouble) M-1, (gdouble) k, (gdouble) l);
      }

      for ( j = 0; j < M; j++) {
      	xi -= p[0][j].x*generic_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	xi -= p[N-1][j].x*generic_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
	yi -= p[0][j].y*generic_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	yi -= p[N-1][j].y*generic_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
	zi -= p[0][j].z*generic_2D_bspline ((gdouble) 0, (gdouble) j, (gdouble) k, (gdouble) l);
      	zi -= p[N-1][j].z*generic_2D_bspline ((gdouble) N-1, (gdouble) j, (gdouble) k, (gdouble) l);
      }

      xi += p[0][0].x*generic_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      xi += p[0][M-1].x*generic_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      xi += p[N-1][M-1].x*generic_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      xi += p[N-1][0].x*generic_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      yi += p[0][0].y*generic_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      yi += p[0][M-1].y*generic_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      yi += p[N-1][M-1].y*generic_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      yi += p[N-1][0].y*generic_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);
      zi += p[0][0].z*generic_2D_bspline ((gdouble) 0, (gdouble) 0, (gdouble) k, (gdouble) l);
      zi += p[0][M-1].z*generic_2D_bspline ((gdouble) 0, (gdouble) M-1, (gdouble) k, (gdouble) l);
      zi += p[N-1][M-1].z*generic_2D_bspline ((gdouble) N-1, (gdouble) M-1, (gdouble) k, (gdouble) l);
      zi += p[N-1][0].z*generic_2D_bspline ((gdouble) N-1, (gdouble) 0, (gdouble) k, (gdouble) l);

      gsl_vector_set (x, (k-1)+(l-1)*(N-2), xi);
      gsl_vector_set (y, (k-1)+(l-1)*(N-2), yi);
      gsl_vector_set (z, (k-1)+(l-1)*(N-2), zi);
    }


  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, n);

  /* Multidimensional fitting usind singular values decomposition method */
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, 1e-6, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, 1e-6, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, 1e-6, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 pre-fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++) {
      p[i][j].x = gsl_vector_get(cx, i-1+(j-1)*(N-2));
      p[i][j].y = gsl_vector_get(cy, i-1+(j-1)*(N-2));
      p[i][j].z = gsl_vector_get(cz, i-1+(j-1)*(N-2));
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
      panel.var[0] = 0.;
      panel.var[1] = 0.;
      panel.var[2] = 0.;
      for ( k = 0; k < N; k++)
      	for ( l = 0; l < M; l++) {
      	  panel.var[0] += p[k][l].x*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
      	  panel.var[1] += p[k][l].y*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
      	  panel.var[2] += p[k][l].z*generic_2D_bspline ((gdouble) k, (gdouble) l, (gdouble) i+0.5, (gdouble) j+0.5);
      	}
      /* if ( j == 0 || i == 0 || i == patch->rows->len-1 || j == patch->rows-1) { */
      /* 	panel.p[0] = p[i+1][j+1]; */
      /* 	panel.p[1] = p[i+1][j]; */
      /* 	panel.p[2] = p[i][j]; */
      /* 	panel.p[3] = p[i][j+1]; */
      /* } */
      g_array_index (row, Panel, i) = panel;
    }
  }
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

  if (hz != NULL) {
    /* Check if we are under the free-surface */
    gdouble dz[4];
    dz[0] = patch_eval (patch, il, jl, 2)
      - hz (patch_eval (patch, il, jl, 0), patch_eval (patch, il, jl, 1), 0, NULL);
    dz[1] = patch_eval (patch, ih, jl, 2)
      - hz (patch_eval (patch, ih, jl, 0), patch_eval (patch, ih, jl, 1), 0, NULL);
    dz[2] = patch_eval (patch, ih, jh, 2)
      - hz (patch_eval (patch, ih, jh, 0), patch_eval (patch, ih, jh, 1), 0, NULL);
    dz[3] = patch_eval (patch, il, jh, 2)
      - hz (patch_eval (patch, il, jh, 0), patch_eval (patch, il, jh, 1), 0, NULL);
    
    if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.))
      return v;
  }

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Unit normal at the center of the small piece of panel */
    /* Vector n = patch_unit_normal (patch, (il+ih)/2., (jl+jh)/2.); */
    Vector n = patch_normal (patch, (il+ih)/2., (jl+jh)/2.);
    
    /* Surface integral of var over the small piece of panel */
    /* NB: Those are already computed by patch unit normal   */
    /* some optimization is possibly needed here             */
    /* gdouble xi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2.,0); */
    /* gdouble xj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2.,0); */
    /* gdouble yi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2.,1); */
    /* gdouble yj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2.,1); */
    /* gdouble J = fabs(xi*yj-xj*yi); */
    gdouble J = vector_norm (n);

    /* /\* Correction for vertical and horizontal walls *\/ */
    /* if ( J == 0)  { */
    /*   gdouble zi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2., 2); */
    /*   gdouble zj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2., 2); */

    /*   if (xi == 0. && xj == 0.) */
    /* 	J = fabs(zi*yj-zj*yi); */
    /*   else */
    /* 	J = fabs(zi*xj-zj*xi);  */
    /* } */

    n = vector_normalise (v);

    gdouble val = patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)*J*fabs(ih-il)*fabs(jh-jl);

    v.x = val*n.x;
    v.y = val*n.y;
    v.z = val*n.z;
    
    return v;
  }

  /* We refine */
  v = vector_sum (v, adaptive_panel_integral (patch, hz, var, il, ihalf, jl, jhalf));
  v = vector_sum (v, adaptive_panel_integral (patch, hz, var, ihalf, ih, jl, jhalf));
  v = vector_sum (v, adaptive_panel_integral (patch, hz, var, il, ihalf, jhalf, jh));
  v = vector_sum (v, adaptive_panel_integral (patch, hz, var, ihalf, ih, jhalf, jh));

  return v;
}

Vector patch_adaptive_panel_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  gint i, j;
  Vector sum;
  sum.x = sum.y = sum.z = 0.;

  for ( j = 0; j < patch->rows->len ; j++) {
    GArray * row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      sum = vector_sum (sum, adaptive_panel_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5));
    }
  }
  
  /* if (patch->t != NULL) */
  /*   return transform_vector (sum, &patch->t); */
  /* else */
    return sum;
}

/**
 * Returns the value of the integral of a given variable times the local unit normal
 * var over the bspline surface of the patch.
 **/
static gdouble adaptive_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var,
					       gdouble il, gdouble ih, gdouble jl, gdouble jh)
{
  gdouble tolerance = 1./* 0.00625 */;
  /* gdouble border_tolerance = 0.001; */
  gdouble v = 0.;

  if (hz != NULL) {
    /* Check if we are under the free-surface */
    gdouble dz[4];
    dz[0] = patch_eval (patch, il, jl, 2)
      - hz (patch_eval (patch, il, jl, 0), patch_eval (patch, il, jl, 1), 0, NULL);
    dz[1] = patch_eval (patch, ih, jl, 2)
      - hz (patch_eval (patch, ih, jl, 0), patch_eval (patch, ih, jl, 1), 0, NULL);
    dz[2] = patch_eval (patch, ih, jh, 2)
      - hz (patch_eval (patch, ih, jh, 0), patch_eval (patch, ih, jh, 1), 0, NULL);
    dz[3] = patch_eval (patch, il, jh, 2)
      - hz (patch_eval (patch, il, jh, 0), patch_eval (patch, il, jh, 1), 0, NULL);
    
    if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.))
      return v;
  }

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Surface integral of var over the small piece of panel */
    /* NB: Those are already computed by patch unit normal   */
    /* some optimization is possibly needed here             */
    /* gdouble xi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2.,0); */
    /* gdouble xj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2.,0); */
    /* gdouble yi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2.,1); */
    /* gdouble yj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2.,1); */
    /* gdouble zi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2., 2); */
    /* gdouble zj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2., 2); */
    /* gdouble J0 = fabs(xi*yj-xj*yi); */

    /* /\* Correction for vertical and horizontal walls *\/ */
    /* if ( J0 == 0)  { */
    /*   /\* gdouble zi = patch_eval_dx (patch, (il+ih)/2., (jl+jh)/2., 2); *\/ */
    /*   /\* gdouble zj = patch_eval_dy (patch, (il+ih)/2., (jl+jh)/2., 2); *\/ */

    /*   if (xi == 0. && xj == 0.) */
    /* 	J0 = fabs(zi*yj-zj*yi); */
    /*   else */
    /* 	J0 = fabs(zi*xj-zj*xi); */
    /* } */
 
    /* fprintf(stdout, "%g \n", patch_local_metric (patch, (il+ih)/2., (jl+jh)/2.)); */

    return /* patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)* */ 
      patch_local_metric (patch, (il+ih)/2., (jl+jh)/2.)*fabs(ih-il)*fabs(jh-jl);
  }

  /* We refine */
  v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jl, jhalf);
  v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jl, jhalf);
  v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jhalf, jh);
  v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jhalf, jh);

  return v;
}

gdouble patch_adaptive_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  gint i, j;
  gdouble sum = 0.;

  for ( j = 0; j < patch->rows->len ; j++) {
    GArray * row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      sum += adaptive_panel_scalar_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5);
    }
  }
  return sum;
}

static void panel_compute_metric (Panel * panel, gpointer data)
{
  Patch * patch = (Patch *) data;
  gdouble xi = patch_eval_dx (patch, PANEL_VAL (panel, 0), PANEL_VAL (panel, 1), 0);
  gdouble xj = patch_eval_dy (patch, PANEL_VAL (panel, 0), PANEL_VAL (panel, 1), 0);
  gdouble yi = patch_eval_dx (patch, PANEL_VAL (panel, 0), PANEL_VAL (panel, 1), 1);
  gdouble yj = patch_eval_dy (patch, PANEL_VAL (panel, 0), PANEL_VAL (panel, 1), 1);
  gdouble zi = patch_eval_dx (patch, PANEL_VAL (panel, 0), PANEL_VAL (panel, 1), 2);
  gdouble zj = patch_eval_dy (patch, PANEL_VAL (panel, 0), PANEL_VAL (panel, 1), 2);
  
  Vector g;
  g.x = yi*zj-yj*zi;
  g.y = zi*xj-zj*xi;
  g.z = xi*yj-xj*yi;

  PANEL_VAL (panel, 9) = vector_norm (g);
}

void patch_compute_metric (Patch * patch)
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

  if (hz != NULL) {
    /* Check if we are under the free-surface */
    gdouble dz[4];
    dz[0] = patch_eval (patch, il, jl, 2)
      - hz (patch_eval (patch, il, jl, 0), patch_eval (patch, il, jl, 1), 0, NULL);
    dz[1] = patch_eval (patch, ih, jl, 2)
      - hz (patch_eval (patch, ih, jl, 0), patch_eval (patch, ih, jl, 1), 0, NULL);
    dz[2] = patch_eval (patch, ih, jh, 2)
      - hz (patch_eval (patch, ih, jh, 0), patch_eval (patch, ih, jh, 1), 0, NULL);
    dz[3] = patch_eval (patch, il, jh, 2)
      - hz (patch_eval (patch, il, jh, 0), patch_eval (patch, il, jh, 1), 0, NULL);
    
    if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.))
      return v;
  }

  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  /* Check if panel size is reached */
  if ( fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance ) {
    /* Surface integral of var over the small piece of panel */
    /* var=3 == metric */
    return patch_eval (patch, (il+ih)/2., (jl+jh)/2., var)
      *patch_eval (patch, (il+ih)/2., (jl+jh)/2., 3)*fabs(ih-il)*fabs(jh-jl);
  }

  /* We refine */
  v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jl, jhalf);
  v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jl, jhalf);
  v += adaptive_panel_scalar_integral (patch, hz, var, il, ihalf, jhalf, jh);
  v += adaptive_panel_scalar_integral (patch, hz, var, ihalf, ih, jhalf, jh);

  return v;
}

gdouble patch_fast_panel_scalar_integral (Patch * patch, HeightCurve hz, gint var)
{
  g_assert (patch != NULL);
  gint i, j;
  gdouble sum = 0.;

  for ( j = 0; j < patch->rows->len ; j++) {
    GArray * row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len; i++) {
      sum += fast_panel_scalar_integral (patch, hz, var, i-0.5, i+0.5, j-0.5, j+0.5);
    }
  }
  return sum;
}

typedef struct {
  Patch * patch;
  HeightCurve hz;
} VTMParams;

static void panel_compute_volume_times_metric (Panel * panel, gpointer data)
{
  VTMParams * par = (VTMParams *) data;

  if (par->hz == NULL)
    PANEL_VAL (panel, 9) = PANEL_VAL (panel, 9);
  else
    PANEL_VAL (panel, 10) = adaptive_panel_scalar_integral (par->patch, par->hz, 3, panel->i-0.5,
							    panel->i+0.5, panel->j-0.5, panel->j+0.5);
}

void patch_compute_volume_times_metric (Patch * patch, HeightCurve hz)
{
  VTMParams par;
  par.patch = patch;
  par.hz = hz;

  patch_forall_panels (patch, panel_compute_volume_times_metric, &par);
}



void patch_store_derivatives (Patch * patch)
{
  gint i, j;
  GArray * row = g_ptr_array_index (patch->rows, 0);
  GArray * row0, * row2;
  
  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 1; i < row->len-1; i++) {
      Panel p = g_array_index (row, Panel, i);
      Panel p2 = g_array_index (row, Panel, i+1);
      Panel p0 = g_array_index (row, Panel, i-1);
      
      p.var[3] = (p2.var[0] - p0.var[0])/2.;
      p.var[5] = (p2.var[1] - p0.var[1])/2.;
      p.var[7] = (p2.var[2] - p0.var[2])/2.;

      g_array_index (row, Panel, i) = p;
    }

    i = 0;
    {
      Panel p = g_array_index (row, Panel, i);
      Panel p1 = g_array_index (row, Panel, i+1);
      Panel p2 = g_array_index (row, Panel, i+2);
      
      /* p.var[3] = p1.var[3] - (p2.var[0] - p.var[0])/2. + 0.5*(p2.var[0] + p.var[0] - 2.*p1.var[0]); */
      /* p.var[5] = p1.var[5] - (p2.var[1] - p.var[1])/2. + 0.5*(p2.var[1] + p.var[1] - 2.*p1.var[1]); */
      /* p.var[7] = p1.var[7] - (p2.var[2] - p.var[2])/2. + 0.5*(p2.var[2] + p.var[2] - 2.*p1.var[2]); */

      p.var[3] = p1.var[3] - (p2.var[3] - p1.var[3]) ;
      p.var[5] = p1.var[5] ;
      p.var[7] = p1.var[7] ;

      g_array_index (row, Panel, i) = p;
    }

    i = row->len-1;
    {
      Panel p = g_array_index (row, Panel, i);
      Panel p1 = g_array_index (row, Panel, i-1);
      Panel p0 = g_array_index (row, Panel, i-2);
      
      p.var[3] = p1.var[3] + (p.var[3] - p0.var[3])/2. /* + 0.5*(p.var[0] + p0.var[0] - 2.*p1.var[0]) */;
      p.var[5] = p1.var[5] + (p.var[1] - p0.var[1])/2. /* + 0.5*(p.var[1] + p0.var[1] - 2.*p1.var[1]) */;
      p.var[7] = p1.var[7] + (p.var[2] - p0.var[2])/2. /* + 0.5*(p.var[2] + p0.var[2] - 2.*p1.var[2]) */;

      g_array_index (row, Panel, i) = p;
    }
  }

  
  for ( j = 1; j < patch->rows->len-1; j++) {
    row = g_ptr_array_index (patch->rows, j);
    row0 = g_ptr_array_index (patch->rows, j-1);
    row2 = g_ptr_array_index (patch->rows, j+1);
    for ( i = 0; i < row->len; i++) {
      Panel p = g_array_index (row, Panel, i);
      Panel p2 = g_array_index (row2, Panel, i);
      Panel p0 = g_array_index (row0, Panel, i);
      
      p.var[4] = (p2.var[0] - p0.var[0])/2.;
      p.var[6] = (p2.var[1] - p0.var[1])/2.;
      p.var[8] = (p2.var[2] - p0.var[2])/2.;

      g_array_index (row, Panel, i) = p;
    }
  }

  row = g_ptr_array_index (patch->rows, 0);
  row2 = g_ptr_array_index (patch->rows, 1);
  for ( i = 0; i < row->len; i++) {
    Panel p = g_array_index (row, Panel, i);
    Panel p2 = g_array_index (row2, Panel, i);
    
    p.var[4] = (p2.var[0] - p.var[0]);
    p.var[6] = (p2.var[1] - p.var[1]);
    p.var[8] = (p2.var[2] - p.var[2]);
    
    g_array_index (row, Panel, i) = p;
  }

  row = g_ptr_array_index (patch->rows, patch->rows->len-1);
  row0 = g_ptr_array_index (patch->rows, patch->rows->len-2);
  for ( i = 0; i < row->len; i++) {
    Panel p = g_array_index (row, Panel, i);
    Panel p0 = g_array_index (row0, Panel, i);
    
    p.var[4] = (p.var[0] - p0.var[0]);
    p.var[6] = (p.var[1] - p0.var[1]);
    p.var[8] = (p.var[2] - p0.var[2]);
    
    g_array_index (row, Panel, i) = p;
  }
}

//****************************************************************************************//
Spline2D * spline2d_fit_geometry2_old (Patch * patch, gint M, gint N, gint order, gint inner, gint outer, gboolean flip)
{
  Spline2D * sp = spline2d_new (M, N, order, inner, outer);

  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint NX = row->len+1;
  gint NY = patch->rows->len+1;
  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  Point p[NX][NY];
  gint i, j, k, l;

  /* Copy data to p */
  for ( j = 0; j <  patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len ; i++) {
      Panel panel = g_array_index (row, Panel, i);
      p[i][j] = panel.p[2];
    }
  }

  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    i = row->len-1;
    Panel panel = g_array_index (row, Panel, i);
    p[row->len][j] = panel.p[1];
  }

  j = patch->rows->len-1;
  row = g_ptr_array_index (patch->rows, j);
  for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      p[i][patch->rows->len] = panel.p[3];
  }

  i = row->len-1;
  Panel panel = g_array_index (row, Panel, i);
  p[row->len][patch->rows->len] = panel.p[0];


  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);

  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;

  cx = gsl_vector_alloc (NU*NV);
  cy = gsl_vector_alloc (NU*NV);
  cz = gsl_vector_alloc (NU*NV);
  cov_x = gsl_matrix_alloc (NU*NV, NU*NV);
  cov_y = gsl_matrix_alloc (NU*NV, NU*NV);
  cov_z = gsl_matrix_alloc (NU*NV, NU*NV);

  X = gsl_matrix_alloc (NX*NY, NU*NV);
  x = gsl_vector_alloc (NX*NY);
  y = gsl_vector_alloc (NX*NY);
  z = gsl_vector_alloc (NX*NY);

  /* Initial guess */
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      gint m = flip ? NX - 1 -i / (NU -1.)*(NX-1.) : i / (NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, i+ j*NU, p[m][n].x);
      gsl_vector_set (cy, i+ j*NU, p[m][n].y);
      gsl_vector_set (cz, i+ j*NU, p[m][n].z);
    }
  
  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      for ( k = 0; k < NX; k++)
  	for ( l = 0; l < NY; l++) {
	  if (flip)
	    gsl_bspline_eval (1.-k/(NX-1.), Bu, sp->w_u);
	  else
	    gsl_bspline_eval (k/(NX-1.), Bu, sp->w_u);
  	  gsl_bspline_eval (l/(NY-1.), Bv, sp->w_v);
  	  gsl_matrix_set (X, k + NX*l, i + j*NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      gsl_vector_set (x, i + j*NX, p[i][j].x);
      gsl_vector_set (y, i + j*NX, p[i][j].y);
      gsl_vector_set (z, i + j*NX, p[i][j].z);
    }
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (NX*NY, NU*NV);

  /* Multidimensional fitting usind singular values decomposition method */
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, 1e-6, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, 1e-6, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, 1e-6, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i+j*NU));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i+j*NU));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i+j*NU));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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
  return sp;
}

Spline2D * spline2d_fit_geometry2 (Patch * patch, gint M, gint N, gint order, gint inner, gint outer, gboolean flip_u, gboolean flip_v, gboolean centripetal_reparam, gboolean swap)
{
  gdouble tolerance_border = 1e-2;
  gdouble tolerance_inside = 1e-3;

  Spline2D * sp = swap ? spline2d_new (N, M, order, inner, outer) : spline2d_new (M, N, order, inner, outer);

  spline2d_init_panels (sp);

  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint NX = swap ? patch->rows->len+1 : row->len+1;
  gint NY = swap ? row->len+1 : patch->rows->len+1;

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  //Point p[NX][NY];
  Point ptmp[NX][NY];
  gint i, j, k, l;

  /* Copy data to p */
  for ( j = 0; j <  patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len ; i++) {
      Panel panel = g_array_index (row, Panel, i);
      if (swap)
	ptmp[j][i] = panel.p[2];
      else
	ptmp[i][j] = panel.p[2];
    }
  }

  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    i = row->len-1;
    Panel panel = g_array_index (row, Panel, i);
    if (swap)
      ptmp[j][row->len] = panel.p[1];
    else
      ptmp[row->len][j] = panel.p[1];
  }

  j = patch->rows->len-1;
  row = g_ptr_array_index (patch->rows, j);
  for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      if (swap)
	ptmp[patch->rows->len][i] = panel.p[3];
      else
	ptmp[i][patch->rows->len] = panel.p[3];
  }

  i = row->len-1;
  Panel panel = g_array_index (row, Panel, i);
  if (swap)
    ptmp[patch->rows->len][row->len] = panel.p[0];
  else
    ptmp[row->len][patch->rows->len] = panel.p[0];

  for (  j = 0; j < NY; j++) {
    
    gboolean water = FALSE;
    for ( i = 0; i < NX; i++) {
      if (ptmp[i][j].z < 0)
	water = TRUE;
    }

    if ( water == FALSE ) {
      int a, b;
      for ( b = j+1; b < NY; b++ ) {
	for ( a = 0; a < NX; a++ ) {
	  ptmp[a][b-1]=ptmp[a][b];
	}
      }
      NY--; j--;
    }
  }


  /***********************************************************************/

  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      fprintf (fp, "%f %f %f %i %i\n", ptmp[i][j].x, ptmp[i][j].y, ptmp[i][j].z, i, j);
    }
  }
  fclose (fp);

  gint NTMP = NX;
  NX = NY; NY = NTMP;
  
  Point p[NX][NY], p2[NX][NY];
  gdouble u[NX][NY], v[NX][NY];

  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      p[i][j] = ptmp[NY-1-j][i];
    }
  }

  for ( i = 0; i < NX; i++) {
    gint ii = flip_u ? NX-1-i : i;
    for ( j = 0; j < NY; j++) {
      gint jj = flip_v ? NY-1-j : j;
      p2[i][j] = p[ii][jj];
    }
  }
  //  g_assert_not_reached ();

  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      p[i][j] = p2[i][j];
      
    }
  }



  /* g_warning ("Tweek to ensure that the hull is perfectly closed\n"); */
  /* for ( j = 0; j < NY; j++) { */
  /*   i = 0; */
  /*   p[i][j].y = 0.; */
  /*   // fprintf (stdout, "%f %f %f \n", p[i][j].x, p[i][j].y, p[i][j].z); */
  /*   i = NX-1; */
  /*   p[i][j].y = 0.; */
  /*   //fprintf (stdout, "%f %f %f \n", p[i][j].x, p[i][j].y, p[i][j].z); */
  /* } */

  /* for ( i = 0; i < NX; i++) { */
  /*   j = NY-1; */
  /*   p[i][j].y = 0.; */
  /*   //fprintf (stdout, "%f %f %f \n", p[i][j].x, p[i][j].y, p[i][j].z); */

  /*   /\* j = 0; *\/ */
  /*   /\* fprintf (stdout, "%f %f %f \n", p[i][j].x, p[i][j].y, p[i][j].z); *\/ */
  /* } */


  // Centripetal reparametrization of the whole surface
  if (centripetal_reparam) {
    for ( i = 0; i < NX; i++) {
      gdouble lv = 0.;
      v[i][0] = lv;
      for ( j = 1; j < NY; j++) {
	lv += point_distance (p[i][j], p[i][j-1]);
	v[i][j] = lv;
      }
      for ( j = 1; j < NY; j++) {
	if (lv != 0.)
		v[i][j] /= lv;
	else
	v[i][j] = j/(NY-1.);
      }
    }

    for ( j = 0; j < NY; j++) {
      gdouble lu = 0.;
      u[0][j] = 0.;
      for ( i = 1; i < NX; i++) {
	lu += point_distance (p[i][j], p[i-1][j]);
	u[i][j] = lu;
      }
      for ( i = 1; i < NX; i++) {
	if (lu != 0.)
	  u[i][j] /= lu;
	else
	  u[i][j] = i/(NX-1.);
      }
    }
  }
  else {
    for ( j = 0; j < NY; j++) {
      for ( i = 0; i < NX; i++) {
	  u[i][j] = i/(NX-1.);
	v[i][j] = j/(NY-1.);
      }
    }
  }



  // Blend of original and centripetal parametrisations
  for ( j = 0; j < NY; j++) {
    for ( i = 0; i < NX; i++) {
      u[i][j] = 0.4*u[i][j] + 0.6*i/(NX-1.);
      v[i][j] = 0.7*v[i][j]+0.3*j/(NY-1.);
    }
  }


  // Set corners
  gint u0 = 0;
  gint u1 = NU-1;
  gint v0 = 0;
  gint v1 = NV-1;

  coeff_assign (sp, u0, v0, 0, p[0][0].x);
  coeff_assign (sp, u0, v0, 1, p[0][0].y);
  coeff_assign (sp, u0, v0, 2, p[0][0].z);
    
  coeff_assign (sp, u1, v0, 0, p[NX-1][0].x);
  coeff_assign (sp, u1, v0, 1, p[NX-1][0].y);
  coeff_assign (sp, u1, v0, 2, p[NX-1][0].z);

  coeff_assign (sp, u1, v1, 0, p[NX-1][NY-1].x);
  coeff_assign (sp, u1, v1, 1, p[NX-1][NY-1].y);
  coeff_assign (sp, u1, v1, 2, p[NX-1][NY-1].z);

  coeff_assign (sp, u0, v1, 0, p[0][NY-1].x);
  coeff_assign (sp, u0, v1, 1, p[0][NY-1].y);
  coeff_assign (sp, u0, v1, 2, p[0][NY-1].z);



  // Fit borders
  

  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;

  cx = gsl_vector_alloc (2*(NU+NV)-8);
  cy = gsl_vector_alloc (2*(NU+NV)-8);
  cz = gsl_vector_alloc (2*(NU+NV)-8);
  cov_x = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_y = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_z = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);

  X = gsl_matrix_alloc (2*(NX+NY)-8, 2*(NU+NV)-8);
  x = gsl_vector_alloc (2*(NX+NY)-8);
  y = gsl_vector_alloc (2*(NX+NY)-8);
  z = gsl_vector_alloc (2*(NX+NY)-8);

  /* Initial guess */ /* Could be improved using a spline fit at first */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    gint m = i/(NU -1.)*(NX-1.);
    gint n = 0;
    gsl_vector_set (cx, i-1, p[m][n].x);
    gsl_vector_set (cy, i-1, p[m][n].y);
    gsl_vector_set (cz, i-1, p[m][n].z);

    j = NV-1;
    n = /* flip_v ? 0 : */ (NY-1.);
    gsl_vector_set (cx, (i-1)+(NU-2), p[m][n].x);
    gsl_vector_set (cy, (i-1)+(NU-2), p[m][n].y);
    gsl_vector_set (cz, (i-1)+(NU-2), p[m][n].z);
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;
    gint m = 0;
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (j-1)+2*(NU-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2), p[m][n].z);

    i = NU-1;
    m = NX - 1;
    gsl_vector_set (cx, (j-1)+2*(NU-2)+(NV-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2)+(NV-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2)+(NV-2), p[m][n].z);
  }

  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++) {

    j = 0;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }

    j = NV-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    
    i = NU-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    j = 0;
    gsl_vector_set (x, (i-1), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    j = NY-1;
    gsl_vector_set (x, (i-1) + (NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1) + (NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1) + (NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  for ( j = 1; j < NY-1; j++) {
    i = 0;
    gsl_vector_set (x, (j-1)+2*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    i = NX-1;
    gsl_vector_set (x, (j-1)+2*(NX-2)+(NY-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2)+(NY-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2)+(NY-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (2.*(NX+NY)-8, 2*(NU+NV)-8);
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, tolerance_border, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_border, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_border, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1));
    
    j = NV-1;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1+(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1+(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1+(NU-2)));
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)));

    i = NU-1;
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)+(NV-2)));
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


  cx = gsl_vector_alloc ((NU-2)*(NV-2));
  cy = gsl_vector_alloc ((NU-2)*(NV-2));
  cz = gsl_vector_alloc ((NU-2)*(NV-2));
  cov_x = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_y = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_z = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));

  X = gsl_matrix_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));
  x = gsl_vector_alloc ((NX-2)*(NY-2));
  y = gsl_vector_alloc ((NX-2)*(NY-2));
  z = gsl_vector_alloc ((NX-2)*(NY-2));

  /* Initial guess */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      gint m = i/(NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, (i-1)+ (j-1)*(NU-2), p[m][n].x);
      gsl_vector_set (cy, (i-1)+ (j-1)*(NU-2), p[m][n].y);
      gsl_vector_set (cz, (i-1)+ (j-1)*(NU-2), p[m][n].z);
    }
  
  /* Set problem */
  Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      for ( k = 1; k < NX-1; k++)
  	for ( l = 1; l < NY-1; l++) {
  	  gsl_bspline_eval (u[k][l], Bu, sp->w_u);
  	  gsl_bspline_eval (v[k][l], Bv, sp->w_v);
  	  gsl_matrix_set (X, (k-1) + (NX-2)*(l-1), (i-1) + (j-1)*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    for ( j = 1; j < NY-1; j++) {
      gsl_vector_set (x, (i-1) + (j-1)*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
      gsl_vector_set (y, (i-1) + (j-1)*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
      gsl_vector_set (z, (i-1) + (j-1)*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
    }
  }

  work = gsl_multifit_linear_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));

  /* Multidimensional fitting usind singular values decomposition method */
  gsl_multifit_linear_svd (X, x, tolerance_inside, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_inside, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_inside, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (i-1)+(j-1)*(NU-2)));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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

  spline2d_reinit_panels_physical_quantities (sp);

  return sp;
}

Spline2D * spline2d_interpolate_geometry (Patch * patch, gint M, gint N, gint order, gint inner, gint outer, gboolean flip_u, gboolean flip_v)
{
  gdouble tolerance_border = 1e-2;
  gdouble tolerance_inside = 1e-3;
  Spline2D * sp = spline2d_new (M, N, order, inner, outer);

  spline2d_init_panels (sp);

  GArray * row = g_ptr_array_index (patch->rows, 0);
  gint NX = row->len+1;
  gint NY = patch->rows->len+1;

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  //Point p[NX][NY];
  Point ptmp[NX][NY];
  gint i, j, k, l;

  /* Copy data to p */
  for ( j = 0; j <  patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    for ( i = 0; i < row->len ; i++) {
      Panel panel = g_array_index (row, Panel, i);
      ptmp[i][j] = panel.p[2];
    }
  }

  for ( j = 0; j < patch->rows->len; j++) {
    row = g_ptr_array_index (patch->rows, j);
    i = row->len-1;
    Panel panel = g_array_index (row, Panel, i);
    ptmp[row->len][j] = panel.p[1];
  }

  j = patch->rows->len-1;
  row = g_ptr_array_index (patch->rows, j);
  for ( i = 0; i < row->len; i++) {
      Panel panel = g_array_index (row, Panel, i);
      ptmp[i][patch->rows->len] = panel.p[3];
  }

  i = row->len-1;
  Panel panel = g_array_index (row, Panel, i);
  ptmp[row->len][patch->rows->len] = panel.p[0];


  /* FILE * fp = fopen ("points2fit.dat","w"); */
  /* for ( i = 0; i < NX; i++) { */
  /*   for ( j = 0; j < NY; j++) { */
  /*     fprintf (fp, "%f %f %f %i %i\n", ptmp[i][j].x, ptmp[i][j].y, ptmp[i][j].z, i, j); */
  /*   } */
  /* } */
  /* fclose (fp); */

  gint NTMP = NX;
  NX = NY; NY = NTMP;
  
  Point p[NX][NY], p2[NX][NY];
  gdouble u[NX][NY], v[NX][NY];

  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      p[i][j] = ptmp[NY-1-j][i];
    }
  }

  for ( i = 0; i < NX; i++) {
    gint ii = flip_u ? NX-1-i : i;
    for ( j = 0; j < NY; j++) {
      gint jj = flip_v ? NY-1-j : j;
      p2[i][j] = p[ii][jj];
    }
  }
  //  g_assert_not_reached ();

  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      p[i][j] = p2[i][j];
    }
  }

  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);

  g_assert_not_reached ();



  /* coeff_assign (sp, u0, v0, 0, p[0][0].x); */
  /* coeff_assign (sp, u0, v0, 1, p[0][0].y); */
  /* coeff_assign (sp, u0, v0, 2, p[0][0].z); */
    
 

  spline2d_reinit_panels_physical_quantities (sp);

  return sp;
}

Spline2D * spline2d_fit_geometry3 (FILE * finput, gint M, gint N, gint order, gint inner, gint outer, gint NX, gint NY)
{
  Spline2D * sp = spline2d_new (M, N, order, inner, outer);

  // GArray * row = g_ptr_array_index (patch->rows, 0);
  //gint NX = row->len+1;
  //gint NY = patch->rows->len+1;

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  Point p[NX][NY];
  gint i, j, k, l;

  float px, py, pz;
  i = j = 0;
  while (fscanf(finput, "%f %f %f\n", &px, &py, &pz) != EOF) {
    p[i][j].x = px;
    p[i][j].y = py;
    p[i][j].z = pz;

    i++;
    if ( i == NX ) {
      i = 0;
      j++;
    }
  }


  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);
  //g_assert_not_reached ();
  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;

  cx = gsl_vector_alloc (NU*NV);
  cy = gsl_vector_alloc (NU*NV);
  cz = gsl_vector_alloc (NU*NV);
  cov_x = gsl_matrix_alloc (NU*NV, NU*NV);
  cov_y = gsl_matrix_alloc (NU*NV, NU*NV);
  cov_z = gsl_matrix_alloc (NU*NV, NU*NV);

  X = gsl_matrix_alloc (NX*NY, NU*NV);
  x = gsl_vector_alloc (NX*NY);
  y = gsl_vector_alloc (NX*NY);
  z = gsl_vector_alloc (NX*NY);

  /* Initial guess */
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      gint m = i / (NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, i+ j*NU, p[m][n].x);
      gsl_vector_set (cy, i+ j*NU, p[m][n].y);
      gsl_vector_set (cz, i+ j*NU, p[m][n].z);
    }
  
  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      for ( k = 0; k < NX; k++)
  	for ( l = 0; l < NY; l++) {
  	  gsl_bspline_eval (k/(NX-1.), Bu, sp->w_u);
  	  gsl_bspline_eval (l/(NY-1.), Bv, sp->w_v);
  	  gsl_matrix_set (X, k + NX*l, i + j*NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      gsl_vector_set (x, i + j*NX, p[i][j].x);
      gsl_vector_set (y, i + j*NX, p[i][j].y);
      gsl_vector_set (z, i + j*NX, p[i][j].z);
    }
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (NX*NY, NU*NV);

  /* Multidimensional fitting usind singular values decomposition method */
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, 1e-6, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, 1e-6, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, 1e-6, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i+j*NU));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i+j*NU));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i+j*NU));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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

  spline2d_init_panels (sp);

  return sp;
}

/**
 * This fitting routine does a best fit of the borders and then a best
 * fitting of the inside of the patch.
 * The quality of the grid is sensitive to the tolerance parameter imposed
 * the svd fitting algorithm.
 * It is also possible to reparametrize the surface using the centripetal
 * approach prior to doing the fit.
 **/
Spline2D * spline2d_fit_geometry4 (FILE * finput, gint M, gint N, gint order, gint inner, gint outer, gint NX, gint NY)
{
  gdouble tolerance_border = 5e-3;
  gdouble tolerance_inside = 1e-3;
  gdouble centripetal_reparam = FALSE;
  Spline2D * sp = spline2d_new (M, N, order, inner, outer);
  spline2d_init_panels (sp);
  coeff_set_var_to_zero (sp, 0);
  coeff_set_var_to_zero (sp, 1);
  coeff_set_var_to_zero (sp, 2);

  // GArray * row = g_ptr_array_index (patch->rows, 0);
  //gint NX = row->len+1;
  //gint NY = patch->rows->len+1;

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  Point p[NX][NY];
  gdouble u[NX][NY], v[NX][NY];
  gint i, j, k, l;

  float px, py, pz;
  i = j = 0;
  while (fscanf(finput, "%f %f %f\n", &px, &py, &pz) != EOF) {
    p[i][j].x = px;
    p[i][j].y = py;
    p[i][j].z = pz;

    i++;
    if ( i == NX ) {
      i = 0;
      j++;
    }
  }

  // Centripetal reparametrization of the whole surface
  if (centripetal_reparam) {
    for ( i = 0; i < NX; i++) {
      gdouble lv = 0.;
      v[i][0] = lv;
      for ( j = 1; j < NY; j++) {
	lv += point_distance (p[i][j], p[i][j-1]);
	v[i][j] = lv;
      }
      for ( j = 1; j < NY; j++) {
	if (lv != 0.)
		v[i][j] /= lv;
	else
	v[i][j] = j/(NY-1.);
      }
    }

    for ( j = 0; j < NY; j++) {
      gdouble lu = 0.;
      u[0][j] = 0.;
      for ( i = 1; i < NX; i++) {
	lu += point_distance (p[i][j], p[i-1][j]);
	u[i][j] = lu;
      }
      for ( i = 1; i < NX; i++) {
	if (lu != 0.)
	  u[i][j] /= lu;
	else
	  u[i][j] = i/(NX-1.);
      }
    }
  }
  else {
    for ( j = 0; j < NY; j++) {
      for ( i = 0; i < NX; i++) {
	u[i][j] = i/(NX-1.);
	v[i][j] = j/(NY-1.);
      }
    }
  }

  
  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);


  // Fit borders
  

  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;

  cx = gsl_vector_alloc (2*(NU+NV)-4);
  cy = gsl_vector_alloc (2*(NU+NV)-4);
  cz = gsl_vector_alloc (2*(NU+NV)-4);
  cov_x = gsl_matrix_alloc (2*(NU+NV)-4, 2*(NU+NV)-4);
  cov_y = gsl_matrix_alloc (2*(NU+NV)-4, 2*(NU+NV)-4);
  cov_z = gsl_matrix_alloc (2*(NU+NV)-4, 2*(NU+NV)-4);

  X = gsl_matrix_alloc (2*(NX+NY)-4, 2*(NU+NV)-4);
  x = gsl_vector_alloc (2*(NX+NY)-4);
  y = gsl_vector_alloc (2*(NX+NY)-4);
  z = gsl_vector_alloc (2*(NX+NY)-4);

  /* Initial guess */
  for ( i = 0; i < NU; i++) {
    j = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, i, p[m][n].x);
    gsl_vector_set (cy, i, p[m][n].y);
    gsl_vector_set (cz, i, p[m][n].z);

    j = NV-1;
    m = i / (NU -1.)*(NX-1.);
    n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, i+NU, p[m][n].x);
    gsl_vector_set (cy, i+NU, p[m][n].y);
    gsl_vector_set (cz, i+NU, p[m][n].z);
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (j-1)+2*NU, p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*NU, p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*NU, p[m][n].z);

    i = NU-1;
    gsl_vector_set (cx, (j-1)+2*NU+(NV-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*NU+(NV-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*NU+(NV-2), p[m][n].z);
  }

  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 0; i < NU; i++) {
    j = 0;

    for ( k = 0; k < NX; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k , i, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k + NX, i, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX, i, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX + NY-2, i, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }

    j = NV-1;

    for ( k = 0; k < NX; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k , i + NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k + NX, i + NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX, i + NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX + (NY-2), i + NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;
    for ( k = 0; k < NX; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k , (j-1)+2*NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k + NX, (j-1)+2*NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX, (j-1)+2*NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX + (NY-2), (j-1)+2*NU, gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    
    i = NU-1;
    for ( k = 0; k < NX; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k , (j-1)+2*NU+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, k + NX, (j-1)+2*NU+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX, (j-1)+2*NU+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*NX + (NY-2), (j-1)+2*NU+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 0; i < NX; i++) {
    j = 0;
    gsl_vector_set (x, i, p[i][j].x);
    gsl_vector_set (y, i, p[i][j].y);
    gsl_vector_set (z, i, p[i][j].z);

    j = NY-1;
    gsl_vector_set (x, i + NX, p[i][j].x);
    gsl_vector_set (y, i + NX, p[i][j].y);
    gsl_vector_set (z, i + NX, p[i][j].z);
  }

  for ( j = 1; j < NY-1; j++) {
    i = 0;
    gsl_vector_set (x, (j-1)+2*NX, p[i][j].x);
    gsl_vector_set (y, (j-1)+2*NX, p[i][j].y);
    gsl_vector_set (z, (j-1)+2*NX, p[i][j].z);

    i = NX-1;
    gsl_vector_set (x, (j-1)+2*NX+(NY-2), p[i][j].x);
    gsl_vector_set (y, (j-1)+2*NX+(NY-2), p[i][j].y);
    gsl_vector_set (z, (j-1)+2*NX+(NY-2), p[i][j].z);
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (2.*(NX+NY)-4, 2*(NU+NV)-4);
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, tolerance_border, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_border, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_border, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 0; i < NU; i++) {
    j = 0;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i));
    
    j = NV-1;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i+NU));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i+NU));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i+NU));
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (j-1)+2*NU));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (j-1)+2*NU));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (j-1)+2*NU));

    i = NU-1;
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (j-1)+2*NU+(NV-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (j-1)+2*NU+(NV-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (j-1)+2*NU+(NV-2)));
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


  cx = gsl_vector_alloc ((NU-2)*(NV-2));
  cy = gsl_vector_alloc ((NU-2)*(NV-2));
  cz = gsl_vector_alloc ((NU-2)*(NV-2));
  cov_x = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_y = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_z = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));

  X = gsl_matrix_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));
  x = gsl_vector_alloc ((NX-2)*(NY-2));
  y = gsl_vector_alloc ((NX-2)*(NY-2));
  z = gsl_vector_alloc ((NX-2)*(NY-2));

  /* Initial guess */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      gint m = i / (NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, (i-1)+ (j-1)*(NU-2), p[m][n].x);
      gsl_vector_set (cy, (i-1)+ (j-1)*(NU-2), p[m][n].y);
      gsl_vector_set (cz, (i-1)+ (j-1)*(NU-2), p[m][n].z);
    }
  
  /* Set problem */
  Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      for ( k = 1; k < NX-1; k++)
  	for ( l = 1; l < NY-1; l++) {
  	  gsl_bspline_eval (u[k][l], Bu, sp->w_u);
  	  gsl_bspline_eval (v[k][l], Bv, sp->w_v);
  	  gsl_matrix_set (X, (k-1) + (NX-2)*(l-1), (i-1) + (j-1)*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    for ( j = 1; j < NY-1; j++) {
      gsl_vector_set (x, (i-1) + (j-1)*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
      gsl_vector_set (y, (i-1) + (j-1)*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
      gsl_vector_set (z, (i-1) + (j-1)*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
    }
  }

  work = gsl_multifit_linear_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));

  /* Multidimensional fitting usind singular values decomposition method */
  gsl_multifit_linear_svd (X, x, tolerance_inside, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_inside, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_inside, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (i-1)+(j-1)*(NU-2)));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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

  spline2d_reinit_panels_physical_quantities (sp);

  return sp;
}

/**
 * This fitting routine enforces the corners of the patch to be exactly that
 * in the input file, then does a best fit of the borders and finally a best
 * fitting of the inside of the patch.
 * The quality of the grid is sensitive to the tolerance parameter imposed
 * the svd fitting algorithm.
 * It is also possible to reparametrize the surface using the centripetal
 * approach prior to doing the fit.
 **/
Spline2D * spline2d_fit_geometry6 (GArray * input, gint M, gint N, gint order, gint inner, gint outer, gint NX, gint NY)
{
  gdouble tolerance_border = 1e-2;
  gdouble tolerance_inside = 1e-6;
  gdouble centripetal_reparam = FALSE;

  Spline2D * sp = spline2d_new (M, N, order, inner, outer);
  spline2d_init_panels (sp);
  coeff_set_var_to_zero (sp, 0);
  coeff_set_var_to_zero (sp, 1);
  coeff_set_var_to_zero (sp, 2);

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  Point p[NX][NY];
  gdouble u[NX][NY], v[NX][NY];
  gint i, j, k, l;

  // Read data
  for ( j = 0; j < NY; j++ ) {
    for ( i = 0; i < NX; i++ ) {
      p[i][j] = g_array_index (input, Point, j + i*NY); 
    }
  } 

  // Centripetal reparametrization of the whole surface
  if (centripetal_reparam) {
    for ( i = 0; i < NX; i++) {
      gdouble lv = 0.;
      v[i][0] = lv;
      for ( j = 1; j < NY; j++) {
	lv += point_distance (p[i][j], p[i][j-1]);
	v[i][j] = lv;
      }
      for ( j = 1; j < NY; j++) {
	if (lv != 0.)
		v[i][j] /= lv;
	else
	v[i][j] = j/(NY-1.);
      }
    }

    for ( j = 0; j < NY; j++) {
      gdouble lu = 0.;
      u[0][j] = 0.;
      for ( i = 1; i < NX; i++) {
	lu += point_distance (p[i][j], p[i-1][j]);
	u[i][j] = lu;
      }
      for ( i = 1; i < NX; i++) {
	if (lu != 0.)
	  u[i][j] /= lu;
	else
	  u[i][j] = i/(NX-1.);
      }
    }
  }
  else {
    for ( j = 0; j < NY; j++) {
      for ( i = 0; i < NX; i++) {
	u[i][j] = i/(NX-1.);
	v[i][j] = j/(NY-1.);
      }
    }
  }

  
  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    /* for ( j = 0; j < NY; j++) */ {
      j = 0;
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);

  // Set corners
  coeff_assign (sp, 0, 0, 0, p[0][0].x);
  coeff_assign (sp, 0, 0, 1, p[0][0].y);
  coeff_assign (sp, 0, 0, 2, p[0][0].z);

  coeff_assign (sp, NU-1, 0, 0, p[NX-1][0].x);
  coeff_assign (sp, NU-1, 0, 1, p[NX-1][0].y);
  coeff_assign (sp, NU-1, 0, 2, p[NX-1][0].z);

  coeff_assign (sp, NU-1, NV-1, 0, p[NX-1][NY-1].x);
  coeff_assign (sp, NU-1, NV-1, 1, p[NX-1][NY-1].y);
  coeff_assign (sp, NU-1, NV-1, 2, p[NX-1][NY-1].z);

  coeff_assign (sp, 0, NV-1, 0, p[0][NY-1].x);
  coeff_assign (sp, 0, NV-1, 1, p[0][NY-1].y);
  coeff_assign (sp, 0, NV-1, 2, p[0][NY-1].z);

  // Fit borders
  

  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;

  cx = gsl_vector_alloc (2*(NU+NV)-8);
  cy = gsl_vector_alloc (2*(NU+NV)-8);
  cz = gsl_vector_alloc (2*(NU+NV)-8);
  cov_x = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_y = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_z = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);

  X = gsl_matrix_alloc (2*(NX+NY)-8, 2*(NU+NV)-8);
  x = gsl_vector_alloc (2*(NX+NY)-8);
  y = gsl_vector_alloc (2*(NX+NY)-8);
  z = gsl_vector_alloc (2*(NX+NY)-8);

  /* Initial guess */ /* Could be improved using a spline fit at first */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, i-1, p[m][n].x);
    gsl_vector_set (cy, i-1, p[m][n].y);
    gsl_vector_set (cz, i-1, p[m][n].z);

    j = NV-1;
    m = i / (NU -1.)*(NX-1.);
    n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (i-1)+(NU-2), p[m][n].x);
    gsl_vector_set (cy, (i-1)+(NU-2), p[m][n].y);
    gsl_vector_set (cz, (i-1)+(NU-2), p[m][n].z);
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (j-1)+2*(NU-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2), p[m][n].z);

    i = NU-1;
    gsl_vector_set (cx, (j-1)+2*(NU-2)+(NV-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2)+(NV-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2)+(NV-2), p[m][n].z);
  }

  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++) {

    j = 0;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }

    j = NV-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    
    i = NU-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    j = 0;
    gsl_vector_set (x, (i-1), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    j = NY-1;
    gsl_vector_set (x, (i-1) + (NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1) + (NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1) + (NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  for ( j = 1; j < NY-1; j++) {
    i = 0;
    gsl_vector_set (x, (j-1)+2*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    i = NX-1;
    gsl_vector_set (x, (j-1)+2*(NX-2)+(NY-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2)+(NY-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2)+(NY-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (2.*(NX+NY)-8, 2*(NU+NV)-8);
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, tolerance_border, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_border, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_border, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1));
    
    j = NV-1;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1+(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1+(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1+(NU-2)));
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)));

    i = NU-1;
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)+(NV-2)));
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


  cx = gsl_vector_alloc ((NU-2)*(NV-2));
  cy = gsl_vector_alloc ((NU-2)*(NV-2));
  cz = gsl_vector_alloc ((NU-2)*(NV-2));
  cov_x = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_y = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_z = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));

  X = gsl_matrix_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));
  x = gsl_vector_alloc ((NX-2)*(NY-2));
  y = gsl_vector_alloc ((NX-2)*(NY-2));
  z = gsl_vector_alloc ((NX-2)*(NY-2));

  /* Initial guess */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      gint m = i / (NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, (i-1)+ (j-1)*(NU-2), p[m][n].x);
      gsl_vector_set (cy, (i-1)+ (j-1)*(NU-2), p[m][n].y);
      gsl_vector_set (cz, (i-1)+ (j-1)*(NU-2), p[m][n].z);
    }
  
  /* Set problem */
  Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      for ( k = 1; k < NX-1; k++)
  	for ( l = 1; l < NY-1; l++) {
  	  gsl_bspline_eval (u[k][l], Bu, sp->w_u);
  	  gsl_bspline_eval (v[k][l], Bv, sp->w_v);
  	  gsl_matrix_set (X, (k-1) + (NX-2)*(l-1), (i-1) + (j-1)*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    for ( j = 1; j < NY-1; j++) {
      gsl_vector_set (x, (i-1) + (j-1)*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
      gsl_vector_set (y, (i-1) + (j-1)*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
      gsl_vector_set (z, (i-1) + (j-1)*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
    }
  }

  work = gsl_multifit_linear_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));

  /* Multidimensional fitting usind singular values decomposition method */
  gsl_multifit_linear_svd (X, x, tolerance_inside, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_inside, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_inside, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (i-1)+(j-1)*(NU-2)));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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

  spline2d_reinit_panels_physical_quantities (sp);

  return sp;
}

Spline2D * spline2d_fit_geometry5 (FILE * finput, gint M, gint N, gint order, gint inner, gint outer, gint NX, gint NY)
{
  gdouble tolerance_border = 1e-2;
  gdouble tolerance_inside = 1e-6;
  gdouble centripetal_reparam = FALSE;

  Spline2D * sp = spline2d_new (M, N, order, inner, outer);
  spline2d_init_panels (sp);
  coeff_set_var_to_zero (sp, 0);
  coeff_set_var_to_zero (sp, 1);
  coeff_set_var_to_zero (sp, 2);

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  Point p[NX][NY];
  gdouble u[NX][NY], v[NX][NY];
  gint i, j, k, l;

  float px, py, pz;
  i = j = 0;
  while (fscanf(finput, "%f %f %f\n", &px, &py, &pz) != EOF) {
    p[i][j].x = px;
    p[i][j].y = py;
    p[i][j].z = pz;

    i++;
    if ( i == NX ) {
      i = 0;
      j++;
    }
  }

  // Centripetal reparametrization of the whole surface
  if (centripetal_reparam) {
    for ( i = 0; i < NX; i++) {
      gdouble lv = 0.;
      v[i][0] = lv;
      for ( j = 1; j < NY; j++) {
	lv += point_distance (p[i][j], p[i][j-1]);
	v[i][j] = lv;
      }
      for ( j = 1; j < NY; j++) {
	if (lv != 0.)
		v[i][j] /= lv;
	else
	v[i][j] = j/(NY-1.);
      }
    }

    for ( j = 0; j < NY; j++) {
      gdouble lu = 0.;
      u[0][j] = 0.;
      for ( i = 1; i < NX; i++) {
	lu += point_distance (p[i][j], p[i-1][j]);
	u[i][j] = lu;
      }
      for ( i = 1; i < NX; i++) {
	if (lu != 0.)
	  u[i][j] /= lu;
	else
	  u[i][j] = i/(NX-1.);
      }
    }
  }
  else {
    for ( j = 0; j < NY; j++) {
      for ( i = 0; i < NX; i++) {
	u[i][j] = i/(NX-1.);
	v[i][j] = j/(NY-1.);
      }
    }
  }

  
  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    /* for ( j = 0; j < NY; j++) */ {
      j = 0;
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);

  // Set corners
  coeff_assign (sp, 0, 0, 0, p[0][0].x);
  coeff_assign (sp, 0, 0, 1, p[0][0].y);
  coeff_assign (sp, 0, 0, 2, p[0][0].z);

  coeff_assign (sp, NU-1, 0, 0, p[NX-1][0].x);
  coeff_assign (sp, NU-1, 0, 1, p[NX-1][0].y);
  coeff_assign (sp, NU-1, 0, 2, p[NX-1][0].z);

  coeff_assign (sp, NU-1, NV-1, 0, p[NX-1][NY-1].x);
  coeff_assign (sp, NU-1, NV-1, 1, p[NX-1][NY-1].y);
  coeff_assign (sp, NU-1, NV-1, 2, p[NX-1][NY-1].z);

  coeff_assign (sp, 0, NV-1, 0, p[0][NY-1].x);
  coeff_assign (sp, 0, NV-1, 1, p[0][NY-1].y);
  coeff_assign (sp, 0, NV-1, 2, p[0][NY-1].z);

  // Fit borders
  

  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;
  double chisqx, chisqy, chisqz;

  cx = gsl_vector_alloc (2*(NU+NV)-8);
  cy = gsl_vector_alloc (2*(NU+NV)-8);
  cz = gsl_vector_alloc (2*(NU+NV)-8);
  cov_x = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_y = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_z = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);

  X = gsl_matrix_alloc (2*(NX+NY)-8, 2*(NU+NV)-8);
  x = gsl_vector_alloc (2*(NX+NY)-8);
  y = gsl_vector_alloc (2*(NX+NY)-8);
  z = gsl_vector_alloc (2*(NX+NY)-8);

  /* Initial guess */ /* Could be improved using a spline fit at first */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, i-1, p[m][n].x);
    gsl_vector_set (cy, i-1, p[m][n].y);
    gsl_vector_set (cz, i-1, p[m][n].z);

    j = NV-1;
    m = i / (NU -1.)*(NX-1.);
    n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (i-1)+(NU-2), p[m][n].x);
    gsl_vector_set (cy, (i-1)+(NU-2), p[m][n].y);
    gsl_vector_set (cz, (i-1)+(NU-2), p[m][n].z);
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (j-1)+2*(NU-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2), p[m][n].z);

    i = NU-1;
    gsl_vector_set (cx, (j-1)+2*(NU-2)+(NV-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2)+(NV-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2)+(NV-2), p[m][n].z);
  }

  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++) {

    j = 0;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }

    j = NV-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    
    i = NU-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    j = 0;
    gsl_vector_set (x, (i-1), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    j = NY-1;
    gsl_vector_set (x, (i-1) + (NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1) + (NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1) + (NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  for ( j = 1; j < NY-1; j++) {
    i = 0;
    gsl_vector_set (x, (j-1)+2*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    i = NX-1;
    gsl_vector_set (x, (j-1)+2*(NX-2)+(NY-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2)+(NY-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2)+(NY-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (2.*(NX+NY)-8, 2*(NU+NV)-8);
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, tolerance_border, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_border, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_border, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1));
    
    j = NV-1;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1+(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1+(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1+(NU-2)));
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)));

    i = NU-1;
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)+(NV-2)));
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


  cx = gsl_vector_alloc ((NU-2)*(NV-2));
  cy = gsl_vector_alloc ((NU-2)*(NV-2));
  cz = gsl_vector_alloc ((NU-2)*(NV-2));
  cov_x = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_y = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_z = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));

  X = gsl_matrix_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));
  x = gsl_vector_alloc ((NX-2)*(NY-2));
  y = gsl_vector_alloc ((NX-2)*(NY-2));
  z = gsl_vector_alloc ((NX-2)*(NY-2));

  /* Initial guess */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      gint m = i / (NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, (i-1)+ (j-1)*(NU-2), p[m][n].x);
      gsl_vector_set (cy, (i-1)+ (j-1)*(NU-2), p[m][n].y);
      gsl_vector_set (cz, (i-1)+ (j-1)*(NU-2), p[m][n].z);
    }
  
  /* Set problem */
  Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      for ( k = 1; k < NX-1; k++)
  	for ( l = 1; l < NY-1; l++) {
  	  gsl_bspline_eval (u[k][l], Bu, sp->w_u);
  	  gsl_bspline_eval (v[k][l], Bv, sp->w_v);
  	  gsl_matrix_set (X, (k-1) + (NX-2)*(l-1), (i-1) + (j-1)*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    for ( j = 1; j < NY-1; j++) {
      gsl_vector_set (x, (i-1) + (j-1)*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
      gsl_vector_set (y, (i-1) + (j-1)*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
      gsl_vector_set (z, (i-1) + (j-1)*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
    }
  }

  work = gsl_multifit_linear_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));

  /* Multidimensional fitting usind singular values decomposition method */
  gsl_multifit_linear_svd (X, x, tolerance_inside, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_inside, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_inside, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (i-1)+(j-1)*(NU-2)));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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

  spline2d_reinit_panels_physical_quantities (sp);

  return sp;
}

Spline2D * spline2d_fit_geometry8 (FILE * finput, gint M, gint N, gint order, gint inner, gint outer, gint NX, gint NY)
{
  gdouble tolerance_border = 1e-5;
  gdouble tolerance_inside = 1e-6;
  gdouble centripetal_reparam = FALSE;

  Spline2D * sp = spline2d_new (M, N, order, inner, outer);
  spline2d_init_panels (sp);
  coeff_set_var_to_zero (sp, 0);
  coeff_set_var_to_zero (sp, 1);
  coeff_set_var_to_zero (sp, 2);

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);

  g_assert ( NX > NU && NY > NV);

  Point p[NX][NY];
  gdouble u[NX][NY], v[NX][NY];
  gint i, j, k, l;

  float px, py, pz;
  i = j = 0;

  /* GArray * ppp = g_array_new (FALSE, FALSE, sizeof(Point)); */
  /* Point ptmp; */
  /* while (fscanf(finput, "%f %f %f\n", &px, &py, &pz) != EOF) { */
  /*   ptmp.x = px; ptmp.y = py; ptmp.z = pz; */
  /*   g_array_append_val (ppp, ptmp); */
  /* } */

  GPtrArray * curves = g_ptr_array_new ();
  DCurve * dc = dcurve_new (0);
  Point p0;
  fscanf(finput, "%f %f %f\n", &px, &py, &pz);
  p0.x = px; p0.y = py; p0.z = pz;
  g_array_append_val (dc->p, p0);
  while (fscanf(finput, "%f %f %f\n", &px, &py, &pz) != EOF) {
    Point p1 /* = g_array_index (ppp, Point, i) */;
    p1.x = px; p1.y = py; p1.z = pz;
    if (p0.z != p1.z) {
      g_ptr_array_add (curves, dc);
      dc = dcurve_new (0);
    }
    g_array_append_val (dc->p, p1);
    //fprintf (stdout, "%f %f %f \n", px, py, pz);
    p0 = p1;
  }
  g_ptr_array_add (curves, dc);


  gint order_1d = 4;
  gint size_1d = 10;

  gsl_bspline_workspace * w_u = gsl_bspline_alloc (order_1d, size_1d+1);
  gsl_bspline_knots_uniform (0., 1., w_u);
  gdouble coeffs[gsl_bspline_ncoeffs (w_u)];

  gsl_vector * cx1d = gsl_vector_alloc (gsl_bspline_ncoeffs (w_u));
  gsl_vector * cy1d = gsl_vector_alloc (gsl_bspline_ncoeffs (w_u));
  gsl_vector * cz1d = gsl_vector_alloc (gsl_bspline_ncoeffs (w_u));
  gsl_matrix * cov_x1d = gsl_matrix_alloc (gsl_bspline_ncoeffs (w_u),
					   gsl_bspline_ncoeffs (w_u));
  gsl_matrix * cov_y1d = gsl_matrix_alloc (gsl_bspline_ncoeffs (w_u),
					   gsl_bspline_ncoeffs (w_u));
  gsl_matrix * cov_z1d = gsl_matrix_alloc (gsl_bspline_ncoeffs (w_u),
					   gsl_bspline_ncoeffs (w_u));
  double chisqx, chisqy, chisqz;

  // Might want to impose boundary points later on
  for ( i = 0; i < /* curves->len */1; i++) {
    dc = g_ptr_array_index (curves, i);
    gint npoints = dc->p->len;
    gsl_vector * x1d = gsl_vector_alloc (npoints);
    gsl_vector * y1d = gsl_vector_alloc (npoints);
    gsl_vector * z1d = gsl_vector_alloc (npoints);
    gsl_matrix * X1d = gsl_matrix_alloc (npoints, gsl_bspline_ncoeffs (w_u));

    /* Initial guess */
    for ( j = 0; j < size_1d; j++) {
      gint m = j/(size_1d)*npoints;
      p0 = g_array_index (dc->p, Point, m);
      gsl_vector_set (cx1d, j, p0.x);
      gsl_vector_set (cy1d, j, p0.y);
      gsl_vector_set (cz1d, j, p0.z);
    }

    /* Set problem */
    gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (w_u));

    p0 = g_array_index (dc->p, Point, 0);
    p0.xi = 0.;
    gdouble l0 = 0.;
    for ( k = 1; k < npoints; k++) {
      Point p1 = g_array_index (dc->p, Point, k);
      l0 += point_distance (p0, p1);
      p1.xi = l0;
      g_array_index (dc->p, Point, k) = p1;
      p0 = p1;
    }

    for (k = 1; k < npoints; k++) {
      Point p1 = g_array_index (dc->p, Point, k);
      p1.xi /= l0;
      g_array_index (dc->p, Point, k) = p1;
    }

    for ( k = 0; k < npoints; k++) {
      p0 = g_array_index (dc->p, Point, k);
      //p0.xi = k/(npoints-1.);
      fprintf (stderr, "xi %f %f \n", p0.xi, p0.x);
      gsl_bspline_eval (p0.xi, Bu, w_u);
      for ( l = 0; l < gsl_bspline_ncoeffs (w_u); l++) {
	gsl_matrix_set (X1d, k , l, gsl_vector_get(Bu, l));
      }
    }

    /* right hand side */
    for ( k = 0; k < npoints; k++) {
      p0 = g_array_index (dc->p, Point, k);
      gsl_vector_set (x1d, k, p0.x);
      gsl_vector_set (y1d, k, p0.y);
      gsl_vector_set (z1d, k, p0.z);
    }

    gsl_multifit_linear_workspace * work_1d = gsl_multifit_linear_alloc (npoints, gsl_bspline_ncoeffs (w_u));
    size_t rank_x_1d, rank_y_1d, rank_z_1d;

    gsl_multifit_linear_svd (X1d, x1d, tolerance_border, &rank_x_1d, cx1d, cov_x1d, &chisqx, work_1d);
    gsl_multifit_linear_svd (X1d, y1d, tolerance_border, &rank_y_1d, cy1d, cov_y1d, &chisqy, work_1d);
    gsl_multifit_linear_svd (X1d, z1d, tolerance_border, &rank_z_1d, cz1d, cov_z1d, &chisqz, work_1d);

    fprintf (stderr, "chisq: %e %e %e\n", chisqx, chisqy, chisqz);

    gdouble u = 0.;
    for ( u = 0.; u < 1.0001; u += 0.0025) {
      Point ptmp;
      ptmp.x = ptmp.y = ptmp.z = 0.;
      gsl_bspline_eval (u, Bu, w_u);
      for ( k = 0; k < gsl_bspline_ncoeffs (w_u); k++) {
	ptmp.x += gsl_vector_get (cx1d, k)*gsl_vector_get (Bu, k);
	ptmp.y += gsl_vector_get (cy1d, k)*gsl_vector_get (Bu, k);
	ptmp.z += gsl_vector_get (cz1d, k)*gsl_vector_get (Bu, k);
      }
      fprintf (stdout, "%f %f %f \n", ptmp.x, ptmp.y, ptmp.z);
    }

    //fprintf (stdout, "%i ", i);
    //dcurve_print (dc, stdout);

    gsl_vector_free (x1d);
    gsl_matrix_free (X1d);
  }
  gsl_matrix_free (cov_x1d);
  gsl_vector_free (cx1d);
  gsl_bspline_free (w_u);
  g_assert_not_reached ();



  // Centripetal fit of the waterlines using 1d b-splines

  // Generates points + other half

  // 2D fit


  // Centripetal reparametrization of the whole surface
  if (centripetal_reparam) {
    for ( i = 0; i < NX; i++) {
      gdouble lv = 0.;
      v[i][0] = lv;
      for ( j = 1; j < NY; j++) {
	lv += point_distance (p[i][j], p[i][j-1]);
	v[i][j] = lv;
      }
      for ( j = 1; j < NY; j++) {
	if (lv != 0.)
		v[i][j] /= lv;
	else
	v[i][j] = j/(NY-1.);
      }
    }

    for ( j = 0; j < NY; j++) {
      gdouble lu = 0.;
      u[0][j] = 0.;
      for ( i = 1; i < NX; i++) {
	lu += point_distance (p[i][j], p[i-1][j]);
	u[i][j] = lu;
      }
      for ( i = 1; i < NX; i++) {
	if (lu != 0.)
	  u[i][j] /= lu;
	else
	  u[i][j] = i/(NX-1.);
      }
    }
  }
  else {
    for ( j = 0; j < NY; j++) {
      for ( i = 0; i < NX; i++) {
	u[i][j] = i/(NX-1.);
	v[i][j] = j/(NY-1.);
      }
    }
  }

  
  FILE * fp = fopen ("points2fit.dat","w");
  for ( i = 0; i < NX; i++) {
    /* for ( j = 0; j < NY; j++) */ {
      j = 0;
      fprintf (fp, "%f %f %f %i %i\n", p[i][j].x, p[i][j].y, p[i][j].z, i, j);
    }
  }
  fclose (fp);

  // Set corners
  coeff_assign (sp, 0, 0, 0, p[0][0].x);
  coeff_assign (sp, 0, 0, 1, p[0][0].y);
  coeff_assign (sp, 0, 0, 2, p[0][0].z);

  coeff_assign (sp, NU-1, 0, 0, p[NX-1][0].x);
  coeff_assign (sp, NU-1, 0, 1, p[NX-1][0].y);
  coeff_assign (sp, NU-1, 0, 2, p[NX-1][0].z);

  coeff_assign (sp, NU-1, NV-1, 0, p[NX-1][NY-1].x);
  coeff_assign (sp, NU-1, NV-1, 1, p[NX-1][NY-1].y);
  coeff_assign (sp, NU-1, NV-1, 2, p[NX-1][NY-1].z);

  coeff_assign (sp, 0, NV-1, 0, p[0][NY-1].x);
  coeff_assign (sp, 0, NV-1, 1, p[0][NY-1].y);
  coeff_assign (sp, 0, NV-1, 2, p[0][NY-1].z);

  // Fit borders
  

  /* Feed in base spline values to the GSL structure*/
  gsl_vector * cx, * cy, *cz, *x, * y, * z;
  gsl_matrix * cov_x, *cov_y, *cov_z, *X;

  cx = gsl_vector_alloc (2*(NU+NV)-8);
  cy = gsl_vector_alloc (2*(NU+NV)-8);
  cz = gsl_vector_alloc (2*(NU+NV)-8);
  cov_x = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_y = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);
  cov_z = gsl_matrix_alloc (2*(NU+NV)-8, 2*(NU+NV)-8);

  X = gsl_matrix_alloc (2*(NX+NY)-8, 2*(NU+NV)-8);
  x = gsl_vector_alloc (2*(NX+NY)-8);
  y = gsl_vector_alloc (2*(NX+NY)-8);
  z = gsl_vector_alloc (2*(NX+NY)-8);

  /* Initial guess */ /* Could be improved using a spline fit at first */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, i-1, p[m][n].x);
    gsl_vector_set (cy, i-1, p[m][n].y);
    gsl_vector_set (cz, i-1, p[m][n].z);

    j = NV-1;
    m = i / (NU -1.)*(NX-1.);
    n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (i-1)+(NU-2), p[m][n].x);
    gsl_vector_set (cy, (i-1)+(NU-2), p[m][n].y);
    gsl_vector_set (cz, (i-1)+(NU-2), p[m][n].z);
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;
    gint m = i / (NU -1.)*(NX-1.);
    gint n = j / (NV -1.)*(NY-1.);
    gsl_vector_set (cx, (j-1)+2*(NU-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2), p[m][n].z);

    i = NU-1;
    gsl_vector_set (cx, (j-1)+2*(NU-2)+(NV-2), p[m][n].x);
    gsl_vector_set (cy, (j-1)+2*(NU-2)+(NV-2), p[m][n].y);
    gsl_vector_set (cz, (j-1)+2*(NU-2)+(NV-2), p[m][n].z);
  }

  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++) {

    j = 0;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }

    j = NV-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (i-1) + (NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    
    i = NU-1;
    for ( k = 1; k < NX-1; k++) {
      l = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) , (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
      
      l = NY-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (k-1) + (NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
    for ( l = 1; l < NY-1; l++) {
      k = 0;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));

      k = NX-1;
      gsl_bspline_eval (u[k][l], Bu, sp->w_u);
      gsl_bspline_eval (v[k][l], Bv, sp->w_v);
      gsl_matrix_set (X, (l-1) + 2*(NX-2) + (NY-2), (j-1)+2*(NU-2)+(NV-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    j = 0;
    gsl_vector_set (x, (i-1), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    j = NY-1;
    gsl_vector_set (x, (i-1) + (NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (i-1) + (NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (i-1) + (NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  for ( j = 1; j < NY-1; j++) {
    i = 0;
    gsl_vector_set (x, (j-1)+2*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));

    i = NX-1;
    gsl_vector_set (x, (j-1)+2*(NX-2)+(NY-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
    gsl_vector_set (y, (j-1)+2*(NX-2)+(NY-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
    gsl_vector_set (z, (j-1)+2*(NX-2)+(NY-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
  }

  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (2.*(NX+NY)-8, 2*(NU+NV)-8);
  size_t rank_x, rank_y, rank_z;
  gsl_multifit_linear_svd (X, x, tolerance_border, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_border, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_border, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++) {
    j = 0;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1));
    
    j = NV-1;
    
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, i-1+(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, i-1+(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, i-1+(NU-2)));
  }

  for ( j = 1; j < NV-1; j++) {
    i = 0;

    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)));

    i = NU-1;
    coeff_assign (sp, i, j, 0, gsl_vector_get(cx, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 1, gsl_vector_get(cy, j-1+2*(NU-2)+(NV-2)));
    coeff_assign (sp, i, j, 2, gsl_vector_get(cz, j-1+2*(NU-2)+(NV-2)));
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


  cx = gsl_vector_alloc ((NU-2)*(NV-2));
  cy = gsl_vector_alloc ((NU-2)*(NV-2));
  cz = gsl_vector_alloc ((NU-2)*(NV-2));
  cov_x = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_y = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));
  cov_z = gsl_matrix_alloc ((NU-2)*(NV-2), (NU-2)*(NV-2));

  X = gsl_matrix_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));
  x = gsl_vector_alloc ((NX-2)*(NY-2));
  y = gsl_vector_alloc ((NX-2)*(NY-2));
  z = gsl_vector_alloc ((NX-2)*(NY-2));

  /* Initial guess */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      gint m = i / (NU -1.)*(NX-1.);
      gint n = j / (NV -1.)*(NY-1.);
      gsl_vector_set (cx, (i-1)+ (j-1)*(NU-2), p[m][n].x);
      gsl_vector_set (cy, (i-1)+ (j-1)*(NU-2), p[m][n].y);
      gsl_vector_set (cz, (i-1)+ (j-1)*(NU-2), p[m][n].z);
    }
  
  /* Set problem */
  Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_u));
  Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (sp->w_v));
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      for ( k = 1; k < NX-1; k++)
  	for ( l = 1; l < NY-1; l++) {
  	  gsl_bspline_eval (u[k][l], Bu, sp->w_u);
  	  gsl_bspline_eval (v[k][l], Bv, sp->w_v);
  	  gsl_matrix_set (X, (k-1) + (NX-2)*(l-1), (i-1) + (j-1)*(NU-2), gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j));
  	}
    }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  /* right-hand side */
  for ( i = 1; i < NX-1; i++) {
    for ( j = 1; j < NY-1; j++) {
      gsl_vector_set (x, (i-1) + (j-1)*(NX-2), p[i][j].x-spline2d_eval (sp, u[i][j], v[i][j], 0));
      gsl_vector_set (y, (i-1) + (j-1)*(NX-2), p[i][j].y-spline2d_eval (sp, u[i][j], v[i][j], 1));
      gsl_vector_set (z, (i-1) + (j-1)*(NX-2), p[i][j].z-spline2d_eval (sp, u[i][j], v[i][j], 2));
    }
  }

  work = gsl_multifit_linear_alloc ((NX-2)*(NY-2), (NU-2)*(NV-2));

  /* Multidimensional fitting usind singular values decomposition method */
  gsl_multifit_linear_svd (X, x, tolerance_inside, &rank_x, cx, cov_x, &chisqx, work);
  gsl_multifit_linear_svd (X, y, tolerance_inside, &rank_y, cy, cov_y, &chisqy, work);
  gsl_multifit_linear_svd (X, z, tolerance_inside, &rank_z, cz, cov_z, &chisqz, work);

  fprintf(stderr, "Initial b2 fit statistics:\n");
  fprintf(stderr, "Chisqx: %g Chisqy: %g Chisqz: %g\n\n", chisqx, chisqy, chisqz);

  gsl_multifit_linear_free (work);

  /* Copy the solution of the linear system */
  for ( i = 1; i < NU-1; i++)
    for ( j = 1; j < NV-1; j++) {
      coeff_assign (sp, i, j, 0, gsl_vector_get(cx, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 1, gsl_vector_get(cy, (i-1)+(j-1)*(NU-2)));
      coeff_assign (sp, i, j, 2, gsl_vector_get(cz, (i-1)+(j-1)*(NU-2)));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
  	  gdouble x = k/(30.-1.); /* True only is the input points are uniformly spaced */
  	  gdouble y = l/(10.-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (sp, x, y, 0), spline2d_eval (sp, x, y, 1), spline2d_eval (sp, x, y, 2));
  	}
  fclose (fp);
  
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

  spline2d_reinit_panels_physical_quantities (sp);

  return sp;
}

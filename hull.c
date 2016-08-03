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
#include "hull.h"


void motion_update_euler_matrix (Motion * m)
{
  m->t.x = m->x[0];
  m->t.y = m->x[1];
  m->t.z = m->x[2];

  m->euler_m.a[0][0] = cos(m->x[4])*cos(m->x[5]);
  m->euler_m.a[1][0] = -cos(m->x[3])*sin(m->x[5]) + sin(m->x[3])*sin(m->x[4])*cos(m->x[5]);
  m->euler_m.a[2][0] = sin(m->x[3])*sin(m->x[5]) + cos(m->x[3])*sin(m->x[4])*cos(m->x[5]);
  m->euler_m.a[0][1] = cos(m->x[4])*sin(m->x[5]);
  m->euler_m.a[1][1] = cos(m->x[3])*cos(m->x[5]) + sin(m->x[3])*sin(m->x[4])*sin(m->x[5]);
  m->euler_m.a[2][1] = -sin(m->x[3])*cos(m->x[5]) + cos(m->x[3])*sin(m->x[4])*sin(m->x[5]);
  m->euler_m.a[0][2] = -sin(m->x[4]);
  m->euler_m.a[1][2] = sin(m->x[3])*cos(m->x[4]);
  m->euler_m.a[2][2] = cos(m->x[3])*cos(m->x[4]);
}

/**
 * Applies a transformation (rotation, translation) to a Point.
 * Returns the transformed Point.
 * This assumes that the Euler matrix was processsed beforehand.
 **/
Point motion_transform_point (Point p, Motion * m)
{
  Vector v;

  /* v.x = p.x - m->xg.x; */
  /* v.y = p.y - m->xg.y; */
  /* v.z = p.z - m->xg.z; */

  v = vector_rotate (&m->euler_m, v);

  /* p.x = m->xg.x + v.x + m->t.x; */
  /* p.y = m->xg.y + v.y + m->t.y; */
  /* p.z = m->xg.z + v.z + m->t.z; */

  return p;
}

/**
 * Applies a transformation (only rotation) to a Vector.
 * Returns the transformed Vector.
 * This assumes that the Euler matrix was processsed beforehand.
 **/
Vector motion_transform_vector (Vector v, Motion * m)
{
  return vector_rotate (&m->euler_m, v);
}

/**
 * Checks whether both points are identical
 **/
static gboolean same_point (Point p1, Point p2)
{
  if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z)
    return TRUE;
  return FALSE;
}

/**
 * Returns TRUE if both panels are neighbors
 * i.e. if they have two points in commun
 **/
static gboolean is_neighbor (Panel * p1, Panel * p2)
{
  gint i, j;
  gint num = 0;

  g_assert ( p1 != NULL && p2 != NULL);

  for ( i = 0; i < 4; i++) {
    for ( j = 0; j < 4; j++) {
      if (same_point (p1->p[i], p2->p[j]))
	num ++;
    }
    if (num > 1)
      return TRUE;
  }
  return FALSE;
}

static gboolean is_left_or_right_neighbor (Panel * p1, Panel * p2)
{

  if ( same_point (p1->p[0], p2->p[3]) && same_point (p1->p[1], p2->p[2]))
    return TRUE;

  if ( same_point (p2->p[0], p1->p[3]) && same_point (p2->p[1], p1->p[2]))
   return TRUE;

  return FALSE;
}

static gboolean is_top_or_bottom_neighbor (Panel * p1, Panel * p2)
{
  if ( same_point (p1->p[0], p2->p[1]) && same_point (p1->p[3], p2->p[2]))
    return TRUE;

  if ( same_point (p2->p[0], p1->p[1]) && same_point (p2->p[3], p1->p[2]))
   return TRUE;

  return FALSE;
}

static gboolean panel_neighbor_to_last_row (Panel * panel,
					    Patch * patch)
{
  gint j;

  GArray * row = g_ptr_array_index (patch->rows, patch->rows->len-1);
  for ( j = row->len-1; j >= 0; j--) {
    Panel p = g_array_index (row, Panel, j);
    if (is_neighbor (panel, &p))
      return TRUE;
  }

  return FALSE;
}

static gboolean panel_belongs_to_row (Panel * panel,
				      GArray * row)
{
  gint i;

  for ( i = row->len-1; i >= 0; i--) {
    Panel p = g_array_index (row, Panel, i);
    if (is_left_or_right_neighbor (panel, &p))
    //if (is_top_or_bottom_neighbor (panel, &p))
  	return TRUE;
  }

  return FALSE;
}

static gboolean panel_belongs_to_last_row (Panel * panel,
					   Patch * patch)
{
  gint i;
  GArray * row = g_ptr_array_index (patch->rows,
				    patch->rows->len - 1);

  for ( i = row->len-1; i >= 0; i--) {
    Panel p = g_array_index (row, Panel, i);
    if (is_left_or_right_neighbor (panel, &p))
    //  if (is_top_or_bottom_neighbor (panel, &p))
	return TRUE;
  }

  return FALSE;
}

/* Hull methods */
Hull * hull_new ()
{
  Hull * new = g_malloc (sizeof(Hull));
  new->patches = NULL;
  new->wet_patches = NULL;
  new->fh = NULL;

  new->mg = 0.;
  new->xg.x = new->xg.y = new->xg.z = 0.;
  new->U.x = new->U.y = new->U.z = 0.;

  gint i, j;

  for ( i = 0; i < 6; i++)
    for ( j = 0; j < 6; j++) {
      new->M[i][j] = 0.;
      new->A[i][j] = 0.;
      new->D[i][j] = 0.;
      new->R[i][j] = 0.;
    }

  for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++)
      new->Ig[i][j] = 0.;

  // Set initial motion to zero
  for ( i = 0; i < 6; i++)
    new->m.x[i] = new->m.v[i] = new->m.u[i] = 0.;

  for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++)
      new->m.euler_m.a[i][j] = new->m.euler_r.a[i][j] = 0.;

  for ( i = 0; i < 3; i++)
    new->m.euler_m.a[i][i] = new->m.euler_r.a[i][i] = 1.;

  new->m.equilibrium.x = new->m.equilibrium.y = new->m.equilibrium.z = 0.;
  new->m.t.x = new->m.t.y = new->m.t.z = 0.;

  new->mass_lu1 = NULL;
  new->mass_lu2 = NULL;

  return new;
}

void hull_destroy (Hull * h)
{
  g_assert (h != NULL);
  GSList * patches = h->patches;
  while (patches) {
    spline2d_destroy (patches->data);
  }
  g_slist_free (h->patches);
  patches = h->wet_patches;
  while (patches) {
    spline2d_destroy (patches->data);
  }
  g_slist_free (h->patches);
  g_free (h);
}

GPtrArray * sort_panels (GArray * all_panels)
{
  gint i;
  GPtrArray * patches = g_ptr_array_new ();
  Patch * patch = patch_new ();

  g_assert (all_panels != NULL);
  Panel p = g_array_index (all_panels, Panel, 0);
  patch_add_row (patch);
  patch_add_panel (patch, p);
  for ( i = 1; i < all_panels->len; i++) {
    p = g_array_index (all_panels, Panel, i);
    if (panel_belongs_to_last_row (&p, patch)) {
      
      patch_add_panel (patch, p);
    }
    else if (panel_neighbor_to_last_row (&p, patch)) {
      patch_add_row (patch);
      patch_add_panel (patch, p);
    }
    else {
      g_ptr_array_add (patches, patch);
      patch = patch_new ();
      patch_add_row (patch);
      patch_add_panel (patch, p);
    }
  }
  g_ptr_array_add (patches, patch);

  patch_tag_borders (patch);
  
  return patches;
}

/**
 * Reads the hull mesh in fp. The hull is made of 
 * quadrilateral panels and is initially stored in GDF format.
 **/
void hull_read (Hull * h, FILE * fp, gint M, gint N, gboolean flip_u, gboolean flip_v, gboolean centripetal_reparam,
		gboolean swap)
{
  Panel p;
  char buffer[100];
  GArray * all_panels = g_array_new (FALSE, FALSE, sizeof(Panel));
  gint i;

  g_assert (h);
  g_assert (fp);

  /* Skip 4 first lines */
  fgets(buffer, 100, fp);
  fgets(buffer, 100, fp);
  fgets(buffer, 100, fp);
  fgets(buffer, 100, fp);
  
  while (panel_read (&p, fp)) {
    g_array_append_val (all_panels, p);
  }

  /* Sorts the panels into patches */
  GPtrArray * patches = sort_panels (all_panels);

  /* Fit the spline surface to the panels */
  GPtrArray * spline_patches = g_ptr_array_new ();
  for ( i = 0; i < patches->len; i++) {
    Patch * patch = g_ptr_array_index (patches, i);
    
 
    g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, M, N, 3, 4, 3, flip_u, flip_v, centripetal_reparam, swap));
    //g_ptr_array_add (spline_patches, spline2d_interpolate_geometry (patch, M, N, 3, 4, 3, flip_u, flip_v));
  }
  
  /* Replace the read patches by their spline represnetation */
  for ( i = 0; i < patches->len; i++) {
    Patch * patch = g_ptr_array_index (patches, i);
    patch_destroy (patch);
    g_ptr_array_free (patches, TRUE);
  }

  /* Linking of the patches of the hull to hull motion */
  for ( i = 0; i < spline_patches->len; i++) {
    Spline2D * sp = g_ptr_array_index (spline_patches, i);
    spline2d_init_panels (sp);
    g_assert (sp != NULL);
    h->patches = g_slist_append (h->patches, sp);
  }

  g_ptr_array_free (spline_patches, FALSE);

  g_array_free (all_panels, TRUE);
}

void hull_read_old (Hull * h, FILE * fp, gint M, gint N, gboolean flip, gboolean centripetal_reparam)
{
  Panel p;
  char buffer[100];
  GArray * all_panels = g_array_new (FALSE, FALSE, sizeof(Panel));
  gint i;

  g_assert (h);
  g_assert (fp);

  /* Skip 4 first lines */
  fgets(buffer, 100, fp);
  fgets(buffer, 100, fp);
  fgets(buffer, 100, fp);
  fgets(buffer, 100, fp);
  
  while (panel_read (&p, fp)) {
    g_array_append_val (all_panels, p);
  }

  /* Sorts the panels into patches */
  GPtrArray * patches = sort_panels (all_panels);

  /* Fit the spline surface to the panels */
  GPtrArray * spline_patches = g_ptr_array_new ();
  for ( i = 0; i < patches->len; i++) {
    Patch * patch = g_ptr_array_index (patches, i);
    
    //  g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, 10, 80, 3, 4, 3));
    // g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, 10, 20, 3, 6, 2));
    //g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, 12, 30, 3, 4, 3));

    //  g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, 20, 20, 3, 4, 3));

    //g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, 20, 15, 3, 3, 3));

    //g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, 15, 15, 3, 4, 3));
    g_ptr_array_add (spline_patches, spline2d_fit_geometry2_old (patch, M, N, 3, 4, 3, flip));
    //g_ptr_array_add (spline_patches, spline2d_fit_geometry2 (patch, M, N, 3, 4, 3, flip, centripetal_reparam));
  }
  
  /* Replace the read patches by their spline represnetation */
  for ( i = 0; i < patches->len; i++) {
    Patch * patch = g_ptr_array_index (patches, i);
    patch_destroy (patch);
    g_ptr_array_free (patches, TRUE);
  }

  /* Linking of the patches of the hull to hull motion */
  for ( i = 0; i < spline_patches->len; i++) {
    Spline2D * sp = g_ptr_array_index (spline_patches, i);
    spline2d_init_panels (sp);
    g_assert (sp != NULL);
    h->patches = g_slist_append (h->patches, sp);
  }

  g_ptr_array_free (spline_patches, FALSE);

  g_array_free (all_panels, TRUE);
}

gdouble hull_normal_velocity (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Hull * hull = (Hull *) data;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);;
  return 0.;
}

Vector hull_normal_time_derivative (Hull * hull, gdouble u, gdouble v)
{
  Vector nt;
  nt.x = nt.y = nt.z = 0.;

 return nt;
}

/**
 * Print the hull in Gnuplot format into a file called hull.out
 **/
void hull_print (Hull * h, FILE * fout)
{
  FILE * fp = fopen ("hull.out","w");

  g_assert ( h != NULL);

  GSList * patches = h->patches;
  while (patches) {
    fprintf(fp, "#New PAtch \n");
    Spline2D * sp = patches->data;
    g_assert (sp != NULL);
    spline2d_print_panels (patches->data, fp);
    patches = patches->next;
  }

  fclose (fp);
}

void spline2d_print_panels_mayavi (Spline2D * splines, FILE * fp)
{
  g_assert (fp != NULL);
  gint i, j;
  
  Spline2D * sp = splines;
  gint M = 0;

  while (sp) {
    M += sp->M+1;
    sp = sp->next;
  }

  sp = splines;
  fprintf (fp, "%i %i\n", M, sp->N+1, 1);
  while (sp) {
    for ( j = 0; j < sp->N+1; j++) {
      gdouble v = gsl_vector_get (sp->w_v->knots, sp->k + j - 1);
      for ( i = 0; i < sp->M+1; i++) {
	gdouble u = gsl_vector_get (sp->wx_u->knots, sp->k + i - 1);
	Point p = spline2d_eval_point (sp, u, v);
	fprintf (fp, "%f %f %f\n", p.x, p.y, p.z);
      }
    }
    sp = sp->next;
  }
}

void hull_print_mayavi (Hull * h, FILE * fout)
{
  

  g_assert ( h != NULL);

  GSList * patches = h->patches;
  gint num = 0;
  while (patches) {
    num++;
    FILE * fp = fopen (g_strdup_printf ("hull_mayavi_patch_%i.out", num),"w");
    Spline2D * sp = patches->data;
    g_assert (sp != NULL);
    spline2d_print_panels_mayavi (patches->data, fp);
    patches = patches->next;
    fclose (fp);
  }
}

void hull_print_gnuplot (Hull * h, FILE * fout, gint var, gdouble t)
{
  FILE * fp = fopen (g_strdup_printf ("hull-color_%5.4f.out", t),"w");

  g_assert ( h != NULL);

  GSList * patches = h->patches;
  while (patches) {
    fprintf(fp, "#New PAtch \n");
    Spline2D * sp = patches->data;
    g_assert (sp != NULL);
    spline2d_print_panels_gnuplot (patches->data, fp, var);
    patches = patches->next;
  }

  fclose (fp);
}

DCurve * hull_intersect_with_free_surface (Hull * h, HeightCurve hz, gdouble t, gpointer data, gint N)
{
  g_assert (h != NULL);
  GSList * inter = NULL;
  GSList * patches = h->patches;

  while (patches) {
    inter = g_slist_concat (inter, spline2d_intersect_with_free_surface (patches->data, hz, t, data));
    patches = patches->next;
  }

  inter = sort_list_by_angle (inter);

  DCurve * dc = NULL;
  if (inter) {
    Point * last = inter->data;
    
    inter = g_slist_append (inter, last);
    
    /* The code crashes in this routine when the size of the object is to small compared to that 
       of the free surface */
    /* FILE * fp = fopen ("inter.tmp","w"); */
    /* GSList * ltmp = inter; */
    /* while (ltmp) { */
    /*   Point * ptmp = ltmp->data; */
    /*   fprintf (fp, "%f %f %f \n", ptmp->x, ptmp->y, ptmp->z); */
    /*   ltmp = ltmp->next; */
    /* } */
    /* fclose (fp); */

    dc = resample_point_list_hybrid (inter, N);
    //dc =  resample_point_list_centripetal (inter, N);
    
    FILE * fp = fopen ("inter.tmp","w");
    dcurve_print (dc, fp);
    fclose (fp);
    
    g_slist_free (inter);
  }
  else {
    dc = dcurve_new (N);
    gint i;
    for ( i = 0; i < N/2; i++) {
      Point p;
      p.x = -1. + 2.*i/(N/2);
      p.y = 0.;
      p.z = hz (p.x, p.y, t, data);
      p.xi = i/(N-1.);
      g_array_append_val (dc->p, p);
    }
    
    gint j;
    for ( j = i; j < N; j++) {
      Point p;
      p.x = 1. - 2.*(j-i)/(N-1.-i);
      p.y = 0.;
      p.z = hz (p.x, p.y, t, data);
      p.xi = j/(N-1.);
      g_array_append_val (dc->p, p);
    }
  }

  return dc;
}

Point spline2d_hull_eval_point (Spline2D * sp, Hull * h, gdouble u, gdouble v)
{
  g_assert_not_reached ();
  // probably obsolete see hull_transformed_eval_point
  return transform_point (spline2d_eval_point (sp, u, v),
			  &h->xg, &h->m.t, &h->m.euler_m);
}

Point spline2d_hull_eval_point_gauss_point (Spline2D * sp,
					    GaussPoints * gp,
					    Hull * h,
					    gint m, gint n)
{
  gint ng = sp->nouter;
  return transform_point (g_array_index (gp->Pi, Point, m + n*ng),
			  &h->xg, &h->m.t, &h->m.euler_m);
}

Point hull_transformed_point (Hull * h, Point p)
{
  return transform_point (p, &h->xg, &h->m.t, &h->m.euler_m);
}

Vector hull_velocity_at_point (Hull * h, Point p)
{
  Motion m = h->m;
  Vector vg, rg, pg;

  vg.x = m.u[0];
  vg.y = m.u[1];
  vg.z = m.u[2];

  rg.x = m.u[3];
  rg.y = m.u[4];
  rg.z = m.u[5];

  pg.x = h->xg.x - p.x;
  pg.y = h->xg.y - p.y;
  pg.z = h->xg.z - p.z;

  // Wisdom for this kind of things comes from Skrew theory (Torseurs)
  // V_P = V_G + PG x  R

  return vector_sum (vg, vector_vector_product (&pg, &rg));
}

Vector hull_velocity_at_point_1 (Hull * h, Point p)
{
  Motion m = h->m;
  Vector vg, rg, pg;

  vg.x = m.u1[0];
  vg.y = m.u1[1];
  vg.z = m.u1[2];

  rg.x = m.u1[3];
  rg.y = m.u1[4];
  rg.z = m.u1[5];

  pg.x = h->xg.x - p.x;
  pg.y = h->xg.y - p.y;
  pg.z = h->xg.z - p.z;

  // Wisdom for this kind of things comes from Skrew theory (Torseurs)
  // V_P = V_G + PG x  R

  return vector_sum (vg, vector_vector_product (&pg, &rg));
}

Point hull_transformed_eval_point (Hull * h, Spline2D * sp,
				   gdouble u, gdouble v)
{
  return hull_transformed_point (h, spline2d_eval_point (sp, u, v));
}

/* Point hull_transformed_eval_point_gauss_point (Hull * h, */
/* 					       Spline2D * sp, */
/* 					       GaussPoints * gp, */
/* 					       gint m, gint n) */
/* { */
/*   gint ng = sp->nouter; */
  
/*   return transform_point (g_array_index (gp->Pi, Point, m + n*ng), */
/* 			  &h->xg, &h->m.t, &h->m.euler_m); */
/* } */

Vector hull_transformed_vector (Hull * h, Vector v)
{
  // This need to be tested properly first
  return transform_vector (v, &h->m.euler_m);
}


/************** Wetted hull ***************************************/


static GSList * intersection (Spline2D * patch, HeightCurve hz, gdouble il,
			      gdouble ih, gdouble jl, gdouble jh, gdouble t, gpointer data)
{
  gdouble tolerance = 1e-3;
  GSList * l = NULL;
  
  /* Check for intersection */
  gdouble dz[4];
  Point c[4];
  gint i;
  c[0] = spline2d_eval_point (patch, il, jl);
  c[1] = spline2d_eval_point (patch, ih, jl);
  c[2] = spline2d_eval_point (patch, ih, jh);
  c[3] = spline2d_eval_point (patch, il, jh);
  
  for ( i = 0; i < 4; i++)
    dz[i] = c[i].z - hz (c[i].x, c[i].y, t, data);
  
  if ( (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.) ||
       (dz[0] < 0. && dz[1] < 0. && dz[2] < 0. && dz[3] < 0.) )
    return l;
  
  /* We refine */
  gdouble ihalf = (il + ih)/2.;
  gdouble jhalf = (jl + jh)/2.;
  if (fabs(il-ih) < tolerance || fabs (jl-jh) < tolerance) {
    Point * p = g_malloc (sizeof(Point));
    *p = spline2d_eval_point (patch, ihalf, jhalf);
    l = g_slist_append (l, p);
    return l;
  }
  
  l = g_slist_concat (l, intersection (patch, hz, il, ihalf, jl, jhalf, t, data));
  l = g_slist_concat (l, intersection (patch, hz, ihalf, ih, jl, jhalf, t, data));
  l = g_slist_concat (l, intersection (patch, hz, il, ihalf, jhalf, jh, t, data));
  l = g_slist_concat (l, intersection (patch, hz, ihalf, ih, jhalf, jh, t, data));

  return l; // ??? Still works when not there
}

GSList * spline2d_intersect_with_free_surface (Spline2D * patch, HeightCurve hz, gdouble t, gpointer data)
{
  GSList * l = NULL;
  gint i, j;

  for ( i = 0; i < patch->M; i++) {
    for ( j = 0; j < patch->N; j++) {
      SPPanel * panel = g_ptr_array_index (patch->panels, i + j*patch->M);
      l = g_slist_concat (l, intersection (patch, hz, panel->u0, panel->u1, panel->v0, panel->v1, t, data));
    }
  }
  return l;
}

GSList * spline2d_intersection_u (Spline2D * sp,
				  gdouble u0, gdouble u1,
				  gdouble v,
				  Hull * h,
				  HeightCurve hz,
				  gdouble t, gpointer data)
{
  GSList * l = NULL;
  gdouble tolerance = 1e-2;
  gdouble dz[2];
  Point c[2];
  gint i;
  /* c[0] = spline2d_eval_point (sp, u0, v); */
  /* c[1] = spline2d_eval_point (sp, u1, v); */
  //c[0] = spline2d_hull_eval_point (sp, h, u0, v);
  c[0] = hull_transformed_eval_point (h, sp, u0, v);
  //c[1] = spline2d_hull_eval_point (sp, h, u1, v);
  c[1] = hull_transformed_eval_point (h, sp, u1, v);
  
  for ( i = 0; i < 2; i++)
    dz[i] = c[i].z - hz (c[i].x, c[i].y, t, data);
  
  if ( (dz[0] > 0. && dz[1] > 0.) || (dz[0] < 0. && dz[1] < 0.) )
    return l;
  
  /* We refine */
  gdouble uhalf = 0.49*u0 + 0.51*u1; // To avoid problems with symmetrical objects
  if (fabs(u1-u0) < tolerance) {
    Point * p = g_malloc (sizeof(Point));
    *p = spline2d_eval_point (sp, uhalf, v);
    p->xi = uhalf;
    l = g_slist_append (l, p);
    return l;
  }

  l = g_slist_concat (l, spline2d_intersection_u (sp, u0, uhalf, v,
						  h, hz, t, data));
  l = g_slist_concat (l, spline2d_intersection_u (sp, uhalf, u1, v,
						  h, hz, t, data));

  return l;
}

GSList * spline2d_intersection_v (Spline2D * sp, gdouble u,
				  gdouble v0, gdouble v1,
				  Hull * h, HeightCurve hz,
				  gdouble t, gpointer data)
{
  GSList * l = NULL;
  gdouble tolerance = 1e-8;
  gdouble dz[2];
  Point c[2];
  gint i;
  /* c[0] = spline2d_eval_point (sp, u, v0); */
  /* c[1] = spline2d_eval_point (sp, u, v1); */
  //c[0] = spline2d_hull_eval_point (sp, h, u, v0);
  c[0] = hull_transformed_eval_point (h, sp, u, v0);
  //c[1] = spline2d_hull_eval_point (sp, h, u, v1);
  c[1] = hull_transformed_eval_point (h, sp, u, v1);
  
  for ( i = 0; i < 2; i++)
    dz[i] = c[i].z - hz (c[i].x, c[i].y, t, data);
  
  if ( (dz[0] > 0. && dz[1] > 0.) || (dz[0] < 0. && dz[1] < 0.) )
    return l;
  
  /* We refine */
  gdouble vhalf = 0.49*v0 + 0.51*v1; // To avoid problem with symmetrical objects
  if (fabs(v1-v0) < tolerance) {
    Point * p = g_malloc (sizeof(Point));
    *p = spline2d_eval_point (sp, u, vhalf);
    p->xi = vhalf;
    l = g_slist_append (l, p);
    return l;
  }

  l = g_slist_concat (l, spline2d_intersection_v (sp, u, v0, vhalf,
						  h, hz, t, data));
  l = g_slist_concat (l, spline2d_intersection_v (sp, u, vhalf, v1,
						  h, hz, t, data));

  return l;
}

/* gboolean sppanel_is_wet (Spline2D * sp, */
/* 			 gdouble u0, gdouble u1, */
/* 			 gdouble v0, gdouble v1, */
/* 			 HeightCurve hz, */
/* 			 gdouble t, gpointer data, */
/* 			 gboolean * intersect) */
/* { */
/*   /\* Check for intersection *\/ */
/*   gdouble dz[4]; */
/*   Point c[4]; */
/*   gint i; */
/*   c[0] = spline2d_eval_point (sp, u0, v0); */
/*   c[1] = spline2d_eval_point (sp, u1, v0); */
/*   c[2] = spline2d_eval_point (sp, u1, v1); */
/*   c[3] = spline2d_eval_point (sp, u0, v1); */
  
/*   for ( i = 0; i < 4; i++) */
/*     dz[i] = c[i].z - hz (c[i].x, c[i].y, t, data); */
  
/*   *intersect = FALSE; */

/*   if ( ! (dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.) ) { */
/*     if ( !(dz[0] < 0. && dz[1] < 0. && dz[2] < 0. && dz[3] < 0.)) { */
/*       * intersect = TRUE; */
/*       /\* fprintf (stderr, "Intersect upstream\n"); *\/ */
/*     } */
/*     return TRUE; */
/*   } */

/*   return FALSE; */
/* } */

static gboolean point_is_wet (Spline2D * sp, gdouble u, gdouble v,
			      Hull * h, HeightCurve hz,
			      gdouble t, gpointer data)
{
  gdouble dz;
  //Point p = spline2d_eval_point (sp, u, v);
  //Point p = spline2d_hull_eval_point (sp, h, u, v);
  Point p = hull_transformed_eval_point (h, sp, u, v);

  dz = p.z - hz (p.x, p.y, t, data);

  if ( dz < 0. )
    return TRUE;

  return FALSE;
}

gboolean sppanel_intersects_curve (SPPanel * spp, HeightCurve hz, gdouble t, gpointer data)
{
  Spline2D * patch = spp->sp;

  /* Check for intersection */
  gdouble dz[4];
  Point c[4];
  gint i;
  c[0] = spline2d_eval_point (patch, spp->u0, spp->v0);
  c[1] = spline2d_eval_point (patch, spp->u1, spp->v0);
  c[2] = spline2d_eval_point (patch, spp->u1, spp->v1);
  c[3] = spline2d_eval_point (patch, spp->u0, spp->v1);
  
  for ( i = 0; i < 4; i++)
    dz[i] = c[i].z - hz (c[i].x, c[i].y, t, data);
  
  if ( !(dz[0] > 0. && dz[1] > 0. && dz[2] > 0. && dz[3] > 0.) &&
	 !(dz[0] < 0. && dz[1] < 0. && dz[2] < 0. && dz[3] < 0.))
    return TRUE;

  return FALSE;
}

Spline2D * spline2d_wet_part (Spline2D * sp, Hull * hull,
			      HeightCurve hz, gdouble t, gpointer data)
{
  gint NU = 3*sp->NU;
  gint NV = 3*sp->NV;

  if (sp->periodic)
    g_warning ("Periodic patch will be turned into standard patch\n");

  Point pp;
  gint i, j, m, n;
  GArray * points = g_array_new (FALSE, FALSE, sizeof(Point));
  gdouble sumdxi = 0.;
  for ( i = 0; i < NU; i++ ) {   
    gdouble u = i/(NU-1.);
    GSList * lp = NULL;

    /* Puts the 0 end in the point list if it is under the surface hz */
    if ( point_is_wet (sp, u, 0., hull, hz, t, data) ) {
      Point * pstart = g_malloc (sizeof(Point));
      *pstart = spline2d_eval_point (sp, u, 0.);
      pstart->xi = 0.;
      lp = g_slist_append (lp, pstart);
    }

    /* Finds values of v for a constant value of u which the patch intersect */
    /* the surface hz */
    lp = g_slist_concat (lp, spline2d_intersection_v (sp, u, 0, 1,
						      hull, hz, t, data));

    /* Puts the 1. end in the point list if it is under the surface hz */
    if ( point_is_wet (sp, u, 1., hull, hz, t, data) ) {
      Point * pstart  = g_malloc (sizeof(Point));
      *pstart = spline2d_eval_point (sp, u, 1.);
      pstart->xi = 1;
      lp = g_slist_append (lp, pstart);
    }

    g_assert ( g_slist_length (lp) == 2 );

    Point * p0 = lp->data;
    Point * p1 = lp->next->data;

    for ( n = 0; n <  NV; n++ ) {
      gdouble v = p0->xi + (gdouble) n / (gdouble) (NV-1.)*(p1->xi-p0->xi);
      pp = spline2d_eval_point (sp, u, v);
      g_array_append_val (points, pp);
    }
    sumdxi += fabs(p1->xi-p0->xi);

    g_free (p0);
    g_free (p1);
    g_slist_free (lp);
  }

  sumdxi /= NU;

  /* if ( sumdxi == 1.) */
  /*   return sp; */

  // Estimates the number of panels in the v direction
  gint new_N = (sumdxi*sp->N)/* +1 */;

  if ( (sp->N*sumdxi - new_N) > 0.2)
    new_N += 1;

  fprintf (stderr, "%i %i %f \n", sp->N, new_N, sumdxi);

  Spline2D * wet_hull = spline2d_fit_geometry6 (points, sp->M, new_N,
						sp->k, sp->ninner,
						sp->nouter, NU, NV);

  NU = gsl_bspline_ncoeffs (wet_hull->w_u);
  NV = gsl_bspline_ncoeffs (wet_hull->w_v);
  /* for ( i = 0; i < NU; i++ ) { */
  /*   for ( j = 0; j < NV; j++ ) { */
  /*     fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /*   } */
  /* } */

  /* j = 0; */
  /* for ( i = 0; i < NU; i++ ) { */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /*   // fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /* } */
  /* //fprintf (stdout, "\n"); */

  /* j = NV-1; */
  /* for ( i = 0; i < NU; i++ ) { */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /*   // fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /* } */
  //fprintf (stdout, "\n");

  /* j = 1; */
  /* for ( i = 0; i < NU; i++ ) { */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /*   // fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /* } */
  /* //fprintf (stdout, "\n"); */

  /* j = NV-2; */
  /* for ( i = 0; i < NU; i++ ) { */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /*   // fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /* } */
  
  /* j = 2; */
  /* for ( i = 0; i < NU; i++ ) { */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /*   // fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /* } */
  /* //fprintf (stdout, "\n"); */

  /* j = NV-3; */
  /* for ( i = 0; i < NU; i++ ) { */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /*   // fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /* } */

  /* i = 0; */
  /* for ( j = 0; j < NV; j++ ) { */
  /*   //fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /* } */
  //fprintf (stdout, "\n");

  /* i = NU-1; */
  /* for ( j = 0; j < NV; j++ ) { */
  /*   //fprintf (stdout, "%i %i %f \n", i, j, coeff (wet_hull, i, j, 1)); */
  /*   coeff_assign (sp, i, j, 1, 0.); */
  /* } */
  /* fprintf (stdout, "\n"); */
  /* g_assert_not_reached (); */

  if ( sumdxi == 1.)
    wet_hull->fully_submerged = TRUE;
					  
  FILE * ffm = fopen ("wet.tmp","w");
  spline2d_print_panels (wet_hull, ffm);
  fclose (ffm);

  g_array_free (points, TRUE);

  return wet_hull;
}

void hull_generate_wet_hull (Hull * hull, HeightCurve hz,
			     gdouble t, gpointer data)
{
  GSList * patches = hull->patches;

  while (patches) {
    Spline2D * sp = patches->data;
    while (sp) {
      hull->wet_patches =
	g_slist_append (hull->wet_patches,
			spline2d_wet_part (sp, hull, hz, t, data));
      sp = sp->next;
    }
    patches = patches->next;
  }
}

void hull_print_wet (Hull * h, FILE * fout)
{

  g_assert ( h != NULL);

  GSList * patches = h->wet_patches;
  
  while (patches) {
    fprintf(fout, "#New PAtch \n");
    Spline2D * sp = patches->data;
    g_assert (sp != NULL);
    spline2d_print_panels (patches->data, fout);
    patches = patches->next;
  }


}

/*****************************************************/
/*                 CLEAN FORMULATION                 */
/*****************************************************/

/* Returns the first order normal from the zero order one */
Vector hull_normal_1 (Hull * hull, Vector * n0)
{
  Vector n1;

  n1.x = -hull->m.x1[5]*n0->y + hull->m.x1[4]*n0->z;
  n1.y = hull->m.x1[5]*n0->x - hull->m.x1[3]*n0->z;
  n1.z = -hull->m.x1[4]*n0->x + hull->m.x1[3]*n0->y;

  n1.x = /* hull->m.Rbi1.a[0][0]*n0->x */
    /* +  */hull->m.Rbi1.a[0][1]*n0->y
    + hull->m.Rbi1.a[0][2]*n0->z;
  n1.y = hull->m.Rbi1.a[1][0]*n0->x
    /* + hull->m.Rbi1.a[1][1]*n0->y */
    + hull->m.Rbi1.a[1][2]*n0->z;
  n1.z = hull->m.Rbi1.a[2][0]*n0->x
    + hull->m.Rbi1.a[2][1]*n0->y
    /* + hull->m.Rbi1.a[2][2]*n0->z */;

  return n1;
}

Vector hull_normal_2 (Hull * hull, Vector * n0)
{
  Vector n2;

  n2.x = hull->m.Rbi2.a[0][0]*n0->x
    +hull->m.Rbi2.a[0][1]*n0->y
    + hull->m.Rbi2.a[0][2]*n0->z;
  n2.y = /* hull->m.Rbi2.a[1][0]*n0->x */
    /* + */ hull->m.Rbi2.a[1][1]*n0->y
    + hull->m.Rbi2.a[1][2]*n0->z;
  n2.z = /* hull->m.Rbi2.a[2][0]*n0->x */
    /* + hull->m.Rbi2.a[2][1]*n0->y */
    /* + */ hull->m.Rbi2.a[2][2]*n0->z;

  return n2;
}

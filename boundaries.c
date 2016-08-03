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
#include "boundaries.h"

static Vector tangent (Point p0, Point p1)
{
  Vector ti;

  ti.x = (p1.x -p0.x)/(p1.xi-p0.xi);
  ti.y = (p1.y -p0.y)/(p1.xi-p0.xi);
  ti.z = (p1.z -p0.z)/(p1.xi-p0.xi);
  return ti;
}

DCurve * hybrid_curve_point_distribution (gint m, BoundaryCurve curve, gdouble time, gpointer * data)
{
  gint N = 10*m;
  /* gdouble lambda_s = 0.2, lambda_k = 0.8; */
  gdouble lambda_s = 1, lambda_k = 0.;
  DCurve * dc = dcurve_new (N);
  DCurve * dc_new = dcurve_new (m);
  /* arc length */
  GArray * s = g_array_new (FALSE, TRUE, sizeof(gdouble));
  Point p0, p1;
  gint i;

  /* High resolution curve sampling */
  for ( i = 0; i < N; i++) {
    p0 = curve (i/(N-1.), time, data);
    g_array_append_val (dc->p, p0);
  }

  /* STEP 1: Init xi to 0 */
  for ( i = 0; i < N; i++) {
    gdouble x0 = 0.;
    g_array_append_val (s, x0);
    p0 = g_array_index (dc->p, Point, i);
    p0.xi = 0.;
    g_array_index (dc->p, Point, i) = p0;
  }

  /* STEP 2: Compute arc length , rescale so that max arc length is 1 */
  /* add to xi weighted by lambda_s */
  p0 = g_array_index (dc->p, Point, 0);
  for ( i = 1; i < N; i++) {
    p1 = g_array_index (dc->p, Point, i);
    gdouble si = g_array_index (s, gdouble, i-1)
      + sqrt((p1.x - p0.x)*(p1.x - p0.x) + (p1.y - p0.y)*(p1.y - p0.y)
  	     + (p1.z - p0.z)*(p1.z - p0.z));
    g_array_index (s, gdouble, i) = si;
    p0 = p1;
  }

  for ( i = 1; i < N; i++) {
    g_array_index (s, gdouble, i) = (m-1.)*g_array_index (s, gdouble, i)
      /g_array_index (s, gdouble, N-1);
    p0 = g_array_index (dc->p, Point, i);
    p0.xi += lambda_s * g_array_index (s, gdouble, i);
    g_array_index (dc->p, Point, i) = p0;
  }

  /* STEP 3: Compute curvature arc length on fine grid. Check if curve has non-    */
  /* trivial amount of curvature. If so, normalize to m, and add into xi, weighted */
  /* by lambda_k. Otherwise, use arc length instead                                */
  GArray * t = g_array_new (FALSE, TRUE, sizeof(Vector));
  GArray * a = g_array_new (FALSE, TRUE, sizeof(gdouble));
  Vector ti;
  ti = tangent (g_array_index (dc->p, Point, 0),
  		g_array_index (dc->p, Point, 1));
  g_array_append_val (t, ti);
  for ( i = 1; i < N-1; i++) {
    ti = tangent (g_array_index (dc->p, Point, i-1),
  		  g_array_index (dc->p, Point, i+1));
    g_array_append_val (t, ti);
  }
  ti = tangent (g_array_index (dc->p, Point, N-2),
  		g_array_index (dc->p, Point, N-1));
  g_array_append_val (t, ti);
  
  for ( i = 0; i < N; i++) {
    g_assert (t->len >= N);
    Vector vtmp = vector_normalise ((Vector) g_array_index (t, Vector, i));
    g_array_index (t, Vector, i) = vtmp;
  }
 
  gdouble ai = 0.;
  Vector t0, t1;
  t0 = g_array_index (t, Vector, 0);
  g_array_append_val (a, ai);
  for ( i = 1; i < N; i++) {
    t1 = g_array_index (t, Vector, i);
    ai = g_array_index (a, gdouble, i-1) + sqrt((t1.x-t0.x)*(t1.x-t0.x) +
  						(t1.y-t0.y)*(t1.y-t0.y) +
  						(t1.z-t0.z)*(t1.z-t0.z));
    g_array_append_val (a, ai);
    t0 = t1;
  }

  if (g_array_index (a, gdouble, N-1) > 0.01) {
    for ( i = 1; i < N; i++) {
      g_array_index (a, gdouble, i) = (m-1.)*g_array_index (a, gdouble, i)
  	/g_array_index (a, gdouble, N-1);
      p0 = g_array_index (dc->p, Point, i);
      p0.xi += lambda_k * g_array_index (a, gdouble, i);
      g_array_index (dc->p, Point, i) = p0;
    }
  }
  else {
    for ( i = 1; i < N; i++) {
      p0 = g_array_index (dc->p, Point, i);
      p0.xi += lambda_k * g_array_index (s, gdouble, i);
      g_array_index (dc->p, Point, i) = p0;
    }
  }

  /* Skip STEP 4: Attractor Points */
  

  /* STEP 5: Obtain point distribution by inverting grid function */
  
  
  p0 = g_array_index (dc->p, Point, N-1);
  p0.xi = m-1.;
  g_array_index (dc->p, Point, N-1) = p0;
  GArray * xi = g_array_sized_new (FALSE, TRUE, sizeof(gdouble), m);
  gdouble x0 = 0.;
  for ( i = 0; i < m; i++)
    g_array_append_val (xi, x0);


  gint j = 1;
  p0 = p1 = g_array_index (dc->p, Point, 0);
  for ( i = 1; i < N; i++) {
    p1 = g_array_index (dc->p, Point, i);
    while ( j <=  p1.xi) {
      g_array_index (xi, gdouble, j) = i/(N-1.) - 1./(N-1.)*(p1.xi-j)/(p1.xi-p0.xi);
      j++;
    }
    p0 = p1;
  }

  for ( i = 0; i < xi->len; i++) {
    p0 = curve (g_array_index (xi, gdouble, i), time, data);
    g_array_append_val (dc_new->p, p0);
  }

  g_array_free (s, TRUE);
  g_array_free (a, TRUE);
  g_array_free (xi, TRUE);
  dcurve_destroy (dc);
  return dc_new;
}

Boundaries * boundaries_new ()
{
  Boundaries * new = g_malloc (sizeof(Boundaries));
  
  new->dcb = NULL;
  new->dct = NULL;
  new->dcl = NULL;
  new->dcr = NULL;

  new->curve_bottom = NULL;
  new->curve_top = NULL;
  new->curve_right = NULL;
  new->curve_left = NULL;

  return new;
}

void boundaries_init (Boundaries * b, gdouble t, gpointer data, gint N, gint M)
{
  g_assert (b != NULL);
  gpointer datum[2];
  datum[0] = b;
  datum[1] = data;

  if (b->curve_bottom != NULL)
    b->dcb = hybrid_curve_point_distribution (N, b->curve_bottom, t, datum);
  if (b->dct == NULL)
    b->dct = hybrid_curve_point_distribution (N, b->curve_top, t, datum);
  if (b->dcl == NULL)
    b->dcl = hybrid_curve_point_distribution (M, b->curve_left, t, datum);
  if (b->dcr == NULL)
    b->dcr = hybrid_curve_point_distribution (M, b->curve_right, t, datum);
}

void boundaries_print (Boundaries * b, FILE * fp)
{
  gint i;

  for (i = 0; i < b->dcb->p->len; i++) {
    Point p0 = g_array_index (b->dcb->p, Point, i);
    fprintf (fp, "%g %g %g\n", p0.x, p0.y, p0.z);
    fprintf (fp, "\n");
  }

  /* for (i = 0; i < b->dcr->p->len; i++) { */
  /*   Point p0 = g_array_index (b->dcr->p, Point, i); */
  /*   fprintf (fp, "%g %g %g\n", p0.x, p0.y, p0.z); */
  /*   fprintf (fp, "\n"); */
  /* } */

  fprintf (fp, "#Top \n");
  for (i = 0; i < b->dct->p->len; i++) {
    Point p0 = g_array_index (b->dct->p, Point, i);
    fprintf (fp, "%g %g %g \n", p0.x, p0.y, p0.z);
    fprintf (fp, "\n");
  }

  /* for (i = 0; i < b->dcl->p->len; i++) { */
  /*   Point p0 = g_array_index (b->dcl->p, Point, i); */
  /*   fprintf (fp, "%g %g %g \n", p0.x, p0.y, p0.z); */
  /*   fprintf (fp, "\n"); */
  /* } */
}


void boundaries_destroy (Boundaries * b)
{
  g_assert (b != NULL);

  dcurve_destroy (b->dcb);
  dcurve_destroy (b->dct);
  dcurve_destroy (b->dcl);
  dcurve_destroy (b->dcr);

  g_free (b);
}

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


gint sign (gdouble x)
{
  return (x > 0) - (x < 0);
}

/* Point methods */

void point_print (Point p, FILE * fp, Transformation * t)
{
  if (t == NULL)
    fprintf(fp, "%g %g %g\n", p.x, p.y ,p.z);
  else {
    Point p2 = transform_point (p, &t->xg, &t->t, &t->euler_m);
    fprintf(fp, "%g %g %g\n", p2.x, p2.y ,p2.z);
  }
}

gdouble point_distance (Point p1, Point p2)
{
  return sqrt((p2.x-p1.x)*(p2.x-p1.x) +
  	      (p2.y-p1.y)*(p2.y-p1.y) +
  	      (p2.z-p1.z)*(p2.z-p1.z));
  /* Vector diff; */
  /* diff.x = p2.x-p1.x; */
  /* diff.y = p2.y-p1.y; */
  /* diff.z = p2.z-p1.z; */
  /* return vector_norm (diff); */
}

gdouble point_distance_squared (Point p1, Point p2)
{
  return ((p2.x-p1.x)*(p2.x-p1.x) +
	  (p2.y-p1.y)*(p2.y-p1.y) +
	  (p2.z-p1.z)*(p2.z-p1.z));
}

/* Vector methods */

gdouble vector_norm (Vector v)
{
  if ( v.x == 0. && v.y == 0. && v.z == 0. )
    return 0.;

  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vector vector_normalise (Vector v)
{
  gdouble n = vector_norm (v);
  Vector v0 = v;
  v0.x /= n;
  v0.y /= n;
  v0.z /= n;
  return v0;
}

gdouble vector_scalar_product (Vector * v0, Vector * v1)
{
  return (v0->x*v1->x+v0->y*v1->y+v0->z*v1->z);
}

gdouble vector_normalized_scalar_product (Vector * v0, Vector * v1)
{
  gdouble l0 = sqrt(v0->x*v0->x+v0->y*v0->y+v0->z*v0->z);
  gdouble l1 = sqrt(v1->x*v1->x+v1->y*v1->y+v1->z*v1->z);

  return (v0->x*v1->x+v0->y*v1->y+v0->z*v1->z)/(l0*l1);
}

Vector vector_vector_product (Vector * v0, Vector * v1)
{
  Vector v2;

  v2.x = v0->y*v1->z-v0->z*v1->y;
  v2.y = v0->z*v1->x-v0->x*v1->z;
  v2.z = v0->x*v1->y-v0->y*v1->x;

  return v2;
}

Vector vector_sum (Vector v0, Vector v1)
{
  Vector v2;

  v2.x = v0.x + v1.x;
  v2.y = v0.y + v1.y;
  v2.z = v0.z + v1.z;

  return v2;
}

Vector vector_times_constant (Vector v, gdouble c)
{
  Vector v0;

  v0.x = c*v.x;
  v0.y = c*v.y;
  v0.z = c*v.z;

  return v0;
}

void vector_set_to_zero (Vector * v)
{
  v->x = 0.;
  v->y = 0.;
  v->z = 0.;
}

/* Matrix3 Methods */

Vector vector_rotate (Matrix3 * m, Vector v0)
{
  Vector v1;

  v1.x = m->a[0][0]*v0.x + m->a[0][1]*v0.y + m->a[0][2]*v0.z;
  v1.y = m->a[1][0]*v0.x + m->a[1][1]*v0.y + m->a[1][2]*v0.z;
  v1.z = m->a[2][0]*v0.x + m->a[2][1]*v0.y + m->a[2][2]*v0.z;

  return v1;
}



void matrix3_inverse (Matrix3 * m)
{
  gdouble a = m->a[0][0];
  gdouble b = m->a[0][1];
  gdouble c = m->a[0][2];
  gdouble d = m->a[1][0];
  gdouble e = m->a[1][1];
  gdouble f = m->a[1][2];
  gdouble g = m->a[2][0];
  gdouble h = m->a[2][1];
  gdouble k = m->a[2][2];

  gdouble det = a*(e*k-f*h)+b*(f*g-k*d)+c*(d*h-e*g);

  g_assert (det != 0.);

  det = 1./det;

  m->a[0][0] = (e*k-f*h)*det;
  m->a[0][1] = (f*g-d*k)*det;
  m->a[0][2] = (d*h-e*g)*det;
  m->a[1][0] = (c*h-b*k)*det;
  m->a[1][1] = (a*k-c*g)*det;
  m->a[1][2] = (g*b-a*h)*det;
  m->a[2][0] = (b*f-c*e)*det;
  m->a[2][1] = (c*d-a*f)*det;
  m->a[2][2] = (a*e-b*d)*det;
}

/* Vector6 method */

void vector6_set_to_zero (Vector6 * v)
{
  v->x[0] = v->x[1] = v->x[2] = v->x[3] = v->x[4] = v->x[5] = 0.;
}

Vector6 vector6_sum (Vector6 v0, Vector6 v1)
{
  Vector6 vsum;
  gint i;
  
  for ( i = 0; i < 6; i++)
    vsum.x[i] = v0.x[i] + v1.x[i];

  return vsum;
}

void vector6_add (Vector6 * v0, Vector6 * v1)
{
  gint i; 
  for ( i = 0; i < 6; i++)
    v0->x[i] += v1->x[i];
}

Vector6 vector6_times_constant (Vector6 v, gdouble c)
{
  Vector6 vsum;
  gint i;
  
  for ( i = 0; i < 6; i++)
    vsum.x[i] =  c*v.x[i];

  return vsum;
}

void vector_multiply_by_constant (Vector6 * v, gdouble c)
{
  gint i;
  
  for ( i = 0; i < 6; i++)
    v->x[i] *=  c;
}

/* Transformation methods */

/**
 * Fills the Euler matrix euler_m of the transfomation
 **/
void transformation_euler_matrix  (Transformation * t)
{
  if (t != NULL) {
    t->euler_m.a[0][0] = cos(t->euler.y)*cos(t->euler.z);
    t->euler_m.a[1][0] = -cos(t->euler.x)*sin(t->euler.z) + sin(t->euler.x)*sin(t->euler.y)*cos(t->euler.z);
    t->euler_m.a[2][0] = sin(t->euler.x)*sin(t->euler.z) + cos(t->euler.x)*sin(t->euler.y)*cos(t->euler.z);
    t->euler_m.a[0][1] = cos(t->euler.y)*sin(t->euler.z);
    t->euler_m.a[1][1] = cos(t->euler.x)*cos(t->euler.z) + sin(t->euler.x)*sin(t->euler.y)*sin(t->euler.z);
    t->euler_m.a[2][1] = -sin(t->euler.x)*cos(t->euler.z) + cos(t->euler.x)*sin(t->euler.y)*sin(t->euler.z);
    t->euler_m.a[0][2] = -sin(t->euler.y);
    t->euler_m.a[1][2] = sin(t->euler.x)*cos(t->euler.y);
    t->euler_m.a[2][2] = cos(t->euler.x)*cos(t->euler.y);
  }
}

/**
 * Applies a transformation (rotation, translation) to a Point.
 * Returns the transformed Point.
 * This assumes that the Euler matrix was processsed beforehand.
 **/
Point transform_point (const Point p, Point * xg, Vector * t, Matrix3 * euler_m)
{
  Vector v;
  Point p0 = p;

  v.x = p.x - xg->x;
  v.y = p.y - xg->y;
  v.z = p.z - xg->z;

  v = vector_rotate (euler_m, v);

  p0.x = xg->x + v.x + t->x;
  p0.y = xg->y + v.y + t->y;
  p0.z = xg->z + v.z + t->z;

  return p0;
}

/**
 * Applies a transformation (only rotation) to a Vector.
 * Returns the transformed Vector.
 * This assumes that the Euler matrix was processsed beforehand.
 **/
Vector transform_vector (Vector v, Matrix3 * euler_m)
{
  return vector_rotate (euler_m, v);
}

/* Panel methods */

Panel * panel_new ()
{
  Panel * new = g_malloc (sizeof(Panel));
  new->border = FALSE;
  return new;
}

void panel_destroy (Panel * p)
{
  g_assert (p);
  g_free (p);
}

gboolean panel_read (Panel * p, FILE * fp)
{
  float x, y, z;
  gint i;
  gdouble px, py, pz;
  px = py = pz = 0.;

  g_assert (fp != NULL);

  for ( i = 0; i < 4; i++) {
    if(fscanf(fp, "%f %f %f\n", &x, &y, &z) != EOF) {
      p->p[i].x = x;
      p->p[i].y = y;
      p->p[i].z = z;
      px += x; py += y; pz += z;
    }
    else
      return FALSE;
  }


  p->var[0] = px/4.;
  p->var[1] = py/4.;
  p->var[2] = pz/4.;

  return TRUE;
}

Vector panel_first_order_normal (Panel * panel)
{
  Vector v0, v1;

  v0.x = panel->p[0].x-panel->var[0];
  v0.y = panel->p[0].y-panel->var[1];
  v0.z = panel->p[0].z-panel->var[2];
  v1.x = panel->p[1].x-panel->var[0];
  v1.y = panel->p[1].y-panel->var[1];
  v1.z = panel->p[1].z-panel->var[2];

  /* vector product should give the direction of the normal of the flat (non-spline) panel */
  /* if Rhino exports panels in a consistant way                                           */
  return vector_vector_product (&v0, &v1);
}

/*******************************/
/*  BSpline panel evaluation   */
/*******************************/

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

gdouble panel_eval (Panel * p, gdouble i, gdouble j, gint var)
{
  /* fprintf(stdout,"Spline: %e %e %e \n", j-p->j, i-p->i, b2 (j-p->j, 1.)*b2 (i-p->i, 1.)); */
  return PANEL_VAL(p,var)*b2 (j-p->j, 1.)*b2 (i-p->i, 1.);
}

/* Unidirectional bspline interpolation for border panels */
gdouble panel_eval_y (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (i-p->i) < 0.5)
    return PANEL_VAL(p,var)*b2 (j-p->j, 1.);
  else
    return 0.;
}

gdouble panel_eval_x (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (j-p->j) < 0.5)
    return PANEL_VAL(p,var)*b2 (i-p->i, 1.);
  else
    return 0.;
}

/*************************************************/
/*  BSpline panel first derivatives evaluation   */
/*************************************************/

static gdouble b2_x (gdouble x, gdouble h)
{
   if (x <= -3.*h/2.)
    return 0.;
  if (x > 3.*h/2.)
    return 0.;

  if (-3.*h/2. < x && x <= -h/2.)
    return 1./(h*h) * (x + 3.*h/2.);

  if (-h/2. < x && h/2. >= x)
    return (-2.*x)/(h*h);

  if (3.*h/2. >= x && x > h/2.)
    return 1./(h*h) * (x - 3.*h/2.);

  g_assert_not_reached ();
}

gdouble panel_eval_dx (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b2 (j-p->j, 1.)*b2_x (i-p->i, 1.);
}

gdouble panel_eval_dy (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b2_x (j-p->j, 1.)*b2 (i-p->i, 1.);
}

/* Unidirectional bspline interpolation for border panels */
gdouble panel_eval_dx_x (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (j-p->j) < 0.5)
    return PANEL_VAL(p,var)*b2_x (i-p->i, 1.);
  else
    return 0.;
}

gdouble panel_eval_dy_y (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (i-p->i) < 0.5)
    return PANEL_VAL(p,var)*b2_x (j-p->j, 1.);
  else
    return 0.;
}

/**************************************************/
/*  BSpline panel second derivatives evaluation   */
/**************************************************/
static gdouble b2_xx (gdouble x, gdouble h)
{
   if (x <= -3.*h/2.)
    return 0.;
  if (x > 3.*h/2.)
    return 0.;

  if (-3.*h/2. < x && x <= -h/2.)
    return 1./(h*h);

  if (-h/2. < x && h/2. >= x)
    return -2./(h*h);

  if (3.*h/2. >= x && x > h/2.)
    return 1./(h*h);

  g_assert_not_reached ();
}

gdouble panel_eval_dxx (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b2 (j-p->j, 1.)*b2_xx (i-p->i, 1.);
}

gdouble panel_eval_dxx_x (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (j-p->j) < 0.5)
    return PANEL_VAL(p,var)*b2_xx (i-p->i, 1.);
  else
    return 0.;
}

gdouble panel_eval_dyy (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b2_xx (j-p->j, 1.)*b2 (i-p->i, 1.);
}

gdouble panel_eval_dyy_y (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (i-p->i) < 0.5)
    return PANEL_VAL(p,var)*b2_xx (j-p->j, 1.);
  else
    return 0.;
}

gdouble panel_eval_dxy (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b2_x (j-p->j, 1.)*b2_x (i-p->i, 1.);
}

/* Panel bspline integration */

static gdouble b2_int (gdouble x, gdouble h)
{
  if (x <= -3.*h/2.)
    return 0.;
  if (x > 3.*h/2.)
    return 0.;

  if (-3.*h/2. < x && x <= -h/2.)
    return h/6.;

  if (-h/2. < x && h/2. >= x)
    return 2.*h/3.;

  if (3.*h/2. >= x && x > h/2.)
    return h/6.;
  g_assert_not_reached ();
}

gdouble panel_eval_int (Panel * p, gdouble i, gdouble j, gint var)
{
  return PANEL_VAL(p,var)*b2_int (j-p->j, 1.)*b2_int (i-p->i, 1.);
}

/* Unidirectional bspline interpolation for border panels */
gdouble panel_eval_int_y (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (i-p->i) < 0.5)
    return PANEL_VAL(p,var)*b2_int (j-p->j, 1.);
  else
    return 0.;
}

gdouble panel_eval_int_x (Panel * p, gdouble i, gdouble j, gint var)
{
  if (fabs (j-p->j) < 0.5)
    return PANEL_VAL(p,var)*b2_int (i-p->i, 1.);
  else
    return 0.;
}

/* DCurve methods */
/**
 * Creates a new @DCurve structure of size N
 */
DCurve * dcurve_new (gint N)
{
  DCurve * new = g_malloc (sizeof(DCurve));
  new->p = g_array_sized_new (FALSE, TRUE, sizeof(Point), N);
  return new;
}

/**
 * Destroys a @DCurve structure
 */
void dcurve_destroy (DCurve * c)
{
  g_assert (c != NULL);
  g_array_free (c->p, FALSE);
  g_free (c);
}

void dcurve_print (DCurve * dc, FILE * fp)
{
  gint i;

  for ( i = 0; i < dc->p->len; i++) {
    Point point = g_array_index (dc->p, Point, i);
    fprintf (fp, "%g %g %g %g\n", point.x, point.y, point.z, point.xi);
  }
}

static Point interpolate_point (GSList * l, gdouble xi)
{
  Point * pp0, * pp1;

  pp0 = l->data;
  while (l) {
    pp1 = l->data;
    
    if (pp1->xi == xi)
      return *pp1;

    if ( pp0->xi <= xi && pp1->xi >= xi) {
      Point p;
      p.x = ((pp1->xi -xi)*pp0->x + (xi - pp0->xi)*pp1->x)/(pp1->xi-pp0->xi);
      p.y = ((pp1->xi -xi)*pp0->y + (xi - pp0->xi)*pp1->y)/(pp1->xi-pp0->xi);
      p.z = ((pp1->xi -xi)*pp0->z + (xi - pp0->xi)*pp1->z)/(pp1->xi-pp0->xi);
      p.xi = xi;
      return p;
    }
    pp0 = pp1;
    l = l->next;
  }

  //g_assert_not_reached ();
}

static Vector tangent (Point * p0, Point * p1)
{
  Vector ti;

  ti.x = (p1->x -p0->x)/(p1->xi-p0->xi);
  ti.y = (p1->y -p0->y)/(p1->xi-p0->xi);
  ti.z = (p1->z -p0->z)/(p1->xi-p0->xi);
  return ti;
}

DCurve * resample_point_list_hybrid (GSList * l, gint m)
{
  /* gdouble lambda_s = 0.2, lambda_k = 0.8; */
  gdouble lambda_s = /* 0.9 */1., lambda_k = /* 0.1 */0.;
  /* DCurve * dc_new = dcurve_new (m); */
  /* arc length */
  GArray * s = g_array_new (FALSE, TRUE, sizeof(gdouble));
  Point * p0, * p1;
  gint i;
  gint N = g_slist_length (l);


  GSList * ll = l;
  /* STEP 1: Compute arc length , rescale so that max arc length is 1 */
  /* add to xi weighted by lambda_s */
  ll = l->next;
  p0 = l->data;
  gdouble length = 0.;
  while (ll) {
    p1 = ll->data;
    length += point_distance (*p0, *p1);
    p1->xi = length;
    p0 = p1;
    ll = ll->next;
  }

  ll = l->next;
  while (ll) {
    p0 = ll->data;
    p0->xi *= lambda_s/length;
    ll = ll->next;
  }

  /* STEP 2: Compute curvature arc length on fine grid. Check if curve has non-    */
  /* trivial amount of curvature. If so, normalize to m, and add into xi, weighted */
  /* by lambda_k. Otherwise, use arc length instead                                */
  GArray * t = g_array_new (FALSE, TRUE, sizeof(Vector));
  GArray * a = g_array_new (FALSE, TRUE, sizeof(gdouble));
  Vector ti;
  
  Point * p2;
  p0 = l->data;
  p1 = l->next->data;
  p2 = g_slist_last (l)->data;
  ll = l->next->next;
  ti = vector_normalise (tangent (p2, p1));
  g_array_append_val (t, ti);
  while (ll) {
    p2 = ll->data;

    ti = vector_normalise (tangent (p0, p2));
    g_array_append_val (t, ti);
    p0 = p1;
    p1 = p2;

    ll = ll->next;
  }
  p2 = l->data;
  ti = vector_normalise (tangent (p0, p2));
  g_array_append_val (t, ti);
  g_assert (t->len == N);

 
  gdouble ai = 0.;
  Vector t0, t1;
  t0 = g_array_index (t, Vector, 0);
  g_array_append_val (a, ai);
  for ( i = 1; i < t->len; i++) {
    t1 = g_array_index (t, Vector, i);
    ai = g_array_index (a, gdouble, i-1) + pow((t1.x-t0.x)*(t1.x-t0.x) +
  						(t1.y-t0.y)*(t1.y-t0.y) +
					       (t1.z-t0.z)*(t1.z-t0.z), 0.5);
    g_array_append_val (a, ai);
    t0 = t1;
  }


  ll = l->next;
  if (g_array_index (a, gdouble, a->len-1) > 0.01) {
    for ( i = 1; i < a->len; i++) {
      p0 = ll->data;
      g_array_index (a, gdouble, i) = g_array_index (a, gdouble, i)
  	/g_array_index (a, gdouble, a->len-1);
      p0->xi += lambda_k * g_array_index (a, gdouble, i);
      ll = ll->next;
    }
  }
  else {
    while (ll) {
      p0 = ll->data;
      p0->xi /= lambda_s;
      ll = ll->next;
    }
  }

  /* Skip STEP 3: Attractor Points */
  

  /* STEP 4: Obtain point distribution by inverting grid function */
  DCurve * dc_new = dcurve_new (m);
  Point * pp0 = l->data;
  Point pstart = *pp0;
  g_array_append_val (dc_new->p, pstart);
  for ( i = 1; i < m-1; i++) {
    Point p0 = interpolate_point (l, i/(m-1.));
    g_array_append_val (dc_new->p, p0);
  }
  g_array_append_val (dc_new->p, pstart);

  return dc_new;

}

DCurve * resample_point_list_centripetal (GSList * l, gint m)
{
  GSList * ll;
  Point p0, p1;
  gint i;

  ll = l;
  Point * pp1 = NULL;
  Point * pp2 = NULL;
  GSList * sharp_points = NULL;
  while (ll) {
    Point * p = ll->data;
    /* p0 = *p; */
 

    if (pp1 != NULL && pp2 != NULL) {
      Vector v1, v2;
      v1.x = p->x - pp1->x;
      v1.y = p->y - pp1->y;
      v1.z = p->z - pp1->z;
      v2.x = pp1->x - pp2->x;
      v2.y = pp1->y - pp2->y;
      v2.z = pp1->z - pp2->z;
      v1 = vector_normalise (v1);
      v2 = vector_normalise (v2);
      
      gdouble angle = fabs(fabs (vector_scalar_product (&v1 , &v2 )) - 1.);
      if ( angle > 0.01) {
    	/* fprintf (stdout, " %f %f %f %f\n", (pp1->x+p->x)/2., (pp1->y+p->y)/2., (pp1->z+p->z)/2., angle); */

    	Point * sharp_point = g_malloc (sizeof(Point));
    	sharp_point->x = /* 1.001* */(pp1->x+p->x)/2.;
    	sharp_point->y = /* 1.001* */(pp1->y+p->y)/2.;
    	sharp_point->z = /* 1.001* */(pp1->z+p->z)/2.;

    	l = g_slist_insert_before (l, ll, sharp_point);
    	sharp_points = g_slist_append (sharp_points, sharp_point);

    	pp1 = NULL;
      	pp2 = NULL;
    	p = NULL;
      }
    }

    /* g_array_append_val (dc->p, p0); */
    pp2 = pp1;
    pp1 = p;
    ll = ll->next;
  }

  /* Compute arc length , rescale so that max arc length is 1 */
  /* add to xi weighted by lambda_s */
  gdouble l0 = 0.;
  Point * pp0;

  ll = l;
  pp0 = l->data;
  pp0->xi = 0.;
  

  ll = l->next;
  while (ll) {
    pp1 = ll->data;
    
    l0 += point_distance (*pp1, *pp0);

    GSList * tmp = sharp_points;
    while (tmp) {
      Point * psharp = tmp->data;
      //l0 += exp(-5*pow(point_distance(*pp1,*psharp),2.));
      //l0 += exp(-6.5*pow(point_distance(*pp1,*psharp),2));
      tmp = tmp->next;
    }

    pp1->xi = l0;

    pp0 = pp1;
    ll = ll->next;
  }


  ll = l->next;
  while (ll) {
    pp0 = ll->data;
    pp0->xi /= l0;
    ll = ll->next;
  }

  g_assert (pp0->xi == 1.);


  ll = sharp_points;
  /* while (ll) { */
    GSList * list = l;
    pp1 = ll->data;
    
    gint ll1 = pp1->xi*(m-1);
    
    if (fabs (pp1->xi- ll1/(m-1)) > fabs (pp1->xi -(ll1+1)/(m-1.)))
      ll1 = ll1+1;

    //fprintf (stdout, "pos: %i of %i %f \n", ll1, m, pp1->xi);

    gdouble oldxi = pp1->xi;
    gboolean past = FALSE;
    while (list) {

      Point * ptmp = list->data;
      if (!past)
      	ptmp->xi *= ll1/(m-1.)/oldxi;
      else {
      	ptmp->xi =  pp1->xi + (1.-pp1->xi)*(ptmp->xi-oldxi)/(1.-oldxi);
	//fprintf (stderr, "PAST %f %f\n", ptmp->xi, oldxi);
      }
      
      if (list->data == ll->data)
      	past = TRUE;

      list = list->next;
    }

    //fprintf (stderr, " %f %f %f \n", pp1->xi, 25./(m-1.), ll1/(m-1.));

    ll = ll->next;
    /* } */


  DCurve * dc_new = dcurve_new (m);
  pp0 = l->data;
  Point pstart = *pp0;
  g_array_append_val (dc_new->p, pstart);
  for ( i = 1; i < m-1; i++) {
    p0 = interpolate_point (l, i/(m-1.));
    g_array_append_val (dc_new->p, p0);
  }
  g_array_append_val (dc_new->p, pstart);

  return dc_new;
}


/* BSpline methods */
/**
 * Creates a BSpline structure and allocates
 * memory for the gsl structures that will allow
 * to evaluate the bspline and its first derivative.
 * It is assumed that the bspline has 20 degrees of
 * freedom.
 * Returns a pointer on a @BSpline structure.
 */
BSpline * bspline_new ()
{
  BSpline * new = g_malloc (sizeof(BSpline));
  
  new->ncoeffs = /* 25 */10;
  /* allocate a cubic bspline workspace (k = 4, dc->x->len-2 break points) */
  new->bi = gsl_vector_alloc(new->ncoeffs);
  new->ci = gsl_vector_alloc(new->ncoeffs);
  new->bdi = gsl_matrix_alloc(new->ncoeffs, 3);
  new->cov = gsl_matrix_alloc(new->ncoeffs, new->ncoeffs);
  new->bw = gsl_bspline_alloc (4, new->ncoeffs - 2);
  new->bdw = gsl_bspline_deriv_alloc (4);

  return new;
}

/**
 * Deallocates memory and destroys a BSpline structure
 */
void bspline_destroy (BSpline * bs)
{
  g_assert (bs != NULL);
  gsl_vector_free (bs->bi);
  gsl_vector_free (bs->ci);
  gsl_matrix_free (bs->bdi);
  gsl_matrix_free (bs->cov);
  gsl_bspline_free (bs->bw);
  gsl_bspline_deriv_free (bs->bdw);
  g_free (bs);
}

/**
 * Evaluate the @BSpline at position xi
 */
gdouble bspline_eval (BSpline * bs, gdouble xi)
{
  gdouble val, err;
  gsl_vector * col = gsl_vector_alloc (bs->ncoeffs);
  
  /* spline */
  gsl_bspline_eval (xi, bs->bi, bs->bw);
  gsl_multifit_linear_est (bs->bi, bs->ci, bs->cov, &val, &err);

  gsl_vector_free (col);
  return val;
}

/**
 * Evaluates the first derivative of a @BSpline at position xi
 */
gdouble bspline_eval_first_derivative (BSpline * bs, gdouble xi)
{
  gdouble val = 0., err;
  gsl_vector * col = gsl_vector_alloc (bs->ncoeffs);

  /* first derivative of spline */
  gsl_bspline_deriv_eval (xi, 1, bs->bdi, bs->bw, bs->bdw);
  gsl_matrix_get_col (col, bs->bdi, 1);
  gsl_multifit_linear_est (col, bs->ci, bs->cov, &val, &err);
  
  gsl_vector_free (col);
  return val;
}

/**
 * Fits a bspline curve of ncoeffs degrees of freedom
 * to the  the curve defined by
 * the coordinate (x,y) respectively stored as @GArray
 * of gdouble.
 * Returns a pointer on a @BSpline
*/
static BSpline * bspline_fit (GArray * x, GArray * y)
{
  g_assert (x != NULL);
  g_assert (y != NULL);
  gint i, j;
  BSpline * bs = bspline_new ();
  
  /* create uniformly spaced knots */
  gsl_bspline_knots_uniform(0., 1., bs->bw);
  
  /* copies values to fit */
  gsl_vector * val_y = gsl_vector_alloc(y->len);
  for ( i = 0; i < y->len; i++) {
    gsl_vector_set (val_y, i, g_array_index (y, gdouble, i));
  }

  /* Linear fitting of the spline */
  gsl_multifit_linear_workspace * mw;
  mw = gsl_multifit_linear_alloc(x->len, bs->ncoeffs);
  gsl_matrix * X = gsl_matrix_alloc(x->len, bs->ncoeffs);
  gdouble chisq;
  
  /* construct the fit matrix X */
  for (i = 0; i < x->len; ++i) {
    /* compute B_j(ti) for all j */
    gsl_bspline_eval (g_array_index (x, gdouble, i), bs->bi, bs->bw);
    /* fill in row i of X */
    for (j = 0; j < bs->ncoeffs; ++j) {
      double Bxj = gsl_vector_get(bs->bi, j);
      gsl_matrix_set(X, i, j, Bxj);
    }
  }
  /* do the fit */
  gsl_multifit_linear (X, val_y, bs->ci, bs->cov, &chisq, mw);
  
  /* free structures */
  gsl_multifit_linear_free (mw);
  gsl_matrix_free (X);

  return bs;
}

Point bcurve_eval (gdouble u, BCurve * bc)
{
  Point p;

  p.xi = u;
  p.x = bspline_eval (bc->bsx, u);
  p.y = bspline_eval (bc->bsy, u);
  p.z = bspline_eval (bc->bsz, u);

  return p;
}

BCurve * bcurve_fit_point_list (GSList * l)
{
  BCurve * bc = g_malloc (sizeof(BCurve));
  gint i;
  GSList * ll;
  gdouble l0 = 0;
  GArray * x = g_array_new (TRUE, TRUE, sizeof(gdouble));
  GArray * y = g_array_new (TRUE, TRUE, sizeof(gdouble));
  GArray * z = g_array_new (TRUE, TRUE, sizeof(gdouble));
  GArray * xi = g_array_new (TRUE, TRUE, sizeof(gdouble));

  ll = l;
  while (ll) {
    Point * p = ll->data;

    g_array_append_val (x, p->x);
    g_array_append_val (y, p->y);
    g_array_append_val (z, p->z);
    l0 += sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
    g_array_append_val (xi, l0);

    ll = ll->next;
  }
  
  for ( i = 0; i < xi->len; i++)
    g_array_index (xi, gdouble, i) = g_array_index (xi, gdouble, i)/l0;

  bc->bsx = bspline_fit (xi, x);
  bc->bsy = bspline_fit (xi, y);
  bc->bsz = bspline_fit (xi, z);

  return bc;
}

/****************************************************************************/


GArray * new_2D_double_array (gint M, gint N)
{
  GArray * new = g_array_new (FALSE, FALSE, sizeof(gdouble));
  gint i;
  gdouble v0 = 0.;
  
  for ( i = 0; i < M*N; i++)
    g_array_append_val (new, v0);

  g_assert (new->len == M*N);
  
  return new;
}

GArray * new_2D_vector_array (gint M, gint N) 
{
  GArray * new = g_array_new (FALSE, FALSE, sizeof(Vector));
  gint i;
  Vector v0;
  v0.x = v0.y = v0.z = 0.;
  
  for ( i = 0; i < M*N; i++)
    g_array_append_val (new, v0);

  g_assert (new->len == M*N);

  return new;
}

/****************************************************************************/

gdouble coeff (Spline2D * s, gint i, gint j, gint var)
{
  SplineCoeffs * sc;

  if (var < 3)
    sc = g_ptr_array_index (s->coeffs, i+ j*s->NXU);
  else {
    if (s->periodic)
      sc = g_ptr_array_index (s->coeffs, i+ j*(s->NU+s->k-1));
    else
      sc = g_ptr_array_index (s->coeffs, i+ j*s->NU);
  }
  
  return sc->v[var];
}

gdouble coeff_x (Spline2D * s, gint i, gint j, gint var)
{
  SplineCoeffs * sc;

  sc = g_ptr_array_index (s->coeffs, i+ j*(s->NXU));
  
  return sc->v[var];
}

void coeff_assign (Spline2D * s, gint i, gint j, gint var, gdouble coeff)
{
  SplineCoeffs * sc;

  if (var < 3)
    sc = g_ptr_array_index (s->coeffs, i + j*s->NXU);
  else
    sc = g_ptr_array_index (s->coeffs, i + j*s->NU);

  sc->v[var] = coeff;
}

void coeff_set_var_to_zero (Spline2D * splines, gint var)
{
  gint i, j;
  Spline2D * s = splines;

  if (var < 3) {
    while (s) {
      for ( i = 0; i < s->NXU; i++) {
	for ( j = 0; j < s->NXV; j++) {
	  SplineCoeffs * sc = g_ptr_array_index (s->coeffs, i + j*s->NXU);
	  sc->v[var] = 0.;
	}
      }
      s = s->next;
    }
  }
  else {
    while (s) {
      for ( i = 0; i < s->NU; i++) {
	for ( j = 0; j < s->NV; j++) {
	  SplineCoeffs * sc = g_ptr_array_index (s->coeffs, i + j*s->NU);
	  sc->v[var] = 0.;
	}
      }
      s = s->next;
    }
  }
}

void coeff_set_var_to_constant (Spline2D * splines, gint var, gdouble val)
{
  gint i, j;
  Spline2D * s = splines;

  while (s) {
    for ( i = 0; i < s->M+s->k-1; i++) {
      for ( j = 0; j < s->N+s->k-1; j++) {
	coeff_assign (s, i, j, var, val);
      }
    }
    s = s->next;
  }
}

void spline2d_copy_var (Spline2D * splines, gint var1, gint var2)
{
  gint i, j;
  Spline2D * sp = splines;

  g_assert ( var1 > 3 && var2 > 3);
  // Otherwise NU is NXU
  
  while (sp) {
    for ( i = 0; i < sp->NU; i++) {
      for ( j = 0; j < sp->NV; j++) {
	SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, i + j*sp->NU);
	sc->v[var2] = sc->v[var1];
      }
    }
    sp = sp->next;
  }
}

void spline2d_list_copy_var (GSList * list, gint var1, gint var2)
{
  while (list) {
    spline2d_copy_var (list->data, var1, var2);
    list = list->next;
  }
}

gint spline2d_size (Spline2D * sp)
{
  return sp->NU*sp->NV;
}

gsl_vector * spline2d_rhs_vector (Spline2D * sp, gsl_vector * rhs)
{
  if (rhs == NULL)
    return gsl_vector_alloc (sp->size(sp));

  if (rhs->size == spline2d_size (sp))
    return rhs;
  
  gsl_vector_free (rhs);
  return gsl_vector_alloc (sp->size(sp));
}

/**
 * N: Number of panels along x
 * M: Number of panels along y
 * k: spline order
 **/
Spline2D * spline2d_new (gint M, gint N, gint k, gint ninner, gint nouter)
{
  Spline2D * s = g_malloc (sizeof(Spline2D));
  
  // Order of the splines
  s->k = k;
  s->N = N;
  s->M = M;
  s->NUT = 0;
  s->ninner = ninner;
  s->nouter = nouter;

  // Periodicity of spline
  s->periodic = FALSE;
  s->noflux = FALSE;
  s->next = NULL;
  s->fs_index = 0;
  s->fs_index_x = 0;

  // Number of knots so that it fits with the notation of Maniar et al,
  s->w_u = gsl_bspline_alloc (k, M+1);
  s->wd_u = gsl_bspline_deriv_alloc (k);

  s->w_v = gsl_bspline_alloc (k, N+1);
  s->wd_v = gsl_bspline_deriv_alloc (k);

  s->wx_u = s->w_u;
  s->wxd_u = s->wd_u;

  // break = M+1/N+1
  // NB: GSL collocate k knots at the start and end of the spline
  gsl_bspline_knots_uniform (0., 1., s->w_u);
  gsl_bspline_knots_uniform (0., 1., s->w_v);

  // Store number of spline coefficients in each direction
  s->NU = gsl_bspline_ncoeffs (s->w_u);
  s->NV = gsl_bspline_ncoeffs (s->w_v);

  s->NXU = s->NU;
  s->NXV = s->NV;

  // Allocates space for the coeffs of the splines.
  s->coeffs = g_ptr_array_new ();
  gint i, j;
  for ( i = 0; i < s->NU; i++) {
    for ( j = 0; j < s->NV; j++) {
      SplineCoeffs * new = g_malloc0 (sizeof(SplineCoeffs));
      g_ptr_array_add (s->coeffs, new);
    }
  }

  // Creates panels
  s->panels = g_ptr_array_new ();
  for ( i = 0; i < M; i++) {
    for ( j = 0; j < N; j++) {
      SPPanel * spp = sppanel_new (k);
      g_ptr_array_add (s->panels, spp);
    }
  }

  // Associated methods
  s->build_fit_matrix = spline2d_build_galerkin_fit_matrix;
  s->build_fit_rhs = build_galerkin_rhs_gauss;
  s->build_fit_noflux_matrix = spline2d_build_galerkin_fit_matrix;
  s->build_fit_noflux_rhs = build_galerkin_rhs_gauss;
  s->copy_fit_solution = spline2d_copy_problem_solution;
  s->add_fit_solution = spline2d_add_problem_solution;
  s->reinit_panels = spline2d_reinit_panels_physical_quantities;
  s->size = spline2d_size;
  s->rhs_vector = spline2d_rhs_vector;

  s->fit = NULL;
  s->fit_noflux = NULL;
  s->gr = NULL;
  s->fully_submerged = FALSE;
  s->hull_patch = NULL;
  s->rhs = NULL;

  s->istart = 999999999; // Will make the code crash if not reassigned

  s->fit_tmp = NULL;

  return s;
}

Spline2D * spline2d_symmetrical_y (Spline2D * sp, gdouble y)
{
  Spline2D * new = spline2d_new (sp->M, sp->N, sp->k, sp->ninner, sp->nouter);
  
  gint i, j;
  for ( i = 0; i < sp->NXU; i++) {
    for ( j = 0; j < sp->NXV; j++) {
      SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, (sp->NXU - 1 - i) + j*sp->NXU);
      SplineCoeffs * scnew = g_ptr_array_index (new->coeffs, i + j*sp->NXU);
      scnew->v[0] = sc->v[0];
      scnew->v[1] = -(sc->v[1]-y) + y;
      scnew->v[2] = sc->v[2];
    }
  }

  spline2d_init_panels (new);

  return new;
}

void spline2d_translate (Spline2D * sp, gdouble x, gdouble y, gdouble z)
{
  gint i, j;
  for ( i = 0; i < sp->NXU; i++) {
    for ( j = 0; j < sp->NXV; j++) {
      SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, i + j*sp->NXU);
      sc->v[0] = sc->v[0] + x;
      sc->v[1] = sc->v[1] + y;
      sc->v[2] = sc->v[2] + z;
    }
  }
  spline2d_reinit_panels_physical_quantities (sp);
}

void spline2d_rescale (Spline2D * sp, gdouble sx, gdouble sy, gdouble sz)
{
  gint i, j;
  for ( i = 0; i < sp->NXU; i++) {
    for ( j = 0; j < sp->NXV; j++) {
      SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, i + j*sp->NXU);
      sc->v[0] = sc->v[0]*sx;
      sc->v[1] = sc->v[1]*sy;
      sc->v[2] = sc->v[2]*sz;
    }
  }
  spline2d_reinit_panels_physical_quantities (sp);
}

void spline2d_rotate (Spline2D * sp, Point xc, gdouble angle)
{
  gdouble ca = cos(angle);
  gdouble sa = sin(angle);

  gint i, j;
  for ( i = 0; i < sp->NXU; i++) {
    for ( j = 0; j < sp->NXV; j++) {
      SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, i + j*sp->NXU);
      gdouble x = sc->v[0]-xc.x;
      gdouble y = sc->v[1]-xc.y;

      sc->v[0] = xc.x + ca*x - sa*y;
      sc->v[1] = xc.y + sa*x + ca*y;
    }
  }
  spline2d_reinit_panels_physical_quantities (sp);
}

void spline2d_destroy (Spline2D * sp)
{
  g_assert (sp != NULL);
  gint i;

  for ( i = 0; i < sp->coeffs->len; i++)
    g_free (g_ptr_array_index (sp->coeffs, i));
  g_ptr_array_free (sp->coeffs, TRUE);

  for ( i = 0; i < sp->panels->len; i++) {
    sppanel_destroy (g_ptr_array_index (sp->panels, i));
  }
  g_ptr_array_free (sp->panels, TRUE);

  gsl_bspline_free (sp->w_u);
  gsl_bspline_deriv_free (sp->wd_u);
  gsl_bspline_free (sp->w_v);
  gsl_bspline_deriv_free (sp->wd_v);

  if (sp->fit != NULL) {
    ccs_problem_destroy (sp->fit);
  }

  g_free (sp);
}

/**
 * Returns the sum of the contributions of all of the splines at location (u,v)
 **/
gdouble spline2d_eval (Spline2D * s, gdouble u, gdouble v, gint var)
{
  g_assert ( s->w_u != NULL);
  g_assert ( s->w_v != NULL);

  gint i, j;
  size_t istart, iend, jstart, jend;
  gsl_vector * Bu = gsl_vector_alloc (s->k);
  gsl_vector * Bv = gsl_vector_alloc (s->k);
  gdouble val = 0.;

  gsl_bspline_eval_nonzero (MAX(1e-12,MIN(1.-1e-12,u)), Bu, &istart, &iend, s->w_u);
  gsl_bspline_eval_nonzero (MAX(1e-12,MIN(1.-1e-12,v)), Bv, &jstart, &jend, s->w_v);

  if (s->periodic)
    istart -= (s->k-1);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_vector_get(Bu, i);
    gint ii = (istart+i);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = gsl_vector_get(Bv, j);
      val += coeff (s, ii, jstart+j, var)*cu*cv;
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);
  return val;
}

gdouble spline2d_eval_x (Spline2D * s, gdouble u, gdouble v, gint var)
{
  g_assert ( s->wx_u != NULL);
  g_assert_not_reached (); // Don't think the routine should be used
  gint i, j;
  size_t istart, iend, jstart, jend;
  gsl_vector * Bu = gsl_vector_alloc (s->k);
  gsl_vector * Bv = gsl_vector_alloc (s->k);
  gdouble val = 0.;

  gsl_bspline_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), Bu, &istart, &iend, s->wx_u);
  gsl_bspline_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), Bv, &jstart, &jend, s->w_v);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_vector_get(Bu, i);
    gint ii = (istart+i);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = gsl_vector_get(Bv, j);
      val += coeff (s, ii,jstart+j,var)*cu*cv;
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);
  return val;
}

gdouble spline2d_eval_gauss_point (Spline2D * s, GaussPoints * gp,
				   gint m, gint n, gint var)
{
  gdouble val = 0.;
  gint i, j;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);
  
  istart = gp->istart;
  if (s->periodic)
    istart -= (s->k-1);
  jstart = gp->jstart;

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = cu*gsl_matrix_get (Bv, j, 0);
      val += coeff (s,istart, (jstart+j),var)*cv;
    }
    istart++;
  }

  return val;
}

gdouble spline2d_eval_greville_point (Spline2D * s, GrevillePoints * gr,
				      gint m, gint n, gint var)
{
  gdouble val = 0.;
  gint i, j;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gr->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gr->Bv, n);
  
  istart = g_array_index (gr->ustart, gint, m);
  jstart = g_array_index (gr->vstart, gint, n);

  if (s->periodic)
    istart -= (s->k-1);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = cu*gsl_matrix_get (Bv, j, 0);
      val += coeff (s,istart, (jstart+j),var)*cv;
    }
    istart++;
  }

  return val;
}

/**
 * Returns coordinates (x,y,z) of the point at location (u,v) in the parametric space
 **/
Point spline2d_eval_point (Spline2D * s, gdouble u, gdouble v)
{
  g_assert ( u >= -1.-1e-8 && u <= 1.+1.e-8);
  if ( v > 1. || v < -1.)
    fprintf(stderr,"%f \n",v);

  g_assert ( v >= -1 && v <= 1.);
  g_assert ( s->wx_u != NULL);
  g_assert ( s->w_v != NULL);

  gint i, j;
  size_t istart, iend, jstart, jend;
  gsl_vector * Bu = gsl_vector_alloc (s->k);
  gsl_vector * Bv = gsl_vector_alloc (s->k);
  Point p;
  p.x = p.y = p.z = p.xi = 0.;

  gsl_bspline_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), Bu, &istart, &iend, s->wx_u);
  gsl_bspline_eval_nonzero (v, Bv, &jstart, &jend, s->w_v);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_vector_get(Bu, i);
    gint ii = istart;
    for ( j = 0; j < s->k; j++) {
      gint jj = (jstart+j);
      gdouble cv = cu*gsl_vector_get(Bv, j);
      p.x += coeff (s, ii, jj, 0)*cv;
      p.y += coeff (s, ii, jj, 1)*cv;
      p.z += coeff (s, ii, jj, 2)*cv;
    }
    istart++;
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);
  return p;
}

Point spline2d_eval_gauss_point_point (Spline2D * s, GaussPoints * gp,
				       gint m, gint n)
{
  gint i, j;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);
  Point p;
  p.x = p.y = p.z = 0.;

  istart = gp->istart;
  jstart = gp->jstart;

  if (s->periodic)
    istart -= (s->k-1);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    gint ii = istart;
    for ( j = 0; j < s->k; j++) {
      gdouble cv = cu*gsl_matrix_get (Bv, j, 0);
      gint jj = (jstart+j);
      p.x += coeff (s, ii, jj,0)*cv;
      p.y += coeff (s, ii, jj,1)*cv;
      p.z += coeff (s, ii, jj,2)*cv;
    }
    istart++;
  }

  return p;
}

/**
 * Returns the contribution of spline (i, j) at point (u,v)
 **/
gdouble spline2d_eval_spline (Spline2D * s, gint i, gint j, gdouble u, gdouble v)
{
  g_assert ( s != NULL);
  g_assert ( u >= -1 && u <= 1.);
  g_assert ( v >= -1 && v <= 1.);
  g_assert ( s->w_u != NULL);
  g_assert ( s->w_v != NULL);

  gdouble val = 0.;
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (s->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (s->w_v));
  gsl_bspline_eval (u, Bu, s->w_u);
  gsl_bspline_eval (v, Bv, s->w_v);
  
  // Can gsl_bspline_eval_nonzero be used instead ??

  val = gsl_vector_get(Bu, i)*gsl_vector_get(Bv, j);

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);
  return val;
}

/**
 * Returns the sum of the contributions of all of the splines at location (u,v) when derived
 * m times in the u direction and n times in the v one.
 **/
gdouble spline2d_derivative_eval (Spline2D * s, gdouble u, gdouble v, gint m, gint n, gint var)
{
  g_assert ( u >= -1e-8 && u <= 1.+1e-8);
  g_assert ( v >= -1e-8 && v <= 1.+1e-8);
  g_assert ( s->w_u != NULL);
  g_assert ( s->w_v != NULL);

  gint  i, j;
  gdouble val = 0.;
  size_t istart, iend, jstart, jend;
  gsl_matrix * Bu = gsl_matrix_alloc (s->k, m+1);
  gsl_matrix * Bv = gsl_matrix_alloc (s->k, n+1);

  if (var > 2) {
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), m, Bu, &istart, &iend, s->w_u, s->wd_u);
   
    if (s->periodic)
      istart -= (s->k-1);
  }
  else
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), m, Bu, &istart, &iend, s->wx_u, s->wxd_u);

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), n, Bv, &jstart, &jend, s->w_v, s->wd_v);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, m);
    gint ii = istart;
    for ( j = 0; j < s->k; j++) {
      gint jj = (jstart+j);
      gdouble cv = gsl_matrix_get (Bv, j, n);
      val += coeff (s, ii, jj,var)*cu*cv;
    }
    istart++;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  return val;
}

gdouble spline2d_derivative_eval_gauss_point (Spline2D * s, GaussPoints * gp, gint m, gint n, gint dm, gint dn, gint var)
{
  gint  i, j;
  gdouble val = 0.;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);
  istart = gp->istart;
  jstart = gp->jstart;

  if (s->periodic)
    istart -= (s->k-1);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, dm);
    gint ii = (istart+i);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = gsl_matrix_get (Bv, j, dn);
      val += coeff (s, ii, (jstart+j),var)*cu*cv;
    }
  }

  return val;
}

gdouble spline2d_derivative_eval_greville_point (Spline2D * s, GrevillePoints * gr, gint m, gint n, gint dm, gint dn, gint var)
{
  gint  i, j;
  gdouble val = 0.;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gr->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gr->Bv, n);
  istart = g_array_index (gr->ustart, gint, m);
  jstart = g_array_index (gr->vstart, gint, n);

  if (s->periodic)
    istart -= (s->k-1);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, dm);
    gint ii = (istart+i);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = gsl_matrix_get (Bv, j, dn);
      val += coeff (s, ii, (jstart+j),var)*cu*cv;
    }
  }

  return val;
}

void spline2d_fit_geometry (Spline2D * s)
{
  gint NX = 101, NY = 101;
  Point p[NX][NY];
  gint i, j, k, l;
  gint NU = s->NU;
  gint NV = s->NV;
  
  /* for ( i = 0; i < NX; i++) { */
  /*   for ( j = 0; j < NY; j++) { */
  /*     p[i][j].x = -0.5*cos(2.*M_PI*i/100); */
  /*     p[i][j].y = -0.5*sin(2.*M_PI*i/100); */
  /*     p[i][j].z = 0.5-1.*(j)/100.; */
  /*   } */
  /* } */

  for ( i = 0; i < NX; i++) {
    for ( j = 0; j < NY; j++) {
      p[i][j].x = /* 0.5* */cos(2.*M_PI*i/100)*sin(M_PI*j/100);
      p[i][j].y = /* 0.5* */sin(2.*M_PI*i/100)*sin(M_PI*j/100);
      p[i][j].z = /* 0.5* */cos(M_PI*j/100);
    }
  }

  /* for ( i = 0; i < NX; i++) { */
  /*   for ( j = 0; j < NY; j++) { */
  /*     p[i][j].x = i*NX; */
  /*     p[i][j].y = j*NY; */
  /*     p[i][j].z = 1.; */
  /*   } */
  /* } */

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
      gsl_vector_set (cx, i+ j*NU, 1.);
      gsl_vector_set (cy, i+ j*NU, 1.);
      gsl_vector_set (cz, i+ j*NU, 1.);
    }
  
  /* Set problem */
  gsl_vector * Bu = gsl_vector_alloc (gsl_bspline_ncoeffs (s->w_u));
  gsl_vector * Bv = gsl_vector_alloc (gsl_bspline_ncoeffs (s->w_v));
  for ( i = 0; i < NU; i++)
    for ( j = 0; j < NV; j++) {
      for ( k = 0; k < NX; k++)
  	for ( l = 0; l < NY; l++) {
	  gsl_bspline_eval (k/(NX-1.), Bu, s->w_u);
	  gsl_bspline_eval (l/(NY-1.), Bv, s->w_v);
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
      coeff_assign (s, i, j, 0, gsl_vector_get(cx, i+j*NU));
      coeff_assign (s, i, j, 1, gsl_vector_get(cy, i+j*NU));
      coeff_assign (s, i, j, 2, gsl_vector_get(cz, i+j*NU));
    }

  fp = fopen("fit.tmp","w");
  for ( k = 0; k < 30; k++)
  	for ( l = 0; l < 10; l++) {
	  gdouble x = k/(30-1.); /* True only is the input points are uniformly spaced */
	  gdouble y = l/(10-1.);
  	  fprintf(fp, "%f %f %f\n", spline2d_eval (s, x, y, 0), spline2d_eval (s, x, y, 1), spline2d_eval (s, x, y, 2));
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
}

S2P * s2p_new (gint k)
{
  S2P * new = g_malloc (sizeof(S2P));

  new->X = g_array_new (FALSE, FALSE, sizeof(Vector));
  new->x = g_array_new (FALSE, FALSE, sizeof(Vector));
  new->x2X = gsl_matrix_alloc (k*k, k*k);
  new->X2x = gsl_matrix_alloc (k*k, k*k);

  gint M = 2*k-2;
  // This was changed from 6 to 3
  new->w = new_2D_double_array (2*M+2, 2*M+2);
  new->n = new_2D_vector_array (2*M+2, 2*M+2);
  new->om = new_2D_vector_array (2*M+2, 2*M+2);

  return new;
}

void s2p_destroy (S2P * s2p)
{
  g_assert (s2p != NULL);
  g_array_free (s2p->X, TRUE);
  g_array_free (s2p->x, TRUE);
  gsl_matrix_free (s2p->x2X);
  gsl_matrix_free (s2p->X2x);

  g_array_free (s2p->w, TRUE);
  g_array_free (s2p->n, TRUE);
  g_array_free (s2p->om, TRUE);
}

static gdouble fact (gint n)
{
  if (n == 0)
    return 1;
  else
      return n*fact(n-1);
}

GaussPoints * gausspoints_new (gint size)
{
  GaussPoints * new = g_malloc (sizeof(GaussPoints));
  gint i;
  gdouble zero = 0.;

  new->ui = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->vj = g_array_new (FALSE, FALSE, sizeof(gdouble));

  new->Ni = g_array_new (FALSE, FALSE, sizeof(Vector));
  new->Pi = g_array_new (FALSE, FALSE, sizeof(Point));
  new->wJij = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->wij = g_array_new (FALSE, FALSE, sizeof(gdouble));

  new->c1 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->c2 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->c3 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->c4 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->c5 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->c6 = g_array_new (FALSE, FALSE, sizeof(gdouble));

  new->Bu = g_ptr_array_new ();
  new->Bv = g_ptr_array_new ();

  new->Bux = g_ptr_array_new ();
 
  for ( i = 0; i < size; i++) {
    g_array_append_val (new->ui, zero);
    g_array_append_val (new->vj, zero);
  }

  new->fsdata = g_array_new (FALSE, FALSE, sizeof(FSData));

  new->interpol = g_ptr_array_new ();

  return new;
}

void gausspoints_destroy (GaussPoints * gp)
{
  g_array_free (gp->ui, TRUE);
  g_array_free (gp->vj, TRUE);
  g_array_free (gp->Ni, TRUE);
  g_array_free (gp->Pi, TRUE);
  g_array_free (gp->wJij, TRUE);
  g_array_free (gp->wij, TRUE);
  g_array_free (gp->fsdata, TRUE);
  g_array_free (gp->c1, TRUE);
  g_array_free (gp->c2, TRUE);
  g_array_free (gp->c3, TRUE);
  g_array_free (gp->c4, TRUE);
  g_array_free (gp->c5, TRUE);
  g_array_free (gp->c6, TRUE);

  gint i;
  for ( i = 0; i < gp->Bu->len; i++)
    gsl_matrix_free (g_ptr_array_index (gp->Bu, i));
  g_ptr_array_free (gp->Bu, TRUE);
  for ( i = 0; i < gp->Bv->len; i++)
    gsl_matrix_free (g_ptr_array_index (gp->Bv, i));
  g_ptr_array_free (gp->Bv, TRUE);
  for ( i = 0; i < gp->Bux->len; i++)
    gsl_matrix_free (g_ptr_array_index (gp->Bux, i));
  g_ptr_array_free (gp->Bux, TRUE);

  for ( i = 0; i < gp->interpol->len; i++)
    g_free (g_ptr_array_index (gp->interpol, i));
  g_ptr_array_free (gp->interpol, TRUE);

  g_free (gp);
}

GrevillePoints * greville_points_new (gint M, gint N)
{
  GrevillePoints * new = g_malloc (sizeof(GrevillePoints));
  gint i;
  gdouble zero = 0.;

  new->ui = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->vj = g_array_new (FALSE, FALSE, sizeof(gdouble));

  new->Ni = g_array_new (FALSE, FALSE, sizeof(Vector));
  new->Pi = g_array_new (FALSE, FALSE, sizeof(Point));
  new->Ji = g_array_new (FALSE, FALSE, sizeof(gdouble));

  new->Bu = g_ptr_array_new ();
  new->Bv = g_ptr_array_new ();
  new->ustart = g_array_new (FALSE, FALSE, sizeof(gint));
  new->vstart = g_array_new (FALSE, FALSE, sizeof(gint));


  for ( i = 0; i < M; i++)
    g_array_append_val (new->ui, zero);

  for ( i = 0; i < N; i++)
    g_array_append_val (new->vj, zero);

  new->fsdata = g_array_new (FALSE, FALSE, sizeof(FSData));

  return new;
}

void greville_points_destroy (GrevillePoints * gr)
{
  g_array_free (gr->ui, TRUE);
  g_array_free (gr->vj, TRUE);
  g_array_free (gr->Ni, TRUE);
  g_array_free (gr->Pi, TRUE);
  g_array_free (gr->Ji, TRUE);
  g_array_free (gr->fsdata, TRUE);
  g_array_free (gr->ustart, TRUE);
  g_array_free (gr->vstart, TRUE);
  gint i;
  for ( i = 0; i < gr->Bu->len; i++)
    gsl_matrix_free (g_ptr_array_index (gr->Bu, i));
  g_ptr_array_free (gr->Bu, TRUE);
  for ( i = 0; i < gr->Bv->len; i++)
    gsl_matrix_free (g_ptr_array_index (gr->Bv, i));
  g_ptr_array_free (gr->Bv, TRUE);
  g_free (gr);
}

SPPanel * sppanel_new (gint k)
{
  SPPanel * new = g_malloc (sizeof(SPPanel));

  new->k = k;
  new->outer = NULL;

  return new;
}

void sppanel_destroy (SPPanel * spp)
{
  gint i;
  g_assert (spp != NULL);

  if (spp->outer)
  gausspoints_destroy (spp->outer);

  g_free (spp);
}

/**
 * Binomial coefficients (k n)
 * with (x+y)^n = \sum_0^k (k n) x^k y^(n-k)
 **/
static gdouble binomial_coeff (gint k, gint n)
{
  return (gdouble) fact(n)/((gdouble) fact(k)*fact(n-k));
}

/**
 * Returns the coefficient of the polynomial expansion of the surface over spp
 * for a different expension point of parametric coordinates (uc, vc)
 **/
GArray * s2p_change_polynomial_origin (S2P * s2p, gint k,
				       gdouble ue, gdouble ve,
				       gdouble uc, gdouble vc)
{
  GArray * xnew = new_2D_vector_array (k, k);
  
  // (u - um)^m = (u - uc + uc - um)^m = \sum_0^n (k m) (u - uc)^k (uc - um)^(m-k)
  gint m, n, k1, k2;
  
  for ( m = 0; m < k; m++) {
    for ( n = 0; n < k; n++) {
      
      for (k1 = 0; k1 <= m; k1++) {
	for (k2 = 0; k2 <= n; k2++) {
	  g_array_index (xnew, Vector, k1 + k2*k) = 
	    vector_sum (g_array_index (xnew, Vector, k1 + k2*k),
			vector_times_constant (g_array_index (s2p->x, Vector, m + n*k),
					       binomial_coeff (k1, m)*binomial_coeff (k2, n)
					       *pow(uc-ue, m - k1)*pow(vc-ve, n - k2)));
	}
      }   
    }
  }
  return xnew;
}

static gdouble s2p_poly_eval (S2P * s2p, SPPanel * spp, gdouble u, gdouble v, gint var)
{
  gdouble sum = 0.;
  gint i, j, m, n;
  gint k = spp->k;
      
  g_assert (s2p != NULL);
  g_assert (spp != NULL);
  g_assert ( u >= spp->u0 && u <= spp->u1);
  g_assert ( v >= spp->v0 && v <= spp->v1);

  for ( m = 0; m < k; m++) {
    for ( n = 0; n < k; n++) {
      
      for ( i = 0; i < k; i++) {
      	for ( j = 0; j < k; j++) {

	  sum += gsl_matrix_get (s2p->X2x, m + n*k, i+ j*k)*coeff(spp->sp, s2p->istart+i, s2p->jstart+j, var)
	    *pow (u - s2p->ue, m)*pow(v - s2p->ve, n);

      	}
      }

    }
  }

  return sum;
}

/**
 *
 *  Returns a SPPanel containing the coefficients X, and x of the spline and polynomial
 *  representation of the surface over the panel (ii, jj) of s. The forward and inverse
 *  transformations between X and x are also returned.
 *  \sum Xij Ni Nj = \sum xij (u-ue)^i*(v-ve)^j
 *
 **/
/* S2P * sppanel_spline_to_poly_matrix (SPPanel * spp, Spline2D * sp, */
/* 				     gdouble ue, gdouble ve) */
/* { */
/*   gint i, j, k, l; */
/*   gboolean TEST = FALSE; */

/*   gsl_matrix * Bdu = gsl_matrix_alloc (sp->k, sp->k+1); */
/*   gsl_matrix * Bdv = gsl_matrix_alloc (sp->k, sp->k+1); */

/*   gsl_matrix * A = gsl_matrix_alloc (sp->k*sp->k, sp->k*sp->k); */
/*   gsl_matrix * Ai = gsl_matrix_alloc (sp->k*sp->k, sp->k*sp->k); */
/*   gsl_matrix * B = gsl_matrix_alloc (sp->k*sp->k, sp->k*sp->k); */

/*   S2P * s2p = s2p_new (sp->k); */
/*   s2p->ue = ue; */
/*   s2p->ve = ve; */

/*   // The coefficient at a point that cannot be a Gauss point */
/*   gdouble uc = spp->u0 + 0.75*(spp->u1-spp->u0); */
/*   gdouble vc = spp->v0 + 0.75*(spp->v1-spp->v0); */

/*   g_assert ( uc != ue && vc != ve); */

/*   // The sp->k derivative of the base splines are calculated */
/*   size_t iend, jend; */
/*   gsl_bspline_deriv_eval_nonzero (uc, sp->k, Bdu, &s2p->istart, &iend, sp->w_u, sp->wd_u); */
/*   gsl_bspline_deriv_eval_nonzero (vc, sp->k, Bdv, &s2p->jstart, &jend, sp->w_v, sp->wd_v); */

/*   // Storage of the X */
/*   for ( j = 0; j < sp->k; j++) { */
/*     for ( i = 0; i < sp->k; i++) { */
/*       Vector X0; */
/*       X0.x = coeff(sp, s2p->istart+i, s2p->jstart+j, 0); */
/*       X0.y = coeff(sp, s2p->istart+i, s2p->jstart+j, 1); */
/*       X0.z = coeff(sp, s2p->istart+i, s2p->jstart+j, 2); */
/*       g_array_append_val (s2p->X, X0); */
/*     } */
/*   } */

/*   // We build A and B such that A*X = Bx */
/*   for ( i = 0; i < sp->k; i++) { */
/*     for ( j = 0; j < sp->k; j++) { */

/*       for ( k = 0; k < sp->k; k++) { */
/*   	for ( l = 0; l < sp->k; l++) { */

/*   	  if ( k <= i && l <= j ) */
/* 	    gsl_matrix_set (B, k+l*sp->k, i+j*sp->k, fact(i)/fact(i-k)* */
/* 			    pow(uc-ue,i-k)*fact(j)/fact(j-l)*pow(vc-ve,j-l)); */
/*   	  else */
/* 	    gsl_matrix_set (B, k+l*sp->k, i+j*sp->k, 0.); */

/* 	  gsl_matrix_set (A, k+l*sp->k, i+j*sp->k, */
/* 			  gsl_matrix_get (Bdu, i, k)*gsl_matrix_get (Bdv, j, l)); */
/*   	} */
/*       } */

/*     } */
/*   } */

/*   // Calculation of the matrix of the transformation X2x = B^(-1)*A */
/*   gint signum; */
/*   gsl_permutation * p = gsl_permutation_alloc (sp->k*sp->k); */

/*   gsl_matrix * Bi = gsl_matrix_alloc (sp->k*sp->k, sp->k*sp->k); */

/*   gsl_linalg_LU_decomp (B, p, &signum); */

/*   gsl_linalg_LU_invert (B, p, Bi); */

/*   for ( j = 0; j < sp->k*sp->k; j++) { */
/*     for ( i = 0; i < sp->k*sp->k; i++) { */
/*       gdouble sum = 0.; */
/*       for ( l = 0; l < sp->k*sp->k; l++) */
/*   	sum += gsl_matrix_get (Bi, i, l)*gsl_matrix_get (A, l, j); */

/*       gsl_matrix_set (s2p->X2x, i, j, sum); */
/*     } */
/*   } */

/*   // Calculation of the x vector of the polynomial representation */
/*   for ( j = 0; j < sp->k*sp->k; j++) { */
/*     Vector sum; */
/*     sum.x = sum.y = sum.z = 0.; */
/*     for ( i = 0; i < sp->k*sp->k; i++) { */
/*       sum.x += gsl_matrix_get (s2p->X2x, j, i)*g_array_index (s2p->X, Vector, i).x; */
/*       sum.y += gsl_matrix_get (s2p->X2x, j, i)*g_array_index (s2p->X, Vector, i).y; */
/*       sum.z += gsl_matrix_get (s2p->X2x, j, i)*g_array_index (s2p->X, Vector, i).z; */
/*     } */
/*     g_array_append_val (s2p->x, sum); */
/*   } */

/*   // Calculation of the matrix of the inverse transformation x2X */
/*   gsl_matrix * tmp = gsl_matrix_alloc (sp->k*sp->k, sp->k*sp->k); */

/*   gsl_matrix_memcpy (tmp, s2p->X2x); */
/*   gsl_linalg_LU_decomp (tmp, p, &signum); */
/*   gsl_linalg_LU_invert (tmp, p, s2p->x2X); */

/*   // Test */
/*   if (TEST) { */
/*     for ( j = 0; j < sp->k*sp->k; j++) { */
/*       Vector sum; */
/*       sum.x = sum.y = sum.z = 0.; */
/*       for ( i = 0; i < sp->k*sp->k; i++) { */
/*     	sum.x += gsl_matrix_get (s2p->x2X, j, i)*g_array_index (s2p->x, Vector, i).x; */
/*     	sum.y += gsl_matrix_get (s2p->x2X, j, i)*g_array_index (s2p->x, Vector, i).y; */
/*     	sum.z += gsl_matrix_get (s2p->x2X, j, i)*g_array_index (s2p->x, Vector, i).z; */
/*       } */
/*       g_array_index (s2p->X, Vector, j) = sum; */
/*     } */
/*     gdouble ucc = spp->u0 + 0.69*(spp->u1-spp->u0); */
/*     gdouble vcc = spp->v0 + 0.37*(spp->v1-spp->v0); */

/*     GArray * shifted = s2p_change_polynomial_origin (s2p, spp->k, spp->ue, spp->ve, ucc, vcc); */
/*     gdouble s2x, s2y, s2z, sx, sy, sz, s3x, s3y, s3z; */
/*     sx = sy =sz = s2x = s2y = s2z = 0.; */
/*     gdouble u, v; */

/*     fprintf(stdout, "# New panel \n"); */
/*     for ( u = spp->u0; u <= spp->u1; u += (spp->u1 - spp->u0)/5.) { */
/*       for ( v = spp->v0; v <= spp->v1; v += (spp->v1 - spp->v0)/5.) { */
/*     	sx = sy =sz = s2x = s2y = s2z = s3x = s3y = s3z = 0.; */
/*     	gsl_bspline_deriv_eval_nonzero (u, sp->k, Bdu, &s2p->istart, &iend, sp->w_u, sp->wd_u); */
/*     	gsl_bspline_deriv_eval_nonzero (v, sp->k, Bdv, &s2p->jstart, &jend, sp->w_v, sp->wd_v); */

/*     	for ( i = 0; i < sp->k; i++) { */
/*     	  for ( j = 0; j < sp->k; j++) { */
/*     	    sx += g_array_index (s2p->X, Vector, i + sp->k*j).x*gsl_matrix_get (Bdu, i, 0)*gsl_matrix_get (Bdv, j, 0); */
/*     	    sy += g_array_index (s2p->X, Vector, i + sp->k*j).y*gsl_matrix_get (Bdu, i, 0)*gsl_matrix_get (Bdv, j, 0); */
/*     	    sz += g_array_index (s2p->X, Vector, i + sp->k*j).z*gsl_matrix_get (Bdu, i, 0)*gsl_matrix_get (Bdv, j, 0); */
/*     	    s2x += g_array_index (s2p->x, Vector, i + sp->k*j).x*pow(u-ue, i)*pow(v-ve, j); */
/*     	    s2y += g_array_index (s2p->x, Vector, i + sp->k*j).y*pow(u-ue, i)*pow(v-ve, j); */
/*     	    s2z += g_array_index (s2p->x, Vector, i + sp->k*j).z*pow(u-ue, i)*pow(v-ve, j); */

/* 	    s3x += g_array_index (shifted, Vector, i + sp->k*j).x*pow(u-ucc, i)*pow(v-vcc, j); */
/* 	    s3y += g_array_index (shifted, Vector, i + sp->k*j).y*pow(u-ucc, i)*pow(v-vcc, j); */
/* 	    s3z += g_array_index (shifted, Vector, i + sp->k*j).z*pow(u-ucc, i)*pow(v-vcc, j); */
/*     	  } */
/*     	} */

/*     	/\* for ( i = 0; i < 2*sp->k-1; i++) { *\/ */
/*     	/\*   for ( j = 0; j <= i /\\* sp->k *\\/; j++) { *\/ */
/*   	/\*     if (j < sp->k && (i-j) < sp->k) { *\/ */
/*   	/\*       s2x += g_array_index (s2p->x, Vector, j + sp->k*(i-j)).x*pow(u-ue, j)*pow(v-ve, i-j); *\/ */
/*   	/\*       s2y += g_array_index (s2p->x, Vector, j + sp->k*(i-j)).y*pow(u-ue, j)*pow(v-ve, i-j); *\/ */
/*   	/\*       s2z += g_array_index (s2p->x, Vector, j + sp->k*(i-j)).z*pow(u-ue, j)*pow(v-ve, i-j); *\/ */
/*   	/\*     } *\/ */
/*     	/\*   } *\/ */
/*     	/\* } *\/ */

/*     	/\* fprintf(stdout,"%f %f %f | %f %f %f | %f %f %f | %f %f %f\n", s3x-s2x, s3y-s2y, s3z-s2z, s2x, s2y, s2z, sx, sy, sz, s3x, s3y, s3z); *\/ */

/* 	fprintf(stdout, "%f %f %f | %f %f %f \n", sx, sy, sz, s2p_poly_eval (s2p, spp, u, v, 0), s2p_poly_eval (s2p, spp, u, v, 1), s2p_poly_eval (s2p, spp, u, v, 2)); */
/*       } */
/*     } */
/*   } */
  
/*   /\* // End test *\/ */
  
/*   gsl_matrix_free (tmp); */
/*   gsl_permutation_free (p); */
/*   gsl_matrix_free (A); */
/*   gsl_matrix_free (Ai); */
/*   gsl_matrix_free (B); */
/*   gsl_matrix_free (Bi); */
/*   gsl_matrix_free (Bdu); */
/*   gsl_matrix_free (Bdv); */

/*   return s2p; */
/* } */

// Some macros to makes things look cleaner
#define q(a,b) g_array_index (sid->q, Vector, a + (b)*(N+1))
#define s(a,b) g_array_index (sid->s, gdouble, a + (b)*(24*max(N,1)))
#define p(a,b) g_array_index (sid->p, gdouble, a + (b)*(11*max(N,1)))
#define g(a,b) g_array_index (sid->g, gdouble, a + (b)*(11*max(N,1)))
#define w(a,b) g_array_index (s2p->w, gdouble, a + (b)*(3*(M-1)+1))
#define n(a,b) g_array_index (s2p->n, Vector, a + (b)*(6*(M-1)+1)) 
#define om(a,b) g_array_index (s2p->om, Vector, a + (b)*(6*(M-1)+1))
#define ui(a) g_array_index (gp->ui, gdouble, a)
#define vj(a) g_array_index (gp->vj, gdouble, a)

static double jacobian (S2P * s2p, gint k, gdouble u, gdouble v)
{
  gdouble J = 0.;
  gint m, n;
  gint M = 2*k-2;
  
  g_assert ( s2p->w != NULL);

  for ( m = 0; m <= 2*(M-1); m++) {
    for ( n = 0; n <= m; n++) {
      J += w(n,m-n)*pow(u-s2p->ue,n)*pow(v-s2p->ve,m-n);
    }
  }

  return J;
}

double spline2d_jacobian (Spline2D * sp, gdouble u, gdouble v)
{
  Vector xu, xv;
  xu.x = xu.y = xu.z = 0.;
  xv.x = xv.y = xv.z = 0.;

  gint  i, j, k = sp->k;
  size_t istart, iend, jstart, jend;
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->wx_u, sp->wxd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);

  for ( i = 0; i < k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    gdouble cdu = gsl_matrix_get (Bu, i, 1);
    gint ii = istart;
    for ( j = 0; j < k; j++) {
      gint jj = (jstart+j);
      gdouble cv = gsl_matrix_get (Bv, j, 0);
      gdouble cudv = cu*gsl_matrix_get (Bv, j, 1);
      gdouble cvdu = cv*cdu;
	
      gdouble v0 = coeff (sp,ii,jj,0);
      gdouble v1 = coeff (sp,ii,jj,1);
      gdouble v2 = coeff (sp,ii,jj,2);
      
      xu.x += v0*cvdu;
      xu.y += v1*cvdu;
      xu.z += v2*cvdu;
      xv.x += v0*cudv;
      xv.y += v1*cudv;
      xv.z += v2*cudv;
    }
    istart++;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  return vector_norm (vector_vector_product (&xu, &xv));
}

static Vector normal (S2P * s2p, gint k, gdouble u, gdouble v)
{
  Vector N;
  N.x = N.y = N.z = 0.;
  gint m, n;
  gint M = 2*k-2;

  for ( m = 0; m <= 2*(M-1); m++) {
    for ( n = 0; n <= m; n++) {
      N = vector_sum (N, vector_times_constant (n(n, m-n), pow(u-s2p->ue,n)*pow(v-s2p->ve,m-n)));
    }
  }

  return N;
}

Vector spline2d_dimensional_normal (Spline2D * sp, gdouble u, gdouble v)
{
  Vector xu, xv;
  gsl_matrix * Bu = gsl_matrix_alloc (sp->k,2);
  gsl_matrix * Bv = gsl_matrix_alloc (sp->k,2);
  size_t ustart, uend, vstart, vend;
  gint m, n;

  xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

  for ( m = 0; m < sp->k; m++) {
    gdouble cu = gsl_matrix_get (Bu, m, 0);
    gdouble cdu = gsl_matrix_get (Bu, m, 1);
    gint ii = (ustart+m);
    for ( n = 0; n < sp->k; n++) {
      gdouble cv = gsl_matrix_get (Bv, n, 0);
      gdouble cudv = cu*gsl_matrix_get (Bv, n, 1);
      gdouble cvdu = cv*cdu;
      gdouble cuv = cu*cv;
      gint jj = (vstart+n);

      SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, ii+ jj*sp->NXU);
      

      xu.x += sc->v[0]*cvdu;
      xu.y += sc->v[1]*cvdu;
      xu.z += sc->v[2]*cvdu;
      xv.x += sc->v[0]*cudv;
      xv.y += sc->v[1]*cudv;
      xv.z += sc->v[2]*cudv;

      /* xu.x += coeff (sp,ii,jj, 0)*cvdu; */
      /* xu.y += coeff (sp,ii,jj, 1)*cvdu; */
      /* xu.z += coeff (sp,ii,jj, 2)*cvdu; */
      /* xv.x += coeff (sp,ii,jj, 0)*cudv; */
      /* xv.y += coeff (sp,ii,jj, 1)*cudv; */
      /* xv.z += coeff (sp,ii,jj, 2)*cudv; */
    }
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  return vector_vector_product (&xu, &xv);
}

Vector spline2d_normal (Spline2D * sp, gdouble u, gdouble v)
{
  return vector_normalise (spline2d_dimensional_normal (sp, u, v));
}

Vector spline2d_xu (Spline2D * sp, gdouble u, gdouble v)
{
  Vector xu;
  gsl_matrix * Bu = gsl_matrix_alloc (sp->k,2);
  gsl_matrix * Bv = gsl_matrix_alloc (sp->k,2);
  size_t ustart, uend, vstart, vend;
  gint m, n;

  xu.x = xu.y = xu.z = 0.;

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
    
  for ( m = 0; m < sp->k; m++) {
    gdouble cu = gsl_matrix_get (Bu, m, 0);
    gdouble cdu = gsl_matrix_get (Bu, m, 1);
    gint ii = (ustart+m);
    for ( n = 0; n < sp->k; n++) {
      gdouble cv = gsl_matrix_get (Bv, n, 0);
      gdouble cudv = cu*gsl_matrix_get (Bv, n, 1);
      gdouble cvdu = cv*cdu;
      gdouble cuv = cu*cv;
      gint jj = (vstart+n);
      xu.x += coeff (sp,ii,jj, 0)*cvdu;
      xu.y += coeff (sp,ii,jj, 1)*cvdu;
      xu.z += coeff (sp,ii,jj, 2)*cvdu;
    }
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  return xu;
}

Vector spline2d_xv (Spline2D * sp, gdouble u, gdouble v)
{
  Vector xv;
  gsl_matrix * Bu = gsl_matrix_alloc (sp->k,2);
  gsl_matrix * Bv = gsl_matrix_alloc (sp->k,2);
  size_t ustart, uend, vstart, vend;
  gint m, n;

  xv.x = xv.y = xv.z = 0.;

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
    
  for ( m = 0; m < sp->k; m++) {
    gdouble cu = gsl_matrix_get (Bu, m, 0);
    gdouble cdu = gsl_matrix_get (Bu, m, 1);
    gint ii = (ustart+m);
    for ( n = 0; n < sp->k; n++) {
      gdouble cv = gsl_matrix_get (Bv, n, 0);
      gdouble cudv = cu*gsl_matrix_get (Bv, n, 1);
      gdouble cvdu = cv*cdu;
      gdouble cuv = cu*cv;
      gint jj = (vstart+n);
      xv.x += coeff (sp,ii,jj, 0)*cudv;
      xv.y += coeff (sp,ii,jj, 1)*cudv;
      xv.z += coeff (sp,ii,jj, 2)*cudv;
    }
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  return xv;
}

static Vector rotational_mode (S2P * s2p, gint k, gdouble u, gdouble v)
{
  Vector OM;
  OM.x = OM.y = OM.z = 0.;
  gint m, n;
  gint M = 2*k-2;
  
  for ( m = 0; m <= 2*(M-1); m++) {
    for ( n = 0; n <= m; n++) {
      OM = vector_sum (OM, vector_times_constant (om(n, m-n), pow(u-s2p->ue,n)*pow(v-s2p->ve,m-n)));
    }
  }

  return OM;
}

static Vector spline2d_rotational_mode (SPPanel * spp, gdouble u, gdouble v)
{
  Vector N = spline2d_normal (spp->sp, u, v);
  Vector X;

  g_assert ( u >= spp->u0 && u <= spp->u1);
  g_assert ( v >= spp->v0 && v <= spp->v1);
  
  X.x = spline2d_eval (spp->sp, u, v, 0);
  X.y = spline2d_eval (spp->sp, u, v, 1);
  X.z = spline2d_eval (spp->sp, u, v, 2);

  return vector_vector_product (&X, &N);
}

/* BE CAREFUL THIS ASSUMES THAT THE EXPENSION POINT IS THE CENTROID OF THE SURFACE */
/* void sppanel_store_jacobian_normal (SPPanel * spp, S2P * s2p) */
/* { */
/*   gint M = 2*spp->k-2;  // I think the definition of M in 5.5.1 is not consistent with its definition in the */
/*                        // previous sections */
/*   gint m, n, mu, nu; */

/*   // Copy of the polynomial coefficient of the surface representation over the panel */
/*   Vector x[2*spp->k-1][2*spp->k-1]; */
/*   Vector V0; */
/*   V0.x = V0.y = V0.z = 0.; */
 
/*   // Shift of the origin of the polynomial expansion if necessary */
/*   if (spp->ue != s2p->ue || spp->ve != s2p->ve) { */
/*     GArray * shifted = s2p_change_polynomial_origin (spp->s2pe, spp->k, spp->ue, spp->ve, s2p->ue, s2p->ve); */
/*     for ( m = 0; m < 2*spp->k-1; m++) */
/*       for ( n = 0; n < 2*spp->k-1; n++) { */
/* 	if ( m < spp->k && n < spp->k) */
/* 	  x[m][n] = g_array_index (shifted, Vector, (m) + (n)*spp->k); */
/* 	else */
/* 	  x[m][n] = V0; */
/*       } */
/*   } */
/*   else { */
/*     for ( m = 0; m < 2*spp->k-1; m++) */
/*       for ( n = 0; n < 2*spp->k-1; n++) { */
/* 	if ( m < spp->k && n < spp->k) */
/* 	  x[m][n] = g_array_index (spp->s2pe->x, Vector, (m) + (n)*spp->k); */
/* 	else */
/* 	  x[m][n] = V0; */
/*       } */
/*   } */

/*   /\* Point pe; *\/ */
/*   /\* for ( m = 0; m < spp->k; m++) *\/ */
/*   /\*   for ( n = 0; n < spp->k; n++) { *\/ */
/*   /\*     Vector xx = g_array_index (shifted, Vector, (m) + (n)*spp->k); *\/ */
/*   /\*     pe.x += xx.x*pow(spp->ue-s2p->ue, m)*pow(spp->ve-s2p->ve, n); *\/ */
/*   /\*     pe.y += xx.y*pow(spp->ue-s2p->ue, m)*pow(spp->ve-s2p->ve, n); *\/ */
/*   /\*     pe.z += xx.z*pow(spp->ue-s2p->ue, m)*pow(spp->ve-s2p->ve, n); *\/ */
/*   /\*   } *\/ */
/*   /\* Point pp = spline2d_eval_point (spp->sp, spp->ue, spp->ve); *\/ */
/*   /\* fprintf(stdout, "%e %e %e \n", pe.x - pp.x, pe.y - pp.y, pe.z - pp.z); *\/ */

/*   Vector W[2*M][2*M-1]; */
/*   for ( m = 0; m < 2*M-1; m++) { */
/*     for ( n = 0; n < 2*M; n++) { */
/*       W[n][m-n].x = W[n][m-n].y = W[n][m-n].z = 0.; */
      
/*       for ( mu = max(0,m-M+1); mu <= min(m, M-1); mu++) { */
/* 	for ( nu = max(0,n-m+mu); nu <= min(n,mu); nu++) { */
/* 	  W[n][m-n] = vector_sum (W[n][m-n], vector_times_constant ( vector_vector_product (&x[nu+1][mu-nu], &x[n-nu][m-mu-n+nu+1]) , (nu+1.)*(m-mu-n+nu+1.))); */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   // Jacobian - p 109 */
/*   gdouble w2[2*M+2][2*M+2]; */
/*   for ( m = 0; m <= 2*(M-1); m++) { */
/*     for ( n = 0; n <= m; n++) { */
/*       w2[n][m-n] = 0.; */
/*       for (mu = max(0,m-2*M+2); mu <= min(m,2*M-2); mu++) { */
/*   	for (nu = max(0,n-m+mu); nu <= min(n,mu); nu++) { */
/*   	  w2[n][m-n] += vector_scalar_product (&W[nu][mu-nu], &W[n-nu][m-mu-n+nu]); */
/*   	} */
/*       } */
/*     } */
/*   } */

/*   // page 110 */
/*   w(0,0) = sqrt(w2[0][0]); */
  
/*   for ( m = 1; m <= 2*(M-1); m++ ) { */
/*     for ( n = 0; n <= m; n++ ) { */

/*   	if ( m <= 4*(M-1) ) */
/*   	  w(n,m-n) = w2[n][m-n]; */
/*   	else */
/*   	  w(n,m-n) = 0.; */
      
/*   	for ( mu = 1; mu <= m-1; mu++ ) { */
/*   	  for ( nu = max(0, n-m+mu); nu <= min(n, mu); nu++ ) { */
/*   	    w(n,m-n) -= w(nu,mu-nu)*w(n-nu,m-mu-n+nu); */
/*   	  } */
/*   	} */
	
/*   	w(n,m-n) /= (2.*w(0,0)); */
/*     } */
/*   } */

/*   // Normal - p 112 */
/*   gdouble wm1[2*M+2][2*M+2]; */
/*   wm1[0][0] = 1./w(0,0); */
/*   for ( m = 1; m <= 2*(M-1); m++ ) { */
/*     for ( n = 0; n <= m; n++ ) { */

/*   	wm1[n][m-n] = 0.; */

/*   	for ( mu = 1; mu <= m; mu++ ) { */
/*   	  for ( nu = max(0, n-m+mu); nu <= min(n, mu); nu++ ) { */
/* 	    wm1[n][m-n] += w(nu,mu-nu)*wm1[n-nu][m-mu-n+nu]; // Corrected from Maniar, 1995 */
/*   	  } */
/*   	} */
	
/*   	wm1[n][m-n] *= -wm1[0][0]; */
/*     } */
/*   } */

/*   // Normal vector - p 112 */
/*   for ( m = 0; m <= 2*(M-1); m++ ) { */
/*     for ( n = 0; n <= m; n++ ) { */

/*       for ( mu = 0; mu <= min(m,2*(M-1)); mu++) { */
/*   	  for ( nu = max(0,n-m+mu); nu <= min(n,mu); nu++) { */
/* 	    n(n,m-n) = vector_sum (n(n,m-n), */
/* 				   vector_times_constant(W[nu][mu-nu], */
/* 							 wm1[n-nu][m-mu-n+nu])); */
/*   	  } */
/*   	} */
/*     } */
/*   } */

/*   /\* Vector N = normal (s2p, spp->k, s2p->ue, s2p->ve); *\/ */
/*   /\* Vector N2 = spline2d_normal (spp, s2p->ue, s2p->ve); *\/ */
/*   /\* fprintf(stdout, "%e %e %e \n", N.x - N2.x, N.y - N2.y, N.z - N2.z); *\/ */



/*   // Rotational modes - p 113 */
/*   for ( m = 0; m <= 2*(M-1); m++ ) { */
/*     for ( n = 0; n <= m; n++ ) { */

/*       for ( mu = 0; mu <= min(m,M); mu++) { */
/*   	for ( nu = max(0, n-m+mu); nu <= min(n,mu); nu++) { */
/* 	  om(n,m-n) = vector_sum (om(n,m-n), vector_vector_product (&x[nu][mu-nu], &n(n-nu,m-mu-n+nu))); */
/*   	} */
/*       } */
/*     } */
/*   } */
/* } */


/* GaussPoints * spline2d_store_gauss_legendre_points (Spline2D * sp, gdouble u0, gdouble u1, gdouble v0, gdouble v1, gint ng) */
/* { */
/*   gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng); */
/*   gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng); */
/*   gdouble ui, vj, wi, wj; */
/*   gint i, j; */
/*   GaussPoints * gp = gausspoints_new (ng); */
  
/*   for ( i = 0; i < ng; i++) { */
/*     gsl_integration_glfixed_point (u0, u1, i, &ui, &wi, itable); */
/*     g_array_index (gp->ui, gdouble, i) = ui; */
/*     gsl_integration_glfixed_point (v0, v1, i, &vj, &wj, jtable); */
/*     g_array_index (gp->vj, gdouble, i) = vj; */

/*     size_t ustart, uend, vstart, vend; */
/*     gsl_matrix * Bu = gsl_matrix_alloc (sp->k, 2); */
/*     gsl_bspline_deriv_eval_nonzero (ui, 1, Bu, &ustart, &uend, sp->w_u, sp->wd_u); */
/*     g_ptr_array_add (gp->Bu, Bu); */

/*     gsl_matrix * Bv = gsl_matrix_alloc (sp->k, 2); */
/*     gsl_bspline_deriv_eval_nonzero (vj, 1, Bv, &vstart, &vend, sp->w_v, sp->wd_u); */
/*     g_ptr_array_add (gp->Bv, Bv); */

/*     if ( i == 0) { */
/*       gp->istart = ustart; */
/*       gp->jstart = vstart; */
/*     } */
/*   } */

/*   // Stores the normal, the physical coordinates and the jacobian at the inner Gauss points */
/*   for ( j = 0; j < ng; j++) { */
/*     gsl_integration_glfixed_point (v0, v1, j, &vj, &wj, jtable); */
/*     for ( i = 0; i < ng; i++) { */
/*       gsl_integration_glfixed_point (u0, u1, i, &ui, &wi, itable); */
      
/*       Vector N = spline2d_normal (sp, ui(i), vj(j)); */
/*       Point P = spline2d_eval_point (sp, ui(i), vj(j)); */

/*       g_array_append_val (gp->Ni, N); */
/*       g_array_append_val (gp->Pi, P); */

/*       gdouble J = wi*wj; */
/*       g_array_append_val (gp->wij, J); */
/*       J *= spline2d_jacobian (sp, ui(i), vj(j)); */
/*       g_array_append_val (gp->wJij, J); */
/*     } */
/*   } */
  
/*   gsl_integration_glfixed_table_free (itable); */
/*   gsl_integration_glfixed_table_free (jtable); */
/*   return gp; */
/* } */

static GrevillePoints * spline2d_store_greville_points (Spline2D * sp)
{
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint i, j;
  GrevillePoints * gr = greville_points_new (NU, NV);
 
  size_t ustart, uend, vstart, vend;
  for ( i = 0; i < NU; i++) {
    gdouble u = g_array_index (gr->ui, gdouble, i) = gsl_bspline_greville_abscissa (i, sp->w_u);
    gsl_matrix * Bu = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->w_u, sp->wd_u);
    g_ptr_array_add (gr->Bu, Bu);
    g_array_append_val (gr->ustart, ustart);
  }
    
  for ( i = 0; i < NV; i++) {
    gdouble v = g_array_index (gr->vj, gdouble, i) = gsl_bspline_greville_abscissa (i, sp->w_v);
    gsl_matrix * Bv = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
    g_ptr_array_add (gr->Bv, Bv);
    g_array_append_val (gr->vstart, vstart);
  }

  // Stores the normal, the physical coordinates and the jacobian at the inner Gauss points
  for ( j = 0; j < NV; j++) {
    gdouble v = g_array_index (gr->vj, gdouble, j);
    for ( i = 0; i < NU; i++) {
      gdouble u = g_array_index (gr->ui, gdouble, i);
      
      Vector N = spline2d_normal (sp, u, v);
      Point P = spline2d_eval_point (sp, u, v);
      gdouble J = spline2d_jacobian (sp, u, v);
      Vector xu = spline2d_xu (sp, u, v);
      Vector xv = spline2d_xv (sp, u, v);

      g_array_append_val (gr->Ni, N);
      g_array_append_val (gr->Pi, P);
      g_array_append_val (gr->Ji, J);
    }
  }

  return gr;
}

void spline2d_store_greville_points_data (Spline2D * sp)
{
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint i, j;

  // Stores the normal, the physical coordinates and the jacobian
  for ( j = 0; j < NV; j++) {
    gdouble v = g_array_index (gr->vj, gdouble, j);
    for ( i = 0; i < NU; i++) {
      gdouble u = g_array_index (gr->ui, gdouble, i);

      g_array_index (gr->Ni, Vector, i + NU*j) =  spline2d_normal (sp, u, v);
      g_array_index (gr->Pi, Point, i + NU*j) =   spline2d_eval_point (sp, u, v);
      g_array_index (gr->Ji, gdouble, i + NU*j) = spline2d_jacobian (sp, u, v);
    }
  }
}

/**
 * Compute and stores Gauss-Legendre points and weights for a given panel
 **/
GaussPoints * sppanel_store_gauss_legendre_points (SPPanel * spp, gint ng)
{
  Spline2D * sp = spp->sp;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui, vj, wi, wj;
  gint i, j;
  GaussPoints * gp = gausspoints_new (ng);
  
  for ( i = 0; i < ng; i++) {
    gsl_integration_glfixed_point (spp->u0, spp->u1, i, &ui, &wi, itable);
    g_array_index (gp->ui, gdouble, i) = ui;
    gsl_integration_glfixed_point (spp->v0, spp->v1, i, &vj, &wj, jtable);
    g_array_index (gp->vj, gdouble, i) = vj;

    size_t ustart, uend, vstart, vend;
    gsl_matrix * Bu = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,ui)), 1, Bu, &ustart, &uend, sp->w_u, sp->wd_u);
    g_ptr_array_add (gp->Bu, Bu);

    gsl_matrix * Bv = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,vj)), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
    g_ptr_array_add (gp->Bv, Bv);
    if ( i == 0) {
      gp->istart = ustart;
      gp->jstart = vstart;
    }

    // Geometry splines
    Bu = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,ui)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
    g_ptr_array_add (gp->Bux, Bu);

    if ( i == 0)
      gp->istart_x = ustart;
  }

  // Stores the normal, the physical coordinates and the jacobian at the inner Gauss points
  for ( j = 0; j < ng; j++) {
    gsl_integration_glfixed_point (spp->v0, spp->v1, j, &vj, &wj, jtable);
    for ( i = 0; i < ng; i++) {
      gsl_integration_glfixed_point (spp->u0, spp->u1, i, &ui, &wi, itable);
      
      Vector N = spline2d_normal (sp, ui(i), vj(j));
      Point P = spline2d_eval_point (sp, ui(i), vj(j));
      Vector xu = spline2d_xu (sp, ui(i), vj(j));
      Vector xv = spline2d_xv (sp, ui(i), vj(j));

      gdouble det = xu.x*(xv.y*N.z-xv.z*N.y)
	- xu.y*(N.z*xv.x-xv.z*N.x)
	+ xu.z*(xv.x*N.y-xv.y*N.x);
      det = 1./det;
      gdouble c1 = (xv.y*N.z-xv.z*N.y)*det;
      gdouble c2 = (xu.z*N.y-xu.y*N.z)*det;
      gdouble c3 = (xv.z*N.x-xv.x*N.z)*det;
      gdouble c4 = (xu.x*N.z-xu.z*N.x)*det;
      gdouble c5 = (xv.x*N.y-xv.y*N.x)*det;
      gdouble c6 = (N.x*xu.y-xu.x*N.y)*det;

      g_array_append_val (gp->Ni, N);
      g_array_append_val (gp->Pi, P);
      g_array_append_val (gp->c1, c1);
      g_array_append_val (gp->c2, c2);
      g_array_append_val (gp->c3, c3);
      g_array_append_val (gp->c4, c4);
      g_array_append_val (gp->c5, c5);
      g_array_append_val (gp->c6, c6);

      gdouble J = wi*wj;
      g_array_append_val (gp->wij, J);
      J *= spline2d_jacobian (sp, ui(i), vj(j));
      g_array_append_val (gp->wJij, J);
    }
  }
  
  gsl_integration_glfixed_table_free (itable);
  gsl_integration_glfixed_table_free (jtable);
  return gp;
}

/**
 * Compute and stores Gauss-Legendre points and weights for a given panel
 **/
void sppanel_store_gauss_legendre_data (GaussPoints * gp, SPPanel * spp)
{
  gint ng = gp->ui->len;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui, vj, wi, wj;
  gint i, j;

  // Stores the normal, the physical coordinates and the jacobian at the inner Gauss points
  for ( j = 0; j < ng; j++) {
    gsl_integration_glfixed_point (spp->v0, spp->v1, j, &vj, &wj, jtable);
    for ( i = 0; i < ng; i++) {
      gsl_integration_glfixed_point (spp->u0, spp->u1, i, &ui, &wi, itable);
      g_assert ( ui(i) <= 1. && ui(i) >= -1.);
      Vector N = spline2d_normal (spp->sp, ui, vj);
      g_array_index (gp->Ni, Vector, i + ng*j) = N;
      g_array_index (gp->Pi, Point, i + ng*j) = spline2d_eval_point (spp->sp, ui, vj);
      g_array_index (gp->wJij, gdouble, i + ng*j) = spline2d_jacobian (spp->sp, ui, vj)*wi*wj;
      Vector xu = spline2d_xu (spp->sp, ui, vj);
      Vector xv = spline2d_xv (spp->sp, ui, vj);
      gdouble det = xu.x*(xv.y*N.z-xv.z*N.y)
	- xu.y*(N.z*xv.x-xv.z*N.x)
	+ xu.z*(xv.x*N.y-xv.y*N.x);
      det = 1./det;
      g_array_index (gp->c1, gdouble, i + ng*j) = (xv.y*N.z-xv.z*N.y)*det;
      g_array_index (gp->c2, gdouble, i + ng*j) = (xu.z*N.y-xu.y*N.z)*det;
      g_array_index (gp->c3, gdouble, i + ng*j) = (xv.z*N.x-xv.x*N.z)*det;
      g_array_index (gp->c4, gdouble, i + ng*j) = (xu.x*N.z-xu.z*N.x)*det;
      g_array_index (gp->c5, gdouble, i + ng*j) = (xv.x*N.y-xv.y*N.x)*det;
      g_array_index (gp->c6, gdouble, i + ng*j) = (N.x*xu.y-xu.x*N.y)*det;
    }
  }
  
  gsl_integration_glfixed_table_free (itable);
  gsl_integration_glfixed_table_free (jtable);
}

/* /\** */
/*  *  Initialise the panels of the spline patch and store */
/*  *  all the quantities that can be precomputed */
/*  **\/ */
/* void spline2d_init_panels (Spline2D * sp) */
/* { */
/*   gint i, j; */

/*   // Loop over all the panels */
/*   for (i = 0; i < sp->M; i++) { */
/*     for (j = 0; j < sp->N; j++) { */
/*       SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M); */
/*       g_assert (spp != NULL); */

/*       // Link back to the spline patch */
/*       spp->sp = sp; */
/*       spp->k = sp->k; */
      
/*       // Panel parametric boundaries */
/*       spp->u0 = gsl_vector_get (sp->w_u->knots,sp->k + i - 1); */
/*       spp->u1 = gsl_vector_get (sp->w_u->knots,sp->k + i); */
/*       spp->v0 = gsl_vector_get (sp->w_v->knots,sp->k + j - 1); */
/*       spp->v1 = gsl_vector_get (sp->w_v->knots,sp->k + j); */
      
/*       // Panel parametric centroid */
/*       spp->ue = (spp->u0 + spp->u1)/2.; */
/*       spp->ve = (spp->v0 + spp->v1)/2.; */
/*       spp->pe = spline2d_eval_point (sp, spp->ue, spp->ve); */

/*       // Compute and store Gauss-Legendre points and weights inner and outer integration */
/*       spp->inner = sppanel_store_gauss_legendre_points (spp, sp->ninner); */
/*       spp->outer = sppanel_store_gauss_legendre_points (spp, sp->nouter); */
/*     } */
/*   } */
/* } */

/**
 *  Initialise the panels of the spline patch and store
 *  all the quantities that can be precomputed
 **/
void spline2d_init_panels (Spline2D * sp)
{
  gint i, j;

  // Loop over all the panels
  for (i = 0; i < sp->M; i++) {
    for (j = 0; j < sp->N; j++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
      g_assert (spp != NULL);

      // Link back to the spline patch
      spp->sp = sp;
      spp->k = sp->k;
      
      // Panel parametric boundaries
      spp->u0 = gsl_vector_get (sp->wx_u->knots, sp->k + i - 1);
      spp->u1 = gsl_vector_get (sp->wx_u->knots, sp->k + i);
      spp->v0 = gsl_vector_get (sp->w_v->knots, sp->k + j - 1);
      spp->v1 = gsl_vector_get (sp->w_v->knots, sp->k + j);
      
      // Panel parametric centroid
      spp->ue = (spp->u0 + spp->u1)/2.;
      spp->ve = (spp->v0 + spp->v1)/2.;
      spp->pe = spline2d_eval_point (sp, spp->ue, spp->ve);

      // Compute and store a polynomial representation of the surface over the panel
      // expended using the centroid of the panel
      /* g_assert (spp->s2pe == NULL); */
      /* spp->s2pe = sppanel_spline_to_poly_matrix (spp, sp, spp->ue, spp->ve); */

      // Store the coefficients of the polynomial expansion of the jacobian, the normal
      // and the rotational modes
      /* sppanel_store_jacobian_normal (spp, spp->s2pe); */

      // Compute and store Gauss-Legendre points and weights inner and outer integration
      spp->outer = sppanel_store_gauss_legendre_points (spp, sp->nouter);

      // Compute and store a polynomial representation of the surface over the panel
      // expended using the each outer Gauss point
      /* gint m, n; */
      /* GaussPoints * gp = spp->outer; */
      /* for ( n = 0; n < sp->nouter; n++) { */
      /* 	for ( m = 0; m < sp->nouter; m++) { */
      /* 	  g_ptr_array_add (spp->s2p, sppanel_spline_to_poly_matrix (spp, sp, ui(m), vj(n))); */
      /* 	  sppanel_store_jacobian_normal (spp, g_ptr_array_index (spp->s2p, m + n*sp->nouter)); */
      /* 	} */
      /* } */

    }
  }

  // Greville Points for collocation approach
  sp->gr = spline2d_store_greville_points (sp);

  g_warning ("Normal orientation needs checking manually before starting any serious run\n");
}

void spline2d_reinit_panels_physical_quantities (Spline2D * sp)
{
  gint i, j;

  // Loop over all the panels
  for (i = 0; i < sp->M; i++) {
    for (j = 0; j < sp->N; j++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
      g_assert (spp != NULL);
      
      // Panel parametric centroid
      spp->pe = spline2d_eval_point (sp, spp->ue, spp->ve);

      // Compute and store a polynomial representation of the surface over the panel
      // expended using the centroid of the panel
      /* s2p_destroy (spp->s2pe); */
      /* spp->s2pe = sppanel_spline_to_poly_matrix (spp, sp, spp->ue, spp->ve); */

      // Store the coefficients of the polynomial expansion of the jacobian, the normal
      // and the rotational modes
      /* sppanel_store_jacobian_normal (spp, spp->s2pe); */

      // Compute and store physical quantities at Gauss-Legendre points
      sppanel_store_gauss_legendre_data (spp->outer, spp);

      // Compute and store a polynomial representation of the surface over the panel
      // expended using the each outer Gauss point
      /* gint m, n; */
      /* GaussPoints * gp = spp->outer; */
      /* for ( n = 0; n < sp->nouter; n++) { */
      /* 	for ( m = 0; m < sp->nouter; m++) { */
      /* 	  s2p_destroy (g_ptr_array_index (spp->s2p, m + n*sp->nouter)); */
      /* 	  g_ptr_array_index (spp->s2p, m + n*sp->nouter) = */
      /* 	    sppanel_spline_to_poly_matrix (spp, sp, ui(m), vj(n)); */
      /* 	  sppanel_store_jacobian_normal (spp, g_ptr_array_index (spp->s2p, m + n*sp->nouter)); */
      /* 	} */
      /* } */

    }
  }

  spline2d_store_greville_points_data (sp);
  g_warning ("Normal orientation needs checking manually before starting any serious run\n");
}

void sppanel_print (SPPanel * spp, FILE * fp)
{
  Spline2D * sp = spp->sp;

  Point p = spline2d_eval_point (sp, spp->u0, spp->v0);
  fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
  p = spline2d_eval_point (sp, spp->u1, spp->v0);
  fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
  p = spline2d_eval_point (sp, spp->u1, spp->v1);
  fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
  p = spline2d_eval_point (sp, spp->u0, spp->v1);
  fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
  p = spline2d_eval_point (sp, spp->u0, spp->v0);
  fprintf(fp, "%f %f %f \n\n\n", p.x, p.y, p.z);
}

void spline2d_print_panels (Spline2D * splines, FILE * fp)
{
  g_assert (fp != NULL);
  gint i, j;
  
  Spline2D * sp = splines;
  while (sp) {
    fprintf(fp, "#New PAtch \n");
    for ( i = 0; i < sp->M; i++) {
      for ( j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
	sppanel_print (spp, fp);
      }
    }
    sp = sp->next;
  }
}

void spline2d_print_transformed_panels (Spline2D * splines, FILE * fp,
					Point * xg, Vector * t,
					Matrix3 * euler_m)
{
  g_assert (fp != NULL);

  gint i, j;
  Spline2D * sp = splines;
  while (sp) {
    for ( i = 0; i < sp->M; i++) {
      for ( j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
	Point p = transform_point (spline2d_eval_point (sp, spp->u0, spp->v0),
				   xg, t, euler_m);
	fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
	p = transform_point (spline2d_eval_point (sp, spp->u1, spp->v0),
			     xg, t, euler_m);
	fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
	p = transform_point (spline2d_eval_point (sp, spp->u1, spp->v1),
			     xg, t, euler_m);
	fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
	p = transform_point (spline2d_eval_point (sp, spp->u0, spp->v1),
			     xg, t, euler_m);
	fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
	p = transform_point (spline2d_eval_point (sp, spp->u0, spp->v0),
			     xg, t, euler_m);
	fprintf(fp, "%f %f %f \n\n\n", p.x, p.y, p.z);
	
      }
    }
    sp = sp->next;
  }
}

/**
 * Plot using:
   set style data pm3d
   set palette rgbformulae 22,13,-31
   splot 'hull-color.out'
 **/
void spline2d_print_panels_gnuplot (Spline2D * sp, FILE * fp, gint var)
{
  g_assert (fp != NULL);

  gdouble u, v;
  gdouble du = 0.00999999999;

  for ( u = 0; u+du < 1.; u+=du) {
    for ( v = 0; v+du < 1.; v+=du) {
      Point p = spline2d_eval_point (sp, u, v);
      fprintf(fp, "%f %f %f %f\n", p.x, p.y, p.z,
	      spline2d_eval (sp, u, v, var));
      p = spline2d_eval_point (sp, u+du, v);
      fprintf(fp, "%f %f %f %f \n\n", p.x, p.y, p.z,
	      spline2d_eval (sp, u+du, v, var));
      
      p = spline2d_eval_point (sp, u, v+du);
      fprintf(fp, "%f %f %f %f \n", p.x, p.y, p.z,
	      spline2d_eval (sp, u, v+du, var));
      p = spline2d_eval_point (sp, u+du, v+du);
      fprintf(fp, "%f %f %f %f \n\n\n", p.x, p.y, p.z,
	      spline2d_eval (sp, u+du, v+du, var));
    }
  }

}

void spline2d_print_normals (Spline2D * splines, FILE * fp)
{
  g_assert (fp != NULL);

  gint i, j;
  Spline2D * sp = splines;
  while (sp) {
    for ( i = 0; i < sp->M; i++) {
      for ( j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
	Point p = spline2d_eval_point (sp, spp->ue, spp->ve);
	Vector n = spline2d_normal (sp, spp->ue, spp->ve);

	fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
	fprintf(fp, "%f %f %f \n\n\n", p.x + n.x, p.y + n.y, p.z + n.z);
      }
    }
    sp = sp->next;
  }
}

void spline2d_print_transformed_normals (Spline2D * sp, FILE * fp,
					 Point * xg, Vector * t,
					 Matrix3 * euler_m)
{
  g_assert (fp != NULL);

  gint i, j;

  for ( i = 0; i < sp->M; i++) {
    for ( j = 0; j < sp->N; j++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
      //    Point p = spline2d_eval_point (sp, spp->ue, spp->ve);
      Point p = transform_point (spline2d_eval_point (sp, spp->ue, spp->ve),
				 xg, t, euler_m);
      Vector n = transform_vector (spline2d_normal (sp, spp->ue, spp->ve), euler_m);

      fprintf(fp, "%f %f %f \n", p.x, p.y, p.z);
      fprintf(fp, "%f %f %f \n\n\n", p.x + n.x, p.y + n.y, p.z + n.z);
    }
  }
}

gdouble spline2d_gauss_integration (Spline2D * sp, GaussFunc func, gpointer data)
{
  gdouble sum = 0.;
  gint ii, m, n;

  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    GaussPoints * gp = spp->outer; // Store Gauss-point and weights
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule

    for ( m = 0; m < ng; m++) {	  
      for ( n = 0; n < ng; n++) {
	sum += g_array_index (gp->wJij, gdouble, m+n*ng) * func(spp, m, n, data);
      }
    }
  }

  return sum;
}

static gdouble adaptive_gauss_integration (SPPanel * spp, 
					   gdouble u1, gdouble u2,
					   gdouble v1, gdouble v2,
					   SPPanelFunc func, gpointer data)
{
  gint ng = spp->outer->ui->len;
  gsl_integration_glfixed_table * mtable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * ntable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble um, vn, wm, wn, sum = 0.;
  gint m, n;

  for ( m = 0; m < ng; m++) {
    gsl_integration_glfixed_point (u1, u2, m, &um, &wm, mtable);
    for ( n = 0; n < ng; n++) {
      gsl_integration_glfixed_point (v1, v2, n, &vn, &wn, ntable);
      
      sum += wm * wn * spline2d_jacobian (spp->sp, um, vn) * func (spp, um, vn, data);
    }
  }

  gsl_integration_glfixed_table_free (mtable);
  gsl_integration_glfixed_table_free (ntable);

  return sum;
}

gdouble spline2d_adaptive_gauss_integration (Spline2D * sp, SPPanelFunc func, gpointer data, gint level)
{
  gint i, j, ii;
  gdouble sum = 0.;

  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    gdouble du = (spp->u1-spp->u0)/pow(2,level);
    gdouble dv = (spp->v1-spp->v0)/pow(2,level);

    for ( i = 0; i < pow(2,level); i++) {
      for ( j = 0; j < pow(2,level); j++) {
	sum += adaptive_gauss_integration (spp,
					   spp->u0 + i*du, /* min( */spp->u0 + (i+1.)*du/* , 1.) */,
					   spp->v0 + j*dv, /* min( */spp->v0 + (j+1.)*dv/* , 1.) */,
					   func, data);
      }
    }
  }
  
  return sum;
}


/************** Self-influence coefficients ************************/

InfluenceCoeffs * influencecoeffs_new (gint k)
{
  InfluenceCoeffs * new = g_malloc (sizeof(InfluenceCoeffs));
  new->psi = gsl_matrix_alloc (k, k);
  gsl_matrix_set_zero (new->psi);
  new->phi = gsl_matrix_alloc (k, k);
  gsl_matrix_set_zero (new->phi);
  return new;
}

void influencecoeffs_destroy (InfluenceCoeffs * ic)
{
  g_assert (ic != NULL);
  g_assert (ic->phi != NULL);
  g_assert (ic->psi != NULL);
  gsl_matrix_free (ic->phi);
  gsl_matrix_free (ic->psi);
  g_free (ic);
}

void influencecoeffs_add (InfluenceCoeffs * ic1, InfluenceCoeffs * ic2)
{
  gsl_matrix_add (ic1->psi, ic2->psi);
  gsl_matrix_add (ic1->phi, ic2->phi);
}

InfluenceCoeffs * influencecoeffs_add_and_destroy (InfluenceCoeffs * ic1,
						   InfluenceCoeffs * ic2)
{
  gsl_matrix_add (ic1->psi, ic2->psi);
  gsl_matrix_add (ic1->phi, ic2->phi);
  influencecoeffs_destroy (ic2);
  return ic1;
}

/* void spline2d_clear_self_influence_coefficients (Spline2D * sp) */
/* { */
/*   gint i; */

/*   g_assert (sp != NULL); */

/*   for (i = 0; i < sp->panels->len; i++) { */
/*     SPPanel * panel = g_ptr_array_index (sp->panels, i); */
    
/*     if (panel->sic == NULL) */
/*       return; */

/*     gint j; */
/*     for ( j = 0; j < panel->sic->len; j++) */
/*       influencecoeffs_destroy (g_ptr_array_index (panel->sic, j)); */

/*     g_ptr_array_free (panel->sic, TRUE); */
/*     panel->sic = NULL; */
/*   } */
/* } */

typedef struct {
  GArray * q;
  gdouble a, b, c;
  GArray * s;
  GArray * p, * g;
} SIData;

static void self_influence_intermediate_calculations  (SIData * sid, gint k, gdouble xie, gdouble eta)
{
  // See section 5.1 of (Maniar,1995)
  gint M = 2*k-1;
  gint N = M-1;
  gint m,n, mu, nu;
  /* gdouble val; */
  /* gdouble delta = 4*sid->a*sid->c-sid->b*sid->b; */

  // page 84
  gdouble q2[2*N+1][2*N+3];
  for ( m = 0; m < 2*N+1; m++) {
    for ( n = 0; n <= m+2; n++) {
      q2[m][n] = 0.;
      for ( mu = max(0,m-N); mu <= min(m, N); mu++) {
  	for ( nu = max(0,n-m+mu-1); nu <= min(n,mu+1); nu++) {
  	  q2[m][n] += vector_scalar_product(&q(mu,nu), &q(m-mu,n-nu));
  	}
      }
    }
  }

  // page 99 + 97
  s(0,0) = 1.;
  s(1,0) = sid->a;
  s(1,1) = sid->b;
  s(1,2) = sid->c;

  for ( m = 2; m < 8*N; m++) {
    for ( n = 0; n <= 2*m; n++) {
      s(m,n) = 0.;
      /* gdouble sum = 0.; */
      for ( nu = max(0, n-2); nu <= min(n, 2*m-2); nu++) {
  	s(m,n) += s(m-1,nu)*s(1,n-nu);
      }
    }
  }

  // page 85
  gdouble h2[2*N+1][6*N+1];
  h2[0][0] = 1.;
  for ( m = 1; m < 2*N+1; m++) {
    for ( n = 0; n <= 3*m; n++) {
      h2[m][n] = 0.;
      for ( nu = max(0,n-2*m+2); nu <= min(n,m+2); nu++) {
  	h2[m][n] += q2[m][nu]*s(m-1,n-nu);
      }
    }
  }


  // page 86
  gdouble hm2[5*N][15*N+3];
  hm2[0][0] = 1.;
  for ( m = 1; m < 5*N; m++) {
    for ( n = 0; n <= 3*m; n++) {
      hm2[m][n] = 0.;
      for ( mu = 1; mu <= min(m,2*N); mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  hm2[m][n] -= h2[mu][nu]*hm2[m-mu][n-nu];
  	}
      }

    }
  }

  // page 87/88
  gdouble hm1[5*N][15*N+3];
  hm1[0][0] = 1.;
  for ( m = 1; m < 5*N; m++) {
    for ( n = 0; n <= 3*m; n++) {
      hm1[m][n] = 0.5*hm2[m][n];
      for ( mu = 1; mu <= m-1; mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  hm1[m][n] -= 0.5*hm1[mu][nu]*hm1[m-mu][n-nu];
  	}
      }
    }
  }
  
  // page 88
  gdouble hm3[3*N+1][9*N+6];
  hm3[0][0] = 1.;
  for ( m = 1; m < 3*N+1; m++) {
    for ( n = 0; n <= 3*m; n++) {
      hm3[m][n] = 0.;
      for ( mu = 0; mu <= m; mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  hm3[m][n] += hm2[mu][nu]*hm1[m-mu][n-nu];
  	}
      }
    }
  }

  // page 92
  Vector w[2*N+1][2*N+1];
  for ( m = 0; m < 2*N+1; m++) {
    for ( n = 0; n <= m; n++) {
      w[m][n].x = w[m][n].y = w[m][n].z = 0.;
      
      for ( mu = max(0, m-N); mu <= min(m,N); mu++) {
  	for ( nu = max(0,n-m+mu); nu <= min(n,mu+1); nu++) {
  	  w[m][n] = vector_sum (w[m][n],vector_times_constant( vector_vector_product(&q(mu,nu),&q(m-mu,n-nu+1)),
  	  						       (mu+1)*(n-nu+1)));
  	}
      }

    }
  }

  // page 92
  gdouble w2[4*N+1][4*N+1];
  for ( m = 0; m < 4*N+1; m++) {
    for ( n = 0; n <= m; n++) {
      w2[m][n] = 0.;

      for ( mu = max(0,m-2*N); mu <= min(m,2*N); mu++) {
  	for ( nu = max(0, n-m+mu); nu <= min(n,mu); nu++) {
  	  w2[m][n] += vector_scalar_product (&w[mu][nu],&w[m-mu][n-nu]);
  	}
      }

    }
  }

  // page 92/92
  /* gdouble w1[8*N+1][8*N+3]; */
  gdouble w1[4*N+2][4*N+3];
  w1[0][0] = sqrt (w2[0][0]);
  for ( m = 1; m < 4*N+2; m++) {
    for ( n = 0; n <= m; n++) {

  	if ( m <= 4*N )
  	  w1[m][n] = w2[m][n];
  	else
  	  w1[m][n] = 0.;

  	for ( mu = 1; mu <= m-1; mu++) {
  	  for ( nu = max(0,n-m+mu); nu <= min(n,mu); nu++) {
  	    w1[m][n] -= w1[mu][nu]*w1[m-mu][n-nu];
  	  }
  	}
	
  	w1[m][n] /= (2.*w1[0][0]); // w1 here is instead of w2. I think this is right and the manuscript is wrong
    }
  }

  // page 93
  Vector cm[2*N][2*N+2];
  for ( m = 0; m < 2*N; m++) {
    for ( n = 0; n <= m+2; n++) {
      cm[m][n].x = cm[m][n].y = cm[m][n].z = 0.;
      
      for ( mu = max(0, m-N); mu <= min(m,N-1); mu++) {
  	for ( nu = max(0,n-m+mu); nu <= min(n,mu+2); nu++) {
  	  cm[m][n] = vector_sum (cm[m][n],vector_times_constant( vector_vector_product(&q(mu+1,nu),&q(m-mu,n-nu+1)),
  								 (mu+1)*(n-nu+1)));
	  
  	}
      }

    }
  }

  // page 94
  gdouble e[3*N+3][3*N+6];
  for ( m = 0; m < 3*N+3; m++) {
    for ( n = 0; n <= m+3; n++) {
      e[m][n] = 0.;

      for ( mu = max(0,m-N); mu <= min(m,2*N-1); mu++) {
  	for ( nu = max(0,n-m+mu-1); nu <= min(n,mu+2); nu++) {
  	  e[m][n] += vector_scalar_product (&cm[mu][nu], &q(m-mu,n-nu));
  	}
      }

    }
  }

  // page 95
  gdouble t[4*N+2][12*N+12];
  for ( m = 0; m < 4*N+2; m++) {
    for ( n = 0; n <= 3*m; n++) {
      t[m][n] = 0.;

      for ( nu = max(0,n-m); nu <= min(n,2*m); nu++) {
  	t[m][n] += s(m,nu)*w1[m][n-nu];
      }

    }
  }

  // page 96
  for ( m = 0; m < 4*N+2; m++) {
    for ( n = 0; n <= 3*m ; n++) {
      p(m,n) = 0.;
      for ( mu = 0; mu <= m; mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  p(m,n) += t[mu][nu]*hm1[m-mu][n-nu];
  	}
      }

    }
  }

  // page 96
  gdouble f[3*N+1][9*(N+1)+1];
  for ( m = 0; m < 3*N+1; m++) {
    for ( n = 0; n <= 9*(N+1)+1; n++) {
      f[m][n] = 0.;

      for ( nu = max(0,n-m-3); nu <= min(n,2*m); nu++) {
  	  f[m][n] += s(m,nu)*e[m][n-nu];
      }

    }
  }

  // page 96
  for ( m = 0; m < 3*N+1; m++) {
    for ( n = 0; n <= 3*m+3; n++) {
      g(m,n) = 0.;
      for ( mu = 0; mu <= min(m,3*N-1); mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu+3); nu++) {
  	  g(m,n) += f[mu][nu]*hm3[m-mu][n-nu];
  	}
      }
      
    }
  }

}

static void self_influence_intermediate_calculations2  (SIData * sid, gint k, gdouble xie, gdouble eta)
{
  // See section 5.1 of (Maniar,1995)
  gint N = 2*k-2;
  gint m,n, mu, nu;

  // page 84
  gdouble q2[2*N+1][2*N+3];

  for ( m = 0; m < 2*N+1; m++) {
    for ( n = 0; n <= m+2; n++) {
      q2[m][n] = 0.;
      for ( mu = max(0,m-N); mu <= min(m, N); mu++) {
  	for ( nu = max(0,n-m+mu-1); nu <= min(n,mu+1); nu++) {
  	  q2[m][n] += vector_scalar_product(&q(mu,nu), &q(m-mu,n-nu));
  	}
      }
    }
  }

  // page 99 + 97
  s(0,0) = 1.;
  s(1,0) = sid->a;
  s(1,1) = sid->b;
  s(1,2) = sid->c;

  for ( m = 2; m < 8*N; m++) {
    for ( n = 0; n <= 2*m; n++) {
      s(m,n) = 0.;
      /* gdouble sum = 0.; */
      for ( nu = max(0, n-2); nu <= min(n, 2*m-2); nu++) {
  	s(m,n) += s(m-1,nu)*s(1,n-nu);
      }
    }
  }

  // page 85
  gdouble h2[2*N+1][6*N+1];
  h2[0][0] = 1.;
  for ( m = 1; m < 2*N+1; m++) {
    for ( n = 0; n <= 3*m; n++) {
      h2[m][n] = 0.;
      for ( nu = max(0,n-2*m+2); nu <= min(n,m+2); nu++) {
  	h2[m][n] += q2[m][nu]*s(m-1,n-nu);
      }
    }
  }


  // page 86
  gdouble hm2[5*N][15*N+3];
  hm2[0][0] = 1.;
  for ( m = 1; m < 5*N; m++) {
    for ( n = 0; n <= 3*m; n++) {
      hm2[m][n] = 0.;
      for ( mu = 1; mu <= min(m,2*N); mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  hm2[m][n] -= h2[mu][nu]*hm2[m-mu][n-nu];
  	}
      }

    }
  }

  // page 87/88
  gdouble hm1[5*N][15*N+3];
  hm1[0][0] = 1.;
  for ( m = 1; m < 5*N; m++) {
    for ( n = 0; n <= 3*m; n++) {
      hm1[m][n] = 0.5*hm2[m][n];
      for ( mu = 1; mu <= m-1; mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  hm1[m][n] -= 0.5*hm1[mu][nu]*hm1[m-mu][n-nu];
  	}
      }
    }
  }
  
  // page 88
  gdouble hm3[3*N+1][9*N+6];
  hm3[0][0] = 1.;
  for ( m = 1; m < 3*N+1; m++) {
    for ( n = 0; n <= 3*m; n++) {
      hm3[m][n] = 0.;
      for ( mu = 0; mu <= m; mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  hm3[m][n] += hm2[mu][nu]*hm1[m-mu][n-nu];
  	}
      }
    }
  }

  // page 92
  Vector w[2*N+1][2*N+1];
  for ( m = 0; m < 2*N+1; m++) {
    for ( n = 0; n <= m; n++) {
      w[m][n].x = w[m][n].y = w[m][n].z = 0.;
      
      for ( mu = max(0, m-N); mu <= min(m,N); mu++) {
  	for ( nu = max(0,n-m+mu); nu <= min(n,mu+1); nu++) {
  	  w[m][n] = vector_sum (w[m][n],vector_times_constant( vector_vector_product(&q(mu,nu),&q(m-mu,n-nu+1)),
  	  						       (mu+1)*(n-nu+1)));
  	}
      }

    }
  }

  // page 92
  gdouble w2[4*N+1][4*N+1];
  for ( m = 0; m < 4*N+1; m++) {
    for ( n = 0; n <= m; n++) {
      w2[m][n] = 0.;

      for ( mu = max(0,m-2*N); mu <= min(m,2*N); mu++) {
  	for ( nu = max(0, n-m+mu); nu <= min(n,mu); nu++) {
  	  w2[m][n] += vector_scalar_product (&w[mu][nu],&w[m-mu][n-nu]);
  	}
      }

    }
  }

  // page 92/92
  /* gdouble w1[8*N+1][8*N+3]; */
  gdouble w1[4*N+2][4*N+3];
  w1[0][0] = sqrt (w2[0][0]);
  for ( m = 1; m < 4*N+2; m++) {
    for ( n = 0; n <= m; n++) {

  	if ( m <= 4*N )
  	  w1[m][n] = w2[m][n];
  	else
  	  w1[m][n] = 0.;

  	for ( mu = 1; mu <= m-1; mu++) {
  	  for ( nu = max(0,n-m+mu); nu <= min(n,mu); nu++) {
  	    w1[m][n] -= w1[mu][nu]*w1[m-mu][n-nu];
  	  }
  	}
	
  	w1[m][n] /= (2.*w1[0][0]); // w1 here is instead of w2. I think this is right and the manuscript is wrong
    }
  }

  // page 93
  Vector cm[2*N][2*N+2];
  for ( m = 0; m < 2*N; m++) {
    for ( n = 0; n <= m+2; n++) {
      cm[m][n].x = cm[m][n].y = cm[m][n].z = 0.;
      
      for ( mu = max(0, m-N); mu <= min(m,N-1); mu++) {
  	for ( nu = max(0,n-m+mu); nu <= min(n,mu+2); nu++) {
  	  cm[m][n] = vector_sum (cm[m][n],vector_times_constant( vector_vector_product(&q(mu+1,nu),&q(m-mu,n-nu+1)),
  								 (mu+1)*(n-nu+1)));
	  
  	}
      }

    }
  }

  // page 94
  gdouble e[3*N+3][3*N+6];
  for ( m = 0; m < 3*N+3; m++) {
    for ( n = 0; n <= m+3; n++) {
      e[m][n] = 0.;

      for ( mu = max(0,m-N); mu <= min(m,2*N-1); mu++) {
  	for ( nu = max(0,n-m+mu-1); nu <= min(n,mu+2); nu++) {
  	  e[m][n] += vector_scalar_product (&cm[mu][nu], &q(m-mu,n-nu));
  	}
      }

    }
  }

  // page 95
  gdouble t[4*N+2][12*N+12];
  for ( m = 0; m < 4*N+2; m++) {
    for ( n = 0; n <= 3*m; n++) {
      t[m][n] = 0.;

      for ( nu = max(0,n-m); nu <= min(n,2*m); nu++) {
  	t[m][n] += s(m,nu)*w1[m][n-nu];
      }

    }
  }

  // page 96
  for ( m = 0; m < 4*N+2; m++) {
    for ( n = 0; n <= 3*m ; n++) {
      p(m,n) = 0.;
      for ( mu = 0; mu <= m; mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu); nu++) {
  	  p(m,n) += t[mu][nu]*hm1[m-mu][n-nu];
  	}
      }

    }
  }

  // page 96
  gdouble f[3*N+1][9*(N+1)+1];
  for ( m = 0; m < 3*N+1; m++) {
    for ( n = 0; n <= 9*(N+1)+1; n++) {
      f[m][n] = 0.;

      for ( nu = max(0,n-m-3); nu <= min(n,2*m); nu++) {
  	  f[m][n] += s(m,nu)*e[m][n-nu];
      }

    }
  }

  // page 96
  for ( m = 0; m < 3*N+1; m++) {
    for ( n = 0; n <= 3*m+3; n++) {
      g(m,n) = 0.;
      for ( mu = 0; mu <= min(m,3*N-1); mu++) {
  	for ( nu = max(0,n-3*m+3*mu); nu <= min(n,3*mu+3); nu++) {
  	  g(m,n) += f[mu][nu]*hm3[m-mu][n-nu];
  	}
      }
      
    }
  }

}

void calculate_triangle_self_influence_coeff (SIData * sid, gint k,
					      InfluenceCoeffs * ic, gdouble xie,
					      gdouble eta1, gdouble eta2, gboolean direct)
{
  gint M = 2*k-1;
  gint N = M-1;
  gint m, n;
  gdouble delta = 4*sid->a*sid->c-sid->b*sid->b;

  // page 99
  gdouble Seta1[8*N+3];
  gdouble Seta2[8*N+3];
  Seta1[0] = Seta2[0] = 1.;
  for ( m = 1; m < 8*N+3; m++) {
    Seta1[m] = Seta2[m] = 0.;
    for ( n = 0; n <= 2*m; n++) {
      Seta1[m] += s(m,n)*pow(eta1,n);
      Seta2[m] += s(m,n)*pow(eta2,n);
    }
  }

  // page 98
  gdouble YY[4*N+2][12*N+4+k];
  
  YY[0][0] = 1./sqrt(sid->c)*(asinh((2.*sid->c*eta2+sid->b)/sqrt(delta))
			    - asinh((2.*sid->c*eta1+sid->b)/sqrt(delta)));

  for ( m = 1; m < 4*N+2; m++) {
    YY[m][0] = 2./((2.*m-1.)*delta)*( (2.*sid->c*eta2 + sid->b)/sqrt(Seta2[2*m-1]) - (2.*sid->c*eta1 + sid->b)/sqrt(Seta1[2*m-1])
    				    + 4.*(m-1.)*sid->c*YY[m-1][0] );
  }

  for ( n = 1; n < 12*N+4+k; n++) {
    YY[0][n] = 1./(n*sid->c)*( pow(eta2,n-1)*sqrt(Seta2[1]) - pow(eta1,n-1)*sqrt(Seta1[1])
  			     - (2*n-1.)*sid->b/2.*YY[0][n-1] - (n-1)*sid->a*YY[0][n-2]);
  }

  for ( m = 1; m < 4*N+2; m++) {
    for ( n = 1; n <= 3*m+3+k; n++) {

      if ( n != 2*m) {
	YY[m][n] = 1./((n-2.*m)*sid->c)*( pow(eta2,n-1)/sqrt(Seta2[2*m-1]) - pow(eta1,n-1)/sqrt(Seta1[2*m-1])
					  - (2*n-2*m-1.)*sid->b/2.*YY[m][n-1] - (n-1)*sid->a*YY[m][n-2]);
      }
      else {
	YY[m][2*m] = -1./sid->c*( pow(eta2,2*m-1)/((2.*m-1.)*sqrt(Seta2[2*m-1]))
				  - pow(eta1,2*m-1)/((2.*m-1.)*sqrt(Seta1[2*m-1]))
				  + sid->b/2.*YY[m][2*m-1] - YY[m-1][2*(m-1)]);
      }

    }
  }


  // page 89
  gdouble XX[4*N+2*k];
  for ( m = 0; m < 4*N+2*k; m++)
    XX[m] = pow(xie, m+1)/(m+1);

  // page 90  OK only for triangles 1 - 2
  gint i, j;

  for ( i = 0; i < k; i++) {
    for ( j = 0; j < k; j++) {

      for ( m = 0; m < 4*N+1; m++) {
  	for ( n = 0; n <= 3*m; n++) {
	  gdouble old = gsl_matrix_get (ic->psi, j, i);
	 
	  if (direct == TRUE)
	    gsl_matrix_set (ic->psi, j, i,  old + XX[m+i+j]*p(m,n)*YY[m][n+j]);
	  else
	    gsl_matrix_set (ic->psi, j, i,  old + XX[m+i+j]*p(m,n)*YY[m][n+j]*pow(-1.,i));
  	}
      }

     

    }
  }

  // page 90
  for ( i = 0; i < k; i++) {
    for ( j = 0; j < k; j++) {

      for ( m = 0; m < 3*N; m++) {
  	for ( n = 0; n <= 3*m+3; n++) {
	  gdouble old = gsl_matrix_get (ic->phi, j, i);

	  if (direct == TRUE)
	    gsl_matrix_set (ic->phi, j, i, old + XX[m+i+j]*g(m,n)*YY[m+1][n+j]);
	  else
	    gsl_matrix_set (ic->phi, j, i, old + XX[m+i+j]*g(m,n)*YY[m+1][n+j]*pow(-1.,i));
  	}
      }
      
    }
  }
}

void calculate_triangle_self_influence_coeff2 (SIData * sid, gint k,
					       InfluenceCoeffs * ic, gdouble xie,
					       gdouble eta1, gdouble eta2, gboolean direct)
{
  gint M = 2*k-1;
  gint N = M-1;
  gint m, n;
  gdouble delta = 4*sid->a*sid->c-sid->b*sid->b;

  // page 99
  gdouble Seta1[8*N+2];
  gdouble Seta2[8*N+2];
  Seta1[0] = Seta2[0] = 1.;
  for ( m = 1; m < 8*N+2; m++) {
    Seta1[m] = Seta2[m] = 0.;
    for ( n = 0; n <= 2*m; n++) {
      Seta1[m] += s(m,n)*pow(eta1,n);
      Seta2[m] += s(m,n)*pow(eta2,n);
    }
  }

  // page 98
  gdouble YY[4*N+1][12*N+4+k];
  
  YY[0][0] = 1./sqrt(sid->c)*(asinh((2.*sid->c*eta2+sid->b)/sqrt(delta))
			    - asinh((2.*sid->c*eta1+sid->b)/sqrt(delta)));

  for ( m = 1; m < 4*N+1; m++) {
    YY[m][0] = 2./((2.*m-1.)*delta)*( (2.*sid->c*eta2 + sid->b)/sqrt(Seta2[2*m-1]) - (2.*sid->c*eta1 + sid->b)/sqrt(Seta1[2*m-1])
    				    + 4.*(m-1.)*sid->c*YY[m-1][0] );
  }

  for ( n = 1; n < 12*N+4+k; n++) {
    YY[0][n] = 1./(n*sid->c)*( pow(eta2,n-1)*sqrt(Seta2[1]) - pow(eta1,n-1)*sqrt(Seta1[1])
  			     - (2*n-1.)*sid->b/2.*YY[0][n-1] - (n-1)*sid->a*YY[0][n-2]);
  }

  for ( m = 1; m < 4*N+1; m++) {
    for ( n = 1; n <= 3*m+3+k; n++) {

      if ( n != 2*m) {
	YY[m][n] = 1./((n-2.*m)*sid->c)*( pow(eta2,n-1)/sqrt(Seta2[2*m-1]) - pow(eta1,n-1)/sqrt(Seta1[2*m-1])
					  - (2*n-2*m-1.)*sid->b/2.*YY[m][n-1] - (n-1)*sid->a*YY[m][n-2]);
      }
      else {
	YY[m][2*m] = -1./sid->c*( pow(eta2,2*m-1)/((2.*m-1.)*sqrt(Seta2[2*m-1]))
				  - pow(eta1,2*m-1)/((2.*m-1.)*sqrt(Seta1[2*m-1]))
				  + sid->b/2.*YY[m][2*m-1] - YY[m-1][2*(m-1)]);
      }

    }
  }


  // page 89
  gdouble XX[4*N-2+2*k];
  for ( m = 0; m < 4*N-2+2*k; m++)
    XX[m] = pow(xie, m+1)/(m+1);

  // page 90  OK only for triangles 1 - 2
  gint i, j;

  for ( i = 0; i < k; i++) {
    for ( j = 0; j < k; j++) {

      for ( m = 0; m < 4*N-1; m++) {
  	for ( n = 0; n <= 3*m; n++) {
	  gdouble old = gsl_matrix_get (ic->psi, j, i);

	  if (direct == TRUE)
	    gsl_matrix_set (ic->psi, j, i,  old + XX[m+i+j]*p(m,n)*YY[m][n+i]);
	  else
	    gsl_matrix_set (ic->psi, j, i,  old + XX[m+i+j]*p(m,n)*YY[m][n+i]*pow(-1.,i));
  	}
      }

    }
  }

  // page 90
  for ( i = 0; i < k; i++) {
    for ( j = 0; j < k; j++) {

      for ( m = 0; m < 3*N-1; m++) {
  	for ( n = 0; n <= 3*m+3; n++) {
	  gdouble old = gsl_matrix_get (ic->phi, j, i);

	  if (direct == TRUE)
	    gsl_matrix_set (ic->phi, j, i, old - XX[m+i+j]*g(m,n)*YY[m+1][n+i]);
	  else
	    gsl_matrix_set (ic->phi, j, i, old - XX[m+i+j]*g(m,n)*YY[m+1][n+i]*pow(-1.,i));
  	}
      }
      
    }
  }
}

static SIData * sidata_new (gint k)
{
  gint N = k-1;
  SIData * new = g_malloc (sizeof(SIData));

  new->q = new_2D_vector_array (N+1,N+2);
  new->s = new_2D_double_array (24*max(N,1), 48*max(N,1));
  new->p = new_2D_double_array (11*max(N,1), 33*max(N,1));
  new->g = new_2D_double_array (11*max(N,1), 33*max(N,1)+10);

  return new;
}

static void sidata_destroy (SIData * sid)
{
  g_assert (sid != NULL);
  g_array_free (sid->q, TRUE);
  g_array_free (sid->s, TRUE);
  g_array_free (sid->p, TRUE);
  g_array_free (sid->g, TRUE);
  g_free (sid);
}

typedef struct {
  S2P * s2p;
  SPPanel * spp;
  gint a, b;
  gdouble u, v;
  Point p;
  gboolean direct;
  gboolean psi;
} tmpStruct;

double f2 (double v, void * params) {
  tmpStruct * tmp = (tmpStruct *) params;
  S2P * s2p = tmp->s2p;
  SPPanel * spp = tmp->spp;
  Point pu = spline2d_eval_point (spp->sp, tmp->u, v);
  /* gdouble r = sqrt(pow(pu.x - tmp->p.x, 2.) + pow(pu.y - tmp->p.y, 2.) + pow(pu.z - tmp->p.z, 2.)); */

  Vector R;
  R.x = pu.x - tmp->p.x;
  R.y = pu.y - tmp->p.y;
  R.z = pu.z - tmp->p.z;
  gdouble r = vector_norm (R);

  Vector N = normal (s2p, spp->k, tmp->u,v);

  /* PSI */
  if (tmp->psi)
    return pow(v-s2p->ve, tmp->b)*jacobian(s2p, spp->k, tmp->u, v)/r;
  else
  /* PHI */
    return pow(v-s2p->ve, tmp->b)*jacobian(s2p, spp->k, tmp->u, v)*vector_scalar_product (&N, &R) / pow(r,3);
}

double f (double u, void * params) {
  tmpStruct * tmp = (tmpStruct *) params;
  S2P * s2p = tmp->s2p;
  SPPanel * spp = tmp->spp;

  tmp->u = u;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
  double result, error;
     
  gsl_function F;
  F.function = &f2;
  F.params = tmp;

  if (tmp->direct) {
    gdouble b1 = s2p->ve + (spp->v0 - s2p->ve)*(u-s2p->ue)/(spp->u1-s2p->ue);
    gdouble b2 = s2p->ve + (spp->v1 - s2p->ve)*(u-s2p->ue)/(spp->u1-s2p->ue);
    
    gsl_integration_qags (&F, b1, b2, 0, 1e-7, 1000,
			  w, &result, &error);
  }
  else {
    gdouble delta = (s2p->ue-u)/(s2p->ue-spp->u0);
    
    gdouble b1 = s2p->ve + (spp->v0 - s2p->ve)*delta;
    gdouble b2 = s2p->ve + (spp->v1 - s2p->ve)*delta;
    
    g_assert ( b1 <= spp->v1 && b1 >= spp->v0);
    g_assert ( b2 <= spp->v1 && b2 >= spp->v0);
    
    gsl_integration_qags (&F, b1, b2, 0, 1e-7, 1000,
			  w, &result, &error);
  }
     
  gsl_integration_workspace_free (w);

  return pow(u-s2p->ue, tmp->a)*result;
}

/* static gdouble direct_integration (SPPanel * spp, gint im, gint in, gint a, gint b, gboolean direct, gboolean psi) */
/* { */
/*   GaussPoints * gp = spp->outer; */
/*   S2P * s2p = g_ptr_array_index (spp->s2p, im + in*gp->ui->len); */

/*   tmpStruct * tmp = g_malloc (sizeof(tmpStruct)); */
/*   tmp->spp = spp; */
/*   tmp->s2p = s2p; */
/*   tmp->a = a; */
/*   tmp->b = b; */
/*   tmp->p = g_array_index (gp->Pi, Point, im + in*gp->ui->len); */
/*   tmp->direct = direct; */
/*   tmp->psi = psi; */

/*   gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000); */

/*   double result, error; */

/*   gsl_function F; */
/*   F.function = &f; */
/*   F.params = tmp; */

/*   if (direct) */
/*     gsl_integration_qags (&F, s2p->ue, spp->u1, 0, 1e-7, 1000, */
/* 			  w, &result, &error); */
/*   else */
/*     gsl_integration_qags (&F, spp->u0, s2p->ue, 0, 1e-7, 1000, */
/* 			  w, &result, &error); */

/*   gsl_integration_workspace_free (w); */

/*   g_free (tmp); */

/*   return result; */
/* } */

/**
 * Returns the self-influence coefficient for a panel spp
 * at a Gauss-Legendre point (um, vn).
 * These are computed following (Maniar,1995)'s method.
 **/
/* InfluenceCoeffs * sppanel_self_influence_coeff (SPPanel * spp, gint im, gint in) */
/* { */
/*   // See section 5.1 of (Maniar,1995) */
/*   gint M = 2*spp->k-1; */
/*   gint N = M-1; */
/*   gdouble eta1, eta2, xie; */
/*   InfluenceCoeffs * ic = influencecoeffs_new (spp->k); */
/*   SIData * sid = sidata_new (M); */
/*   gint m,n; */

/*   // Get the coefficient of the polynomial expansion of panel surface */
/*   // with the Gauss-Legendre point as expansion point */
/*   GaussPoints * gp = spp->outer; */
/*   S2P * s2p = g_ptr_array_index (spp->s2p, im + in*gp->ui->len);   */
  
/*   // For T2 Xm and Ymn only need to be recomputed with eta = -eta */
/*   // evaluating pmn, gmn is done twice, one for each pair (T1T2/T3T4) and Xm Ymn 4 times */
/*   // The analysis for the second pair of triangles can be directly obtained from the first */
/*   // pair of triangles if we simply "fill-in" the matrix qmn with xn,m-n+1 instead xm-n+1,n. */
/*   /\* fprintf(stderr, "B\n"); *\/ */
/*   // triangle 1 or 2 */
/*   // page 97 */

/*   // page 83 */
/*   Vector V0; */
/*   V0.x = V0.y = V0.z = 0.; */
/*   for ( n = 0; n < N+2; n++) */
/*     for ( m = 0; m < N+1; m++) { */
/*       if (m-n+1 < spp->k && n < spp->k) */
/*   	q(m,n) = g_array_index (s2p->x, Vector, (m-n+1) + n*spp->k); */
/*       else */
/*   	q(m,n) = V0; */
/*     } */

/*   // Middle page 83 */
/*   sid->a = vector_scalar_product(&g_array_index (s2p->x, Vector, 1),&g_array_index (s2p->x, Vector, 1)); */
/*   sid->b = 2.*vector_scalar_product(&g_array_index (s2p->x, Vector, 1),&g_array_index (s2p->x, Vector, spp->k)); */
/*   sid->c = vector_scalar_product(&g_array_index (s2p->x, Vector, spp->k),&g_array_index (s2p->x, Vector, spp->k)); */

/*   xie = eta1 = 0.; // Dummy only used for debugging purposes by self_influence_intermediate_calculations2 */
/*   self_influence_intermediate_calculations2 (sid, spp->k, xie, eta1); */

/*   // T1 */
/*   xie = spp->u1 - ui(im); */
/*   eta1 = (spp->v0 - vj(in))/xie; */
/*   eta2 = (spp->v1 - vj(in))/xie; */

/*   calculate_triangle_self_influence_coeff (sid, spp->k, ic, xie, eta1, eta2, TRUE); */

/*   // T2 */
/*   xie = ui(im) - spp->u0; */
/*   eta1 = (spp->v0 - vj(in))/xie; */
/*   eta2 = (spp->v1 - vj(in))/xie; */
  
/*   calculate_triangle_self_influence_coeff (sid, spp->k, ic, xie, eta1, eta2, FALSE); */



/*   // triangle 3 or 4 */
/*   // page 97 */
/*   for ( n = 0; n < N+2; n++) */
/*     for ( m = 0; m < N+1; m++) { */
/*       if (n < spp->k && m-n+1 < spp->k && m-n+1 >= 0) */
/*   	q(m,n) = g_array_index (s2p->x, Vector, n + (m-n+1)*spp->k); */
/*       else */
/*   	q(m,n) = V0; */
/*     } */

/*   // Middle page 83 */
/*   sid->a = vector_scalar_product(&g_array_index (s2p->x, Vector, spp->k),&g_array_index (s2p->x, Vector, spp->k)); */
/*   sid->b = 2.*vector_scalar_product(&g_array_index (s2p->x, Vector, 1),&g_array_index (s2p->x, Vector, spp->k)); */
/*   sid->c = vector_scalar_product(&g_array_index (s2p->x, Vector, 1),&g_array_index (s2p->x, Vector, 1)); */

/*   self_influence_intermediate_calculations (sid, spp->k, xie, eta1); */


/*   // T3 */
/*   xie = spp->v1 - vj(in); */
/*   eta1 = (spp->u0 - ui(im))/xie; */
/*   eta2 = (spp->u1 - ui(im))/xie; */

/*   calculate_triangle_self_influence_coeff2 (sid, spp->k, ic, xie, eta1, eta2, TRUE); */


/*   // T4 */
/*   xie = spp->v0 - vj(in); */
/*   eta2 = -(spp->u0 - ui(im))/xie; */
/*   eta1 = -(spp->u1 - ui(im))/xie; */

/*   calculate_triangle_self_influence_coeff2 (sid, spp->k, ic, xie, eta1, eta2, FALSE); */

/*   sidata_destroy (sid); */

/*   g_assert (ic); */
/*   return ic; */
/* } */

typedef struct
{
  gdouble s;
  gdouble d;
} Source;




gboolean need_refining (Spline2D * sp, Point p,
			 gdouble u0, gdouble u1,
			 gdouble v0, gdouble v1,
			 gint count)
{
  if (count > 5)
    return FALSE;

  Point pc = spline2d_eval_point (sp, (u0 + u1)/2., (v0 + v1)/2.);
  gdouble d = point_distance (pc, p);
  gdouble L = 0.;

  L = max (L, point_distance (p, spline2d_eval_point (sp, u0, v0)));
  L = max (L, point_distance (p, spline2d_eval_point (sp, u0, v1)));
  L = max (L, point_distance (p, spline2d_eval_point (sp, u1, v0)));
  L = max (L, point_distance (p, spline2d_eval_point (sp, u1, v1)));

  if (L/d > 0.66666)
    return FALSE;
  else 
    return TRUE;
}

/*** Structure for faster integration ***/
static void gpcell_fill (GPCell * gpc, Spline2D * sp)
{
  gint i, j, ng = gpc->ng;
  gint k = sp->k, NU = sp->NXU;
  gdouble di[ng][ng], d = 0., dtmp;

   /*  Init */
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui[ng], wi[ng], vj[ng], wj[ng];
  size_t ustart[ng], uend, vstart[ng], vend;
  size_t ustart_x[ng], uend_x;

  gpc->Ni = g_array_new (FALSE, FALSE, sizeof(Vector));
  gpc->Pi = g_array_new (FALSE, FALSE, sizeof(Point));
  gpc->wJij = g_array_new (FALSE, FALSE, sizeof(gdouble));
  
  gpc->Bu = g_ptr_array_new ();
  gpc->Bv = g_ptr_array_new ();
  gpc->Bux = g_ptr_array_new ();

  for ( i = 0; i < ng; i++) {
    gsl_integration_glfixed_point (gpc->u0, gpc->u1, i, &ui[i], &wi[i], itable);
    gsl_integration_glfixed_point (gpc->v0, gpc->v1, i, &vj[i], &wj[i], jtable);

    gsl_vector * Bu = gsl_vector_alloc (k);
    gsl_bspline_eval_nonzero (ui[i], Bu, &ustart[i], &uend, sp->w_u);
    g_ptr_array_add (gpc->Bu, Bu);

    gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
    gsl_bspline_deriv_eval_nonzero (vj[i], 1, Bv, &vstart[i], &vend, sp->w_v, sp->wd_v);
    g_ptr_array_add (gpc->Bv, Bv);

    gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
    gsl_bspline_deriv_eval_nonzero (ui[i], 1, Bux, &ustart_x[i], &uend_x, sp->wx_u, sp->wxd_u);
    g_ptr_array_add (gpc->Bux, Bux);
  }
  gsl_integration_glfixed_table_free (itable);
  gsl_integration_glfixed_table_free (jtable);

  // Stores the normal, the physical coordinates and the jacobian at the inner Gauss points
  for ( j = 0; j < ng; j++) {
    gsl_matrix * Bv = g_ptr_array_index (gpc->Bv, j);
    for ( i = 0; i < ng; i++) {
      gsl_matrix * Bux = g_ptr_array_index (gpc->Bux, i);

      Vector N;
      N.x = N.y = N.z = 0.;
      Point P;
      P.x = P.y = P.z = 0.;

      gint m, n;
      Vector xu, xv;
      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      
      for ( m = 0; m < k; m++) {
  	gdouble cu = gsl_matrix_get (Bux, m, 0);
  	gdouble cdu = gsl_matrix_get (Bux, m, 1);
	gint ii = (ustart_x[i]+m);
  	for ( n = 0; n < k; n++) {
  	  gdouble cv = gsl_matrix_get (Bv, n, 0);
  	  gdouble cudv = cu*gsl_matrix_get (Bv, n, 1);
  	  gdouble cvdu = cv*cdu;
  	  gdouble cuv = cu*cv;
	  gint jj = (vstart[j]+n);

	  SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, ii+ jj*NU);
	  gdouble v0 = sc->v[0];
	  gdouble v1 = sc->v[1];
	  gdouble v2 = sc->v[2];

  	  xu.x += v0*cvdu;
  	  xu.y += v1*cvdu;
  	  xu.z += v2*cvdu;
  	  xv.x += v0*cudv;
  	  xv.y += v1*cudv;
  	  xv.z += v2*cudv;
  	  P.x += v0*cuv;
  	  P.y += v1*cuv;
  	  P.z += v2*cuv;
  	}
      }

      // Calculate distance to each of the Gauss points
      // Compute the inner characteristic distance (minimum distance)
      dtmp = point_distance_squared (gpc->pe, P);
      if ( d < dtmp)
	d = dtmp;

      N = vector_vector_product (&xu, &xv);

      gdouble J = vector_norm (N);

      N.x /= J; N.y /= J; N.z /= J;
      
      g_array_append_val (gpc->Ni, N);
      g_array_append_val (gpc->Pi, P);

      J *= wi[i]*wj[j];
      g_array_append_val (gpc->wJij, J);
    }
  }

  gpc->d = d;
}

static GPCell * gpcell_new (Spline2D * sp, gdouble u0, gdouble u1, gdouble v0, gdouble v1, gint ng, gint level)
{
  GPCell * new = g_malloc (sizeof (GPCell));
  
  new->u0 = u0;
  new->u1 = u1;
  new->v0 = v0;
  new->v1 = v1;
  new->level = level;
  new->ue = (u0+u1)/2.;
  new->ve = (v0+v1)/2.;
  new->pe = spline2d_eval_point (sp, new->ue, new->ve);
  new->k = sp->k;
  new->ng = ng;

  new->Pi = NULL;
  new->Ni = NULL;
  new->wJij = NULL;
  new->Bu = NULL;
  new->Bv = NULL;
  new->Bux = NULL;

  new->children[0] = NULL;
  new->children[1] = NULL;
  new->children[2] = NULL;
  new->children[3] = NULL;

  gpcell_fill (new, sp);

  return new;
}

static void gpcell_destroy (GPCell * gpc)
{
  if (gpc->Ni != NULL) {
    g_array_free (gpc->Ni, TRUE);
    g_array_free (gpc->Pi, TRUE);
    g_array_free (gpc->wJij, TRUE);
    
    gint i;
    for ( i = 0; i < gpc->Bu->len; i++)
      gsl_vector_free (g_ptr_array_index (gpc->Bu, i));
    g_ptr_array_free (gpc->Bu, TRUE);
    for ( i = 0; i < gpc->Bv->len; i++)
      gsl_matrix_free (g_ptr_array_index (gpc->Bv, i));
    g_ptr_array_free (gpc->Bv, TRUE);
    for ( i = 0; i < gpc->Bux->len; i++)
      gsl_matrix_free (g_ptr_array_index (gpc->Bux, i));
    g_ptr_array_free (gpc->Bux, TRUE);
  }

  g_free (gpc);
}

GPCell * gpcell_tree_new (SPPanel * spp)
{
  GPCell * gpc = gpcell_new (spp->sp, spp->u0, spp->u1,
  			     spp->v0, spp->v1,
  			     spp->sp->ninner, 0);

  return gpc;
}

static void gpcell_destroy_recursive (GPCell * gpc)
{
  if (gpc->children[0] != NULL)
    gpcell_destroy_recursive ((GPCell *) gpc->children[0]);
  gpc->children[0] = NULL;

  if (gpc->children[1] != NULL)
    gpcell_destroy_recursive ((GPCell *) gpc->children[1]);
  gpc->children[1] = NULL;

  if (gpc->children[2] != NULL)
    gpcell_destroy_recursive ((GPCell *) gpc->children[2]);
  gpc->children[2] = NULL;

  if (gpc->children[3] != NULL)
    gpcell_destroy_recursive ((GPCell *) gpc->children[3]);
  gpc->children[3] = NULL;

  gpcell_destroy (gpc);
}

void gpcell_tree_destroy (GPCell * gpc)
{
  gpcell_destroy_recursive (gpc);

  /* if (gpc->children[0] != NULL) */
  /*   gpcell_destroy_recursive ((GPCell *) gpc->children[0]); */
  /* gpc->children[0] = NULL; */

  /* if (gpc->children[1] != NULL) */
  /*   gpcell_destroy_recursive ((GPCell *) gpc->children[1]); */
  /* gpc->children[1] = NULL; */

  /* if (gpc->children[2] != NULL) */
  /*   gpcell_destroy_recursive ((GPCell *) gpc->children[2]); */
  /* gpc->children[2] = NULL; */

  /* if (gpc->children[3] != NULL) */
  /*   gpcell_destroy_recursive ((GPCell *) gpc->children[3]); */
  /* gpc->children[3] = NULL; */

  /* g_free (gpc); */
}

/**
 * Recursive routine for the adaptive calculation of near-influence coeffcients
 **/
void spline_near_field_influence_coeff_recursive (SPPanel * spp,
						  GPCell * gpc,
						  Point p,
						  gdouble * psi,
						  gdouble * phi,
						  gboolean edge)
{
  // Here we use the adaptive subdivision method (see (Maniar,1995) section 5.2)
  gdouble Cr;
  gint ng = gpc->ng;
  gint i, j;
  
  // Calculate Cr for refinement criteria
  Cr = point_distance_squared (p, gpc->pe);

  if ( (gpc->d/Cr <= 0.3 || gpc->level == 6)&&(!edge || gpc->level > 0) ) { //Optimum according to (Maniar,1995) p101

    //if ( 3.*gpc->d < 2.*Cr || gpc->level == 6 ) { //Optimum according to (Maniar,1995) p101

    gint m, n;
    gint k = gpc->k;

    gdouble tmpBu[k][ng];
    gdouble tmpBv[k][ng];
    for ( i = 0; i < ng; i++ ) {
      gsl_vector * Bu = g_ptr_array_index (gpc->Bu, i);
      gsl_matrix * Bv = g_ptr_array_index (gpc->Bv, i);
      for ( j = 0; j < k; j++) {
	tmpBu[j][i] = gsl_vector_get (Bu, j);
	tmpBv[j][i] = gsl_matrix_get (Bv, j, 0);
      }
    }

    // Loop over Gauss-Legendre points
    gint ii = 0;
    for ( i = 0; i < ng; i++) {
      for ( j = 0; j < ng; j++) {
	ii = i+j*ng;

	// Physical coordinates of Gauss-Legendre Point
	Point pg = g_array_index (gpc->Pi, Point, ii);

	Vector R;
	R.x = pg.x - p.x;
	R.y = pg.y - p.y;
	R.z = pg.z - p.z;
	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z; 

	// Normal at Gauss-Legendre point
	Vector N = g_array_index (gpc->Ni, Vector, ii);
	  
	// Jacobian at Gauss-Legendre point
	gdouble c1 = g_array_index (gpc->wJij, gdouble, ii)/sqrt(r2);
	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

#if 1
	R.z = pg.z+p.z+2.*15.45;
	r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	gdouble c1_tmp = g_array_index (gpc->wJij, gdouble, ii)/sqrt(r2);
	c2 += c1_tmp*vector_scalar_product (&N, &R)/r2;
	c1 += c1_tmp;
#endif

	// Loop over the splines included in spp	
	gint mm = 0;
	for ( n = 0; n < k; n++) {
	  gdouble c3 = c1*tmpBv[n][j];
	  gdouble c4 = c2*tmpBv[n][j];
	  for ( m = 0; m < k; m++) {
	    psi[mm] += c3*tmpBu[m][i];
	    phi[mm] += c4*tmpBu[m][i];
	    mm++;
	  }
	}
      }
    }
    return;
  }

  // Subdivide the surface
  if ( gpc->children[0] == NULL ) {
    gpc->children[0] = gpcell_new (spp->sp, gpc->u0, gpc->ue, gpc->v0, gpc->ve, gpc->ng, gpc->level+1);
    gpc->children[1] = gpcell_new (spp->sp, gpc->u0, gpc->ue, gpc->ve, gpc->v1, gpc->ng, gpc->level+1);
    gpc->children[2] = gpcell_new (spp->sp, gpc->ue, gpc->u1, gpc->v0, gpc->ve, gpc->ng, gpc->level+1);
    gpc->children[3] = gpcell_new (spp->sp, gpc->ue, gpc->u1, gpc->ve, gpc->v1, gpc->ng, gpc->level+1);
  }

  spline_near_field_influence_coeff_recursive (spp, gpc->children[0], p, psi, phi, edge);
  spline_near_field_influence_coeff_recursive (spp, gpc->children[1], p, psi, phi, edge);
  spline_near_field_influence_coeff_recursive (spp, gpc->children[2], p, psi, phi, edge);
  spline_near_field_influence_coeff_recursive (spp, gpc->children[3], p, psi, phi, edge);
}

void spline_near_field_influence_coeff_recursive2 (SPPanel * spp,
						   GPCell * gpc,
						   GaussPoints * gp,
						   gdouble * psi,
						   gdouble * phi)
{
  // Here we use the adaptive subdivision method (see (Maniar,1995) section 5.2)
  gdouble Cr = 0./* G_MAXDOUBLE */;
  gint ng = gpc->ng;
  gint i, j;
  GArray * p = gp->Pi;
  
  // Calculate Cr for refinement criteria
  for ( i = 0; i < p->len; i++) {
    gdouble d = point_distance (g_array_index (p, Point, i), gpc->pe);
    if ( d > Cr )
      Cr = d;
  }

  if ( pow(gpc->d/Cr,2.) <= 0.4  || gpc->level == 3) { //Optimum according to (Maniar,1995) p101
    gint m, n, a;
    gint k = gpc->k;

    gdouble tmpBu[k][ng];
    gdouble tmpBv[k][ng];
    for ( i = 0; i < ng; i++ ) {
      gsl_matrix * Bu = g_ptr_array_index (gpc->Bu, i);
      gsl_matrix * Bv = g_ptr_array_index (gpc->Bv, i);
      for ( j = 0; j < k; j++) {
	tmpBu[j][i] = gsl_matrix_get (Bu, j, 0);
	tmpBv[j][i] = gsl_matrix_get (Bv, j, 0);
      }
    }

    // Loop over Gauss-Legendre points
    for ( i = 0; i < ng; i++) {
 
      for ( j = 0; j < ng; j++) {

	// Physical coordinates of Gauss-Legendre Point
	Point pg = g_array_index (gpc->Pi, Point, i+j*ng);

	// Normal at Gauss-Legendre point
	Vector N = g_array_index (gpc->Ni, Vector, i+j*ng);

	// Weight at Gauss-Legendre point
	gdouble wJij = g_array_index (gpc->wJij, gdouble, i+j*ng);

	Vector R;
	gdouble r2;
	for ( a = 0; a < p->len; a++) {
	  Point pp = g_array_index (p, Point, a);
	  R.x = pg.x - pp.x;
	  R.y = pg.y - pp.y;
	  R.z = pg.z - pp.z;
	  r2 = R.x*R.x + R.y*R.y + R.z*R.z;

	  gdouble c1 = wJij/sqrt(r2);
	  gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

	  // Loop over the splines included in spp
	  for ( m = 0; m < k; m++) {
	    gdouble c3 = c1*tmpBu[m][i];
	    gdouble c4 = c2*tmpBu[m][i];
	    for ( n = 0; n < k; n++) {
	      psi[m + n*k + a*k*k] += c3*tmpBv[n][j];
	      phi[m + n*k + a*k*k] += c4*tmpBv[n][j];
	    }
	  }
	}
      }
    }
    return;
  }

  // Subdivide the surface
  if ( gpc->children[0] == NULL ) {
    gpc->children[0] = gpcell_new (spp->sp, gpc->u0, gpc->ue, gpc->v0, gpc->ve, gpc->ng, gpc->level+1);
    gpc->children[1] = gpcell_new (spp->sp, gpc->u0, gpc->ue, gpc->ve, gpc->v1, gpc->ng, gpc->level+1);
    gpc->children[2] = gpcell_new (spp->sp, gpc->ue, gpc->u1, gpc->v0, gpc->ve, gpc->ng, gpc->level+1);
    gpc->children[3] = gpcell_new (spp->sp, gpc->ue, gpc->u1, gpc->ve, gpc->v1, gpc->ng, gpc->level+1);
  }

  spline_near_field_influence_coeff_recursive2 (spp, gpc->children[0], gp, psi, phi);
  spline_near_field_influence_coeff_recursive2 (spp, gpc->children[1], gp, psi, phi);
  spline_near_field_influence_coeff_recursive2 (spp, gpc->children[2], gp, psi, phi);
  spline_near_field_influence_coeff_recursive2 (spp, gpc->children[3], gp, psi, phi);
}

/* static InfluenceCoeffs *  spline_near_field_adaptive_subdivision (SPPanel * spp, */
/* 								  gdouble u1, gdouble u2, */
/* 								  gdouble v1, gdouble v2, */
/* 								  Point p) */
/* { */
/*   Spline2D * sp = spp->sp; */
/*   gdouble d, Cr; */
/*   gdouble alpha = 0.4; //Optimum according to (Maniar,1995) p101 */
/*   gint ng = sp->ninner; // Optimum number of Gauss points according to (Maniar,1995) p101 */
/*   Point pg[ng][ng]; */
/*   gdouble wi[ng][ng]; */
/*   gdouble ui[ng]; */
/*   gdouble vi[ng]; */
/*   gdouble di[ng][ng]; */
/*   gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng); */
/*   gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng); */
/*   InfluenceCoeffs * ic = influencecoeffs_new (spp->k); */
/*   gint i, j, m, n; */
/*   gdouble psi[spp->k][spp->k]; */
/*   gdouble phi[spp->k][spp->k]; */

/*   for ( m = 0; m < spp->k; m++) { */
/*     for ( n = 0; n < spp->k; n++) { */
/*       psi[m][n] = 0.; */
/*       phi[m][n] = 0.; */
/*     } */
/*   } */

/*   // Centroid of the surface or sub-surface */
/*   gdouble um = (u1+u2)/2.; */
/*   gdouble vm = (v1+v2)/2.; */
/*   Point pm = spline2d_eval_point (spp->sp, um, vm); */

/*   // Calculate d, Cr */
/*   d = 0.; */
/*   Cr = point_distance (p, pm); */
  
/*   // Fill in the Gauss points */
/*   gdouble wui, wvi; */
/*   for ( i = 0; i < ng; i++) { */
/*     gsl_integration_glfixed_point (u1, u2, i, &ui[i], &wui, itable); */
/*     for ( j = 0; j < ng; j++) { */
/*       gsl_integration_glfixed_point (v1, v2, j, &vi[j], &wvi, jtable); */
/*       pg[i][j] = spline2d_eval_point (spp->sp, ui[i], vi[j]); */
/*       di[i][j] = point_distance (pm, pg[i][j]); */
/*       wi[i][j] = wui*wvi; */
/*     } */
/*   } */

/*   // Compute the inner characteristic distance */
/*   for ( i = 0; i < ng; i ++) { */
/*     for ( j = 0; j < ng; j ++) { */
/*       if ( d < di[i][j]) */
/* 	d = di[i][j]; */
/*     } */
/*   } */

/*   if ( pow(d/Cr,2.) <= alpha ) { */
/*     gsl_matrix * Bu = gsl_matrix_alloc (spp->k, 2); */
/*     gsl_matrix * Bv = gsl_matrix_alloc (spp->k, 2); */
/*     size_t istart, iend, jstart, jend; */
/*     // Loop over Gauss-Legendre points */
/*     for ( i = 0; i < ng; i++) { */
/*       gsl_bspline_deriv_eval_nonzero (ui[i], 1, Bu, &istart, &iend, sp->w_u, sp->wd_u); */
/*       if (sp->periodic) */
/* 	istart -= (sp->k-1); */
/*       for ( j = 0; j < ng; j++) { */
/* 	gsl_bspline_deriv_eval_nonzero (vi[j], 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v); */
/* 	Vector R; */
/* 	R.x = pg[i][j].x - p.x; */
/* 	R.y = pg[i][j].y - p.y; */
/* 	R.z = pg[i][j].z - p.z; */

/* 	//gdouble r = vector_norm (R); */
/* 	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z; */

/* 	Vector xu, xv; */
/* 	gint  k, l; */
/* 	xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.; */

/* 	for ( k = 0; k < sp->k; k++) { */
/* 	  gdouble cu = gsl_matrix_get (Bu, k, 0); */
/* 	  gdouble cdu = gsl_matrix_get (Bu, k, 1); */
/* 	  gint ii = (istart+k); */
/* 	  for ( l = 0; l < sp->k; l++) { */
/* 	    gdouble cudv = cu*gsl_matrix_get (Bv, l, 1); */
/* 	    gdouble cvdu = gsl_matrix_get (Bv, l, 0)*cdu; */

/* 	    gdouble v0 = coeff (sp, ii, (jstart+l), 0); */
/* 	    gdouble v1 = coeff (sp, ii, (jstart+l), 1); */
/* 	    gdouble v2 = coeff (sp, ii, (jstart+l), 2); */

/* 	    xu.x += v0*cvdu; */
/* 	    xu.y += v1*cvdu; */
/* 	    xu.z += v2*cvdu; */
/* 	    xv.x += v0*cudv; */
/* 	    xv.y += v1*cudv; */
/* 	    xv.z += v2*cudv; */
/* 	  } */
/* 	} */

/* 	Vector N = vector_vector_product (&xu, &xv); */

/* 	// Jacobian */
/* 	gdouble J = vector_norm (N); */

/* 	gdouble c1 = wi[i][j]/sqrt(r2)*J; */
	
/* 	// Normal normalise */
/* 	N.x /= J; N.y /= J; N.z /= J; */
	
/* 	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2; */

/* 	// Loop over the splines included in spp */
/* 	for ( m = 0; m < spp->k; m++) { */
/* 	  gdouble c3 = c1*gsl_matrix_get(Bu, m, 0); */
/* 	  gdouble c4 = c2*gsl_matrix_get(Bu, m, 0); */
/* 	  for ( n = 0; n < spp->k; n++) { */
/* 	    psi[m][n] += c3*gsl_matrix_get(Bv, n, 0); */
/* 	    phi[m][n] += c4*gsl_matrix_get(Bv, n, 0); */
/* 	  } */
/* 	} */
/*       } */
/*     } */

/*     for ( m = 0; m < spp->k; m++) { */
/*       for ( n = 0; n < spp->k; n++) { */
/* 	gsl_matrix_set (ic->psi, n, m, psi[m][n]); */
/* 	gsl_matrix_set (ic->phi, n, m, phi[m][n]); */
/*       } */
/*     } */
    
/*     gsl_matrix_free (Bu); */
/*     gsl_matrix_free (Bv); */

/*     gsl_integration_glfixed_table_free (itable); */
/*     gsl_integration_glfixed_table_free (jtable); */

/*     return ic; */
/*   } */

/*   // Subdivide the surface */
/*   ic = influencecoeffs_add_and_destroy (ic, */
/* 					spline_near_field_adaptive_subdivision (spp, u1, um, v1, vm, p)); */
/*   ic = influencecoeffs_add_and_destroy (ic, */
/* 					spline_near_field_adaptive_subdivision (spp, u1, um, vm, v2, p)); */
/*   ic = influencecoeffs_add_and_destroy (ic, */
/* 					spline_near_field_adaptive_subdivision (spp, um, u2, v1, vm, p)); */
/*   ic = influencecoeffs_add_and_destroy (ic, */
/* 					spline_near_field_adaptive_subdivision (spp, um, u2, vm, v2, p)); */

/*   return ic; */
/* } */

gdouble spp_characteristic_length (SPPanel * spp, Spline2D * sp)
{
  gdouble L = 0.;
  L = max (L, point_distance (spp->pe, spline2d_eval_point (sp, spp->u0, spp->v0)));
  L = max (L, point_distance (spp->pe, spline2d_eval_point (sp, spp->u0, spp->v1)));
  L = max (L, point_distance (spp->pe, spline2d_eval_point (sp, spp->u1, spp->v0)));
  L = max (L, point_distance (spp->pe, spline2d_eval_point (sp, spp->u1, spp->v1)));

  return L;
}

/******* LACHAT-WATSON Method for self-influence coefficients ********/

static gdouble xi (gdouble xi, gdouble eta, gdouble xip, gint se)
{
  switch (se) {
  case 1: return 0.5*((1.+xi)*eta + (1.-xi)*xip);
  case 2: return 0.5*((1.+xi) + (1-xi)*xip);
  case 3: return 0.5*(-(1.+xi)*eta + (1.-xi)*xip);
  case 4: return 0.5*(-(1.+xi) + (1.-xi)*xip);
  default: g_assert_not_reached ();
  }
}

static gdouble eta (gdouble xi, gdouble eta, gdouble etap, gint se)
{
  switch (se) {
  case 1: return 0.5*(-(1.+xi) + (1-xi)*etap);
  case 2: return 0.5*((1.+xi)*eta + (1-xi)*etap);
  case 3: return 0.5*((1.+xi) + (1-xi)*etap);
  case 4: return 0.5*(-(1.+xi)*eta + (1-xi)*etap);
  default: g_assert_not_reached ();
  }
}

static gdouble J_se (gdouble xi, gdouble eta, gdouble xip, gdouble etap, gint se)
{
  switch (se) {
  case 1: return 0.25*(1.+xi)*(1.+etap);
  case 2: return 0.25*(1.+xi)*(1.-xip);
  case 3: return 0.25*(1.+xi)*(1.-etap);
  case 4: return 0.25*(1.+xi)*(1.+xip);
  default: g_assert_not_reached ();
  }
}

static void LW_coeffs (gdouble xi, gdouble eta, gdouble xip, gdouble etap, gdouble * xi_se, gdouble * eta_se, gdouble * Jse)
{
  gdouble xip1 = 0.5*(xi+1.); // is not case dependent
  gdouble xim1 = 0.5*(1.-xi); // is not case dependent
  gdouble half_xip1 = 0.5*xip1; // is not case dependent
  gdouble xip1_eta = xip1*eta; // is not case dependent

  gdouble xim1_xip = 1.+xim1*xip;
  gdouble xim1_etap = 1.+xim1*etap;


  *(xi_se) = xip1_eta + xim1_xip;
  *(xi_se+1) = xip1 + xim1_xip;
  *(xi_se+2) = -xip1_eta + xim1_xip;
  *(xi_se+3) = -xip1 + xim1_xip;

  *(eta_se) = -xip1 + xim1_etap;
  *(eta_se+1) = xip1_eta + xim1_etap;
  *(eta_se+2) = xip1 + xim1_etap;
  *(eta_se+3) = -xip1_eta + xim1_etap;

  *(Jse) = half_xip1*(1.+etap);
  *(Jse+1) = half_xip1*(1.-xip);
  *(Jse+2) = half_xip1*(1.-etap);
  *(Jse+3) = half_xip1*(1.+xip);
}

/**
 * Calculate the spline self-influence coefficients usind of method of
 * (Lachat and Watson, 1976), Effective numerical treatment of boundary
 * integral equations: A formulation for three dimensional elastostatics.
 **/
static const double ui_35[35] =   /* abscissea for 35 points Gauss-Lengendre rule */
{
  -0.9977065690996,  
  -0.9879357644439,  
  -0.9704376160392,  
  -0.9453451482078,  
  -0.9128542613593,  
  -0.8732191250252,  
  -0.8267498990922,  
  -0.7738102522869,  
  -0.7148145015566,  
  -0.6502243646659,  
  -0.5805453447498,  
  -0.5063227732415,  
  -0.4281375415178,  
  -0.3466015544308,  
  -0.2623529412093,  
  -0.1760510611660,  
  -0.0883713432757,  
  -0.0000000000000,  
  0.0883713432757,  
  0.1760510611660,  
  0.2623529412093,  
  0.3466015544308,  
  0.4281375415178,  
  0.5063227732415,  
  0.5805453447498,  
  0.6502243646659,  
  0.7148145015566,  
  0.7738102522869,  
  0.8267498990922,  
  0.8732191250252,  
  0.9128542613593,  
  0.9453451482078,  
  0.9704376160392,  
  0.9879357644439,  
  0.9977065690996
};

static const double wi_35[35] =   /* weights for 35 points Gauss-Lengendre rule */
{
  0.0058834334204,  
  0.0136508283484,  
  0.0213229799115,  
  0.0288292601089,  
  0.0361101158635,  
  0.0431084223262,  
  0.0497693704014,  
  0.0560408162124,  
  0.0618736719661,  
  0.0672222852691,  
  0.0720447947726,  
  0.0763034571554,  
  0.0799649422423,  
  0.0830005937289,  
  0.0853866533921,  
  0.0871044469972,  
  0.0881405304303,  
  0.0884867949071,  
  0.0881405304303,  
  0.0871044469972,  
  0.0853866533921,  
  0.0830005937289,  
  0.0799649422423,  
  0.0763034571554,  
  0.0720447947726,  
  0.0672222852691,  
  0.0618736719661,  
  0.0560408162124,  
  0.0497693704014,  
  0.0431084223262,  
  0.0361101158635,  
  0.0288292601089,  
  0.0213229799115,  
  0.0136508283484,  
  0.0058834334204 
};

static const double ui_10[10] =   /* abscissea for 10 points Gauss-Lengendre rule */
{
  -0.9739065285172,
  -0.8650633666890,
  -0.6794095682990,
  -0.4333953941292,
  -0.1488743389816,
  0.1488743389816,
  0.4333953941292,
  0.6794095682990,
  0.8650633666890,
  0.9739065285172
};

static const double wi_10[10] =   /* weights for 10 points Gauss-Lengendre rule */
{
  0.0666713443087, 
  0.1494513491506, 
  0.2190863625160, 
  0.2692667193100, 
  0.2955242247148, 
  0.2955242247148, 
  0.2692667193100, 
  0.2190863625160, 
  0.1494513491506, 
  0.0666713443087
};

static const double ui_16[16] =   /* abscissea for 10 points Gauss-Lengendre rule */
{
  -0.989400934991,
  -0.944575023073,
  -0.865631202387,
  -0.755404408355,
  -0.617876244402,
  -0.458016777657,
  -0.281603550779,
  -0.095012509837,
  0.0950125098376,
  0.281603550779,
  0.458016777657,
  0.617876244402,
  0.75540440835,
  0.86563120238,
  0.94457502307,
  0.98940093499
};

static const double wi_16[16] =   /* weights for 10 points Gauss-Lengendre rule */
{
  0.0271524594117,
  0.0622535239386,
  0.0951585116824,
  0.1246289712555,
  0.1495959888165,
  0.1691565193950,
  0.1826034150449,
  0.1894506104550,
  0.1894506104550,
  0.1826034150449,
  0.1691565193950,
  0.1495959888165,
  0.1246289712555,
  0.0951585116824,
  0.0622535239386,
  0.0271524594117
};

void lachat_watson_self_influence_coefficients (SPPanel * spp,
						gdouble up, gdouble vp,
						Point p,
						gdouble * psi,
						gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 10, i, j, m, n, se;
  gdouble du = 0.5*(spp->u1-spp->u0);
  gdouble dv = 0.5*(spp->v1-spp->v0);
  gdouble xip = -1. + (up-spp->u0)/du;
  gdouble etap = -1. + (vp-spp->v0)/dv;
  gint k = spp->k, NU = sp->NXU;
  
  /* g_assert_not_reached (); */

  gsl_vector * Bu = gsl_vector_alloc (k);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  
  size_t istart, iend, jstart, jend;
  size_t istart_x, iend_x;

  // Loop over the Gauss-Point
  for ( i = 0; i < ng; i++) {
    for ( j = 0; j < ng; j++) {
      gdouble Jse[4], xi_se[4], eta_se[4];
      LW_coeffs (ui_10[i], ui_10[j], xip, etap, &xi_se[0], &eta_se[0], &Jse[0]);

      // Loop over the triangles
      for ( se = 0; se < 4; se ++) {

	gdouble u = spp->u0 + xi_se[se]*du;
	gdouble v = spp->v0 + eta_se[se]*dv;
	gsl_bspline_eval_nonzero (u, Bu, &istart, &iend, sp->w_u);
	gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
	gsl_bspline_deriv_eval_nonzero (u, 1, Bux, &istart_x, &iend_x, sp->wx_u, sp->wxd_u);

	Point pg;
	Vector xu, xv;
	gint  a, b;

	xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
	pg.x = pg.y = pg.z = 0.;

	gint ii = istart_x;
	for ( a = 0; a < k; a++) {
	  gdouble cu = gsl_matrix_get (Bux, a, 0);
	  gdouble cdu = gsl_matrix_get (Bux, a, 1);
	  gint jj = jstart;
	  for ( b = 0; b < k; b++) {
	    gdouble cv = gsl_matrix_get (Bv, b, 0);
	    gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	    gdouble cvdu = cv*cdu;
	    gdouble cuv = cu*cv;

	    SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, ii+ jj*NU);
	    gdouble v0 = sc->v[0];
	    gdouble v1 = sc->v[1];
	    gdouble v2 = sc->v[2];

	    xu.x += v0*cvdu;
	    xu.y += v1*cvdu;
	    xu.z += v2*cvdu;
	    xv.x += v0*cudv;
	    xv.y += v1*cudv;
	    xv.z += v2*cudv;
	    pg.x += v0*cuv;
	    pg.y += v1*cuv;
	    pg.z += v2*cuv;
	    jj++;
	  }
	  ii++;
	}
	
	// Normal at Gauss-Legendre point
	Vector N = vector_vector_product (&xu, &xv);

	// Physical distance to Gauss-Legendre Point
	Vector R;
	R.x = pg.x - p.x;
	R.y = pg.y - p.y;
	R.z = pg.z - p.z;
	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
	// Jacobian at Gauss-Legendre point = norm of N
	gdouble J = vector_norm (N);
	gdouble c1 = wi_10[i]*wi_10[j]*Jse[se]/sqrt(r2);
	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

#if 1
	R.z = pg.z+p.z+2.*(15.45);
	r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	gdouble c1_tmp = wi_10[i]*wi_10[j]*Jse[se]/sqrt(r2);
	c2 += c1_tmp*vector_scalar_product (&N, &R)/r2;
	c1 += c1_tmp;
#endif

	c1 *= J;

	// Loop over the splines included in spp
	gint mm = 0;
	for ( n = 0; n < k; n++) {
	  gdouble c3 = c1*gsl_matrix_get (Bv, n, 0);
	  gdouble c4 = c2*gsl_matrix_get (Bv, n, 0);
	  for ( m = 0; m < k; m++) {
	    psi[mm] += c3*gsl_vector_get (Bu, m);
	    phi[mm] += c4*gsl_vector_get (Bu, m);
	    mm++;
	  }
	}

      }
    }
  }

  gdouble J = du*dv;
  for ( m = 0; m < k*k; m++) {
    psi[m] *= J;
    phi[m] *= J;
  }

  gsl_vector_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);
}

void lachat_watson_ye_self_influence_coefficients (SPPanel * spp,
						   gdouble up, gdouble vp,
						   Point p,
						   gdouble * psi,
						   gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 16, i, j, m, n, se;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui[ng], wi[ng];
  gdouble vj[ng], wj[ng];
  gdouble du = (spp->u1-spp->u0)/2.;
  gdouble dv = (spp->v1-spp->v0)/2.;
  gdouble xip = -1 + (up-spp->u0)/du;
  gdouble etap = -1 + (vp-spp->v0)/dv;
  gint k = spp->k;

  for ( i = 0; i < ng; i++) {
    gsl_integration_glfixed_point (/* -1. */0., 1., i, &ui[i], &wi[i], itable);
    gsl_integration_glfixed_point (-1., 1., i, &vj[i], &wj[i], jtable);
  } // Could be precalculated
  
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;

  // Loop over the triangles
  for ( se = 1; se <= 1/* 4 */; se ++) {

    // Loop over the Gauss-Point
    for ( i = 0; i < ng; i++) {
      /* gdouble uuuu = ui[i]; */
      
      /* gdouble uuu = 2*uuuu-1.; */

      /* gdouble tmp = (ui[i]+1.)/2.; */
      /* gdouble uuu = 1./(5.-1.)*((1.-pow(tmp,5))/(1.-tmp)-1.); */

      /* gdouble J_u = 1./(5.-1.)*0.5*( (1.-pow(tmp,5))/(1.-pow(1-tmp,2.)) - 5.*pow(tmp,5.-1.)/(1.-tmp)); */

      gdouble uuu = 1./(5.-1.)*((1.-pow(ui[i],5))/(1.-ui[i])-1.);

      gdouble J_u = 1./(5.-1.)*( (1.-pow(ui[i],5))/(pow(1-ui[i],2.)) - 5.*pow(ui[i],5.-1.)/(1.-ui[i]));

      /* gdouble uu = uuu*2.-1.; */

      /* uuu = ui[i]*2-1.; */
      /* J_u = 1.; */

      for ( j = 0; j < ng; j++) {
	/* gdouble uuu = (xi (ui[i], vj[j], xip, se)+1.)/2.; */
	/* gdouble uu = 1./(5.-1.)*((1.-pow(uuu,5))/(1.-uuu)-1.); */
	/* uu = uu*2.-1.; */
	/* 	J_u = 1./(5.-1.)*( (1.-pow(uuu,5))/(1.-pow(1-uuu,2.)) - 5.*pow(uuu,5.-1.)/(1.-uuu)); */


	gdouble uu = xi (/* ui[i] */uuu*2.-1., vj[j], xip, se);
	//	fprintf (stdout, "%f \n",uu);
	gdouble vv = eta (/* ui[i] */uuu*2.-1., vj[j], etap, se);
	gdouble u = spp->u0 + (uu+1.)*du;
	gdouble v = spp->v0 + (vv+1.)*dv;
	
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
		
	Point pg;
	Vector xu, xv;
	gint  a, b;

	xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
	pg.x = pg.y = pg.z = 0.;

	for ( a = 0; a < k; a++) {
	  gdouble cu = gsl_matrix_get (Bux, a, 0);
	  gdouble cdu = gsl_matrix_get (Bux, a, 1);
	  gint ii = (istart+a);
	  for ( b = 0; b < k; b++) {
	    gdouble cv = gsl_matrix_get (Bv, b, 0);
	    gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	    gdouble cvdu = cv*cdu;
	    gdouble cuv = cu*cv;
	    gint jj = (jstart+b);
	    gdouble v0 = coeff (sp, ii, jj, 0);
	    gdouble v1 = coeff (sp, ii, jj, 1);
	    gdouble v2 = coeff (sp, ii, jj, 2);

	    xu.x += v0*cvdu;
	    xu.y += v1*cvdu;
	    xu.z += v2*cvdu;
	    xv.x += v0*cudv;
	    xv.y += v1*cudv;
	    xv.z += v2*cudv;
	    pg.x += v0*cuv;
	    pg.y += v1*cuv;
	    pg.z += v2*cuv;
	  }
	}
	
	// Normal at Gauss-Legendre point
	Vector N = vector_vector_product (&xu, &xv);

	// Physical distance to Gauss-Legendre Point
	Vector R;
	R.x = pg.x - p.x;
	R.y = pg.y - p.y;
	R.z = pg.z - p.z;
	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
	// Jacobian at Gauss-Legendre point = norm of N
	gdouble J = vector_norm (N);
	gdouble c1 = wi[i]*wj[j]*J_se (/* ui[i] */uuu*2.-1., vj[j], xip, etap, se)*J/sqrt(r2)*J_u;
	N.x /= J; N.y /= J; N.z /= J;
	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

	// Loop over the splines included in spp
	for ( m = 0; m < k; m++) {
	  gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	  gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	  for ( n = 0; n < k; n++) {
	    psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	    phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	  }
	}
      }
    }
  }

  gdouble J = du*dv/* /2. */;
  for ( m = 0; m < k*k; m++) {
    psi[m] *= J;
    phi[m] *= J;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

  gsl_integration_glfixed_table_free (itable);
  gsl_integration_glfixed_table_free (jtable);
}

void lachat_watson_self_influence_coefficients_adaptive (SPPanel * spp,
							 gdouble up, gdouble vp,
							 Point p,
							 gdouble * psi,
							 gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 3, i, j, m, n, se;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui, wi, vj, wj;
  gdouble du = (spp->u1-spp->u0)/2.;
  gdouble dv = (spp->v1-spp->v0)/2.;
  gdouble xip = -1 + (up-spp->u0)/du;
  gdouble etap = -1 + (vp-spp->v0)/dv;
  gint k = spp->k;
  
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;

  gint nsub = 20;
  gint mu, nu;
  gdouble u0 = 0., v0 = 0., u1 = -1., v1 = -1.;

  // Loop over the triangles
  for ( se = 1; se <= 4; se ++) {
    u1 = -1.;
    for ( mu = 0; mu < nsub; mu++) {
      /* u0 = -1. + (2.*mu)/nsub; */
      /* u1 = -1. + 2.*(mu+1.)/nsub; */

      u0 = u1;
      u1 = -1. + 2./pow(2., nsub-mu-1);
      

      // Loop over the Gauss-Point
      for ( i = 0; i < ng; i++) {
	gsl_integration_glfixed_point (u0, u1, i, &ui, &wi, itable);

	for ( nu = 0; nu < nsub; nu++) {
	  v0 = -1. + (2.*nu)/nsub;
	  v1 = -1. + 2.*(nu+1.)/nsub;
	
	  for ( j = 0; j < ng; j++) {
	    gsl_integration_glfixed_point (v0, v1, j, &vj, &wj, jtable);
	    gdouble uu = xi (ui, vj, xip, se);
	    gdouble vv = eta (ui, vj, etap, se);
	    gdouble u = spp->u0 + (uu+1.)*du;
	    gdouble v = spp->v0 + (vv+1.)*dv;

	    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
	    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
	     gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);
	
	    if (sp->periodic)
	      istart -= (k-1);
	
	    Point pg;
	    Vector xu, xv;
	    gint  a, b;

	    xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
	    pg.x = pg.y = pg.z = 0.;

	    for ( a = 0; a < k; a++) {
	      gdouble cu = gsl_matrix_get (Bux, a, 0);
	      gdouble cdu = gsl_matrix_get (Bux, a, 1);
	      gint ii = (istart+a);
	      for ( b = 0; b < k; b++) {
		gdouble cv = gsl_matrix_get (Bv, b, 0);
		gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
		gdouble cvdu = cv*cdu;
		gdouble cuv = cu*cv;
		gint jj = (jstart+b);
		gdouble v0 = coeff (sp, ii, jj, 0);
		gdouble v1 = coeff (sp, ii, jj, 1);
		gdouble v2 = coeff (sp, ii, jj, 2);

		xu.x += v0*cvdu;
		xu.y += v1*cvdu;
		xu.z += v2*cvdu;
		xv.x += v0*cudv;
		xv.y += v1*cudv;
		xv.z += v2*cudv;
		pg.x += v0*cuv;
		pg.y += v1*cuv;
		pg.z += v2*cuv;
	      }
	    }
	
	    // Normal at Gauss-Legendre point
	    Vector N = vector_vector_product (&xu, &xv);

	    // Physical distance to Gauss-Legendre Point
	    Vector R;
	    R.x = pg.x - p.x;
	    R.y = pg.y - p.y;
	    R.z = pg.z - p.z;
	    gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
	    // Jacobian at Gauss-Legendre point = norm of N
	    gdouble J = vector_norm (N);
	    gdouble c1 = wi*wj*J_se (ui, vj, xip, etap, se)*vector_norm (N)/sqrt(r2);
	    N.x /= J; N.y /= J; N.z /= J;
	    gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

	    // Loop over the splines included in spp
	    for ( m = 0; m < k; m++) {
	      gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	      gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	      for ( n = 0; n < k; n++) {
		psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
		phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	      }
	    }
	  }
	}

      }
    }

  }

  gdouble J = du*dv/* /(nsub*nsub) */;
  for ( m = 0; m < k*k; m++) {
    psi[m] *= J;
    phi[m] *= J;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

  gsl_integration_glfixed_table_free (itable);
  gsl_integration_glfixed_table_free (jtable);
}

/**
 * Quadrature-based self-influence coefficient integration method adapted from
 * function intkt22 kindly provided by Junjie Rong from the College of
 * Astronautics, Northwestern Polytechnical University, Xi'an, China.
 *
 * The details of the method can be found in the case of triangular elements
 * in the following paper:
 * Junjie Rong, Lihua Wen and Jinyou Xiao. Efficiency improvement of the polar
 * coordinate transformation for evaluating BEM singular integrals on curved
 * elements. Engineering Analysis with Boundary Elements. 83-93, 38, 2014.
 **/

static const double X_6[6] =
  {
    0.619309593042,
    0.830604693233,
    0.966234757102,
    0.380690406958,
    0.169395306767,
    0.033765242898
  };

static const double W_6[6] =
  {
    0.233956967286,
    0.180380786524,
    0.085662246190,
    0.233956967286,
    0.180380786524,
    0.085662246190
  };

static const double X_8[8] =
  {
    0.591717321248,
    0.762766204958,
    0.898333238707,
    0.980144928249,
    0.408282678752,
    0.237233795042,
    0.101666761293,
    0.019855071751
  };

static const double W_8[8] =
  {
    0.181341891689,
    0.156853322939,
    0.111190517227,
    0.050614268145,
    0.181341891689,
    0.156853322939,
    0.111190517227,
    0.050614268145
  };

static const double X_10[10] =
  {
    0.574437169491,
    0.716697697065,
    0.839704784150,
    0.932531683344,
    0.986953264259,
    0.425562830509,
    0.283302302935,
    0.160295215850,
    0.067468316656,
    0.013046735741
  };

static const double W_10[10] =
  {
    0.147762112357,
    0.134633359655,
    0.109543181258,
    0.074725674575,
    0.033335672154,
    0.147762112357,
    0.134633359655,
    0.109543181258,
    0.074725674575,
    0.033335672154
  };

static const double X_12[12] =
  {
    0.562616704256,
    0.683915749499,
    0.793658977143,
    0.884951337097,
    0.952058628185,
    0.990780317123,
    0.437383295744,
    0.316084250501,
    0.206341022857,
    0.115048662903,
    0.047941371815,
    0.009219682877
  };

static const double W_12[12] =
  {
    0.124573522907,
    0.116746268269,
    0.101583713362,
    0.080039164272,
    0.053469662998,
    0.023587668193,
    0.124573522907,
    0.116746268269,
    0.101583713362,
    0.080039164272,
    0.053469662998,
    0.023587668193
  };

static const double X_14[14] =
  {
    0.554027474354,
    0.659556184464,
    0.757624318179,
    0.843646452406,
    0.913600657535,
    0.964217441832,
    0.993141904348,
    0.445972525646,
    0.340443815536,
    0.242375681821,
    0.156353547594,
    0.086399342465,
    0.035782558168,
    0.006858095652
  };

static const double W_14[14] =
  {
    0.107631926732,
    0.102599231861,
    0.092769198739,
    0.078601583579,
    0.060759285344,
    0.040079043580,
    0.017559730166,
    0.107631926732,
    0.102599231861,
    0.092769198739,
    0.078601583579,
    0.060759285344,
    0.040079043580,
    0.017559730166
  };

static const double X_16[16] =
  {
    0.547506254919,
    0.640801775390,
    0.729008388829,
    0.808938122201,
    0.877702204178,
    0.932815601194,
    0.972287511537,
    0.994700467496,
    0.452493745081,
    0.359198224610,
    0.270991611171,
    0.191061877799,
    0.122297795822,
    0.067184398806,
    0.027712488463,
    0.005299532504
  };

static const double W_16[16] =
  {
    0.094725305228,
    0.091301707522,
    0.084578259698,
    0.074797994408,
    0.062314485628,
    0.047579255841,
    0.031126761969,
    0.013576229706,
    0.094725305228,
    0.091301707522,
    0.084578259698,
    0.074797994408,
    0.062314485628,
    0.047579255841,
    0.031126761969,
    0.013576229706
  };

static const double X_18[18] =
  {
    0.542387506521,
    0.625943112846,
    0.705875580731,
    0.779885415537,
    0.845843521530,
    0.901852479486,
    0.946301233249,
    0.977911974786,
    0.995782584210,
    0.457612493479,
    0.374056887154,
    0.294124419269,
    0.220114584463,
    0.154156478470,
    0.098147520514,
    0.053698766751,
    0.022088025214,
    0.004217415790
  };

static const double W_18[18] =
  {
    0.084571191482,
    0.082138241873,
    0.077342337563,
    0.070321457335,
    0.061277603356,
    0.050471022053,
    0.038212865127,
    0.024857274447,
    0.010808006763,
    0.084571191482,
    0.082138241873,
    0.077342337563,
    0.070321457335,
    0.061277603356,
    0.050471022053,
    0.038212865127,
    0.024857274447,
    0.010808006763
  };

static const double X_20[20] =
  {
    0.538263260567,
    0.613892925571,
    0.686853044358,
    0.755433500975,
    0.818026840363,
    0.873165953230,
    0.919558485911,
    0.956117214126,
    0.981985963639,
    0.996564299593,
    0.461736739433,
    0.386107074429,
    0.313146955642,
    0.244566499025,
    0.181973159637,
    0.126834046770,
    0.080441514089,
    0.043882785874,
    0.018014036361,
    0.003435700407
  };

static const double W_20[20] =
  {
    0.076376693565,
    0.074586493236,
    0.071048054659,
    0.065844319225,
    0.059097265981,
    0.050965059909,
    0.041638370788,
    0.031336024167,
    0.020300714900,
    0.008807003570,
    0.076376693565,
    0.074586493236,
    0.071048054659,
    0.065844319225,
    0.059097265981,
    0.050965059909,
    0.041638370788,
    0.031336024167,
    0.020300714900,
    0.008807003570
  };

static const double X_132[132] =
  {
    0.505927352263,
    0.517778724764,
    0.529620103059,
    0.541444830589,
    0.553246260153,
    0.565017757647,
    0.576752705793,
    0.588444507860,
    0.600086591370,
    0.611672411795,
    0.623195456235,
    0.634649247078,
    0.646027345644,
    0.657323355802,
    0.668530927567,
    0.679643760668,
    0.690655608093,
    0.701560279596,
    0.712351645180,
    0.723023638545,
    0.733570260491,
    0.743985582297,
    0.754263749053,
    0.764398982947,
    0.774385586518,
    0.784217945856,
    0.793890533758,
    0.803397912836,
    0.812734738573,
    0.821895762329,
    0.830875834287,
    0.839669906354,
    0.848273034994,
    0.856680384010,
    0.864887227260,
    0.872888951315,
    0.880681058055,
    0.888259167193,
    0.895619018739,
    0.902756475397,
    0.909667524889,
    0.916348282209,
    0.922794991810,
    0.929004029712,
    0.934971905542,
    0.940695264496,
    0.946170889222,
    0.951395701631,
    0.956366764629,
    0.961081283763,
    0.965536608799,
    0.969730235207,
    0.973659805569,
    0.977323110911,
    0.980718091940,
    0.983842840205,
    0.986695599177,
    0.989274765240,
    0.991578888607,
    0.993606674167,
    0.995356982274,
    0.996828829564,
    0.998021390017,
    0.998933997355,
    0.999566155774,
    0.999917650413,
    0.494072647737,
    0.482221275236,
    0.470379896941,
    0.458555169411,
    0.446753739847,
    0.434982242353,
    0.423247294207,
    0.411555492140,
    0.399913408630,
    0.388327588205,
    0.376804543765,
    0.365350752922,
    0.353972654356,
    0.342676644198,
    0.331469072433,
    0.320356239332,
    0.309344391907,
    0.298439720404,
    0.287648354820,
    0.276976361455,
    0.266429739509,
    0.256014417703,
    0.245736250947,
    0.235601017053,
    0.225614413482,
    0.215782054144,
    0.206109466242,
    0.196602087164,
    0.187265261427,
    0.178104237671,
    0.169124165713,
    0.160330093646,
    0.151726965006,
    0.143319615990,
    0.135112772740,
    0.127111048685,
    0.119318941945,
    0.111740832807,
    0.104380981261,
    0.097243524603,
    0.090332475111,
    0.083651717791,
    0.077205008190,
    0.070995970288,
    0.065028094458,
    0.059304735504,
    0.053829110778,
    0.048604298369,
    0.043633235371,
    0.038918716237,
    0.034463391201,
    0.030269764793,
    0.026340194431,
    0.022676889089,
    0.019281908060,
    0.016157159795,
    0.013304400823,
    0.010725234760,
    0.008421111393,
    0.006393325833,
    0.004643017726,
    0.003171170436,
    0.001978609983,
    0.001066002645,
    0.000433844226,
    0.000082349587
  };

static const double W_132[132] =
  {
    0.011854149158,
    0.011847485418,
    0.011834161684,
    0.011814185447,
    0.011787567935,
    0.011754324111,
    0.011714472664,
    0.011668035995,
    0.011615040209,
    0.011555515097,
    0.011489494120,
    0.011417014393,
    0.011338116658,
    0.011252845268,
    0.011161248158,
    0.011063376818,
    0.010959286266,
    0.010849035016,
    0.010732685046,
    0.010610301760,
    0.010481953956,
    0.010347713784,
    0.010207656705,
    0.010061861453,
    0.009910409985,
    0.009753387438,
    0.009590882083,
    0.009422985270,
    0.009249791382,
    0.009071397778,
    0.008887904742,
    0.008699415423,
    0.008506035779,
    0.008307874518,
    0.008105043034,
    0.007897655348,
    0.007685828043,
    0.007469680194,
    0.007249333310,
    0.007024911256,
    0.006796540190,
    0.006564348490,
    0.006328466681,
    0.006089027362,
    0.005846165133,
    0.005600016519,
    0.005350719889,
    0.005098415386,
    0.004843244841,
    0.004585351697,
    0.004324880929,
    0.004061978960,
    0.003796793581,
    0.003529473868,
    0.003260170096,
    0.002989033663,
    0.002716216997,
    0.002441873484,
    0.002166157386,
    0.001889223779,
    0.001611228525,
    0.001332328354,
    0.001052681395,
    0.000772449927,
    0.000491820334,
    0.000211329833,
    0.011854149158,
    0.011847485418,
    0.011834161684,
    0.011814185447,
    0.011787567935,
    0.011754324111,
    0.011714472664,
    0.011668035995,
    0.011615040209,
    0.011555515097,
    0.011489494120,
    0.011417014393,
    0.011338116658,
    0.011252845268,
    0.011161248158,
    0.011063376818,
    0.010959286266,
    0.010849035016,
    0.010732685046,
    0.010610301760,
    0.010481953956,
    0.010347713784,
    0.010207656705,
    0.010061861453,
    0.009910409985,
    0.009753387438,
    0.009590882083,
    0.009422985270,
    0.009249791382,
    0.009071397778,
    0.008887904742,
    0.008699415423,
    0.008506035779,
    0.008307874518,
    0.008105043034,
    0.007897655348,
    0.007685828043,
    0.007469680194,
    0.007249333310,
    0.007024911256,
    0.006796540190,
    0.006564348490,
    0.006328466681,
    0.006089027362,
    0.005846165133,
    0.005600016519,
    0.005350719889,
    0.005098415386,
    0.004843244841,
    0.004585351697,
    0.004324880929,
    0.004061978960,
    0.003796793581,
    0.003529473868,
    0.003260170096,
    0.002989033663,
    0.002716216997,
    0.002441873484,
    0.002166157386,
    0.001889223779,
    0.001611228525,
    0.001332328354,
    0.001052681395,
    0.000772449927,
    0.000491820334,
    0.000211329833
  };

static const double X_60[60] =
  {
    0.512979886151,
    0.538904666975,
    0.564724567698,
    0.590369982437,
    0.615771775688,
    0.640861468712,
    0.665571424134,
    0.689835028288,
    0.713586870792,
    0.736762920881,
    0.759300700029,
    0.781139450377,
    0.802220298524,
    0.822486414245,
    0.841883163691,
    0.860358256678,
    0.877861887653,
    0.894346869966,
    0.909768763081,
    0.924085992393,
    0.937259961323,
    0.949255155405,
    0.960039238089,
    0.969583138058,
    0.977861127920,
    0.984850894383,
    0.990533600876,
    0.994893947611,
    0.997920262559,
    0.999605061614,
    0.487020113849,
    0.461095333025,
    0.435275432302,
    0.409630017563,
    0.384228224312,
    0.359138531288,
    0.334428575866,
    0.310164971712,
    0.286413129208,
    0.263237079119,
    0.240699299971,
    0.218860549623,
    0.197779701476,
    0.177513585755,
    0.158116836309,
    0.139641743322,
    0.122138112347,
    0.105653130034,
    0.090231236919,
    0.075914007607,
    0.062740038677,
    0.050744844595,
    0.039960761911,
    0.030416861942,
    0.022138872080,
    0.015149105617,
    0.009466399124,
    0.005106052389,
    0.002079737441,
    0.000394938386
  };

static const double W_60[60] =
  {
    0.025953938816,
    0.025883971587,
    0.025744225750,
    0.025535078035,
    0.025257092266,
    0.024911017845,
    0.024497787728,
    0.024018515910,
    0.023474494424,
    0.022867189858,
    0.022198239398,
    0.021469446418,
    0.020682775618,
    0.019840347726,
    0.018944433785,
    0.017997449026,
    0.017001946362,
    0.015960609510,
    0.014876245750,
    0.013751778375,
    0.012590238811,
    0.011394758472,
    0.010168560365,
    0.008914950507,
    0.007637309298,
    0.006339083238,
    0.005023778591,
    0.003694965582,
    0.002356364963,
    0.001013405984,
    0.025953938816,
    0.025883971587,
    0.025744225750,
    0.025535078035,
    0.025257092266,
    0.024911017845,
    0.024497787728,
    0.024018515910,
    0.023474494424,
    0.022867189858,
    0.022198239398,
    0.021469446418,
    0.020682775618,
    0.019840347726,
    0.018944433785,
    0.017997449026,
    0.017001946362,
    0.015960609510,
    0.014876245750,
    0.013751778375,
    0.012590238811,
    0.011394758472,
    0.010168560365,
    0.008914950507,
    0.007637309298,
    0.006339083238,
    0.005023778591,
    0.003694965582,
    0.002356364963,
    0.001013405984
  };

static const double X_30[30] =
  {
    0.525735921278,
    0.576934956804,
    0.627318463084,
    0.676352362765,
    0.723516884769,
    0.768312074071,
    0.810263091495,
    0.848925247397,
    0.883888716052,
    0.914782881191,
    0.941280267896,
    0.963100023715,
    0.980010932484,
    0.991834061640,
    0.998446742037,
    0.474264078722,
    0.423065043196,
    0.372681536916,
    0.323647637235,
    0.276483115231,
    0.231687925929,
    0.189736908505,
    0.151074752603,
    0.116111283948,
    0.085217118809,
    0.058719732104,
    0.036899976285,
    0.019989067516,
    0.008165938360,
    0.001553257963
  };

static const double W_30[30] =
  {
    0.051426326447,
    0.050881194874,
    0.049796710293,
    0.048184368587,
    0.046061261119,
    0.043449893601,
    0.040377947615,
    0.036877987369,
    0.032987114941,
    0.028746578109,
    0.024201336415,
    0.019399596285,
    0.014392353942,
    0.009233234156,
    0.003984096248,
    0.051426326447,
    0.050881194874,
    0.049796710293,
    0.048184368587,
    0.046061261119,
    0.043449893601,
    0.040377947615,
    0.036877987369,
    0.032987114941,
    0.028746578109,
    0.024201336415,
    0.019399596285,
    0.014392353942,
    0.009233234156,
    0.003984096248
  };

static const double X_40[40] =
  {
    0.519386208753,
    0.558042035338,
    0.596348790351,
    0.634076092504,
    0.670997045413,
    0.706889602186,
    0.741537900843,
    0.774733562548,
    0.806276944834,
    0.835978342307,
    0.863659127595,
    0.889152825713,
    0.912306115417,
    0.932979751606,
    0.951049403484,
    0.966406404139,
    0.978958409607,
    0.988629974992,
    0.995363119350,
    0.999118854855,
    0.480613791247,
    0.441957964662,
    0.403651209649,
    0.365923907496,
    0.329002954587,
    0.293110397814,
    0.258462099157,
    0.225266437452,
    0.193723055166,
    0.164021657693,
    0.136340872405,
    0.110847174287,
    0.087693884583,
    0.067020248394,
    0.048950596516,
    0.033593595861,
    0.021041590393,
    0.011370025008,
    0.004636880650,
    0.000881145145
  };

static const double W_40[40] =
  {
    0.038752973989,
    0.038519909082,
    0.038055180950,
    0.037361584529,
    0.036443291198,
    0.035305823696,
    0.033956022908,
    0.032402006728,
    0.030653121246,
    0.028719884550,
    0.026613923492,
    0.024347903818,
    0.021935454093,
    0.019391083987,
    0.016730097641,
    0.013968503490,
    0.011122924597,
    0.008210529191,
    0.005249142266,
    0.002260638549,
    0.038752973989,
    0.038519909082,
    0.038055180950,
    0.037361584529,
    0.036443291198,
    0.035305823696,
    0.033956022908,
    0.032402006728,
    0.030653121246,
    0.028719884550,
    0.026613923492,
    0.024347903818,
    0.021935454093,
    0.019391083987,
    0.016730097641,
    0.013968503490,
    0.011122924597,
    0.008210529191,
    0.005249142266,
    0.002260638549
  };

static const double X_50[50] =
  {
    0.515549169164,
    0.546587350780,
    0.577445294999,
    0.608003618438,
    0.638144096890,
    0.667750122710,
    0.696707155949,
    0.724903167487,
    0.752229072454,
    0.778579152257,
    0.803851463592,
    0.827948232843,
    0.850776234353,
    0.872247151113,
    0.892277916450,
    0.910791035430,
    0.927714884715,
    0.942983989762,
    0.956539278328,
    0.968328309472,
    0.978305477621,
    0.986432192553,
    0.992677042024,
    0.997015984716,
    0.999433202210,
    0.484450830836,
    0.453412649220,
    0.422554705001,
    0.391996381562,
    0.361855903110,
    0.332249877290,
    0.303292844051,
    0.275096832513,
    0.247770927546,
    0.221420847743,
    0.196148536408,
    0.172051767157,
    0.149223765647,
    0.127752848887,
    0.107722083550,
    0.089208964570,
    0.072285115285,
    0.057016010238,
    0.043460721672,
    0.031671690528,
    0.021694522379,
    0.013567807447,
    0.007322957976,
    0.002984015284,
    0.000566797790
  };

static const double W_50[50] =
  {
    0.031088308328,
    0.030968033710,
    0.030727949795,
    0.030368985421,
    0.029892529352,
    0.029300424907,
    0.028594962824,
    0.027778872403,
    0.026855310944,
    0.025827851535,
    0.024700469225,
    0.023477525652,
    0.022163752169,
    0.020764231545,
    0.019284378306,
    0.017729917808,
    0.016106864112,
    0.014421496790,
    0.012680336785,
    0.010890121585,
    0.009057780357,
    0.007190411381,
    0.005295274192,
    0.003379899598,
    0.001454311277,
    0.031088308328,
    0.030968033710,
    0.030727949795,
    0.030368985421,
    0.029892529352,
    0.029300424907,
    0.028594962824,
    0.027778872403,
    0.026855310944,
    0.025827851535,
    0.024700469225,
    0.023477525652,
    0.022163752169,
    0.020764231545,
    0.019284378306,
    0.017729917808,
    0.016106864112,
    0.014421496790,
    0.012680336785,
    0.010890121585,
    0.009057780357,
    0.007190411381,
    0.005295274192,
    0.003379899598,
    0.001454311277
  };

//!calculate the distance between two points (two-dimensional)
double dstc2(double *p, double *q) {
	return sqrt((p[0]-q[0])*(p[0]-q[0])+(p[1]-q[1])*(p[1]-q[1]));
}

//!calculate the dot product between two vector (p,q0) an (p,q1) (two-dimensional)
double Dot_Pdt2(double *p, double *q0, double *q1) {
	return (q0[0]-p[0])*(q1[0]-p[0]) + (q0[1]-p[1])*(q1[1]-p[1]);
}

void rong_self_influence_coefficients (SPPanel * spp,
				       gdouble up, gdouble vp,
				       Point p,
				       gdouble * psi,
				       gdouble * phi)
{


  Spline2D * sp = spp->sp;
  gint i, j, m, n, kk;
  
  int SingAngN = 20; // Number of Gauss points in the angular direction
  int SingRadN = 10; // Number of Gauss points in the radial direction

  gdouble du = (spp->u1-spp->u0); // Mapping on [0:1]x[0:1] square
  gdouble dv = (spp->v1-spp->v0);

  gint k = spp->k, NU = sp->NXU;

  gsl_vector * Bu = gsl_vector_alloc (k);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  
  size_t istart, iend, jstart, jend;
  size_t istart_x, iend_x;
  
  // Metric and normal at singularity
  Point ps = p;
  Vector xu = spline2d_xu (sp, up, vp);
  Vector xv = spline2d_xv (sp, up, vp);
  Vector Ns = vector_vector_product (&xu, &xv);
  gdouble Js_inv = 1./vector_norm (Ns);
  Ns.x *= Js_inv; Ns.y *= Js_inv; Ns.z *= Js_inv;

  // Base for conformal mapping
  double u1, u2, u1dotu2, lambda_inv, sing, cosg;
  u1 = vector_norm (xu);
  u2 = vector_norm (xv);
  u1dotu2 = vector_scalar_product (&xu, &xv);
  lambda_inv = u2/u1;
  cosg = u1dotu2/(u1*u2);
  sing = sqrt(1.0-cosg*cosg);

  // Coordinate of element corners in (eta1-eta2) coordinate system
  double etav[2*4] /* = {0.0} */;
  etav[0] = 0.; etav[1] = 0.;
  etav[2] = cosg*lambda_inv; etav[3] = sing*lambda_inv;
  etav[4] = 1 + etav[2]; etav[5] = etav[3];
  etav[6] = 1.; etav[7] = 0.;

  // Coordinate of singularity in (eta1-eta2) coordinate system
  double etas[2];
  etas[0] = (up-spp->u0)/du + etav[2]*(vp-spp->v0)/dv;
  etas[1] = etav[3]*(vp-spp->v0)/dv;

  //calculate the distance from singular point to corners
  double d[4], ddv[4];
  ddv[0] = dstc2(etas,etav+6);
  ddv[1] = dstc2(etas,etav+4);
  ddv[2] = dstc2(etas,etav+2);
  ddv[3] = dstc2(etas,etav+0);
  d[0] = dstc2(etav+2*2,etav+3*2);
  d[1] = dstc2(etav+1*2,etav+2*2);
  d[2] = dstc2(etav,etav+1*2);
  d[3] = 1.0;

  //d(s)/d(eta), transformation matrix T as in equation (27) reverse transformation
  double dse[4];
  dse[0] = du;
  dse[1] = 0.0;
  dse[2] = -etav[2]*du/etav[3];
  dse[3] = dv/etav[3];

  //jacobian cofficent of the conformal transformation, i.e. normal of transformation matrix
  double J2 = dse[3]*dse[0];

  //perform integral
  //------------------------------------------------------------------------------

  double ab0, ab1;
  double theta, rho; //quadrature points in polar coordinates
  double s, t;	   //quadrature points
  double eta0, eta1; //quadrature points in eta0-eta1 system
  double rr2 = 0;    //up bound of rho
  double cta, sta;

  //calculate the angles
  double theta0;
  double dp = Dot_Pdt2(etav+6,etas,etav); // (etav+6,etas) an (etav+6,etav)
  theta0 = -acos(dp/ddv[0]);
  double SumTheta = 0.0;

  int no;//index of sub-triangle
  for (no = 0; no < 4; no++) {//for each sub-triangles
    
    //calculate the angles, i.e. deltatheta, alpha and h
    double cgmma, deltatheta, alpha, cost, h;
    if (no != 3) {
      kk = 3-no; // sideno - no - 1;
      dp = Dot_Pdt2(etas,etav+kk*2,etav+(kk-1)*2);
      cost = dp/(ddv[no]*ddv[no+1]);
      deltatheta = acos(cost); // Inner angle

      dp = Dot_Pdt2(etav+kk*2,etas,etav+(kk-1)*2);
      cgmma = dp/(ddv[no]*d[no]);
      alpha = 0.5*M_PI-acos(cgmma);
      h = ddv[no]*cos(alpha); // Altitude of triangle
      SumTheta += deltatheta;
    }
    else {
      deltatheta = 2.0*M_PI-SumTheta; // Inner angle
      alpha = 1.5*M_PI-theta0;
      h = etas[1]; // Altitude of triangle
    }

    // Sigmoidal transformation, change to new variable z.
    // The up and lower bound of z is calculated here
    double thetabar0, thetabar1, z, z0, z1, powt;

    thetabar0 = -alpha/M_PI+0.5; // lower end of angular range of integration + transformation (33)
    thetabar1 = (deltatheta-alpha)/M_PI+0.5; // upper end of angular range of integration + transformation (33)


    powt = pow(thetabar0/(1.0-thetabar0),1.0/3.0/* m */); // Looks close to sigmoidal transfo (31)
    z0 = powt/(1.0+powt);
    powt = pow(thetabar1/(1.0-thetabar1),1.0/3.0/* m */);
    z1 = powt/(1.0+powt);

    ab0 = z1 - z0; // Some kind of dz
    ab1 = z0; // Start of interval


    for (i = 0; i < SingAngN; i++) { // Loop over angular quadrature points

      z = ab1 + ab0*X_20[i];//integral variable z

      //do sigmoidal transformation
      double Jc, powt2, pt, temp2;
      powt = pow(z,3.-1); // Here 3 is optimum value for this type of singularity
      powt2 = pow((1.0-z),3.-1);
      temp2 = powt*z+powt2*(1.0-z);
      pt = powt*z/temp2; // (31)

      theta = M_PI*(pt-0.5) + alpha; // Should not be negative

      // upper bound of rho for angle theta
      rr2 = h/cos(theta-alpha);

      //d(theta)/dz Jacobian of sigmoidal transformation times 3 other Jacobian/Gauss-Weights coefficients
      Jc = M_PI*3.*powt*powt2/(temp2*temp2)*ab0*rr2*W_20[i];


      theta += theta0; //Radial coordinate of integration point

      cta = cos(theta);
      sta = sin(theta);

      //compute the line integral of the equation # 27 in M. Guiggiani's paper

      //compute the double integral
      for (j = 0; j < SingRadN; j++) { // Loop over radial quadrature points

	rho = rr2*X_10[j]; // Radial coordinate of intagration point in (eta1, eta2)

	//source point (gauss integral point) in the parameter coordinate
	eta0 = etas[0] + rho*cta;
	eta1 = etas[1] + rho*sta;

	gdouble u = dse[0]*eta0 + dse[2]*eta1 + spp->u0; //conformal transformation (reverse)
	gdouble v = /* dse[1]*eta0 + */ dse[3]*eta1 + spp->v0;

	gsl_bspline_eval_nonzero (u, Bu, &istart, &iend, sp->w_u);
	gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
	gsl_bspline_deriv_eval_nonzero (u, 1, Bux, &istart_x, &iend_x, sp->wx_u, sp->wxd_u);

	Point pg;
	gint a, b;

	xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
	pg.x = pg.y = pg.z = 0.;

	gint ii = istart_x;
	for ( a = 0; a < k; a++) {
	  gdouble cu = gsl_matrix_get (Bux, a, 0);
	  gdouble cdu = gsl_matrix_get (Bux, a, 1);
	  gint jj = jstart;
	  for ( b = 0; b < k; b++) {
	    gdouble cv = gsl_matrix_get (Bv, b, 0);
	    gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	    gdouble cvdu = cv*cdu;
	    gdouble cuv = cu*cv;

	    SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, ii+ jj*NU);
	    gdouble v0 = sc->v[0];
	    gdouble v1 = sc->v[1];
	    gdouble v2 = sc->v[2];

	    xu.x += v0*cvdu;
	    xu.y += v1*cvdu;
	    xu.z += v2*cvdu;
	    xv.x += v0*cudv;
	    xv.y += v1*cudv;
	    xv.z += v2*cudv;
	    pg.x += v0*cuv;
	    pg.y += v1*cuv;
	    pg.z += v2*cuv;
	    jj++;
	  }
	  ii++;
	}

	// Normal at Gauss-Legendre point
	Vector N = vector_vector_product (&xu, &xv);

	// Physical distance to Gauss-Legendre Point
	Vector R;
	R.x = pg.x - p.x;
	R.y = pg.y - p.y;
	R.z = pg.z - p.z;
	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
	// Jacobian at Gauss-Legendre point = norm of N
	gdouble J = vector_norm (N);
	gdouble c1 = Jc*rho*W_10[j]/sqrt(r2);
	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

#if 1
	R.z = pg.z+p.z+2.*(15.45);
	r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	gdouble c1_tmp = Jc*rho*W_10[j]/sqrt(r2);
	c2 += c1_tmp*vector_scalar_product (&N, &R)/r2;
	c1 += c1_tmp;
#endif

	c1 *= J;

	// Loop over the splines included in spp
	gint mm = 0;
	for ( n = 0; n < k; n++) {
	  gdouble c3 = c1*gsl_matrix_get (Bv, n, 0);
	  gdouble c4 = c2*gsl_matrix_get (Bv, n, 0); 
	  for ( m = 0; m < k; m++) {
	    psi[mm] += c3*gsl_vector_get (Bu, m);
	    phi[mm] += c4*gsl_vector_get (Bu, m);
	    mm++;
	  }
	}

      }

    }
    theta0 += deltatheta;
  }
  

  for ( m = 0; m < k*k; m++) {
    psi[m] *= J2;
    phi[m] *= J2;
  }

  gsl_vector_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

}

void  wamit_self_influence_coefficients (SPPanel * spp,
					 gdouble up, gdouble vp,
					 Point p,
					 gdouble * psi,
					 gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 40, i, j, m, n, se;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble u, v, uu, vv, J;
  gint k = spp->k;
  gdouble A = 0.881373587019543; // asinh (1)
  
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;
  gdouble J2;

  gdouble ui[ng], wi[ng];
  for ( i = 0; i < ng; i++)
    gsl_integration_glfixed_point (-1., 1., i, &ui[i], &wi[i], itable);
  


  for ( i = 0; i < ng; i++) {
    for ( j = 0; j < ng; j++) {

      /* Term 1 */
      uu = ui[i];
      vv = ui[i]*sinh(A*ui[j]);
      
      if ( uu < 0 ) {
	if ( vv < 0 ) {
	  u = spp->u0 + (up-spp->u0)*(1.+uu);
	  v = spp->v0 + (vp-spp->v0)*(1.+vv);
	  J2 = (up-spp->u0)*(vp-spp->v0);
	}
	else {
	  u = spp->u0 + (up-spp->u0)*(1.+uu);
	  v = vp + (spp->v1-vp)*vv;
	  J2 = (up-spp->u0)*(spp->v1-vp);
	}
      }
      else {
	if ( vv < 0 ) {
	  u = up + (spp->u1-up)*uu;
	  v = spp->v0 + (vp-spp->v0)*(1.+vv);
	  J2 = (spp->u1-up)*(vp-spp->v0);
	}
	else {
	  u = up + (spp->u1-up)*uu;
	  v = vp + (spp->v1-vp)*vv;
	  J2 = (spp->u1-up)*(spp->v1-vp);
	}
      }


      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);


      Point pg;
      Vector xu, xv;
      gint  a, b;

      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      pg.x = pg.y = pg.z = 0.;

      for ( a = 0; a < k; a++) {
	gdouble cu = gsl_matrix_get (Bux, a, 0);
	gdouble cdu = gsl_matrix_get (Bux, a, 1);
	gint ii = (istart+a);
	for ( b = 0; b < k; b++) {
	  gdouble cv = gsl_matrix_get (Bv, b, 0);
	  gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	  gdouble cvdu = cv*cdu;
	  gdouble cuv = cu*cv;
	  gint jj = (jstart+b);
	  gdouble v0 = coeff (sp, ii, jj, 0);
	  gdouble v1 = coeff (sp, ii, jj, 1);
	  gdouble v2 = coeff (sp, ii, jj, 2);

	  xu.x += v0*cvdu;
	  xu.y += v1*cvdu;
	  xu.z += v2*cvdu;
	  xv.x += v0*cudv;
	  xv.y += v1*cudv;
	  xv.z += v2*cudv;
	  pg.x += v0*cuv;
	  pg.y += v1*cuv;
	  pg.z += v2*cuv;
	}
      }
	
      // Normal at Gauss-Legendre point
      Vector N = vector_vector_product (&xu, &xv);

      // Physical distance to Gauss-Legendre Point
      Vector R;
      R.x = pg.x - p.x;
      R.y = pg.y - p.y;
      R.z = pg.z - p.z;
      gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
      // Jacobian at Gauss-Legendre point = norm of N
      gdouble J = vector_norm (N);
      gdouble c1 = wi[i]*wi[j]*A*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
      N.x /= J; N.y /= J; N.z /= J;
      gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

      // Loop over the splines included in spp
      for ( m = 0; m < k; m++) {
	gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	for ( n = 0; n < k; n++) {
	  psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	  phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	}
      }


      /* Term 2 */
      uu = vv;
      vv = ui[i];

      if ( uu < 0 ) {
	if ( vv < 0 ) {
	  u = spp->u0 + (up-spp->u0)*(1.+uu);
	  v = spp->v0 + (vp-spp->v0)*(1.+vv);
	  J2 = (up-spp->u0)*(vp-spp->v0);
	}
	else {
	  u = spp->u0 + (up-spp->u0)*(1.+uu);
	  v = vp + (spp->v1-vp)*vv;
	  J2 = (up-spp->u0)*(spp->v1-vp);
	}
      }
      else {
	if ( vv < 0 ) {
	  u = up + (spp->u1-up)*uu;
	  v = spp->v0 + (vp-spp->v0)*(1.+vv);
	  J2 = (spp->u1-up)*(vp-spp->v0);
	}
	else {
	  u = up + (spp->u1-up)*uu;
	  v = vp + (spp->v1-vp)*vv;
	  J2 = (spp->u1-up)*(spp->v1-vp);
	}
      }

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);
      
      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      pg.x = pg.y = pg.z = 0.;

      for ( a = 0; a < k; a++) {
	gdouble cu = gsl_matrix_get (Bux, a, 0);
	gdouble cdu = gsl_matrix_get (Bux, a, 1);
	gint ii = (istart+a);
	for ( b = 0; b < k; b++) {
	  gdouble cv = gsl_matrix_get (Bv, b, 0);
	  gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	  gdouble cvdu = cv*cdu;
	  gdouble cuv = cu*cv;
	  gint jj = (jstart+b);
	  gdouble v0 = coeff (sp, ii, jj, 0);
	  gdouble v1 = coeff (sp, ii, jj, 1);
	  gdouble v2 = coeff (sp, ii, jj, 2);

	  xu.x += v0*cvdu;
	  xu.y += v1*cvdu;
	  xu.z += v2*cvdu;
	  xv.x += v0*cudv;
	  xv.y += v1*cudv;
	  xv.z += v2*cudv;
	  pg.x += v0*cuv;
	  pg.y += v1*cuv;
	  pg.z += v2*cuv;
	}
      }
	
      // Normal at Gauss-Legendre point
      N = vector_vector_product (&xu, &xv);

      // Physical distance to Gauss-Legendre Point
      R.x = pg.x - p.x;
      R.y = pg.y - p.y;
      R.z = pg.z - p.z;
      r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
      // Jacobian at Gauss-Legendre point = norm of N
      J = vector_norm (N);
      c1 = wi[i]*wi[j]*A*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
      N.x /= J; N.y /= J; N.z /= J;
      c2 = c1*vector_scalar_product (&N, &R)/r2;


      // Loop over the splines included in spp
      for ( m = 0; m < k; m++) {
	gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	for ( n = 0; n < k; n++) {
	  psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	  phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	}
      }
    }
  }
  
  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

  gsl_integration_glfixed_table_free (itable);
}

static void wamit_self_influence_point_contribution (SPPanel * spp,
						     gdouble up, gdouble vp,
						     Point p,
						     gdouble * psi,
						     gdouble * phi,
						     gdouble uu,
						     gdouble vv,
						     gdouble weight)
{
  Spline2D * sp = spp->sp;
  gdouble u, v, J2;
  gint a, b, m, n, k = spp->k;
  size_t istart, iend, jstart, jend;

  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);

  if ( uu < 0 ) {
    if ( vv < 0 ) {
      u = spp->u0 + (up-spp->u0)*(1.+uu);
      v = spp->v0 + (vp-spp->v0)*(1.+vv);
      J2 = (up-spp->u0)*(vp-spp->v0);
    }
    else {
      u = spp->u0 + (up-spp->u0)*(1.+uu);
      v = vp + (spp->v1-vp)*vv;
      J2 = (up-spp->u0)*(spp->v1-vp);
    }
  }
  else {
    if ( vv < 0 ) {
      u = up + (spp->u1-up)*uu;
      v = spp->v0 + (vp-spp->v0)*(1.+vv);
      J2 = (spp->u1-up)*(vp-spp->v0);
    }
    else {
      u = up + (spp->u1-up)*uu;
      v = vp + (spp->v1-vp)*vv;
      J2 = (spp->u1-up)*(spp->v1-vp);
    }
  }


  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
  
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);

  Point pg;
  Vector xu, xv;

  xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
  pg.x = pg.y = pg.z = 0.;

  for ( a = 0; a < k; a++) {
    gdouble cu = gsl_matrix_get (Bux, a, 0);
    gdouble cdu = gsl_matrix_get (Bux, a, 1);
    gint ii = (istart+a);
    for ( b = 0; b < k; b++) {
      gdouble cv = gsl_matrix_get (Bv, b, 0);
      gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
      gdouble cvdu = cv*cdu;
      gdouble cuv = cu*cv;
      gint jj = (jstart+b);
      gdouble v0 = coeff (sp, ii, jj, 0);
      gdouble v1 = coeff (sp, ii, jj, 1);
      gdouble v2 = coeff (sp, ii, jj, 2);

      xu.x += v0*cvdu;
      xu.y += v1*cvdu;
      xu.z += v2*cvdu;
      xv.x += v0*cudv;
      xv.y += v1*cudv;
      xv.z += v2*cudv;
      pg.x += v0*cuv;
      pg.y += v1*cuv;
      pg.z += v2*cuv;
    }
  }
	
  // Normal at Gauss-Legendre point
  Vector N = vector_vector_product (&xu, &xv);

  // Physical distance to Gauss-Legendre Point
  Vector R;
  R.x = pg.x - p.x;
  R.y = pg.y - p.y;
  R.z = pg.z - p.z;
  gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
  // Jacobian at Gauss-Legendre point = norm of N
  gdouble J = vector_norm (N);
  gdouble c1 = weight*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
  N.x /= J; N.y /= J; N.z /= J;
  gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

  // Loop over the splines included in spp
  for ( m = 0; m < k; m++) {
    gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
    gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
    for ( n = 0; n < k; n++) {
      psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
      phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
    }
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);
}

void  wamit_self_influence_coefficients_cauchy (SPPanel * spp,
						gdouble up, gdouble vp,
						Point p,
						gdouble * psi,
						gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 6, i, j;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble uu, vv;
  gint k = spp->k;
  gdouble A = 0.881373587019543; // asinh (1)

  gdouble J2;

  gdouble epsilon = 1e-4;

  gdouble ui[ng], wi[ng];
  for ( i = 0; i < ng; i++)
    gsl_integration_glfixed_point (epsilon, 1., i, &ui[i], &wi[i], itable);
  


  for ( i = 0; i < ng; i++) {
    for ( j = 0; j < ng; j++) {

      gdouble weight = A*wi[i]*wi[j];

      /* Term 1 */
      uu = ui[i];
      vv = ui[i]*sinh(A*ui[j]);
      
      
      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);

      /* Term 2 */
      uu = vv;
      vv = ui[i];

      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);


      /* Term 1 */
      uu = ui[i];
      vv = ui[i]*sinh(A*-ui[j]);
      
      
      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);

      /* Term 2 */
      uu = vv;
      vv = ui[i];

      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);


      
      /* Term 1 */
      uu = -ui[i];
      vv = -ui[i]*sinh(A*ui[j]);
      
      
      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);

      /* Term 2 */
      uu = vv;
      vv = -ui[i];

      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);


      /* Term 1 */
      uu = -ui[i];
      vv = -ui[i]*sinh(A*-ui[j]);
      
      
      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);

      /* Term 2 */
      uu = vv;
      vv = -ui[i];

      wamit_self_influence_point_contribution (spp, up, vp, p, psi, phi, uu, vv, weight);
    }
  }

  gsl_integration_glfixed_table_free (itable);
}

static void wamit_rectangle_integral (SPPanel * spp,
				      gdouble u0, gdouble u1,
				      gdouble v0, gdouble v1,
				      gdouble up, gdouble vp,
				      Point p,
				      gdouble * psi,
				      gdouble * phi)
{
  GPCell * panel_tree = gpcell_new (spp->sp, u0, u1,
				    v0, v1,
				    6, 0);

  spline_near_field_influence_coeff_recursive (spp, panel_tree, p, psi, phi, FALSE);

  gpcell_tree_destroy (panel_tree);
}

static void wamit_singularity_integral (SPPanel * spp,
					gdouble u0, gdouble u1,
					gdouble v0, gdouble v1,
					gdouble up, gdouble vp,
					Point p,
					gdouble * psi,
					gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 8, i, j, m, n, se;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble u, v, uu, vv, J;
  gint k = spp->k;
  gdouble A = 0.881373587019543; // asinh (1)
  
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;
  gdouble du = (u1-u0)/2.;
  gdouble dv = (v1-v0)/2.;
  gdouble J2 = A*du*dv;
  

  gdouble ui[ng], wi[ng];
  for ( i = 0; i < ng; i++)
    gsl_integration_glfixed_point (-1., 1., i, &ui[i], &wi[i], itable);

  for ( i = 0; i < ng; i++) {
    for ( j = 0; j < ng; j++) {

      /* Term 1 */
      uu = ui[i];
      vv = ui[i]*sinh(A*ui[j]);

      if ( uu < 0 )
	u = u0 + du*(1.+uu);
      else
	u = up + du*uu;

      if ( vv < 0 )
	v = v0 + dv*(1.+vv);
      else
	v = vp + dv*vv;

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
	
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);

      Point pg;
      Vector xu, xv;
      gint  a, b;

      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      pg.x = pg.y = pg.z = 0.;

      for ( a = 0; a < k; a++) {
	gdouble cu = gsl_matrix_get (Bux, a, 0);
	gdouble cdu = gsl_matrix_get (Bux, a, 1);
	gint ii = (istart+a);
	for ( b = 0; b < k; b++) {
	  gdouble cv = gsl_matrix_get (Bv, b, 0);
	  gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	  gdouble cvdu = cv*cdu;
	  gdouble cuv = cu*cv;
	  gint jj = (jstart+b);
	  gdouble v0 = coeff (sp, ii, jj, 0);
	  gdouble v1 = coeff (sp, ii, jj, 1);
	  gdouble v2 = coeff (sp, ii, jj, 2);

	  xu.x += v0*cvdu;
	  xu.y += v1*cvdu;
	  xu.z += v2*cvdu;
	  xv.x += v0*cudv;
	  xv.y += v1*cudv;
	  xv.z += v2*cudv;
	  pg.x += v0*cuv;
	  pg.y += v1*cuv;
	  pg.z += v2*cuv;
	}
      }
	
      // Normal at Gauss-Legendre point
      Vector N = vector_vector_product (&xu, &xv);

      // Physical distance to Gauss-Legendre Point
      Vector R;
      R.x = pg.x - p.x;
      R.y = pg.y - p.y;
      R.z = pg.z - p.z;
      gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
      // Jacobian at Gauss-Legendre point = norm of N
      gdouble J = vector_norm (N);
      gdouble c1 = wi[i]*wi[j]*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
      N.x /= J; N.y /= J; N.z /= J;
      gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

      // Loop over the splines included in spp
      for ( m = 0; m < k; m++) {
	gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	for ( n = 0; n < k; n++) {
	  psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	  phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	}
      }


      /* Term 2 */
      uu = vv;
      vv = ui[i];

      if ( uu < 0 )
	u = u0 + du*(1.+uu);
      else
	u = up + du*uu;

      if ( vv < 0 )
	v = v0 + dv*(1.+vv);
      else
	v = vp + dv*vv;

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);
      

      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      pg.x = pg.y = pg.z = 0.;

      for ( a = 0; a < k; a++) {
	gdouble cu = gsl_matrix_get (Bux, a, 0);
	gdouble cdu = gsl_matrix_get (Bux, a, 1);
	gint ii = (istart+a);
	for ( b = 0; b < k; b++) {
	  gdouble cv = gsl_matrix_get (Bv, b, 0);
	  gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	  gdouble cvdu = cv*cdu;
	  gdouble cuv = cu*cv;
	  gint jj = (jstart+b);
	  gdouble v0 = coeff (sp, ii, jj, 0);
	  gdouble v1 = coeff (sp, ii, jj, 1);
	  gdouble v2 = coeff (sp, ii, jj, 2);

	  xu.x += v0*cvdu;
	  xu.y += v1*cvdu;
	  xu.z += v2*cvdu;
	  xv.x += v0*cudv;
	  xv.y += v1*cudv;
	  xv.z += v2*cudv;
	  pg.x += v0*cuv;
	  pg.y += v1*cuv;
	  pg.z += v2*cuv;
	}
      }
	
      // Normal at Gauss-Legendre point
      N = vector_vector_product (&xu, &xv);

      // Physical distance to Gauss-Legendre Point
      R.x = pg.x - p.x;
      R.y = pg.y - p.y;
      R.z = pg.z - p.z;
      r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
      // Jacobian at Gauss-Legendre point = norm of N
      J = vector_norm (N);
      c1 = wi[i]*wi[j]*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
      N.x /= J; N.y /= J; N.z /= J;
      c2 = c1*vector_scalar_product (&N, &R)/r2;


      // Loop over the splines included in spp
      for ( m = 0; m < k; m++) {
	gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	for ( n = 0; n < k; n++) {
	  psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	  phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	}
      }
    }
  }
  
  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

  gsl_integration_glfixed_table_free (itable);
}

void  centered_wamit_self_influence_coefficients (SPPanel * spp,
						  gdouble up, gdouble vp,
						  Point p,
						  gdouble * psi,
						  gdouble * phi)
{


  if ( fabs (up-spp->u0) > fabs (spp->u1-up) ) {
    gdouble um = spp->u1 - 2.*(spp->u1-up);
    if ( fabs (vp-spp->v0) > fabs (spp->v1-vp) ) {    
      gdouble vm = spp->v1 - 2.*(spp->v1-vp);

      wamit_singularity_integral (spp, um, spp->u1, vm, spp->v1, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, spp->u0, um, vm, spp->v1, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, um, spp->u1, spp->v0, vm, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, spp->u0, um, spp->v0, vm, up, vp, p, psi, phi);
    }
    else {
      gdouble vm = spp->v0 + 2.*(vp-spp->v0);

      wamit_singularity_integral (spp, um, spp->u1, spp->v0, vm, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, spp->u0, um, spp->v0, vm, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, um, spp->u1, vm, spp->v1, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, spp->u0, um, vm, spp->v1, up, vp, p, psi, phi);
    }
  }
  else {
    gdouble um = spp->u0 + 2.*(up-spp->u0);
    if ( fabs (vp-spp->v0) > fabs (spp->v1-vp) ) {
      gdouble vm = spp->v1 - 2.*(spp->v1-vp);

      wamit_singularity_integral (spp, spp->u0, um, vm, spp->v1, up, vp, p, psi, phi);  
      wamit_rectangle_integral (spp, um, spp->u1, vm, spp->v1, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, spp->u0, um, spp->v0, vm, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, um, spp->u1, spp->v0, vm, up, vp, p, psi, phi);
      
    }
    else {
      gdouble vm = spp->v0 + 2.*(vp-spp->v0);

      wamit_singularity_integral (spp, spp->u0, um, spp->v0, vm, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, um, spp->u1, spp->v0, vm, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, spp->u0, um, vm, spp->v1, up, vp, p, psi, phi);
      wamit_rectangle_integral (spp, um, spp->u1, vm, spp->v1, up, vp, p, psi, phi);
    }
  }
}

void  wamit_self_influence_coefficients_adaptive (SPPanel * spp,
						  gdouble up, gdouble vp,
						  Point p,
						  gdouble * psi,
						  gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = 16, i, j, m, n, se;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui, wi, vj, wj, u, v, uu, vv, J;
  gint k = spp->k;
  gdouble A = 0.881373587019543; // asinh (1)
  
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;
  gdouble J2;

  for ( i = 0; i < ng; i++) {
    gsl_integration_glfixed_point (-1., 1., i, &ui, &wi, itable);

    for ( j = 0; j < ng; j++) {
      gsl_integration_glfixed_point (-1., 1., j, &vj, &wj, jtable);

      /* Term 1 */
      uu = ui;
      vv = ui*sinh(A*vj);

      if ( uu < 0 && vv < 0 ) {
	u = spp->u0 + (up-spp->u0)*(1.+uu);
	v = spp->v0 + (vp-spp->v0)*(1.+vv);
	J2 = (up-spp->u0)*(vp-spp->v0);
      }

      if ( uu > 0 && vv < 0 ) {
	u = up + (spp->u1-up)*uu;
	v = spp->v0 + (vp-spp->v0)*(1.+vv);
	J2 = (spp->u1-up)*(vp-spp->v0);
      }
      
      if ( uu < 0 && vv > 0) {
	u = spp->u0 + (up-spp->u0)*(1.+uu);
	v = vp + (spp->v1-vp)*vv;
	J2 = (up-spp->u0)*(spp->v1-vp);
      }

      if ( uu > 0 && vv > 0 ) {
	u = up + (spp->u1-up)*uu;
	v = vp + (spp->v1-vp)*vv;
	J2 = (spp->u1-up)*(spp->v1-vp);
      }
      
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);
      

      Point pg;
      Vector xu, xv;
      gint  a, b;

      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      pg.x = pg.y = pg.z = 0.;

      for ( a = 0; a < k; a++) {
	gdouble cu = gsl_matrix_get (Bux, a, 0);
	gdouble cdu = gsl_matrix_get (Bux, a, 1);
	gint ii = (istart+a);
	for ( b = 0; b < k; b++) {
	  gdouble cv = gsl_matrix_get (Bv, b, 0);
	  gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	  gdouble cvdu = cv*cdu;
	  gdouble cuv = cu*cv;
	  gint jj = (jstart+b);
	  gdouble v0 = coeff (sp, ii, jj, 0);
	  gdouble v1 = coeff (sp, ii, jj, 1);
	  gdouble v2 = coeff (sp, ii, jj, 2);

	  xu.x += v0*cvdu;
	  xu.y += v1*cvdu;
	  xu.z += v2*cvdu;
	  xv.x += v0*cudv;
	  xv.y += v1*cudv;
	  xv.z += v2*cudv;
	  pg.x += v0*cuv;
	  pg.y += v1*cuv;
	  pg.z += v2*cuv;
	}
      }
	
      // Normal at Gauss-Legendre point
      Vector N = vector_vector_product (&xu, &xv);

      // Physical distance to Gauss-Legendre Point
      Vector R;
      R.x = pg.x - p.x;
      R.y = pg.y - p.y;
      R.z = pg.z - p.z;
      gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
      // Jacobian at Gauss-Legendre point = norm of N
      gdouble J = vector_norm (N);
      gdouble c1 = wi*wj*A*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
      N.x /= J; N.y /= J; N.z /= J;
      gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

      // Loop over the splines included in spp
      for ( m = 0; m < k; m++) {
	gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	for ( n = 0; n < k; n++) {
	  psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	  phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	}
      }


      /* Term 2 */
      uu = ui*sinh(A*vj);
      vv = ui;

      if ( uu < 0 && vv < 0 ) {
	u = spp->u0 + (up-spp->u0)*(1.+uu);
	v = spp->v0 + (vp-spp->v0)*(1.+vv);
	J2 = (up-spp->u0)*(vp-spp->v0);
      }

      if ( uu > 0 && vv < 0 ) {
	u = up + (spp->u1-up)*uu;
	v = spp->v0 + (vp-spp->v0)*(1.+vv);
	J2 = (spp->u1-up)*(vp-spp->v0);
      }
      
      if ( uu < 0 && vv > 0) {
	u = spp->u0 + (up-spp->u0)*(1.+uu);
	v = vp + (spp->v1-vp)*vv;
	J2 = (up-spp->u0)*(spp->v1-vp);
      }

      if ( uu > 0 && vv > 0 ) {
	u = up + (spp->u1-up)*uu;
	v = vp + (spp->v1-vp)*vv;
	J2 = (spp->u1-up)*(spp->v1-vp);
      }

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
       gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);


      xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
      pg.x = pg.y = pg.z = 0.;

      for ( a = 0; a < k; a++) {
	gdouble cu = gsl_matrix_get (Bux, a, 0);
	gdouble cdu = gsl_matrix_get (Bux, a, 1);
	gint ii = (istart+a);
	for ( b = 0; b < k; b++) {
	  gdouble cv = gsl_matrix_get (Bv, b, 0);
	  gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	  gdouble cvdu = cv*cdu;
	  gdouble cuv = cu*cv;
	  gint jj = (jstart+b);
	  gdouble v0 = coeff (sp, ii, jj, 0);
	  gdouble v1 = coeff (sp, ii, jj, 1);
	  gdouble v2 = coeff (sp, ii, jj, 2);

	  xu.x += v0*cvdu;
	  xu.y += v1*cvdu;
	  xu.z += v2*cvdu;
	  xv.x += v0*cudv;
	  xv.y += v1*cudv;
	  xv.z += v2*cudv;
	  pg.x += v0*cuv;
	  pg.y += v1*cuv;
	  pg.z += v2*cuv;
	}
      }
	
      // Normal at Gauss-Legendre point
      N = vector_vector_product (&xu, &xv);

      // Physical distance to Gauss-Legendre Point
      R.x = pg.x - p.x;
      R.y = pg.y - p.y;
      R.z = pg.z - p.z;
      r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
      // Jacobian at Gauss-Legendre point = norm of N
      J = vector_norm (N);
      c1 = wi*wj*A*J2*sqrt(uu*uu+vv*vv)*J/sqrt(r2);
      N.x /= J; N.y /= J; N.z /= J;
      c2 = c1*vector_scalar_product (&N, &R)/r2;


      // Loop over the splines included in spp
      for ( m = 0; m < k; m++) {
	gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	for ( n = 0; n < k; n++) {
	  psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	  phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	}
      }
    }
  }
  
  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

  gsl_integration_glfixed_table_free (itable);
  gsl_integration_glfixed_table_free (jtable);
}

typedef struct {
  gdouble u;
  SPPanel * spp;
  gint se;
  gdouble xip, etap;
  gdouble du, dv;
  Point p;
  gint m, n;
  gdouble * psi, * phi;
  gdouble up, vp;
  gint k;
  gsl_matrix * Bu, * Bv;
  gsl_matrix * Bux;
} Qags_struct;

double dipole_v (double v, void * params) {
  Qags_struct * qs = (Qags_struct *) params;
  gdouble u = qs->u;

  return u*v;
}

double dipole_u (double u, void * params) {
  Qags_struct * qs = (Qags_struct *) params;

  qs->u = u;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
  double result, error;
     
  gsl_function F;
  F.function = &dipole_v;
  F.params = params;

  gsl_integration_qags (&F, -1, 1, 0, 1e-7, 1000,
			w, &result, &error);
     
  gsl_integration_workspace_free (w);

  return result;
}

double dipole ()
{
  Qags_struct * qs = g_malloc (sizeof(Qags_struct));
  

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &dipole_u;
  F.params = qs;
  
  double result, error;

  gsl_integration_qags (&F, -1., 1., 0, 1e-6, 1000,
			w, &result, &error);

  gsl_integration_workspace_free (w);

  return result;
}

void lachat_watson_self_influence_coefficients_qags (SPPanel * spp,
						     gdouble up, gdouble vp,
						     Point p,
						     gdouble * psi,
						     gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint ng = sp->ninner, i, j, m, n, se;
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  gdouble ui[ng], wi[ng];
  gdouble du = (spp->u1-spp->u0)/2.;
  gdouble dv = (spp->v1-spp->v0)/2.;
  gdouble xip = -1 + (up-spp->u0)/du;
  gdouble etap = -1 + (vp-spp->v0)/dv;
  gint k = spp->k;

  for ( i = 0; i < ng; i++) {
    gsl_integration_glfixed_point (-1., 1., i, &ui[i], &wi[i], itable);
  }
  
  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bux = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;

  // Loop over the triangles
  for ( se = 1; se <= 4; se ++) {

    // Loop over the Gauss-Point
    for ( i = 0; i < ng; i++) {
      for ( j = 0; j < ng; j++) {
	gdouble uu = xi (ui[i], ui[j], xip, se);
	gdouble vv = eta (ui[i], ui[j], etap, se);
	gdouble u = spp->u0 + (uu+1.)*du;
	gdouble v = spp->v0 + (vv+1.)*dv;
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);
	
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bux, &istart, &iend, sp->wx_u, sp->wxd_u);

	Point pg;
	Vector xu, xv;
	gint  a, b;

	xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
	pg.x = pg.y = pg.z = 0.;

	for ( a = 0; a < k; a++) {
	  gdouble cu = gsl_matrix_get (Bux, a, 0);
	  gdouble cdu = gsl_matrix_get (Bux, a, 1);
	  gint ii = (istart+a);
	  for ( b = 0; b < k; b++) {
	    gdouble cv = gsl_matrix_get (Bv, b, 0);
	    gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	    gdouble cvdu = cv*cdu;
	    gdouble cuv = cu*cv;
	    gint jj = (jstart+b);
	    gdouble v0 = coeff (sp, ii, jj, 0);
	    gdouble v1 = coeff (sp, ii, jj, 1);
	    gdouble v2 = coeff (sp, ii, jj, 2);

	    xu.x += v0*cvdu;
	    xu.y += v1*cvdu;
	    xu.z += v2*cvdu;
	    xv.x += v0*cudv;
	    xv.y += v1*cudv;
	    xv.z += v2*cudv;
	    pg.x += v0*cuv;
	    pg.y += v1*cuv;
	    pg.z += v2*cuv;
	  }
	}
	
	// Normal at Gauss-Legendre point
	Vector N = vector_vector_product (&xu, &xv);

	// Physical distance to Gauss-Legendre Point
	Vector R;
	R.x = pg.x - p.x;
	R.y = pg.y - p.y;
	R.z = pg.z - p.z;
	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
	// Jacobian at Gauss-Legendre point = norm of N
	gdouble J = vector_norm (N);
	gdouble c1 = wi[i]*wi[j]*J_se (ui[i], ui[j], xip, etap, se)*vector_norm (N)/sqrt(r2);
	N.x /= J; N.y /= J; N.z /= J;
	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

	// Loop over the splines included in spp
	for ( m = 0; m < k; m++) {
	  gdouble c3 = c1*gsl_matrix_get (Bu, m, 0);
	  gdouble c4 = c2*gsl_matrix_get (Bu, m, 0);
	  for ( n = 0; n < k; n++) {
	    psi[m + n*k] += c3*gsl_matrix_get (Bv, n, 0);
	    phi[m + n*k] += c4*gsl_matrix_get (Bv, n, 0);
	  }
	}
      }
    }
  }

  gdouble J = du*dv;
  for ( m = 0; m < k*k; m++) {
    psi[m] *= J;
    phi[m] *= J;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  gsl_matrix_free (Bux);

  gsl_integration_glfixed_table_free (itable);
}

typedef struct
{
  SPPanel * spp;
  Point p;
} ClosePointData;

/* int */
/*      rosenbrock_f (const gsl_vector * x, void *params, */
/*                    gsl_vector * f) */
/*      { */
/*        double a = ((struct rparams *) params)->a; */
/*        double b = ((struct rparams *) params)->b; */
     
/*        const double x0 = gsl_vector_get (x, 0); */
/*        const double x1 = gsl_vector_get (x, 1); */
     
/*        const double y0 = a * (1 - x0); */
/*        const double y1 = b * (x1 - x0 * x0); */
     
/*        gsl_vector_set (f, 0, y0); */
/*        gsl_vector_set (f, 1, y1); */
     
/*        return GSL_SUCCESS; */
/*      } */

static gdouble distance_function (const gsl_vector * x, void *params)
{

  ClosePointData * cp = (ClosePointData *) params;
  gdouble u = gsl_vector_get (x, 0);
  gdouble v = gsl_vector_get (x, 1);

  Point ptmp = spline2d_eval_point (cp->spp->sp, u, v);

  return point_distance (ptmp, cp->p);
}

static void distance_function_fdf (const gsl_vector * x, void * params,
				   double * f, gsl_vector * g)
{

  ClosePointData * cp = (ClosePointData *) params;
  Spline2D * sp = cp->spp->sp;
  gdouble u = gsl_vector_get (x, 0);
  gdouble v = gsl_vector_get (x, 1);
  

  gint i, j;
  size_t ustart, uend, vstart, vend;
  gsl_matrix * Bu = gsl_matrix_alloc (sp->k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (sp->k, 2);

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

  gdouble vx = 0., vxdu = 0., vxdv = 0.;
  gdouble vy = 0., vydu = 0., vydv = 0.;
  gdouble vz = 0., vzdu = 0., vzdv = 0.;

  for ( i = 0; i < sp->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    gdouble cdu = gsl_matrix_get (Bu, i, 1);
    gint ii = ustart;
    for ( j = 0; j < sp->k; j++) {
      gint jj = (ustart+j);
      gdouble cv = gsl_matrix_get (Bv, j, 0);
      gdouble cuv = cu*cv;
      gdouble cudv = cu*gsl_matrix_get (Bv, j, 1);
      gdouble cvdu = cu*cv;

      gdouble cx = coeff (sp, ii, jj, 0);
      gdouble cy = coeff (sp, ii, jj, 1);
      gdouble cz = coeff (sp, ii, jj, 2);

      vx += cx*cuv;
      vxdu += cx*cvdu;
      vxdv += cx*cudv;

      vy += cy*cuv;
      vydu += cy*cvdu;
      vydv += cy*cudv;

      vz += cz*cuv;
      vzdu += cz*cvdu;
      vzdv += cz*cudv;
    }
    ustart++;
  }

  gdouble dx = (vx-cp->p.x);
  gdouble dy = (vy-cp->p.y);
  gdouble dz = (vz-cp->p.z);

  gsl_vector_set (g, 0, 0.5*(dx*vxdu+dy*vydu+dz*vzdu));
  gsl_vector_set (g, 1, 0.5*(dx*vxdv+dy*vydv+dz*vzdv));
  
  *f = dx*dx + dy*dy + dz*dz;
}

/* static void closest_point_on_panel (SPPanel * spp, Point p) */
/* { */
/*   ClosePointData cp; */
/*   cp.spp = spp; */
/*   cp.p = p; */
     
/*   const gsl_multimin_fdfminimizer_type *T = */
/*     gsl_multimin_fdfminimizer_vector_bfgs2; */

/*   //    gsl_multimin_fminimizer_nmsimplex2; */
/*   gsl_multimin_fdfminimizer *s = NULL; */
/*   gsl_vector *ss, *x; */
/*   gsl_multimin_function_fdf minex_func; */
     
/*   size_t iter = 0; */
/*   int status; */
/*   double size; */
     
/*   /\* Starting point *\/ */
/*   x = gsl_vector_alloc (2); */
/*   gsl_vector_set (x, 0, spp->ue); */
/*   gsl_vector_set (x, 1, spp->ve); */
     
/*   /\* Set initial step sizes *\/ */
/*   ss = gsl_vector_alloc (2); */
/*   gsl_vector_set (ss, 0, (spp->u1-spp->u0)/2.); */
/*   gsl_vector_set (ss, 1, (spp->v1-spp->v0)/2.); */
     
/*   /\* Initialize method and iterate *\/ */
/*   minex_func.n = 2; */
/*   minex_func.fdf = distance_function_fdf; */
/*   minex_func.params = &cp; */
     
/*   s = gsl_multimin_fdfminimizer_alloc (T, 2); */
/*   gsl_multimin_fdfminimizer_set (s, &minex_func, x, ss, 1e-4); */
     
/*   /\* do *\/ */
/*   /\*   { *\/ */
/*   /\*     iter++; *\/ */
/*   /\*     status = gsl_multimin_fminimizer_iterate(s); *\/ */
           
/*   /\*     if (status) *\/ */
/*   /\* 	break; *\/ */
     
/*   /\*     size = gsl_multimin_fminimizer_size (s); *\/ */
/*   /\*     status = gsl_multimin_test_size (size, 1e-2); *\/ */
     
/*   /\*     if (status == GSL_SUCCESS) *\/ */
/*   /\* 	{ *\/ */
/*   /\* 	  printf ("converged to minimum at\n"); *\/ */
/*   /\* 	} *\/ */
     
/*   /\*     printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", *\/ */
/*   /\* 	      iter, *\/ */
/*   /\* 	      gsl_vector_get (s->x, 0), *\/ */
/*   /\* 	      gsl_vector_get (s->x, 1), *\/ */
/*   /\* 	      s->fval, size); *\/ */
/*   /\*   } *\/ */
/*   /\* while (status == GSL_CONTINUE && iter < 100); *\/ */
       
/*   /\* gsl_vector_free(x); *\/ */
/*   /\* gsl_vector_free(ss); *\/ */
/*   /\* gsl_multimin_fminimizer_free (s); *\/ */
     
/*   // return status; */
/* } */
      

/******** END LACHAT-WATSON method for self-influence coefficients ********/

/**
 * Analytical solution for the potential flow past a translating sphere
 **/
static gdouble analytical_potential (Point P, gdouble r0)
{
  gdouble r = sqrt (P.x*P.x + P.y*P.y + P.z*P.z);

  return -pow(r0,3)*P.x/(2.*pow(r,3));
}

static gdouble analytical_derivative (Point P, gdouble r0)
{
  gdouble r = sqrt (P.x*P.x + P.y*P.y + P.z*P.z);

  return -pow(r0,3)*P.x/(pow(r,4));
}

static gdouble analytical_dz (Point P, gdouble r0)
{
  gdouble r = sqrt (P.x*P.x + P.y*P.y + P.z*P.z);

  return 3./2.*pow(r0,3)*P.x*P.z/(pow(r,5));
}

static gdouble analytical_dzdz (Point P, gdouble r0)
{
  gdouble r = sqrt (P.x*P.x + P.y*P.y + P.z*P.z);

  return 3./2.*pow(r0,3)*P.x/(pow(r,5)) - 15./2.*pow(r0,3)*P.x*P.z*P.z/(pow(r,7));
}

static gdouble potential_at_point (GSList * list, Point p)
{
  gint m, n, ii, jj;
  GSList * patches = list;
  gdouble potential = 0.;

  while (patches) {
    Spline2D * patch = patches->data;
    g_assert (patch != NULL);
    
    // Loop over the panels of the patch
    for ( ii = 0; ii < patch->M; ii++) {
      for ( jj = 0; jj < patch->N ; jj++) { 
	SPPanel * panel = g_ptr_array_index (patch->panels, ii + jj*patch->M);
	g_assert (panel != NULL);
	
	// Gauss inner-integration
	GaussPoints * gp = panel->outer;
	gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
	for ( m = 0; m < ng; m++) {
	  for ( n = 0; n < ng; n++) {
	    gdouble wJ = g_array_index (gp->wJij, gdouble, m + n*ng);
	    Point pg = g_array_index (gp->Pi, Point, m + n*ng);
	    Vector N = g_array_index (gp->Ni, Vector, m+n*ng);
	    Vector R;
	    R.x = pg.x - p.x;
	    R.y = pg.y - p.y;
	    R.z = pg.z - p.z;
	    gdouble r = vector_norm (R);

	    // NB: The coefficient is 4*pi for the inside of the fluid domain
	    // it would have been 2*pi on the boundary (see Lee&Newman, Computation of wave effects using the panel method)
	    potential += 1./(4.*M_PI)*( N.x/r
					+ spline2d_eval (patch, ui(m), vj(n), 3)
					*vector_scalar_product (&N,&R)/pow(r,3.))*wJ;
	  }
	} // End Gauss integration

      }
    }  // End loop over panels

    patches = patches->next;
  }

  return potential;
}

void print_potential_sphere (GSList * list)
{
  FILE * fp = fopen("potential_sphere.tmp","w");
  gint i, j;
  gdouble NP = 100;
  gdouble r0 = 1.;

  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      Point p;
      p.x = -2.5 + 5./NP*i;
      p.y = 0.;
      p.z = -2.5 + 5./NP*j;
      gdouble r = sqrt (p.x*p.x + p.y*p.y + p.z*p.z);
      if (  r >= r0 ) {
	gdouble potential = potential_at_point (list, p);
	fprintf(fp, "%f %f %f %f %f\n", p.x, p.z,
		potential, analytical_potential (p, r0), analytical_potential (p, r0)/potential);
      }
      else {
      	fprintf(fp, "%f %f %f %f %f \n", p.x, p.z, 0., 0., 1.);
      }
    }
    fprintf(fp, "\n");
  }

  fclose (fp);
}

void print_potential_spheroid (GSList * list)
{
  FILE * fp = fopen("potential_spheroid.tmp","w");
  gint i, j;
  gdouble NP = 100;
  gdouble r0 = 1.;

  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      Point p;
      p.x = -1 + 2./NP*i;
      p.y = 0.;
      p.z = -1.5 + 2./NP*j;
      gdouble r = sqrt (p.x*p.x + p.y*p.y + p.z*p.z);
 
	gdouble potential = potential_at_point (list, p);
	fprintf(fp, "%f %f %f %f %f\n", p.x, p.z,
		potential, analytical_potential (p, r0), analytical_potential (p, r0)/potential);

    }
    fprintf(fp, "\n");
  }

  fclose (fp);
}

void print_potential_sphere_wall (GSList * list)
{
  FILE * fp = fopen("potential_sphere_wall.tmp","w");
  gint i, j;
  gdouble NP = 100;
  gdouble r0 = 0.5;

  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      Point p;
      p.x = -2. + 4./NP*i;
      p.y = 0.;
      p.z = -0.2 + 2.3/NP*j;
      gdouble r = sqrt (p.x*p.x + p.y*p.y + (p.z-1)*(p.z-1));
      if (  r >= r0 ) {
	gdouble potential = potential_at_point (list, p);
	fprintf(fp, "%f %f %f\n", p.x, p.z, potential);
      }
      else {
	fprintf(fp, "%f %f %f \n", p.x, p.z, 0.);
      }
    }
    fprintf(fp, "\n");
  }

  fclose (fp);
}

/**
 * Outputs the potential on each face of the cube for the test case of the
 * translating unit cube.
 **/
void print_potential_square (GSList * list)
{
  FILE * fp = fopen("potential_square-1.tmp","w");
  gint i, j;
  gdouble NP = 60;
  GSList * plist = list;
  Spline2D * sp = plist->data;

  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      gdouble ui = i/(NP);
      gdouble vj = j/(NP);
      
      fprintf (fp, "%f %f %f %f\n", spline2d_eval (sp, ui, vj, 0),
	       spline2d_eval (sp, ui, vj, 1),
	       spline2d_eval (sp, ui, vj, 2),
	       spline2d_eval (sp, ui, vj, 3));

    }
    fprintf(fp, "\n");
  }
  fclose (fp);
  plist = plist->next;
  

  sp = plist->data;
  fp = fopen("potential_square-2.tmp","w");
  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      gdouble ui = i/(NP);
      gdouble vj = j/(NP);

      fprintf (fp, "%f %f %f %f\n", spline2d_eval (sp, ui, vj, 0),
	       spline2d_eval (sp, ui, vj, 1),
	       spline2d_eval (sp, ui, vj, 2),
	       spline2d_eval (sp, ui, vj, 3));

    }
    fprintf(fp, "\n");
  }
  fclose (fp);
  plist = plist->next;

  sp = plist->data;
  fp = fopen("potential_square-3.tmp","w");
  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      gdouble ui = i/(NP);
      gdouble vj = j/(NP);

      fprintf (fp, "%f %f %f %f\n", spline2d_eval (sp, ui, vj, 0),
	       spline2d_eval (sp, ui, vj, 1),
	       spline2d_eval (sp, ui, vj, 2),
	       spline2d_eval (sp, ui, vj, 3));

    }
    fprintf(fp, "\n");
  }
  fclose (fp);
  plist = plist->next;

  sp = plist->data;
  fp = fopen("potential_square-4.tmp","w");
  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      gdouble ui = i/(NP);
      gdouble vj = j/(NP);

      fprintf (fp, "%f %f %f %f\n", spline2d_eval (sp, ui, vj, 0),
	       spline2d_eval (sp, ui, vj, 1),
	       spline2d_eval (sp, ui, vj, 2),
	       spline2d_eval (sp, ui, vj, 3));

    }
    fprintf(fp, "\n");
  }
  fclose (fp);
  plist = plist->next;

  sp = plist->data;
  fp = fopen("potential_square-5.tmp","w");
  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      gdouble ui = i/(NP);
      gdouble vj = j/(NP);

      fprintf (fp, "%f %f %f %f\n", spline2d_eval (sp, ui, vj, 0),
	       spline2d_eval (sp, ui, vj, 1),
	       spline2d_eval (sp, ui, vj, 2),
	       spline2d_eval (sp, ui, vj, 3));

    }
    fprintf(fp, "\n");
  }
  fclose (fp);
  plist = plist->next;

  sp = plist->data;
  fp = fopen("potential_square-6.tmp","w");
  for ( i = 0; i <= NP; i++) {
    for (j = 0; j <= NP; j++) {
      gdouble ui = i/(NP);
      gdouble vj = j/(NP);

      fprintf (fp, "%f %f %f %f\n", spline2d_eval (sp, ui, vj, 0),
	       spline2d_eval (sp, ui, vj, 1),
	       spline2d_eval (sp, ui, vj, 2),
	       spline2d_eval (sp, ui, vj, 3));

    }
    fprintf(fp, "\n");
  }
  fclose (fp);
}

/**
 * Outputs the potential profile for the cylinder test case
 * as descirbed in (Maniar,1995) page 136.
 **/
void print_cylinder_profile (GSList * list)
{
  FILE * fp = fopen("cylinder_profile.tmp","w");
  gint i;
  gdouble NP = 100;
  GSList * plist = list;

   
  Spline2D * sp = plist->data;
  for ( i = 0; i <= NP; i++) {
    fprintf(fp, "%f %f \n",
	    -spline2d_eval (sp, 0., 0.5*((gdouble) i)/NP, 0),
	    spline2d_eval (sp, 0., 0.5*((gdouble) i)/NP, 3));
  }
  plist = plist->next;


  sp = plist->data;
  for ( i = 0; i <= NP; i++) {
    fprintf(fp, "%f %f\n",
	    2.+spline2d_eval (sp, 0., ((gdouble) i)/NP, 1),
	    spline2d_eval (sp, 0., ((gdouble) i)/NP, 3));
  }
}

/* static gdouble potential_on_surface (GSList * list, SPPanel * spp, Point p, gint ui, gint vi) */
/* { */
/*   gint m, n, ii, jj; */
/*   GSList * patches = list; */
/*   gdouble potential = 0.; */

/*   while (patches) { */
/*     Spline2D * patch = patches->data; */
/*     g_assert (patch != NULL); */
    
/*     // Loop over the panels of the patch */
/*     for ( ii = 0; ii < patch->M; ii++) { */
/*       for ( jj = 0; jj < patch->N ; jj++) { */
/*   	SPPanel * panel = g_ptr_array_index (patch->panels, ii + jj*patch->M); */
/* 	gint k = panel->k; */
/*   	g_assert (panel != NULL); */
	
/*   	if (panel == spp) { // Self-influence */

/*   	  InfluenceCoeffs * ic = NULL; */
/*   	  GaussPoints * gp = spp->outer; */
/*   	  gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule */
/*   	  gint a, b, c, d; */

/*   	  // Outer Gauss-Legendre point parametric coodinates */
/*   	  /\* gdouble um = g_array_index (gp->ui, gdouble, ui); *\/ */
/*   	  /\* gdouble vn = g_array_index (gp->vj, gdouble, vi); *\/ */
	  
/*   	  S2P * s2p = g_ptr_array_index (spp->s2p, ui + vi*ng); */
/*   	  ic = sppanel_self_influence_coeff (spp, ui, vi); */

/* 	  g_assert (s2p != NULL); */
/* 	  g_assert (ic != NULL); */
	  
/*   	  for ( c = 0; c < k; c++) { */
/*   	    for ( d = 0; d < k; d++) { */
/* 	      g_assert (s2p->x2X != NULL); */
/* 	      g_assert (ic->psi != NULL); */
/*   	      for ( a = 0; a < k; a++) { */
/*   		for ( b = 0; b < k; b++) { */
		  
/* 		  potential += gsl_matrix_get (s2p->X2x, c + d*k, a+ b*k)*coeff(spp->sp, s2p->istart+a, s2p->jstart+b, 3) */
/* 		    *gsl_matrix_get (ic->psi, d, c)/2.; */
/*   		} */
/*   	      } */
	      
/*   	    } */
/*   	  } */

/*   	  influencecoeffs_destroy (ic); */
/*   	} */
/*   	else {	// Gauss inner-integration */
/*   	  GaussPoints * gp = panel->inner; */
/*   	  gint ng = gp->ui->len; // Order of inner Gauss-Legendre rule */
/*   	  for ( m = 0; m < ng; m++) { */
/*   	    for ( n = 0; n < ng; n++) { */
/*   	      gdouble wJ = g_array_index (gp->wJij, gdouble, m + n*ng); */
/*   	      Point pg = g_array_index (gp->Pi, Point, m + n*ng); */
/*   	      Vector R; */
/*   	      R.x = pg.x - p.x; */
/*   	      R.y = pg.y - p.y; */
/*   	      R.z = pg.z - p.z; */
/*   	      gdouble r = vector_norm (R); */
	      
/*   	      potential += spline2d_eval (patch, ui(m), vj(n), 3)*wJ/r; */
/*   	    } */
/*   	  } */
/*   	} // End Gauss integration */

/*       } */
/*     }  // End loop over panels */

/*     patches = patches->next; */
/*   } */

/*   return potential; */
/* } */

static gdouble spp_added_mass (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;

  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  
  return N.x*spline2d_eval_gauss_point (spp->sp, gp, m, n, 3);
}

gdouble calculate_added_mass (GSList * list)
{
  g_test_timer_start ();
  GSList * plist = list;
  gdouble sum = 0.;

  // Go over the list of patches
  while (plist != NULL) {
    Spline2D * sp = plist->data;
    g_assert (sp != NULL);

    while (sp) {
      sum += spline2d_gauss_integration (sp, spp_added_mass, NULL);
      sp = sp->next;
    }

    plist = plist->next;
  }

  fprintf (stdout, "   Calculate_added_mass: %f \n", g_test_timer_elapsed());
  return sum;
}

void print_potential_gauss (GSList * list)
{
  gint m, n, ii, jj;
  gint istart = 0;
  GSList * plist = list;
  FILE * fp = fopen ("potential_gauss.tmp","w");

  fprintf (fp, "#1:Theta1 2:Theta2 3:Phi 4:Phi_ref 5:Phin 6:Phin_ref 7:Phi_dz 8:Phi_dz_ref 9:Phi_dzz  10:Phidzz_ref \n");

  // Go over the list of patches
  while (plist != NULL) {
    Spline2D * sp = plist->data;
    g_assert (sp != NULL);

    while (sp) {
      // Loop over the panels of the patch
      for ( ii = 0; ii < sp->M; ii++) {
	for ( jj = 0; jj < sp->N ; jj++) { 
	  SPPanel * spp = g_ptr_array_index (sp->panels, ii + jj*sp->M);
	  g_assert (spp != NULL);
	
	  // Gauss inner-integration
	  GaussPoints * gp = spp->outer;
	  gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
	  for ( m = 0; m < ng; m++) {
	    for ( n = 0; n < ng; n++) {
	      Point P = g_array_index (gp->Pi, Point, m + n*ng);
	      gdouble r = sqrt (P.x*P.x + P.y*P.y + P.z*P.z);

	      Vector grad = potential_gradient_on_surface_gauss_point (sp, gp, m, n, 14);

	      fprintf (fp, "%f %f %f %f %f %f %f %f %f %f\n", 
		       atan(P.y/P.x),
		       acos(P.z/r),
		       spline2d_eval_gauss_point (sp, gp, m, n, 3),
		       analytical_potential (P, 0.5),
		       spline2d_eval_gauss_point (sp, gp, m, n, 4),
		       analytical_derivative (P, 0.5),
		       spline2d_eval_gauss_point (sp, gp, m, n, 14),
		       analytical_dz  (P, 0.5),
		       grad.z,
		       analytical_dzdz  (P, 0.5));
	    }
	  } // End Gauss integration

	}
      }  // End loop over panels

      sp = sp->next;
    }

   // shift index for next patch
   // istart += sp->NU*sp->NV;
    plist = plist->next;
  }
  fclose (fp);
}

gdouble zero_potential (SPPanel * spp, gint m, gint n, gpointer data)
{
  return 0.;
}

Vector zero_velocity (Spline2D * sp, gdouble u, gdouble v, gpointer data)
{
  Vector v0;
  v0.x = v0.y = v0.z = 0.;
  return v0;
}

Vector uniform_velocity (Spline2D * sp, gdouble u, gdouble v, gpointer data)
{
  Vector * U = data;

  return *U;
}

gdouble uniform_normal_velocity (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  Vector * U = data;
  //fprintf (stdout, "1:%f \n", vector_scalar_product (U, &N));
  return vector_scalar_product (U, &N);
}

Vector m_terms_dirichlet_bc (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  Vector * U = data;
  Vector grad = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  grad.x = U->x-grad.x;
  grad.y = U->x-grad.y;
  grad.z = U->x-grad.z;

  return grad;
}

gdouble m1_term_dirichlet_bc (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  Vector * U = data;
  Vector grad = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  return U->x-grad.x;
}

gdouble m2_term_dirichlet_bc (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  Vector * U = data;
  Vector grad = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  return U->y-grad.y;
}

gdouble m3_term_dirichlet_bc (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  Vector * U = data;
  Vector grad = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  return U->z-grad.z;
}

gdouble mode_forcing_1 (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);

  return N.x;
}

gdouble mode_forcing_2 (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  
  return N.y;
}

gdouble mode_forcing_3 (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);

  return N.z;
}

gdouble mode_forcing_4 (SPPanel * spp, gint m, gint n, gpointer data)
{
  Point * xg = (Point *) data;
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  Point P = g_array_index (gp->Pi, Point, m + n*ng);

  return (P.y-xg->y)*N.z-(P.z-xg->z)*N.y;
}

gdouble mode_forcing_5 (SPPanel * spp, gint m, gint n, gpointer data)
{
  Point * xg = (Point *) data;
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  Point P = g_array_index (gp->Pi, Point, m + n*ng);

  return (P.z-xg->z)*N.x-(P.x-xg->x)*N.z;
}

gdouble mode_forcing_6 (SPPanel * spp, gint m, gint n, gpointer data)
{
  Point * xg = (Point *) data;
  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;
  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  Point P = g_array_index (gp->Pi, Point, m + n*ng);

  return (P.x-xg->x)*N.y-(P.y-xg->y)*N.x;
}

gdouble zero_normal_velocity (SPPanel * spp, gint m, gint n, gpointer data)
{
  return 0.;
}

Vector potential_gradient_on_surface (Spline2D * sp, gdouble u, gdouble v, gint var)
{
  Vector grad;
  Vector xu, xv;
  
  xu.x = spline2d_derivative_eval (sp, u, v, 1, 0, 0);
  xu.y = spline2d_derivative_eval (sp, u, v, 1, 0, 1);
  xu.z = spline2d_derivative_eval (sp, u, v, 1, 0, 2);
  xv.x = spline2d_derivative_eval (sp, u, v, 0, 1, 0);
  xv.y = spline2d_derivative_eval (sp, u, v, 0, 1, 1);
  xv.z = spline2d_derivative_eval (sp, u, v, 0, 1, 2);

  /* Derivatives of the variables on the surface with respect to parametric variables */
  gdouble dvardu = spline2d_derivative_eval (sp, u, v, 1, 0, var);
  gdouble dvardv = spline2d_derivative_eval (sp, u, v, 0, 1, var);

  /* Normal to the surface */
  Vector N = vector_normalise (vector_vector_product (&xu, &xv));

  /* Normal derivative from boundary conditions */
  gdouble Vn = spline2d_eval (sp, u, v, var + 1);

  gdouble det = xu.x*(xv.y*N.z-xv.z*N.y)
    - xu.y*(N.z*xv.x-xv.z*N.x)
    + xu.z*(xv.x*N.y-xv.y*N.x);

  det = 1./det;

  /* (xv.y*N.z-xv.z*N.y)*dvardu + (xu.z*N.y-xu.y*N.z)*dvardv */
  /*   + (xu.y*xv.z-xu.z*xv.y)*dvardw; */
  /* (xv.z*N.x-xv.x*N.z)*dvardu + (xu.x*N.z-xu.z*N.x)*dvardv */
  /*   + (xu.z*xv.x-xu.x*xv.z)*dvardw; */
  /* (xv.x*N.y-xv.y*N.x)*dvardu + (N.x*xu.y-xu.x*N.y)*dvardv */
  /*   + (xu.x*xv.y-xu.y*xv.x)*dvardw; */

  grad.x = ((xv.y*N.z-xv.z*N.y)*dvardu + (xu.z*N.y-xu.y*N.z)*dvardv)*det
  	    + Vn*N.x;
  grad.y = ((xv.z*N.x-xv.x*N.z)*dvardu + (xu.x*N.z-xu.z*N.x)*dvardv)*det
  	    + Vn*N.y;
  grad.z = ((xv.x*N.y-xv.y*N.x)*dvardu + (N.x*xu.y-xu.x*N.y)*dvardv)*det
  	    + Vn*N.z;

  return grad;
}

Vector potential_gradient_on_surface_gauss_point (Spline2D * sp, GaussPoints * gp, gint m, gint n, gint var)
{
  Vector grad;
  gdouble dvardu = 0., dvardv = 0., Vn = 0.;

  gint  i, j, ng = sp->nouter;
  gint gindex = m+n*ng;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);

  istart = gp->istart;
  jstart = gp->jstart;

  if (sp->periodic)
    istart -= (sp->k-1);

  for ( i = 0; i < sp->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    gdouble cdu = gsl_matrix_get (Bu, i, 1);
    gint ii = (istart+i);
    for ( j = 0; j < sp->k; j++) {
      gdouble cv = gsl_matrix_get (Bv, j, 0);
      gint jj = (jstart+j);

      SplineCoeffs * sc;
      if (sp->periodic)
	sc = g_ptr_array_index (sp->coeffs, ii+ jj*(sp->NU+sp->k-1));
      else
	sc = g_ptr_array_index (sp->coeffs, ii+ jj*(sp->NU));

      dvardu += sc->v[var]*cv*cdu;
      dvardv += sc->v[var]*cu*gsl_matrix_get (Bv, j, 1);
      Vn += sc->v[var+1]*cu*cv;
    }
  }

  /* Normal to the surface */
  Vector N = g_array_index (gp->Ni, Vector, gindex);

  /* Metric coefficients for general gradient */
  gdouble c1 = g_array_index (gp->c1, gdouble, gindex);
  gdouble c2 = g_array_index (gp->c2, gdouble, gindex);
  gdouble c3 = g_array_index (gp->c3, gdouble, gindex);
  gdouble c4 = g_array_index (gp->c4, gdouble, gindex);
  gdouble c5 = g_array_index (gp->c5, gdouble, gindex);
  gdouble c6 = g_array_index (gp->c6, gdouble, gindex);

  grad.x = c1*dvardu + c2*dvardv + Vn*N.x;
  grad.y = c3*dvardu + c4*dvardv + Vn*N.y;
  grad.z = c5*dvardu + c6*dvardv + Vn*N.z;

  return grad;
}

Vector potential_gradient2d_on_surface_gauss_point (Spline2D * sp,
						    GaussPoints * gp,
						    gint m, gint n,
						    gint var)
{
  Vector grad;
  gdouble dvardu = 0., dvardv = 0., Vn = 0.;

  gint  i, j, ng = sp->nouter;
  gint gindex = m+n*ng;
  size_t istart, jstart;
  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);

  istart = gp->istart;
  jstart = gp->jstart;

  if (sp->periodic)
    istart -= (sp->k-1);

  for ( i = 0; i < sp->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, 0);
    gdouble cdu = gsl_matrix_get (Bu, i, 1);
    gint ii = (istart+i);
    for ( j = 0; j < sp->k; j++) {
      gdouble cv = gsl_matrix_get (Bv, j, 0);
      gint jj = (jstart+j);

      SplineCoeffs * sc;
      if (sp->periodic)
	sc = g_ptr_array_index (sp->coeffs, ii+ jj*(sp->NU+sp->k-1));
      else
	sc = g_ptr_array_index (sp->coeffs, ii+ jj*(sp->NU));

      dvardu += sc->v[var]*cv*cdu;
      dvardv += sc->v[var]*cu*gsl_matrix_get (Bv, j, 1);
    }
  }

  /* Metric coefficients for general gradient */
  gdouble c1 = g_array_index (gp->c1, gdouble, gindex);
  gdouble c2 = g_array_index (gp->c2, gdouble, gindex);
  gdouble c3 = g_array_index (gp->c3, gdouble, gindex);
  gdouble c4 = g_array_index (gp->c4, gdouble, gindex);

  grad.x = c1*dvardu + c2*dvardv;
  grad.y = c3*dvardu + c4*dvardv;

  return grad;
}

Vector gradient_at_point (GSList * list, Point p, gint var)
{
  gint m, n, ii, jj;
  GSList * patches = list;
  Vector grad;
  grad.x = grad.y = grad.z = 0.;

  while (patches) {
    Spline2D * patch = patches->data;
    g_assert (patch != NULL);
    
    // Loop over the panels of the patch
    for ( ii = 0; ii < patch->M; ii++) {
      for ( jj = 0; jj < patch->N ; jj++) { 
	SPPanel * panel = g_ptr_array_index (patch->panels, ii + jj*patch->M);
	g_assert (panel != NULL);
	
	// Gauss inner-integration
	GaussPoints * gp = panel->outer;
	gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
	for ( m = 0; m < ng; m++) {
	  for ( n = 0; n < ng; n++) {
	    gdouble wJ = g_array_index (gp->wJij, gdouble, m + n*ng);
	    Point pg = g_array_index (gp->Pi, Point, m + n*ng);
	    Vector N = g_array_index (gp->Ni, Vector, m+n*ng);
	    Vector R;
	    R.x = pg.x - p.x;
	    R.y = pg.y - p.y;
	    R.z = pg.z - p.z;
	    gdouble r = vector_norm (R);

	    // NB: The coefficient is 4*pi for the inside of the fluid domain
	    // it would have been 2*pi on the boundary (see Lee&Newman, Computation of wave effects using the panel method)
	    // potential += 1./(4.*M_PI)*( N.x/r + spline2d_eval (patch, ui(m), vj(n), 3)*vector_scalar_product (&N,&R)/pow(r,3.))*wJ;

	    

	    /* grad.x = 1./(2.*M_PI)*(spline2d_eval (patch, ui(m), vj(n), var)*(N.x*r*r-R.x*R.x)/pow(r,6.) + N.x*R.x/pow(r,3)); */
	    /* grad.y = 1./(2.*M_PI)*(spline2d_eval (patch, ui(m), vj(n), var)*(N.y*r*r-R.y*R.y)/pow(r,6.) + N.y*R.y/pow(r,3)); */
	    /* grad.z = 1./(2.*M_PI)*(spline2d_eval (patch, ui(m), vj(n), var)*(N.z*r*r-R.z*R.z)/pow(r,6.) + N.z*R.z/pow(r,3)); */


	    grad.x += 1./(2.*M_PI)*(spline2d_eval (patch, ui(m), vj(n), var)*(N.x + 3.*R.x*vector_scalar_product (&N,&R)*r*r)/pow(r,3.) + N.x*R.x/pow(r,3))*wJ;
	    grad.y += 1./(2.*M_PI)*(spline2d_eval (patch, ui(m), vj(n), var)*(N.y + 3.*R.y*vector_scalar_product (&N,&R)*r*r)/pow(r,3.) + N.y*R.y/pow(r,3))*wJ;
	    grad.z += 1./(2.*M_PI)*(spline2d_eval (patch, ui(m), vj(n), var)*(N.z + 3.*R.z*vector_scalar_product (&N,&R)*r*r)/pow(r,3.) + N.z*R.z/pow(r,3))*wJ;


	  }
	} // End Gauss integration

      }
    }  // End loop over panels

    patches = patches->next;
  }

  return grad;
}

void print_velocity_on_surface (GSList * list, NeumannFunc Vn)
{
  gint ii, jj;
  GSList * patches = list;
  FILE * fp = fopen ("velocity.tmp","w");

  while (patches) {
    Spline2D * patch = patches->data;
    g_assert (patch != NULL);
    
    // Loop over the panels of the patch
    for ( ii = 0; ii < patch->M; ii++) {
      for ( jj = 0; jj < patch->N ; jj++) {
  	SPPanel * panel = g_ptr_array_index (patch->panels, ii + jj*patch->M);
	Point p = spline2d_eval_point (patch, panel->ue, panel->ve);

	/* Vector vel = gradient_on_surface (patch, panel->ue, panel->ve, */
	/* 				  Vn, 3); */

	Vector vel = gradient_at_point (list, p, 3);
	
	fprintf (fp, "%f %f %f %f %f %f\n", p.x, p.y, p.z, panel->ue, panel->ve, vector_norm (vel));
	fprintf (fp, "%f %f %f %f %f\n\n\n", p.x+ 0.25*vel.x, p.y+0.25*vel.y, p.z+0.25*vel.z, panel->ue, panel->ve);

      }
    }  // End loop over panels

    patches = patches->next;
  }


  fclose (fp);
}

/***** Optimized QAG integration *****/

static inline
void initialise (gsl_integration_workspace ** workspace, double a, double b, gint size)
{
  gint i;
  
  for ( i = 0; i < size; i++) {
    workspace[i]->size = 0;
    workspace[i]->nrmax = 0;
    workspace[i]->i = 0;
    workspace[i]->alist[0] = a;
    workspace[i]->blist[0] = b;
    workspace[i]->rlist[0] = 0.0;
    workspace[i]->elist[0] = 0.0;
    workspace[i]->order[0] = 0;
    workspace[i]->level[0] = 0;

    workspace[i]->maximum_level = 0;
  }
}

static inline
void set_initial_result (gsl_integration_workspace ** workspace,
                         double * result, double * error, gint size)
{
  gint i;

  for ( i = 0; i < size; i++) {
    workspace[i]->size = 1;
    workspace[i]->rlist[0] = result[i];
  }
  workspace[size-1]->elist[0] = *error;
}

static double rescale_error (double err, const double result_abs,
			     const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
        double scale = pow((200 * err / result_asc), 1.5) ;
        
        if (scale < 1)
          {
            err = result_asc * scale ;
          }
        else 
          {
            err = result_asc ;
          }
      }
  if (result_abs > GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))
    {
      double min_err = 50 * GSL_DBL_EPSILON * result_abs ;

      if (min_err > err) 
        {
          err = min_err ;
        }
    }
  
  return err ;
}

static inline void
retrieve (gsl_integration_workspace ** workspace, 
          double * a, double * b, double * r, double * e, gint size)
{
  gint j;

  for ( j = 0; j < size; j++) {
    const size_t i = workspace[j]->i;
    double * rlist = workspace[j]->rlist;
    r[j] = rlist[i];
  }

  size_t i = workspace[size-1]->i;
  double * elist = workspace[size-1]->elist;
  double * alist = workspace[size-1]->alist;
  double * blist = workspace[size-1]->blist;
  *e = elist[i] ;
  *a = alist[i] ;
  *b = blist[i] ;
}

static inline int
subinterval_too_small (double a1, double a2, double b2)
{
  const double e = GSL_DBL_EPSILON;
  const double u = GSL_DBL_MIN;

  double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);

  int status = fabs (a1) <= tmp && fabs (b2) <= tmp;

  return status;
}

static inline
void qpsrt (gsl_integration_workspace * workspace)
{
  const size_t last = workspace->size - 1;
  const size_t limit = workspace->limit;

  double * elist = workspace->elist;
  size_t * order = workspace->order;

  double errmax ;
  double errmin ;
  int i, k, top;

  size_t i_nrmax = workspace->nrmax;
  size_t i_maxerr = order[i_nrmax] ;
  
  /* Check whether the list contains more than two error estimates */

  if (last < 2) 
    {
      order[0] = 0 ;
      order[1] = 1 ;
      workspace->i = i_maxerr ;
      return ;
    }

  errmax = elist[i_maxerr] ;

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
    {
      order[i_nrmax] = order[i_nrmax - 1] ;
      i_nrmax-- ;
    } 

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */
  
  if(last < (limit/2 + 2)) 
    {
      top = last ;
    }
  else
    {
      top = limit - last + 1;
    }
  
  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). */
  
  i = i_nrmax + 1 ;
  
  /* The order of the tests in the following line is important to
     prevent a segmentation fault */

  while (i < top && errmax < elist[order[i]])
    {
      order[i-1] = order[i] ;
      i++ ;
    }
  
  order[i-1] = i_maxerr ;
  
  /* Insert errmin by traversing the list bottom-up */
  
  errmin = elist[last] ;
  
  k = top - 1 ;
  
  while (k > i - 2 && errmin >= elist[order[k]])
    {
      order[k+1] = order[k] ;
      k-- ;
    }
  
  order[k+1] = last ;

  /* Set i_max and e_max */

  i_maxerr = order[i_nrmax] ;
  
  workspace->i = i_maxerr ;
  workspace->nrmax = i_nrmax ;
}

static inline
void update (gsl_integration_workspace ** workspace,
             double a1, double b1, double * area1, double * error1,
             double a2, double b2, double * area2, double * error2,
	     gint size)
{
  gint i;

  for ( i = 0; i < size; i++ ) {

    double * alist = workspace[i]->alist ;
    double * blist = workspace[i]->blist ;
    double * rlist = workspace[i]->rlist ;
    double * elist = workspace[i]->elist ;
    size_t * level = workspace[i]->level ;

    const size_t i_max = workspace[i]->i ;
    const size_t i_new = workspace[i]->size ;

    const size_t new_level = workspace[i]->level[i_max] + 1;

    /* append the newly-created intervals to the list */
  
    if (error2 > error1)
      {
	alist[i_max] = a2;      /* blist[maxerr] is already == b2 */
	rlist[i_max] = area2[i];
	elist[i_max] = *error2;
	level[i_max] = new_level;
      
	alist[i_new] = a1;
	blist[i_new] = b1;
	rlist[i_new] = area1[i];
	elist[i_new] = *error1;
	level[i_new] = new_level;
      }
    else
      {
	blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
	rlist[i_max] = area1[i];
	elist[i_max] = *error1;
	level[i_max] = new_level;
      
	alist[i_new] = a2;
	blist[i_new] = b2;
	rlist[i_new] = area2[i];
	elist[i_new] = *error2;
	level[i_new] = new_level;
      }
  
    workspace[i]->size++;

    if (new_level > workspace[i]->maximum_level)
      {
	workspace[i]->maximum_level = new_level;
      }

    qpsrt (workspace[i]) ;
  }
}

static inline double
sum_results (gsl_integration_workspace * workspace)
{
  const double * const rlist = workspace->rlist ;
  const size_t n = workspace->size;

  size_t k;
  double result_sum = 0;

  for (k = 0; k < n; k++)
    {
      result_sum += rlist[k];
    }
  
  return result_sum;
}

/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */

static const double xgk15[8] =    /* abscissae of the 15-point kronrod rule */
{
  0.991455371120812639206854697526329,
  0.949107912342758524526189684047851,
  0.864864423359769072789712788640926,
  0.741531185599394439863864773280788,
  0.586087235467691130294144838258730,
  0.405845151377397166906606412076961,
  0.207784955007898467600689403773245,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */

static const double wg15[4] =     /* weights of the 7-point gauss rule */
{
  0.129484966168869693270611432679082,
  0.279705391489276667901467771423780,
  0.381830050505118944950369775488975,
  0.417959183673469387755102040816327
};

static const double wgk15[8] =    /* weights of the 15-point kronrod rule */
{
  0.022935322010529224963732008058970,
  0.063092092629978553290700663189204,
  0.104790010322250183839876322541518,
  0.140653259715525918745189590510238,
  0.169004726639267902826583426598550,
  0.190350578064785409913256402421014,
  0.204432940075298892414161999234649,
  0.209482141084727828012999174891714
};


static const double xgk7[4] =    /* abscissae of the 7-point kronrod rule */
{
  0.96049126870802028342350709262908,
  0.7745966692414833770358530799564799,
  0.4342437493468025580020715028446278,
  0.0000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */

static const double wg7[2] =     /* weights of the 7-point gauss rule */
{
  0.55555555555555555555555555555556,
  0.8888888888888888888888888888888889
};

static const double wgk7[4] =    /* weights of the 7-point kronrod rule */
{
  0.104656226026467265193823857192073,
  0.26848808986833344072856928066671,
  0.401397414775962222905051818618432,
  0.450916538658474142345110087045571
};


static const double xgk5[3] =    /* abscissae of the 5-point kronrod rule */
{
  0.9258200997725514615665667765839995,
  0.5773502691896257645091487805019575,
  0.0000000000000000000000000000000000
};

static const double wg5[1] =     /* weights of the 5-point gauss rule */
{
  1.
};

static const double wgk5[3] =    /* weights of the 5-point kronrod rule */
{
  0.197979797979797979797979797979798,
  0.4909090909090909090909090909090909,
  0.6222222222222222222222222222222222
};

static const double xgk9[5] =    /* abscissae of the 5-point kronrod rule */
{
  0.9765602507375731115345053593699196,
  0.8611363115940525752239464888928095,
  0.6402862174963099824046890231574920,
  0.3399810435848562648026657591032447,
  0.0000000000000000000000000000000000
};

static const double wg9[2] =     /* weights of the 5-point gauss rule */
{
  0.34785484513745385737306394922199941,
  0.65214515486254614262693605077800059
};

static const double wgk9[5] =    /* weights of the 5-point kronrod rule */
{
  0.06297737366547301476549248855281868,
  0.17005360533572272680273885329620659,
  0.26679834045228444803277062841785567,
  0.32694918960145162955845946561731919,
  0.34644298189013636168107712823159978
};

static const double xgk11[6] =    /* abscissae of the 5-point kronrod rule */
{
  0.9840853600948424644961729346361395,
  0.9061798459386639927976268782993930,
  0.7541667265708492204408171669461159,
  0.5384693101056830910363144207002088,
  0.2796304131617831934134665227489774,
  0.0000000000000000000000000000000000
};

static const double wg11[3] =     /* weights of the 5-point gauss rule */
{
  0.23692688505618908751426404071991736,
  0.47862867049936646804129151483563819,
  0.56888888888888888888888888888888889
};

static const double wgk11[6] =    /* weights of the 5-point kronrod rule */
{
  0.04258203675108183286450945084767009,
  0.11523331662247339402462684588057354,
  0.18680079655649265746780002687848597,
  0.24104033922864758669994261122326211,
  0.27284980191255892234099326448445552,
  0.28298741785749121320425560137110554
};

static const double xgk13[7] =    /* abscissae of the 5-point kronrod rule */
{
  0.9887032026126788575046459517121851,
  0.9324695142031520278123015544939946,
  0.8213733408650279400456498342439503,
  0.6612093864662645136613995950199053,
  0.4631182124753046121567583640191766,
  0.2386191860831969086305017216807119,
  0.0000000000000000000000000000000000
};

static const double wg13[3] =     /* weights of the 5-point gauss rule */
{
  0.1713244923791703450402961421727329,
  0.3607615730481386075698335138377161,
  0.4679139345726910473898703439895510
};

static const double wgk13[7] =    /* weights of the 5-point kronrod rule */
{
  0.03039615411981976885196454467602788,
  0.08369444044690662613284560348241111,
  0.13732060463444692308714987253378181,
  0.18107199432313761518699209331551194,
  0.21320965227196227916289416351688930,
  0.23377086411699440662283572598899837,
  0.24107258017346476191063599297275917
};

struct dipole_function_struct
{
  void (* function) (double x, void * params, gdouble * res);
  void * params;
};

typedef struct dipole_function_struct dipole_function ;

#define DIPOLE_FN_EVAL(F,x,res) (*((F)->function))(x,(F)->params,res)

void gsl_integration_qk_new (const int n, 
			     const double xgk[], const double wg[], const double wgk[],
			     double fv1[], double fv2[],
			     const dipole_function * f,
			     double a, double b,
			     double *result, double *abserr,
			     double *resabs, double *resasc, gint size)
{
  gint i;

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);

  // That or the table is passed as argument ??
  double f_center [size]; 
  DIPOLE_FN_EVAL (f, center, f_center);

  gdouble result_gauss = 0.;
  gdouble result_kronrod [size];
  gdouble result_abs;
  gdouble result_asc = 0.;
  gdouble mean;
  gdouble err;

  for ( i = 0; i < size; i++)
    result_kronrod[i] = f_center[i] * wgk[n - 1];

  result_abs = fabs (result_kronrod[size-1]);

  int j;

  if (n % 2 == 0)
    result_gauss = f_center[size-1] * wg[n / 2 - 1];

  for (j = 0; j < (n - 1) / 2; j++) {
    const int jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
    const double abscissa = half_length * xgk[jtw];
    gdouble fval1 [size];
    gdouble fval2 [size];
    DIPOLE_FN_EVAL (f, center - abscissa, fval1);
    DIPOLE_FN_EVAL (f, center + abscissa, fval2);

    for ( i = 0; i < size; i++)
      result_kronrod[i] += wgk[jtw] * (fval1[i] + fval2[i]);

    fv1[jtw] = fval1[size-1];
    fv2[jtw] = fval2[size-1];
    result_abs += wgk[jtw] * (fabs (fval1[size-1]) + fabs (fval2[size-1]));
    result_gauss += wg[j] * (fval1[size-1] + fval2[size-1]);
  }

  for (j = 0; j < n / 2; j++) {
    int jtwm1 = j * 2;
    const double abscissa = half_length * xgk[jtwm1];
    
    gdouble fval1 [size];
    gdouble fval2 [size];
    DIPOLE_FN_EVAL (f, center - abscissa, fval1);
    DIPOLE_FN_EVAL (f, center + abscissa, fval2);

    for ( i = 0; i < size; i++)	
      result_kronrod[i] += wgk[jtwm1] * (fval1[i] + fval2[i]);
    
    fv1[jtwm1] = fval1[size-1];
    fv2[jtwm1] = fval2[size-1];
    result_abs += wgk[jtwm1] * (fabs (fval1[size-1]) + fabs (fval2[size-1]));
  }

  mean = result_kronrod[size-1] * 0.5;

  result_asc = wgk[n - 1] * fabs (f_center[size-1] - mean);

  for (j = 0; j < n - 1; j++)
    result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));

  /* scale by the width of the integration region */
  for ( i = 0; i < size; i++)
    result[i] = result_kronrod[i]*half_length;

  *resabs = result_abs*abs_half_length;
  *resasc = result_asc*abs_half_length;

  err = (result_kronrod[size-1] - result_gauss)* half_length;

  *abserr = rescale_error (err, result_abs, result_asc);
}

void gsl_integration_qk15_new (const dipole_function * f,
			       double a, double b,
			       double *result, double *abserr,
			       double *resabs, double *resasc, gint k)
{
  double fv1[8], fv2[8];
  gsl_integration_qk_new (8, xgk15, wg15, wgk15, fv1, fv2, f, a, b,
			  result, abserr, resabs, resasc, k);
}

void gsl_integration_qk13_new (const dipole_function * f,
			       double a, double b,
			       double *result, double *abserr,
			       double *resabs, double *resasc, gint k)
{
  double fv1[7], fv2[7];
  gsl_integration_qk_new (7, xgk13, wg13, wgk13, fv1, fv2, f, a, b,
			  result, abserr, resabs, resasc, k);
}

void gsl_integration_qk11_new (const dipole_function * f,
			       double a, double b,
			       double *result, double *abserr,
			       double *resabs, double *resasc, gint size)
{
  double fv1[6], fv2[6];
  gsl_integration_qk_new (6, xgk11, wg11, wgk11, fv1, fv2, f, a, b,
			  result, abserr, resabs, resasc, size);
}

void gsl_integration_qk9_new (const dipole_function * f,
			      double a, double b,
			      double *result, double *abserr,
			      double *resabs, double *resasc, gint k)
{
  double fv1[5], fv2[5];
  gsl_integration_qk_new (5, xgk9, wg9, wgk9, fv1, fv2, f, a, b,
			  result, abserr, resabs, resasc, k);
}

void gsl_integration_qk7_new (const dipole_function * f,
			      double a, double b,
			      double *result, double *abserr,
			      double *resabs, double *resasc, gint k)
{
  double fv1[4], fv2[4];
  gsl_integration_qk_new (4, xgk7, wg7, wgk7, fv1, fv2, f, a, b,
			  result, abserr, resabs, resasc, k);
}

void gsl_integration_qk5_new (const dipole_function * f,
			      double a, double b,
			      double *result, double *abserr,
			      double *resabs, double *resasc, gint k)
{
  double fv1[3], fv2[3];
  gsl_integration_qk_new (3, xgk5, wg5, wgk5, fv1, fv2, f, a, b,
			  result, abserr, resabs, resasc, k);
}

void dipole_qag_lw_v (double v, void * params, gdouble * res) {
  Qags_struct * qs = (Qags_struct *) params;

  Spline2D * sp = qs->spp->sp;
  gint k = sp->k;

  g_assert_not_reached ();
  // Needs to be updated to new splines structure

  gdouble uu = xi (qs->u, v, qs->xip, qs->se);
  gdouble vv = eta (qs->u, v, qs->etap, qs->se);
  gdouble uu2 = qs->spp->u0 + (uu+1.)*qs->du;
  gdouble vv2 = qs->spp->v0 + (vv+1.)*qs->dv;

  size_t istart, iend, jstart, jend;

  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,uu2)), 1, qs->Bu, &istart, &iend, sp->wx_u, sp->wxd_u);
  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,vv2)), 1, qs->Bv, &jstart, &jend, sp->w_v, sp->wd_v);

  if (sp->periodic)
    istart -= (k-1);
  
  Point pg;
  Vector xu, xv;
  gint a, b;

  xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
  pg.x = pg.y = pg.z = 0.;

  for ( a = 0; a < k; a++) {
    gdouble cu = gsl_matrix_get (qs->Bu, a, 0);
    gdouble cdu = gsl_matrix_get (qs->Bu, a, 1);
    gint ii = (istart+a);
    for ( b = 0; b < k; b++) {
      gdouble cv = gsl_matrix_get (qs->Bv, b, 0);
      gdouble cudv = cu*gsl_matrix_get (qs->Bv, b, 1);
      gdouble cvdu = cv*cdu;
      gdouble cuv = cu*cv;
      gint jj = (jstart+b);
      gdouble v0 = coeff (sp, ii, jj, 0);
      gdouble v1 = coeff (sp, ii, jj, 1);
      gdouble v2 = coeff (sp, ii, jj, 2);

      xu.x += v0*cvdu;
      xu.y += v1*cvdu;
      xu.z += v2*cvdu;
      xv.x += v0*cudv;
      xv.y += v1*cudv;
      xv.z += v2*cudv;
      pg.x += v0*cuv;
      pg.y += v1*cuv;
      pg.z += v2*cuv;
    }
  }
	
  // Normal at Gauss-Legendre point
  Vector N = vector_vector_product (&xu, &xv);

  // Physical distance to Gauss-Legendre Point
  Vector R;
  R.x = pg.x - qs->p.x;
  R.y = pg.y - qs->p.y;
  R.z = pg.z - qs->p.z;
  gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
  // Jacobian at Gauss-Legendre point = norm of N
  gdouble J = vector_norm (N);
  gdouble c1 = J_se (qs->u, v, qs->xip, qs->etap, qs->se)/sqrt(r2);

  gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

  c1 *= J;

  for ( a = 0; a < k; a++) {
    gdouble cc1 = c1*gsl_matrix_get (qs->Bu, a, 0);
    gdouble cc2 = c2*gsl_matrix_get (qs->Bu, a, 0);
    for ( b = 0; b < k; b++) {
      res[a+b*k] = cc1*gsl_matrix_get (qs->Bv, b, 0);
      res[a+b*k+k*k] = cc2*gsl_matrix_get (qs->Bv, b, 0);
    }
  }
}

static int qag_new (const dipole_function * f,
		    const double a, const double b,
		    const double epsabs, const double epsrel,
		    const size_t limit,
		    double *result, double *abserr,		    
		    gint k)
{
  gint size = 2*k*k;
  double area[size], errsum;
  double result0[size], abserr0, resabs0, resasc0;
  gsl_integration_workspace * workspaces[size];
  double tolerance;
  size_t iteration = 0;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;
  gint i;

  for ( i = 0; i < size; i++)
    workspaces[i] = gsl_integration_workspace_alloc (limit);

  double round_off;

  /* Initialize results */

  initialise (workspaces, a, b, size);

  for ( i = 0; i < size; i++) {
    result[i] = 0.;
  }
  *abserr = 0.;

  if (limit > workspaces[size-1]->limit)
    {
      GSL_ERROR ("iteration limit exceeds available workspace", GSL_EINVAL) ;
    }

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
      GSL_ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
                 GSL_EBADTOL);
    }

  /* perform the first integration */


  // gsl_integration_qk15_new (f, a, b, result0, &abserr0, &resabs0, &resasc0, size);
  /* gsl_integration_qk13_new (f, a, b, result0, &abserr0, &resabs0, &resasc0, size); */
  gsl_integration_qk11_new (f, a, b, result0, &abserr0, &resabs0, &resasc0, size);
  /* gsl_integration_qk9_new (f, a, b, result0, &abserr0, &resabs0, &resasc0, size); */
  /* gsl_integration_qk7_new (f, a, b, result0, &abserr0, &resabs0, &resasc0, size); */
  /* gsl_integration_qk5_new (f, a, b, result0, &abserr0, &resabs0, &resasc0, size); */

  set_initial_result (workspaces, result0, &abserr0, size);

  /* Test on accuracy */

  tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result0[size-1]));

  /* need IEEE rounding here to match original quadpack behavior */

  round_off = gsl_coerce_double (50 * GSL_DBL_EPSILON * resabs0);

  if (abserr0 <= round_off && abserr0 > tolerance)
    {
      
      for ( i = 0;i < size; i++) {
  	result[i] = result0[i];
      }
      *abserr = abserr0;

      GSL_ERROR ("cannot reach tolerance because of roundoff error "
                 "on first attempt", GSL_EROUND);
    }
  else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
    {

      for ( i = 0;i < size; i++)
  	result[i] = result0[i];
      *abserr = abserr0;

      for ( i = 0; i < size; i++)
	gsl_integration_workspace_free (workspaces[i]);
      return GSL_SUCCESS;
    }
  else if (limit == 1)
    {

      for ( i = 0;i < size; i++)
  	result[i] = result0[i];
      *abserr = abserr0; 

      GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
    }

  for ( i = 0;i < size; i++)
    area[i] = result0[i];
  errsum = abserr0;

  iteration = 1;

  do
    {
      double a1, b1, a2, b2;
      double a_i, b_i, r_i[size], e_i;
      double area1[size], area2[size], area12[size];
      double error1 = 0., error2 = 0., error12 = 0.;
      double resasc1, resasc2;
      double resabs1, resabs2;

      for ( i = 0; i < size; i++)
  	area1[i] = area2[i] = area12[i] = 0.;

      /* Bisect the subinterval with the largest error estimate */

      retrieve (workspaces, &a_i, &b_i, r_i, &e_i, size);

      a1 = a_i;
      b1 = 0.5 * (a_i + b_i);
      a2 = b1;
      b2 = b_i;

      /* gsl_integration_qk15_new (f, a1, b1, area1, &error1, &resabs1, &resasc1, size); */
      /* gsl_integration_qk15_new (f, a2, b2, area2, &error2, &resabs2, &resasc2, size); */

      /* gsl_integration_qk13_new (f, a1, b1, area1, &error1, &resabs1, &resasc1, size); */
      /* gsl_integration_qk13_new (f, a2, b2, area2, &error2, &resabs2, &resasc2, size); */

      gsl_integration_qk11_new (f, a1, b1, area1, &error1, &resabs1, &resasc1, size);
      gsl_integration_qk11_new (f, a2, b2, area2, &error2, &resabs2, &resasc2, size);

      /* gsl_integration_qk9_new (f, a1, b1, area1, &error1, &resabs1, &resasc1, size); */
      /* gsl_integration_qk9_new (f, a2, b2, area2, &error2, &resabs2, &resasc2, size); */

      /* gsl_integration_qk7_new (f, a1, b1, area1, &error1, &resabs1, &resasc1, size); */
      /* gsl_integration_qk7_new (f, a2, b2, area2, &error2, &resabs2, &resasc2, size); */

      /* gsl_integration_qk5_new (f, a1, b1, area1, &error1, &resabs1, &resasc1, size); */
      /* gsl_integration_qk5_new (f, a2, b2, area2, &error2, &resabs2, &resasc2, size);  */


      error12 = error1 + error2;
      errsum += (error12 - e_i);
      for ( i = 0; i < size; i++) {
  	area12[i] = area1[i] + area2[i];
  	area[i] += area12[i] - r_i[i];
      }

      if (resasc1 != error1 && resasc2 != error2) {
	double delta = r_i[size-1] - area12[size-1];

	if (fabs (delta) <= 1.0e-5 * fabs (area12[size-1]) && error12 >= 0.99 * e_i)
	  {
	    roundoff_type1++;
	  }
	if (iteration >= 10 && error12 > e_i) {
	  roundoff_type2++;
	}
      }

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area[size-1]));

      if (errsum > tolerance)
        {
          if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
            {
              error_type = 2;   /* round off error */
            }

          /* set error flag in the case of bad integrand behaviour at
             a point of the integration range */

          if (subinterval_too_small (a1, a2, b2))
            {
              error_type = 3;
            }
        }

      update (workspaces, a1, b1, area1, &error1, a2, b2, area2, &error2, size);

      retrieve (workspaces, &a_i, &b_i, r_i, &e_i, size);

      iteration++;

    }
  while (iteration < limit && !error_type && errsum > tolerance);

  for ( i = 0; i < size; i++) {
    result[i] = sum_results (workspaces[i]);
  }
  *abserr = errsum;

  /* if (errsum[size-1] <= tolerance) */
  /*   { */
  /*     return GSL_SUCCESS; */
  /*   } */
  /* else if (error_type == 2) */
  /*   { */
  /*     GSL_ERROR ("roundoff error prevents tolerance from being achieved", */
  /*                GSL_EROUND); */
  /*   } */
  /* else if (error_type == 3) */
  /*   { */
  /*     GSL_ERROR ("bad integrand behavior found in the integration interval", */
  /*                GSL_ESING); */
  /*   } */
  /* else if (iteration == limit) */
  /*   { */
  /*     GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER); */
  /*   } */
  /* else */
  /*   { */
  /*     GSL_ERROR ("could not integrate function", GSL_EFAILED); */
  /*   } */

  for ( i = 0; i < size; i++)
    gsl_integration_workspace_free (workspaces[i]);
}

void dipole_qag_lw_u (double u, void * params, gdouble * result)
{
  Qags_struct * qs = (Qags_struct *) params;
  qs->u = u;
  gdouble error;

  dipole_function F;
  F.function = &dipole_qag_lw_v;
  F.params = params;

  qag_new (&F, -1., 1., 0., 1e-6, 1000,
	   result, &error, qs->k);

}

void dipole_qag_lw (SPPanel * spp, gint se,
		    gdouble xip, gdouble etap,
		    gdouble du, gdouble dv, Point p,
		    gdouble * psi,
		    gdouble * phi,
		    gint k)
{
  Qags_struct * qs = g_malloc (sizeof(Qags_struct));
  qs->spp = spp;
  qs->se = se;
  qs->xip = xip;
  qs->etap = etap;
  qs->du = du;
  qs->dv = dv;
  qs->p = p;
  qs->k = k;
  qs->Bu = gsl_matrix_alloc (k, 2);
  qs->Bv = gsl_matrix_alloc (k, 2);

  dipole_function F;
  F.function = &dipole_qag_lw_u;
  F.params = qs;
  
  double result[2*k*k], error;

  qag_new (&F, -1., 1., 0., 1e-6, 1000,
	   result, &error, k);

  gsl_matrix_free (qs->Bu);
  gsl_matrix_free (qs->Bv);

  gint i;

  for ( i = 0; i < k*k; i++) {
    psi[i] += result[i];
    phi[i] += result[i+k*k];
  }
}

void lachat_watson_self_influence_coefficients_qag (SPPanel * spp,
						    gdouble up, gdouble vp,
						    Point p,
						    gdouble * psi,
						    gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gdouble du = (spp->u1-spp->u0)/2.;
  gdouble dv = (spp->v1-spp->v0)/2.;
  gdouble xip = -1. + (up-spp->u0)/du;
  gdouble etap = -1. + (vp-spp->v0)/dv;
  gint k = spp->k, m, n, se;

  for ( se = 1; se <= 4; se ++)
    dipole_qag_lw (spp, se, xip, etap,
		   du, dv, p, psi, phi, k);

  gdouble J = du*dv;
  for ( m = 0; m < k*k; m++) {
    psi[m] *= J;
    phi[m] *= J;
  }
}

/***** END - Optimised QAG integration *****/

static gdouble xi_ye (gdouble xi, gdouble eta, gdouble xip, gint se)
{
  switch (se) {
  case 1: return 0.5*((1.+xi)*eta + (1.-xi)*xip);
  case 2: return 0.5*((1.+xi) + (1-xi)*xip);
  case 3: return 0.5*(-(1.+xi)*eta + (1.-xi)*xip);
  case 4: return 0.5*(-(1.+xi) + (1.-xi)*xip);
  default: g_assert_not_reached ();
  }
}



void ye_transformation (SPPanel * spp)
{
  Point xs; /* Eval point */

  /* Find the closest point */
  Point xp;
  gdouble up, vp;

  /* Ends of the triangle */
  Point p1, p2;
  gdouble u1, v1, u2, v2;

  gdouble theta, rho;

  gdouble d = point_distance (xp,xs);
  gdouble c2 = pow(cos(theta),2.);
  gdouble s2 = pow(sin(theta),2.);

  gdouble eta1 = rho*c2;
  gdouble eta2 = rho*s2;

  Vector vp1, vp2;
  vp1.x = xp.x - p1.x;
  vp1.y = xp.y - p1.y;
  vp1.z = xp.z - p1.z;

  vp2.x = xp.x - p2.x;
  vp2.y = xp.y - p2.y;
  vp2.z = xp.z - p2.z;

  Vector tmp1, tmp2;
  tmp1.x = (xs.x-xp.x)/d;
  tmp1.y = (xs.y-xp.y)/d;
  tmp1.z = (xs.z-xp.z)/d;

  tmp2.x = c2*vp1.x+s2*vp2.x;
  tmp2.y = c2*vp1.y+s2*vp2.y;
  tmp2.z = c2*vp1.z+s2*vp2.z;

  gdouble f1 = vector_scalar_product (&tmp1, &tmp2);
  gdouble f2 = c2*c2*vector_scalar_product (&vp1, &vp1)
    + 2.*c2*s2*vector_scalar_product (&vp1, &vp2)
    + s2*s2*vector_scalar_product (&vp2, &vp2);
  

  gdouble u = eta1*u1 + eta2*u2 + (1.-eta1-eta2)*up;
  // or
  u = eta1*(u1-up) + eta2*(u2-up) + up;

  gdouble v = eta1*v1 + eta2*v2 + (1.-eta1-eta2)*vp;
  // or
  v = eta1*(v1-vp) + eta2*(v2-vp) + vp;
  
  
}

/* void qin_self_influence_coefficients (SPPanel * spp, */
/* 				      gdouble up, gdouble vp, */
/* 				      Point p, */
/* 				      gdouble * psi, */
/* 				      gdouble * phi) */
/* { */
/*   Spline2D * sp = spp->sp; */


/*   /\*  *\/ */
/*   gint i, j, m, n, se, ii, jj; */
/*   gint k = spp->k; */
/*   gint ng = 10.; */
/*   gdouble ui[k], wi[k]; */
/*   gdouble vj[k], wj[k]; */
/*   gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng); */
/*   gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng); */

/*   for ( i = 0; i < ng; i++) */
/*     gsl_integration_glfixed_point (0., 1., i, &ui[i], &wi[i], itable); */

/*   gsl_matrix * Bu = gsl_matrix_alloc (k, 2); */
/*   gsl_matrix * Bv = gsl_matrix_alloc (k, 2); */
/*   size_t istart, iend, jstart, jend; */

/*   /\* 10 Gauss-points over [0:1] *\/ */

/*   /\* NB: G is the Jacobian of the surface *\/ */

/*   /\* *\/ */
/*   gdouble uc, vc;  */
/*   Point pc = spline2d_eval_point (sp, uc, vc); */

/*   // For now  */
/*   pc = p; */
/*   uc = up; */
/*   vc = vp; */

/*   gdouble d[3]; */
/*   d[0] = pc.x - p.x; */
/*   d[1] = pc.y - p.y; */
/*   d[2] = pc.z - p.z; */

/*   gdouble r0_2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; */

/*   Vector xu, xv; */

/*   xu.x = spline2d_derivative_eval (sp, uc, vc, 1, 0, 0); */
/*   xu.y = spline2d_derivative_eval (sp, uc, vc, 1, 0, 1); */
/*   xu.z = spline2d_derivative_eval (sp, uc, vc, 1, 0, 2); */

/*   xv.x = spline2d_derivative_eval (sp, uc, vc, 0, 1, 0); */
/*   xv.y = spline2d_derivative_eval (sp, uc, vc, 0, 1, 1); */
/*   xv.z = spline2d_derivative_eval (sp, uc, vc, 0, 1, 2); */


/*   for (se = 1 ; se <= 4; se++) { */
/*     /\* Other ends of triangle *\/ */
/*     gdouble u1, v1; */
/*     gdouble u2, v2; */

/*     switch (se) { */
/*     case 1: u1 = spp->u0; v1 = spp->v0; u2 = spp->u1; v2 = spp->v0; break; */
/*     case 2: u1 = spp->u1; v1 = spp->v0; u2 = spp->u1; v2 = spp->v1;break; */
/*     case 3: u1 = spp->u1; v1 = spp->v1; u2 = spp->u0; v2 = spp->v1;break; */
/*     case 4: u1 = spp->u0; v1 = spp->v1; u2 = spp->u0; v2 = spp->v0;break; */
/*     default: g_assert_not_reached (); */
/*     } */

/*     gdouble kappa = fabs(u1*v2 + u2*vc + uc*v1 - u2*v1 - uc*v2 - u1*vc); */

/*     for ( m = 0; m < k; m++ ) { */
/*       gdouble beta = ui[m]; */



/*       gdouble A[3]; */
/*       A[0] = xu.x*(u1-uc + (u2-u1)*beta) + xv.x*(v1-vc + (v2-v1)*beta); */
/*       A[1] = xu.y*(u1-uc + (u2-u1)*beta) + xv.y*(v1-vc + (v2-v1)*beta); */
/*       A[2] = xu.z*(u1-uc + (u2-u1)*beta) + xv.z*(v1-vc + (v2-v1)*beta); */

/*       gdouble a = 0.; */
/*       for (i = 0; i < 3; i++) */
/* 	a += A[i]*A[i]; */

/*       gdouble b = 0.; */
/*       for (i = 0; i < 3; i++) */
/* 	b += 2.*d[i]*A[i]; */

/*       gdouble b_2a = b/(2.*a); */

/*       gdouble delta_2 = r0_2/a - b_2a*b_2a; */


/*       /\* gdouble alpha0 = log (0 + b_2a + g (0, beta)); *\/ */
/*       /\* gdouble alpha1 = log (1 + b_2a + g (1, beta)); *\/ */

/*       gdouble eta0 = log (0. + b_2a + sqrt (pow(0. + b_2a,2.) + delta_2) ); */
/*       gdouble eta1 = log (1. + b_2a + sqrt (pow(1. + b_2a,2.) + delta_2) ); */

/*       for ( i = 0; i < ng; i++) */
/* 	gsl_integration_glfixed_point (eta0, eta1, i, &vj[i], &wj[i], jtable); */

/*       for ( n = 0; n < ng; n++ ) { */

/* 	gdouble eta = vj[n]; */

/* 	gdouble lambda = 0.5*(exp(eta) - delta_2*exp(-eta)) - b_2a; */

/* 	gdouble g = sqrt (pow(eta + b_2a,2.) + delta_2); */

/* 	Vector N = spline2d_dimensional_normal (sp, eta, beta); */
/* 	// Physical distance to Gauss-Legendre Point */
/* 	Point pg = spline2d_eval_point (sp, eta, beta); */
/* 	Vector R; */
/* 	R.x = pg.x - p.x; */
/* 	R.y = pg.y - p.y; */
/* 	R.z = pg.z - p.z; */
/* 	gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z; */
	  
/* 	// Jacobian at Gauss-Legendre point = norm of N */
/* 	gdouble J = vector_norm (N); */
/* 	gdouble c1 = wi[i]*wj[j]/\* *J *\//sqrt(r2)*g*kappa*lambda; */


/* 	/\* N.x /= J; N.y /= J; N.z /= J; *\/ */

/* 	gdouble c2 = c1*vector_scalar_product (&N, &R)/r2; */

/* 	c1 *= J; */

/* 	// Loop over the splines included in spp */
/* 	for ( ii = 0; ii < k; ii++) { */
/* 	  gdouble c3 = c1*gsl_matrix_get (Bu, ii, 0); */
/* 	  gdouble c4 = c2*gsl_matrix_get (Bu, ii, 0); */
/* 	  for ( jj = 0; jj < k; jj++) { */
/* 	    psi[ii + jj*k] += c3*gsl_matrix_get (Bv, jj, 0); */
/* 	    phi[ii + jj*k] += c4*gsl_matrix_get (Bv, jj, 0); */
/* 	  } */
/* 	} */

/*       } */
/*     } */
/*   } */

/*   gsl_matrix_free (Bu); */
/*   gsl_matrix_free (Bv); */

/*   gsl_integration_glfixed_table_free (itable); */
/*   gsl_integration_glfixed_table_free (jtable); */


/* } */



/* void qin_self_influence_coefficients (SPPanel * spp, */
/* 				      gdouble up, gdouble vp, */
/* 				      Point p, */
/* 				      gdouble * psi, */
/* 				      gdouble * phi) */
/* { */
/*   Spline2D * sp = spp->sp; */
/*   gint i, j, m, n, se, ii, jj; */
/*   gint k = spp->k; */
/*   gint ng = 16.; */
/*   gdouble ui[ng], wi[ng]; */
/*   gdouble vj[ng], wj[ng]; */
/*   gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng); */
/*   gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng); */

/*   for ( i = 0; i < ng; i++) */
/*     gsl_integration_glfixed_point (0., 1., i, &ui[i], &wi[i], itable); */

/*   gsl_matrix * Bu = gsl_matrix_alloc (k, 2); */
/*   gsl_matrix * Bv = gsl_matrix_alloc (k, 2); */
/*   size_t istart, iend, jstart, jend; */

/*   /\* 10 Gauss-points over [0:1] *\/ */

/*   /\* NB: G is the Jacobian of the surface *\/ */

/*   /\* *\/ */
/*   gdouble uc, vc;  */
/*   //Point pc = spline2d_eval_point (sp, uc, vc); */

/*   // For now  */
/*   Point pc = p; */
/*   uc = up; */
/*   vc = vp; */

/*   gdouble d[3]; */
/*   d[0] = pc.x - p.x; */
/*   d[1] = pc.y - p.y; */
/*   d[2] = pc.z - p.z/\* +1e-6 *\/; */
  

/*   gdouble r0_2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; */

/*   Vector xu, xv; */

/*   xu.x = spline2d_derivative_eval (sp, uc, vc, 1, 0, 0); */
/*   xu.y = spline2d_derivative_eval (sp, uc, vc, 1, 0, 1); */
/*   xu.z = spline2d_derivative_eval (sp, uc, vc, 1, 0, 2); */

/*   xv.x = spline2d_derivative_eval (sp, uc, vc, 0, 1, 0); */
/*   xv.y = spline2d_derivative_eval (sp, uc, vc, 0, 1, 1); */
/*   xv.z = spline2d_derivative_eval (sp, uc, vc, 0, 1, 2); */


/*   for (se = 1 ; se <= 4; se++) { */
/*     /\* Other ends of triangle *\/ */
/*     gdouble u1, v1; */
/*     gdouble u2, v2; */

/*     switch (se) { */
/*     case 1: u1 = spp->u0; v1 = spp->v0; u2 = spp->u1; v2 = spp->v0; break; */
/*     case 2: u1 = spp->u1; v1 = spp->v0; u2 = spp->u1; v2 = spp->v1;break; */
/*     case 3: u1 = spp->u1; v1 = spp->v1; u2 = spp->u0; v2 = spp->v1;break; */
/*     case 4: u1 = spp->u0; v1 = spp->v1; u2 = spp->u0; v2 = spp->v0;break; */
/*     default: g_assert_not_reached (); */
/*     } */

/*     if ( (u1 != uc || u2 != uc) && (v1 != vc || v2 != vc ) ) { */

/*       /\* gdouble kappa = fabs(u1*v2 + u2*vc + uc*v1 - u2*v1 - uc*v2 - u1*vc); *\/ */
/*       gdouble kappa = ((u1-uc)*(v2-v1) - (u2-u1)*(v1-vc)); */

/*       for ( m = 0; m < ng; m++ ) { */
/* 	gdouble beta = ui[m]; */


/* 	gdouble A[3]; */
/* 	A[0] = xu.x*(u1-uc + (u2-u1)*beta) + xv.x*(v1-vc + (v2-v1)*beta); */
/* 	A[1] = xu.y*(u1-uc + (u2-u1)*beta) + xv.y*(v1-vc + (v2-v1)*beta); */
/* 	A[2] = xu.z*(u1-uc + (u2-u1)*beta) + xv.z*(v1-vc + (v2-v1)*beta); */

/* 	gdouble a = 0.; */
/* 	for (i = 0; i < 3; i++) */
/* 	  a += A[i]*A[i]; */
	
/* 	gdouble b = 0.; */
/* 	for (i = 0; i < 3; i++) */
/* 	  b += 2.*d[i]*A[i]; */

/* 	gdouble b_2a = b/(2.*a); */

/* 	gdouble delta_2 = r0_2/a - b_2a*b_2a; */

/* 	gdouble eta0 = /\* log (0. + b_2a + sqrt (pow(0. ,2) + delta_2 )) *\/0.; */
/* 	gdouble eta1 = /\* log (1. + b_2a + sqrt (pow(1. ,2) + delta_2 )) *\/1.; */

/* 	for ( i = 0; i < ng; i++) */
/* 	  gsl_integration_glfixed_point (eta0, eta1, i, &vj[i], &wj[i], itable); */

/* 	for ( n = 0; n < ng; n++ ) { */
/* 	  gdouble alpha = /\* ui[n] *\/vj[n]; */

/* 	  gdouble lambda = 0.5*(exp(alpha) - delta_2*exp(-alpha) ) - b_2a; */
/* 	  gdouble g = sqrt ( pow(alpha + b_2a, 2.) + delta_2); */

/* 	  gdouble u = uc + (u1-uc)*alpha + (u2-u1)*alpha*beta; */
/* 	  gdouble v = vc + (v1-vc)*alpha + (v2-v1)*alpha*beta; */

/* 	   /\* fprintf (stdout, "%f %f \n", u, v); *\/ */

/* 	  gsl_bspline_deriv_eval_nonzero (u, 1, Bu, &istart, &iend, sp->w_u, sp->wd_u); */
/* 	  gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v); */

/* 	  Vector N = spline2d_dimensional_normal (sp, u, v); */
/* 	  // Physical distance to Gauss-Legendre Point */
/* 	  Point pg = spline2d_eval_point (sp, u, v); */
/* 	  Vector R; */
/* 	  R.x = pg.x - p.x; */
/* 	  R.y = pg.y - p.y; */
/* 	  R.z = pg.z - p.z; */
/* 	  gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z; */
	  
/* 	  // Jacobian at Gauss-Legendre point = norm of N */
/* 	  gdouble J = vector_norm (N); */
/* 	  /\* gdouble c1 = wi[m]*wi[n]*J/sqrt(r2)*alpha*kappa; *\/ */
/* 	  gdouble c1 = wi[m]*wj[n]/\* *J *\//sqrt(r2)*g*lambda*kappa; */

/* 	  gdouble c2 = c1*vector_scalar_product (&N, &R)/r2; */

/* 	  c1 *= J; */

/* 	  // Loop over the splines included in spp */
/* 	  for ( ii = 0; ii < k; ii++) { */
/* 	    gdouble c3 = c1*gsl_matrix_get (Bu, ii, 0); */
/* 	    gdouble c4 = c2*gsl_matrix_get (Bu, ii, 0); */
/* 	    for ( jj = 0; jj < k; jj++) { */
/* 	      psi[ii + jj*k] += c3*gsl_matrix_get (Bv, jj, 0); */
/* 	      phi[ii + jj*k] += c4*gsl_matrix_get (Bv, jj, 0); */
/* 	    } */
/* 	  } */

/* 	} */
/*       } */
/*     } */
/*   } */

/*   gsl_matrix_free (Bu); */
/*   gsl_matrix_free (Bv); */

/*   gsl_integration_glfixed_table_free (itable); */
/*   gsl_integration_glfixed_table_free (jtable); */
/*   /\* g_assert_not_reached (); *\/ */


/*    /\*  for ( m = 0; m < k; m++) { *\/ */
/*   /\*   for ( n = 0; n < k; n++) { *\/ */
/*   /\*     fprintf (stdout, "%e %e \n", psi[m + n*k], phi[m + n*k]); *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */

/*   /\* g_assert_not_reached (); *\/ */

/* } */


void qin_self_influence_coefficients (SPPanel * spp,
				      gdouble up, gdouble vp,
				      Point p,
				      gdouble * psi,
				      gdouble * phi)
{
  Spline2D * sp = spp->sp;
  gint i, j, m, n, se, ii, jj;
  gint k = spp->k;
  gint ng = 35;
  gdouble ui[ng], wi[ng];
  /* gdouble vj[ng], wj[ng]; */
  gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng);
  /* gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng); */

  for ( i = 0; i < ng; i++)
    gsl_integration_glfixed_point (0., 1., i, &ui[i], &wi[i], itable);

  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;

  /* 10 Gauss-points over [0:1] */

  /* NB: G is the Jacobian of the surface */

  /* */
  gdouble uc, vc;
  //Point pc = spline2d_eval_point (sp, uc, vc);

  // For now
  Point pc = p;
  uc = up;
  vc = vp;

  /* gdouble d[3]; */
  /* d[0] = pc.x - p.x; */
  /* d[1] = pc.y - p.y; */
  /* d[2] = pc.z - p.z; */

  /* gdouble r0_2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; */

  /* Vector xu, xv; */

  /* xu.x = spline2d_derivative_eval (sp, uc, vc, 1, 0, 0); */
  /* xu.y = spline2d_derivative_eval (sp, uc, vc, 1, 0, 1); */
  /* xu.z = spline2d_derivative_eval (sp, uc, vc, 1, 0, 2); */

  /* xv.x = spline2d_derivative_eval (sp, uc, vc, 0, 1, 0); */
  /* xv.y = spline2d_derivative_eval (sp, uc, vc, 0, 1, 1); */
  /* xv.z = spline2d_derivative_eval (sp, uc, vc, 0, 1, 2); */


  for (se = 1 ; se <= 4; se++) {
    /* Other ends of triangle */
    gdouble u1, v1;
    gdouble u2, v2;

    switch (se) {
    case 1: u1 = spp->u0; v1 = spp->v0; u2 = spp->u1; v2 = spp->v0; break;
    case 2: u1 = spp->u1; v1 = spp->v0; u2 = spp->u1; v2 = spp->v1;break;
    case 3: u1 = spp->u1; v1 = spp->v1; u2 = spp->u0; v2 = spp->v1;break;
    case 4: u1 = spp->u0; v1 = spp->v1; u2 = spp->u0; v2 = spp->v0;break;
    default: g_assert_not_reached ();
    }

    if ( (u1 != uc || u2 != uc) && (v1 != vc || v2 != vc ) ) {

      /* gdouble kappa = fabs(u1*v2 + u2*vc + uc*v1 - u2*v1 - uc*v2 - u1*vc); */
      gdouble kappa = ((u1-uc)*(v2-v1) - (u2-u1)*(v1-vc));

      for ( m = 0; m < ng; m++ ) {
	gdouble beta = ui[m];

	for ( n = 0; n < ng; n++ ) {
	  gdouble mm = /* 110. */2.;
	  //gdouble alpha = ui[n];
	  gdouble alpha = 1./(mm-1.)*( (1. - pow(ui[n],mm))/(1. - ui[n]) - 1. );

	  gdouble Jm = 1./(mm-1.)*( (1. - pow(ui[n],mm))/pow(1-ui[n],2.) - mm*pow(ui[n],mm-1)/(1.-ui[n]));

	  gdouble u = uc + (u1-uc)*alpha + (u2-u1)*alpha*beta;
	  gdouble v = vc + (v1-vc)*alpha + (v2-v1)*alpha*beta;

	  /* fprintf (stdout, "%f %f \n", u, v); */

	  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &istart, &iend, sp->w_u, sp->wd_u);
	  gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);

	  Vector N = spline2d_dimensional_normal (sp, u, v);
	  // Physical distance to Gauss-Legendre Point
	  Point pg = spline2d_eval_point (sp, u, v);
	  Vector R;
	  R.x = pg.x - p.x;
	  R.y = pg.y - p.y;
	  R.z = pg.z - p.z;
	  gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z;
	  
	  // Jacobian at Gauss-Legendre point = norm of N
	  gdouble J = vector_norm (N);
	  gdouble c1 = wi[m]*wi[n]/* *J *//sqrt(r2)*alpha*kappa*Jm;


	  /* N.x /= J; N.y /= J; N.z /= J; */

	  gdouble c2 = c1*vector_scalar_product (&N, &R)/r2;

	  c1 *= J;

	  // Loop over the splines included in spp
	  for ( ii = 0; ii < k; ii++) {
	    gdouble c3 = c1*gsl_matrix_get (Bu, ii, 0);
	    gdouble c4 = c2*gsl_matrix_get (Bu, ii, 0);
	    for ( jj = 0; jj < k; jj++) {
	      psi[ii + jj*k] += c3*gsl_matrix_get (Bv, jj, 0);
	      phi[ii + jj*k] += c4*gsl_matrix_get (Bv, jj, 0);
	    }
	  }

	}
      }
    }
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  gsl_integration_glfixed_table_free (itable);
  /* gsl_integration_glfixed_table_free (jtable); */
  /* g_assert_not_reached (); */


   /*  for ( m = 0; m < k; m++) { */
  /*   for ( n = 0; n < k; n++) { */
  /*     fprintf (stdout, "%e %e \n", psi[m + n*k], phi[m + n*k]); */
  /*   } */
  /* } */

  /* g_assert_not_reached (); */

}


/* void qin_self_influence_coefficients_LW (SPPanel * spp, */
/* 					 gdouble up, gdouble vp, */
/* 					 Point p, */
/* 					 gdouble * psi, */
/* 					 gdouble * phi) */
/* { */
/*   Spline2D * sp = spp->sp; */
/*   gint i, j, m, n, se, ii, jj; */
/*   gint k = spp->k; */
/*   gint ng = 16.; */
/*   gdouble ui[ng], wi[ng]; */
/*   /\* gdouble vj[ng], wj[ng]; *\/ */
/*   gsl_integration_glfixed_table * itable =  gsl_integration_glfixed_table_alloc (ng); */
/*   /\* gsl_integration_glfixed_table * jtable =  gsl_integration_glfixed_table_alloc (ng); *\/ */

/*   for ( i = 0; i < ng; i++) */
/*     gsl_integration_glfixed_point (0., 1., i, &ui[i], &wi[i], itable); */

/*   gsl_matrix * Bu = gsl_matrix_alloc (k, 2); */
/*   gsl_matrix * Bv = gsl_matrix_alloc (k, 2); */
/*   size_t istart, iend, jstart, jend; */

/*   /\* 10 Gauss-points over [0:1] *\/ */

/*   /\* NB: G is the Jacobian of the surface *\/ */

/*   /\* *\/ */
/*   gdouble uc, vc; */
/*   //Point pc = spline2d_eval_point (sp, uc, vc); */

/*   // For now */
/*   Point pc = p; */
/*   uc = up; */
/*   vc = vp; */

/*   /\* gdouble d[3]; *\/ */
/*   /\* d[0] = pc.x - p.x; *\/ */
/*   /\* d[1] = pc.y - p.y; *\/ */
/*   /\* d[2] = pc.z - p.z; *\/ */

/*   /\* gdouble r0_2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; *\/ */

/*   /\* Vector xu, xv; *\/ */

/*   /\* xu.x = spline2d_derivative_eval (sp, uc, vc, 1, 0, 0); *\/ */
/*   /\* xu.y = spline2d_derivative_eval (sp, uc, vc, 1, 0, 1); *\/ */
/*   /\* xu.z = spline2d_derivative_eval (sp, uc, vc, 1, 0, 2); *\/ */

/*   /\* xv.x = spline2d_derivative_eval (sp, uc, vc, 0, 1, 0); *\/ */
/*   /\* xv.y = spline2d_derivative_eval (sp, uc, vc, 0, 1, 1); *\/ */
/*   /\* xv.z = spline2d_derivative_eval (sp, uc, vc, 0, 1, 2); *\/ */


/*   for (se = 1 ; se <= 4; se++) { */
/*     /\* Other ends of triangle *\/ */
/*     gdouble u1, v1; */
/*     gdouble u2, v2; */

/*     switch (se) { */
/*     case 1: u1 = spp->u0; v1 = spp->v0; u2 = spp->u1; v2 = spp->v0; break; */
/*     case 2: u1 = spp->u1; v1 = spp->v0; u2 = spp->u1; v2 = spp->v1;break; */
/*     case 3: u1 = spp->u1; v1 = spp->v1; u2 = spp->u0; v2 = spp->v1;break; */
/*     case 4: u1 = spp->u0; v1 = spp->v1; u2 = spp->u0; v2 = spp->v0;break; */
/*     default: g_assert_not_reached (); */
/*     } */

/*     if ( (u1 != uc || u2 != uc) && (v1 != vc || v2 != vc ) ) { */

/*       /\* gdouble kappa = fabs(u1*v2 + u2*vc + uc*v1 - u2*v1 - uc*v2 - u1*vc); *\/ */
/*       gdouble kappa = ((u1-uc)*(v2-v1) - (u2-u1)*(v1-vc)); */

/*       for ( m = 0; m < ng; m++ ) { */
/* 	gdouble beta = ui[m]; */

/* 	for ( n = 0; n < ng; n++ ) { */
/* 	  gdouble alpha = ui[n]; */

/* 	  gdouble u = uc + (u1-uc)*alpha + (u2-u1)*alpha*beta; */
/* 	  gdouble v = vc + (v1-vc)*alpha + (v2-v1)*alpha*beta; */

/* 	  gsl_bspline_deriv_eval_nonzero (u, 1, Bu, &istart, &iend, sp->w_u, sp->wd_u); */
/* 	  gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v); */

/* 	  Vector N = spline2d_dimensional_normal (sp, u, v); */
/* 	  // Physical distance to Gauss-Legendre Point */
/* 	  Point pg = spline2d_eval_point (sp, u, v); */
/* 	  Vector R; */
/* 	  R.x = pg.x - p.x; */
/* 	  R.y = pg.y - p.y; */
/* 	  R.z = pg.z - p.z; */
/* 	  gdouble r2 = R.x*R.x + R.y*R.y + R.z*R.z; */
	  
/* 	  // Jacobian at Gauss-Legendre point = norm of N */
/* 	  gdouble J = vector_norm (N); */
/* 	  gdouble c1 = wi[m]*wi[n]*J/sqrt(r2)*alpha*kappa; */


/* 	  N.x /= J; N.y /= J; N.z /= J; */

/* 	  gdouble c2 = c1*vector_scalar_product (&N, &R)/r2; */

/* 	  /\* c1 *= J; *\/ */

/* 	  // Loop over the splines included in spp */
/* 	  for ( ii = 0; ii < k; ii++) { */
/* 	    gdouble c3 = c1*gsl_matrix_get (Bu, ii, 0); */
/* 	    gdouble c4 = c2*gsl_matrix_get (Bu, ii, 0); */
/* 	    for ( jj = 0; jj < k; jj++) { */
/* 	      psi[ii + jj*k] += c3*gsl_matrix_get (Bv, jj, 0); */
/* 	      phi[ii + jj*k] += c4*gsl_matrix_get (Bv, jj, 0); */
/* 	    } */
/* 	  } */

/* 	} */
/*       } */
/*     } */
/*   } */

/*   gsl_matrix_free (Bu); */
/*   gsl_matrix_free (Bv); */

/*   gsl_integration_glfixed_table_free (itable); */
/*   /\* gsl_integration_glfixed_table_free (jtable); *\/ */
/*   /\* g_assert_not_reached (); *\/ */


/*    /\*  for ( m = 0; m < k; m++) { *\/ */
/*   /\*   for ( n = 0; n < k; n++) { *\/ */
/*   /\*     fprintf (stdout, "%e %e \n", psi[m + n*k], phi[m + n*k]); *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */

/*   /\* g_assert_not_reached (); *\/ */

/* } */

/**
 * Return the closest point in the panel
 **/
Point newton_raphson (SPPanel * spp, Point p)
{
  Spline2D * sp = spp->sp;
  gint k = sp->k;

  Point pc;
  gint i, j;

  FILE * fp = fopen ("panel.tmp","w");
  
  for ( i = 0; i < 20; i++) {
    gdouble upp = spp->u0 + i/(19.)*(spp->u1-spp->u0);
    for ( j = 0; j < 20; j++) {
      gdouble vpp = spp->v0 + j/(19.)*(spp->v1-spp->v0);
      Point pp = spline2d_eval_point (sp, upp, vpp);
      fprintf (fp, "%f %f %f \n", pp.x, pp.y, pp.z);
    }
  }
  
  fclose (fp);

  gdouble uc, vc;

  /* Evaluates gradient at point */
  //uc += uc 

  gsl_matrix * Bu = gsl_matrix_alloc (k, 2);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 2);
  size_t istart, iend, jstart, jend;

  gdouble u0, v0, u1, v1;

  /* First guess */
  u1 = u0 = 0.5*(spp->u0+spp->u1);
  v1 = v0 = 0.5*(spp->v0+spp->v1);

  /* We do 5 iterations */

  for ( i = 0; i < 150; i++) {
    Point p0;
    Vector xu, xv;
    gint  a, b;

    u0 = u1;
    v0 = v1;

    xu.x = xu.y = xu.z = xv.x = xv.y = xv.z = 0.;
    p0.x = p0.y = p0.z = 0.;

    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u0)), 1, Bu, &istart, &iend, sp->wx_u, sp->wxd_u);
    gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v0)), 1, Bv, &jstart, &jend, sp->w_v, sp->wd_v);

    if (sp->periodic)
      istart -= (k-1);

    for ( a = 0; a < k; a++) {
      gdouble cu = gsl_matrix_get (Bu, a, 0);
      gdouble cdu = gsl_matrix_get (Bu, a, 1);
      gint ii = (istart+a);
      for ( b = 0; b < k; b++) {
	gdouble cv = gsl_matrix_get (Bv, b, 0);
	gdouble cudv = cu*gsl_matrix_get (Bv, b, 1);
	gdouble cvdu = cv*cdu;
	gdouble cuv = cu*cv;
	gint jj = (jstart+b);
	gdouble vx = coeff (sp, ii, jj, 0);
	gdouble vy = coeff (sp, ii, jj, 1);
	gdouble vz = coeff (sp, ii, jj, 2);
      
	xu.x += vx*cvdu;
	xu.y += vy*cvdu;
	xu.z += vz*cvdu;
	xv.x += vx*cudv;
	xv.y += vy*cudv;
	xv.z += vz*cudv;
	p0.x += vx*cuv;
	p0.y += vy*cuv;
	p0.z += vz*cuv;
      }
    }

    gdouble dx = (p0.x-p.x);
    gdouble dy = (p0.y-p.y);
    gdouble dz = (p0.z-p.z);

    gdouble f = dx*dx+dy*dy+dz*dz;
    gdouble dfdu = 2.*(xu.x*dx + xu.y*dy + xu.z*dz);
    gdouble dfdv = 2.*(xv.x*dx + xv.y*dy + xv.z*dz);

    /* u1 = u0 - f(u0,v0)/d f/du; */
    /* v1 = v0 - f(u0,v0)/d f/dv; */
    u1 = u0 - 0.25*f/dfdu;
    v1 = v0 - 0.25*f/dfdv;

    /* Make sure we stay within the panel */
    if ( u1 > spp->u1 )
      u1 = spp->u1;
    if ( u1 < spp->u0 )
      u1 = spp->u0;

    if ( v1 > spp->v1 )
      v1 = spp->v1;
    if ( v1 < spp->v0 )
      v1 = spp->v0;

    pc = spline2d_eval_point (sp, u1,v1);
    fprintf (stdout, "%e \n", (pc.x-p.x)*(pc.x-p.x) + (pc.y-p.y)*(pc.y-p.y) + (pc.z-p.z)*(pc.z-p.z));
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  //pc = spline2d_eval_point (sp, u1,v1);
  //fprintf (stdout, "%e \n", (pc.x-p.x)*(pc.x-p.x) + (pc.y-p.y)*(pc.y-p.y) + (pc.z-p.z)*(pc.z-p.z));
  fprintf (stdout, "%f %f %f \n", p.x, p.y, p.z);
  fprintf (stdout, "%e %e %e \n", pc.x, pc.y, pc.z);

  g_assert_not_reached ();

  return pc;
}

/*****************************************************/
/*                 CLEAN FORMULATION                 */
/*****************************************************/

gdouble zero_spline2d_func (Spline2D * sp,
			    gdouble u, gdouble v,
			    gpointer data)
{
  return 0.;
}

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
#include "boundaries.h"
#include "surface.h"

/*    Mesh generation */
#define TEST_2D

static gdouble ramp_function (WaveParams * wp, gdouble t)
{
  gint np = 3.;

  if (t < np*2.*M_PI/wp->w)
    return sin (0.25*t/(np/wp->w));
  else
    return 1.;

  return MIN(1.,t/(2.*M_PI/wp->w));
  return MIN(1.,t/(10.*M_PI/wp->w));
}

static gdouble ramp_function_dt (WaveParams * wp, gdouble t)
{
  if (t < 3.*2.*M_PI/wp->w)
    return 0.25/(3./wp->w)*cos (0.25*t/(3./wp->w));
  else
    return 0.;
}

/**
 * h = water depth
 * t = time
 * g = gravity
 * w = frequency
 **/

gdouble finite_depth_wave_potential (WaveParams * wp, Point p, gdouble t)
{
  /* complex double phi; */

  /* /\* Complex representation *\/ */
  /* phi = wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* return creal(phi); */

  return MIN(1.,t/(10.*M_PI/wp->w))*wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*
    sin (p.x*wp->k*wp->cosb+wp->k*p.y*wp->sinb-wp->w*t);
}

gdouble finite_depth_wave_potential_dt (WaveParams * wp, Point p, gdouble t)
{
  /* complex double dphidt; */

  /* /\* Complex representation *\/ */
  /* dphidt = I*wp->w*wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* return creal(dphidt); */

  //return -wp->g*wp->A*cos (p.x*wp->k-wp->w*t);

  //-wp->g*wp->A*cos(wp->k*(x*wp->cosb+y*wp->sinb-wp->w*t);

  //  return -wp->g*wp->A*cos(wp->k*p.x-wp->w*t) -wp->g*wp->A*cos(2.115565e-2*p.x-2.*M_PI/60.*t);

  return MIN(1.,t/(10.*M_PI/wp->w))*-wp->g*wp->A*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*
    cos (p.x*wp->k*wp->cosb+wp->k*p.y*wp->sinb-wp->w*t);
}

Vector finite_depth_wave_potential_gradient (WaveParams * wp, Point p, gdouble t)
{
  /* complex double phi = wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  Vector grad;

  complex double c1 = wp->k*wp->g*wp->A/wp->w/cosh(wp->k*wp->h)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  grad.x = creal (-I*c1*wp->cosb*cosh(wp->k*(p.z+wp->h)));
  grad.y = creal (-I*c1*wp->sinb*cosh(wp->k*(p.z+wp->h)));
  grad.z = creal (c1*sinh(wp->k*(p.z+wp->h)));

  /* grad.x = creal (-I*wp->k*wp->cosb*phi); */
  /* grad.y = creal (-I*wp->k*wp->sinb*phi); */
  /* grad.z = creal (wp->g*wp->A/wp->w*wp->k*sinh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t)); */

  grad.x = MIN(1.,t/(10.*M_PI/wp->w))*wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*wp->k*cos (p.x*wp->k*wp->cosb+wp->k*p.y*wp->sinb-wp->w*t)*wp->cosb;
  grad.y = MIN(1.,t/(10.*M_PI/wp->w))*wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*wp->k*cos (p.x*wp->k*wp->cosb+wp->k*p.y*wp->sinb-wp->w*t)*wp->sinb;
  grad.z = MIN(1.,t/(10.*M_PI/wp->w))*wp->g*wp->A/wp->w*wp->k*sinh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*
    sin (p.x*wp->k-wp->w*t);

  return grad;
}

gdouble finite_depth_wave_potential_dz_dt (WaveParams * wp, Point p, gdouble t)
{
  complex double c1 = wp->k*wp->g*wp->A/wp->w/cosh(wp->k*wp->h)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  return MIN(1.,t/(10.*M_PI/wp->w))*creal (I*wp->w*c1*sinh(wp->k*(p.z+wp->h)));
}

Vector finite_depth_wave_potential_z_derivative_gradient (WaveParams * wp,
							  Point p, gdouble t)
{
  Vector grad;

  complex double c1 = wp->k*wp->k*wp->g*wp->A/wp->w/cosh(wp->k*wp->h)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  grad.x = MIN(1.,t/(10.*M_PI/wp->w))*creal (-I*c1*wp->cosb*sinh(wp->k*(p.z+wp->h)));
  grad.y = MIN(1.,t/(10.*M_PI/wp->w))*creal (-I*c1*wp->sinb*cosh(wp->k*(p.z+wp->h)));
  grad.z = MIN(1.,t/(10.*M_PI/wp->w))*creal (c1*cosh(wp->k*(p.z+wp->h)));

  return grad;
}

gdouble finite_depth_wave_elevation (gdouble x, gdouble y, gdouble t, /* WaveParams * wp */gpointer data)
{
  WaveParams * wp = (WaveParams *) data;

  /* g_assert_not_reached (); */

  //return creal(wp->g*wp->A/(wp->w*wp->w)*sinh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
  /* return creal (wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t)); */

  /* return wp->A*cos(/\* wp->k*x- *\/wp->w*t); */

  

  return MIN(1.,t/(10.*M_PI/wp->w))*wp->A*cos(wp->k*x*wp->cosb+wp->k*y*wp->sinb-wp->w*t);

  //return wp->A*cos(wp->k*(x*wp->cosb+y*wp->sinb-wp->w*t);

  return wp->A*cos(wp->k*x-wp->w*t) + wp->A*cos(2.115565e-2*x-2.*M_PI/60.*t);


  return creal (-I*wp->A*cexp(-I*wp->k*x*wp->cosb - I*wp->k*y*wp->sinb) * cexp(I*wp->w*t));
}

Vector finite_depth_wave_elevation_gradient (WaveParams * wp, Point p, gdouble t)
{
  Vector grad;

  grad.x = creal (-wp->k*wp->cosb*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
  grad.y = creal (-wp->k*wp->sinb*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
  grad.z = 0;

  grad.x = MIN(1.,t/(10.*M_PI/wp->w))*-wp->k*wp->A*sin(wp->k*p.x-wp->k*t);;
  grad.y = MIN(1.,t/(10.*M_PI/wp->w))*0.;
  grad.z = MIN(1.,t/(10.*M_PI/wp->w))*0.;

  return grad;
}

gdouble finite_depth_wave_elevation_time_derivative (WaveParams * wp, Point p, gdouble t)
{
  return MIN(1.,t/(10.*M_PI/wp->w))*wp->A*wp->k*sin(wp->k*p.x-wp->k*t);
  return creal (wp->w*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
}

Vector finite_depth_wave_normal_time_derivative (WaveParams * wp, Point p, gdouble t)
{
  complex double B = wp->k*wp->A*cexp( I*wp->w*t - I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb);
  complex double c = I*wp->w*B/cpow(1.+B*B, 3./2.);
  Vector nt;

  nt.x = MIN(1.,t/(10.*M_PI/wp->w))*creal (c*wp->cosb);
  nt.y = MIN(1.,t/(10.*M_PI/wp->w))*creal (c*wp->sinb);
  nt.z = MIN(1.,t/(10.*M_PI/wp->w))*creal (-B*c);

  return nt;
}

gdouble finite_depth_wave_potential2 (WaveParams * wp, Point p, gdouble t)
{
  gdouble phi0_alpha2 = wp->A*wp->A*wp->w*3./8.*cosh(wp->k*(p.z+wp->h))/pow(sinh(wp->k*wp->h),4.);

  return phi0_alpha2*cosh(2.*wp->k*(p.z+wp->h))*sin(2.*(p.x*wp->k-wp->w*t));

  g_assert_not_reached ();
  /* complex double phi; */

  /* /\* Complex representation *\/ */
  /* phi = wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* return creal(phi); */

  return wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*
    sin (p.x*wp->k-wp->w*t);
}

gdouble finite_depth_wave_potential2_dt (WaveParams * wp, Point p, gdouble t)
{
  gdouble phi0_alpha2 = wp->A*wp->A*wp->w*3./8.*cosh(wp->k*(p.z+wp->h))/pow(sinh(wp->k*wp->h),4.);

  return -2.*wp->w*phi0_alpha2*cosh(2*wp->k*(p.z+wp->h))*sin(2.*(p.x*wp->k-wp->w*t));

  g_assert_not_reached ();
  /* complex double dphidt; */

  /* /\* Complex representation *\/ */
  /* dphidt = I*wp->w*wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* return creal(dphidt); */

  //return -wp->g*wp->A*cos (p.x*wp->k-wp->w*t);

  //-wp->g*wp->A*cos(wp->k*(x*wp->cosb+y*wp->sinb-wp->w*t);

  //  return -wp->g*wp->A*cos(wp->k*p.x-wp->w*t) -wp->g*wp->A*cos(2.115565e-2*p.x-2.*M_PI/60.*t);

  return -wp->g*wp->A*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*
    cos (p.x*wp->k-wp->w*t);
}

Vector finite_depth_wave_potential2_gradient (WaveParams * wp, Point p, gdouble t)
{
  Vector grad;

  gdouble phi0_alpha2 = wp->A*wp->A*wp->w*3./8.*cosh(wp->k*(p.z+wp->h))/pow(sinh(wp->k*wp->h),4.);

  grad.x = 2*wp->k*phi0_alpha2*cosh(wp->k*(p.z+wp->h))*cos(p.x*wp->k-wp->w*t);
  grad.y = 0.;
  grad.z = 2*wp->k*phi0_alpha2*sinh(wp->k*(p.z+wp->h))*sin(p.x*wp->k-wp->w*t);

  return grad;
  g_assert_not_reached ();
  /* complex double phi = wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  //Vector grad;

  complex double c1 = wp->k*wp->g*wp->A/wp->w/cosh(wp->k*wp->h)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  grad.x = creal (-I*c1*wp->cosb*cosh(wp->k*(p.z+wp->h)));
  grad.y = creal (-I*c1*wp->sinb*cosh(wp->k*(p.z+wp->h)));
  grad.z = creal (c1*sinh(wp->k*(p.z+wp->h)));

  /* grad.x = creal (-I*wp->k*wp->cosb*phi); */
  /* grad.y = creal (-I*wp->k*wp->sinb*phi); */
  /* grad.z = creal (wp->g*wp->A/wp->w*wp->k*sinh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t)); */

  grad.x = wp->g*wp->A/wp->w*cosh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*wp->k*
    cos (p.x*wp->k-wp->w*t);
  grad.y = 0.;
  grad.z = wp->g*wp->A/wp->w*wp->k*sinh(wp->k*(p.z+wp->h))/cosh(wp->k*wp->h)*
    sin (p.x*wp->k-wp->w*t);

  return grad;
}

gdouble finite_depth_wave_potential2_dz_dt (WaveParams * wp, Point p, gdouble t)
{
  g_assert_not_reached ();
  complex double c1 = wp->k*wp->g*wp->A/wp->w/cosh(wp->k*wp->h)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  return creal (I*wp->w*c1*sinh(wp->k*(p.z+wp->h)));
}

Vector finite_depth_wave_potential2_z_derivative_gradient (WaveParams * wp,
							  Point p, gdouble t)
{
  Vector grad;
  g_assert_not_reached ();
  complex double c1 = wp->k*wp->k*wp->g*wp->A/wp->w/cosh(wp->k*wp->h)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  grad.x = creal (-I*c1*wp->cosb*sinh(wp->k*(p.z+wp->h)));
  grad.y = creal (-I*c1*wp->sinb*cosh(wp->k*(p.z+wp->h)));
  grad.z = creal (c1*cosh(wp->k*(p.z+wp->h)));

  return grad;
}

gdouble finite_depth_wave_elevation2 (gdouble x, gdouble y, gdouble t, /* WaveParams * wp */gpointer data)
{
  WaveParams * wp = (WaveParams *) data;
  gdouble alpha = tanh(wp->k*wp->h);

  /* (Shao, 2010) p 123 */

  return 0.25*(3.-alpha*alpha)/(alpha*alpha*alpha)*wp->k*wp->A*wp->A*cos(2.*(wp->k*x-wp->w*t));
}

Vector finite_depth_wave_elevation2_gradient (WaveParams * wp, Point p, gdouble t)
{
  Vector grad;
  g_assert_not_reached ();
  grad.x = creal (-wp->k*wp->cosb*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
  grad.y = creal (-wp->k*wp->sinb*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
  grad.z = 0;

  grad.x = -wp->k*wp->A*sin(wp->k*p.x-wp->k*t);;
  grad.y = 0.;
  grad.z = 0.;

  return grad;
}

gdouble finite_depth_wave_elevation2_time_derivative (WaveParams * wp, Point p, gdouble t)
{
  gdouble alpha = tanh(wp->k*wp->h);

  /* (Shao, 2010) p 123 */

  return 2.*wp->w*0.25*(3.-alpha*alpha)/(alpha*alpha*alpha)*wp->k*wp->A*wp->A*sin(2.*(wp->k*p.x-wp->w*t));

  g_assert_not_reached ();
  return wp->A*wp->k*sin(wp->k*p.x-wp->k*t);
  return creal (wp->w*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
}

Vector finite_depth_wave_normal2_time_derivative (WaveParams * wp, Point p, gdouble t)
{
  g_assert_not_reached ();
  complex double B = wp->k*wp->A*cexp( I*wp->w*t - I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb);
  complex double c = I*wp->w*B/cpow(1.+B*B, 3./2.);
  Vector nt;

  nt.x = creal (c*wp->cosb);
  nt.y = creal (c*wp->sinb);
  nt.z = creal (-B*c);

  return nt;
}

gdouble infinite_depth_wave_potential (WaveParams * wp, Point p, gdouble t)
{
  /* complex double phi; */

  /* /\* Complex representation *\/ */
  /* phi = wp->g*wp->A/wp->w*exp(wp->k*p.z)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* return creal(phi); */
  return ramp_function (wp, t)*wp->g*wp->A/wp->w*exp(wp->k*p.z)*cos(wp->w*t - wp->k*p.x*wp->cosb - wp->k*p.y*wp->sinb);
}

gdouble infinite_depth_wave_potential_dt (WaveParams * wp, Point p, gdouble t)
{
  /* complex double dtphi; */

  /* dtphi = I*wp->w*wp->g*wp->A/wp->w*exp(wp->k*p.z)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* return creal(dtphi); */
  return -ramp_function (wp, t)*wp->g*wp->A*exp(wp->k*p.z)*sin(wp->w*t - wp->k*p.x*wp->cosb - wp->k*p.y*wp->sinb);
}

Vector infinite_depth_wave_potential_gradient (WaveParams * wp, Point p, gdouble t)
{
  Vector grad;

  /* complex double c1 = wp->k*wp->g*wp->A/wp->w*exp(wp->k*p.z)* */
  /*   cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t); */

  /* grad.x = creal (-I*c1*wp->cosb); */
  /* grad.y = creal (-I*c1*wp->sinb); */
  /* grad.z = creal (c1); */

  gdouble c0 =  ramp_function (wp, t)*wp->k*wp->g*wp->A/wp->w*exp(wp->k*p.z);
  gdouble theta = wp->w*t-wp->k*p.x*wp->cosb-wp->k*p.y*wp->sinb;

  grad.x = wp->cosb*c0*sin(theta);
  grad.y = wp->sinb*c0*sin(theta);
  grad.z = c0*cos(theta);

  return grad;
}

gdouble infinite_depth_wave_potential_dz_dt (WaveParams * wp, Point p, gdouble t)
{
  complex double phi_dz_dt;
  g_assert_not_reached ();
  phi_dz_dt = ramp_function (wp, t)*I*wp->w*wp->k*wp->g*wp->A/wp->w*exp(wp->k*p.z)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  return creal (phi_dz_dt);
}

Vector infinite_depth_wave_potential_z_derivative_gradient (WaveParams * wp,
							    Point p, gdouble t)
{
  Vector grad;
  g_assert_not_reached ();
  complex double c1 = ramp_function (wp, t)*wp->k*wp->k*wp->g*wp->A/wp->w*exp(wp->k*p.z)*
    cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  grad.x = creal (-I*c1*wp->cosb);
  grad.y = creal (-I*c1*wp->sinb);
  grad.z = creal (c1);

  return grad;
}

//gdouble infinite_depth_wave_elevation (WaveParams * wp, Point p, gdouble t)
gdouble infinite_depth_wave_elevation (gdouble x, gdouble y, gdouble t, gpointer data)
{
  WaveParams * wp = (WaveParams *) data;
  //return creal (-I*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t));
  return ramp_function (wp, t)*wp->A*sin(wp->w*t-wp->k*x*wp->cosb-wp->k*y*wp->sinb);
}

Vector infinite_depth_wave_elevation_gradient (WaveParams * wp, Point p, gdouble t)
{
  Vector grad;
  g_assert_not_reached ();
  complex double c1 = -ramp_function (wp, t)*I*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);

  grad.x = -I*wp->k*wp->cosb*c1;
  grad.y = - I*wp->k*wp->sinb*c1;
  grad.z = 0.;

  return grad;
}

gdouble infinite_depth_wave_elevation_time_derivative (WaveParams * wp, Point p, gdouble t)
{
  complex double dteta;
  g_assert_not_reached ();
  dteta = wp->w*wp->A*cexp(-I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb) * cexp(I*wp->w*t);
  g_assert_not_reached ();
  return creal (dteta);
}

Vector infinite_depth_wave_normal_time_derivative (WaveParams * wp, Point p, gdouble t)
{
  complex double B = wp->k*wp->A*cexp( I*wp->w*t - I*wp->k*p.x*wp->cosb - I*wp->k*p.y*wp->sinb);
  complex double c = ramp_function (wp, t)*I*wp->w*B/cpow(1.+B*B, 3./2.);
  Vector nt;
  g_assert_not_reached ();
  nt.x = creal (c*wp->cosb);
  nt.y = creal (c*wp->sinb);
  nt.z = creal (-B*c);

  return nt;
}


static gdouble f_dxdu (gdouble u, gdouble v)
{
  return 1.;
}

static gdouble f_dxdv (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dxdudu (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dxdvdv (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dxdudv (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dydu (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dydv (gdouble u, gdouble v)
{
  return 1.;
}

static gdouble f_dydudu (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dydvdv (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_dydudv (gdouble u, gdouble v)
{
  return 0.;
}

static gdouble f_z (gdouble u, gdouble v, gdouble t, gpointer data)
{
  return 0.;
  return cos(u/5.);
#ifdef TEST_2D
  return 0.;
#else
  //  return sin(M_PI*u)*sin(M_PI*v);
  return -0.015 + 0.005*cos(5.*u)*sin(10.*v);
  return 1.0 - 0.5*(u*u*(3.-2.*u)*(1.-sin(M_PI*(0.5-v))));
#endif
}

static gdouble f_dzdu (gdouble u, gdouble v)
{
#ifdef TEST_2D
  return 0.;
#else
  // return M_PI*cos(M_PI*u)*sin(M_PI*v);
  return -0.025*sin(5.*u)*sin(10.*v);
  return 3.*(u*u - u)*(1.+sin(M_PI*(0.5-v)));
#endif
}

static gdouble f_dzdv (gdouble u, gdouble v)
{
#ifdef TEST_2D
  return 0.;
#else
  //return M_PI*sin(M_PI*u)*cos(M_PI*v);
  return 0.05*cos(5.*u)*cos(10.*v);
  return -0.5*(u*u*(3.-2.*u))*M_PI*cos(M_PI*(0.5-v));
#endif
}

static gdouble f_dzdudu (gdouble u, gdouble v)
{
#ifdef TEST_2D
  return 0.;
#else
  //return -M_PI*M_PI*sin(M_PI*u)*sin(M_PI*v);
  return -5.*0.025*cos(5.*u)*sin(10.*v);
  return 3.*(2.*u - 1.)*(1.+sin(M_PI*(0.5-v)));
#endif
}

static gdouble f_dzdvdv (gdouble u, gdouble v)
{
#ifdef TEST_2D
  return 0.;
#else
  //return -M_PI*M_PI*sin(M_PI*u)*sin(M_PI*v);
  return -0.05*cos(5.*u)*sin(10.*v);
  return 0.5*(u*u*(3.-2.*u))*M_PI*M_PI*sin(M_PI*(0.5-v));
#endif
}

static gdouble f_dzdudv (gdouble u, gdouble v)
{
#ifdef TEST_2D
  return 0.;
#else
  //return -M_PI*M_PI*cos(M_PI*u)*cos(M_PI*v);
  return -5.*0.05*sin(5.*u)*cos(10.*v);
  return 3.*(u*u - u)*(M_PI*cos(M_PI*(0.5-v)));
#endif
}

/********************************************************************/

/* Methods for Surface */

Surface * surface_new ()
{
  Surface * new = g_malloc (sizeof(Surface));
  new->b = boundaries_new ();
  new->patches = NULL;
  new->hz = NULL;
  return new;
}

void surface_destroy (Surface * s)
{
  g_assert (s != NULL);
  GSList * patches = s->patches;
  while (patches) {
    spline2d_destroy (patches->data);
  }
  g_slist_free (s->patches);
  boundaries_destroy (s->b);
  g_free (s);
}

gdouble surface_eval (Surface * s, gdouble x, gdouble y)
{
  gdouble sum = 0.;

  GSList * patches = s->patches;
  while (patches) {
    sum += patch_eval (patches->data, x, y, 2);
    patches = patches->next;
  }

  return sum;
}

void surface_print_grid (Surface * s, FILE * fp)
{
  patch_print (s->patches->data, fp);
}

void spline2d_surface_print_grid (Surface * s, FILE * fp)
{
  spline2d_print_panels (s->patches->data, fp);
}

static Point freesurface_top (gdouble u, gdouble t, gpointer * data)
{
  Point p;
  Boundaries * b = data[0];
  WaveParams * wp = data[1];

  g_assert ( u >= 0. && u <= 1. );

  p.x = -wp->r1*cos((u)*2.*M_PI);
  p.y = wp->r2*sin((u)*2.*M_PI);
  p.z = wp->wave_elevation (p.x, p.y, t, wp);
  p.xi = u;

  return p;
}

static Point freesurface_left (gdouble v, gdouble t, gpointer * data)
{
  Point p;
  Boundaries * b = data[0];
  Point p0 = g_array_index (b->dcb->p, Point, 0);
  Point p1 = g_array_index (b->dct->p, Point, 0);
  WaveParams * wp = data[1];

  g_assert ( v >= 0. && v <= 1. );

  p.x = (p1.x - p0.x)*v + p0.x;
  p.y = v;
  p.z = wp->wave_elevation (p.x, p.y, t, data[1]);
  p.xi = v;

  return p;
}

static Point freesurface_right (gdouble v, gdouble t, gpointer * data)
{
  Point p;
  Boundaries * b = data[0];
  Point p0 = g_array_index (b->dcb->p, Point, 0);
  Point p1 = g_array_index (b->dct->p, Point, 0);
  WaveParams * wp = data[1];

  g_assert ( v >= 0. && v <= 1. );

  p.x = (p1.x - p0.x)*v + p0.x;
  p.y = v;
  p.z = wp->wave_elevation (p.x, p.y, t, data[1]);
  p.xi = v;

  return p;
}

/* Methods for FreeSurface */

FreeSurface * freesurface_new ()
{
  FreeSurface * new = g_malloc (sizeof(FreeSurface));
  new->s = surface_new ();
  return new;
}

void freesurface_init (FreeSurface * f, WaveParams * wp)
{
  Surface * s = f->s;

  s->b->curve_top = freesurface_top;
  s->b->curve_right = freesurface_right;
  s->b->curve_left = freesurface_left;

  s->hz = wp->wave_elevation;
}

void freesurface_destroy (FreeSurface * f)
{
  if (f) {
    surface_destroy (f->s);
    g_free (f);
  }
}

static gdouble b_z (gdouble x, gdouble y, gdouble t, gpointer data)
{
  WaveParams * wp = data;
  return 0.;
  return -wp->h;
}

static Point bathymetry_top (gdouble u, gpointer * data)
{
  Point p;
  Boundaries * b = data[0];
  WaveParams * wp = data[1];

  g_assert ( u >= 0. && u <= 1. );

  p.x = -wp->r1*cos((u)*2.*M_PI);
  p.y = wp->r2*sin((u)*2.*M_PI);
  p.z = b_z (p.x, p.y, 0, wp);
  p.xi = u;

  return p;
}

static Point bathymetry_left (gdouble v, gpointer * data)
{
  Point p;
  Boundaries * b = data[0];
  Point p0 = g_array_index (b->dcb->p, Point, 0);
  Point p1 = g_array_index (b->dct->p, Point, 0);

  g_assert ( v >= 0. && v <= 1. );

  p.x = (p1.x - p0.x)*v + p0.x;
  p.y = v;
  p.z = b_z (p.x, p.y, 0, data[1]);
  p.xi = v;

  return p;
}

static Point bathymetry_right (gdouble v, gpointer * data)
{
  Point p;
  Boundaries * b = data[0];
  Point p0 = g_array_index (b->dcb->p, Point, 0);
  Point p1 = g_array_index (b->dct->p, Point, 0);

  g_assert ( v >= 0. && v <= 1. );

  p.x = (p1.x - p0.x)*v + p0.x;
  p.y = v;
  p.z = b_z (p.x, p.y, 0, data[1]);
  p.xi = v;

  return p;
}

/* Methods for Bathymetry */

Bathymetry * bathymetry_new ()
{
  
  Bathymetry * new = g_malloc (sizeof(Bathymetry));
  new->s = surface_new ();
  return new;
}

static gdouble dc_max_x (DCurve * dc)
{
  gdouble max = -G_MAXDOUBLE;
  gint i;

  for ( i = 0; i < dc->p->len; i++) {
    Point p = g_array_index (dc->p, Point, i);
    if ( p.x > max )
      max = p.x;
  }
  return max;
}

static gdouble dc_min_x (DCurve * dc)
{
  gdouble min = G_MAXDOUBLE;
  gint i;

  for ( i = 0; i < dc->p->len; i++) {
    Point p = g_array_index (dc->p, Point, i);
    if ( p.x < min )
      min = p.x;
  }
  return min;
}

void bathymetry_init (Bathymetry * b, DCurve * dcb, WaveParams * wp)
{
  Surface * s = b->s;
  DCurve * bottom = dcurve_new (dcb->p->len);
  gint i;

  gdouble max_x = dc_max_x (dcb);
  gdouble min_x = dc_min_x (dcb);

  for ( i = 0; i < dcb->p->len/2; i++) {
    Point p;
    p.x = min_x + (max_x-min_x)*i/(dcb->p->len/2);
    p.y = 0.;
    p.z = b_z (p.x, p.y, 0, wp);
    g_array_append_val (bottom->p, p);
  }

  for ( i = dcb->p->len/2; i < dcb->p->len; i++) {
    Point p;
    p.x = max_x + (min_x-max_x)*(i-dcb->p->len/2)/(dcb->p->len-dcb->p->len/2-1);
    p.y = 0.;
    p.z = b_z (p.x, p.y, 0, wp);
    g_array_append_val (bottom->p, p);
  }

  s->b->dcb = bottom;
  s->b->curve_top = bathymetry_top;
  s->b->curve_right = bathymetry_right;
  s->b->curve_left = bathymetry_left;

  s->hz = b_z;
}

void bathymetry_destroy (Bathymetry * b)
{
  if (b) {
    surface_destroy (b->s);
    g_free (b);
  }
}

static void generate_grid_by_linear_transfinite_interpolation (Boundaries * b,
							       LinearProblem * lpx,
							       LinearProblem * lpy)
{
  gint i, j;
  gdouble u, v;
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;

  /* Copy boundary points */
  Point p;
  for ( i = 0; i < N; i++) {
    j = 0;
    p = g_array_index (b->dcb->p, Point, i);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
    j = M-1;
    p = g_array_index (b->dct->p, Point, i);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
  }
  for ( j = 0; j < M; j++) {
    i = 0;
    p = g_array_index (b->dcl->p, Point, j);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
    i = N-1;
    p = g_array_index (b->dcr->p, Point, j);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
  }

  Point p00 = g_array_index (b->dcb->p, Point, 0);
  Point pN0 = g_array_index (b->dcb->p, Point, b->dcb->p->len-1);
  Point pNN = g_array_index (b->dct->p, Point, b->dcb->p->len-1);
  Point p0N = g_array_index (b->dct->p, Point, 0);

  for ( i = 1; i < N-1; i++) {
    Point pi0 = g_array_index (b->dcb->p, Point, i);
    Point piN = g_array_index (b->dct->p, Point, i);
    for ( j = 1; j < M-1; j++) {
       Point p0j = g_array_index (b->dcl->p, Point, j);
       Point pNj = g_array_index (b->dcr->p, Point, j);
       gdouble si = i/(N-1.);
       gdouble tj = j/(M-1.);

       u = (1.-si)*p0j.x + si*pNj.x
	 + (1.-tj)*pi0.x + tj*piN.x
	 - (1.-si)*((1.-tj)*p00.x + tj*p0N.x)
	 - si*((1.-tj)*pN0.x + tj*pNN.x);

       v = (1.-si)*p0j.y + si*pNj.y
	 + (1.-tj)*pi0.y + tj*piN.y
	 - (1.-si)*((1.-tj)*p00.y + tj*p0N.y)
	 - si*((1.-tj)*pN0.y + tj*pNN.y);

       g_array_index (lpx->lhs, gdouble, i + j*N) = u;
       g_array_index (lpy->lhs, gdouble, i + j*N) = v;
    }
  }
}

/********************************************************************/

/* Some macros and functions to calculate the face coefficients f */
/* as described in (Eca, 1996) */
#define u(a,b) g_array_index (lpx->lhs, gdouble, a + (b)*N)
#define v(a,b) g_array_index (lpy->lhs, gdouble, a + (b)*N)
/* #define u(a,b) g_array_index (lp->lhs, gdouble, a + (b)*N) */
/* #define v(a,b) g_array_index (lp->lhs, gdouble, N*M + a + (b)*N) */
#define z(a,b) (f_z (u(a,b), v(a,b), 0, NULL))

#define g11_bar (f_dxdu(u(i,j),v(i,j))*f_dxdu(u(i,j),v(i,j)) + f_dydu(u(i,j),v(i,j))*f_dydu(u(i,j),v(i,j)) + f_dzdu(u(i,j),v(i,j))*f_dzdu(u(i,j),v(i,j)))
#define g12_bar (f_dxdu(u(i,j),v(i,j))*f_dxdv(u(i,j),v(i,j)) + f_dydu(u(i,j),v(i,j))*f_dydv(u(i,j),v(i,j)) + f_dzdu(u(i,j),v(i,j))*f_dzdv(u(i,j),v(i,j)))
#define g22_bar (f_dxdv(u(i,j),v(i,j))*f_dxdv(u(i,j),v(i,j)) + f_dydv(u(i,j),v(i,j))*f_dydv(u(i,j),v(i,j)) + f_dzdv(u(i,j),v(i,j))*f_dzdv(u(i,j),v(i,j)))

/*********************************************************************/

static gdouble H00 (gdouble t)
{
  return (t-1.)*(t-1.)*(2.*t+1.);
}

static gdouble H01 (gdouble t)
{
  return t*t*(3.-2.*t);
}

static gdouble H11 (gdouble t)
{
  return (t-1.)*(t-1.)*t;
}

static gdouble H10 (gdouble t)
{
  return t*t*(t-1.);
}

static gdouble u_x (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		    gint i, gint j)
{
  gint N = b->dcb->p->len;
  /* gint M = b->dcl->p->len; */

  if ( i == 0 ) {
    return (u(1,j)-u(0,j));
  }
  if ( i == N-1)  {
    return u(N-1,j)-u(N-2,j);
  }
  return (u(i+1,j)-u(i-1,j))/2.;
}

static gdouble v_x (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		    gint i, gint j)
{
  /* gint M = b->dcl->p->len; */
  gint N = b->dcb->p->len;
  return -g12_bar/g22_bar*u_x(b, lpx, lpy, i, j);
}

static gdouble v_e (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		    gint i, gint j)
{
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;

  if ( j == 0 ) {
    return v(i,1)-v(i,0);
  }
  if ( j == M-1)  {
    return v(i,M-1)-v(i,M-2);
  }
  return (v(i,j+1)-v(i,j-1))/2.;
}

static gdouble u_e (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		    gint i, gint j)
{
  gint N = b->dcb->p->len;
  /* gint M = b->dcl->p->len; */

  return -g12_bar/g11_bar*v_e(b, lpx, lpy, i, j);
}

#define us(i,j) (u_x (b,lpx,lpy,i,j))
#define vs(i,j) (v_x (b,lpx,lpy,i,j))
#define ut(i,j) (u_e (b,lpx,lpy,i,j))
#define vt(i,j) (v_e (b,lpx,lpy,i,j))

static gdouble u_ex (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		     gint i, gint j)
{
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;

  if (j == 0)
    return us(i,1)-us(i,0);
  if (j == M-1)
    return us(i,N-1)-us(i,N-2);

  return (us(i,j+1)-us(i,j-1))/2.;
}

static gdouble v_ex (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		     gint i, gint j)
{
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;

  if (j == 0)
    return vs(i,1)-vs(i,0);
  if (j == M-1)
    return vs(i,N-1)-vs(i,N-2);

  return (vs(i,j+1)-vs(i,j-1))/2.;
}

#define ust(i,j) (u_ex (b,lpx,lpy,i,j))
#define vst(i,j) (v_ex (b,lpx,lpy,i,j))

static void generate_grid_by_hermite_transfinite_interpolation (Boundaries * b,
								LinearProblem * lpx,
								LinearProblem * lpy)
{
  gint i, j;
  gdouble u, v;
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;

  /* Copy boundary points */
  Point p;
  for ( i = 0; i < N; i++) {
    j = 0;
    p = g_array_index (b->dcb->p, Point, i);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
    j = M-1;
    p = g_array_index (b->dct->p, Point, i);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
  }
  for ( j = 0; j < M; j++) {
    i = 0;
    p = g_array_index (b->dcl->p, Point, j);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
    i = N-1;
    p = g_array_index (b->dcr->p, Point, j);
    g_array_index (lpx->lhs, gdouble, i + j*N) = p.x;
    g_array_index (lpy->lhs, gdouble, i + j*N) = p.y;
  }

  Point p00 = g_array_index (b->dcb->p, Point, 0);
  Point pN0 = g_array_index (b->dcb->p, Point, b->dcb->p->len-1);
  Point pNN = g_array_index (b->dct->p, Point, b->dcb->p->len-1);
  Point p0N = g_array_index (b->dct->p, Point, 0);

  for ( i = 1; i < N-1; i++) {
    Point pi0 = g_array_index (b->dcb->p, Point, i);
    Point piN = g_array_index (b->dct->p, Point, i);
    for ( j = 1; j < M-1; j++) {
       Point p0j = g_array_index (b->dcl->p, Point, j);
       Point pNj = g_array_index (b->dcr->p, Point, j);
       gdouble si = i/(N-1.);
       gdouble tj = j/(M-1.);

       u = H00(si)*p0j.x + H10(si)*us(0,j) + H11(si)*us(N-1,j) + H01(si)*pNj.x
       	 + H00(tj)*pi0.x + H10(tj)*ut(i,0) + H11(tj)*us(i,M-1) + H01(tj)*piN.x
       	 - H00(si)*(p00.x*H00(tj) + ut(0,0)*H10(tj) + ut(0,M-1)*H11(tj) + p0N.x*H01(tj))
       	 - H10(si)*(us(0,0)*H00(tj) + ust(0,0)*H10(tj) + ust(0,M-1)*H11(tj) + us(0,M-1)*H01(tj))
       	 - H11(si)*(us(N-1,0)*H00(tj) + ust(N-1,0)*H10(tj) + ust(N-1,M-1)*H11(tj)
       	 	    + us(N-1,M-1)*H00(tj))
       	 - H01(si)*(pN0.x*H00(tj) + ut(N-1,0)*H10(tj) + ut(N-1,M-1)*H11(tj) + pNN.x*H01(tj));

       v = H00(si)*p0j.y + H10(si)*vs(0,j) + H11(si)*vs(N-1,j) + H01(si)*pNj.y
       	 + H00(tj)*pi0.y + H10(tj)*vt(i,0) + H11(tj)*vs(i,M-1) + H01(tj)*piN.y
       	 - H00(si)*(p00.y*H00(tj) + vt(0,0)*H10(tj) + vt(0,M-1)*H11(tj) + p0N.y*H01(tj))
       	 - H10(si)*(vs(0,0)*H00(tj) + vst(0,0)*H10(tj) + vst(0,M-1)*H11(tj) + vs(0,M-1)*H01(tj))
       	 - H11(si)*(vs(N-1,0)*H00(tj) + vst(N-1,0)*H10(tj) + vst(N-1,M-1)*H11(tj)
       	 	    + vs(N-1,M-1)*H01(tj))
       	 - H01(si)*(pN0.y*H00(tj) + vt(N-1,0)*H10(tj) + vt(N-1,M-1)*H11(tj) + pNN.y*H01(tj));

       g_array_index (lpx->lhs, gdouble, i + j*N) = u;
       g_array_index (lpy->lhs, gdouble, i + j*N) = v;
    }
  }

  FILE * fp = fopen("tt.tmp","w");
  for (i = 0; i < N; i++) {
    p = g_array_index (b->dcb->p, Point, i);
    fprintf(fp, "%g %g %g %g\n", p.xi,u(i,0),v(i,0), us(i,0)*(N-1));
    // u(a,b) g_array_index (lpx->lhs, gdouble, a + (b)*M)
  }
  fclose (fp);
}

/**
 * Prints the grid in a format that can be plotted using Gnuplot.
 * The command for plotting the grid is splot 'test.tmp' w lp
 */
static void print_grid (LinearProblem * lpx, LinearProblem * lpy, gint N, gint M)
{
  FILE * fout = fopen("grid.tmp", "w");
  gint i, j;
  
  for ( i = 0; i < N; i++ ) {
    for ( j = 0; j < M; j++ )
      fprintf(fout, "%f %f %f\n", g_array_index(lpx->lhs,gdouble,i+j*N),
	      g_array_index(lpy->lhs,gdouble,i+j*N), f_z (g_array_index(lpx->lhs,gdouble,i+j*N),
							  g_array_index(lpy->lhs,gdouble,i+j*N), 0, NULL));
    fprintf (fout, "\n");
  }

  for ( j = 0; j < M; j++ ) {
    for ( i = 0; i < N; i++ )
      fprintf(fout, "%f %f %f\n", g_array_index(lpx->lhs,gdouble,i+j*N),
	      g_array_index(lpy->lhs,gdouble,i+j*N), f_z (g_array_index(lpx->lhs,gdouble,i+j*N),
							  g_array_index(lpy->lhs,gdouble,i+j*N), 0, NULL));
    fprintf (fout, "\n");
  }

  fclose (fout);
}

/**
 * Calculates the maximum deviation from orthogonality
 **/
static gdouble m_do (LinearProblem * lpx,
		     LinearProblem * lpy,
		     gint N, gint M)
{
  gint i, j;
  gdouble mdo = 0.;
  gdouble max = 0.;
  /* FILE * fa = fopen ("angle.tmp","w"); */

  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++) {
      gdouble v1x = u(i+1,j) - u(i-1,j);
      gdouble v1y = v(i+1,j) - v(i-1,j);
      gdouble v1z = z(i+1,j) - z(i-1,j);
      gdouble v2x = u(i,j+1) - u(i,j-1);
      gdouble v2y = v(i,j+1) - v(i,j-1);
      gdouble v2z = z(i,j+1) - z(i,j-1);
      gdouble tmp = fabs(v1x*v2x + v1y*v2y +v1z*v2z);
      tmp /= (sqrt(v1x*v1x+v1y*v1y+v1z*v1z)*sqrt(v2x*v2x+v2y*v2y+v2z*v2z));

      tmp = fabs(acos(tmp)*180./M_PI);

      if ( tmp > 90.)
	tmp = 180. - tmp;

      tmp = 90. - tmp;

      mdo += tmp;
      
      /* fprintf(fa, "%g %g %g \n", u(i,j), v(i,j), tmp); */

      if (max < tmp)
	max = tmp;
    }
  fprintf(stdout, "%g %g\n",mdo/(N*M),max);
  /* fclose (fa); */
  return mdo/(N*N)/* max */;
}

static void wild_grid_generation (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy)
{
  gint i, j;
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;

  for (i = 0; i < N; i++) {
    Point p0 = g_array_index (b->dcb->p, Point, i);
    Point pN = g_array_index (b->dct->p, Point, i);
    for (j = 0; j < M; j++) {
      gdouble v = j/(M-1.);
      v = v < 0.7 ? 0.35/0.7*v : 0.35 + 0.65*(v-0.7)/0.3;
      u(i,j) = p0.x*(1.-v) + v*pN.x;
      v(i,j) = p0.y*(1.-v) + v*pN.y;
      /* u(i,j) = p0.x*(1.-j/(M-1.)) + j/(M-1.)*pN.x; */
      /* v(i,j) = p0.y*(1.-j/(M-1.)) + j/(M-1.)*pN.y; */
    }
  }
}

/* page 9 */
#define g11u_bar (2.*(f_dxdu(u(i,j),v(i,j))*f_dxdudu(u(i,j),v(i,j)) + f_dydu(u(i,j),v(i,j))*f_dydudu(u(i,j),v(i,j)) + f_dzdu(u(i,j),v(i,j))*f_dzdudu(u(i,j),v(i,j))))
#define g11v_bar (2.*(f_dxdu(u(i,j),v(i,j))*f_dxdudv(u(i,j),v(i,j)) + f_dydu(u(i,j),v(i,j))*f_dydudv(u(i,j),v(i,j)) + f_dzdu(u(i,j),v(i,j))*f_dzdudv(u(i,j),v(i,j))))
#define g22u_bar (2.*(f_dxdv(u(i,j),v(i,j))*f_dxdudv(u(i,j),v(i,j)) + f_dydv(u(i,j),v(i,j))*f_dydudv(u(i,j),v(i,j)) + f_dzdv(u(i,j),v(i,j))*f_dzdudv(u(i,j),v(i,j))))
#define g22v_bar (2.*(f_dxdv(u(i,j),v(i,j))*f_dxdvdv(u(i,j),v(i,j)) + f_dydv(u(i,j),v(i,j))*f_dydvdv(u(i,j),v(i,j)) + f_dzdv(u(i,j),v(i,j))*f_dzdvdv(u(i,j),v(i,j))))
#define g12u_bar (( f_dxdu(u(i,j),v(i,j))*f_dxdudv(u(i,j),v(i,j)) + f_dydu(u(i,j),v(i,j))*f_dydudv(u(i,j),v(i,j)) + f_dzdu(u(i,j),v(i,j))*f_dzdudv(u(i,j),v(i,j))) + (f_dxdv(u(i,j),v(i,j))*f_dxdudu(u(i,j),v(i,j)) + f_dydv(u(i,j),v(i,j))*f_dydudu(u(i,j),v(i,j)) + f_dzdv(u(i,j),v(i,j))*f_dzdudu(u(i,j),v(i,j))))
#define g12v_bar ((f_dxdu(u(i,j),v(i,j))*f_dxdvdv(u(i,j),v(i,j)) + f_dydu(u(i,j),v(i,j))*f_dydvdv(u(i,j),v(i,j)) + f_dzdu(u(i,j),v(i,j))*f_dzdvdv(u(i,j),v(i,j))) + (f_dxdv(u(i,j),v(i,j))*f_dxdudv(u(i,j),v(i,j)) + f_dydv(u(i,j),v(i,j))*f_dydudv(u(i,j),v(i,j)) + f_dzdv(u(i,j),v(i,j))*f_dzdudv(u(i,j),v(i,j))))


/* page  4 */
#define J_bar (sqrt(g11_bar*g22_bar - g12_bar*g12_bar))

/* page 8 */
#define J_u_bar (1./(2.*J_bar)*( g11_bar*g22u_bar + g22_bar*g11u_bar - 2.*g12_bar*g12u_bar ))
#define J_v_bar (1./(2.*J_bar)*( g11_bar*g22v_bar + g22_bar*g11v_bar - 2.*g12_bar*g12v_bar ))

/* page 7 */
#define dudxi ((u(i+1,j)-u(i-1,j))/2.)
#define dvdxi ((v(i+1,j)-v(i-1,j))/2.)
#define dudeta ((u(i,j+1)-u(i,j-1))/2.)
#define dvdeta ((v(i,j+1)-v(i,j-1))/2.)

/* page 6 */
/* #define g11 (g11_bar*dudxi*dudxi + 2.*g12_bar*dudxi*dvdxi + g22_bar*dvdxi*dvdxi) */
/* #define g22 (g11_bar*dudeta*dudeta + 2.*g12_bar*dudeta*dvdeta + g22_bar*dvdeta*dvdeta) */
/* #define g12 (g11_bar*dudxi*dudeta + g12_bar*(dudxi*dvdeta + dudeta*dvdxi) + g22_bar*dvdxi*dvdeta) */
#define g11 (dudxi*dudxi + dvdxi*dvdxi)
#define g22 (dudeta*dudeta + dvdeta*dvdeta)
#define g12 (dudxi*dudeta + dvdxi*dvdeta)

/* page 4 */
#define J (dudxi*dvdeta - dudeta*dvdxi)


#define P 0.
#define Q (-100./(50.*50.)*exp(-0.5*(j-1)))
//#define Q 0.
//#define P g_array_index (p, gdouble, i + N*(j))
//#define Q g_array_index (q, gdouble, i + N*(j))

/* page 8 */
#ifdef TEST_2D
#define Delta2u 0.
#define Delta2v 0.
#else
#define Delta2u (1./J_bar*( J_bar*(g22u_bar - g12v_bar) - (g22_bar*J_u_bar - g12_bar*J_v_bar) ))
#define Delta2v (1./J_bar*( J_bar*(g11v_bar - g12u_bar) - (g11_bar*J_v_bar - g12_bar*J_u_bar) ))
#endif

static LinearProblem * build_x_elliptic_linear_system (LinearProblem * lpx,
						       LinearProblem * lpy,
						       GArray * p,
						       GArray * q,
						       gint N,
						       gint M)
{
  LinearProblem * lpx_new = linear_problem_new ();
  gint i, j;

  linear_problem_init_size (lpx_new, N, M);

  /* Boundary points should remain unchanged */
  for ( i = 0; i < N; i++ ) {
    j = 0.;
    linear_problem_add_stencil (lpx_new, i+N*j, i+N*j, 1.);
    g_array_index (lpx_new->rhs, gdouble, i+N*j) =
      g_array_index (lpx->lhs, gdouble, i+N*j);
    j = M-1.;
    linear_problem_add_stencil (lpx_new, i+N*j, i+N*j, 1.);
    g_array_index (lpx_new->rhs, gdouble, i+N*j) =
      g_array_index (lpx->lhs, gdouble, i+N*j);
  }
  
  /* Boundary points should remain unchanged */
  for ( j = 1; j < M-1; j++ ) {
    /* i = 0.; */
    /* linear_problem_add_stencil (lpx_new, i+N*j, i+N*j, 1.); */
    /* g_array_index (lpx_new->rhs, gdouble, i+N*j) = */
    /*   g_array_index (lpx->lhs, gdouble, i+N*j); */
    i = N-1.;
    linear_problem_add_stencil (lpx_new, i+N*j, i+N*j, 1.);
    g_array_index (lpx_new->rhs, gdouble, i+N*j) =
      g_array_index (lpx->lhs, gdouble, 0+N*j);
  }


  for ( j = 1; j < M-1; j++ ) {
    i = 0.;
    linear_problem_add_stencil (lpx_new, i+N*j, (i+1) + N*(j+1), -g12/2.);
    linear_problem_add_stencil (lpx_new, i+N*j, (i+1) + N*(j), g22+g22*J*P/2.);
    linear_problem_add_stencil (lpx_new, i+N*j, (i+1) + N*(j-1), g12/2.);
    
    linear_problem_add_stencil (lpx_new, i+N*j, (i) + N*(j+1), g11+g11*J*Q/2.);
    linear_problem_add_stencil (lpx_new, i+N*j, (i) + N*(j), -2.*g22-2.*g11);
    linear_problem_add_stencil (lpx_new, i+N*j, (i) + N*(j-1), g11-g22*J*Q/2.);
      
    linear_problem_add_stencil (lpx_new, i+N*j, (N-2) + N*(j+1), g12/2.);
    linear_problem_add_stencil (lpx_new, i+N*j, (N-2) + N*(j), g22-g11*J*P/2.);
    linear_problem_add_stencil (lpx_new, i+N*j, (N-2) + N*(j-1), -g12/2.);
  }
  
  /* Inner domain */
  for ( i = 1; i < N-1; i++ )
    for ( j = 1; j < M-1; j++ ) {
      linear_problem_add_stencil (lpx_new, i+N*j, (i+1) + N*(j+1), -g12/2.);
      linear_problem_add_stencil (lpx_new, i+N*j, (i+1) + N*(j), g22+g22*J*P/2.);
      linear_problem_add_stencil (lpx_new, i+N*j, (i+1) + N*(j-1), g12/2.);
      
      linear_problem_add_stencil (lpx_new, i+N*j, (i) + N*(j+1), g11+g11*J*Q/2.);
      linear_problem_add_stencil (lpx_new, i+N*j, (i) + N*(j), -2.*g22-2.*g11);
      linear_problem_add_stencil (lpx_new, i+N*j, (i) + N*(j-1), g11-g22*J*Q/2.);
      
      linear_problem_add_stencil (lpx_new, i+N*j, (i-1) + N*(j+1), g12/2.);
      linear_problem_add_stencil (lpx_new, i+N*j, (i-1) + N*(j), g22-g11*J*P/2.);
      linear_problem_add_stencil (lpx_new, i+N*j, (i-1) + N*(j-1), -g12/2.);
    }

  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++)
      g_array_index (lpx_new->rhs, gdouble, i+j*N) = 0.;

  /* Left hand side */
  for ( i = 0; i < N; i++)
    for ( j = 0; j < M; j++)
      g_array_index (lpx_new->lhs, gdouble, i+j*N) =
  	g_array_index (lpx->lhs, gdouble, i+j*N);

  /* Ghost cells */
  for ( i = 0; i < N; i++) {
    g_array_index (lpx_new->gxi0, gdouble, i)
      = g_array_index (lpx->gxi0, gdouble, i);
    g_array_index (lpx_new->gxiN, gdouble, i)
      = g_array_index (lpx->gxiN, gdouble, i);
  }
  for ( i = 0; i < M; i++) {
    g_array_index (lpx_new->geta0, gdouble, i)
      = g_array_index (lpx->geta0, gdouble, i);
    g_array_index (lpx_new->getaN, gdouble, i)
      = g_array_index (lpx->getaN, gdouble, i);
  }
  return lpx_new;
}

static LinearProblem * build_y_elliptic_linear_system (LinearProblem * lpx,
						       LinearProblem * lpy,
						       GArray * p,
						       GArray * q,
						       gint N,
						       gint M)
{
  LinearProblem * lpy_new = linear_problem_new ();
  gint i, j;

  linear_problem_init_size (lpy_new, N, M);

  /* Boundary points should remain unchanged */
  for ( i = 0; i < N; i++ ) {
    j = 0.;
    linear_problem_add_stencil (lpy_new, i+N*j, i+N*j, 1.);
    g_array_index (lpy_new->rhs, gdouble, i+N*j) =
      g_array_index (lpy->lhs, gdouble, i+N*j);
    j = M-1.;
    linear_problem_add_stencil (lpy_new, i+N*j, i+N*j, 1.);
    g_array_index (lpy_new->rhs, gdouble, i+N*j) =
      g_array_index (lpy->lhs, gdouble, i+N*j);
  }
  
  /* Boundary points should remain unchanged */
  for ( j = 1; j < M-1; j++ ) {
    /* i = 0.; */
    /* linear_problem_add_stencil (lpy_new, i+N*j, i+N*j, 1.); */
    /* g_array_index (lpy_new->rhs, gdouble, i+N*j) = */
    /*   g_array_index (lpy->lhs, gdouble, i+N*j); */
    i = N-1.;
    linear_problem_add_stencil (lpy_new, i+N*j, i+N*j, 1.);
    g_array_index (lpy_new->rhs, gdouble, i+N*j) =
      g_array_index (lpy->lhs, gdouble, 0+N*j);
  }


  for ( j = 1; j < M-1; j++ ) {
    i = 0.;
    linear_problem_add_stencil (lpy_new, i+N*j, (i+1) + N*(j+1), -g12/2.);
    linear_problem_add_stencil (lpy_new, i+N*j, (i+1) + N*(j), g22+g22*J*P/2.);
    linear_problem_add_stencil (lpy_new, i+N*j, (i+1) + N*(j-1), g12/2.);
    
    linear_problem_add_stencil (lpy_new, i+N*j, (i) + N*(j+1), g11+g11*J*Q/2.);
    linear_problem_add_stencil (lpy_new, i+N*j, (i) + N*(j), -2.*g22-2.*g11);
    linear_problem_add_stencil (lpy_new, i+N*j, (i) + N*(j-1), g11-g22*J*Q/2.);
      
    linear_problem_add_stencil (lpy_new, i+N*j, (N-2) + N*(j+1), g12/2.);
    linear_problem_add_stencil (lpy_new, i+N*j, (N-2) + N*(j), g22-g11*J*P/2.);
    linear_problem_add_stencil (lpy_new, i+N*j, (N-2) + N*(j-1), -g12/2.);
  }
  
  /* Inner domain */
  for ( i = 1; i < N-1; i++ )
    for ( j = 1; j < M-1; j++ ) {
      linear_problem_add_stencil (lpy_new, i+N*j, (i+1) + N*(j+1), -g12/2.);
      linear_problem_add_stencil (lpy_new, i+N*j, (i+1) + N*(j), g22+g22*J*P/2.);
      linear_problem_add_stencil (lpy_new, i+N*j, (i+1) + N*(j-1), g12/2.);
      
      linear_problem_add_stencil (lpy_new, i+N*j, (i) + N*(j+1), g11+g11*J*Q/2.);
      linear_problem_add_stencil (lpy_new, i+N*j, (i) + N*(j), -2.*g22-2.*g11);
      linear_problem_add_stencil (lpy_new, i+N*j, (i) + N*(j-1), g11-g22*J*Q/2.);
      
      linear_problem_add_stencil (lpy_new, i+N*j, (i-1) + N*(j+1), g12/2.);
      linear_problem_add_stencil (lpy_new, i+N*j, (i-1) + N*(j), g22-g11*J*P/2.);
      linear_problem_add_stencil (lpy_new, i+N*j, (i-1) + N*(j-1), -g12/2.);
    }

  for ( i = 1; i < N-1; i++)
    for ( j = 1; j < M-1; j++)
      g_array_index (lpy_new->rhs, gdouble, i+j*N) = 0.;

  /* Left hand side */
  for ( i = 0; i < N; i++)
    for ( j = 0; j < M; j++)
      g_array_index (lpy_new->lhs, gdouble, i+j*N) =
  	g_array_index (lpy->lhs, gdouble, i+j*N);

  /* Ghost cells */
  for ( i = 0; i < N; i++) {
    g_array_index (lpy_new->gxi0, gdouble, i)
      = g_array_index (lpy->gxi0, gdouble, i);
    g_array_index (lpy_new->gxiN, gdouble, i)
      = g_array_index (lpy->gxiN, gdouble, i);
  }
  for ( i = 0; i < M; i++) {
    g_array_index (lpy_new->geta0, gdouble, i)
      = g_array_index (lpy->geta0, gdouble, i);
    g_array_index (lpy_new->getaN, gdouble, i)
      = g_array_index (lpy->getaN, gdouble, i);
  }
  return lpy_new;
}

void elliptic_smoothing_periodic (Spline2D * sp)
{
  gint i, j, k = sp->k, a, b;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble A[size][size];
  size_t ustart, vstart, uend, vend;
  
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A[i][j] = 0.;
    }
  }

  gsl_vector * rhs_x = gsl_vector_alloc (size);
  gsl_vector * rhs_y = gsl_vector_alloc (size);
  gsl_vector_set_zero (rhs_x);
  gsl_vector_set_zero (rhs_y);

  GrevillePoints * gr = sp->gr;
  for ( i = 0; i < sp->NU; i++ ) {
    gdouble u = g_array_index (gr->ui, gdouble, i);
    gsl_matrix * Bu = gsl_matrix_alloc (k, 3);
    gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,u), 2, Bu, &ustart, &uend, sp->w_u, sp->wd_u);

    for ( j = 1; j < NV-1; j++ ) {
      gdouble v = g_array_index (gr->vj, gdouble, j);
      gsl_matrix * Bv = gsl_matrix_alloc (k, 3);
      gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
      gint indexi = i + j*NU;

       gdouble dxdeta = spline2d_derivative_eval_greville_point (sp, gr, i, j, 0, 1, 0);
       gdouble dydeta = spline2d_derivative_eval_greville_point (sp, gr, i, j, 0, 1, 1);
       gdouble dxdxi = spline2d_derivative_eval_greville_point (sp, gr, i, j, 1, 0, 0);
       gdouble dydxi = spline2d_derivative_eval_greville_point (sp, gr, i, j, 1, 0, 1);

       gdouble alpha = dxdeta*dxdeta + dydeta*dydeta;
       gdouble beta = dxdxi*dxdxi + dydxi*dydxi;
       gdouble gamma = dxdxi*dxdeta + dydxi*dydeta;
       gdouble JJ = dxdxi*dydeta - dxdeta*dydxi;
       gdouble PP = 0.;
       gdouble QQ = -100./pow(10., 2.)*exp (-0.5*(i-2.));

       for ( a = ustart; a < ustart+k; a++ ) {
	 gint aa = a-ustart;
	 gint atmp = (a-k+1)%NU;
	 //for ( b = MAX(1, vstart); b < MIN(vstart+k, NV-1); b++) {
	 for ( b = vstart; b < vstart+k; b++) {
	   gint bb = b-vstart;

	   A[atmp+(b/* -1 */)*NU][indexi] += 
	     alpha*gsl_matrix_get (Bu, aa, 2)*gsl_matrix_get (Bv, bb, 0)
	     - 2.*gamma*gsl_matrix_get (Bu, aa, 1)*gsl_matrix_get (Bv, bb, 1)
	     + beta*gsl_matrix_get (Bu, aa, 0)*gsl_matrix_get (Bv, bb, 2)
	     + JJ*(PP*gsl_matrix_get (Bu, aa, 1)*gsl_matrix_get (Bv, bb, 0)+
		   QQ*gsl_matrix_get (Bu, aa, 0)*gsl_matrix_get (Bv, bb, 1));

	 }
       }
    }

    j = 0;
    gint indexi = i + j*NU;
    gdouble v = g_array_index (gr->vj, gdouble, j);
    gsl_matrix * Bv = gsl_matrix_alloc (k, 3);
    gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
    for ( a = ustart; a < ustart+k; a++ ) {
      gint aa = a-ustart;
      gint atmp = (a-k+1)%NU;
      //for ( b = MAX(1, vstart); b < MIN(vstart+k, NV-1); b++) {
      for ( b = vstart; b < vstart+k; b++) {
	gint bb = b-vstart;
	A[atmp+(b/* -1 */)*NU][indexi] += gsl_matrix_get (Bu, aa, 0)*gsl_matrix_get (Bv, bb, 0);
      }
    }
    gsl_vector_set (rhs_x, indexi, spline2d_eval_greville_point (sp, gr, i, j, 0));
    gsl_vector_set (rhs_y, indexi, spline2d_eval_greville_point (sp, gr, i, j, 1));

    j = NV-1;
    indexi = i + j*NU;
    v = g_array_index (gr->vj, gdouble, j);
    Bv = gsl_matrix_alloc (k, 3);
    gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
    for ( a = ustart; a < ustart+k; a++ ) {
      gint aa = a-ustart;
      gint atmp = (a-k+1)%NU;
      //for ( b = MAX(1, vstart); b < MIN(vstart+k, NV-1); b++) {
      for ( b = vstart; b < vstart+k; b++) {
	gint bb = b-vstart;
	A[atmp+(b/* -1 */)*NU][indexi] += gsl_matrix_get (Bu, aa, 0)*gsl_matrix_get (Bv, bb, 0);
      }
    }
    gsl_vector_set (rhs_x, indexi, spline2d_eval_greville_point (sp, gr, i, j, 0));
    gsl_vector_set (rhs_y, indexi, spline2d_eval_greville_point (sp, gr, i, j, 1));
  }
  fprintf (stderr,"Toto %i %i %i \n", NU, NV, size);
  /* Storage in Compressed Column Storage (CCS) format for superlu */
  FILE * fp = fopen ("elliptic.tmp","w");
  fprintf (stderr,"Toto %i \n", size);
  CCSProblem * css = ccs_problem_new ();
  gint count = 0;
  size = NU*NV;
  for ( i = 0; i < size; i++) {
    g_array_append_val (css->index, count);
    for ( j = 0; j < size; j++) {
      fprintf (fp, "%i %i %e\n",i,j,A[i][j]);

      if (A[i][j] != 0.) {
  	g_array_append_val (css->matrix, A[i][j]);
  	g_array_append_val (css->column, j);
  	count++;
      }
      
    }
    fprintf (fp, "\n");
  }
  g_array_append_val (css->index, count);

  fclose (fp);

  ccs_problem_lu_solve (css, rhs_x);
  ccs_problem_lu_solve (css, rhs_y);

  for ( i = 0; i < NU+k-1; i++ ) {
    gint indexi = (i + k - 1)%(NU+k-1);
    for ( j = 0; j < NV; j++) {
      coeff_assign (sp, indexi, j, 0, gsl_vector_get (rhs_x, i%NU+j*NU));
      coeff_assign (sp, indexi, j, 1, gsl_vector_get (rhs_y, i%NU+j*NU));
    }
  }

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  ccs_problem_destroy (css);

}


/**
 * Returns the value of x in (i,j). If i or j equals -1 or N, then
 * the function return the ghost cell value as defined in
 * chapter 6 of the Handbook of grid generation
 */
static gdouble bx (LinearProblem * lpx, LinearProblem * lpy,
		   gint N, gint M, gint i, gint j)
{
  if ( i == -1 && (j == -1 || j == M))
    g_assert_not_reached ();
  if ( i == N && (j == -1 || j == M))
    g_assert_not_reached ();

  /* u(-1,j) = 2*u(0,j)-u(1,j)*/
  if (i == -1)
    return g_array_index (lpx->geta0, gdouble, j);

  /* u(N,j) = 2*u(N-1,j)-u(N-2,j) */
  if (i == N)
    return g_array_index (lpx->getaN, gdouble, j);

  /* u(i,-1) = u(i,0) + g12_bar/g11_bar*(v(i,1)-v(i,0)) */
  if (j == -1)
    return g_array_index (lpx->gxi0, gdouble, i);

  /* u (i,N) = u(i,N-1) - g12_bar/g11_bar*(v(i,N-1) - v(i,N-2) ) */
  if (j == M)
    return g_array_index (lpx->gxiN, gdouble, i);

  return u(i,j);
}

/**
 * Returns the value of y in (i,j). If i or j equals -1 or N, then
 * the function return the ghost cell value as defined in
 * chapter 6 of the Handbook of grid generation
 */
static gdouble by (LinearProblem * lpx, LinearProblem * lpy,
		   gint N, gint M, gint i, gint j)
{
  if ( i == -1 && (j == -1 || j == M))
    g_assert_not_reached ();
  if ( i == N && (j == -1 || j == M))
    g_assert_not_reached ();
  
  /* v(-1,j) = v(0,j) + g12_bar/g22_bar*(u(1,j)-u(0,j))*/
  if (i == -1)
    return g_array_index (lpy->geta0, gdouble, j)/* 1./(N-1.) */;

  /* v(N,j) = v(n-1,j) - g12_bar/g22_bar (u(N-1,j) - u(N-2,j)) */
  if (i == N)
    return g_array_index (lpy->getaN, gdouble, j);

  /* v(i,-1) = 2*v(i,0)-v(i,1) */
  if (j == -1)
    return g_array_index (lpy->gxi0, gdouble, i);

  /* v(i,N) = 2*v(i,N-1) - v(i,N-2) */
  if (j == M)
    return g_array_index (lpy->gxiN, gdouble, i);

  return v(i,j);
}

/* Macros for values of x and y (including ghost cells */
#define bvx(a,b) (bx (lpx,lpy,N,M,a,b))
#define bvy(a,b) (by (lpx,lpy,N,M,a,b))

static void print_ghost_cells (LinearProblem * lpx, LinearProblem * lpy, gint N, gint M)
{
  FILE * fg = fopen("ghost.tmp","w");
  gint i, j;

  for ( i = 0; i < N; i++ ) {
    j = -1;
    fprintf(fg, "%f %f %f\n", /* i/(N-1.) */bvx(i,j), /* j/(N-1.) */bvy(i,j), f_z (bvx(i,j),bvy(i,j), 0, NULL));
  }

  for ( i = 0; i < N; i++ ) {
    j = M;
    fprintf(fg, "%f %f %f\n", /* i/(N-1.) */bvx(i,j), /* j/(N-1.) */bvy(i,j), f_z (bvx(i,j),bvy(i,j), 0, NULL));
  }

  for ( j = 0; j < M; j++ ) {
    i = -1;
    fprintf(fg, "%f %f %f\n", /* i/(N-1.) */bvx(i,j), /* j/(N-1.) */bvy(i,j), f_z (bvx(i,j),bvy(i,j), 0, NULL));
  }

  for ( j = 0; j < M; j++ ) {
    i = N;
    fprintf(fg, "%f %f %f\n", /* i/(N-1.) */bvx(i,j), /* j/(N-1.) */bvy(i,j), f_z (bvx(i,j),bvy(i,j), 0, NULL));
  }

  fclose (fg);
}

/**
 * Ghost cells can be set in an attempt to impose boundary othogonality
 **/
static void set_ghost_cells (LinearProblem * lpx, LinearProblem * lpy,
			     gint N, gint M)
{
  g_assert (lpx != NULL);
  g_assert (lpy != NULL);
  gint i, j;

  /* xi boundaries */
  for ( i = 0; i < N; i++) {
    /* j = -1 */;
    j = 0;
    /* u(i,-1) = u(i,0) + g12_bar/g11_bar*(v(i,1)-v(i,0)) */
    g_array_index (lpx->gxi0, gdouble, i)
      = 2.*u(i,0) - u(i,1);
    /* v(i,-1) = 2*v(i,0)-v(i,1) */
    g_array_index (lpy->gxi0, gdouble, i)
      = 2*v(i,0)-v(i,1);
    /* j = N */
    j = M-1;
    /* u (i,N) = u(i,N-1) - g12_bar/g11_bar*(v(i,N-1) - v(i,N-2) ) */
    g_array_index (lpx->gxiN, gdouble, i)
      = 2*u(i,M-1)  - u(i,M-2);
    /* v(i,N) = 2*v(i,N-1) - v(i,N-2) */
    g_array_index (lpy->gxiN, gdouble, i)
      = 2*v(i,M-1) - v(i,M-2);
  }

  /* eta boundaries */
  for ( j = 0; j < M; j++) {
    /* i = -1 */
    i = 0;
    /* u(-1,j) = 2*u(0,j)-u(1,j)*/
    g_array_index (lpx->geta0, gdouble, j)
      = 2*u(0,j)-u(1,j);
    /* v(-1,j) = v(0,j) + g12_bar/g22_bar*(u(1,j)-u(0,j))*/
    g_array_index (lpy->geta0, gdouble, j)
      = 2*v(0,j)-v(1,j);
    /* i = N */
    i = N-1;
    /* u(N,j) = 2*u(N-1,j)-u(N-2,j) */
    g_array_index (lpx->getaN, gdouble, j)
      = 2*u(N-1,j)-u(N-2,j);
    /* v(N,j) = v(n-1,j) - g12_bar/g22_bar (u(N-1,j) - u(N-2,j)) */
    g_array_index (lpy->getaN, gdouble, j)
      = 2*v(N-1,j)-v(N-2,j);
  }
}

/* Macros for second derivatives */
#define uxx(a,b) (bvx(a+1,b)-2.*bvx(a,b)+bvx(a-1,b))
#define vxx(a,b) (bvy(a+1,b)-2.*bvy(a,b)+bvy(a-1,b))
#define uee(a,b) (bvx(a,b+1)-2.*bvx(a,b)+bvx(a,b-1))
#define vee(a,b) (bvy(a,b+1)-2.*bvy(a,b)+bvy(a,b-1))

#define uxe(a,b) (0.25*(bvx(a+1,b+1)+bvx(a-1,b-1)-bvx(a+1,b-1)-bvx(a-1,b+1)))
#define vxe(a,b) (0.25*(bvy(a+1,b+1)+bvy(a-1,b-1)-bvy(a+1,b-1)-bvy(a-1,b+1)))

/* Macros for first derivatives */
#define ux(a,b) (u_x (lpx,lpy,N,M,a,b))
#define vx(a,b) (v_x (lpx,lpy,N,M,a,b))
#define ue(a,b) (u_e (lpx,lpy,N,M,a,b))
#define ve(a,b) (v_e (lpx,lpy,N,M,a,b))

#define G11 (g11_bar*us(i,j)*us(i,j) + 2.*g12_bar*us(i,j)*vs(i,j) + g22_bar*vs(i,j)*vs(i,j))
#define G22 (g11_bar*us(i,j)*ut(i,j) + 2.*g12_bar*ut(i,j)*vt(i,j) + g22_bar*vt(i,j)*vt(i,j))
#define G12 (g11_bar*us(i,j)*ut(i,j) + g12_bar*(us(i,j)*vt(i,j) + ut(i,j)*vs(i,j)) + g22_bar*vs(i,j)*vt(i,j))

/**
 * Returns the value of the control function P at a given point
 * of the boundary.
 */
static double pb (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		  gint N, gint M, gint i, gint j)
{
  /* gdouble J2 = (us(i,j)*vt(i,j) - ut(i,j)*vs(i,j)); */

#ifdef TEST_2D
  return - (g11_bar*uee(i,j)*us(i,j) + g22_bar*vee(i,j)*vs(i,j))/G22
    - (g11_bar*uxx(i,j)*us(i,j) + g22_bar*vxx(i,j)*vs(i,j))/G11;
#else
  return J2*J2*( g11_bar*Delta2u*us(i,j) + g22_bar*Delta2v*vs(i,j)
		 + g12_bar*(Delta2u*vs(i,j) + Delta2v*us(i,j)))/(G11*G22)
    - (g11_bar*uee(i,j)*us(i,j) + g22_bar*vee(i,j)*vs(i,j)
       + g12_bar*(uee(i,j)*vs(i,j) + vee(i,j)*us(i,j)))/G22
    - (g11_bar*uxx(i,j)*us(i,j) + g22_bar*vxx(i,j)*vs(i,j)
       + g12_bar*(uxx(i,j)*vs(i,j) + vxx(i,j)*us(i,j)))/G11;
#endif
}

/**
 * Returns the value of the control function Q at a given point
 * of the boundary.
 */
static double qb (Boundaries * b, LinearProblem * lpx, LinearProblem * lpy,
		  gint N, gint M, gint i, gint j)
{
  /* gdouble J2 = (us(i,j)*vt(i,j) - ut(i,j)*vs(i,j)); */

#ifdef TEST_2D
  return - (g11_bar*uee(i,j)*ut(i,j) + g22_bar*vee(i,j)*vt(i,j))/G22
    - (g11_bar*uxx(i,j)*ut(i,j) + g22_bar*vxx(i,j)*vt(i,j))/G11;
#else
  return J2*J2*( g11_bar*Delta2u*ut(i,j) + g22_bar*Delta2v*vt(i,j)
		 + g12_bar*(Delta2u*vt(i,j) + Delta2v*ut(i,j)))/(G11*G22)
    - (g11_bar*uee(i,j)*ut(i,j) + g22_bar*vee(i,j)*vt(i,j)
       + g12_bar*(uee(i,j)*vt(i,j) + vee(i,j)*ut(i,j)))/G22
    - (g11_bar*uxx(i,j)*ut(i,j) + g22_bar*vxx(i,j)*vt(i,j)
       + g12_bar*(uxx(i,j)*vt(i,j) + vxx(i,j)*ut(i,j)))/G11;
#endif
}

#define p(a,b) g_array_index (p, gdouble, a + (b)*N)
#define q(a,b) g_array_index (q, gdouble, a + (b)*N)

static GArray * smooth_p (GArray * p, gint N, gint M)
{
  gint i, j;
  GArray * sp = g_array_sized_new (FALSE, TRUE, sizeof(gdouble), N*M);

  for ( i = 0; i < N; i++) {
    for ( j = 0; j < M; j++) {
      if (j == 0.) {
	gdouble p0 = (p(i,j+1)+p(i,j))/2.;
	g_array_append_val (sp, p0);
      }
      else if (j == M-1) {
	gdouble p0 = (p(i,j)+p(i,j-1))/2.;
	g_array_append_val (sp, p0);
      }
      
      /* if (i == 0 || i == N-1 || j == 0 || j == M-1) { */
      
      /* 	g_array_append_val (sp, p(i,j)); */
      /* } */
      else {
	gdouble p0 = (p(i,j+1)+p(i,j-1))/2.;
	g_array_append_val (sp, p0);
      }
    }
  }

  for ( i = 0; i < N; i++)
    for ( j = 0; j < M; j++)
      p(i,j) = g_array_index (sp, gdouble, i + (j)*N);

  g_array_free (sp, FALSE);

  return p;
}

static GArray * smooth_q (GArray * q, gint N, gint M)
{
  gint i, j;
  GArray * sq = g_array_sized_new (FALSE, TRUE, sizeof(gdouble), N*M);

  for ( i = 0; i < N; i++) {
    for ( j = 0; j < M; j++) {
      if (i == 0.) {
	gdouble q0 = (q(i+1,j)+q(i,j))/2.;
	g_array_append_val (sq, q0);
      }
      else if (i == N-1) {
	gdouble q0 = (q(i,j)+q(i-1,j))/2.;
	g_array_append_val (sq, q0);
      }


      /* if (i == 0 || i == N-1 || j == 0 || j == M-1) { */
      /* 	g_array_append_val (sq, q(i,j)); */
      /* } */
      else {
	gdouble q0 = (q(i+1,j)+q(i-1,j))/2.;
	g_array_append_val (sq, q0);
      }
    }
  }

  for ( i = 0; i < N; i++)
    for ( j = 0; j < M; j++)
      q(i,j) = g_array_index (sq, gdouble, i + (j)*N);

  g_array_free (sq, FALSE);

  return q;
}

/**
 * Computes the control functions P and Q as defined in Chapter 6 of the 
 * handbook of grid generation on the boundaries and then interpolate them
 * to the rest of the domain using linear transfinite interpolation.
 */
static void update_p (Boundaries * b, LinearProblem * lpx,
		      LinearProblem * lpy,
		      GArray * p, GArray * q, gint N, gint M)
{
  gint i, j;

  /* xi boundaries */
  for (i = 0; i < N; i++) {
    j = 0;
    g_array_index (p, gdouble, i + j*M) = pb (b, lpx, lpy, N, M, i, j);
    g_array_index (q, gdouble, i + j*M) = qb (b, lpx, lpy, N, M, i, j);
    j = M-1;
    g_array_index (p, gdouble, i + j*M) = pb (b, lpx, lpy, N, M, i, j);
    g_array_index (q, gdouble, i + j*M) = qb (b, lpx, lpy, N, M, i, j);
  }

  /* eta boundaries */
  for (j = 1; j < M-1; j++) {
    i = 0;
    g_array_index (p, gdouble, i + j*M) = pb (b, lpx, lpy, N, M, i, j);
    g_array_index (q, gdouble, i + j*M) = qb (b, lpx, lpy, N, M, i, j);
    i = N-1;
    g_array_index (p, gdouble, i + j*M) = pb (b, lpx, lpy, N, M, i, j);
    g_array_index (q, gdouble, i + j*M) = qb (b, lpx, lpy, N, M, i, j);
  }

  /* Transfinite linear interpolation to the rest of the domain */
  gdouble p00 = g_array_index (p, gdouble, 0);
  gdouble pN0 = g_array_index (p, gdouble, N-1);
  gdouble pNN = g_array_index (p, gdouble, N-1 + (M-1)*N);
  gdouble p0N = g_array_index (p, gdouble, (M-1)*N);

  gdouble q00 = g_array_index (q, gdouble, 0);
  gdouble qN0 = g_array_index (q, gdouble, N-1);
  gdouble qNN = g_array_index (q, gdouble, N-1 + (M-1)*N);
  gdouble q0N = g_array_index (q, gdouble, (M-1)*N);

  for (i = 1; i < N-1; i++) {
    gdouble xi = i/(N-1.);
    gdouble pxi0 = g_array_index (p, gdouble, i);
    gdouble pxiN = g_array_index (p, gdouble, i + (M-1)*N);
    gdouble qxi0 = g_array_index (q, gdouble, i);
    gdouble qxiN = g_array_index (q, gdouble, i + (M-1)*N);

    for (j = 1; j < M-1; j++) {
      gdouble eta = j/(M-1.);
      gdouble peta0 = g_array_index (p, gdouble, 0 + j*N);
      gdouble petaN = g_array_index (p, gdouble, (N-1) + j*N);
      gdouble qeta0 = g_array_index (q, gdouble, 0 + j*N);
      gdouble qetaN = g_array_index (q, gdouble, (N-1) + j*N);

      g_array_index (p, gdouble, i + j*N) = (1.-eta)*pxi0 + eta*pxiN
      	+ (1.-xi)*peta0 + xi*petaN- p00*(1.-xi)*(1.-eta) - pN0*xi*(1-eta)
      	- p0N*(1-xi)*eta - pNN*xi*eta;
      g_array_index (q, gdouble, i + j*N) = (1.-eta)*qxi0 + eta*qxiN
      	+ (1.-xi)*qeta0 + xi*qetaN - q00*(1.-xi)*(1.-eta) - qN0*xi*(1-eta)
      	- q0N*(1-xi)*eta - qNN*xi*eta;
    }
  }
  
  /* /\* 5 iterations of smoothing *\/ */
  /* for ( i = 0; i < 4; i++ ){ */
  /*   p = smooth_p (p, N, M); */
  /*   q = smooth_q (q, N, M); */
  /* } */

  FILE * fp = fopen("p.tmp","w");
  for (i=0;i<N;i++)
    for (j=0;j<M;j++) {
      fprintf(fp,"%i %i %g %g\n",i,j,g_array_index (p, gdouble, i + j*N),g_array_index (q, gdouble, i + j*N));
    }
  fclose(fp);
}

static void store_grid (Surface * s, LinearProblem * lpx, LinearProblem * lpy, gint N, gint M)
{
  Patch * patch = patch_new ();
  gint i, j, k;
  
  for ( j = 1; j < M; j++) {
    patch_add_row (patch);
    for ( i = 1; i < N; i++) {
      Panel * p = panel_new ();
      
      p->p[0].x = u(i-1, j-1);
      p->p[0].y = v(i-1, j-1);
      p->p[0].z = s->hz (p->p[0].x, p->p[0].y, 0., NULL);
      p->p[1].x = u(i, j-1);
      p->p[1].y = v(i, j-1);
      p->p[1].z = s->hz (p->p[1].x, p->p[1].y, 0., NULL);
      p->p[2].x = u(i, j);
      p->p[2].y = v(i, j);
      p->p[2].z = s->hz (p->p[2].x, p->p[2].y, 0., NULL);
      p->p[3].x = u(i-1, j);
      p->p[3].y = v(i-1, j);
      p->p[3].z = s->hz (p->p[3].x, p->p[3].y, 0., NULL);

      p->var[0] = p->var[1] = p->var[2] = 0.;
      for (k = 0; k < 4; k++) {
  	p->var[0] += p->p[k].x;
  	p->var[1] += p->p[k].y;
  	p->var[2] += p->p[k].z;
      }
      p->var[0] /= 4.; p->var[1] /= 4.; p->var[2] /= 4.;

      if (j == 1 || j == M-1 || i == 1 || i == N-1)
  	p->border = TRUE;
      
      patch_add_panel (patch, *p);
    }
  }

  s->patches = g_slist_append (s->patches, patch);
}

static gdouble eval_grid_x (Spline2D * grid, gdouble u, gdouble v, gpointer data)
{
  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval_x (sp, u, v, 0);
}

static gdouble eval_grid_y (Spline2D * grid, gdouble u, gdouble v, gpointer data)
{
  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval_x (sp, u, v, 1);
}

/* static gdouble eval_grid_x_gauss (Spline2D * grid, */
/* 				  gdouble u, gdouble v, */
/* 				  gpointer data) */
/* { */
/*   Spline2D * sp = (Spline2D *) data; */

/*   return spline2d_eval (sp, u, v, 0); */
/* } */

static gdouble eval_grid_x_gauss (SPPanel * spp,
				  gint m, gint n,
				  gpointer data)
{
  GaussPoints * gp = spp->outer;
  gdouble u = g_array_index (gp->ui, gdouble, m);
  gdouble v = g_array_index (gp->vj, gdouble, n);
  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval_x (sp, u, v, 0);
}

static gdouble eval_grid_y_gauss (SPPanel * spp,
				  gint m, gint n,
				  gpointer data)
{
  GaussPoints * gp = spp->outer;
  gdouble u = g_array_index (gp->ui, gdouble, m);
  gdouble v = g_array_index (gp->vj, gdouble, n);
  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval_x (sp, u, v, 1);
}

void free_surface_galerkin_fit (Spline2D * sp, Spline2DFunc func, gpointer data, gint var)
{
  gint i, j, m, n, ii, a, b;
  size_t ustart, vstart;
  gint k = sp->k;
  gint size = /* gsl_bspline_ncoeffs (sp->w_u) */(3*k+sp->M-1)*gsl_bspline_ncoeffs (sp->w_v);
  gdouble A0[size][size];
  gdouble RHS[size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
    RHS[i] = 0.;
  }
  
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
    ustart = gp->istart < k-1+sp->M ? gp->istart - (k-1) : gp->istart - sp->M - (k-1);
    vstart = gp->jstart;
    
    for ( m = 0; m < ng; m++) {
      gdouble um = g_array_index (gp->ui, gdouble, m);
      gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
  	gdouble vn = g_array_index (gp->vj, gdouble, n);
  	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
  	gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);

  	gdouble rhsmn = func (sp, um, vn, data);

  	// Loop over the splines whose support is included in the panel
  	for ( i = ustart; i < ustart + sp->k; i++) {
  	  gdouble wmni = wmn*gsl_matrix_get(Bu, i-ustart, 0);
  	  for ( j = vstart; j < vstart + sp->k; j++) {
  	    gdouble wmnij = wmni*gsl_matrix_get(Bv, j-vstart, 0);
  	    gint indexi = i + j*gsl_bspline_ncoeffs (sp->w_u);

  	    for ( a = 0; a < sp->k; a++) {
  	      gdouble wmnija = wmnij*gsl_matrix_get(Bu, a, 0);
  	      for ( b = 0; b < sp->k; b++) {
  		A0[indexi][(ustart+a) + (vstart+b)*gsl_bspline_ncoeffs (sp->w_u)] +=  wmnija*gsl_matrix_get(Bv, b, 0);
  	      }
  	    }
  	    RHS[indexi] += wmnij*rhsmn;
  	  }
  	}
	
      }
    }
  }

  FILE * fp = fopen ("galerkin.tmp","w");
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      
      fprintf (fp, "%i %i %f \n", i, j, A0[i][j]);
    }
    fprintf (fp, "\n");
    //fprintf (fp, "%i %f \n", i, RHS[i]);
  }
  fclose (fp);

  gsl_vector * gsl_rhs = gsl_vector_alloc (size);

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * fit = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (fit->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (fit->matrix, A0[i][j]);
  	g_array_append_val (fit->column, j);
  	count++;
      }
    }
    gsl_vector_set (gsl_rhs, i, RHS[i]);
  }
  g_array_append_val (fit->index, count);


  ccs_problem_lu_solve (fit, gsl_rhs);

  for ( i = 0; i < gsl_bspline_ncoeffs (sp->w_u); i++) {
    for ( j = 0; j < gsl_bspline_ncoeffs (sp->w_v); j++) {
      coeff_assign (sp, i+k-1, j, var, gsl_vector_get (gsl_rhs, i+j*gsl_bspline_ncoeffs (sp->w_u)));
      if ( i < k-1)
	coeff_assign (sp, i+k-1+sp->M, j, var, gsl_vector_get (gsl_rhs, i+j*gsl_bspline_ncoeffs (sp->w_u)));
    }
  }

  gsl_vector_free (gsl_rhs);
  ccs_problem_destroy (fit);
}

static void spline2d_store_grid (Surface * s, LinearProblem * lpx, LinearProblem * lpy,
				 gint N, gint M, gboolean flip)
{
  gint k = 3;
  Spline2D * sp =  periodic_fs_new (N-1, M-k+1, k, 4, 3);
  periodic_fs_init_panels (sp);

  gint i, j, a, b;
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gsl_vector * rhs_x = gsl_vector_alloc (size);
  gsl_vector * rhs_y = gsl_vector_alloc (size);
  gdouble A0[size][size];

  for ( i = 0; i < size; i++)
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;
  
  /* gint k = sp->k; */
  gsl_vector * Bv = gsl_vector_alloc (k);
  for ( i = 0; i < NU; i++) {
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i);
    gint ustart = g_array_index (gr->ustart, gint, i);

    for ( j = 0; j < NV; j++) {
      gdouble v = j/(NV-1.);
      size_t vstart, vend;      
      gsl_bspline_eval_nonzero (MIN(1.-1e-12,v), Bv, &vstart, &vend, sp->w_v);
      gint indexi = i+j*NU;
      
      for ( a = 0; a < k; a++) {
	gint atmp = (ustart+a-k+1)%NU;
	for ( b = 0; b < k; b++) {
	  A0[atmp + (vstart+b)*NU][indexi] = gsl_matrix_get (Bu, a, 0)*gsl_vector_get (Bv, b);
	}
      }

      if (flip) {
  	gsl_vector_set (rhs_x, indexi, u(N-(i+1),j));
  	gsl_vector_set (rhs_y, indexi, v(N-(i+1),j));
      }
      else {
  	gsl_vector_set (rhs_x, indexi,  u(i,j));
  	gsl_vector_set (rhs_y, indexi,  v(i,j));
      }
    }
  } 
  gsl_vector_free (Bv);

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * css = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (css->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (css->matrix, A0[i][j]);
  	g_array_append_val (css->column, j);
  	count++;
      }
    }
  }
  g_array_append_val (css->index, count);

  ccs_problem_lu_solve (css, rhs_x);
  ccs_problem_lu_solve (css, rhs_y);

  for ( i = 0; i < NU+k-1; i++ ) {
    gint indexi = (i + k - 1)%(NU+k-1);
    for ( j = 0; j < NV; j++) {     
      coeff_assign (sp, indexi, j, 0, gsl_vector_get (rhs_x, i%NU+j*NU));
      coeff_assign (sp, indexi, j, 1, gsl_vector_get (rhs_y, i%NU+j*NU));
    }
  }

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  ccs_problem_destroy (css);

  periodic_fs_reinit_panels_physical_quantities (sp);

  s->patches = g_slist_append (s->patches, sp);
}

typedef struct {
  LinearProblem * lpx, * lpy;
  gint M, N;
} GridData;

/* gdouble test_u (Spline2D * sp, gdouble u, gdouble v, gpointer data) */
/* { */
  
/*   gdouble * A = (gdouble *) data; */
  
/*   gint m = u*sp->M; */


/*   return; */
/* } */

/* static void spline2d_store_grid (Surface * s, LinearProblem * lpx, LinearProblem * lpy, */
/* 				 gint N, gint M, gboolean flip) */
/* { */
/*   gint k = 3; */
/*   Spline2D * patch = spline2d_new (N, M, k, 4, 3); */
/*   /\* Spline2D * grid =  periodic_fs_new (N, M-k+1, k, 4, 3); *\/ */
/*   /\* periodic_fs_init_panels (grid); *\/ */

/*   Spline2D * grid =  spline2d_new (N-k+1, M-k+1, k, 4, 3); */
/*   spline2d_init_panels (grid); */


/*   gint i, j, m, n, ii, a, b; */
/*   gint NU = grid->NU; */
/*   gint NV = grid->NV; */
/*   g_assert (NU == N && NV == M); */
/*   gint size = NU*NV; */
/*   gdouble A0[size][size]; */
/*   gsl_vector * rhs_x = gsl_vector_alloc (size); */
/*   gsl_vector * rhs_y = gsl_vector_alloc (size); */
  
/*   fprintf (stderr, "%i %i %i %i \n",N, M, patch->NU, patch->NV);  */

/*   for ( i = 0; i < size; i++) */
/*     for ( j = 0; j < size; j++) */
/*       A0[i][j] = 0.; */

/*   gsl_vector * Bu = gsl_vector_alloc (k); */
/*   gsl_vector * Bv = gsl_vector_alloc (k); */
/*   size_t istart, iend, jstart, jend; */

/*   for ( i = 0; i < NU; i++) { */
/*     gdouble u = i/(NU-1.); */
/*     gsl_bspline_eval_nonzero (u, Bu, &istart, &iend, grid->w_u); */

/*     for ( j = 0; j < NV; j++) { */
/*       gdouble v = j/(NV-1.); */
/*       gsl_bspline_eval_nonzero (v, Bv, &jstart, &jend, grid->w_v); */

/*       /\* gsl_matrix * Bv = g_ptr_array_index (gr->Bv, j+1); *\/ */
/*       /\* gint jstart = g_array_index (gr->vstart, gint, j+1); *\/ */

/*       /\* for ( a = MAX(1, istart); a < MIN(istart + sp->k, NU-1); a++) { *\/ */
/*       /\* 	for ( b = MAX(1, jstart); b < MIN(jstart + sp->k, NV-1); b++) { *\/ */
/*       for ( a = istart; a < istart + grid->k; a++) { */
/*       	for ( b = jstart; b < jstart + grid->k; b++) { */
/* 	  A0[i+j*NU][a+b*NU] = gsl_vector_get (Bu, a-istart)*gsl_vector_get (Bv, b-jstart); */
/*       	} */
/*       } */

/*       /\* gsl_vector_set (rhs, i+j*(NU-2), func (sp, u, v, data) - spline2d_eval_greville_point (sp, gr, i+1, j+1, var)); *\/ */
/*       if (flip) { */
/* 	gsl_vector_set (rhs_x, i+j*NU, u(N-(i+1),j)); */
/* 	gsl_vector_set (rhs_x, i+j*NU, v(N-(i+1),j)); */
/*       } */
/*       else { */
/* 	gsl_vector_set (rhs_x, i+j*NU,  u(i,j)); */
/* 	gsl_vector_set (rhs_y, i+j*NU,  v(i,j)); */
/*       } */

/*     } */
/*   } */
  
/*   /\* Storage in Compressed Column Storage (CCS) format for superlu *\/ */
/*   CCSProblem * css = ccs_problem_new (); */
/*   gint count = 0; */
/*   for ( i = 0; i < size; i++) { */
/*     g_array_append_val (css->index, count); */
/*     for ( j = 0; j < size; j++) { */
/*       if (A0[i][j] != 0.) { */
/*   	g_array_append_val (css->matrix, A0[i][j]); */
/*   	g_array_append_val (css->column, j); */
/*   	count++; */
/*       } */
/*     } */
/*   } */
/*   g_array_append_val (css->index, count); */

/*   ccs_problem_lu_solve (css, rhs_x); */
/*   ccs_problem_lu_solve (css, rhs_y); */

/*   for ( i = 0; i < NU; i++) { */
/*     for ( j = 0; j < NV; j++) { */
/*       coeff_assign (grid, i, j, 0, gsl_vector_get (rhs_x, i + j*NU)); */
/*       coeff_assign (grid, i, j, 1, gsl_vector_get (rhs_y, i + j*NU)); */
/*     } */
/*   } */

/*   ccs_problem_destroy (css); */
/*   gsl_vector_free (rhs_x); */
/*   gsl_vector_free (rhs_y); */




/*   /\* GridData gd; *\/ */
/*   /\* gd.lpx = lpx; *\/ */
/*   /\* gd.lpy = lpy; *\/ */
/*   /\* gd.M = M; *\/ */
/*   /\* gd.N = N; *\/ */

/*   /\* gint NU = gsl_bspline_ncoeffs (patch->w_u); *\/ */
/*   /\* gint NV = gsl_bspline_ncoeffs (patch->w_v); *\/ */

/*   /\* fprintf(stderr, "SIZE %i %i %i %i\n",M,N,NU,NV); *\/ */

/*   /\* if (flip) { *\/ */
/*   /\*   for ( i = 0; i < N; i++) { *\/ */
/*   /\*     for ( j = 0; j < M; j++) { *\/ */
/*   /\* 	coeff_assign (patch, i, j, 0, u(N-(i+1),j)); *\/ */
/*   /\* 	coeff_assign (patch, i, j, 1, v(N-(i+1),j)); *\/ */
/*   /\* 	coeff_assign (patch, i, j, 2, 0./\\* s->hz (u(N-(i+1),j), v(N-(i+1),j), NULL) *\\/); *\/ */
/*   /\* 	fprintf (stdout, "%f %f %f \n", u(N-(i+1),j), v(N-(i+1),j), 0.); *\/ */
/*   /\*     } *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */
/*   /\* else { *\/ */
/*   /\*   for ( i = 0; i < N; i++) { *\/ */
/*   /\*     for ( j = 0; j < M; j++) { *\/ */
/*   /\* 	coeff_assign (patch, i, j, 0, u(i/\\* N-(i+1) *\\/,j)); *\/ */
/*   /\* 	coeff_assign (patch, i, j, 1, v(i/\\* N-(i+1) *\\/,j)); *\/ */
/*   /\* 	coeff_assign (patch, i, j, 2, 0./\\* s->hz (u(i/\\\* N-(i+1) *\\\/,j), v(i/\\\* N-(i+1) *\\\/,j), NULL) *\\/); *\/ */
/*   /\* 	fprintf (stdout, "%f %f %f \n", u(i,j), v(i,j), 0.); *\/ */
/*   /\*     } *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */
/*   //g_assert_not_reached (); */

/*   // Spline2D * grid = spline2d_new (N-1, M-1, 2, 3, 3); */
/*   //Spline2D * grid =  periodic_fs_new (N-1, M-1, 3, 4, 3); */

/*   //spline2d_init_panels (grid); */
  

/*   //spline2d_fit_greville_border (grid, eval_grid_x, patch, 0); */
/*   //spline2d_fit_greville_border (grid, eval_grid_y, patch, 1); */

/*   /\* spline2d_fit_galerkin (grid, eval_grid_x_gauss, patch, 0); *\/ */
/*   /\* spline2d_fit_galerkin (grid, eval_grid_y_gauss, patch, 1); *\/ */

  

/*   /\* spline2d_destroy (patch); *\/ */

 
/*   //spline2d_reinit_panels_physical_quantities (grid); */
/*   periodic_fs_reinit_panels_physical_quantities (grid); */

/*   s->patches = g_slist_append (s->patches, grid); */
/* } */

static gsl_vector * finite_depth_free_surface_rhs (Spline2D * sp, WaveParams * wp, gdouble t)
{
  gint i, j, m, n, ii;
  size_t ustart, uend, vstart, vend;
  gint size = gsl_bspline_ncoeffs (sp->w_u)*gsl_bspline_ncoeffs (sp->w_v);
  gdouble RHS[size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_vector * Bu = gsl_vector_alloc (sp->k);
  gsl_vector * Bv = gsl_vector_alloc (sp->k);
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
    for ( m = 0; m < ng; m++) {
      gdouble um = g_array_index (gp->ui, gdouble, m);
      gsl_bspline_eval_nonzero (um, Bu, &ustart, &uend, sp->w_u); // Maybe store that
	  
      for ( n = 0; n < ng; n++) {
	gdouble vn = g_array_index (gp->vj, gdouble, n);
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	gsl_bspline_eval_nonzero (vn, Bv, &vstart, &vend, sp->w_v); // Maybe store that
	Point p = g_array_index (gp->Pi, Point, m+n*ng);
	
	gdouble fmn = wp->wave_elevation (p.x, p.y, t, wp) ;

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_vector_get(Bu, i-ustart);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_vector_get(Bv, j-vstart);
	    gint indexi = i + j*gsl_bspline_ncoeffs (sp->w_u);
	 
	    RHS[indexi] += wmnij*fmn;
	  }
	}
	
      }
    }
    
  }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  gint NU = gsl_bspline_ncoeffs (sp->w_u);
  gint NV = gsl_bspline_ncoeffs (sp->w_v);
  // Continuity of variable on boundary
  for ( i = 0; i < NV; i++) {
    gint index1 = (i+1)*NU - 1;
    gint index2 = i*NU;

    RHS[index1] += RHS[index2];
    RHS[index2] = 0.;
  }

  // Copy lhs and rhs to gsl structures
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, RHS[i]);

  return rhs;
}

static gdouble flat_bathymetry_rhs_gauss (SPPanel * spp, gint m, gint n, gpointer data)
{
  return -20.;
}

static gsl_vector * flat_bathymetry_rhs (Spline2D * sp, WaveParams * wp, gdouble t)
{
  gint i, j, m, n, ii;
  size_t ustart, uend, vstart, vend;
  gint size = gsl_bspline_ncoeffs (sp->w_u)*gsl_bspline_ncoeffs (sp->w_v);
  gdouble RHS[size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_vector * Bu = gsl_vector_alloc (sp->k);
  gsl_vector * Bv = gsl_vector_alloc (sp->k);
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
    for ( m = 0; m < ng; m++) {
      gdouble um = g_array_index (gp->ui, gdouble, m);
      gsl_bspline_eval_nonzero (um, Bu, &ustart, &uend, sp->w_u); // Maybe store that
	  
      for ( n = 0; n < ng; n++) {
	gdouble vn = g_array_index (gp->vj, gdouble, n);
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	gsl_bspline_eval_nonzero (vn, Bv, &vstart, &vend, sp->w_v); // Maybe store that
	Point p = g_array_index (gp->Pi, Point, m+n*ng);

	gdouble fmn = b_z (p.x, p.y, t, wp);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_vector_get(Bu, i-ustart);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_vector_get(Bv, j-vstart);
	    gint indexi = i + j*gsl_bspline_ncoeffs (sp->w_u);
	 
	    RHS[indexi] += wmnij*fmn;
	  }
	}
	
      }
    }
    
  }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  // Copy lhs and rhs to gsl structures
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, RHS[i]);

  return rhs;
}

static gdouble discretize_free_surface_rhs_gauss (SPPanel * spp,
						  gint m, gint n,
						  gpointer data)
{
  //  Simulation * sim = (Simulation *) data;
  WaveParams * wp = (WaveParams *) data;
  Spline2D * sp = spp->sp;
  GaussPoints * gp = spp->outer;
  /* gdouble u = g_array_index (gp->ui, gdouble, m); */
  /* gdouble v = g_array_index (gp->vj, gdouble, n); */
  Point p = g_array_index (gp->Pi, Point, m + n*sp->nouter);

  return wp->wave_elevation (p.x, p.y, 0./* sim->time.t */, wp) ; 
}

void spline2d_discretize_free_surface (Spline2D * sp, WaveParams * wp, gdouble t)
{

  if (!sp->fit)
    sp->fit = sp->build_fit_matrix (sp);

  gsl_vector * rhs = sp->build_fit_rhs (sp, discretize_free_surface_rhs_gauss, wp, NULL, NULL, rhs); // Need to pass the time as well

  ccs_problem_lu_solve (sp->fit, rhs);

   // Copy results to patch
  sp->copy_fit_solution (sp, rhs, 2);

  sp->reinit_panels (sp);// spline2d_reinit_panels_physical_quantities (sp);

  gsl_vector_free (rhs);

  
  //  spline2d_build_freesurface_noflux_galerkin_fit_matrix (sp);
}

void spline2d_discretize_bathymetry (Spline2D * sp, WaveParams * wp, gdouble t)
{
  if (!sp->fit)
    sp->fit = sp->build_fit_matrix (sp);

  //  gsl_vector * rhs = flat_bathymetry_rhs (sp, wp, t);
  gsl_vector * rhs = sp->build_fit_rhs (sp, flat_bathymetry_rhs_gauss, NULL, NULL, NULL, rhs);

  ccs_problem_lu_solve (sp->fit, rhs);

   // Copy results to patch
  sp->copy_fit_solution (sp, rhs, 2);

  gsl_vector_free (rhs);
}

void surface_generate_grid (Surface * s, gboolean flip)
{
  Boundaries * b = s->b;
  gint N = b->dcb->p->len;
  gint M = b->dcl->p->len;
  gint i;

  /* Linear problem in sparse matrix format */
  LinearProblem * lpx = linear_problem_new ();
  LinearProblem * lpy = linear_problem_new ();
  linear_problem_init_size (lpx, N, M);
  linear_problem_init_size (lpy, N, M);

  wild_grid_generation (b, lpx, lpy);
  //m_do (lpx, lpy, N, M);
  //print_grid (lpx, lpy, N, M);

#if 1
  /* Control functions */
  GArray * p = g_array_new (FALSE, FALSE, sizeof(double));
  GArray * q = g_array_new (FALSE, FALSE, sizeof(double));

  for ( i = 0; i < N*M; i++) {
    gdouble x0 = 0.;
    g_array_append_val (p, x0);
    g_array_append_val (q, x0);
  }

  /* Convergence boolean */
  gboolean converged = FALSE;  
  
  LinearProblem * lpx_new = build_x_elliptic_linear_system (lpx, lpy, p, q, N, M);
  LinearProblem * lpy_new = build_y_elliptic_linear_system (lpx, lpy, p, q, N, M);
  
  linear_problem_destroy (lpx);
  linear_problem_destroy (lpy);
  
  lpx = lpx_new;
  lpy = lpy_new;
  
  GArray * x0 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  GArray * y0 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  for (i = 0; i < lpx->lhs->len; i++) {
    g_array_append_val (x0, g_array_index (lpx->lhs, gdouble, i));
    g_array_append_val (y0, g_array_index (lpy->lhs, gdouble, i));
  }

  linear_problem_SOR_iteration (lpx);
  linear_problem_SOR_iteration (lpy);
 
  //  m_do (lpx, lpy, N, M);
  //print_grid (lpx, lpy, N, M);
  
  /* update_p (b, lpx, lpy, p, q, N, M); */
  /*   g_assert_not_reached (); */
  /* Start of the loop */
  gint niter = 0;
  while (!converged) {
    //sleep (1);
    for (i = 0; i < lpx->lhs->len; i++) {
      g_array_index (x0, gdouble, i) = g_array_index (lpx->lhs, gdouble, i);
      g_array_index (y0, gdouble, i) = g_array_index (lpy->lhs, gdouble, i);
    }

    /* update_p (b, lpx, lpy, p, q, N, M); */
    /* g_assert_not_reached (); */
    lpx_new = build_x_elliptic_linear_system (lpx, lpy, p, q, N, M);
    lpy_new = build_y_elliptic_linear_system (lpx, lpy, p, q, N, M);

    linear_problem_destroy (lpx);
    linear_problem_destroy (lpy);
    lpx = lpx_new;
    lpy = lpy_new;

    linear_problem_SOR_iteration (lpx);
    linear_problem_SOR_iteration (lpy);

    /* print_grid (lpx, lpy, N, M); */
   

    gdouble d = 0.;
    for (i = 0; i < lpx->lhs->len; i++) {
      gdouble x1 = g_array_index (x0, gdouble, i);
      gdouble y1 = g_array_index (y0, gdouble, i);
      gdouble x2 = g_array_index (lpx->lhs, gdouble, i);
      gdouble y2 = g_array_index (lpy->lhs, gdouble, i);

      if (sqrt ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) > d) {
  	d = sqrt ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
      }
    }
    /* fprintf (stdout,"Delta: %i %g\n", niter, d); */
    /* For now only one iteration */
    
    /* gdouble mdo = m_do (lpx, lpy, N, M); */
    /* fprintf (stdout,"Delta: %i %g %g\n", niter, d, mdo); */

    niter++;
    if (d < 5.e-3)
      converged = TRUE;
  }
#endif
  //print_grid (lpx, lpy, N, M);
  /* store_grid (s, lpx, lpy, N, M); */
  spline2d_store_grid (s, lpx, lpy, N, M, flip);

  /* FILE * fp = fopen ("grid.tmp", "w"); */
  /* Patch * patch = g_ptr_array_index (s->patches, 0); */
  /* patch_print (patch, fp); */
  /* fclose(fp); */
  

  linear_problem_destroy (lpx);
  linear_problem_destroy (lpy);
#if 0
  lpx_new = lpy_new = NULL;
  g_array_free (p, TRUE);
  g_array_free (q, TRUE);
  g_array_free (x0, TRUE);
  g_array_free (y0, TRUE);
#endif
}




struct dispersion_params
{
  double omega, H, g;
};

static double dispersion (double k, void *params)
{
  struct dispersion_params * p 
    = (struct dispersion_params *) params;

  return k*tanh(k*p->H) - pow(p->omega,2)/p->g;
}
     
gdouble solve_dispersion_relation (WaveParams * wp)
{
  int status;
  int iter = 0, max_iter = 2000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0.;
  double x_lo = 0.0, x_hi = 100.0;
  
  gsl_function F;
  struct dispersion_params params = {wp->w, wp->h, wp->g};
  
  F.function = &dispersion;
  F.params = &params;
     
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
         
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi,
				     0, 0.001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
     
  g_assert (status != GSL_CONTINUE);
  
  gsl_root_fsolver_free (s);

  return r;
}

gdouble spline2d_eval2 (Spline2D * s, gdouble u, gdouble v, gint var, gint num)
{
  /* g_assert ( u >= -1 && u <= 1.); */
  /* if ( v > 1 || v < -1) */
  /*   fprintf(stderr,"%f \n",v); */

  /* g_assert ( v >= -1 && v <= 1.); */
  g_assert ( s->w_u != NULL);
  g_assert ( s->w_v != NULL);

  gint i, j;
  size_t istart, iend, jstart, jend;
  gsl_vector * Bu = gsl_vector_alloc (s->k);
  gsl_vector * Bv = gsl_vector_alloc (s->k);
  gdouble val = 0.;

  gsl_bspline_eval_nonzero (MIN(1.-1e-12,u), Bu, &istart, &iend, s->w_u);
  gsl_bspline_eval_nonzero (MIN(1.-1e-12,v), Bv, &jstart, &jend, s->w_v);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_vector_get(Bu, i);
    for ( j = 0; j < s->k; j++) {
      gdouble cv = gsl_vector_get(Bv, j);
      if (num ==  (istart+i) % (s->M+s->k-1))
	val += coeff (s, (istart+i) % (s->M+s->k-1),(jstart+j) % (s->N+s->k-1),var)*cu*cv;
    }
  }

  gsl_vector_free (Bu);
  gsl_vector_free (Bv);
  return val;
}

gdouble spline2d_derivative_eval2 (Spline2D * s, gdouble u, gdouble v, gint m, gint n, gint var, gint num)
{
  g_assert ( u >= -1e-8 && u <= 1.+1e-8);
  g_assert ( v >= -1e-8 && v <= 1.+1e-8);
  /* g_assert ( u >= 0. && u <= 1.); */
  /* g_assert ( v >= 0. && v <= 1.); */
  g_assert ( s->w_u != NULL);
  g_assert ( s->w_v != NULL);

  gint  i, j;
  gdouble val = 0.;
  size_t istart, iend, jstart, jend;
  gsl_matrix * Bu = gsl_matrix_alloc (s->k, m+1);
  gsl_matrix * Bv = gsl_matrix_alloc (s->k, n+1);

  gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,u), m, Bu, &istart, &iend, s->w_u, s->wd_u);
  gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), n, Bv, &jstart, &jend, s->w_v, s->wd_v);

  //fprintf (stderr, "index: %i\n", istart, istart+s->k-1);

  for ( i = 0; i < s->k; i++) {
    gdouble cu = gsl_matrix_get (Bu, i, m);
    gint ii = istart % (s->M+s->k-1);
    for ( j = 0; j < s->k; j++) {
      gint jj = (jstart+j) % (s->N+s->k-1);
      gdouble cv = gsl_matrix_get (Bv, j, n);
      if (num == ii)
	val += coeff (s, ii, jj,var)*cu*cv;
    }
    istart++;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);
  return val;
}

gdouble tmp_rhs (SPPanel * spp, gint m, gint n, gpointer data)
{
  Spline2D * sp = spp->sp;
  GaussPoints * gp = spp->outer;

  return spline2d_eval_gauss_point (sp, gp, m, n, 3);
}

void print_free_surface_potential (GSList * list, WaveParams * wp, gdouble t)
{
  FILE * fp = fopen ("potential.tmp","w");
  Spline2D * sp = list->data;
  
  gint var = 3/* 15 *//* 9 */; //9 = zeta
  gint var2 = 4/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      
  sp->noflux = FALSE;

  /* CCSProblem * fit = sp->build_fit_matrix (sp); */
  /* gsl_vector * rhs = sp->build_fit_rhs (sp, tmp_rhs, NULL, NULL, NULL); */
  /* ccs_problem_lu_solve (fit, rhs); */
  /* sp->copy_fit_solution (sp, rhs, 3); */
  /* sp->noflux = TRUE; */

  /* ccs_problem_destroy (fit); */
  /* gsl_vector_free (rhs); */

  gint m, n;
  m = 1;
  n = 0;



  gint i, j;
  while (sp) {
    for ( i = 0; i < sp->M; i++) {
      for ( j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
	
	Point p = spline2d_eval_point (sp, spp->u0, spp->v0);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v0, var), spline2d_eval (sp, spp->u0, spp->v0, var2), spline2d_eval (sp, spp->u0, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
  
	p = spline2d_eval_point (sp, spp->u1, spp->v0);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u1, spp->v0, var), spline2d_eval (sp, spp->u1, spp->v0, var2), spline2d_eval (sp, spp->u1, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
  
	p = spline2d_eval_point (sp, spp->u1, spp->v1);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u1, spp->v1, var), spline2d_eval (sp, spp->u1, spp->v1, var2), spline2d_eval (sp, spp->u1, spp->v1, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	
	p = spline2d_eval_point (sp, spp->u0, spp->v1);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v1, var), spline2d_eval (sp, spp->u0, spp->v1, var2), spline2d_eval (sp, spp->u0, spp->v1, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
      
	p = spline2d_eval_point (sp, spp->u0, spp->v0);
	fprintf(fp, "%f %f %f %f %f %f\n\n\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v0, var), spline2d_eval (sp, spp->u0, spp->v0, var2), spline2d_eval (sp, spp->u0, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	
      }
    }
    sp = sp->next;
  }

  fclose (fp);
}

void print_hull_potential (GSList * list, WaveParams * wp, gdouble t)
{
  FILE * fp = fopen ("hull_potential.tmp","w");
  /* Spline2D * sp = list->data; */
  
  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      
  /* sp->noflux = FALSE; */

  /* CCSProblem * fit = sp->build_fit_matrix (sp); */
  /* gsl_vector * rhs = sp->build_fit_rhs (sp, tmp_rhs, NULL, NULL, NULL); */
  /* ccs_problem_lu_solve (fit, rhs); */
  /* sp->copy_fit_solution (sp, rhs, 3); */
  /* sp->noflux = TRUE; */

  /* ccs_problem_destroy (fit); */
  /* gsl_vector_free (rhs); */

  gint m, n;
  m = 1;
  n = 0;



  gint i, j;

  GSList * splines = list;
  while (splines) {
    Spline2D * sp = splines->data;
    for ( i = 0; i < sp->M; i++) {
      for ( j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
	
	Point p = spline2d_eval_point (sp, spp->u0, spp->v0);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v0, var), spline2d_eval (sp, spp->u0, spp->v0, var2), spline2d_eval (sp, spp->u0, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
  
	p = spline2d_eval_point (sp, spp->u1, spp->v0);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u1, spp->v0, var), spline2d_eval (sp, spp->u1, spp->v0, var2), spline2d_eval (sp, spp->u1, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
  
	p = spline2d_eval_point (sp, spp->u1, spp->v1);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u1, spp->v1, var), spline2d_eval (sp, spp->u1, spp->v1, var2), spline2d_eval (sp, spp->u1, spp->v1, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	
	p = spline2d_eval_point (sp, spp->u0, spp->v1);
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v1, var), spline2d_eval (sp, spp->u0, spp->v1, var2), spline2d_eval (sp, spp->u0, spp->v1, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
      
	p = spline2d_eval_point (sp, spp->u0, spp->v0);
	fprintf(fp, "%f %f %f %f %f %f\n\n\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v0, var), spline2d_eval (sp, spp->u0, spp->v0, var2), spline2d_eval (sp, spp->u0, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	
      }
    }
    splines = splines->next;
  }

  fclose (fp);
}


void print_free_surface (GSList * list, WaveParams * wp, gdouble t)
{
  //FILE * fp = fopen ("freesurface.tmp","w");
  FILE * fp = fopen (g_strdup_printf ("freesurface_%5.4f.tmp", t),"w");
  Spline2D * sp = list->data;
  
  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      

  

  
  gint i, j;
  while (sp) {
    for ( i = 0; i < sp->M; i++) {
      for ( j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);

	Point p = spline2d_eval_point (sp, spp->u0, spp->v0);
	fprintf(fp, "%e %e %e %e %e %e\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v0, var), spline2d_eval (sp, spp->u0, spp->v0, var2), spline2d_eval (sp, spp->u0, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
  
	p = spline2d_eval_point (sp, spp->u1, spp->v0);
	fprintf(fp, "%e %e %e %e %e %e\n", p.x, p.y,
		spline2d_eval (sp, spp->u1, spp->v0, var), spline2d_eval (sp, spp->u1, spp->v0, var2), spline2d_eval (sp, spp->u1, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
      
	p = spline2d_eval_point (sp, spp->u1, spp->v1);
	fprintf(fp, "%e %e %e %e %e %e\n", p.x, p.y,
		spline2d_eval (sp, spp->u1, spp->v1, var), spline2d_eval (sp, spp->u1, spp->v1, var2), spline2d_eval (sp, spp->u1, spp->v1, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	      
	p = spline2d_eval_point (sp, spp->u0, spp->v1);
	fprintf(fp, "%e %e %e %e %e %e\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v1, var), spline2d_eval (sp, spp->u0, spp->v1, var2), spline2d_eval (sp, spp->u0, spp->v1, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	      
	p = spline2d_eval_point (sp, spp->u0, spp->v0);
	fprintf(fp, "%e %e %e %e %e %e\n\n\n", p.x, p.y,
		spline2d_eval (sp, spp->u0, spp->v0, var), spline2d_eval (sp, spp->u0, spp->v0, var2), spline2d_eval (sp, spp->u0, spp->v0, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
      }
    }
    sp = sp->next; 
  }

  fclose (fp);
  /* fp = fopen ("color.tmp","w"); */
  /* for ( u = 0; u < 1; u += du) { */
  /*   for ( v = 0; v < 1; v += du) { */
  /*     Point p = spline2d_eval_point (sp, u, v); */
  /*     fprintf(fp, "%f %f %f \n", p.x, p.y, finite_depth_wave_elevation (p.x, p.y, t, wp) + */
  /* 	      spline2d_eval (sp, u, v, 9)); */
  /*   } */
  /*   fprintf (fp, "\n"); */
  /* } */
  /* fclose (fp); */
}

void print_free_surface_mayavi (GSList * list, WaveParams * wp, gdouble t)
{
  //FILE * fp = fopen ("freesurface.tmp","w");
  FILE * fp = fopen (g_strdup_printf ("fs_mayavi_%5.4f.tmp", t),"w");
  Spline2D * sp = list->data;
  
  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2  

  
  gint i, j;
  gdouble u, v, du = 0.01, dv = 0.02;

  gint NU = 0, NV;

  sp = list->data;
  while (sp) {
    NU += (gint) 1./du;
    sp = sp->next;
  }
  NU += 1; // Close the patch
  NV = 1./dv + 1;

  fprintf (fp,"# NX %i NY %i NZ 1 \n", NU, NV);
  
  for ( v = 0.; v <= 1.0001; v+= 0.02) {
    sp = list->data;
    while (sp) {
      for ( u = 0.; u < 1.0001; u += 0.01 ) {
	Point p =spline2d_eval_point (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.));
	fprintf (fp, "%f %f %f %f \n", p.x, p.y, spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var),
		 spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var2));
      }
      if (sp->next == NULL) {
	Point p =spline2d_eval_point (sp, 1, MIN(MAX(v,0.),1.));
	fprintf (fp, "%f %f %f %f \n", p.x, p.y, spline2d_eval (sp, 1, MIN(MAX(v,0.),1.), var),
		 spline2d_eval (sp, 1, MIN(MAX(v,0.),1.), var2));
      }
      sp = sp->next;
    }
  }

  fclose (fp);
}

void print_waterline (GSList * list, WaveParams * wp, gdouble t)
{
  //FILE * fp = fopen ("freesurface.tmp","w");
  FILE * fp = fopen (g_strdup_printf ("waterline_%5.4f.tmp", t),"w");
  Spline2D * sp = list->data;
  
  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      
  gint i, j;
  while (sp) {
    gdouble u, v = 0.;

    for ( u = 0.; u <= 1.; u += 0.01) {
      Point p = spline2d_eval_point (sp, u, v);
      fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
	      spline2d_eval (sp, u, v, var), spline2d_eval (sp, u, v, var2), spline2d_eval (sp, u, v, var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
    }

    sp = sp->next;
  }

  fclose (fp);
}

void print_free_surface_tecplot (GSList * list, WaveParams * wp, gdouble t)
{
  FILE * fp = fopen (g_strdup_printf ("freesurface_%5.4f.tp", t),"w");
  Spline2D * sp = list->data;
  gint MM = 0;
  
  while (sp) {
    MM += (sp->M+1);
    sp = sp->next;
  }

  sp = list->data;
  fprintf (fp, "TITLE = \"Free-surface 2D Plot for paraview\"\n");
  fprintf (fp, "VARIABLES = \"X\", \"Y\", \"Eta\", \"Phi2\", \"Phin\", \"eta0\"\n");
  fprintf (fp, "ZONE T=\"BIG ZONE\", I=%i, J=%i, F=POINT\n",sp->N+1,MM);

  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      
  
  gint i, j;  
  while (sp) {
    
    for ( i = -1; i < sp->M; i++ ) {
      gdouble u = gsl_vector_get (sp->wx_u->knots, sp->k + i);
      for ( j = -1; j < sp->N; j++ ) {
	gdouble v = gsl_vector_get (sp->w_v->knots, sp->k + j);
	Point p = spline2d_eval_point (sp, u, v);
	
	fprintf (fp, "%f %f %f %f %f %f\n", p.x, p.y,
		 spline2d_eval (sp, u, v, var),
		 spline2d_eval (sp, u, v, var2),
		 spline2d_eval (sp, u, v, var3),
		 finite_depth_wave_elevation (p.x, p.y, t, wp));
      }
    }
    sp = sp->next;
  }
  fclose (fp);
}

void print_free_surface_tmp (GSList * list, WaveParams * wp, gdouble t)
{
  FILE * fp = fopen ("freesurface2.tmp","w");
  Spline2D * sp = list->data;
  
  gint var = 3; //9 = zeta
  gint var2 = 15/* 4 *//* 4 */; //7 = Phi2
      
  sp->noflux = FALSE;

  /* CCSProblem * fit = sp->build_fit_matrix (sp); */
  /* gsl_vector * rhs = sp->build_fit_rhs (sp, tmp_rhs, NULL); */
  /* ccs_problem_lu_solve (fit, rhs); */
  /* sp->copy_fit_solution (sp, rhs, 3); */
  /* sp->noflux = TRUE; */

  /* ccs_problem_destroy (fit); */
  /* gsl_vector_free (rhs); */


  
  gint m, n;
  m = 1;
  n = 0;

  gint i, j;
  while (sp) {
    gdouble u, v, du = 0.01;
    for ( u = 0; u < 0.9999; u += du ) {
      gdouble u0 = MIN(MAX(u,0.),1.);
      gdouble u1 = MIN(MAX(u+du,0.),1.);
      for ( v = 0; v < 0.999; v += du ) {
	gdouble v0 = MIN(MAX(v,0.),1.);
	gdouble v1 = MIN(MAX(v+du,0.),1.);
	
	gdouble der = spline2d_derivative_eval (sp, u0, v0, 1, 0, var);
	gdouble der2 = spline2d_derivative_eval (sp, u0, v0, 0, 1, var);
	Vector grad =  potential_gradient_on_surface (sp, u0, v0, var);
	Point p = spline2d_eval_point (sp, u0, v0);
	//gdouble der = spline2d_eval (sp, spp->u0, spp->v0, var+1);
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f \n", p.x, p.y, p.z, spline2d_eval (sp, u0, v0, var), spline2d_eval (sp, u0, v0, var2),
		der, der2, u0, v0, grad.x, grad.y, grad.z);
  
	der = spline2d_derivative_eval (sp, u1, v0, 1, 0, var);
	der2 = spline2d_derivative_eval (sp, u1, v0, 0, 1, var);
	grad =  potential_gradient_on_surface (sp, u1, v0, var);
	p = spline2d_eval_point (sp, u1, v0);
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f \n", 
		p.x, p.y, p.z,
		spline2d_eval (sp, u1, v0, var), spline2d_eval (sp, u1, v0, var2),
		der, der2, u1, v0, grad.x, grad.y, grad.z);
      
	der = spline2d_derivative_eval (sp, u1, v1, 1, 0, var);
	der2 = spline2d_derivative_eval (sp, u1, v1, 0, 1, var);
	grad =  potential_gradient_on_surface (sp, u1, v1, var);
	p = spline2d_eval_point (sp, u1, v1);
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f \n", 
		p.x, p.y, p.z,
		spline2d_eval (sp, u1, v1, var), spline2d_eval (sp, u1, v1, var2),
		der, der2, u1, v1, grad.x, grad.y, grad.z);
      
	der = spline2d_derivative_eval (sp, u0, v1, 1, 0, var);
	der2 = spline2d_derivative_eval (sp, u0, v1, 0, 1, var);
	grad =  potential_gradient_on_surface (sp, u0, v1, var);
	p = spline2d_eval_point (sp, u0, v1);
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f \n", 
		p.x, p.y, p.z,
		spline2d_eval (sp, u0, v1, var), spline2d_eval (sp, u0, v1, var2),
		der, der2, u0, v1, grad.x, grad.y, grad.z);
      
	der = spline2d_derivative_eval (sp, u0, v0, 1, 0, var);
	der2 = spline2d_derivative_eval (sp, u0, v0, 0, 1, var);
	grad =  potential_gradient_on_surface (sp, u0, v0, var);
	p = spline2d_eval_point (sp, u0, v0);
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f \n\n\n",
		p.x, p.y, p.z,
		spline2d_eval (sp, u0, v0, var), spline2d_eval (sp, u0, v0, var2),
		der, der2, u0, v0, grad.x, grad.y, grad.z);
      
      }
    }
    sp = sp->next;
  }

  fclose (fp);

}

gdouble zero_scalar_wave_func (WaveParams * wp,
			       Point p, gdouble t)
{
  return 0.;
}

Vector  zero_vector_wave_func (WaveParams * wp,
			       Point p, gdouble t)
{
  Vector v;
  v.x = v.y = v.z = 0.;
  return v;
}

gdouble zero_wave_elevation (gdouble x, gdouble y,
			     gdouble t,
			     gpointer data
			     /* WaveParams * wp */)
{
  return 0.;
}

void print_free_surface_hr (GSList * list, WaveParams * wp, gdouble t)
{
  //FILE * fp = fopen ("freesurface.tmp","w");
  FILE * fp = fopen (g_strdup_printf ("freesurface_hr_%5.4f.tmp", t),"w");
  Spline2D * sp = list->data;
  
  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      

  

  
  gint i, j;
  while (sp) {
    gdouble u, v, du = 0.01;
    for ( u = 0; u < 0.9999; u += du ) {
      for ( v = 0; v < 0.999; v += du ) {
	Point p = spline2d_eval_point (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.));
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var), spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var2), spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
  
	p = spline2d_eval_point (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v,0.),1.));
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v,0.),1.), var), spline2d_eval (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v,0.),1.), var2), spline2d_eval (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v,0.),1.), var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
      
	p = spline2d_eval_point (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v+du,0.),1.));
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v+du,0.),1.), var), spline2d_eval (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v+du,0.),1.), var2), spline2d_eval (sp, MIN(MAX(u+du,0.),1.), MIN(MAX(v+du,0.),1.), var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	      
	p = spline2d_eval_point (sp, MIN(MAX(u,0.),1.), MIN(MAX(v+du,0.),1.));
	fprintf(fp, "%f %f %f %f %f %f\n", p.x, p.y,
		spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v+du,0.),1.), var), spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v+du,0.),1.), var2), spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v+du,0.),1.), var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
	      
	p = spline2d_eval_point (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.));
	fprintf(fp, "%f %f %f %f %f %f\n\n\n", p.x, p.y,
		spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var), spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var2), spline2d_eval (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.), var3), finite_depth_wave_elevation (p.x, p.y, t, wp));
      }
    }
    sp = sp->next; 
  }

  fclose (fp);
}

void print_free_surface_tecplot_hr (GSList * list,
				    WaveParams * wp, gdouble t)
{
  FILE * fp = fopen (g_strdup_printf ("freesurface_hr_%5.4f.tp", t),"w");
  Spline2D * sp = list->data;
  gint MM = 0;
  gdouble u, v, du = 0.01;
  
  while (sp) {
    MM += /* (sp->M+1) */ (gint) (1./du+1);
    sp = sp->next;
  }
  gint N = (gint) (1./du+1);

  sp = list->data;
  fprintf (fp, "TITLE = \"Free-surface 2D Plot for paraview\"\n");
  fprintf (fp, "VARIABLES = \"X\", \"Y\", \"Eta\", \"Phi2\", \"Phin\", \"eta0\"\n");
  fprintf (fp, "ZONE T=\"BIG ZONE\", I=%i, J=%i, F=POINT\n",N,MM);

  gint var = 9/* 15 *//* 9 */; //9 = zeta
  gint var2 = 7/* 14 *//* 7 */; //7 = Phi2
  gint var3 = 8/* 8 *//* 14 *//* 8 */;
      
  
  gint i, j;  
  while (sp) {
    
    for ( u = 0; u < 1.000001; u += du ) {
      for ( v = 0; v < 1.00001; v += du ) {

    /* for ( i = -1; i < sp->M; i++ ) { */
    /*   gdouble u = gsl_vector_get (sp->wx_u->knots, sp->k + i); */
    /*   for ( j = -1; j < sp->N; j++ ) { */
    /* 	gdouble v = gsl_vector_get (sp->w_v->knots, sp->k + j); */
	Point p = spline2d_eval_point (sp, MIN(MAX(u,0.),1.), MIN(MAX(v,0.),1.));
	
	fprintf (fp, "%f %f %f %f %f %f\n", p.x, p.y,
		 spline2d_eval (sp, MIN(MAX(u,0.),1.),
				MIN(MAX(v,0.),1.), var),
		 spline2d_eval (sp, MIN(MAX(u,0.),1.),
				MIN(MAX(v,0.),1.), var2),
		 spline2d_eval (sp, MIN(MAX(u,0.),1.),
				MIN(MAX(v,0.),1.), var3),
		 finite_depth_wave_elevation (p.x, p.y, t, wp));
      }
    }
    sp = sp->next;
  }
  fclose (fp);
}

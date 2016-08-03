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
#include "hull.h"
#include "motion.h"


static Forces forces (gdouble t, gdouble t0, gdouble dt, gdouble v, gdouble x)
{
  Forces f;

  gdouble xi = 2;
  gdouble omega = M_PI;


  /* Test for harmonic oscillator */
  /* Type of formula we will be using */
  f.fv = -2.*xi*v - x + (t0+dt-t)/dt*cos(omega*t0) + (t-t0)/dt*cos(omega*(t0+dt));
  f.fx = v;

  /* Standard analytical */
  /* f.fv = -2.*xi*v - x + cos(omega*t); */
  /* f.fx = v; */

  return f;
}

Forces forces_sum (Forces f1, Forces f2)
{
  Forces fsum;

  // Forces
  fsum.Fm = vector_sum (f1.Fm, f2.Fm);
  fsum.Fl = vector_sum (f1.Fl, f2.Fl);
  fsum.Fh = vector_sum (f1.Fh, f2.Fh);

  // Moments
  fsum.Mm = vector_sum (f1.Mm, f2.Mm);
  fsum.Ml = vector_sum (f1.Ml, f2.Ml);
  fsum.Mh = vector_sum (f1.Mh, f2.Mh);

  gint i;
  for ( i = 0; i < 6; i++) {
    fsum.forces_m[i] = f1.forces_m[i] + f2.forces_m[i];
    fsum.forces_l[i] = f1.forces_l[i] + f2.forces_l[i];
    fsum.forces_h[i] = f1.forces_h[i] + f2.forces_h[i];
    fsum.forces_fk[i] = f1.forces_fk[i] + f2.forces_fk[i];
    fsum.forces_r[i] = f1.forces_r[i] + f2.forces_r[i];
    fsum.forces_e[i] = f1.forces_e[i] + f2.forces_e[i];
    fsum.forces_ext[i] = f1.forces_ext[i] + f2.forces_ext[i];
    fsum.phi2[i] += f1.phi2[i] + f2.phi2[i];
  }

  return fsum;
}

void forces_set_to_zero (Forces * f)
{
  vector_set_to_zero (&f->Fm);
  vector_set_to_zero (&f->Fl);
  vector_set_to_zero (&f->Fh);
  
  vector_set_to_zero (&f->Mm);
  vector_set_to_zero (&f->Ml);
  vector_set_to_zero (&f->Mh);

  f->fv = 0.;
  f->fx = 0.;

  gint i;
  for ( i = 0; i < 6; i++) {
    f->v[i] = f->u[i] = f->forces_m[i] = f->forces_l[i] = f->forces_h[i] = 0.;
    f->forces_fk[i] = f->forces_r[i] = f->forces_e[i] = f->forces_ext[i] = 0.;
    f->x[i] = f->a[i] = f->phi2[i] = f->phi1[i] = 0.;
    f->forces_1[i] = f->forces_fk1[i] = f->forces_ext1[i] = f->forces_h1[i] = f->forces_v1[i] = f->x1[i] = f->u1[i] = f->gradphi1[i].x = f->gradphi1[i].y = f->gradphi1[i].z = f->v1[i] = 0.;
    f->forces_2[i] = f->forces_fk2[i] = f->forces_ext2[i] = f->forces_h2[i] = f->x2[i] = f->forces_v2[i] = f->u2[i] = f->v2[i] = 0.;
  }
}

void forces_print (Forces * f, FILE * fp)
{
  fprintf (fp, "Pressure forces: \n");
  fprintf (fp, "  memory:        %e %e %e\n", f->forces_m[0], f->forces_m[1], f->forces_m[2]);
  fprintf (fp, "  local:         %e %e %e\n", f->forces_l[0], f->forces_l[1],f->forces_l[2] );
  fprintf (fp, "  hydrostatic:   %e %e %e\n", f->forces_h[0], f->forces_h[1], f->forces_h[2]);
  fprintf (fp, "  froude-krylov: %e %e %e\n", f->forces_fk[0], f->forces_fk[1],f->forces_fk[2] );
  fprintf (fp, "  restoring:     %e %e %e\n", f->forces_r[0], f->forces_r[1], f->forces_r[2]);
  fprintf (fp, "  entrainement:  %e %e %e\n", f->forces_e[0], f->forces_e[1], f->forces_e[2]);
  fprintf (fp, "  external:      %e %e %e\n", f->forces_ext[0], f->forces_ext[1], f->forces_ext[2]);
  fprintf (fp, "Pressure moment: \n");
  fprintf (fp, "  memory:        %e %e %e\n", f->forces_m[3], f->forces_m[4], f->forces_m[5]);
  fprintf (fp, "  local:         %e %e %e\n", f->forces_l[3], f->forces_l[4], f->forces_l[5]);
  fprintf (fp, "  hydrostatic:   %e %e %e\n", f->forces_h[3], f->forces_h[4], f->forces_h[5]);
  fprintf (fp, "  froude-krylov: %e %e %e\n", f->forces_fk[3], f->forces_fk[4],f->forces_fk[5] );
  fprintf (fp, "  restoring:     %e %e %e\n", f->forces_r[3], f->forces_r[4], f->forces_r[5]);
  fprintf (fp, "  entrainement:  %e %e %e\n", f->forces_e[3], f->forces_e[4], f->forces_e[5]);
  fprintf (fp, "  external:      %e %e %e\n", f->forces_ext[3], f->forces_ext[4], f->forces_ext[5]);
}

void forces_total_print (Forces * f, FILE * fp, Time * t)
{
  fprintf (fp, " %f %i %e %e %e %e %e %e\n",
	   t->t,
	   t->itime,
	   f->forces_m[0]+f->forces_l[0]+f->forces_h[0]+f->forces_fk[0]
	   +f->forces_r[0]+f->forces_e[0]+f->forces_ext[0],
	   f->forces_m[1]+f->forces_l[1]+f->forces_h[1]+f->forces_fk[1]
	   +f->forces_r[1]+f->forces_e[1]+f->forces_ext[1],
	   f->forces_m[2]+f->forces_l[2]+f->forces_h[2]+f->forces_fk[2]
	   +f->forces_r[2]+f->forces_e[2]+f->forces_ext[2],
	   f->forces_m[3]+f->forces_l[3]+f->forces_h[3]+f->forces_fk[3]
	   +f->forces_r[3]+f->forces_e[3]+f->forces_ext[3],
	   f->forces_m[4]+f->forces_l[4]+f->forces_h[4]+f->forces_fk[4]
	   +f->forces_r[4]+f->forces_e[4]+f->forces_ext[4],
	   f->forces_m[5]+f->forces_l[5]+f->forces_h[5]+f->forces_fk[5]
	   +f->forces_r[5]+f->forces_e[5]+f->forces_ext[5]);
  fflush (fp);
}

void forces_add_contribution_forces (Forces * f, Vector * d,
				     gdouble wm, gdouble wl, gdouble wh,
				     gdouble wfk, gdouble wr)
{
  f->forces_m[0] += wm*d->x;
  f->forces_m[1] += wm*d->y;
  f->forces_m[2] += wm*d->z;

  f->forces_l[0] += wl*d->x;
  f->forces_l[1] += wl*d->y;
  f->forces_l[2] += wl*d->z;

  f->forces_h[0] += wh*d->x;
  f->forces_h[1] += wh*d->y;
  f->forces_h[2] += wh*d->z;

  f->forces_fk[0] += wfk*d->x;
  f->forces_fk[1] += wfk*d->y;
  f->forces_fk[2] += wfk*d->z;

  f->forces_r[0] += wr*d->x;
  f->forces_r[1] += wr*d->y;
  f->forces_r[2] += wr*d->z;
}

void forces_add_contribution_moments (Forces * f, Vector * d,
				      gdouble wm, gdouble wl, gdouble wh,
				      gdouble wfk, gdouble wr)
{
  f->forces_m[3] += wm*d->x;
  f->forces_m[4] += wm*d->y;
  f->forces_m[5] += wm*d->z;

  f->forces_l[3] += wl*d->x;
  f->forces_l[4] += wl*d->y;
  f->forces_l[5] += wl*d->z;

  f->forces_h[3] += wh*d->x;
  f->forces_h[4] += wh*d->y;
  f->forces_h[5] += wh*d->z;

  f->forces_fk[3] += wfk*d->x;
  f->forces_fk[4] += wfk*d->y;
  f->forces_fk[5] += wfk*d->z;

  f->forces_r[3] += wr*d->x;
  f->forces_r[4] += wr*d->y;
  f->forces_r[5] += wr*d->z;
}

void add_force_to_history (ForcesHistory * fh, Forces f)
{
  if (g_slist_length (fh->f) < 4) {
    Forces * new = g_malloc (sizeof(Forces));
    *new = f;
    fh->f = g_slist_append (fh->f, new);
  }
  else {
    Forces * new = fh->f->data;
    fh->f = g_slist_remove (fh->f, new);
    *new = f;
    fh->f = g_slist_append (fh->f, new);
  }
}

void forces_history_destroy (ForcesHistory * fh)
{
  GSList * list = fh->f;

  while (list) {
    Forces * f = list->data;
    g_free (f);
    list = list->next;
  }

  g_slist_free (fh->f);
}

static Motion RK4 (gdouble t, gdouble dt, Motion m0, Forces f0, Forces f1)
{
  Motion m1, k1, k2, k3, k4;
  Forces ftmp;

  // k1 = dt*f(t0, y0)
  ftmp = forces (t, t, dt, m0.v0, m0.x0);
  k1.v0 = dt*ftmp.fv;
  k1.x0 = dt*ftmp.fx;
  
  // k2 = dt*f(t0+dt/2, y0)
  ftmp = forces (t+dt/2, t, dt, m0.v0, m0.x0);
  k2.v0 = dt*ftmp.fv;
  k2.x0 = dt*ftmp.fx;
  
  // k3 = dt*f(t0+dt/2, y0+k2/2)
  ftmp = forces (t+dt/2, t, dt, m0.v0+k2.v0, m0.x0+k2.x0);
  k3.v0 = dt*ftmp.fv;
  k3.x0 = dt*ftmp.fx;

  // k4 = dt*f(t0+dt, y0+k3)
  ftmp = forces (t+dt, t, dt, m0.v0+k3.v0, m0.x0+k3.x0);
  k4.v0 = dt*ftmp.fv;
  k4.x0 = dt*ftmp.fx;

  // y1 = y0 + k1/6 + k2/3 + k3/3 +k4/6
  m1.v0 = m0.v0 + k1.v0/6. + k2.v0/3 + k3.v0/3 + k4.v0/6.;
  m1.x0 = m0.x0 + k1.x0/6. + k2.x0/3 + k3.x0/3 + k4.x0/6.;

  return m1;
}

static Forces forces_interpolation (gdouble t, gdouble t0,
				    gdouble dt, Forces f0, Forces f1,
				    gdouble e[6], gdouble v[6])
{
  Forces f;

  gdouble xi = 2;
  gdouble omega = M_PI;
  gint i;

  for ( i = 0; i < 6; i++ ) {
    f.forces_m[i] = (t0+dt-t)/dt*f0.forces_m[i] + (t-t0)/dt*f1.forces_m[i];
    f.forces_l[i] = (t0+dt-t)/dt*f0.forces_l[i] + (t-t0)/dt*f1.forces_l[i];
    f.forces_h[i] = (t0+dt-t)/dt*f0.forces_h[i] + (t-t0)/dt*f1.forces_h[i];
    f.forces_fk[i] = (t0+dt-t)/dt*f0.forces_fk[i] + (t-t0)/dt*f1.forces_fk[i];
    f.forces_r[i] = (t0+dt-t)/dt*f0.forces_r[i] + (t-t0)/dt*f1.forces_r[i];
    //   f.forces_e[i] = (t0+dt-t)/dt*f0.forces_e[i] + (t-t0)/dt*f1.forces_e[i];
    // WHY ?????
    f.forces_ext[i] = (t0+dt-t)/dt*f0.forces_ext[i] + (t-t0)/dt*f1.forces_ext[i];
    f.v[i] = v[i];
  }

  return f;
}

static Motion AMB4 (gdouble t, gdouble dt, Motion m0, ForcesHistory * fh)
{
  gint i = 0;
  Forces f[4];

  GSList * fl = fh->f;;
  while (fl) {
    f[i] = *(Forces *) fl->data;
    fl = fl->next;
    i++;
  }

  Motion m1;
  /* Prediction */
  m1.x0 = m0.x0 + dt/24.*(-9.*f[0].fx + 37.*f[1].fx - 59.*f[2].fx + 55*f[3].fx);
  m1.v0 = m0.v0 + dt/24.*(-9.*f[0].fv + 37.*f[1].fv - 59.*f[2].fv + 55*f[3].fv);

  /* Correction */
  Forces f4 = forces (t+dt, t, dt, m1.v0, m1.x0);
  m1.x0 = m0.x0 + dt/24.*(f[1].fx - 5.*f[2].fx + 19.*f[3].fx + 9.*f4.fx);
  m1.v0 = m0.v0 + dt/24.*(f[1].fv - 5.*f[2].fv + 19.*f[3].fv + 9.*f4.fv);

  return m1;
}

Motion motion_ABM4_time_integration (gdouble t, gdouble dt, Motion m0, ForcesHistory * fh, Forces f1)
{
  gint i = 0;
  Forces f[4];

  GSList * fl = fh->f;;
  while (fl) {
    f[i] = *(Forces *) fl->data;
    fl = fl->next;
    i++;
  }

  Motion m1;
  /* Prediction */
  for ( i = 0; i < 6; i++) {
    m1.x[i] = m0.x[i]
      + dt/24.*( -9.*f[0].v[i] + 37.*f[1].v[i] - 59.*f[2].v[i] + 55.*f[3].v[i]);
    m1.v[i] = m0.v[i] 
      + dt/24.*( -9.*(f[0].forces_m[i] + f[0].forces_l[i] + f[0].forces_h[i] +
		      f[0].forces_fk[i] + f[0].forces_r[i] + f[0].forces_ext[i])
		 + 37.*(f[1].forces_m[i] + f[1].forces_l[i] + f[1].forces_h[i] +
			f[1].forces_fk[i] + f[1].forces_r[i] + f[1].forces_ext[i])
		 - 59.*(f[2].forces_m[i] + f[2].forces_l[i] + f[2].forces_h[i] +
			f[2].forces_fk[i] + f[2].forces_r[i] + f[2].forces_ext[i])
		 + 55.*(f[3].forces_m[i] + f[3].forces_l[i] + f[3].forces_h[i] +
			f[3].forces_fk[i] + f[3].forces_r[i] + f[3].forces_ext[i]));
  }



  /* Correction */
  Forces f4 = forces_interpolation (t+dt, t, dt, f[3], f1, m1.x, m1.v);

  for ( i = 0; i < 6; i++) {
    m1.x[i] = m0.x[i]
      + dt/24.*( f[1].v[i] - 5.*f[2].v[i] + 19.*f[3].v[i] + 9.*f4.v[i]);
    
    m1.v[i] = m0.v[i]
      + dt/24.*( (f[1].forces_m[i] + f[1].forces_l[i] + f[1].forces_h[i] +
		  f[1].forces_fk[i] + f[1].forces_r[i] + f[1].forces_ext[i])
    		 - 5.*(f[2].forces_m[i] + f[2].forces_l[i] + f[2].forces_h[i] +
		       f[2].forces_fk[i] + f[2].forces_r[i] + f[2].forces_ext[i])
    		 + 19.*(f[3].forces_m[i] + f[3].forces_l[i] + f[3].forces_h[i] +
			f[3].forces_fk[i] + f[3].forces_r[i] + f[3].forces_ext[i])
    		 + 9.*(f4.forces_m[i] + f4.forces_l[i] + f4.forces_h[i] +
		       f4.forces_fk[i] + f4.forces_r[i] + f4.forces_ext[i]));
  }

  return m1;
}

/* Solution test of harmonic oscillator */
static gdouble m_x (gdouble t)
{
  gdouble xi = 2.;
  gdouble w = M_PI;

  gdouble phi = atan(2*xi*w/(w*w-1));

  return -phi/fabs(phi)*1./sqrt(pow(1-w*w, 2.) + pow(2*xi*w,2.))*cos(w*t+phi)
    + exp(-xi*t)*(exp(t*sqrt(xi*xi-1)) + exp(-t*sqrt(xi*xi-1)));
}

void test_rk ()
{
  gdouble t = 0;
  gdouble dt = 0.001;
  FILE * fp = fopen ("rk.tmp","w");
  Motion m;
  m.x0 = 2;
  m.v0 = -4.;

  for ( t = 0; t < 100; t+=dt) {
    m = RK4 (t, dt, m, forces (t, t, dt, m.v0, m.x0), forces (t+dt, t, dt, m.v0, m.x0));
    fprintf (fp, "%f %f %f \n", t+dt, m.x0, m_x(t+dt));
  }
  fprintf (stdout, "Error: %e \n", fabs(m.x0-m_x(t)));

  fclose (fp);
}

void test_abm ()
{
  gdouble t = 0;
  gdouble dt = 0.001;
  FILE * fp = fopen ("abm.tmp","w");
  Motion m;
  m.x0 = 2.;
  m.v0 = -4.;

  // Initialise with 3 steps of RK4
  gint i;
  ForcesHistory * fh = g_malloc (sizeof(ForcesHistory));
  
  for (i = 0; i < 3; i++) {
    add_force_to_history (fh, forces (t, t, dt, m.v0, m.x0));
    m = RK4 (t, dt, m, forces (t, t, dt, m.v0, m.x0), forces (t+dt, t, dt, m.v0, m.x0));
    t += dt;
  }

  // Adam-Bashford-Moulton
  for ( t = t; t < 100; t+=dt) {
    add_force_to_history (fh, forces (t, t, dt, m.v0, m.x0));
    m = AMB4 (t, dt, m, fh);

    fprintf (fp, "%f %f %f \n", t+dt, m.x0, /* m.v0 */m_x(t+dt));
  }
  fprintf (stdout, "Error: %e \n", fabs(m.x0-m_x(t)));

  fclose (fp);
  forces_history_destroy (fh);
}



/* void solve_equation_of_motion_new (Hull * hull, */
/* 				   Time * time, */
/* 				   Simulation * sim) */
/* { */
/*   gdouble t = time->t; */
/*   gdouble dt = time->dt; */

  

/*   // Test */
/*   Forces * f1 = g_malloc (sizeof(Forces)); */

/*   // Calculates forces on the hull at present time */
/*   /\*  *f1 = wet_hull_pressure_force_integration (sim->hull, *\/ */
/*   /\* 						   sim, *\/ */
/*   /\* 						   sim->fs->s->hz, *\/ */
/*   /\* 						   t, &sim->wp); *\/ */

/*   gdouble xi = 2; */
/*   gdouble omega = M_PI; */
/*   forces_set_to_zero (f1); */
/*   //  f1.forces_m[1] = -2.*xi*v - x + cos(omega*(t+dt)); */
/*   f1->forces_m[4] = -2.*xi*hull->m.u[4] - hull->m.x[4] + cos(omega*(t+dt)); */
/*   f1->u[4] = hull->m.u[4]; */


/*   // Initialise with 3 steps of RK4 */
/*   g_assert (hull->fh); */
/*   if (time->itime < 3) { */
/*     // Finds the latest forces in history records */
/*     Forces * f0 = g_slist_last (hull->fh->f)->data; */
/*     hull->m = motion_RK4_time_integration_new  (hull, time, hull->m, f0, f1); */
/*   } */
/*   else { */
/*     // Adam-Bashford-Moulton */
/*     hull->m = motion_ABM4_time_integration_new (hull, time, hull->m, hull->fh, f1); */
/*   } */

/*   // The integration method requires storage of the first derivatives */
/*   gint i; */
/*   for ( i = 0; i < 6; i++) { */
/*     f1->u[i] = hull->m.u[i]; */
/*     f1->x[i] = hull->m.x[i]; */
/*   } */
/*   add_force_to_history (hull->fh, *f1); */
/*   g_free (f1); */
/* } */

static void motion_motion (Hull * hull, Simulation * sim)
{
  Vector6 u;
  gint i, j;

  // Build mass matrix
  for ( i = 0; i < 6; i++) {
    for ( j = 0; j < 6; j++) {
      hull->M[i][j] = 0.;
    }
  }

  for ( i = 0; i < 3; i++)
    hull->M[i][i] = hull->mg;

  for ( i = 0; i < 3; i++) {
    for ( j = 0; j < 3; j++) {
      hull->M[i+3][j+3] = hull->Ig[i][j];
    }
  }

  // !!!! Here we add added_mass matrix instead of time-local
  // forces
  for ( i = 0; i < 6; i++) {
    for ( j = 0; j < 6; j++) {
      hull->M[i+3][j+3] += hull->A[i][j];
    }
  }

  // Build matrix of change of coordinates  
  update_rotation_matrix (&hull->m, hull->m.x, hull->m.u);

  // m omega x u
  gdouble tmp1[3];
  tmp1[0] = hull->mg*(u.x[4]*u.x[2]-u.x[5]*u.x[1]);
  tmp1[1] = hull->mg*(u.x[5]*u.x[0]-u.x[3]*u.x[2]);
  tmp1[2] = hull->mg*(u.x[3]*u.x[1]-u.x[4]*u.x[0]);

  // omega x Ig omega
  gdouble tmp2[3];
  tmp2[0] = tmp2[1] = tmp2[2] = 0.;
  
  for ( i = 0; i < 3; i++) {
    for ( j = 0; j < 3; j++) {
      tmp2[i] += hull->Ig[i][j]*u.x[3+j];
    }
  }

  tmp2[0] = u.x[4]*tmp2[2]-u.x[5]*tmp2[1];
  tmp2[1] = u.x[5]*tmp2[0]-u.x[3]*tmp2[2];
  tmp2[2] = u.x[3]*tmp2[1]-u.x[4]*tmp2[0];

  // Total forces
  Vector6 ft;
  Forces f;

  ft.x[0] = f.forces_m[0]+f.forces_l[0]+f.forces_h[0]
    +f.forces_fk[0]+f.forces_r[0];
  ft.x[1] = f.forces_m[1]+f.forces_l[1]+f.forces_h[1]
    +f.forces_fk[1]+f.forces_r[1];
  ft.x[2] = f.forces_m[2]+f.forces_l[2]+f.forces_h[2]
    +f.forces_fk[2]+f.forces_r[2];
  ft.x[3] = f.forces_m[3]+f.forces_l[3]+f.forces_h[3]
    +f.forces_fk[3]+f.forces_r[3];
  ft.x[4] = f.forces_m[4]+f.forces_l[4]+f.forces_h[4]
    +f.forces_fk[4]+f.forces_r[4];
  ft.x[5] = f.forces_m[5]+f.forces_l[5]+f.forces_h[5]
    +f.forces_fk[5]+f.forces_r[5];

  // Entrainement forces are functions of y(t)

  // (M+A)*du/dt= F + f(u(t))

  // k1.x = dt*B(u(t))*u(t)
  // k1.u = (M+A)^-1 * (F+f(u(t)))


  // Apply RK4 and ABM4 on that

  // Store entrainement forces at the end !!

}

/**
 * Boundary condition on the hull for the time-local flow
 * as stated page 31 of (Huang, 1997). This is for weak-scatter formulation (?).
 **/
gdouble hull_motion_time_local_bc_huang (SPPanel * spp, gint m, gint n, gpointer data)
{
  Hull * hull = (Hull *) data;
  Motion m1 = hull->m;

  GaussPoints * gp = spp->outer;
  gint ng = spp->sp->nouter;

  Vector N = g_array_index (gp->Ni, Vector, m + n*ng);
  Point P = g_array_index (gp->Pi, Point, m + n*ng);

  Vector delta;

  // delta = d/dt ( xi_t + xi_r x xp )
  // In the reference frame of the boat xp is constant (time derivative = 0)
  // + invariance of scalar-product when working on orthonormal bases
  delta.x = m1.v[0] + m1.v[4]*P.z - m1.v[5]*P.y;
  delta.y = m1.v[1] + m1.v[5]*P.x - m1.v[3]*P.z;
  delta.z = m1.v[2] + m1.v[3]*P.y - m1.v[4]*P.x;

  return vector_scalar_product (&delta, &N);
}

/**
 * Boundary condition on the hull for the time-local flow
 * as stated page 40 of (Kim, 2011). This is for the linear formulation.
 **/
/* gdouble hull_motion_time_local_bc_linear_kim (SPPanel * spp, gint m, gint n, gpointer data) */
/* { */
/*   Hull * hull = (Hull *) data; */
/*   Motion m1 = hull->m; */
/*   Vector U = hull->U; */

/*   Spline2D * sp = spp->sp; */
/*   GaussPoints * gp = spp->outer; */
/*   gint ng = spp->sp->nouter; */

/*   Vector Ni = g_array_index (gp->Ni, Vector, m + n*ng); */
/*   Point P = g_array_index (gp->Pi, Point, m + n*ng); */

/*   gdouble N[6], M[6]; */
/*   N[0] = Ni.x; */
/*   N[1] = Ni.y; */
/*   N[2] = Ni.z; */
/*   N[3] = P.y*Ni.z - P.z*Ni.y; */
/*   N[4] = P.z*Ni.x - P.x*Ni.z; */
/*   N[5] = P.x*Ni.y - P.y*Ni.x; */

/*   // m-terms from (Wu, 1991) */
/*   gint i, j; */
/*   size_t istart = gp->istart, jstart = gp->jstart; */
/*   gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m); */
/*   gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n); */
/*   gdouble Phi_x = 0., Phi_y = 0., Phi_z = 0.; */

/*   for ( i = 0; i < 6; i++) */
/*     M[i] = 0.; */

/*   if (sp->periodic) */
/*     istart -= (sp->k-1); */

/*   for ( i = 0; i < sp->k; i++) { */
/*     gdouble cu = gsl_matrix_get (Bu, i, 0); */
/*     gint ii = istart; */
/*     for ( j = 0; j < sp->k; j++) { */
/*       gint jj = (jstart+j); */
/*       gdouble cuv = cu*gsl_matrix_get (Bv, j, 0); */
      
/*       M[0] += coeff (sp, ii, jj, 25)*cuv; */
/*       M[1] += coeff (sp, ii, jj, 27)*cuv; */
/*       M[2] += coeff (sp, ii, jj, 29)*cuv; */

/*       // These are the values used as bc to find m1, m2 and m3 */
/*       // Alternatively, the values from the spline derivatives could be used. */
/*       Phi_x += coeff (sp, ii, jj, 24)*cuv; */
/*       Phi_y += coeff (sp, ii, jj, 26)*cuv; */
/*       Phi_z += coeff (sp, ii, jj, 28)*cuv; */
/*     } */
/*     istart++; */
/*   } */

/*   M[3] = M[2]*P.y - Ni.y*(Phi_z-U.z) - M[1]*P.z + Ni.z*(Phi_y-U.y); */
/*   M[4] = M[0]*P.z - Ni.z*(Phi_x-U.x) - M[2]*P.x + Ni.x*(Phi_z-U.z); */
/*   M[5] = M[1]*P.x - Ni.x*(Phi_y-U.y) - M[0]*P.y + Ni.y*(Phi_x-U.x); */

/*   // bc = sum xi_t * n + xi * m */
/*   gdouble val = 0.; */
/*   for ( i = 0; i < 6; i++) */
/*     val += m1.v[i]*N[i] + m1.x[i]*M[i]; */

/*   return val; */
/* } */

/*****************************************************/
/*                 CLEAN FORMULATION                 */
/*****************************************************/

/* Updates the first order rotation matrix according to x1 */
void motion_update_rotation_matrix_1 (Motion * m, gdouble x1[6])
{
  gint i, j;

  m->Rbi1.a[0][0] = 0.; m->Rbi1.a[0][1] = -x1[5]; m->Rbi1.a[0][2] = x1[4];
  m->Rbi1.a[1][0] = x1[5]; m->Rbi1.a[1][1] = 0.; m->Rbi1.a[1][2] = -x1[3];
  m->Rbi1.a[2][0] = -x1[4]; m->Rbi1.a[2][1] = x1[3]; m->Rbi1.a[2][2] = 0.;

  for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++)
      m->Rib1.a[j][i] = m->Rbi1.a[i][j];

  gdouble c1 = cos(x1[3]), c2 = cos(x1[4]), c3 = cos(x1[5]);
  gdouble s1 = sin(x1[3]), s2 = sin(x1[4]), s3 = sin(x1[5]);
  gdouble t2 = tan(x1[4]);

  // For small angles could be replaced by identity
  // This is complex and given by the euler angles this is the inverse of the 
  // Matrix defined by formula (2.25) of (Shao, 2010)
  // Transports angular momentum from boat to inertial frame of reference
  m->euler_r.a[0][0] = 1.; m->euler_r.a[0][1] = s1*t2; m->euler_r.a[0][2] = c1*t2;
  m->euler_r.a[1][0] = 0.; m->euler_r.a[1][1] = c1;    m->euler_r.a[1][2] = -s1;
  m->euler_r.a[2][0] = 0.; m->euler_r.a[2][1] = s1/c2; m->euler_r.a[2][2] = c1/c2;

  // Could be linearized for small angles
  m->euler_m.a[0][0] = c2*c3; m->euler_m.a[0][1] = s1*s2*c3-c1*s3; m->euler_m.a[0][2] = c1*s2*c3+s1*s3;
  m->euler_m.a[1][0] = c2*s3; m->euler_m.a[1][1] = s1*s2*s3+c1*c3; m->euler_m.a[1][2] = c1*s2*s3-s1*c3;
  m->euler_m.a[2][0] = -s2;   m->euler_m.a[2][1] = s1*c2;          m->euler_m.a[2][2] = c1*c2;

  m->t.x = x1[0];
  m->t.y = x1[1];
  m->t.z = x1[2];
}

void motion_update_rotation_matrix_2 (Motion * m,
				      gdouble x1[6], gdouble x2[6])
{
  gint i, j;

  m->Rbi2.a[0][0] = -0.5*(x1[4]*x1[4]+x1[5]*x1[5]);
  m->Rbi2.a[0][1] = -x2[5] + x1[3]*x1[4];
  m->Rbi2.a[0][2] = x2[4] + x1[3]*x1[5];
  m->Rbi2.a[1][0] = x2[5];
  m->Rbi2.a[1][1] = -0.5*(x1[3]*x1[3]+x1[5]*x1[5]);
  m->Rbi2.a[1][2] = -x2[3] + x1[4]*x1[5];
  m->Rbi2.a[2][0] = -x2[4];
  m->Rbi2.a[2][1] = x2[3];
  m->Rbi2.a[2][2] = -0.5*(x1[3]*x1[3]+x1[4]*x1[4]);

  for ( i = 0; i < 3; i++)
    for ( j = 0; j < 3; j++)
      m->Rib2.a[j][i] = m->Rbi2.a[i][j];
}

void hull_initialise_motion (Hull * h)
{
  //Hull * h = sim->hull;

  gint i, j;
  for ( i = 0; i < 6; i++) {
    h->m.x1[i] = h->m.x2[i] = 0.;
    h->m.u1[i] = h->m.u2[i] = 0.;
  }

  motion_update_rotation_matrix_1 (&h->m, h->m.x1);
  motion_update_rotation_matrix_2 (&h->m, h->m.x1, h->m.x2);


  //Build the mass matrix
  for ( i = 0; i < 3; i++) {
    for ( j = 0; j < 3; j++) {
      h->M2[i][j] = 0.;
      h->M2[i+3][j] = h->M2[i][j+3] = 0.;
      h->M2[i+3][j+3] = h->Ig[i][j];
    }
  }

  for ( i = 0; i < 3; i++)
    h->M2[i][i] = h->mg;

  /* if (h->mass_lu2 == NULL) { */
  /*   gsl_matrix * M = gsl_matrix_alloc (6,6); */
  /*   for ( i = 0; i < 6; i++ ) */
  /*     for ( j = 0; j < 6; j++) */
  /* 	gsl_matrix_set (M, i, j, h->M[j][i]); */
  /*   //sim->mass_lu2 = sim->lu_factorise (M); */
  /*   g_assert_not_reached (); */
  /* } */

  /* // Add added mass matrix */
  /* for ( i = 0; i < 6; i++) { */
  /*   for ( j = 0; j < 6; j++) { */
  /*     h->M1[i][j] = h->M2[i][j] + h->A[i][j]; */
  /*   } */
  /* } */

  // Fix added mass matrix for symmetry
  for ( i = 0; i < 6; i++)
    for ( j = 0; j < 6; j++)
      if ( (i+j)%2 == 0)
  	h->M1[i][j] = h->M2[i][j] + (h->A[i][j]+h->A[j][i])/2.;

  /* if (h->mass_lu1 == NULL) { */
  /*   gsl_matrix * M = gsl_matrix_alloc (6,6); */
  /*   for ( i = 0; i < 6; i++ ) */
  /*     for ( j = 0; j < 6; j++) */
  /* 	gsl_matrix_set (M, i, j, h->M1[j][i]); */
  /*   h->mass_lu1 = sim->lu_factorise (M); */
  /*   //g_assert_not_reached (); */
  /* } */

  // Allocates memory to store Forces History required by ABM4
  /* h->fh1 = g_malloc (sizeof(ForcesHistory)); */
  /* h->fh2 = g_malloc (sizeof(ForcesHistory)); */
  h->fh = g_malloc (sizeof(ForcesHistory));

  // Creates an "empty vector" and add forces to history
  //h->fh1->f = NULL; // List initially empty
  //h->fh2->f = NULL; // List initially empty
  h->fh->f = NULL; // List initially empty
  for ( i = 0; i < 4; i++) {
    Forces * f0 = g_malloc (sizeof(Forces));
    forces_set_to_zero (f0);
    h->fh->f = g_slist_append (h->fh->f, f0);
    /* f0 = g_malloc (sizeof(Forces)); */
    /* forces_set_to_zero (f0); */
    /* h->fh1->f = g_slist_append (h->fh1->f, f0); */
    /* f0 = g_malloc (sizeof(Forces)); */
    /* forces_set_to_zero (f0); */
    /* h->fh2->f = g_slist_append (h->fh2->f, f0); */
  }
}

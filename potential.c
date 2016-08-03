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
#include "linearproblem.h"
#include "patch.h"
#include "b4.h"
#include "boundaries.h"
#include "surface.h"
#include "hull.h"
#include "potential.h"

typedef struct {
  gint i;
  Point p;
} RankineParams;

void panel_rankine_source_set (Panel * p, gpointer data)
{
  RankineParams * par = (RankineParams *) data;
  gdouble dx, dy, dz;
  dx = PANEL_VAL(p,0) - par->p.x;
  dy = PANEL_VAL(p,1) - par->p.y;
  dy = PANEL_VAL(p,2) - par->p.z;

  p->var[par->i] = 2.*M_PI/sqrt(dx*dx+dy*dy+dz*dz);
  /* fprintf (stdout,"%i %g \n",par->i , p->var[par->i]); */
}

void potential_problem_build (Hull * hull, FreeSurface * fs, Bathymetry * bathy)
{
  LinearProblem * lp = linear_problem_new ();

  //linear_problem_init_size (lp, gint N, gint M);

  gint i, j, k;
  FILE * fp = fopen("potential.tmp","w");

  RankineParams params;
  params.i = 12;
  gint N = 30;
  gint M = 60;

  for ( i = 0; i <= M; i++) {
    for ( j = 0; j <= N; j++) {
      params.p.x =  0;
      params.p.y = -50 + i/M*100;
      params.p.z = -10 + 15*j/N;

      for ( k = 0; k < hull->patches->len; k++) {
	Patch * patch = g_ptr_array_index (hull->patches, k);
	patch_forall_panels (patch, panel_rankine_source_set, &params);
      }

      /* for ( k = 0; k < fs->s->patches->len; k++) { */
      /* 	Patch * patch = g_ptr_array_index (fs->s->patches, k); */
      /* 	patch_forall_panels (patch, panel_rankine_source_set, &params); */
      /* } */

      /* for ( k = 0; k < bathy->s->patches->len; k++) { */
      /* 	Patch * patch = g_ptr_array_index (bathy->s->patches, k); */
      /* 	patch_forall_panels (patch, panel_rankine_source_set, &params); */
      /* } */

      gdouble sum = 0.;
      for ( k = 0; k < hull->patches->len; k++) {
	Patch * patch = g_ptr_array_index (hull->patches, k);
	//sum += b4_patch_adaptive_panel_scalar_integral (patch, fs->s->hz, params.i);
	sum += b4_patch_fast_panel_scalar_integral (patch, fs->s->hz, params.i);
      }

      /* for ( k = 0; k < fs->s->patches->len; k++) { */
      /* 	Patch * patch = g_ptr_array_index (fs->s->patches, k); */
      /* 	sum += patch_adaptive_panel_scalar_integral (patch, NULL, params.i); */
      /* } */

      /* for ( k = 0; k < bathy->s->patches->len; k++) { */
      /* 	Patch * patch = g_ptr_array_index (bathy->s->patches, k); */
      /* 	sum += patch_adaptive_panel_scalar_integral (patch, params.i); */
      /* } */
      fprintf(stderr, "%i %i %e\n",i,j, sum);
      fprintf (fp, "%g %g %g %g \n", params.p.x, params.p.y, params.p.z, sum);
    }
    fprintf (fp, "\n");
  }

  

  
  /* fprintf (fp, "%g %g %g %g \n", params.px) */
  fclose (fp);

  linear_problem_destroy (lp);
}

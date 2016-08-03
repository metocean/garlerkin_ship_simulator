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

/* LinearProblem methods */

LinearProblem * linear_problem_new ()
{
  LinearProblem * new = g_malloc (sizeof(LinearProblem));
  new->A = g_ptr_array_new ();
  new->Ai = g_ptr_array_new ();
  new->A2 = g_ptr_array_new ();
  new->Ai2 = g_ptr_array_new ();
  new->rhs = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->rhs2 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->lhs = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->lhs2 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->gxi0 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->gxiN = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->geta0 = g_array_new (FALSE, FALSE, sizeof(gdouble));
  new->getaN = g_array_new (FALSE, FALSE, sizeof(gdouble));
  return new;
}

void linear_problem_destroy (LinearProblem * lp)
{
  g_assert (lp != NULL);
  g_assert (lp->A->len == lp->Ai->len);
  gint i;
  
  if (lp->A != NULL) {
    for ( i = 0; i < lp->A->len; i++)
      g_array_free ((GArray *) g_ptr_array_index(lp->A, i), TRUE);
    g_ptr_array_free (lp->A, TRUE);
  }

  if (lp->Ai != NULL) {
    for ( i = 0; i < lp->Ai->len; i++)
      g_array_free ((GArray *) g_ptr_array_index(lp->Ai, i), TRUE);
    g_ptr_array_free (lp->Ai, TRUE);
  }

  if (lp->A2 != NULL) {
    for ( i = 0; i < lp->A2->len; i++)
      g_array_free ((GArray *) g_ptr_array_index(lp->A2, i), TRUE);
    g_ptr_array_free (lp->A2, TRUE);
  }

  if (lp->Ai2 != NULL) {
    for ( i = 0; i < lp->Ai2->len; i++)
      g_array_free ((GArray *) g_ptr_array_index(lp->Ai2, i), TRUE);
    g_ptr_array_free (lp->Ai2, TRUE);
  }
  
  g_array_free (lp->rhs, TRUE);
  g_array_free (lp->rhs2, TRUE);
  g_array_free (lp->lhs, TRUE);
  g_array_free (lp->lhs2, TRUE);
  g_array_free (lp->gxi0, TRUE);
  g_array_free (lp->gxiN, TRUE);
  g_array_free (lp->geta0, TRUE);
  g_array_free (lp->getaN, TRUE);

  g_free (lp);
}

void linear_problem_init_size (LinearProblem * lp, gint N, gint M)
{
  gint i;

  for ( i = 0; i < N*M; i++) {
    GArray * new = g_array_new (FALSE, TRUE, sizeof(gdouble));
    g_ptr_array_add (lp->A, new);
    new = g_array_new (FALSE, TRUE, sizeof(gdouble));
    g_ptr_array_add (lp->A2, new);
    new = g_array_new (FALSE, TRUE, sizeof(gint));
    g_ptr_array_add (lp->Ai, new);
    new = g_array_new (FALSE, TRUE, sizeof(gint));
    g_ptr_array_add (lp->Ai2, new);
    gdouble val = 0.;
    g_array_append_val (lp->rhs, val);
    g_array_append_val (lp->rhs2, val);
    /* complex double val0 = 0.; */
    /* g_array_append_val (lp->rhsi, val0); */
    val = 1.;
    g_array_append_val (lp->lhs, val);
    g_array_append_val (lp->lhs2, val);
  }
  for ( i = 0; i < N; i++) {
    gdouble val = 0.;
    g_array_append_val (lp->gxi0, val);
    g_array_append_val (lp->gxiN, val);
  }
  for ( i = 0; i < M; i++) {
    gdouble val = 0.;
    g_array_append_val (lp->geta0, val);
    g_array_append_val (lp->getaN, val);
  }
}

void linear_problem_add_stencil (LinearProblem * lp, gint i,
				 gint j, gdouble coeff)
{
  g_assert (lp != NULL);
  g_assert (g_ptr_array_index(lp->A,i) != NULL);
  g_assert (g_ptr_array_index(lp->Ai,i) != NULL);
  /* if ( i > 300 || j > 300) */
  /*   fprintf(stdout, "Phi: %i %i \n", i, j); */
  gint k;
  // Finds the row of the linear problem
  GArray * coeffs = g_ptr_array_index(lp->A,i);
  GArray * index = g_ptr_array_index(lp->Ai,i);

  // Finds whether there already is a contribution for the same column
  for ( k = 0; k < index->len; k++) {
    if ( g_array_index (index, gint, k) == j ) {
      g_array_index (coeffs, gdouble, k) += coeff;
      return;
    }
  }

  g_array_append_val (coeffs, coeff);
  g_array_append_val (index, j);
}

/* Solve the linear problem using the Successive Over-Relaxation method */
gboolean linear_problem_solve_using_SOR (LinearProblem * lp, gdouble tolerance,
					 gboolean verbose)
{
  gdouble err = 1;
  gdouble w = 1.;
  gint i, j, niter = 0;
  GArray * x0 = g_array_new (FALSE, FALSE, sizeof(gdouble));

  for ( i = 0; i < lp->lhs->len; i++)
    g_array_append_val (x0, g_array_index(lp->rhs ,gdouble, i));
  
  while (err > tolerance) {
    err = 0;
    for ( i = 0; i < lp->A->len; i++) {
      gdouble sum = g_array_index (lp->rhs, gdouble, i);
      GArray * coeff = (GArray * ) g_ptr_array_index (lp->A, i);
      GArray * index = g_ptr_array_index (lp->Ai, i);
      gdouble aii = 0.;
      
      for (j = 0; j < index->len ; j++) {
    	if (g_array_index(index ,gint ,j) == i)
    	  aii = g_array_index(coeff ,gdouble , j);
    	else if (g_array_index(index ,gint ,j) > i)
    	  sum -= g_array_index(coeff ,gdouble , j)
    	    *g_array_index(lp->lhs ,gdouble, g_array_index(index ,gint ,j));
    	else
    	  sum -= g_array_index(coeff ,gdouble , j)
    	    *g_array_index(x0 ,gdouble, g_array_index(index ,gint ,j));
      }
      g_array_index (x0, gdouble, i) = (1. - w)*g_array_index(lp->lhs ,gdouble, i)
    	+ w/aii*sum;
    }
    
    for ( i = 0; i < lp->A->len; i++) {
      if (fabs(g_array_index(lp->lhs ,gdouble, i)
    	       - g_array_index(x0 ,gdouble, i)) > err)
    	err = fabs(g_array_index(lp->lhs ,gdouble, i)
    		   - g_array_index(x0 ,gdouble, i));
    }

    for ( i = 0; i < lp->lhs->len; i++)
      g_array_index (lp->lhs, gdouble, i) = g_array_index (x0, gdouble, i);
    if (verbose)
      fprintf(stdout, "SOR error %e after %i iterations \n", err, niter);
    niter++;
  }

  g_array_free (x0, FALSE);

  /* fprintf(stdout, "SOR error %e after %i iterations \n", err, niter); */

  return TRUE;
}

gboolean linear_problem_SOR_iteration (LinearProblem * lp)
{
  gdouble err = 0;
  gdouble w = 1.;
  gint i, j;
  GArray * x0 = g_array_new (FALSE, FALSE, sizeof(gdouble));

  for ( i = 0; i < lp->lhs->len; i++)
    g_array_append_val (x0, g_array_index(lp->rhs ,gdouble, i));
  
  for ( i = 0; i < lp->A->len; i++) {
    gdouble sum = g_array_index (lp->rhs, gdouble, i);
    GArray * coeff = (GArray * ) g_ptr_array_index (lp->A, i);
    GArray * index = g_ptr_array_index (lp->Ai, i);
    gdouble aii = 0.;
    
    for (j = 0; j < index->len ; j++) {
      if (g_array_index(index ,gint ,j) == i)
	aii = g_array_index(coeff ,gdouble , j);
      else if (g_array_index(index ,gint ,j) > i)
	sum -= g_array_index(coeff ,gdouble , j)
	  *g_array_index(lp->lhs ,gdouble, g_array_index(index ,gint ,j));
      else
	sum -= g_array_index(coeff ,gdouble , j)
	  *g_array_index(x0 ,gdouble, g_array_index(index ,gint ,j));
    }
    g_array_index (x0, gdouble, i) = (1. - w)*g_array_index(lp->lhs ,gdouble, i)
      + w/aii*sum;
  }
  
  for ( i = 0; i < lp->A->len; i++) {
    /* gdouble error =  */
    if (fabs(g_array_index(lp->lhs ,gdouble, i)
	     - g_array_index(x0 ,gdouble, i)) > err)
      err = fabs(g_array_index(lp->lhs ,gdouble, i)
		 - g_array_index(x0 ,gdouble, i));
  }
  
  for ( i = 0; i < lp->lhs->len; i++)
    g_array_index (lp->lhs, gdouble, i) = g_array_index (x0, gdouble, i);

  g_array_free (x0, FALSE);

  /* fprintf(stdout, "SOR error %g\n", err); */

  return TRUE;
}

/* Normal from polynomial expension */
/* static Vector normal (S2P * s2p, gint k, gdouble u, gdouble v) */
/* { */
/*   Vector N; */
/*   N.x = N.y = N.z = 0.; */
/*   gint m, n; */
/*   gint M = 2*k-2; */

/*   for ( m = 0; m <= 2*(M-1); m++) { */
/*     for ( n = 0; n <= m; n++) { */
/*       N = vector_sum (N, vector_times_constant (g_array_index (s2p->n, Vector, n + (m-n)*(6*(M-1)+1)), */
/* 						pow(u-s2p->ue,n)*pow(v-s2p->ve,m-n))); */
/*     } */
/*   } */

/*   return N; */
/* } */

gint spline_numbering (GSList * list)
{
  GSList * plist = list;
  gint SIZE = 0;

  // Here the numbering is tricky as we need to make sure
  // that the free-surface is numbered properly is order to
  // ensure periodicity
  // It would be easier if by convention the free-surface came first
  // in the list of patches

  while (plist != NULL) {
    Spline2D * sp = plist->data;
    sp->istart = SIZE;
    if (sp->periodic) {
      g_assert (sp->NUT != 0);
      SIZE += sp->NUT*sp->NV;
    }
    else {
      SIZE += sp->NU*sp->NV;
    }
    plist = plist->next;
  }

  return SIZE;
}

gint problem_size (GSList * list)
{
  GSList * plist = list;
  gint SIZE = 0;

  while (plist != NULL) {
    Spline2D * sp = plist->data;
    if (sp->periodic) {
      g_assert (sp->NUT != 0);
      SIZE += sp->NUT*sp->NV;
    }
    else {
      SIZE += sp->NU*sp->NV;
    }
    plist = plist->next;
  }

  return SIZE;
}

/* void compute_self_influence_coefficients (Spline2D * sp) */
/* { */
/*   gint i, m, n; */

/*   // Loop over the panels of the patch */
/*   for ( i = 0; i < sp->M*sp->N; i++) { */
/*     SPPanel * spp = g_ptr_array_index (sp->panels, i); */

/*     if (spp->sic != NULL) */
/*       return; */

/*     spp->sic = g_ptr_array_new (); */
     	
/*     // Gauss outer-integration */
/*     GaussPoints * gp = spp->outer; */
/*     gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule */
    	  
/*     for ( n = 0; n < ng; n++) { */
/*       for ( m = 0; m < ng; m++) { */
/* 	InfluenceCoeffs * ic =  sppanel_self_influence_coeff (spp, m, n); */
/* 	g_ptr_array_add (spp->sic, ic); */
/*       } */
/*     } */
/*   } */
/* } */

/**
 * Compute and stores the influence coefficients thats represent the influence of
 * the panel @spp on each of the outer gauss nodes of the patch sp.
 **/
/* void compute_influence_coefficients (SPPanel * spp, Spline2D * sp) */
/* { */
/*   gint ii, m, n; */

/*   /\* GPtrArray * level1 = g_ptr_array_new (); *\/ */
/*   /\* GPtrArray * level1 = g_ptr_array_new (); *\/ */
/*   /\* GPtrArray * level1 = g_ptr_array_new (); *\/ */
/*   /\* GPtrArray * level1 = g_ptr_array_new (); *\/ */
/*   /\* GPtrArray * level1 = g_ptr_array_new (); *\/ */

/*   for ( ii = 0; ii < sp->panels->len; ii++ ) { */
/*     SPPanel * panel = g_ptr_array_index (sp->panels, ii); */
/*     GPtrArray * sic = g_ptr_array_new (); */
/*     GaussPoints * gp = panel->outer; */
/*     gint ng = gp->ui->len; */
    
/*     // Loop over the outer Gauss-Points */
/*     for ( m = 0; m < ng; m++ ) { */
/*       for ( n = 0; n < ng; n++ ) { */
	
/*       } */
/*     } */


/*     panel->sic = sic; */
/*   } */


/*   /\* gint i, m, n; *\/ */

/*   /\* // Loop over the panels of the patch *\/ */
/*   /\* for ( i = 0; i < sp->M*sp->N; i++) { *\/ */
/*   /\*   SPPanel * spp = g_ptr_array_index (sp->panels, i); *\/ */

/*   /\*   if (spp->sic != NULL) *\/ */
/*   /\*     return; *\/ */

/*   /\*   spp->sic = g_ptr_array_new (); *\/ */
     	
/*   /\*   // Gauss outer-integration *\/ */
/*   /\*   GaussPoints * gp = spp->outer; *\/ */
/*   /\*   gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule *\/ */
    	  
/*   /\*   for ( n = 0; n < ng; n++) { *\/ */
/*   /\*     for ( m = 0; m < ng; m++) { *\/ */
/*   /\* 	InfluenceCoeffs * ic =  sppanel_self_influence_coeff (spp, m, n); *\/ */
/*   /\* 	g_ptr_array_add (spp->sic, ic); *\/ */
/*   /\*     } *\/ */
/*   /\*   } *\/ */
/*   /\* } *\/ */
/* } */

gboolean linear_problem_solve_using_gauss (LinearProblem * lp)
{
  int n = lp->lhs->len;
  gdouble a[lp->lhs->len+1][lp->lhs->len];

  // Copy Sparse problem to full matrix
  gint i, j, k, maxrow;
  double tmp;

  for (i = 0; i < lp->lhs->len; i++)
    for (j = 0; j < lp->lhs->len+1; j++)
      a[i][j] = 0.;

  for (i = 0; i < lp->A->len; i++) {
    GArray * coeffs = g_ptr_array_index (lp->A, i);
    GArray * index = g_ptr_array_index (lp->Ai, i);

    for (j = 0; j < coeffs->len; j++) {
      a[i][g_array_index (index, gint, j)] = g_array_index (coeffs, gdouble, j);
    }

    a[n][i] = g_array_index (lp->rhs, gdouble, i);
  }

 
  for ( i = 0; i < n; i++) {

      /* Find the row with the largest first value */
      maxrow = i;
      for ( j = i+1; j < n; j++) {
         if (fabs(a[i][j]) > fabs(a[i][maxrow]))
            maxrow = j;
      }

      /* Swap the maxrow and ith row */
      for ( k = i; k < n+1; k++) {
  	tmp = a[k][i];
  	a[k][i] = a[k][maxrow];
  	a[k][maxrow] = tmp;
      }

      /* Singular matrix? */
      if (fabs(a[i][i]) < 1e-12)
         return(FALSE);

      /* Eliminate the ith element of the jth row */
      for ( j = i+1; j < n; j++) {
         for ( k = n; k >= i;k--) {
            a[k][j] -= a[k][i] * a[i][j] / a[i][i];
         }
      }
  }

   /* Do the back substitution */
  for ( j = n-1; j >= 0; j--) {
      tmp = 0;
      for ( k = j+1; k < n; k++)
  	tmp += a[k][j] * g_array_index (lp->lhs, gdouble, k);
      g_array_index (lp->lhs, gdouble, j) = (a[n][j] - tmp) / a[j][j];
  }
  
  return TRUE;
}

/* void copy_solution_to_patches (GSList * list, LinearProblem * lp) */
/* { */
/*   gint i, j; */
/*   gint istart = 0; */
/*   GSList * plist = list; */

/*   // Go over the list of patches */
/*   while (plist != NULL) { */
/*     Spline2D * sp = plist->data; */
/*     g_assert (sp != NULL); */
    
/*     // Loop over the b-splines */
/*     for ( i = 0; i < gsl_bspline_ncoeffs (sp->w_u); i++) { */
/*       for ( j = 0; j < gsl_bspline_ncoeffs (sp->w_v); j++) { */
/* 	coeff_assign (sp, i, j, 3, */
/* 		      g_array_index (lp->lhs, gdouble, istart + i + j*gsl_bspline_ncoeffs (sp->w_u))); */
	
/*       } */
/*     } */
    
/*     // shift index for next patch */
/*     istart += gsl_bspline_ncoeffs (sp->w_u)*gsl_bspline_ncoeffs (sp->w_v); */
/*     plist = plist->next; */
/*   } */
/* } */

/************************************************************************/

/* BoundaryProblem methods */
BoundaryProblem * boundary_problem_new (gint N, gint M)
{
  BoundaryProblem * new = g_malloc (sizeof(BoundaryProblem));

  new->A = gsl_matrix_alloc (N, M);
  new->rhs = gsl_vector_alloc (N);

  return new;
}

void boundary_problem_destroy (BoundaryProblem * bp)
{
  g_assert (bp != NULL);
  if (bp->A)
    gsl_matrix_free (bp->A);
  gsl_vector_free (bp->rhs);
  g_free (bp);
}

BoundarySubProblem * boundary_subproblem_new (gint N, gint M)
{
  BoundarySubProblem * new = g_malloc (sizeof(BoundarySubProblem));

  new->neumann = boundary_problem_new (N, M);
  new->dirichlet = boundary_problem_new (N, M);

  return new;
}

void boundary_subproblem_destroy (BoundarySubProblem * b)
{
  boundary_problem_destroy (b->neumann);
  boundary_problem_destroy (b->dirichlet);
  g_free (b);
}

/* BoundarySubProblem *  build_boundary_subproblem (Spline2D * sp1, Spline2D * sp2) */
/* { */
/*   gint i, j, m, n, ii, a, b, mu, nu; */
/*   size_t ustart, vstart; */
/*   gdouble M_PI2 = 2.*M_PI; */
/*   gint N1 = gsl_bspline_ncoeffs (sp1->w_u)*gsl_bspline_ncoeffs (sp1->w_v); */
/*   gint N2 = gsl_bspline_ncoeffs (sp2->w_u)*gsl_bspline_ncoeffs (sp2->w_v); */

/*   // Initializes the coefficients to 0 */
/*   gdouble A[N1][N2]; */
/*   gdouble B[N1][N2]; */

/*   for ( i = 0; i < N1; i++) */
/*     for ( j = 0; j < N2; j++) { */
/*       A[i][j] = 0.; */
/*       B[i][j] = 0.; */
/*     } */
/*   // !! Change limit of the stack : ulimit -s 32768 */
/*   // !! To see limit : ulimit -s */

/*   g_assert (sp1 != NULL && sp2 != NULL); */
    
/*   // Loop over the panels of the patch */
/*   for ( ii = 0; ii < sp1->M*sp1->N; ii++) { */
/*     SPPanel * spp = g_ptr_array_index (sp1->panels, ii); */
/*     g_assert (spp != NULL); */
     	
/*     // Gauss outer-integration */
/*     GaussPoints * gp = spp->outer; */
/*     gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule */
/*     for ( m = 0; m < ng; m++) { */
/*       gdouble um = g_array_index (gp->ui, gdouble, m); */
/*       gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m); */
/*       ustart = g_array_index (gp->ustart, gint, m); */
	  
/*       for ( n = 0; n < ng; n++) { */
/* 	gdouble vn = g_array_index (gp->vj, gdouble, n);	     */
/* 	gdouble wmn = M_PI2*g_array_index (gp->wJij, gdouble, m+n*ng); */
/* 	gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n); */
/* 	vstart = g_array_index (gp->vstart, gint, n); */
/* 	Point p = g_array_index (gp->Pi, Point, m+n*ng); */

/* 	// Loop over the splines whose support is included in the panel */
/* 	for ( i = ustart; i < ustart + sp1->k; i++) { */
/* 	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0); */
/* 	  for ( j = vstart; j < vstart + sp1->k; j++) { */
/* 	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0); */
/* 	    gint indexi = i + j*gsl_bspline_ncoeffs (sp1->w_u); */
	      
/* 	    // First term of equation (2*pi*phi) */
/* 	    if (sp1 == sp2) */
/* 	      for ( a = 0; a < sp1->k; a++) { */
/* 		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0); */
/* 		for ( b = 0; b < sp1->k; b++) { */
/* 		  A[indexi][(ustart+a) + (vstart+b)*gsl_bspline_ncoeffs (sp1->w_u)] +=  wmnija*gsl_matrix_get (Bv, b, 0); */
/* 		} */
/* 	      } */
/* 	  } */
/* 	}	     */
	    
/* 	wmn = g_array_index (gp->wJij, gdouble, m+n*ng); */
/* 	// Second term */
/* 	// -> Here it would come */
/* 	// GPtrArray * ics = sppanel_influence_coefficients (panel, sp2) ?? */
/* 	// Loop over all the panels of the patch */
/* 	for ( mu = 0; mu < sp2->M; mu++) { */
/* 	  for ( nu = 0; nu < sp2->N; nu++) { */
/* 	    SPPanel * panel = g_ptr_array_index (sp2->panels, mu + nu*sp2->M); */
	     
/* 	    // Spline coefficients */
/* 	    InfluenceCoeffs * ic = NULL; */

/* 	    if ( spp == panel ) */
/* 	      ic = lachat_watson_self_influence_coefficients (panel, um, vn, p); */
/* 	    else */
/* 	      ic = spline_near_field_influence_coeff (panel, p); */

/* 	    // Loop over the splines whose support is included in the panel */
/* 	    for ( i = ustart; i < ustart + sp1->k; i++) { */
/* 	      gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0); */
/* 	      for ( j = vstart; j < vstart + sp1->k; j++) { */
/* 	    	gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0); */
/* 	    	gint indexi = i + j*gsl_bspline_ncoeffs (sp1->w_u); */
		      
/* 		for ( a = 0; a < sp2->k; a++) { */
/* 		  for ( b = 0; b < sp2->k; b++) { */
/* 		    gint indexa = (mu+a) + (nu+b)*gsl_bspline_ncoeffs (sp2->w_u); */

/* 		    A[indexi][indexa] -= wmnij*gsl_matrix_get (ic->phi, b, a); */
/* 		    B[indexi][indexa] += wmnij*gsl_matrix_get (ic->psi, b, a); */
		    
/* 		  } */
/* 		} */
	
/* 	      } */
/* 	    } */

/* 	    influencecoeffs_destroy (ic); */
/* 	  } */
/* 	} */
	    	
/* 	// end 2nd term loop     */
/*       } */
/*     } // End Gauss integration */
      
/*   } // End loop over panels */

/*   // Copy the problem back to the old linear problem structure */
/*   BoundarySubProblem * bsp = boundary_subproblem_new (N1, N2); */

/*   bsp->neumann->istart = bsp->dirichlet->istart = sp2->istart; */
/*   bsp->neumann->jstart = bsp->dirichlet->jstart = sp1->istart; */

/*   /\* FILE * fp = fopen ("check2.tmp","w"); *\/ */
  
/*   for ( j = 0; j < N1; j++) { */
/*     for ( i = 0; i < N2; i++) { */
/*       /\* fprintf (fp, "%i %i %e %e \n", i, j, A[j][i], B[j][i]); *\/ */
/*       gsl_matrix_set (bsp->neumann->A, j, i, A[j][i]); */
/*       gsl_matrix_set (bsp->dirichlet->A, j, i, B[j][i]); */
/*     } */
/*     /\* fprintf (fp, "\n"); *\/ */
/*   } */
/*   /\* fclose (fp); *\/ */

/*   return bsp; */
/* } */

BoundarySubProblem *  build_boundary_subproblem_galerkin_old (Spline2D * sp1, Spline2D * sp2,
							  SelfInfluenceFunc self_influence_coeffs)
{
  gint i, j, m, n, ii, a, b, mu, nu;
  size_t ustart, vstart;
  gdouble M_PI2 = 2.*M_PI;
  gint N1 = sp1->NU*sp1->NV;
  gint N2 = sp2->NU*sp2->NV;
  gint k1 = sp1->k, k2 = sp2->k;
  gint ng = sp1->nouter;

  // Initializes the coefficients to 0
  gdouble A[N1][N2];
  gdouble B[N1][N2];

  for ( i = 0; i < N1; i++)
    for ( j = 0; j < N2; j++) {
      A[i][j] = 0.;
      B[i][j] = 0.;
    }
  // !! Change limit of the stack : ulimit -s 32768
  // !! To see limit : ulimit -s

  g_assert (sp1 != NULL && sp2 != NULL);

  gint NU1 = sp1->NU;
  gint NU2 = sp2->NU;

  // Loop over the panels of the patch
  if (sp1 == sp2) {
    for ( ii = 0; ii < sp1->M*sp1->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp1->panels, ii);
      g_assert (spp != NULL);

      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      /* ustart = sp1->periodic ? ( gp->istart - k1 + 1) : gp->istart; */
      ustart = sp1->periodic ? ( gp->istart - k1 + 1) : gp->istart;
      gint NUT = sp1->periodic ? sp1->NUT: 100000000;
      // Eventual shift for index for periodic free-surface
      vstart = gp->jstart;

      for ( m = 0; m < ng; m++) {
  	gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
  	for ( n = 0; n < ng; n++) {
  	  gdouble wmn = M_PI2*g_array_index (gp->wJij, gdouble, m+n*ng);
  	  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);
	
  	  // Loop over the splines whose support is included in the panel
  	  for ( i = 0; i < k1; i++) {
  	    gdouble wmni = wmn*gsl_matrix_get (Bu, i, 0);
	    gint itmp = (ustart + i)%NU1;
  	    for ( j = 0; j < k1; j++) {
  	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j, 0);
  	      gint indexi =  itmp + (vstart+j)*NU1;
	      
  	      // First term of equation (2*pi*phi)
  	      for ( a = 0; a < k1; a++) {
  		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		gint atmp = (ustart+a)%NU1;
  		for ( b = 0; b < k1; b++) {
		  A[indexi][atmp + (vstart+b)*NU1] +=  wmnija*gsl_matrix_get (Bv, b, 0);
  		}
  	      }
  	    }
  	  }
  	}
      }
    }
  }


  gdouble * psi = g_malloc (k2*k2*sizeof(gdouble));
  gdouble * phi = g_malloc (k2*k2*sizeof(gdouble));

  for ( mu = 0; mu < sp2->M; mu++) {
    for ( nu = 0; nu < sp2->N; nu++) {
      SPPanel * panel = g_ptr_array_index (sp2->panels, mu + nu*sp2->M);

      GPCell * panel_tree = gpcell_tree_new (panel);
      for ( ii = 0; ii < sp1->M*sp1->N; ii++) {
	SPPanel * spp = g_ptr_array_index (sp1->panels, ii);
   	
	// Gauss outer-integration
	GaussPoints * gp = spp->outer;
	ustart = gp->istart;
	vstart = gp->jstart;
	
	gint pustart = sp1->periodic ? (ustart - k1 + 1) : ustart; // Eventual shift for index for periodic free-surface

	// Temporary storage
	gdouble tmpBv[k1][ng];
 	for ( i = 0; i < ng; i++ ) {
	  gsl_matrix * Bv = g_ptr_array_index (gp->Bv, i);
	  for ( j = 0; j < k1; j++) {
	    tmpBv[j][i] = gsl_matrix_get (Bv, j, 0);
	  }
	}

	for ( m = 0; m < ng; m++) {
	  gdouble um = g_array_index (gp->ui, gdouble, m);
	  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
	  for ( n = 0; n < ng; n++) {
	    gdouble vn = g_array_index (gp->vj, gdouble, n);	    
	    gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	    Point p = g_array_index (gp->Pi, Point, m+n*ng);

	     
	    // Spline coefficients
	    gint a;
	    for ( a = 0; a < k2*k2; a++)
	        psi[a] = phi[a] = 0.;

	    if ( spp == panel )
	       self_influence_coeffs (panel, um, vn, p, psi, phi);
	    else {
	      
	      spline_near_field_influence_coeff_recursive (panel, panel_tree, p, psi, phi, FALSE);
	    }

	    // Loop over the splines whose support is included in the panel
	    for ( i = 0; i < k1; i++) {
	      gdouble wmni = wmn*gsl_matrix_get (Bu, i, 0);
	      gint itmp = (pustart + i) % NU1;
	      for ( j = 0; j < k1; j++) {
	    	gdouble wmnij = wmni*tmpBv[j][n];
	    	gint indexi = itmp + (vstart+j)*NU1;

		for ( a = 0; a < k2; a++) {
		  gint atmp = (mu + a) % NU2;
		  for ( b = 0; b < k2; b++) {
		    gint indexa = atmp + (nu+b)*NU2;

		    A[indexi][indexa] -= wmnij*phi[a + b*k2];
		    B[indexi][indexa] += wmnij*psi[a + b*k2];
		    
      		  }
		}

	      }
	    }

	  }
	}
	    	   
      }
      
      gpcell_tree_destroy (panel_tree);
    }
      
  }
  g_free (psi);
  g_free (phi);

  // Copy the problem back to the old linear problem structure
  BoundarySubProblem * bsp = boundary_subproblem_new (N1, N2);

  bsp->neumann->istart = bsp->dirichlet->istart = sp2->istart;
  bsp->neumann->jstart = bsp->dirichlet->jstart = sp1->istart;

  /* FILE * fp = fopen ("check2.tmp","w"); */
  
  for ( j = 0; j < N1; j++) {
    for ( i = 0; i < N2; i++) {
      /* fprintf (fp, "%i %i %e %e \n", i, j, A[j][i], B[j][i]); */
      gsl_matrix_set (bsp->neumann->A, j, i, A[j][i]);
      gsl_matrix_set (bsp->dirichlet->A, j, i, B[j][i]);
    }
    /* fprintf (fp, "\n"); */
  }
  /* fclose (fp); */

  return bsp;
}

BoundarySubProblem *  build_boundary_subproblem_galerkin (Spline2D * sp1, Spline2D * sp2,
							  SelfInfluenceFunc self_influence_coeffs)
{
  //  fprintf (stderr, "Start \n");
  
  gint i, j, m, n, ii, a, b, mu, nu;
  size_t ustart, vstart;
  gdouble M_PI2 = 2.*M_PI;
  gint N1 = sp1->periodic ? sp1->NUT*sp1->NV : sp1->NU*sp1->NV;
  gint N2 = sp2->periodic ? sp2->NUT*sp2->NV : sp2->NU*sp2->NV;
  gint k1 = sp1->k, k2 = sp2->k;
  gint ng = sp1->nouter;
  //fprintf (stderr, "size: %i %i\n",N1,N2);
  //g_assert_not_reached ();
  // Initializes the coefficients to 0
  gdouble A[N1][N2];
  gdouble B[N1][N2];

  for ( i = 0; i < N1; i++)
    for ( j = 0; j < N2; j++) {
      A[i][j] = 0.;
      B[i][j] = 0.;
    }
  // !! Change limit of the stack : ulimit -s 32768
  // !! To see limit : ulimit -s

  g_assert (sp1 != NULL && sp2 != NULL);

  gint NU1 = sp1->periodic ? sp1->NUT : sp1->NU;
  gint NU2 = sp2->periodic ? sp2->NUT : sp2->NU;

  // Loop over the panels of the patch
  if (sp1 == sp2) {

    Spline2D * sp11 = sp1;
    
    while (sp11) {
      for ( ii = 0; ii < sp11->M*sp11->N; ii++) {
	SPPanel * spp = g_ptr_array_index (sp11->panels, ii);
	g_assert (spp != NULL);

	// Gauss outer-integration
	GaussPoints * gp = spp->outer;
	ustart = sp1->periodic ? ( gp->istart + sp11->fs_index - k1 + 1) : gp->istart;
	// Eventual shift for index for periodic free-surface
	vstart = gp->jstart;

	for ( m = 0; m < ng; m++) {
	  gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
	  for ( n = 0; n < ng; n++) {
	    gdouble wmn = M_PI2*g_array_index (gp->wJij, gdouble, m+n*ng);
	    gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);
	
	    // Loop over the splines whose support is included in the panel
	    for ( i = 0; i < k1; i++) {
	      gdouble wmni = wmn*gsl_matrix_get (Bu, i, 0);
	      gint itmp = (ustart + i)%NU1;
	      for ( j = 0; j < k1; j++) {
		gdouble wmnij = wmni*gsl_matrix_get (Bv, j, 0);
		gint indexi =  itmp + (vstart+j)*NU1;
		
		// First term of equation (2*pi*phi)
		for ( a = 0; a < k1; a++) {
		  gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		  gint atmp = (ustart+a)%NU1;
		  for ( b = 0; b < k1; b++) {
		    A[indexi][atmp + (vstart+b)*NU1] +=  wmnija*gsl_matrix_get (Bv, b, 0);
		  }
		}
	      }
	    }
	  }
	}
      }
      
      sp11 = sp11->next;
    }
  }


  gdouble * psi = g_malloc (k2*k2*sizeof(gdouble));
  gdouble * phi = g_malloc (k2*k2*sizeof(gdouble));

  Spline2D * sp22 = sp2;
  while (sp22) {
    Spline2D * sp11 = sp1;
    while (sp11) {

      for ( mu = 0; mu < sp22->M; mu++) {
	gint mui = mu + sp22->fs_index;
  	for ( nu = 0; nu < sp22->N; nu++) {
	  gboolean edge = (mu == 0 || nu == 0 || mu == sp22->M-1 || nu == sp22->N-1) ? TRUE : FALSE;
	  gint nui = nu*NU2;
  	  SPPanel * panel = g_ptr_array_index (sp22->panels, mu + nu*sp22->M);

  	  GPCell * panel_tree = gpcell_tree_new (panel);	  
  	  for ( ii = 0; ii < sp11->M*sp11->N; ii++) {
  	    SPPanel * spp = g_ptr_array_index (sp11->panels, ii);
   	
  	    // Gauss outer-integration
  	    GaussPoints * gp = spp->outer;
  	    vstart = gp->jstart;
	    ustart = sp1->periodic ? ( gp->istart + sp11->fs_index - k1 + 1) : gp->istart;

  	    // Temporary storage
  	    gdouble tmpBv[k1][ng];
  	    for ( i = 0; i < ng; i++ ) {
  	      gsl_matrix * Bv = g_ptr_array_index (gp->Bv, i);
  	      for ( j = 0; j < k1; j++) {
  		tmpBv[j][i] = gsl_matrix_get (Bv, j, 0);
  	      }
  	    }

  	    for ( m = 0; m < ng; m++) {
  	      gdouble um = g_array_index (gp->ui, gdouble, m);
  	      gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
  	      for ( n = 0; n < ng; n++) {
  		gdouble vn = g_array_index (gp->vj, gdouble, n);
  		gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
  		Point p = g_array_index (gp->Pi, Point, m+n*ng);

	     
  		// Spline coefficients
  		gint a;
  		for ( a = 0; a < k2*k2; a++)
  		  psi[a] = phi[a] = 0.;

  		if ( spp == panel )
  		  self_influence_coeffs (panel, um, vn, p, psi, phi);
  		else
  		  spline_near_field_influence_coeff_recursive (panel, panel_tree, p, psi, phi, edge);

  		// Loop over the splines whose support is included in the panel
  		for ( i = 0; i < k1; i++) {
  		  gdouble wmni = wmn*gsl_matrix_get (Bu, i, 0);
  		  gint indexi = (ustart + i) % NU1;
  		  for ( j = 0; j < k1; j++) {
  		    gdouble wmnij = wmni*tmpBv[j][n];
  		    gint indexij = indexi + (vstart+j)*NU1;

		    if ( (mui + k2) > NU2 ) {

		      for ( a = 0; a < k2; a++) {
			gint indexa = (mui + a) % NU2 + nui;
			gint aa = a;
			for ( b = 0; b < k2; b++) {
			  A[indexij][indexa] -= wmnij*phi[aa];
			  B[indexij][indexa] += wmnij*psi[aa];
			  indexa += NU2;
			  aa += k2;
			}
		      }

		    }
		    else {
		      
		      for ( a = 0; a < k2; a++) {
			gint indexa = mui + a + nui;
			gint aa = a;
			for ( b = 0; b < k2; b++) {
			  A[indexij][indexa] -= wmnij*phi[aa];
			  B[indexij][indexa] += wmnij*psi[aa];
			  indexa += NU2;
			  aa += k2;
			}
		      }

		    }

  		  }
  		}

  	      }
  	    }
	    	   
  	  }
      
  	  gpcell_tree_destroy (panel_tree);
  	}
      
      }
      sp11 = sp11->next;
    }
    sp22 = sp22->next;
  }
  g_free (psi);
  g_free (phi);

  // Copy the problem back to the old linear problem structure
  BoundarySubProblem * bsp = boundary_subproblem_new (N1, N2);

  bsp->neumann->istart = bsp->dirichlet->istart = sp2->istart;
  bsp->neumann->jstart = bsp->dirichlet->jstart = sp1->istart;

  FILE * fp = fopen ("check2.tmp","w");
  
  for ( j = 0; j < N1; j++) {
    for ( i = 0; i < N2; i++) {
      fprintf (fp, "%i %i %e %e \n", i, j, A[j][i], B[j][i]);
      gsl_matrix_set (bsp->neumann->A, j, i, A[j][i]);
      gsl_matrix_set (bsp->dirichlet->A, j, i, B[j][i]);
    }
    fprintf (fp, "\n");
  }
  fclose (fp);

  //fprintf (stderr, "start %i %i \n", bsp->neumann->istart, bsp->neumann->jstart);

  //g_assert_not_reached ();

  return bsp;
}

BoundarySubProblem *  build_boundary_subproblem_collocation (Spline2D * sp1, Spline2D * sp2,
							     SelfInfluenceFunc self_influence_coeffs)
{
  gint i, j, m, n, ii, a, b, mu, nu;
  /* size_t ustart, vstart; */
  gdouble M_PI2 = 2.*M_PI;
  gint N1 = sp1->NU*sp1->NV;
  gint N2 = sp2->NU*sp2->NV;
  gint k1 = sp1->k, k2 = sp2->k;
  gint ng = sp1->nouter;
  
  // Initializes the coefficients to 0
  gdouble A[N1][N2];
  gdouble B[N1][N2];

  for ( i = 0; i < N1; i++)
    for ( j = 0; j < N2; j++) {
      A[i][j] = 0.;
      B[i][j] = 0.;
    }
  // !! Change limit of the stack : ulimit -s 32768
  // !! To see limit : ulimit -s

  gint NU1 = sp1->NU;
  gint NU2 = sp2->NU;
  gint NV1 = sp1->NV;
  gint NV2 = sp2->NV;

  GrevillePoints * gr = sp1->gr;
  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;

  // Loop over the panels of the patch
  if (sp1 == sp2) {
    for ( i = 0; i < NU1; i++) {
      Bu = g_ptr_array_index (gr->Bu, i);
      ustart = sp1->periodic ? ( g_array_index (gr->ustart, gint, i) - k1 + 1) : g_array_index (gr->ustart, gint, i);

      for ( j = 0; j < NV1; j++) {
	Bv = g_ptr_array_index (gr->Bv, j);
	vstart = g_array_index (gr->vstart, gint, j);
	gint indexi = i + j*NU1;

	for ( a = 0; a < k1; a++) {
	  gint atmp = (ustart + a)%NU1;
	  gdouble val = M_PI2*gsl_matrix_get (Bu, a, 0);
	  for ( b = 0; b < k1; b++) {
	    A[indexi][atmp+ (b+vstart)*NU1] = val*gsl_matrix_get (Bv, b, 0);
	  }
	}
      }
    }
  }

  gdouble * psi = g_malloc (k2*k2*sizeof(gdouble));
  gdouble * phi = g_malloc (k2*k2*sizeof(gdouble));

  for ( mu = 0; mu < sp2->M; mu++) {
    for ( nu = 0; nu < sp2->N; nu++) {
      SPPanel * panel = g_ptr_array_index (sp2->panels, mu + nu*sp2->M);

      GPCell * panel_tree = gpcell_tree_new (panel);
      for ( i = 0; i < NU1; i++) {
	gdouble u = g_array_index (gr->ui, gdouble, i);
	for ( j = 0; j < NV1; j++) {
	  gdouble v = g_array_index (gr->vj, gdouble, j);
	  gint indexi = i + j*NU1;

	  Point p = g_array_index (gr->Pi, Point, i + j*NU1);
	     
	  // Spline coefficients
	  gint a;
	  for ( a = 0; a < k2*k2; a++)
	    psi[a] = phi[a] = 0.;

	  if ( sp1 == sp2 && u > panel->u0 && panel->u1 > u && v > panel->v0 && panel->v1 > v )
	    self_influence_coeffs (panel, u, v, p, psi, phi);
	  else {
	    spline_near_field_influence_coeff_recursive (panel, panel_tree, p, psi, phi, FALSE);
	  }

	  // Loop over the splines whose support is included in the panel
	  for ( a = 0; a < k2; a++) {
	    gint atmp = (mu + a) % NU2;
	    for ( b = 0; b < k2; b++) {
	      gint indexa = atmp + (nu+b)*NU2;

	      A[indexi][indexa] -= phi[a + b*k2];
	      B[indexi][indexa] += psi[a + b*k2];	  
	    }
	  }
	
	}
      }

      gpcell_tree_destroy (panel_tree);
    }
  } 
         
  g_free (psi);
  g_free (phi);

  // Copy the problem back to the old linear problem structure
  BoundarySubProblem * bsp = boundary_subproblem_new (N1, N2);

  bsp->neumann->istart = bsp->dirichlet->istart = sp2->istart;
  bsp->neumann->jstart = bsp->dirichlet->jstart = sp1->istart;

  /* FILE * fp = fopen ("check2.tmp","w"); */
  
  for ( j = 0; j < N1; j++) {
    for ( i = 0; i < N2; i++) {
      /* fprintf (fp, "%i %i %e %e \n", i, j, A[j][i], B[j][i]); */
      gsl_matrix_set (bsp->neumann->A, j, i, A[j][i]);
      gsl_matrix_set (bsp->dirichlet->A, j, i, B[j][i]);
    }
    /* fprintf (fp, "\n"); */
  }
  /* fclose (fp); */

  return bsp;
}

BoundarySubProblem *  build_boundary_subproblem_collocation_legendre (Spline2D * sp1,
								      Spline2D * sp2,
								      SelfInfluenceFunc self_influence_coeffs)
{
  gint i, j, k, l, m, n, ii, a, b, mu, nu;
  /* size_t ustart, vstart; */
  gdouble M_PI2 = 2.*M_PI;
  gint k1 = sp1->k, k2 = sp2->k;
  gint ng = sp1->nouter;
  gint N1 = sp1->NU*sp1->NV*ng*ng;
  gint N2 = sp2->NU*sp2->NV;
  
  // Initializes the coefficients to 0
  gdouble A[N1][N2];
  gdouble B[N1][N2];

  for ( i = 0; i < N1; i++)
    for ( j = 0; j < N2; j++) {
      A[i][j] = 0.;
      B[i][j] = 0.;
    }
  // !! Change limit of the stack : ulimit -s 32768
  // !! To see limit : ulimit -s

  gint NU1 = sp1->NU;
  gint NU2 = sp2->NU;
  gint NV1 = sp1->NV;
  gint NV2 = sp2->NV;

  GrevillePoints * gr = sp1->gr;
  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;

  // Loop over the panels
  if (sp1 == sp2) {
    Spline2D * sp = sp1;
    gint index = 0;
    while (sp) {

        for ( i = 0; i < sp->M; i++) {
	  for ( j = 0; j < sp->N; j++) {
	    SPPanel * panel = g_ptr_array_index (sp->panels, i + j*sp->M);
	    GaussPoints * gp = panel->outer;

	    for ( m = 0; m < ng; m++ ) {
	      Bu = g_ptr_array_index (gp->Bu, m);
	      ustart = sp->periodic ? (gp->istart - k1 + 1) : gp->istart;

	      for ( n = 0; n < ng; n++ ) {
		gint indexi = ng*i+m + (ng*j+n)*NU1*ng;
		Bv = g_ptr_array_index (gp->Bv, n);
		vstart = gp->jstart;

		
		for ( a = 0; a < k1; a++) {
		  gint atmp = (ustart + a)%NU1;
		  gdouble val = M_PI2*gsl_matrix_get (Bu, a, 0);
		  for ( b = 0; b < k1; b++) {
		    A[indexi][atmp+ (b+vstart)*NU1] = val*gsl_matrix_get (Bv, b, 0);
		  }
		}

	      }
	    }

	  }
	}

	index += ng*sp->M;
	sp =sp->next;
    }
  }

  gdouble * psi = g_malloc (k2*k2*sizeof(gdouble));
  gdouble * phi = g_malloc (k2*k2*sizeof(gdouble));

  Spline2D * sp_2 = sp2;
  while (sp_2) {
    for ( mu = 0; mu < sp_2->M; mu++) {
      for ( nu = 0; nu < sp_2->N; nu++) {
	SPPanel * panel = g_ptr_array_index (sp_2->panels, mu + nu*sp_2->M);

	GPCell * panel_tree = gpcell_tree_new (panel);
	Spline2D * sp_1 = sp1;
	while (sp_1) {
	  for ( i = 0; i < sp1->M; i++) {
	    // gdouble u = g_array_index (gr->ui, gdouble, i);
	    for ( j = 0; j < sp1->N; j++) {
	      SPPanel * spp = g_ptr_array_index (sp_1->panels, i + j*sp_1->M);
	      GaussPoints * gp = spp->outer;

	      for ( m = 0; m < ng; m++ ) {
		gdouble u = g_array_index (gp->ui, gdouble, m);
		//Bu = g_ptr_array_index (gp->Bu, m);
		//ustart = sp->periodic ? (gp->istart - k1 + 1) : gp->istart;

		for ( n = 0; n < ng; n++ ) {
		  gint indexi = ng*i+m + (ng*j+n)*NU1*ng;
		  //Bv = g_ptr_array_index (gp->Bv, n);
		  //vstart = gp->jstart;
		  gdouble v = g_array_index (gp->vj, gdouble, n);

		  //gdouble v = g_array_index (gr->vj, gdouble, j);
		  // gint indexi = i + j*NU1;

		  Point p = g_array_index (gr->Pi, Point, i + j*NU1);
	     
		  // Spline coefficients
		  gint a;
		  for ( a = 0; a < k2*k2; a++)
		    psi[a] = phi[a] = 0.;

		  if ( sp1 == sp_2 && u > panel->u0 && panel->u1 > u && v > panel->v0 && panel->v1 > v )
		    self_influence_coeffs (panel, u, v, p, psi, phi);
		  else {
		    spline_near_field_influence_coeff_recursive (panel, panel_tree, p, psi, phi, FALSE);
		  }

		  // Loop over the splines whose support is included in the panel
		  for ( a = 0; a < k2; a++) {
		    gint atmp = (mu + a) % NU2;
		    for ( b = 0; b < k2; b++) {
		      gint indexa = atmp + (nu+b)*NU2;

		      A[indexi][indexa] -= phi[a + b*k2];
		      B[indexi][indexa] += psi[a + b*k2];	  
		    }
		  }
	
		}
	      }

	    }
	  }
	  sp_1 = sp_1->next;
	}

	gpcell_tree_destroy (panel_tree);
      }
    }
    sp_2 = sp_2->next;
  }
         
  g_free (psi);
  g_free (phi);

  // Copy the problem back to the old linear problem structure
  BoundarySubProblem * bsp = boundary_subproblem_new (N1, N2);

  bsp->neumann->istart = bsp->dirichlet->istart = sp2->istart;
  bsp->neumann->jstart = bsp->dirichlet->jstart = sp1->istart;

  /* FILE * fp = fopen ("check2.tmp","w"); */
  
  for ( j = 0; j < N1; j++) {
    for ( i = 0; i < N2; i++) {
      /* fprintf (fp, "%i %i %e %e \n", i, j, A[j][i], B[j][i]); */
      gsl_matrix_set (bsp->neumann->A, j, i, A[j][i]);
      gsl_matrix_set (bsp->dirichlet->A, j, i, B[j][i]);
    }
    /* fprintf (fp, "\n"); */
  }
  /* fclose (fp); */

  return bsp;
}

void boundary_problem_assemble_neumann_rhs (BoundaryProblem * bp,
					    GPtrArray * sub_problems)
{
  gint i, p;

  gsl_vector_set_zero (bp->rhs);

  for ( p = 0; p < sub_problems->len; p++) {
    BoundarySubProblem * bsp = g_ptr_array_index (sub_problems, p);
    BoundaryProblem * sub_problem = bsp->neumann;

    for ( i = 0; i < sub_problem->rhs->size; i++) {
      gdouble tmp = gsl_vector_get (bp->rhs, sub_problem->jstart + i) + gsl_vector_get (sub_problem->rhs, i);

      gsl_vector_set (bp->rhs, sub_problem->jstart + i, tmp);
    }
  }

  FILE * fp = fopen ("neumann-rhs.tmp","w");
  for ( i = 0; i < bp->rhs->size; i++)
    fprintf (fp, "%i %e\n", i, gsl_vector_get (bp->rhs, i));
  fclose (fp);
}

void boundary_problem_assemble_neumann (BoundaryProblem * bp,
					 GPtrArray * sub_problems)
{
  gint i, j, p;

  for ( p = 0; p < sub_problems->len; p++) {
    BoundarySubProblem * sub_problem = g_ptr_array_index (sub_problems, p);
    BoundaryProblem * neumann = sub_problem->neumann;

    for ( i = 0; i < neumann->A->size1; i++)
      for ( j = 0; j < neumann->A->size2; j++)
	gsl_matrix_set (bp->A, neumann->istart + j, neumann->jstart + i,
			gsl_matrix_get (neumann->A, i, j));
  }

#if DEBUG
    FILE * fp = fopen ("neumann-bp.tmp","w");
    for ( i = 0; i < bp->A->size1; i++) {
      for ( j = 0; j < bp->A->size2; j++)
  	fprintf(fp, "%i %i %e  \n", i, j, gsl_matrix_get (bp->A, i, j));
      fprintf(fp,"\n");
    }
    fclose (fp);
#endif
}

void boundary_problem_assemble_dirichlet_rhs (BoundaryProblem * bp,
					      GPtrArray * sub_problems)
{
  gint i, p;

  gsl_vector_set_zero (bp->rhs);

  for ( p = 0; p < sub_problems->len; p++) {
    BoundarySubProblem * bsp = g_ptr_array_index (sub_problems, p);
    BoundaryProblem * sub_problem = bsp->dirichlet;
    
    for ( i = 0; i < sub_problem->rhs->size; i++) {
      gdouble tmp = gsl_vector_get (bp->rhs, sub_problem->jstart + i) + gsl_vector_get (sub_problem->rhs, i);

      gsl_vector_set (bp->rhs, sub_problem->jstart + i, tmp);
    }
  }
  
  /* FILE * fp = fopen ("dirichlet-rhs.tmp","w"); */
  /* for ( i = 0; i < bp->rhs->size; i++) */
  /*   fprintf (fp, "%i %f\n", i, gsl_vector_get (bp->rhs, i)); */
  /* fclose (fp); */
}

void boundary_problem_assemble_dirichlet (BoundaryProblem * bp,
					   GPtrArray * sub_problems)
{
  gint i, j, p;

  for ( p = 0; p < sub_problems->len; p++) {
    BoundarySubProblem * sub_problem = g_ptr_array_index (sub_problems, p);
    BoundaryProblem * dirichlet = sub_problem->dirichlet;

    for ( i = 0; i < dirichlet->A->size1; i++)
      for ( j = 0; j < dirichlet->A->size2; j++)
	gsl_matrix_set (bp->A, dirichlet->istart + j, dirichlet->jstart + i,
			gsl_matrix_get (dirichlet->A, i, j));
  }

#if DEBUG
  FILE * fp = fopen ("dirichlet-bp.tmp","w");
  for ( i = 0; i < bp->A->size1; i++) {
      for ( j = 0; j < bp->A->size2; j++)
  	fprintf(fp, "%i %i %e  \n", i, j, gsl_matrix_get (bp->A, i, j));
      fprintf(fp,"\n");
  }
  fclose (fp);
#endif
}

CCSProblem * spline2d_build_galerkin_fit_matrix (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    ustart = gp->istart;
    vstart = gp->jstart;
    gint ng = sp->nouter; // Order of outer Gauss-Legendre rule

    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get(Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get(Bv, j-vstart, 0);
	    gint indexi = i + j*NU;

	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get(Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*NU] +=  wmnija*gsl_matrix_get(Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * fit = ccs_problem_new ();
  gint count = 0;
  FILE * fp = fopen ("galerkin.tmp","w");
  for ( i = 0; i < size; i++) {
    g_array_append_val (fit->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (fit->matrix, A0[i][j]);
  	g_array_append_val (fit->column, j);
  	count++;
      }
      fprintf (fp, "%i %i %f\n", i,j,A0[i][j]);
    }
    fprintf (fp, "\n");
  }
  g_array_append_val (fit->index, count);
  fclose (fp);
  return fit;
}

CCSProblem * spline2d_build_galerkin_fit_matrix_bc (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  gsl_matrix * Bu, * Bv;
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ustart = gp->istart;
    gint vstart = gp->jstart;
    gint ng = sp->nouter; // Order of outer Gauss-Legendre rule

    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get(Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get(Bv, j-vstart, 0);
	    gint indexi = i + j*NU;

	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get(Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*NU] +=  wmnija*gsl_matrix_get(Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

#if 1
  // Dirichlet on edges
  gsl_matrix * Bu2 = gsl_matrix_alloc (sp->k, 3);
  gsl_matrix * Bv2 = gsl_matrix_alloc (sp->k, 3);
  GrevillePoints * gr = sp->gr;
  size_t ustart, uend, vstart, vend;
  for ( i = 0; i < sp->NU; i++) {
    gint index1 = i;
    gdouble u = g_array_index (gr->ui, gdouble, i);
    gdouble v = 0.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_deriv_eval_nonzero (u, 2, Bu2, &ustart, &uend, sp->w_u, sp->wd_u);
    gsl_bspline_deriv_eval_nonzero (v, 2, Bv2, &vstart, &vend, sp->w_v, sp->wd_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 0);
      for ( n = 0; n < sp->k; n++) {
	gdouble cv = gsl_matrix_get (Bv2, n, 0);
	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }

    index1 = i + sp->NU;

    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 0);
      for ( n = 0; n < sp->k; n++) {
	gdouble cv = gsl_matrix_get (Bv2, n, 1);
	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }
    
    /* index1 = i + 2.*sp->NU; */

    /* for ( j = 0; j < size; j++) */
    /*   A0[j][index1] = 0.; */

    /* for ( m = 0; m < sp->k; m++) { */
    /*   gdouble cu = gsl_matrix_get (Bu2, m, 0); */
    /*   for ( n = 0; n < sp->k; n++) { */
    /* 	gdouble cv = gsl_matrix_get (Bv2, n, 2); */
    /* 	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv; */
    /*   } */
    /* } */

    /*****************************/
    index1 = i + NU*(NV-1);
    v = 1.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_deriv_eval_nonzero (u, 2, Bu2, &ustart, &uend, sp->w_u, sp->wd_u);
    gsl_bspline_deriv_eval_nonzero (v, 2, Bv2, &vstart, &vend, sp->w_v, sp->wd_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 0);
      for ( n = 0; n < sp->k; n++) {
    	gdouble cv = gsl_matrix_get (Bv2, n, 0);
    	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }

    index1 = i + NU*(NV-2);

    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 0);
      for ( n = 0; n < sp->k; n++) {
	gdouble cv = gsl_matrix_get (Bv2, n, 1);
	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }

    /* index1 = i + NU*(NV-3); */

    /* for ( j = 0; j < size; j++) */
    /*   A0[j][index1] = 0.; */

    /* for ( m = 0; m < sp->k; m++) { */
    /*   gdouble cu = gsl_matrix_get (Bu2, m, 0); */
    /*   for ( n = 0; n < sp->k; n++) { */
    /* 	gdouble cv = gsl_matrix_get (Bv2, n, 2); */
    /* 	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv; */
    /*   } */
    /* } */
  }

  for ( i = 0; i < sp->NV; i++) {
    gint index1 = i*sp->NU;
    gdouble v = g_array_index (gr->vj, gdouble, i);
    gdouble u = 0.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_deriv_eval_nonzero (u, 2, Bu2, &ustart, &uend, sp->w_u, sp->wd_u);
    gsl_bspline_deriv_eval_nonzero (v, 2, Bv2, &vstart, &vend, sp->w_v, sp->wd_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 0);
      for ( n = 0; n < sp->k; n++) {
  	gdouble cv = gsl_matrix_get (Bv2, n, 0);
  	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }

    index1 = i*sp->NU + 1;

    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 1);
      for ( n = 0; n < sp->k; n++) {
  	gdouble cv = gsl_matrix_get (Bv2, n, 0);
  	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }
    
    /* index1 = i*sp->NU + 2; */

    /* for ( j = 0; j < size; j++) */
    /*   A0[j][index1] = 0.; */

    /* for ( m = 0; m < sp->k; m++) { */
    /*   gdouble cu = gsl_matrix_get (Bu2, m, 2); */
    /*   for ( n = 0; n < sp->k; n++) { */
    /* 	gdouble cv = gsl_matrix_get (Bv2, n, 0); */
    /* 	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv; */
    /*   } */
    /* } */

    /**************************/
    index1 = NU-1 + i*sp->NU;
    u = 1.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_deriv_eval_nonzero (u, 2, Bu2, &ustart, &uend, sp->w_u, sp->wd_u);
    gsl_bspline_deriv_eval_nonzero (v, 2, Bv2, &vstart, &vend, sp->w_v, sp->wd_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 0);
      for ( n = 0; n < sp->k; n++) {
    	gdouble cv = gsl_matrix_get (Bv2, n, 0);
    	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }

    index1 = NU-2 + i*sp->NU;

    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_matrix_get (Bu2, m, 1);
      for ( n = 0; n < sp->k; n++) {
  	gdouble cv = gsl_matrix_get (Bv2, n, 0);
  	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }

    /* index1 = NU-3 + i*sp->NU; */

    /* for ( j = 0; j < size; j++) */
    /*   A0[j][index1] = 0.; */

    /* for ( m = 0; m < sp->k; m++) { */
    /*   gdouble cu = gsl_matrix_get (Bu2, m, 2); */
    /*   for ( n = 0; n < sp->k; n++) { */
    /* 	gdouble cv = gsl_matrix_get (Bv2, n, 0); */
    /* 	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv; */
    /*   } */
    /* } */
  }
#endif

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * fit = ccs_problem_new ();
  gint count = 0;
  FILE * fp = fopen ("galerkin.tmp","w");
  for ( i = 0; i < size; i++) {
    g_array_append_val (fit->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (fit->matrix, A0[i][j]);
  	g_array_append_val (fit->column, j);
  	count++;
      }
      fprintf (fp, "%i %i %f\n", i,j,A0[i][j]);
    }
    fprintf (fp, "\n");
  }
  g_array_append_val (fit->index, count);
  fclose (fp);
  return fit;
}

CCSProblem * spline2d_build_greville_fit_matrix (Spline2D * sp)
{
  gint i, j, a, b;
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble A0[size][size];

  for ( i = 0; i < size; i++)
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;
  
  for ( i = 0; i < NU; i++) {
    gdouble u = g_array_index (gr->ui, gdouble, i);
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i);
    gint istart = g_array_index (gr->ustart, gint, i);
    for ( j = 0; j < NV; j++) {
      gdouble v = g_array_index (gr->vj, gdouble, j);
      gsl_matrix * Bv = g_ptr_array_index (gr->Bv, j);
      gint jstart = g_array_index (gr->vstart, gint, j);

      for ( a = istart; a < istart + sp->k; a++) {
      	for ( b = jstart; b < jstart + sp->k; b++) {
      	  A0[i+j*NU][a+b*NU] = gsl_matrix_get (Bu, a-istart, 0)*gsl_matrix_get (Bv, b-jstart, 0);
      	}
      }
    }
  }

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
  }
  g_array_append_val (fit->index, count);

  return fit;
}

void spline2d_copy_problem_solution (Spline2D * sp, gsl_vector * lhs, gint var)
{
  gint i, j;

  if (var > 2) {
    for ( i = 0; i < sp->NU; i++) {
      for ( j = 0; j < sp->NV; j++) {
	SplineCoeffs * sc =
	  g_ptr_array_index (sp->coeffs, i + j*sp->NU);
	sc->v[var] = gsl_vector_get (lhs, i + j*sp->NU);
      }
    }
  }
  else { // x,y,z
    for ( i = 0; i < sp->NXU; i++) {
      for ( j = 0; j < sp->NXV; j++) {
	SplineCoeffs * sc =
	  g_ptr_array_index (sp->coeffs, i + j*sp->NXU);
	sc->v[var] = gsl_vector_get (lhs, i + j*sp->NXU);
      }
    }
  }
}

void spline2d_add_problem_solution (Spline2D * sp, gsl_vector * lhs, gint var, gint var0)
{
  gint i, j;

  if (var > 2) {
    for ( i = 0; i < sp->NU; i++) {
      for ( j = 0; j < sp->NV; j++) {
	SplineCoeffs * sc =
	  g_ptr_array_index (sp->coeffs, i + j*sp->NU);
	sc->v[var] = sc->v[var0] + gsl_vector_get (lhs, i + j*sp->NU);
      }
    }
  }
  else { // x,y,z
    for ( i = 0; i < sp->NXU; i++) {
      for ( j = 0; j < sp->NXV; j++) {
	SplineCoeffs * sc =
	  g_ptr_array_index (sp->coeffs, i + j*sp->NXU);
	sc->v[var] = sc->v[var0] + gsl_vector_get (lhs, i + j*sp->NXU);
      }
    }
  }
}

void spline2d_build_freesurface_galerkin_fit_matrix (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ustart = gp->istart;
    gint vstart = gp->jstart;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule

    for ( m = 0; m < ng; m++) {
      gdouble um = g_array_index (gp->ui, gdouble, m);
      gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble vn = g_array_index (gp->vj, gdouble, n);
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get(Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get(Bv, j-vstart, 0);
	    gint indexi = i + j*NU;

	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get(Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*NU] +=  wmnija*gsl_matrix_get(Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

  // Continuity of variable on boundary
  for ( i = 0; i < NV; i++) {
    gint index1 = (i+1)*NU - 1;
    gint index2 = i*NU;
    
    for ( j = 0; j < size; j++)
      A0[j][index1] += A0[j][index2];

    for ( j = 0; j < size; j++)
      A0[j][index2] = 0.;

    A0[index1][index2] = 1.;
    A0[index2][index2] = -1.;
  }

  /* FILE * fp = fopen ("galerkin.tmp","w"); */
  /* for ( i = 0; i < size; i++) { */
  /*   for ( j = 0; j < size; j++) { */
  /*     fprintf (fp, "%i %i %f \n", i, j, A0[i][j]); */
  /*   } */
  /*   fprintf (fp, "\n"); */
  /* } */
  /* fclose (fp); */

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  sp->fit = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (sp->fit->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (sp->fit->matrix, A0[i][j]);
  	g_array_append_val (sp->fit->column, j);
  	count++;
      }
    }
  }
  g_array_append_val (sp->fit->index, count);
}

void spline2d_build_freesurface_noflux_galerkin_fit_matrix (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ustart = gp->istart;
    gint vstart = gp->jstart;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule

    for ( m = 0; m < ng; m++) {
      gdouble um = g_array_index (gp->ui, gdouble, m);
      gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble vn = g_array_index (gp->vj, gdouble, n);
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get(Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get(Bv, j-vstart, 0);
	    gint indexi = i + j*NU;

	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get(Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*NU] +=  wmnija*gsl_matrix_get(Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

/* #ifdef NO_FLUX */
  // On the solid boundary dn phi = - dn phi0
  // For now, dn phi = 0 i.e. dv phi = 0
  size_t ustart, uend, vstart, vend;
  gsl_vector * Bu = gsl_vector_alloc (sp->k);
  gsl_matrix * Bv = gsl_matrix_alloc (sp->k, 2);
  for ( i = 0; i < NU; i++) {
    gint index1 = i;
    gdouble u = gsl_bspline_greville_abscissa (i, sp->w_u);
    gdouble v = 0.;

    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_eval_nonzero (MIN(1.-1e-12,u), Bu, &ustart, &uend, sp->w_u);
    gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_vector_get (Bu, m);
      for ( n = 0; n < sp->k; n++) {
    	gdouble cv = gsl_matrix_get (Bv, n, 1);
	A0[ustart+m + (vstart+n)*NU][index1] = cu*cv;
      }
    }
  }
  gsl_vector_free (Bu);
  gsl_matrix_free (Bv);
/* #endif 0 */

  // Continuity of variable on boundary
  for ( i = 0; i < NV; i++) {
    gint index1 = (i+1)*NU - 1;
    gint index2 = i*NU;
    
    for ( j = 0; j < size; j++)
      A0[j][index1] += A0[j][index2];

    for ( j = 0; j < size; j++)
      A0[j][index2] = 0.;

    A0[index1][index2] = 1.;
    A0[index2][index2] = -1.;
  }

  /* FILE * fp = fopen ("galerkin.tmp","w"); */
  /* for ( i = 0; i < size; i++) { */
  /*   for ( j = 0; j < size; j++) { */
  /*     fprintf (fp, "%i %i %f \n", i, j, A0[i][j]); */
  /*   } */
  /*   fprintf (fp, "\n"); */
  /* } */
  /* fclose (fp); */

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  sp->fit = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (sp->fit->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (sp->fit->matrix, A0[i][j]);
  	g_array_append_val (sp->fit->column, j);
  	count++;
      }
    }
  }
  g_array_append_val (sp->fit->index, count);
}

void spline2d_build_boundary_galerkin_fit_matrix (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  gint size = sp->NU*sp->NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
    for ( m = 0; m < ng; m++) {
      gdouble um = g_array_index (gp->ui, gdouble, m);
      gsl_matrix * Bu = g_ptr_array_index (gp->Bu, m);
      gint ustart = gp->istart;/* g_array_index (gp->ustart, gint, m); */
	  
      for ( n = 0; n < ng; n++) {
	gdouble vn = g_array_index (gp->vj, gdouble, n);
	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	/* gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng); */
	gsl_matrix * Bv = g_ptr_array_index (gp->Bv, n);
	gint vstart = gp->jstart;/* g_array_index (gp->vstart, gint, n); */

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get(Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get(Bv, j-vstart, 0);
	    gint indexi = i + j*sp->NU;

	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get(Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*sp->NU] +=  wmnija*gsl_matrix_get (Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

  // We replace one "line" of the matrix by the boundary conditions i.e. no flux through the boundary
  /* for ( i = 0; i < sp->N; i++ ) { */
  /*   gdouble u = i == 0 ? 0 : i == sp->N-1 ? 1 : gsl_bspline_greville_abscissa (i, sp->w_v); */
  /* } */

  // We need to "close" the free-surface and impose that the boundary nodes where both
  // ends of the circle meet reconnect
  for ( i = 0; i < sp->N; i++ ) {
    j = sp->M-1;
    gint index1 = i*sp->N;
    gint index2 = (i+1)*sp->N-1;

    A0[index1][index1] = 1.;
    A0[index1][index2] = -1.;
    // rhs
  }


  FILE * fp = fopen ("galerkin.tmp","w");
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      fprintf (fp, "%i %i %f \n", i, j, A0[i][j]);
    }
    fprintf (fp, "\n");
  }
  fclose (fp);

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  sp->fit = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (sp->fit->index, count);
    for ( j = 0; j < size; j++) {
      if (A0[i][j] != 0.) {
  	g_array_append_val (sp->fit->matrix, A0[i][j]);
  	g_array_append_val (sp->fit->column, j);
  	count++;
      }
    }
  }
  g_array_append_val (sp->fit->index, count);
}

void spline2d_fit_galerkin (Spline2D * sp, GaussFunc func, gpointer data, gint var)
{

  if (sp->fit == NULL)
    sp->fit = sp->build_fit_matrix (sp);

  sp->rhs = sp->build_fit_rhs (sp, func, data, NULL, NULL, sp->rhs);

  ccs_problem_lu_solve (sp->fit, sp->rhs);

  sp->copy_fit_solution (sp, sp->rhs, var);

  // gsl_vector_free (rhs);
}

static void set_values_at_boundaries (Spline2D * sp, Spline2DFunc func, gpointer data, gint var)
{
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = 2*NU+2*NV-4;
  gdouble A[size][size];
  gint jj[size];
  gint ii[size];
  gsl_vector * rhs = gsl_vector_alloc (size);
  gint i, j;

  for ( i = 0; i < size; i++ ) {
    ii[i] = 0;
    jj[i] = 0;
    for ( j = 0; j < size; j++ )
      A[i][j] = 0.;
  }

  for ( i = 0; i < NU; i++) {
    ii[i] = i;
    ii[i+NU] = i;
    jj[i] = 0;
    jj[i+NU] = NV-1;
  }

  for ( j = 1; j < NV-1; j++) {
    ii[2*NU+j-1] = 0;
    jj[2*NU+j-1] = j;
    ii[2*NU+NV-2+j-1] = NU-1;
    jj[2*NU+NV-2+j-1] = j;
  }

  for ( i = 0; i < size; i++) {
    gdouble u = ii[i] == 0 ? 0 : ii[i] == NU-1 ? 1 : /* g_array_index (gr->ui, gdouble, ii[i]) */gsl_bspline_greville_abscissa (ii[i], sp->w_u);
    gdouble v = jj[i] == 0 ? 0 : jj[i] == NV-1 ? 1 : /* g_array_index (gr->ui, gdouble, jj[i]) */gsl_bspline_greville_abscissa (jj[i], sp->w_v);
   
    for ( j = 0; j < size; j++) {
      A[j][i] = spline2d_eval_spline (sp, ii[j], jj[j], u, v);
    }
    gsl_vector_set (rhs, i, func (sp, u, v, data));
  }

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * css = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (css->index, count);
    for ( j = 0; j < size; j++) {
      if (A[i][j] != 0.) {
  	g_array_append_val (css->matrix, A[i][j]);
  	g_array_append_val (css->column, j);
  	count++;
      }
    }
  }
  g_array_append_val (css->index, count);

  ccs_problem_lu_solve (css, rhs);

  for ( i = 0; i < size; i++)
    coeff_assign (sp, ii[i], jj[i], var, gsl_vector_get (rhs, i));

  ccs_problem_destroy (css);
  gsl_vector_free (rhs);
}

void spline2d_fit_greville_border (Spline2D * sp, Spline2DFunc func, gpointer data, gint var)
{
  gint i, j, a, b;
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = (NU-2)*(NV-2);
  gsl_vector * rhs = gsl_vector_alloc (size);
  gdouble A0[size][size];
  
  set_values_at_boundaries (sp, func, data, var);

  for ( i = 0; i < size; i++)
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;

  for ( i = 0; i < NU-2; i++) {
    gdouble u = g_array_index (gr->ui, gdouble, i+1);
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i+1);
    gint istart = g_array_index (gr->ustart, gint, i+1);
    for ( j = 0; j < NV-2; j++) {
      gdouble v = g_array_index (gr->vj, gdouble, j+1);
      gsl_matrix * Bv = g_ptr_array_index (gr->Bv, j+1);
      gint jstart = g_array_index (gr->vstart, gint, j+1);

      for ( a = MAX(1, istart); a < MIN(istart + sp->k, NU-1); a++) {
      	for ( b = MAX(1, jstart); b < MIN(jstart + sp->k, NV-1); b++) {
	    A0[i+j*(NU-2)][(a-1)+(b-1)*(NU-2)] = gsl_matrix_get (Bu, a-istart, 0)*gsl_matrix_get (Bv, b-jstart, 0);
      	}
      }

      gsl_vector_set (rhs, i+j*(NU-2), func (sp, u, v, data) - spline2d_eval_greville_point (sp, gr, i+1, j+1, var));
    }
  }

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

  ccs_problem_lu_solve (css, rhs);

  for ( i = 0; i < NU-2; i++) {
    for ( j = 0; j < NV-2; j++) {
      coeff_assign (sp, i+1, j+1, var, gsl_vector_get (rhs, i+j*(NU-2)));
    }
  }

  gsl_vector_free (rhs);
  ccs_problem_destroy (css);
}

void spline2d_fit_greville (Spline2D * sp, GrevilleFunc func, gpointer data, gint var)
{
  gint i, j, a, b;
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gsl_vector * rhs = gsl_vector_alloc (size);
  gdouble A0[size][size];

  for ( i = 0; i < size; i++)
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;
  
  for ( i = 0; i < NU; i++) {
    gdouble u = g_array_index (gr->ui, gdouble, i);
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i);
    gint istart = g_array_index (gr->ustart, gint, i);
    for ( j = 0; j < NV; j++) {
      gdouble v = g_array_index (gr->vj, gdouble, j);
      gsl_matrix * Bv = g_ptr_array_index (gr->Bv, j);
      gint jstart = g_array_index (gr->vstart, gint, j);

      for ( a = istart; a < istart + sp->k; a++) {
      	for ( b = jstart; b < jstart + sp->k; b++) {
      	  A0[i+j*NU][a+b*NU] = gsl_matrix_get (Bu, a-istart, 0)*gsl_matrix_get (Bv, b-jstart, 0);
      	}
      }
      gsl_vector_set (rhs, i+j*NU, func (sp, i, j, data));
    }
  }

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

  ccs_problem_lu_solve (css, rhs);

  for ( i = 0; i < sp->NU; i++) {
    for ( j = 0; j < sp->NV; j++) {
      coeff_assign (sp, i, j, var, gsl_vector_get (rhs, i+j*sp->NU));
    }
  }

  gsl_vector_free (rhs);
}

void spline2d_fit_greville_periodic (Spline2D * sp, GrevilleFunc func, gpointer data, gint var)
{
  gint i, j, a, b;
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gsl_vector * rhs = gsl_vector_alloc (size);
  gdouble A0[size][size];

  for ( i = 0; i < size; i++)
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;
  
  gint k = sp->k;
  for ( i = 0; i < NU; i++) {
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i);
    gint ustart = g_array_index (gr->ustart, gint, i);
    for ( j = 0; j < NV; j++) {
      gint vstart = g_array_index (gr->vstart, gint, j);
      gsl_matrix * Bv = g_ptr_array_index (gr->Bv, j);
      gint indexi = i+j*NU;
      
      for ( a = 0; a < k; a++) {
	gint atmp = (ustart+a-k+1)%NU;
	for ( b = 0; b < k; b++) {
	  A0[atmp + (vstart+b)*NU][indexi] = gsl_matrix_get (Bu, a, 0)*gsl_matrix_get (Bv, b, 0);
	}
      }

      gsl_vector_set (rhs, indexi, func (sp, i, j, data));
    }
  } 

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

  ccs_problem_lu_solve (css, rhs);

  for ( i = 0; i < NU+k-1; i++ ) {
    gint indexi = (i + k - 1)%(NU+k-1);
    for ( j = 0; j < NV; j++) {     
      coeff_assign (sp, indexi, j, var, gsl_vector_get (rhs, i%NU+j*NU));
    }
  }

  gsl_vector_free (rhs);
  ccs_problem_destroy (css);
}

static void set_values_at_border_periodic (Spline2D * sp,
					   GrevilleFunc func,
					   gpointer data, gint var)
{
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = 2*NU;
  gdouble A[size][size];
  gint k = sp->k;

  gsl_vector * rhs = gsl_vector_alloc (size);
  gint i, j, a, b;

  for ( i = 0; i < size; i++ ) {
    for ( j = 0; j < size; j++ )
      A[i][j] = 0.;
  }

  for ( i = 0; i < NU; i++) {
    gint ustart = g_array_index (gr->ustart, gint, i);
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i);
    gsl_matrix * Bv0 = g_ptr_array_index (gr->Bv, 0);
    gsl_matrix * Bv1 = g_ptr_array_index (gr->Bv, NV-1);


    for ( a = 0; a < k; a++) {
      gint indexa = (ustart+a-k+1)%NU;
      A[indexa][i] +=  gsl_matrix_get (Bu, a, 0) * gsl_matrix_get (Bv0, 0 ,0);
      A[indexa+NU][i+NU] += gsl_matrix_get (Bu, a, 0) * gsl_matrix_get (Bv1, k-1 ,0) ;
    }

    gsl_vector_set (rhs, i, func (sp, i, 0, data));
    gsl_vector_set (rhs, i+NU, func (sp, i, NV-1, data));
  }

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * css = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (css->index, count);
    for ( j = 0; j < size; j++) {
      if (A[i][j] != 0.) {
  	g_array_append_val (css->matrix, A[i][j]);
  	g_array_append_val (css->column, j);
  	count++;
      }
    }
  }
  g_array_append_val (css->index, count);

  ccs_problem_lu_solve (css, rhs);

  for ( i = 0; i < NU+k-1; i++ ) {
    gint indexi = (i + k - 1)%(NU+k-1);
    coeff_assign (sp, indexi, 0, var, gsl_vector_get (rhs, i%NU));
    coeff_assign (sp, indexi, NV-1, var, gsl_vector_get (rhs, i%NU+NU));
  }

  ccs_problem_destroy (css);
  gsl_vector_free (rhs);
}

void spline2d_fit_greville_periodic_border (Spline2D * sp,
					    GrevilleFunc func,
					    gpointer data, gint var)
{
  gint i, j, a, b;
  GrevillePoints * gr = sp->gr;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*(NV-2);
  gsl_vector * rhs = gsl_vector_alloc (size);
  gdouble A0[size][size];

  coeff_set_var_to_zero (sp, var);

  set_values_at_border_periodic (sp, func, data, var);

  for ( i = 0; i < size; i++)
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;
  
  gint k = sp->k;
  for ( i = 0; i < NU; i++) {
    gsl_matrix * Bu = g_ptr_array_index (gr->Bu, i);
    gint ustart = g_array_index (gr->ustart, gint, i);
    for ( j = 0; j < NV-2; j++) {
      gint vstart = g_array_index (gr->vstart, gint, j+1);
      gsl_matrix * Bv = g_ptr_array_index (gr->Bv, j+1);
      gint indexi = i+j*NU;
      
      for ( a = ustart; a < ustart + k; a++) {
	gint atmp = (a-k+1)%NU;
	for ( b = MAX(1, vstart); b < MIN(vstart + k, NV-1); b++) {
	  A0[atmp + (b-1)*NU][indexi] = gsl_matrix_get (Bu, a-ustart, 0)*gsl_matrix_get (Bv, b-vstart, 0);
	}
      }

      gdouble u = g_array_index (gr->ui, gdouble, i);
      gdouble v = g_array_index (gr->vj, gdouble, j+1);
      gsl_vector_set (rhs, indexi, func (sp, i, j, data) - spline2d_eval_greville_point (sp, gr, i, j+1, var) /* spline2d_eval (sp, u, v, var) */);
    }
  } 

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

  ccs_problem_lu_solve (css, rhs);

  for ( i = 0; i < NU+k-1; i++ ) {
    gint indexi = (i + k - 1)%(NU+k-1);
    for ( j = 0; j < NV-2; j++) {
      coeff_assign (sp, indexi, j+1, var, gsl_vector_get (rhs, i%NU+j*NU));
    }
  }

  gsl_vector_free (rhs);
  ccs_problem_destroy (css);
}

CCSProblem * spline2d_build_galerkin_fit_matrix_no_metric (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  size_t ustart, uend, vstart, vend;
  gint NU = sp->NXU;
  gint NV = sp->NXV;
  gint size = NU*NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  gsl_matrix * Bu, * Bv;

  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
    ustart = gp->istart_x;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bux, m);
	  
      for ( n = 0; n < ng; n++) {	    
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	    gint indexi = i + j*NU;
	      
	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*NU] +=  wmnija*gsl_matrix_get (Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

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
  }
  g_array_append_val (fit->index, count);

  return fit;
}

CCSProblem * spline2d_build_galerkin_fit_matrix_no_metric_bc_dirichlet (Spline2D * sp)
{
  gint i, j, m, n, ii, a, b;
  size_t ustart, uend, vstart, vend;
  gint NU = sp->NXU;
  gint NV = sp->NXV;
  gint size = NU*NV;
  gdouble A0[size][size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  gsl_matrix * Bu, * Bv;

  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->M*sp->N; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
    ustart = gp->istart_x;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bux, m);
	  
      for ( n = 0; n < ng; n++) {	    
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	    gint indexi = i + j*NU;
	      
	    for ( a = 0; a < sp->k; a++) {
	      gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
	      for ( b = 0; b < sp->k; b++) {
		A0[indexi][(ustart+a) + (vstart+b)*NU] +=  wmnija*gsl_matrix_get (Bv, b, 0);
	      }
	    }
	  }
	}
	
      }
    }
    
  }

#if 1
  // Dirichlet on edges
  gsl_vector * Bu2 = gsl_vector_alloc (sp->k);
  gsl_vector * Bv2 = gsl_vector_alloc (sp->k);
  GrevillePoints * gr = sp->gr;
  for ( i = 0; i < sp->NU; i++) {
    gint index1 = i;
    gdouble u = g_array_index (gr->ui, gdouble, i);
    gdouble v = 0.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_eval_nonzero (u, Bu2, &ustart, &uend, sp->w_u);
    gsl_bspline_eval_nonzero (v, Bv2, &vstart, &vend, sp->w_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_vector_get (Bu2, m);
      for ( n = 0; n < sp->k; n++) {
	gdouble cv = gsl_vector_get (Bv2, n);
	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }
    

    index1 = i + NU*(NV-1);
    v = 1.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_eval_nonzero (u, Bu2, &ustart, &uend, sp->w_u);
    gsl_bspline_eval_nonzero (v, Bv2, &vstart, &vend, sp->w_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_vector_get (Bu2, m);
      for ( n = 0; n < sp->k; n++) {
    	gdouble cv = gsl_vector_get (Bv2, n);
    	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }
  }

  for ( i = 0; i < sp->NV; i++) {
    gint index1 = i*sp->NU;
    gdouble v = g_array_index (gr->vj, gdouble, i);
    gdouble u = 0.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_eval_nonzero (u, Bu2, &ustart, &uend, sp->w_u);
    gsl_bspline_eval_nonzero (v, Bv2, &vstart, &vend, sp->w_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_vector_get (Bu2, m);
      for ( n = 0; n < sp->k; n++) {
  	gdouble cv = gsl_vector_get (Bv2, n);
  	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }
    
    index1 = NU-1 + i*sp->NU;
    u = 1.;
      
    for ( j = 0; j < size; j++)
      A0[j][index1] = 0.;

    gsl_bspline_eval_nonzero (u, Bu2, &ustart, &uend, sp->w_u);
    gsl_bspline_eval_nonzero (v, Bv2, &vstart, &vend, sp->w_v);

    for ( m = 0; m < sp->k; m++) {
      gdouble cu = gsl_vector_get (Bu2, m);
      for ( n = 0; n < sp->k; n++) {
    	gdouble cv = gsl_vector_get (Bv2, n);
    	A0[(ustart+m) + (vstart+n)*NU] [index1] = cu*cv;
      }
    }
  }
#endif

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
  }
  g_array_append_val (fit->index, count);

  return fit;
}

static gdouble bathy_x (gdouble u, gdouble v)
{
  return -4 + 6.*u;
  return -6 + 9.*u;
  return -10. + 15.*u;
  //return -15. + 20.*u;
}

static gdouble bathy_y (gdouble u, gdouble v)
{
  // return -3. + v*6;

  return v < 0.15 ? -3. + v/0.15*1.5 : v < 0.85 ? -1.5 + (v-0.15)/0.7*3. : 1.5 + (v-0.85)/0.15*1.5;
  
  //return -5. + 10*v;
  return -4. + 8*v;
}

static gsl_vector * rectangle_bathy_rhs_x (Spline2D * sp)
{
  gint i, j, m, n, ii;
  size_t ustart, uend, vstart, vend;
  gint size = sp->NU*sp->NV;

  gsl_vector * rhs_x = gsl_vector_alloc (size);

  gsl_vector_set_zero (rhs_x);

  g_assert (sp != NULL);
  
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
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	gsl_bspline_eval_nonzero (vn, Bv, &vstart, &vend, sp->w_v); // Maybe store that
	
	gdouble fmnx = bathy_x (um, vn);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_vector_get(Bu, i-ustart);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_vector_get(Bv, j-vstart);
	    gint indexi = i + j*sp->NU;
	 
	    gdouble tmpx = wmnij*fmnx + gsl_vector_get (rhs_x, indexi);
	    gsl_vector_set (rhs_x, indexi, tmpx);
	  }
	}
	
      }
    }
    
  }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  return rhs_x;
}

static gsl_vector * rectangle_bathy_rhs_y (Spline2D * sp)
{
  gint i, j, m, n, ii;
  size_t ustart, uend, vstart, vend;
  gint size = sp->NU*sp->NV;

  gsl_vector * rhs_y = gsl_vector_alloc (size);

  gsl_vector_set_zero (rhs_y);

  g_assert (sp != NULL);

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
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	gsl_bspline_eval_nonzero (vn, Bv, &vstart, &vend, sp->w_v); // Maybe store that
	
	gdouble fmnx = bathy_x (um, vn);
	gdouble fmny = bathy_y (um, vn);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_vector_get(Bu, i-ustart);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_vector_get(Bv, j-vstart);
	    gint indexi = i + j*sp->NU;
	 
	    gdouble tmpy = wmnij*fmny + gsl_vector_get (rhs_y, indexi);

	    gsl_vector_set (rhs_y, indexi, tmpy);
	  }
	}
	
      }
    }
    
  }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  return rhs_y;
}

Spline2D * rectangular_grid (gint M, gint N)
{
  Spline2D * grid =  spline2d_new (N-1, M-1, 3, 3, 3);

  coeff_set_var_to_zero (grid, 0);
  coeff_set_var_to_zero (grid, 1);
  coeff_set_var_to_zero (grid, 2);

  spline2d_init_panels (grid);

  gint size = gsl_bspline_ncoeffs (grid->w_u)*gsl_bspline_ncoeffs (grid->w_v);
  CCSProblem * A = spline2d_build_galerkin_fit_matrix_no_metric (grid);

  gsl_vector * rhs_x = rectangle_bathy_rhs_x (grid);
  gsl_vector * rhs_y = rectangle_bathy_rhs_y (grid);
  
  // LU decomposition
  ccs_problem_lu_solve (A, rhs_x);

  gint i, j;
  for ( i = 0; i < gsl_bspline_ncoeffs (grid->w_u); i++) {
    for ( j = 0; j < gsl_bspline_ncoeffs (grid->w_v); j++) {
      coeff_assign (grid, i, j, 0, gsl_vector_get (rhs_x, i + j*gsl_bspline_ncoeffs (grid->w_u)));
    }
  }

  ccs_problem_lu_solve (A, rhs_y);

  for ( i = 0; i < gsl_bspline_ncoeffs (grid->w_u); i++) {
    for ( j = 0; j < gsl_bspline_ncoeffs (grid->w_v); j++) {
      coeff_assign (grid, i, j, 1, gsl_vector_get (rhs_y, i + j*gsl_bspline_ncoeffs (grid->w_u)));
    }
  }

  spline2d_reinit_panels_physical_quantities (grid);

  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_x);
  ccs_problem_destroy (A);

  return grid;
}

static gsl_vector * parametric_bathy_rhs (Spline2D * sp, Spline2DFunc func, gpointer data)
{
  gint i, j, m, n, ii;
  size_t ustart, uend, vstart, vend;
  gint size = sp->NU*sp->NV;

  gsl_vector * rhs = gsl_vector_alloc (size);

  gsl_vector_set_zero (rhs);

  g_assert (sp != NULL);

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
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	gsl_bspline_eval_nonzero (vn, Bv, &vstart, &vend, sp->w_v); // Maybe store that
	
	gdouble fmn = func (sp, um, vn, data);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_vector_get(Bu, i-ustart);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_vector_get(Bv, j-vstart);
	    gint indexi = i + j*sp->NU;
	 
	    gdouble tmpy = wmnij*fmn + gsl_vector_get (rhs, indexi);

	    gsl_vector_set (rhs, indexi, tmpy);
	  }
	}
	
      }
    }
    
  }
  gsl_vector_free (Bu);
  gsl_vector_free (Bv);

  return rhs;
}

Spline2D * parametric_grid (gint M, gint N, GaussFunc func_x, GaussFunc func_y, gpointer data)
{
  Spline2D * grid =  spline2d_new (N-1, M-1, 3, 3, 3);

  grid->build_fit_matrix = spline2d_build_galerkin_fit_matrix_bc;
  grid->build_fit_rhs = build_galerkin_rhs_gauss_bc;

  coeff_set_var_to_zero (grid, 0);
  coeff_set_var_to_zero (grid, 1);
  coeff_set_var_to_zero (grid, 2);

  spline2d_init_panels (grid);

  gint size = grid->NU*grid->NV;

  CCSProblem * fit = spline2d_build_galerkin_fit_matrix_no_metric (grid);

  gsl_vector * rhs_x = build_galerkin_rhs_gauss_no_metric (grid, func_x, data, NULL, NULL);
  gsl_vector * rhs_y = build_galerkin_rhs_gauss_no_metric (grid, func_y, data, NULL, NULL);  
  
  // LU decomposition
  ccs_problem_lu_solve (fit, rhs_x);
  ccs_problem_lu_solve (fit, rhs_y);

  gint i, j;
  grid->copy_fit_solution (grid, rhs_x, 0);
  grid->copy_fit_solution (grid, rhs_y, 1);

  spline2d_reinit_panels_physical_quantities (grid);

  grid->fit = grid->build_fit_matrix (grid);

  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_x);
  ccs_problem_destroy (fit);

  return grid;
}

Spline2D * spline2d_parametric_patch (gint M, gint N, GaussFunc func_x,
				      GaussFunc func_y, GaussFunc func_z,
				      gpointer data, gint k, gint ninner,
				      gint nouter)
{
  Spline2D * grid =  spline2d_new (N, M, k, ninner, nouter);

  coeff_set_var_to_zero (grid, 0);
  coeff_set_var_to_zero (grid, 1);
  coeff_set_var_to_zero (grid, 2);

  spline2d_init_panels (grid);

  gint size = grid->NU*grid->NV;

  CCSProblem * fit = spline2d_build_galerkin_fit_matrix_no_metric (grid);

  gsl_vector * rhs_x = build_galerkin_rhs_gauss_no_metric (grid, func_x, data, NULL, NULL);
  gsl_vector * rhs_y = build_galerkin_rhs_gauss_no_metric (grid, func_y, data, NULL, NULL);
  gsl_vector * rhs_z = build_galerkin_rhs_gauss_no_metric (grid, func_z, data, NULL, NULL);
  
  // LU decomposition
  ccs_problem_lu_solve (fit, rhs_x);
  ccs_problem_lu_solve (fit, rhs_y);
  ccs_problem_lu_solve (fit, rhs_z);

  gint i, j;
  grid->copy_fit_solution (grid, rhs_x, 0);
  grid->copy_fit_solution (grid, rhs_y, 1);
  grid->copy_fit_solution (grid, rhs_z, 2);

  spline2d_reinit_panels_physical_quantities (grid);

  grid->fit = grid->build_fit_matrix (grid);

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_z);

  ccs_problem_destroy (fit);

  return grid;
}

Spline2D * spline2d_better_parametric_patch (gint M, gint N, GaussFunc func_x,
					     GaussFunc func_y, GaussFunc func_z,
					     gpointer data, gint k, gint ninner,
					     gint nouter)
{
  Spline2D * grid =  spline2d_new (N, M, k, ninner, nouter);

  // Maniar, 1995 appendix A p 194
  // Centripetal reparam
  
  // Set k first and last nodes

  // Fit

  coeff_set_var_to_zero (grid, 0);
  coeff_set_var_to_zero (grid, 1);
  coeff_set_var_to_zero (grid, 2);

  spline2d_init_panels (grid);

  gint size = grid->NU*grid->NV;

  CCSProblem * fit = spline2d_build_galerkin_fit_matrix_no_metric (grid);

  gsl_vector * rhs_x = build_galerkin_rhs_gauss_no_metric (grid, func_x, data, NULL, NULL);
  gsl_vector * rhs_y = build_galerkin_rhs_gauss_no_metric (grid, func_y, data, NULL, NULL);
  gsl_vector * rhs_z = build_galerkin_rhs_gauss_no_metric (grid, func_z, data, NULL, NULL);
  
  // LU decomposition
  ccs_problem_lu_solve (fit, rhs_x);
  ccs_problem_lu_solve (fit, rhs_y);
  ccs_problem_lu_solve (fit, rhs_z);

  gint i, j;
  grid->copy_fit_solution (grid, rhs_x, 0);
  grid->copy_fit_solution (grid, rhs_y, 1);
  grid->copy_fit_solution (grid, rhs_z, 2);

  spline2d_reinit_panels_physical_quantities (grid);

  grid->fit = grid->build_fit_matrix (grid);

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_z);

  ccs_problem_destroy (fit);

  return grid;
}

Spline2D * spline2d_parametric_periodic_patch (gint M, gint N, GaussFunc func_x,
					       GaussFunc func_y, GaussFunc func_z,
					       gpointer data)
{
  Spline2D * grid =  periodic_fs_new (N, M, 3, 4, 3);

  coeff_set_var_to_zero (grid, 0);
  coeff_set_var_to_zero (grid, 1);
  coeff_set_var_to_zero (grid, 2);

  grid->noflux = FALSE;
  //spline2d_init_panels (grid);
  periodic_fs_init_panels (grid);
  //grid->hull_patch = sim->hull->patches->data;

  gint size = grid->NU*grid->NV;

  CCSProblem * fit = periodic_fs_build_galerkin_fit_matrix_no_metric (grid);

  gsl_vector * rhs_x = periodic_fs_build_galerkin_rhs_gauss_no_metric (grid, func_x, data, NULL, NULL);
  gsl_vector * rhs_y = periodic_fs_build_galerkin_rhs_gauss_no_metric (grid, func_y, data, NULL, NULL);
  gsl_vector * rhs_z = periodic_fs_build_galerkin_rhs_gauss_no_metric (grid, func_z, data, NULL, NULL);
  
  // LU decomposition
  ccs_problem_lu_solve (fit, rhs_x);
  ccs_problem_lu_solve (fit, rhs_y);
  ccs_problem_lu_solve (fit, rhs_z);

  gint i, j;
  periodic_fs_copy_problem_solution_no_metric (grid, rhs_x, 0);
  periodic_fs_copy_problem_solution_no_metric (grid, rhs_y, 1);
  periodic_fs_copy_problem_solution_no_metric (grid, rhs_z, 2);

  //spline2d_reinit_panels_physical_quantities (grid);
  periodic_fs_reinit_panels_physical_quantities (grid);

  grid->fit = grid->build_fit_matrix (grid);

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_z);
  ccs_problem_destroy (fit);

  periodic_fs_numbering (grid);

  return grid;
}

Spline2D * parametric_grid2 (gint M, gint N, GaussFunc func_x, GaussFunc func_y, gpointer data)
{
  /* Spline2D * grid =  spline2d_new (N-1, M-1, 3, 3, 3); */
  Spline2D * grid =  periodic_fs_new (N-1, M-1, 3, 4, 3);

  /* coeff_set_var_to_zero (grid, 0); */
  /* coeff_set_var_to_zero (grid, 1); */
  /* coeff_set_var_to_zero (grid, 2); */

  /* spline2d_init_panels (grid); */
  periodic_fs_init_panels (grid);

  gint size = grid->NU*grid->NV;
  /* CCSProblem * A = spline2d_build_galerkin_fit_matrix2 (grid); */

  CCSProblem * fit = periodic_fs_build_galerkin_fit_matrix_no_metric (grid);
  gsl_vector * rhs_x = periodic_fs_build_galerkin_rhs_gauss_no_metric (grid, func_x, data, NULL, NULL);
  gsl_vector * rhs_y = periodic_fs_build_galerkin_rhs_gauss_no_metric (grid, func_y, data, NULL, NULL);
  // LU decomposition
  ccs_problem_lu_solve (fit, rhs_x);

  gint i, j;
  periodic_fs_copy_problem_solution_no_metric (grid, rhs_x, 0);

  ccs_problem_lu_solve (fit, rhs_y);

  periodic_fs_copy_problem_solution_no_metric (grid, rhs_y, 1);

  periodic_fs_reinit_panels_physical_quantities (grid);

  periodic_fs_numbering (grid);

  ccs_problem_destroy (fit);
  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_x);

  return grid;
}

Spline2D * parametric_grid3 (gint M, gint N, GrevilleFunc func_x, GrevilleFunc func_y, gpointer data)
{
  /* Spline2D * grid =  spline2d_new (N-1, M-1, 3, 3, 3); */
  Spline2D * grid =  periodic_fs_new (N, M, 3, 4, 3);

  /* coeff_set_var_to_zero (grid, 0); */
  /* coeff_set_var_to_zero (grid, 1); */
  /* coeff_set_var_to_zero (grid, 2); */

  /* spline2d_init_panels (grid); */
  periodic_fs_init_panels (grid);

  gint size = grid->NU*grid->NV;
  /* CCSProblem * A = spline2d_build_galerkin_fit_matrix2 (grid); */

  spline2d_fit_greville_periodic/* _border */ (grid, func_x, data, 0);
  spline2d_fit_greville_periodic/* _border */ (grid, func_y, data, 1);
  
  /* grid->noflux = FALSE; */
  /* CCSProblem * fit = grid->build_fit_matrix (grid); */
  /* gsl_vector * rhs_x = grid->build_fit_rhs (grid, func_x, data); */
  /* gsl_vector * rhs_y = grid->build_fit_rhs (grid, func_y, data); */
  /* grid->noflux = TRUE; */
  /* // LU decomposition */
  /* ccs_problem_lu_solve (fit, rhs_x); */

  /* gint i, j; */
  /* grid->copy_fit_solution (grid, rhs_x, 0); */

  /* ccs_problem_lu_solve (fit, rhs_y); */

  /* grid->copy_fit_solution (grid, rhs_y, 1); */

  periodic_fs_reinit_panels_physical_quantities (grid);

  /* ccs_problem_destroy (fit); */
  /* gsl_vector_free (rhs_y); */
  /* gsl_vector_free (rhs_x); */
  

  return grid;
}

gsl_vector * build_galerkin_rhs_gauss (Spline2D * sp,
				       GaussFunc func, gpointer data,
				       Spline2DFunc bc_func, gpointer bc_data,
				       gsl_vector * rhs)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble RHS[size];

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
	
	gdouble fmn = func (spp, m, n, data);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	    gint indexi = i + j*NU;
	 
	    RHS[indexi] += wmnij*fmn;
	  }
	}
	
      }
    }
  }

  // Copy lhs and rhs to gsl structures
  /* FILE * fp = fopen ("rhs-galerkin.tmp","w"); */
  //  gsl_vector * rhs = gsl_vector_alloc (size);
  /* gsl_vector *  */rhs = sp->rhs_vector (sp, /* sp-> */rhs);
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs, i, RHS[i]);
    /* fprintf (fp, "%i %f \n", i, RHS[i]); */
  }
  /* fclose (fp); */

  return rhs;
}

gsl_vector * build_galerkin_rhs_gauss_bc (Spline2D * sp,
					  GaussFunc func, gpointer data,
					  Spline2DFunc bc_func, gpointer bc_data,
					  gsl_vector * rhs)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint size = NU*NV;
  gdouble RHS[size];

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
	
	gdouble fmn = func (spp, m, n, data);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	    gint indexi = i + j*NU;
	 
	    RHS[indexi] += wmnij*fmn;
	  }
	}
	
      }
    }
  }

#if 1
  // Edges
  GrevillePoints * gr = sp->gr;
  for ( i = 0; i < sp->NU; i++) {
      gint index1 = i;
      gdouble u = g_array_index (gr->ui, gdouble, i);
      gdouble v = 0.;
      RHS[index1] = /* bc_func (sp, u, v, bc_data) */0.;
      RHS[i + sp->NU] = 0.;
      //RHS[i + 2*sp->NU] = 0.;

      index1 = i + sp->NU*(sp->NV-1);
      u = g_array_index (gr->ui, gdouble, i);
      v = 1.;
      RHS[index1] = /* bc_func (sp, u, v, bc_data) */0.;
      RHS[i + sp->NU*(sp->NV-2)] = 0.;
      //RHS[i + sp->NU*(sp->NV-3)] = 0.;
  }

  for ( i = 0; i < sp->NV; i++) {
    gint index1 = i*sp->NU;
    gdouble u = 0.;
    gdouble v = g_array_index (gr->vj, gdouble, i);
    RHS[index1] = /* bc_func (sp, u, v, bc_data) */0.;
    RHS[i*sp->NU + 1] = 0.;
    //RHS[i*sp->NU + 2] = 0.;

    index1 = sp->NU-1 + i*sp->NU;
    u = 0.;
    v = g_array_index (gr->vj, gdouble, i);
    RHS[index1] = /* bc_func (sp, u, v, bc_data) */0.;
    RHS[NU-2 + i*sp->NU] = 0.;
    //RHS[NU-3 + i*sp->NU] = 0.;
  }
#endif  

  // Copy lhs and rhs to gsl structures
  /* FILE * fp = fopen ("rhs-galerkin.tmp","w"); */
  //  gsl_vector * rhs = gsl_vector_alloc (size);
  /* gsl_vector *  */rhs = sp->rhs_vector (sp, /* sp-> */rhs);
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs, i, RHS[i]);
    /* fprintf (fp, "%i %f \n", i, RHS[i]); */
  }
  /* fclose (fp); */

  return rhs;
}

/* gsl_vector * build_greville_rhs (Spline2D * sp, */
/* 				 GaussFunc func, gpointer data, */
/* 				 Spline2DFunc bc_func, gpointer bc_data) */
/* { */
/*   gint i, j, m, n, ii; */
/*   size_t ustart, vstart; */
/*   gint NU = sp->NU; */
/*   gint NV = sp->NV; */
/*   gint size = NU*NV; */
/*   gdouble RHS[size]; */

/*   g_assert (sp != NULL); */

/*   // Initializes the coefficients to 0 */
/*   for ( i = 0; i < size; i++) { */
/*     RHS[i] = 0.; */
/*   } */
  
/*   gsl_matrix * Bu, * Bv; */
/*   // Loop over the panels of the patch */
/*   for ( ii = 0; ii < sp->panels->len; ii++) { */
/*     SPPanel * spp = g_ptr_array_index (sp->panels, ii); */
/*     g_assert (spp != NULL); */
     	
/*     // Gauss outer-integration */
/*     GaussPoints * gp = spp->outer; */
/*     gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule */
/*     ustart = gp->istart; */
/*     vstart = gp->jstart; */
/*     for ( m = 0; m < ng; m++) { */
/*       Bu = g_ptr_array_index (gp->Bu, m); */
	  
/*       for ( n = 0; n < ng; n++) { */
/* 	gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng); */
/* 	Bv = g_ptr_array_index (gp->Bv, n); */
	
/* 	gdouble fmn = func (spp, m, n, data); */

/* 	// Loop over the splines whose support is included in the panel */
/* 	for ( i = ustart; i < ustart + sp->k; i++) { */
/* 	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0); */
/* 	  for ( j = vstart; j < vstart + sp->k; j++) { */
/* 	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0); */
/* 	    gint indexi = i + j*NU; */
	 
/* 	    RHS[indexi] += wmnij*fmn; */
/* 	  } */
/* 	} */
	
/*       } */
/*     } */
/*   } */

/*   // Copy lhs and rhs to gsl structures */
/*   gsl_vector * rhs = sp->rhs_vector (sp); */
/*   for ( i = 0; i < size; i++) { */
/*     gsl_vector_set (rhs, i, RHS[i]); */
/*   } */

/*   return rhs; */
/* } */

gsl_vector * build_galerkin_rhs_gauss_no_metric (Spline2D * sp,
						 GaussFunc func,
						 gpointer data,
						 Spline2DFunc bc_func,
						 gpointer bc_data)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NXU;
  gint NV = sp->NXV;
  gint size = NU*NV;
  gdouble RHS[size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_matrix * Bu;
  gsl_matrix * Bv;   
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    ustart = gp->istart_x;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bux, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);
	
	gdouble fmn = func (spp, m, n, data);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	    gint indexi = i + j*NU;
	 
	    RHS[indexi] += wmnij*fmn;
	  }
	}
	
      }
    }
  }

  // Copy lhs and rhs to gsl structures
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, RHS[i]);

  return rhs;
}

gsl_vector * build_galerkin_rhs_gauss_no_metric_bc_dirichlet (Spline2D * sp,
							      GaussFunc func,
							      gpointer data,
							      Spline2DFunc bc_func,
							      gpointer bc_data)
{
  gint i, j, m, n, ii;
  size_t ustart, vstart;
  gint NU = sp->NXU;
  gint NV = sp->NXV;
  gint size = NU*NV;
  gdouble RHS[size];

  g_assert (sp != NULL);

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_matrix * Bu;
  gsl_matrix * Bv;   
  // Loop over the panels of the patch
  for ( ii = 0; ii < sp->panels->len; ii++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, ii);
    g_assert (spp != NULL);
     	
    // Gauss outer-integration
    GaussPoints * gp = spp->outer;
    gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
    ustart = gp->istart_x;
    vstart = gp->jstart;
    for ( m = 0; m < ng; m++) {
      Bu = g_ptr_array_index (gp->Bux, m);
	  
      for ( n = 0; n < ng; n++) {
	gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	Bv = g_ptr_array_index (gp->Bv, n);
	
	gdouble fmn = func (spp, m, n, data);

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	    gint indexi = i + j*NU;
	 
	    RHS[indexi] += wmnij*fmn;
	  }
	}
	
      }
    }
  }

#if 1
  // Edges
  GrevillePoints * gr = sp->gr;
  for ( i = 0; i < sp->NU; i++) {
      gint index1 = i;
      gdouble u = g_array_index (gr->ui, gdouble, i);
      gdouble v = 0.;
      RHS[index1] = bc_func (sp, u, v, bc_data);

      index1 = i + sp->NU*(sp->NV-1);
      u = g_array_index (gr->ui, gdouble, i);
      v = 1.;
      RHS[index1] = bc_func (sp, u, v, bc_data);
  }

  for ( i = 0; i < sp->NV; i++) {
    gint index1 = i*sp->NU;
    gdouble u = 0.;
    gdouble v = g_array_index (gr->vj, gdouble, i);
    RHS[index1] = bc_func (sp, u, v, bc_data);

    index1 = sp->NU-1 + i*sp->NU;
    u = 0.;
    v = g_array_index (gr->vj, gdouble, i);
    RHS[index1] = bc_func (sp, u, v, bc_data);
  }
#endif  

  // Copy lhs and rhs to gsl structures
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, RHS[i]);

  return rhs;
}

void apply_dirichlet_conditions (Spline2D * sp, GaussFunc func, gpointer data, gint var)
{
  gint i, j;

  g_assert (sp != NULL);

  // Get matrix if not already stored
  if (!sp->fit)
    sp->fit = sp->build_fit_matrix (sp);

  // Get rhs
  sp->rhs = sp->build_fit_rhs (sp, func, data, NULL, NULL, sp->rhs);

  // Solve problem
  ccs_problem_lu_solve (sp->fit, sp->rhs);

  /* FILE * fp = fopen ("sol.tmp","w"); */
  /* for ( i = 0; i < rhs->size; i++) */
  /*   fprintf (fp, "%i %f\n", i, gsl_vector_get (rhs, i)); */
  /* fclose (fp); */

  // Copy results to patch
  sp->copy_fit_solution (sp, sp->rhs, var);

  //  gsl_vector_free (rhs);
}

static gsl_vector * apply_neumann_conditions_rhs (Spline2D * sp, NeumannFunc func, gpointer data)
{
  gint i, j, m, n, ii;
  size_t ustart, uend, vstart, vend;
  gint size = sp->NU*sp->NV;
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

	// Loop over the splines whose support is included in the panel
	for ( i = ustart; i < ustart + sp->k; i++) {
	  gdouble wmni = wmn*gsl_vector_get(Bu, i-ustart);
	  for ( j = vstart; j < vstart + sp->k; j++) {
	    gdouble wmnij = wmni*gsl_vector_get(Bv, j-vstart);
	    gint indexi = i + j*sp->NU;
	 
	    Vector N = spline2d_normal (sp, um, vn);
	    Vector V = func (sp, um, vn, data);
	    RHS[indexi] += wmnij*vector_scalar_product (&V, &N);
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

void apply_neumann_conditions (Spline2D * sp, NeumannFunc func, gpointer data, gint var)
{
  gint size = sp->NU*sp->NV;
  gint i, j;

  g_assert (sp != NULL);

  if (!sp->fit)
    sp->fit = sp->build_fit_matrix (sp);

  // Get rhs
  gsl_vector * rhs = apply_neumann_conditions_rhs (sp, func, data);

  // Solve the linear system
  ccs_problem_lu_solve (sp->fit, rhs);
  
  // Copy results to patch
  sp->copy_fit_solution (sp, rhs, var);

  gsl_vector_free (rhs);
}

void spline_set_var_to_constant (Spline2D * sp, gint var, gdouble val)
{
  gint i, j;

  for ( i = 0; i < gsl_bspline_ncoeffs (sp->w_u); i++) {
    for ( j = 0; j < gsl_bspline_ncoeffs (sp->w_v); j++) {
      coeff_assign (sp, i, j, var, val);
    }
  }
}

/**
 * Build the rhs for the neumann boundary problem
 * using the normal derivative stored in var.
 **/
void boundary_subproblem_build_rhs_neumann (Spline2D * sp1, Spline2D * sp2,
					    BoundarySubProblem * bsp,
					    gint var)
{
  gint i, j;
  BoundaryProblem * dirichlet = bsp->dirichlet;
  BoundaryProblem * neumann = bsp->neumann;

  gint NU1 = sp1->periodic ? sp1->NUT:sp1->NU;
  gint NU2 = sp2->periodic ? sp2->NUT:sp2->NU;

  gdouble tmp[NU2*sp2->NV];
  gdouble rhs[NU1*sp1->NV];

  gsl_vector_set_zero (neumann->rhs);

  for ( i = 0; i < NU1*sp1->NV; i++)
    rhs[i] = 0.;

  // Loop over the b-splines
  if (!sp2->periodic) {
    for ( i = 0; i < sp2->NU; i++) {
      for ( j = 0; j < sp2->NV; j++) {
	tmp[i + j*sp2->NU] = coeff (sp2, i, j, var);
      }
    }
  }
  else {
    Spline2D * sp22 = sp2;
    gint NUT = sp22->NUT;
    while (sp22) {
      for ( i = 0; i < sp22->NU; i++) {
	gint indexi = (sp22->fs_index + i)%NUT;
      	for ( j = 0; j < sp22->NV; j++) {
      	  tmp[indexi+j*NUT] = coeff (sp22, i, j, var);
      	}
      }
      sp22 = sp22->next;
    }
  }

  for ( i = 0; i < dirichlet->A->size1; i++) {
    for (j = 0; j < dirichlet->A->size2; j++) {
      rhs[i] += gsl_matrix_get (dirichlet->A, i, j)*tmp[j];
    }
  } 

  for ( i = 0; i < NU1*sp1->NV; i++)
    gsl_vector_set (neumann->rhs, i, rhs[i]);
}

void boundary_subproblem_build_rhs_dirichlet (Spline2D * sp1, Spline2D * sp2,
					      BoundarySubProblem * bsp,
					      gint var)
{
  gint i, j;
  BoundaryProblem * dirichlet = bsp->dirichlet;
  BoundaryProblem * neumann = bsp->neumann;

  gint NU1 = sp1->periodic ? sp1->NUT:sp1->NU;
  gint NU2 = sp2->periodic ? sp2->NUT:sp2->NU;

  gdouble tmp[NU2*sp2->NV];
  gdouble rhs[NU1*sp1->NV];

  gsl_vector_set_zero (dirichlet->rhs);

  for ( i = 0; i < NU1*sp1->NV; i++)
    rhs[i] = 0.;

  // Loop over the b-splines
  if (!sp2->periodic) {
    for ( i = 0; i < sp2->NU; i++) {
      for ( j = 0; j < sp2->NV; j++) {
	tmp[i + j*sp2->NU] = coeff (sp2, i, j, var);
      }
    }
  }
  else {
    Spline2D * sp22 = sp2;
    gint NUT = sp22->NUT;
    while (sp22) {
      for ( i = 0; i < sp22->NU; i++) {
	gint indexi = (sp22->fs_index + i)%NUT;
      	for ( j = 0; j < sp22->NV; j++) {
      	  tmp[indexi+j*NUT] = coeff (sp22, i, j, var);
      	}
      }
      sp22 = sp22->next;
    }
  }

  for ( i = 0; i < neumann->A->size1; i++) {
    for (j = 0; j < neumann->A->size2; j++) {
      rhs[i] += gsl_matrix_get (neumann->A, i, j)*tmp[j];
    }
  }

  for ( i = 0; i < NU1*sp1->NV; i++)
    gsl_vector_set (dirichlet->rhs, i, rhs[i]);
}

/***** Free Surface *******/

//  0: x
//  1: y
//  2: z
//  3: Phi (basis flow potential)
//  4: Phin
//  5: phi (time-local flow potential)
//  6: phin
//  7: Phi2 (disturbance flow potential)
//  8: Phi2n
//  9: zeta (surface elevation perturbation)
//  10: zetan
//  11: Phidt (Phi at previous time step)
//  12: Temporary variable
//  13: 
//  14: 
//  15: 
//  16: 

gboolean gsl_problem_solve_using_SOR (gsl_matrix * Aij,
				      gsl_vector * rhsi,
				      gsl_vector * lhsi,
				      gdouble tolerance,
				      gboolean verbose)
{
  gdouble err = 1;
  gdouble w = 1.;
  gint i, j;
  gdouble rhs[rhsi->size];
  gdouble lhs[lhsi->size];
  gdouble A[Aij->size1][Aij->size2];

  for ( i = 0; i < rhsi->size; i++)
    rhs[i] = gsl_vector_get (rhsi, i);

  for ( i = 0; i < lhsi->size; i++)
    lhs[i] = gsl_vector_get (lhsi, i);

  for ( i = 0; i < Aij->size1; i++)
    for ( j = 0; j < Aij->size2; j++)
      A[j][i] = gsl_matrix_get (Aij, i, j);

  gdouble tmp[lhsi->size];

  gdouble error = 2*tolerance;
  gint niter = 0;

  while (error > tolerance && niter < 1000) {
    error = 0.;
    for ( i = 0; i < lhsi->size; i++) {
      tmp[i] = rhs[i];
      for ( j = 0; j < Aij->size2; j++) {
  	if ( i != j )
  	  tmp[i] = -A[i][j]*lhs[i];
      }
      tmp[i] /= A[i][i];
    }
    
    for ( i = 0; i < lhsi->size; i++) {
      if (fabs(lhs[i]-tmp[i]) > error)
    	error = fabs(lhs[i] - tmp[i]);
    }

    for ( i = 0; i < lhsi->size; i++)
      lhs[i] = tmp[i];

    niter++;
    if (verbose)
      fprintf(stdout, "Jacobi error %e after %i iterations \n", error, niter);
  }

 


  /* while (error > tolerance) { */
  /*   error = 0; */
  /*   for ( i = 0; i < rhsi->size; i++) { */
  /*     gdouble sum = rhs[i]; */
      
  /*     for (j = 0; j < rhsi->size; j++) { */
  /*   	if (i != j) { */
  /*   	  if (j > i) */
  /* 	    sum -= A[i][j]*lhs[j]; */
  /* 	  else */
  /* 	    sum -= A[i][j]*tmp[j]; */
  /* 	} */
  /*     } */
  /*     tmp[i] = (1. - w)*lhs[i] + w/A[i][i]*sum; */
  /*   } */
    
  /*   for ( i = 0; i < lhsi->size; i++) { */
  /*     if (fabs(lhs[i] - tmp[i]) > err) */
  /*   	error = fabs(lhs[i] - tmp[i]); */
  /*   } */

  /*   for ( i = 0; i < lhsi->size; i++) */
  /*     lhs[i] = tmp[i]; */
  /*   if (verbose) */
  /*     fprintf(stdout, "SOR error %e after %i iterations \n", error, niter); */
  /*   niter++; */
  /* } */

  /* g_array_free (x0, FALSE); */

  /* fprintf(stdout, "SOR error %e after %i iterations \n", err, niter); */

  for ( i = 0; i < lhsi->size; i++)
    gsl_vector_set (lhsi, i, lhs[i]);

  return TRUE;
}

void boundary_problem_direct_sparse_solve (gdouble ** A,
					   gdouble * RHS,
					   gint size)
{
  /* Storage in Compressed Column Storage (CCS) format for superli */
  GArray * matrix = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray * column = g_array_new (FALSE, FALSE, sizeof (gint));
  GArray * index = g_array_new (FALSE, FALSE, sizeof (gint));
  gint i, j, count = 0;

  for ( i = 0; i < size; i++) {
 
   g_array_append_val (index, count);
    for ( j = 0; j < size; j++) {
      if (A[i][j] != 0.) {
  	g_array_append_val (matrix, A[i][j]);
  	g_array_append_val (column, j);
  	count++;
      }
    }
  }
  g_array_append_val (index, count);

  SuperMatrix A2, L2, U2, B2;
  gint * perm_r = g_malloc (size*sizeof(gint));
  gint * perm_c = g_malloc (size*sizeof(gint));
  gint info;
  superlu_options_t options;
  SuperLUStat_t stat;

  // Create superlu matrix
  dCreate_CompCol_Matrix (&A2, size, size, matrix->len, (double *) matrix->data, (int *) column->data, (int *) index->data, SLU_NC, SLU_D, SLU_GE);

  // Create superlu rhs
  dCreate_Dense_Matrix (&B2, size, 1, &RHS[0], size, SLU_DN, SLU_D, SLU_GE);

  // Solver options are defaults
  set_default_options (&options);
  options.ColPerm = NATURAL;

  // Initialize stats
  StatInit(&stat);

  // Solve system using LU decomposition
  dgssv (&options, &A2, perm_c, perm_r, &L2, &U2, &B2, &stat, &info);
  
  // Check for success
  g_assert (info == 0);

  // Here is the solution
  double *sol = (double*) ((DNformat*) B2.Store)->nzval;

  for ( i = 0; i < size; i++)
    RHS[i] = sol[i];
  
  /* // Copy the solution to the patch */
  /* for ( i = 0; i < gsl_bspline_ncoeffs (sp->w_u); i++) { */
  /*   for ( j = 0; j < gsl_bspline_ncoeffs (sp->w_v); j++) { */
  /*     coeff_assign (sp, i, j, 7, sol [ i + j*gsl_bspline_ncoeffs (sp->w_u)]); */
  /*   } */
  /* } */
  
  // Free superlu structures
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix (&A2);
  Destroy_SuperMatrix_Store (&B2);
  Destroy_SuperNode_Matrix (&L2);
  Destroy_CompCol_Matrix (&U2);
  StatFree (&stat);
}

/**** CUDA ****/

#if CUDA
void cuda_init (CUdevice * dev, CUcontext * context)
{
  /* CUdevice  dev; */
  /* CUcontext context; */
  if( CUDA_SUCCESS != cuInit( 0 ) ) {
    fprintf(stderr, "CUDA: Not initialized\n" ); exit(-1);
  }
  if( CUDA_SUCCESS != cuDeviceGet( dev, 0 ) ) {
    fprintf(stderr, "CUDA: Cannot get the device\n"); exit(-1);
  }
  if( CUDA_SUCCESS != cuCtxCreate( context, 0, *dev ) ) {
    fprintf(stderr, "CUDA: Cannot create the context\n"); exit(-1);
  }
  if( CUBLAS_STATUS_SUCCESS != cublasInit( ) ) {
    fprintf(stderr, "CUBLAS: Not initialized\n"); exit(-1);
  }
  printout_devices( );
}

void cuda_finalize (CUcontext * context)
{
  cuCtxDetach( *context );
  cublasShutdown();
}


CudaLUProblem * cuda_lu_problem_new (gint M, gint N)
{
  CudaLUProblem * new = g_malloc ( sizeof (CudaLUProblem));

  g_assert ( cudaSuccess == cudaMalloc ( (void **) &new->A, (M*N) * sizeof (double)));
  g_assert ( cudaSuccess == cudaMalloc ( (void **) &new->rhs, (M) * sizeof (double)));
  new->ipiv = g_malloc ( min(M,N) * sizeof (gint));

  return new;
}

void cuda_lu_problem_destroy (CudaLUProblem * cp)
{
  cudaFree (cp->A);
  cudaFree (cp->rhs);
}

CudaLUProblem * cuda_lu_factorise (gsl_matrix * A)
{
  /* Matrix size */
  magma_int_t M = A->size1, N = A->size2;
  CudaLUProblem * cp = cuda_lu_problem_new (M, N);
  magma_int_t info;
  
  // Copy matrix A on the GPU
  g_assert (CUBLAS_STATUS_SUCCESS == cublasSetMatrix( M, N, sizeof (double), A->data, A->size1, cp->A, M));

  // Computes LU factorisation
  magma_dgetrf_gpu( M, N, cp->A, M, cp->ipiv, &info);
  g_assert (info == MAGMA_SUCCESS);

  return cp;
}

void cuda_solve_factorised_lu (CudaLUProblem * cp, gsl_vector * rhs)
{
  magma_int_t info;

  // Copy rhs on the GPU
  g_assert (CUBLAS_STATUS_SUCCESS == cublasSetMatrix( rhs->size, 1, sizeof (double), rhs->data, rhs->size, cp->rhs, rhs->size));

  // Solves A x = b using the partial pivoting
  magma_dgetrs_gpu (MagmaNoTrans, rhs->size, 1, cp->A, rhs->size, cp->ipiv, cp->rhs, rhs->size, &info);
  g_assert (info == MAGMA_SUCCESS);

  // Get results from GPU
  g_assert (CUBLAS_STATUS_SUCCESS == cublasGetMatrix (rhs->size, 1, sizeof (double), cp->rhs, rhs->size, rhs->data, rhs->size));  
}

void cuda_lu_solve (gsl_matrix * A, gsl_vector * rhs)
{
  CudaLUProblem * cp = cuda_lu_factorise (A);

  cuda_solve_factorised_lu (cp, rhs);

  cuda_lu_problem_destroy (cp);
}

CudaSVDProblem * cuda_svd_problem_new (gint M, gint N)
{
  CudaSVDProblem * new = g_malloc ( sizeof (CudaSVDProblem));

  new->s = g_malloc (min(M,N)*sizeof(double));

  g_assert ( cudaSuccess ==  cudaMalloc( (void **) &new->drhs,
					 M*sizeof(double) ) );
  g_assert ( cudaSuccess ==  cudaMalloc( (void **) &new->du,
					 M*N*sizeof(double) ) );
  g_assert ( cudaSuccess ==  cudaMalloc( (void **) &new->dvt,
					 M*N*sizeof(double) ) );
  g_assert ( cudaSuccess ==  cudaMalloc( (void **) &new->dw,
					 M*sizeof(double) ) );

  new->M = M;
  new->N = N;

  new->w = gsl_vector_alloc (M);

  return new;
}

void cuda_svd_problem_destroy (CudaSVDProblem * cp)
{
  gsl_vector_free (cp->w);
  g_free (cp->u);
  g_free (cp->vt);
  g_free (cp->s);
  cudaFree (cp->du);
  cudaFree (cp->dvt);
  cudaFree (cp->drhs);
  cudaFree (cp->dw);
  g_free (cp);
}

CudaSVDProblem * cuda_svd_factorise (gsl_matrix * A)
{
  /* Matrix size */
  magma_int_t M = A->size1, N = A->size2;
  magma_int_t min_mn = min(M, N);
  magma_int_t nb = magma_get_dgebrd_nb (N);
  magma_int_t lwork = max(5*min_mn, (3*min_mn + max(M,N)))*nb;
  CudaSVDProblem * cp = cuda_svd_problem_new (M, N);
  magma_int_t info;
  double * hwork;
  double * u = g_malloc (M*M*sizeof(double));
  double * vt = g_malloc (N*N*sizeof(double));

  // Allocates workspace on GPU
  g_assert ( cudaSuccess == cudaMallocHost( (void**) &hwork, lwork*sizeof (double)));

  // Computes svd decomposition
  magma_dgesvd ('A', 'A', M, N,
  		A->data, M, cp->s,
  		u, M, vt, N,
  		hwork, lwork,
  		&info );

  // Copies results back on GPU for future use
  cublasSetMatrix (cp->M, cp->N, sizeof (double), u, cp->M,
		   cp->du, cp->M);
  cublasSetMatrix (cp->M, cp->N, sizeof (double), vt, cp->M,
		   cp->dvt, cp->M);

  // Free work-space
  cudaFree (hwork);
  g_free (u);
  g_free (vt);

  g_assert (info == MAGMA_SUCCESS);

  return cp;
}

/**  Solves the system A x = b using the SVD factorization
 *
 *  A = U S V^T
 *
 *  to obtain x. For M x N systems it finds the solution in the least
 *  squares sense.
 **/
void cuda_solve_factorised_svd (CudaSVDProblem * cp, gsl_vector * rhs)
{
   gint i;

   gsl_vector_set_zero (cp->w);
   
   // Copy rhs and w to GPU
   cublasSetVector (rhs->size, sizeof (double),
		    rhs->data, (int) rhs->stride,
		    cp->drhs, (int) rhs->stride);
   cublasSetVector( cp->w->size, sizeof (double),
		    cp->w->data, (int) cp->w->stride,
		    cp->dw, (int) cp->w->stride);


   // w = 1.0*U^T*rhs + 0*w
   cublasDgemv ('T', rhs->size, rhs->size, 1.0, cp->du, rhs->size,
		cp->drhs, (int) rhs->stride,
		0.0, cp->dw, (int) cp->w->stride);
   
   // Get result from GPU
   cublasGetVector (cp->w->size, sizeof (double),
		    cp->dw, (int) cp->w->stride,
		    cp->w->data, (int) cp->w->stride);
   
   // Singular values (diagonal matrix)
   for (i = 0; i < rhs->size; i++) {
     double wi = gsl_vector_get (cp->w, i);
     double alpha = *(cp->s + i);
     if (alpha != 0.)
       alpha = 1.0 / alpha;
     gsl_vector_set (cp->w, i, alpha * wi);
   }

   // Copies w to GPU
   cublasSetVector (cp->w->size, sizeof (double),
		    cp->w->data, (int) cp->w->stride,
		    cp->dw, (int) cp->w->stride);

   // rhs = 1.0*vt*w + 0*rhs
   cublasDgemv ('T', rhs->size, rhs->size, 1.0, cp->dvt, rhs->size,
		cp->dw, (int) cp->w->stride,
		0.0, cp->drhs, (int) rhs->stride );
    
   // Get results from GPU
   cublasGetVector (rhs->size, sizeof (double),
		    cp->drhs, (int) rhs->stride,
		    rhs->data, (int) rhs->stride);
}

void cuda_svd_solve (gsl_matrix * A, gsl_vector * rhs)
{
  CudaSVDProblem * cp = cuda_svd_factorise (A);

  cuda_solve_factorised_svd (cp, rhs);

  cuda_svd_problem_destroy (cp);
}

#endif

/**** OpenCL ****/

#if OPENCL
void opencl_init (/* CUdevice * dev, CUcontext * context */)
{
  magma_init ();
}

void opencl_finalize ()
{
  magma_finalize ();
}

OpenCLLUProblem * opencl_lu_problem_new (gint M, gint N)
{
  magma_err_t err;
  OpenCLLUProblem * new = g_malloc ( sizeof (OpenCLLUProblem));

  err = magma_get_devices (&new->device, 1, &new->num );
  if ( err != 0 || new->num < 1 ) {
    fprintf( stderr, "magma_get_devices failed: %d\n", err );
    exit(-1);
  }

  err = magma_queue_create( new->device, &new->queue );
  if ( err != 0 ) {
    fprintf( stderr, "magma_queue_create failed: %d\n", err );
    exit(-1);
  }

  g_assert ( MAGMA_SUCCESS == magma_malloc((magma_ptr *) &new->A, (M*N) * sizeof (double)));
  g_assert ( MAGMA_SUCCESS == magma_malloc((magma_ptr *) &new->rhs, (M) * sizeof (double)));
  new->ipiv = g_malloc ( min(M,N) * sizeof (gint));

  return new;
}

void opencl_lu_problem_destroy (OpenCLLUProblem * cp)
{
  magma_free (cp->A);
  magma_free (cp->rhs);
  magma_queue_destroy (cp->queue);
  free (cp);
}

OpenCLLUProblem * opencl_lu_factorise (gsl_matrix * A)
{
  /* Matrix size */
  magma_int_t M = A->size1, N = A->size2;
  OpenCLLUProblem * cp = opencl_lu_problem_new (M, N);
  magma_int_t info;
  
  // Copy matrix A on the GPU
  g_assert (MAGMA_SUCCESS == magma_dsetmatrix( M, N, A->data, 0, A->size1, cp->A, 0, M, cp->queue));

  // Computes LU factorisation
  magma_dgetrf_gpu( M, N, cp->A, 0, M, cp->ipiv, &info, cp->queue);
  g_assert (info == MAGMA_SUCCESS);

  return cp;
}

void opencl_solve_factorised_lu (OpenCLLUProblem * cp, gsl_vector * rhs)
{
  magma_int_t info;

  // Copy rhs on the GPU
  //  g_assert (CUBLAS_STATUS_SUCCESS == cublasSetMatrix( rhs->size, 1, sizeof (double), rhs->data, rhs->size, cp->rhs, rhs->size));
  g_assert (MAGMA_SUCCESS == magma_dsetmatrix( rhs->size, 1, rhs->data, 0, rhs->size, cp->rhs, 0, rhs->size, cp->queue));

  // Solves A x = b using the partial pivoting
  magma_dgetrs_gpu (MagmaNoTrans, rhs->size, 1, cp->A, 0, rhs->size, cp->ipiv, cp->rhs, 0, rhs->size, &info, cp->queue);
  g_assert (info == MAGMA_SUCCESS);

  // Get results from GPU
  /* g_assert (CUBLAS_STATUS_SUCCESS == cublasGetMatrix (rhs->size, 1, sizeof (double), cp->rhs, rhs->size, rhs->data, rhs->size)); */
  g_assert (MAGMA_SUCCESS == magma_dgetmatrix (rhs->size, 1, cp->rhs, 0, rhs->size, rhs->data, 0, rhs->size, cp->queue));
  
}

void opencl_lu_solve (gsl_matrix * A, gsl_vector * rhs)
{
  OpenCLLUProblem * cp = opencl_lu_factorise (A);

  opencl_solve_factorised_lu (cp, rhs);

  opencl_lu_problem_destroy (cp);
}

OpenCLSVDProblem * opencl_svd_problem_new (gint M, gint N)
{
  magma_err_t err;
  OpenCLSVDProblem * new = g_malloc ( sizeof (OpenCLSVDProblem));

  err = magma_get_devices (&new->device, 1, &new->num );
  if ( err != 0 || new->num < 1 ) {
    fprintf( stderr, "magma_get_devices failed: %d\n", err );
    exit(-1);
  }

  err = magma_queue_create( new->device, &new->queue );
  if ( err != 0 ) {
    fprintf( stderr, "magma_queue_create failed: %d\n", err );
    exit(-1);
  }

  new->s = g_malloc (min(M,N)*sizeof(double));

  g_assert ( MAGMA_SUCCESS == magma_malloc((magma_ptr *) &new->drhs, M* sizeof (double)));
  g_assert ( MAGMA_SUCCESS == magma_malloc((magma_ptr *) &new->du, M*N*sizeof (double)));
  g_assert ( MAGMA_SUCCESS == magma_malloc((magma_ptr *) &new->dvt, M*N*sizeof (double)));
  g_assert ( MAGMA_SUCCESS == magma_malloc((magma_ptr *) &new->dw, M*sizeof (double)));

  new->M = M;
  new->N = N;

  new->w = gsl_vector_alloc (M);

  return new;
}

void opencl_svd_problem_destroy (OpenCLSVDProblem * cp)
{
  gsl_vector_free (cp->w);
  g_free (cp->u);
  g_free (cp->vt);
  g_free (cp->s);
  magma_free (cp->du);
  magma_free (cp->dvt);
  magma_free (cp->drhs);
  magma_free (cp->dw);
  magma_queue_destroy (cp->queue);
  g_free (cp);
}

OpenCLSVDProblem * opencl_svd_factorise (gsl_matrix * A)
{
  /* Matrix size */
  magma_int_t M = A->size1, N = A->size2;
  magma_int_t min_mn = min(M, N);
  magma_int_t nb = magma_get_dgebrd_nb (N);
  magma_int_t lwork = max(5*min_mn, (3*min_mn + max(M,N)))*nb;
  OpenCLSVDProblem * cp = opencl_svd_problem_new (M, N);
  magma_int_t info;
  double * hwork;
  double * u = g_malloc (M*M*sizeof(double));
  double * vt = g_malloc (N*N*sizeof(double));

  // Allocates workspace on GPU
  g_assert ( MAGMA_SUCCESS == magma_malloc_host ((void**) &hwork, lwork*sizeof (double)));

  // Computes svd decomposition
  magma_dgesvd ('A', 'A', M, N,
  		A->data, M, cp->s,
  		u, M, vt, N,
  		hwork, lwork,
  		&info, cp->queue );

  // Copies results back on GPU for future use
  magma_dsetmatrix(cp->M, cp->N, u, 0, cp->M, cp->du, 0, cp->M,
		   cp->queue)
  magma_dsetmatrix(cp->M, cp->N, vt, 0, cp->M, cp->dvt, 0, cp->M,
		   cp->queue)

  // Free work-space
  magma_free (hwork);
  g_free (u);
  g_free (vt);

  g_assert (info == MAGMA_SUCCESS);

  return cp;
}

/**  Solves the system A x = b using the SVD factorization
 *
 *  A = U S V^T
 *
 *  to obtain x. For M x N systems it finds the solution in the least
 *  squares sense.
 **/
void cuda_solve_factorised_svd (CudaSVDProblem * cp, gsl_vector * rhs)
{
   gint i;

   gsl_vector_set_zero (cp->w);
   
   // Copy rhs and w to GPU
   magma_dsetvector (rhs->size, rhs->data, 0, (int) rhs->stride,
		     cp->drhs, 0, (int) rhs->stride, cp->queue);
   magma_dsetvector (cp->w->size, cp->w->data, 0, (int) cp->w->stride,
		     cp->dw, 0, (int) cp->w->stride, cp->queue);


   // w = 1.0*U^T*rhs + 0*w
   magmablas_dgmev ('T', rhs->size, rhs->size, 1.0, cp->du, 0,
		    rhs->size, cp->drhs, 0, (int) rhs->stride,
		    0.0, cp->dw, 0, (int) cp->w->stride, cp->queue);
   
   // Get result from GPU
   magma_dsetvector (cp->w->size, cp->dw, 0, (int) cp->w->stride,
		     cp->w->data, 0, (int) cp->w->stride, cp->queue);
   
   // Singular values (diagonal matrix)
   for (i = 0; i < rhs->size; i++) {
     double wi = gsl_vector_get (cp->w, i);
     double alpha = *(cp->s + i);
     if (alpha != 0.)
       alpha = 1.0 / alpha;
     gsl_vector_set (cp->w, i, alpha * wi);
   }

   // Copies w to GPU
   magma_dsetvector (cp->w->size, cp->w->data, 0, (int) cp->w->stride,
		     cp->dw, 0, (int) cp->w->stride, cp->queue);

   // rhs = 1.0*vt*w + 0*rhs
   magmablas_dgmev ('T', rhs->size, rhs->size, 1.0, cp->dvt, 0,
		    rhs->size, cp->dw, 0, (int) cp->w->stride,
		    0.0, cp->drhs, 0, (int) rhs->stride, cp->queue);
    
   // Get results from GPU
   magma_dgetvector (rhs->size, cp->drhs, 0, (int) rhs->stride,
		     rhs->data, 0, (int) rhs->stride, cp->queue);
}

void opencl_svd_solve (gsl_matrix * A, gsl_vector * rhs)
{
  OpenCLSVDProblem * cp = opencl_svd_factorise (A);

  opencl_solve_factorised_svd (cp, rhs);

  opencl_svd_problem_destroy (cp);
}

#endif

void cuda_boundary_problem_copy_solution_to_patches (GSList * list,
						     BoundaryProblem * bp,
						     gint var)
{
  gint i, j;
  gint istart = 0;
  GSList * plist = list;

#if DEBUG
  FILE * fp = fopen ("sol.tmp","w");
  for ( i = 0; i < bp->rhs->size; i++)
    fprintf (fp, "%i %f\n", i, gsl_vector_get (bp->rhs, i));
  fclose (fp);
#endif

  // Go over the list of patches
  while (plist != NULL) {
    Spline2D * sp = plist->data;
    g_assert (sp != NULL);
    
    if (!sp->periodic) {
      // Loop over the b-splines
      for ( i = 0; i < sp->NU; i++) {
    	for ( j = 0; j < sp->NV; j++) {
    	  coeff_assign (sp, i, j, var,
    			gsl_vector_get (bp->rhs, istart + i + j*sp->NU));
	  
    	}
      }
      // shift index for next patch
      istart += sp->NU*sp->NV;
    }
    else {
      // Loop over the b-splines
      Spline2D * splines = sp;
      gint NUT = splines->NUT;
      gint k = splines->k;
      while (splines) {
	gint fs_index = splines->fs_index;
	gint M = splines->M+k-1;
	for ( i = 0; i < M; i++) {
	  for ( j = 0; j < splines->NV; j++ ) {
	    SplineCoeffs * sc = g_ptr_array_index (splines->coeffs, i + j*M);
	    sc->v[var] = gsl_vector_get (bp->rhs, istart + (fs_index + i)%NUT + j*NUT);
	  }
	}
	splines = splines->next;
      }
      // shift index for next patch
      istart += sp->NUT*sp->NV;
    }
    
    plist = plist->next;
  }
}

/**** SUPERLU ****/

CCSProblem * ccs_problem_new ()
{
  CCSProblem * new = g_malloc ( sizeof (CCSProblem));

  new->matrix = g_array_new (FALSE, FALSE, sizeof (gdouble));
  new->column = g_array_new (FALSE, FALSE, sizeof (gint));
  new->index = g_array_new (FALSE, FALSE, sizeof (gint));
  new->factorised = FALSE;

  return new;
}

void ccs_problem_destroy (CCSProblem * ccs)
{
  g_array_free (ccs->matrix, TRUE);
  g_array_free (ccs->column, TRUE);
  g_array_free (ccs->index, TRUE);
  if ( ccs->factorised ) {
    // Free superlu structures
    SUPERLU_FREE (ccs->perm_r);
    SUPERLU_FREE (ccs->perm_c);
    Destroy_SuperNode_Matrix (&ccs->L);
    Destroy_CompCol_Matrix (&ccs->U);
    StatFree (&ccs->stat);
  }
  g_free (ccs);
  ccs = NULL;
}

static void ccs_problem_lu_factorise (CCSProblem * ccs)
{
  SuperMatrix A;
  gint size = ccs->index->len-1;
  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2);
  double   drop_tol = 0.;
  int lwork = 0, * etree;
  SuperMatrix AC; /* Matrix postmultiplied by Pc */
  int permc_spec;
  superlu_options_t options;

  // Create superlu matrix
  dCreate_CompCol_Matrix (&A, size, size, ccs->matrix->len,
			  (double *) ccs->matrix->data,
			  (int *) ccs->column->data, (int *) ccs->index->data,
			  SLU_NC, SLU_D, SLU_GE);



  // Solver options are defaults
  set_default_options (&options);
  options.ColPerm = NATURAL;
  options.Trans = TRANS;

  // Initialize stats
  StatInit(&ccs->stat);

  ccs->info = 0;

  ccs->perm_r = g_malloc (size*sizeof(gint));
  ccs->perm_c = g_malloc (size*sizeof(gint));
  permc_spec = options.ColPerm;
  get_perm_c(permc_spec, &A, ccs->perm_c);

  etree = intMalloc(A.ncol);

  sp_preorder(&options, &A, ccs->perm_c, etree, &AC);

  dgstrf(&options, &AC, /*drop_tol, */relax, panel_size,
	 etree, NULL, lwork, ccs->perm_c, ccs->perm_r, &ccs->L, &ccs->U, &ccs->stat, &ccs->info);

  SUPERLU_FREE (etree);
  Destroy_CompCol_Permuted(&AC);
  ccs->factorised = TRUE;
}

static void ccs_problem_solve_factorised_LU (CCSProblem * ccs, gsl_vector * rhs)
{
  SuperMatrix B;
  trans_t  trans = NOTRANS;
  gint size = rhs->size;
  gint i;

  g_assert (ccs->factorised == TRUE);

  // Create superlu rhs
  dCreate_Dense_Matrix (&B, size, 1, rhs->data, size, SLU_DN, SLU_D, SLU_GE);

  dgstrs (trans, &ccs->L, &ccs->U, ccs->perm_c, ccs->perm_r, &B, &ccs->stat, &ccs->info);
  
  // Check for success
  g_assert (ccs->info == 0);

  // Here is the solution
  double *sol = (double*) ((DNformat*) B.Store)->nzval;
  
  // Copy the solution to the patch
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, sol [i]);
}

void ccs_problem_lu_solve (CCSProblem * ccs, gsl_vector * rhs)
{
  if (!ccs->factorised)
    ccs_problem_lu_factorise (ccs);
  ccs_problem_solve_factorised_LU (ccs, rhs);
}

/** Plasma **/


PlasmaLUProblem * plasma_lu_problem_new (gint M, gint N)
{
  PlasmaLUProblem * new = g_malloc ( sizeof (PlasmaLUProblem));

  //new->A = g_malloc ((M*N) * sizeof (double));
  //new->rhs = g_malloc (M * sizeof (double));
  //new->ipiv = g_malloc ( min(M,N) * sizeof (gint));

  return new;
}

void plasma_lu_problem_destroy (PlasmaLUProblem * pp)
{
  /* g_free (pp->A); */
  /* g_free (pp->rhs); */
  g_free (pp->L);
  g_free (pp->ipiv);
}

PlasmaLUProblem * plasma_lu_factorise (gsl_matrix * A)
{
  /* Matrix size */
  gint M = A->size1, N = A->size2;
  PlasmaLUProblem * pp = plasma_lu_problem_new (M, N);
  gint info;
  
  pp->A = A->data;

  /* Allocate L and IPIV */
  g_assert (PLASMA_Alloc_Workspace_dgetrf_incpiv (N, N, &pp->L, &pp->ipiv) == 0);

  // Computes LU factorisation
  g_assert (PLASMA_dgetrf_incpiv (M, N, pp->A, M, pp->L, pp->ipiv) == 0);

  return pp;
}

void plasma_solve_factorised_lu (PlasmaLUProblem * pp, gsl_vector * rhs)
{
  // Copy rhs on Plasma framework
  pp->rhs = rhs->data;

  // Solves A x = b using the partial pivoting
  g_assert (PLASMA_dgetrs_incpiv (PlasmaNoTrans, rhs->size, 1, pp->A, rhs->size, pp->L, pp->ipiv, pp->rhs, rhs->size) == 0);
}

void plasma_lu_solve (gsl_matrix * A, gsl_vector * rhs)
{
  PlasmaLUProblem * pp = plasma_lu_factorise (A);

  plasma_solve_factorised_lu (pp, rhs);

  plasma_lu_problem_destroy (pp);
}

void plasma_boundary_problem_copy_solution_to_patches (GSList * list,
						       BoundaryProblem * bp,
						       gint var)
{
  gint i, j;
  gint istart = 0;
  GSList * plist = list;

  // Go over the list of patches
  while (plist != NULL) {
    Spline2D * sp = plist->data;
    g_assert (sp != NULL);
    
    if (!sp->periodic) {
      // Loop over the b-splines
      for ( i = 0; i < sp->NU; i++) {
	for ( j = 0; j < sp->NV; j++) {
	  coeff_assign (sp, i, j, var,
			gsl_vector_get (bp->rhs, istart + i + j*sp->NU));
	
	}
      }
    }
    else {
      // Loop over the b-splines
      for ( i = sp->k-1; i < sp->NU+2*sp->k-2; i++) {
	for ( j = 0; j < sp->NV; j++) {
	  gint indexi = (i - sp->k + 1 ) % sp->NU + j*sp->NU;

	  coeff_assign (sp, i% (sp->M+sp->k-1), j/* % (sp->N+sp->k-1) */, var,
			gsl_vector_get (bp->rhs, istart + indexi));
	}
      }
    }
    
    // shift index for next patch
    istart += sp->NU*sp->NV;
    plist = plist->next;
  }
}


/****** Periodic spline ******/

gint periodic_fs_numbering (Spline2D * fs)
{
  Spline2D * sp = fs;

  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }

  sp = fs;
  while (sp) {
    sp->NUT = NUT;
    sp = sp->next;
  }

  return NUT*fs->NV;
}

void periodic_fs_copy_problem_solution (Spline2D * splines, gsl_vector * lhs, gint var)
{
  gint i, j;
  Spline2D * sp = splines;

  g_assert (var > 2);

  gint NUT = sp->NUT;
  gint k = sp->k;

  while (sp) {
    gint istart = sp->fs_index;
    gint M = sp->M+k-1;
    for ( i = 0; i < sp->NU+k-1; i++) {
      gint ii = (istart + i)%NUT;
      for ( j = 0; j < sp->NV; j++ ) {
	SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, i + j*M);
	sc->v[var] = gsl_vector_get (lhs, ii + j*NUT);
      }
    }
    sp = sp->next;
  }
}

void periodic_fs_add_problem_solution (Spline2D * splines, gsl_vector * lhs, gint var, gint var0)
{
  gint i, j;
  Spline2D * sp = splines;

  g_assert (var > 2);

  gint NUT = sp->NUT;
  gint k = sp->k;

  while (sp) {
    gint istart = sp->fs_index;
    gint M = sp->M+k-1;
    for ( i = 0; i < sp->NU+k-1; i++) {
      gint ii = (istart + i)%NUT;
      for ( j = 0; j < sp->NV; j++ ) {
	SplineCoeffs * sc = g_ptr_array_index (sp->coeffs, i + j*M);
	sc->v[var] = sc->v[var0] + gsl_vector_get (lhs, ii + j*NUT);
      }
    }
    sp = sp->next;
  }
}

void periodic_fs_copy_problem_solution_no_metric (Spline2D * splines,
						  gsl_vector * lhs,
						  gint var)
{
  gint i, j;
  Spline2D * sp = splines;

  gint NUT = 0;
  while (sp) {
    NUT += sp->NXU;
    sp = sp->next;
  }

  sp = splines;
  while (sp) {

    for ( i = 0; i < sp->NXU; i++) {
      for ( j = 0; j < sp->NXV; j++) {
	gint indexi = (sp->fs_index_x + i) + j*NUT;
	coeff_assign (sp, i, j, var, gsl_vector_get (lhs, indexi));
      }
    }

    sp = sp->next;
  }
}

CCSProblem * periodic_fs_build_galerkin_fit_matrix (Spline2D * splines)
{
  gint i, j, m, n, ii, a, b;
  gint NV = splines->NV; // Should be the same for all the splines

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }
  
  gint size = NUT*NV;
  gdouble A0[size][size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }

  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;
  sp = splines;
  while (sp) {// Loop over the spline patches
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = sp->nouter; // Order of outer Gauss-Legendre rule
      ustart = sp->periodic ? 
	gp->istart + sp->fs_index - sp->k + 1 :
	gp->istart;
      vstart = gp->jstart;
      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    gint indexi = i % NUT;
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      gint indexij = indexi + j*NUT;
	      
	      for ( a = 0; a < sp->k; a++) {
		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		gint indexa = (ustart + a)%NUT;
		for ( b = 0; b < sp->k; b++) {
		  A0[indexij][indexa + (vstart+b)*NUT] +=  wmnija*gsl_matrix_get(Bv, b, 0);
		}
	      }

	    }
	  }
	
	}
      }
    }
    sp = sp->next;
  }

  
  //FILE * fp = fopen ("fit.tmp","w");
  CCSProblem * fit = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (fit->index, count);
    for ( j = 0; j < size; j++) {
      //  fprintf(fp,"%i %i %f\n", i, j, A0[i][j]);
      if (A0[i][j] != 0.) {
  	g_array_append_val (fit->matrix, A0[i][j]);
  	g_array_append_val (fit->column, j);
  	count++;
      }
    }
    //fprintf (fp, "\n");
  }
  g_array_append_val (fit->index, count);
  //fclose (fp);

  return fit;
}

CCSProblem * periodic_fs_build_galerkin_fit_matrix_bc (Spline2D * splines)
{
  gint i, j, k, m, n, ii, a, b;
  gint NV = splines->NV; // Should be the same for all the splines

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }
  
  gint size = NUT*NV;
  gdouble A0[size][size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }

  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;
  sp = splines;
  while (sp) {// Loop over the spline patches
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = sp->nouter; // Order of outer Gauss-Legendre rule
      ustart = sp->periodic ? 
	gp->istart + sp->fs_index - sp->k + 1 :
	gp->istart;
      vstart = gp->jstart;
      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    gint indexi = i % NUT;
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      gint indexij = indexi + j*NUT;
	      
	      for ( a = 0; a < sp->k; a++) {
		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		gint indexa = (ustart + a)%NUT;
		for ( b = 0; b < sp->k; b++) {
		  A0[indexij][indexa + (vstart+b)*NUT] +=  wmnija*gsl_matrix_get(Bv, b, 0);
		}
	      }

	    }
	  }
	
	}
      }
    }
    sp = sp->next;
  }

#if 1
  sp = splines;
  while (sp) {
    // Dirichlet on edges
    gsl_vector * Bu2 = gsl_vector_alloc (sp->k);
    gsl_matrix * Bv2 = gsl_matrix_alloc (sp->k, sp->k);
    GrevillePoints * gr = sp->gr;
    size_t ustart, uend, vstart, vend;
    for ( i = 0; i < sp->NU; i++) {
      gdouble u = g_array_index (gr->ui, gdouble, i);
      gdouble v = 1.;

      gsl_bspline_eval_nonzero (u, Bu2, &ustart, &uend, sp->w_u);
      gsl_bspline_deriv_eval_nonzero (v, sp->k-1, Bv2, &vstart, &vend,
				      sp->w_v, sp->wd_v);

      for ( k = 0; k < sp->k; k++) {
	gint index1 = sp->fs_index + i + NUT*(NV-1-k);      
	for ( j = 0; j < size; j++)
	  A0[j][index1] = 0.;

	for ( m = 0; m < sp->k; m++) {
	  gdouble cu = gsl_vector_get (Bu2, m);
	  for ( n = 0; n < sp->k; n++) {
	    gdouble cv = gsl_matrix_get (Bv2, n, k);
	    A0[(sp->fs_index + ustart + m)%NUT + (vstart+n)*NUT] [index1] = cu*cv;
	  }
	}

      }
      /* index1 = sp->fs_index + i + NUT*(NV-2); */
      
      /* for ( j = 0; j < size; j++) */
      /* 	A0[j][index1] = 0.; */

      /* for ( m = 0; m < sp->k; m++) { */
      /* 	gdouble cu = gsl_matrix_get (Bu2, m, 0); */
      /* 	for ( n = 0; n < sp->k; n++) { */
      /* 	  gdouble cv = gsl_matrix_get (Bv2, n, 1); */
      /* 	  A0[(sp->fs_index + ustart + m)%NUT + (vstart+n)*NUT] [index1] = cu*cv; */
      /* 	} */
      /* } */


    }
    sp = sp->next;
    gsl_vector_free (Bu2);
    gsl_matrix_free (Bv2);
  }
#endif

  //FILE * fp = fopen ("fit.tmp","w");
  CCSProblem * fit = ccs_problem_new ();
  gint count = 0;
  for ( i = 0; i < size; i++) {
    g_array_append_val (fit->index, count);
    for ( j = 0; j < size; j++) {
      //  fprintf(fp,"%i %i %f\n", i, j, A0[i][j]);
      if (A0[i][j] != 0.) {
  	g_array_append_val (fit->matrix, A0[i][j]);
  	g_array_append_val (fit->column, j);
  	count++;
      }
    }
    //fprintf (fp, "\n");
  }
  g_array_append_val (fit->index, count);
  //fclose (fp);

  return fit;
}


CCSProblem * periodic_fs_build_galerkin_fit_matrix_no_metric (Spline2D * splines)
{
  gint i, j, m, n, ii, a, b;
  gint NV = splines->NXV; // Should be the same for all the splines

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index_x = NUT;
    NUT += sp->NXU;
    sp = sp->next;
  }
  gint size = NUT*NV;
  gdouble A0[size][size];


  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }
  
  gsl_matrix * Bu, * Bv;
  size_t ustart, uend, vstart, vend;
  sp = splines;
  while (sp) {
    gint fs_index = sp->fs_index_x;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = sp->nouter; // Order of outer Gauss-Legendre rule
      ustart = gp->istart_x + fs_index;
      vstart = gp->jstart;
      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bux, m);
	  
	for ( n = 0; n < ng; n++) {	    
	  gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	  Bv = g_ptr_array_index (gp->Bv, n);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      gint indexi = i + j*NUT;
	      
	      for ( a = 0; a < sp->k; a++) {
		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		for ( b = 0; b < sp->k; b++) {
		  A0[indexi][(ustart+a) + (vstart+b)*NUT] +=  wmnija*gsl_matrix_get (Bv, b, 0);
		}
	      }
	    }
	  }
	
	}
      }
    
    }
    
    // Continuity condition
    Spline2D * sp2 = sp->next == NULL ? splines: sp->next;

    GrevillePoints * gr = sp->gr;
    gsl_vector * Bu2 = gsl_vector_alloc (sp->k);
    gsl_vector * Bv2 = gsl_vector_alloc (sp->k);
    for ( j = 0; j < sp->NXV; j++) {
      gint indexi = (sp->fs_index_x + sp->NXU - 1) + j*NUT;

      for ( i = 0; i < size; i++)
    	A0[i][indexi] = 0.;

      gdouble u = 1.;
      gdouble v = g_array_index (gr->vj, gdouble, j);
      gsl_bspline_eval_nonzero (MIN(1.-1e-12,u), Bu2, &ustart, &uend, sp->wx_u);
      gsl_bspline_eval_nonzero (MIN(1.-1e-12,v), Bv2, &vstart, &vend, sp->w_v);

      for ( m = 0; m < sp->k; m++) {
    	gdouble cu = gsl_vector_get (Bu2, m);
    	for ( n = 0; n < sp->k; n++) {
    	  A0[sp->fs_index_x + ustart + m + (vstart+n)*NUT][indexi] =
    	    cu*gsl_vector_get (Bv2, n);
    	}
      }

      u = 0.;
      gsl_bspline_eval_nonzero (MIN(1.-1e-12,u), Bu2, &ustart, &uend, sp2->wx_u);
      gsl_bspline_eval_nonzero (MIN(1.-1e-12,v), Bv2, &vstart, &vend, sp2->w_v);
      
      for ( m = 0; m < sp2->k; m++) {
    	gdouble cu = gsl_vector_get (Bu2, m);
    	for ( n = 0; n < sp2->k; n++) {
    	  A0[sp2->fs_index_x + ustart + m + (vstart+n)*NUT][indexi] -=
    	    cu*gsl_vector_get (Bv2, n);
    	}
      }
    }
    gsl_vector_free (Bu2);
    gsl_vector_free (Bv2);

    sp = sp->next;
  }

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
  }
  g_array_append_val (fit->index, count);

  return fit;
}

/**
 * This imposes noflux using a Neumann condition on the edge of the freesurface
 **/
CCSProblem * periodic_fs_build_galerkin_fit_noflux_matrix/* _neumann */ (Spline2D * splines)
{
  gint i, j, m, n, ii, a, b;
  gint NV = splines->NV; // Should be the same for all the splines

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }
  
  gint size = NUT*NV;
  gdouble A0[size][size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }

  gsl_matrix * Bu, * Bv;
  sp = splines;
  while (sp) {
    // Loop over the spline patches
    gint fs_index = sp->fs_index;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = sp->nouter; // Order of outer Gauss-Legendre rule
      gint ustart = gp->istart + fs_index - sp->k + 1;
      gint vstart = gp->jstart;
      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      gint indexi = i % NUT + j*NUT;
	      
	      for ( a = 0; a < sp->k; a++) {
		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		gint indexa = (ustart + a)%NUT;
		for ( b = 0; b < sp->k; b++) {
		  A0[indexi][indexa + (vstart+b)*NUT] +=  wmnija*gsl_matrix_get(Bv, b, 0);
		}
	      }

	    }
	  }
	
	}
      }
    }
    sp = sp->next;
  }

  // On the solid boundary dn phi = - dn phi0
  // For now, dn phi = 0 i.e. dv phi = 0
  sp = splines;
  while (sp) {
    size_t ustart, uend, vstart, vend;
    gsl_vector * Bu2 = gsl_vector_alloc (sp->k);
    Bv = gsl_matrix_alloc (sp->k, 2);
    GrevillePoints * gr = sp->gr;
    for ( i = 0; i < sp->NU; i++) {
      gint index1 = sp->fs_index + i;
      gdouble u = g_array_index (gr->ui, gdouble, i);
      gdouble v = 0.;
      
      for ( j = 0; j < size; j++)
	A0[j][index1] = 0.;

      gsl_bspline_eval_nonzero (MIN(1.-1e-12,u), Bu2, &ustart, &uend, sp->w_u);
      gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

      for ( m = 0; m < sp->k; m++) {
	gdouble cu = gsl_vector_get (Bu2, m);
	for ( n = 0; n < sp->k; n++) {
	  gdouble cv = gsl_matrix_get (Bv, n, 1);
	  A0[(sp->fs_index+ustart+m - sp->k +1)%NUT + (vstart+n)*NUT][index1] = cu*cv;
	}
      }
    }
    gsl_vector_free (Bu2);
    gsl_matrix_free (Bv);
    sp = sp->next;
  }

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
  }
  g_array_append_val (fit->index, count);

  return fit;
}


CCSProblem * periodic_fs_build_galerkin_fit_noflux_matrix_dirichlet (Spline2D * splines)
{
  gint i, j, m, n, ii, a, b;
  gint NV = splines->NV; // Should be the same for all the splines

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }
  
  gint size = NUT*NV;
  gdouble A0[size][size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    for ( j = 0; j < size; j++) {
      A0[i][j] = 0.;
    }
  }

  gsl_matrix * Bu, * Bv;
  sp = splines;
  while (sp) {
    // Loop over the spline patches
    gint fs_index = sp->fs_index;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = sp->nouter; // Order of outer Gauss-Legendre rule
      gint ustart = gp->istart + fs_index - sp->k + 1;
      gint vstart = gp->jstart;
      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      gint indexi = i % NUT + j*NUT;
	      
	      for ( a = 0; a < sp->k; a++) {
		gdouble wmnija = wmnij*gsl_matrix_get (Bu, a, 0);
		for ( b = 0; b < sp->k; b++) {
		  A0[indexi][(( ustart + a)%NUT) + (vstart+b)*NUT] +=  wmnija*gsl_matrix_get(Bv, b, 0);
		}
	      }

	    }
	  }
	
	}
      }
    }
    sp = sp->next;
  }

  // On the solid boundary dn phi = - dn phi0
  // For now, dn phi = 0 i.e. dv phi = 0
  sp = splines;
  while (sp) {
    size_t ustart, uend, vstart, vend;
    gsl_matrix * Bu2 = gsl_matrix_alloc (sp->k, 2);
    gsl_matrix * Bv2 = gsl_matrix_alloc (sp->k, 2);
    GrevillePoints * gr = sp->gr;
    for ( i = 0; i < sp->NU; i++) {
      gint index1 = sp->fs_index + i;
      gdouble u = g_array_index (gr->ui, gdouble, i);
      gdouble v = 0.;
      
      for ( j = 0; j < size; j++)
	A0[j][index1] = 0.;

      gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,u), 1, Bu2, &ustart, &uend, sp->w_u, sp->wd_u);
      gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 1, Bv2, &vstart, &vend, sp->w_v, sp->wd_v);

      for ( m = 0; m < sp->k; m++) {
	gdouble cu = gsl_matrix_get (Bu2, m, 0);
	for ( n = 0; n < sp->k; n++) {
	  gdouble cv = gsl_matrix_get (Bv2, n, 0);
	  A0[(sp->fs_index + ustart + m - sp->k +1)%NUT + (vstart+n)*NUT][index1] = cu*cv;
	}
      }
    }
    gsl_matrix_free (Bu2);
    gsl_matrix_free (Bv2);
    sp = sp->next;
  }

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
  }
  g_array_append_val (fit->index, count);

  return fit;
}

gsl_vector * periodic_fs_build_galerkin_rhs_gauss (Spline2D * splines,
						   GaussFunc func,
						   gpointer data,
						   Spline2DFunc bc_func,
						   gpointer bc_data,
						   gsl_vector * rhs)
{
  gint i, j, m, n, ii;

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }

  gint NV = splines->NV;
  gint size = NUT*NV;
  gdouble RHS[size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++)
    RHS[i] = 0.;

  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;
  sp = splines;
  while (sp) {
    gint fs_index = sp->fs_index;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
      ustart = gp->istart + fs_index - sp->k + 1 ;
      vstart = gp->jstart;

      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	  gdouble fmn = func (spp, m, n, data);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    gint indexi = i % NUT;
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      RHS[indexi + j*NUT] += wmni*gsl_matrix_get (Bv, j-vstart, 0)*fmn;
	    }
	  }

	}
      }
    }
    sp = sp->next;
  }

  

  // Copy RHS
  /* FILE * fp = fopen ("fit_rhs.tmp","w"); */
  rhs = splines->rhs_vector (splines, rhs);
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs, i, RHS[i]);
    /* fprintf (fp, "%i %f\n", i, RHS[i]); */
  }
  /* fclose (fp); */

  return rhs;
}

gsl_vector * periodic_fs_build_galerkin_rhs_gauss_bc (Spline2D * splines,
						      GaussFunc func,
						      gpointer data,
						      Spline2DFunc bc_func,
						      gpointer bc_data,
						      gsl_vector * rhs)
{
  gint i, j, k, m, n, ii;

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }

  gint NV = splines->NV;
  gint size = NUT*NV;
  gdouble RHS[size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++)
    RHS[i] = 0.;

  gsl_matrix * Bu, * Bv;
  gint ustart, vstart;
  sp = splines;
  while (sp) {
    gint fs_index = sp->fs_index;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
      ustart = gp->istart + fs_index - sp->k + 1 ;
      vstart = gp->jstart;

      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wJij, gdouble, m+n*ng);
	  gdouble fmn = func (spp, m, n, data);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    gint indexi = i % NUT;
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      RHS[indexi + j*NUT] += wmni*gsl_matrix_get (Bv, j-vstart, 0)*fmn;
	    }
	  }

	}
      }
    }
    sp = sp->next;
  }

#if 1
  // Edges
  sp = splines;
  while (sp) {
    GrevillePoints * gr = sp->gr;
    for ( i = 0; i < sp->NU; i++) {
      for ( k = 0; k < sp->k; k++) {
	gint index1 = sp->fs_index + i + NUT*(NV-1-k);
	gdouble u = g_array_index (gr->ui, gdouble, i);
	gdouble v = 1.;
	RHS[index1] = 0.;
	// RHS[sp->fs_index + i + NUT*(NV-2)] = 0.;
      }
    }
    sp = sp->next;
  }
#endif    

  // Copy RHS
  /* FILE * fp = fopen ("fit_rhs.tmp","w"); */
  rhs = splines->rhs_vector (splines, rhs);
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs, i, RHS[i]);
    /* fprintf (fp, "%i %f\n", i, RHS[i]); */
  }
  /* fclose (fp); */

  return rhs;
}

gsl_vector * periodic_fs_build_galerkin_rhs_gauss_no_metric (Spline2D * splines,
							     GaussFunc func,
							     gpointer data,
							     Spline2DFunc bc_func,
							     gpointer bc_data)
{
  gint i, j, m, n, ii, a, b;
  gint NV = splines->NXV; // Should be the same for all the splines

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index_x = NUT;
    NUT += sp->NXU;
    sp = sp->next;
  }
  gint size = NUT*NV;
  gdouble RHS[size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++) {
    RHS[i] = 0.;
  }
  
  gsl_matrix * Bu, * Bv;
  size_t ustart, uend, vstart, vend;
  sp = splines;
  while (sp) {
    gint fs_index = sp->fs_index_x;
    gint NU = sp->NXU;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->panels->len; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = gp->ui->len;// Order of outer Gauss-Legendre rule
      ustart = gp->istart_x + fs_index;
      vstart = gp->jstart;
      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bux, m);
	  
	for ( n = 0; n < ng; n++) {
	  gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);
	  Bv = g_ptr_array_index (gp->Bv, n);
	
	  gdouble fmn = func (spp, m, n, data);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      gint indexi = i + j*NUT;
	 
	      RHS[indexi] += wmnij*fmn;
	    }
	  }
	
	}
      }
    }

    // Continuity condition
    for ( j = 0; j < sp->NXV; j++) {
      gint indexi = (sp->fs_index_x + sp->NXU - 1) + j*NUT;
      RHS[indexi] = 0.;
    }

    sp = sp->next;
  }

  // Copy lhs and rhs to gsl structures
  gsl_vector * rhs = gsl_vector_alloc (size);
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, RHS[i]);

  return rhs;
}


gsl_vector * periodic_fs_build_galerkin_noflux_rhs_gauss (Spline2D * splines, 
							  GaussFunc func, 
							  gpointer data, 
							  Spline2DFunc bc_func, 
							  gpointer bc_data,
							  gsl_vector * rhs)
{
  gint i, j, m, n, ii;

  Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }

  gint NV = splines->NV;
  gint size = NUT*NV;
  gdouble RHS[size];

  // Initializes the coefficients to 0
  for ( i = 0; i < size; i++)
    RHS[i] = 0.;

  gsl_matrix * Bu, * Bv;
  sp = splines;
  while (sp) {
    gint fs_index = sp->fs_index;
    // Loop over the panels of the patch
    for ( ii = 0; ii < sp->M*sp->N; ii++) {
      SPPanel * spp = g_ptr_array_index (sp->panels, ii);
      g_assert (spp != NULL);
     	
      // Gauss outer-integration
      GaussPoints * gp = spp->outer;
      gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
      gint ustart = gp->istart + fs_index - sp->k + 1;
      gint vstart = gp->jstart;

      for ( m = 0; m < ng; m++) {
	Bu = g_ptr_array_index (gp->Bu, m);
	  
	for ( n = 0; n < ng; n++) {
	  Bv = g_ptr_array_index (gp->Bv, n);
	  gdouble wmn = g_array_index (gp->wij, gdouble, m+n*ng);

	  gdouble fmn = func (spp, m, n, data);

	  // Loop over the splines whose support is included in the panel
	  for ( i = ustart; i < ustart + sp->k; i++) {
	    gdouble wmni = wmn*gsl_matrix_get (Bu, i-ustart, 0);
	    gint indexi = i % NUT;
	    for ( j = vstart; j < vstart + sp->k; j++) {
	      gdouble wmnij = wmni*gsl_matrix_get (Bv, j-vstart, 0);
	      

	      RHS[indexi + j*NUT] += wmnij*fmn;
	    }
	  }
	
	}
      }
    }
    sp = sp->next;
  }

  sp = splines;
  while (sp) {
    GrevillePoints * gr = sp->gr;
    for ( i = 0; i < sp->NU; i++) {
      gint index1 = i;
      gdouble u = g_array_index (gr->ui, gdouble, sp->fs_index + i);;
      gdouble v = 0.;

      RHS[i] = bc_func (sp, u, v, bc_data);
    }
    sp = sp->next;
  }

  // Copy RHS
  rhs = splines->rhs_vector (splines, rhs);
  for ( i = 0; i < size; i++)
    gsl_vector_set (rhs, i, RHS[i]);

  return rhs;
}

gint periodic_fs_size (Spline2D * sp)
{
  gint NV = sp->NV;
  gint NUT = 0;

  while (sp) {
    sp->fs_index = NUT;
    NUT += sp->NU;
    sp = sp->next;
  }
  
  return NUT*NV;
}

Spline2D * periodic_fs_new (gint M, gint N, gint k, gint ninner, gint nouter)
{
  Spline2D * s = g_malloc (sizeof(Spline2D));
  
  // Order of the splines
  s->k = k;
  s->N = N;
  s->M = M;
  s->NUT = 0;
  s->ninner = ninner;
  s->nouter = nouter;

  // Neighbours (periodic)
  s->next = NULL;

  // Periodicity
  s->periodic = TRUE;
  s->noflux = TRUE;

  // Number of knots so that it fits with the notation of Maliar et al,
  s->w_u = gsl_bspline_alloc (k, M+2*k-1);
  s->wd_u = gsl_bspline_deriv_alloc (k);

  s->w_v = gsl_bspline_alloc (k, N+1);
  s->wd_v = gsl_bspline_deriv_alloc (k);

  s->wx_u = gsl_bspline_alloc (k, M+1);
  s->wxd_u = gsl_bspline_deriv_alloc (k);

  // break = M+1/N+1
  // NB: GSL collocate k knots at the start and end of the spline
  gsl_bspline_knots_uniform (0. - (k-1.)/M, 1. + (k-1.)/M, s->w_u);
  gsl_bspline_knots_uniform (0. , 1., s->w_v);

  gsl_bspline_knots_uniform (0., 1., s->wx_u);
  
  // Store number of spline coefficients in each direction
  s->NU = s->M;
  s->NV = gsl_bspline_ncoeffs (s->w_v);

  s->NXU = gsl_bspline_ncoeffs (s->wx_u);
  s->NXV = s->NV;

  // Allocates space for the coeffs of the splines.
  s->coeffs = g_ptr_array_new ();
  gint i, j;
  for ( i = 0; i < M+k-1; i++) { //
    for ( j = 0; j < N+ k-1; j++) {
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

  //  s->build_fit_matrix = periodic_fs_build_galerkin_fit_matrix;
  s->build_fit_matrix = periodic_fs_build_galerkin_fit_matrix_bc;
  //s->build_fit_rhs = periodic_fs_build_galerkin_rhs_gauss;
  s->build_fit_rhs = periodic_fs_build_galerkin_rhs_gauss_bc;
  s->build_fit_noflux_matrix = periodic_fs_build_galerkin_fit_noflux_matrix;
  s->build_fit_noflux_rhs = periodic_fs_build_galerkin_noflux_rhs_gauss;
  s->copy_fit_solution = periodic_fs_copy_problem_solution;
  s->add_fit_solution = periodic_fs_add_problem_solution;
  s->reinit_panels = periodic_fs_reinit_panels_physical_quantities;
  s->size = periodic_fs_size;
  s->rhs_vector = spline2d_rhs_vector;

  s->fit = NULL;
  s->fit_noflux = NULL;
  s->fully_submerged = FALSE;
  s->hull_patch = NULL;
  s->rhs = NULL;

  coeff_set_var_to_zero (s, 0);
  coeff_set_var_to_zero (s, 1);
  coeff_set_var_to_zero (s, 2);

  return s;
}

static GrevillePoints * periodic_fs_store_greville_points (Spline2D * sp)
{
  gint NU = sp->NU;
  gint NV = sp->NV;
  gint i, j;
  GrevillePoints * gr = greville_points_new (NU, NV);
 
  size_t ustart, uend, vstart, vend;
  for ( i = 0; i < NU; i++) {
    SPPanel * spp = g_ptr_array_index (sp->panels, i);

    gdouble u = g_array_index (gr->ui, gdouble, i) = spp->ue;
    g_assert ( u >= 0. && u <= 1.);
    gsl_matrix * Bu = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,u), 1, Bu, &ustart, &uend, sp->w_u, sp->wd_u);
    g_ptr_array_add (gr->Bu, Bu);
    g_array_append_val (gr->ustart, ustart);
  }
    
  for ( i = 0; i < NV; i++) {
    gdouble v = g_array_index (gr->vj, gdouble, i) = gsl_bspline_greville_abscissa (i, sp->w_v);
    gsl_matrix * Bv = gsl_matrix_alloc (sp->k, 2);
    gsl_bspline_deriv_eval_nonzero (MIN(1.-1e-12,v), 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
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

      g_array_append_val (gr->Ni, N);
      g_array_append_val (gr->Pi, P);
      g_array_append_val (gr->Ji, J);
    }
  }

  return gr;
}

void periodic_fs_init_panels (Spline2D * sp)
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

      // Compute and store Gauss-Legendre points and weights inner and outer integration
      spp->outer = sppanel_store_gauss_legendre_points (spp, sp->nouter);

    }
  }

  // Greville Points for collocation approach
  sp->gr = periodic_fs_store_greville_points (sp);
}

void periodic_fs_reinit_panels_physical_quantities (Spline2D * splines)
{
  gint i, j;
  Spline2D * sp = splines;

  // Make sure that numbering is OK
  periodic_fs_numbering (splines);

  // Loop over all the panels
  while (sp) {
    for (i = 0; i < sp->M; i++) {
      for (j = 0; j < sp->N; j++) {
	SPPanel * spp = g_ptr_array_index (sp->panels, i + j*sp->M);
	g_assert (spp != NULL);
	
	// Panel parametric centroid
	spp->pe = spline2d_eval_point (sp, spp->ue, spp->ve);
	
	// Compute and store physical quantities at Gauss-Legendre points
	sppanel_store_gauss_legendre_data (spp->outer, spp);
      }
    }
    spline2d_store_greville_points_data (sp);

    sp = sp->next;
  }

  g_warning ("Normal orientation needs checking manually before starting any serious run\n");
}

/**
 * Applies an iteration of elliptic grid smoother as described in the Annexes of (Huang, 1995)
 * and in (Kho et al., 2011) page 855.
 **/
void periodic_fs_elliptic_smoothing (Spline2D * splines)
{
  gint i, j, a, b, m, n;
  g_assert (splines != NULL);
  
   Spline2D * sp = splines;
  // Find the total number of coeffs in direction u
  // and number the spline patches
  gint NUT = 0;
  while (sp) {
    sp->fs_index_x = NUT;
    NUT += sp->NXU;
    sp = sp->next;
  }

  gint NV = splines->NXV;
  gint size = NUT*NV;
  gdouble A0[size][size];
  gdouble RHSX[size];
  gdouble RHSY[size];

  for ( i = 0; i < size; i++) {
    RHSX[i] = 0.;
    RHSY[i] = 0.;
    for ( j = 0; j < size; j++)
      A0[i][j] = 0.;
  }
  gint k = splines->k;
  gsl_matrix * Bu = gsl_matrix_alloc (k, 3);
  gsl_matrix * Bv = gsl_matrix_alloc (k, 3);
  size_t ustart, uend, vstart, vend;

  sp = splines;
  while (sp) {
    gint fs_index = sp->fs_index_x;
    for ( i = 0; i < sp->NXU; i++) {
      gdouble u = gsl_bspline_greville_abscissa (i, sp->wx_u);
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 2, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
      ustart += fs_index;

      for ( j = 1; j < NV-1; j++) { // Inner domain
	gdouble v = gsl_bspline_greville_abscissa (j, sp->w_v);
	gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,v)), 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v);
      
	gdouble dxdu = 0., dydu = 0., dxdv = 0., dydv = 0.;

	for ( m = 0; m < sp->k; m++) {
	  gdouble cu = gsl_matrix_get (Bu, m, 0);
	  gdouble cdu = gsl_matrix_get (Bu, m, 1);
	  for ( n = 0; n < sp->k; n++) {
	    gdouble cvdu = cdu*gsl_matrix_get (Bv, n, 0);
	    gdouble cudv = cu*gsl_matrix_get (Bv, n, 1);
	    gdouble v0 = coeff (sp, ustart+m, vstart+n, 0);
	    gdouble v1 = coeff (sp, ustart+m, vstart+n, 1);

	    dxdu += v0*cvdu;
	    dxdv += v0*cudv;
	    dydu += v1*cvdu;
	    dydv += v1*cudv;
	  }
	}


	gdouble J = dxdu*dydv - dydu*dxdv;
	gdouble alpha = dxdv*dxdv + dydv*dydv;
	gdouble beta = dxdu*dxdv + dydu*dydv; // Typo in (Kho et al., 2011) beta and gamma are swapped
	gdouble gamma = dxdu*dxdu + dydu*dydu;
	
	gdouble lpp = /* 70 */10;
	gdouble Q = lpp == 0. ? 0 : 100./(lpp*lpp)*exp(-0.5*(j-1)); // Grid control functions (Kho et al., 2011) page 856
	gdouble P = 0.;

	gint indexi = fs_index + i + j*NUT; // no %NU as it is the index of the point and not that of the spline
      
	for ( a = 0; a < k; a++) {
	  gint atmp = (ustart + a);
	  for ( b = 0; b < k; b++) {
	    A0[atmp + (vstart+b)*NUT][indexi] +=
	      alpha*gsl_matrix_get (Bu, a, 2)*gsl_matrix_get (Bv, b, 0)
	      -2.*beta*gsl_matrix_get (Bu, a, 1)*gsl_matrix_get (Bv, b, 1)
	      +gamma*gsl_matrix_get (Bu, a, 0)*gsl_matrix_get (Bv, b, 2)
	      + J*J*(P*gsl_matrix_get (Bu, a, 1)*gsl_matrix_get (Bv, b, 0)
		     + Q*gsl_matrix_get (Bu, a, 0)*gsl_matrix_get (Bv, b, 1));
	  }
	}
	// NB: RHS is zero
      }
    }
    
    // Keep boundary points at v = 0 and v = 1 where they are
    for ( i = 0; i < sp->NXU; i++) {
      gdouble u = gsl_bspline_greville_abscissa (i, sp->wx_u);

      gint index1 = sp->fs_index_x + i; // Inner boundary
      gdouble v = 0.;
      
      for ( j = 0; j < size; j++)
    	A0[j][index1] = 0.;

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 2, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
      gsl_bspline_deriv_eval_nonzero (v, 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

      for ( a = 0; a < sp->k; a++) {
    	gdouble cu = gsl_matrix_get (Bu, a, 0);
    	for ( b = 0; b < sp->k; b++) {
    	  gdouble cv = gsl_matrix_get (Bv, b, 0);
    	  A0[(sp->fs_index_x + ustart + a) + (vstart+b)*NUT][index1] = cu*cv;
    	}
      }

      Point p = spline2d_eval_point (sp, u, v);
      RHSX[index1] = p.x;
      RHSY[index1] = p.y;

 
      index1 = sp->fs_index_x + i + (sp->NXV-1)*NUT; // Outer boundary
      v = 1.;
      
      for ( j = 0; j < size; j++)
    	A0[j][index1] = 0.;

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 2, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
      gsl_bspline_deriv_eval_nonzero (v, 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

      for ( a = 0; a < sp->k; a++) {
    	gdouble cu = gsl_matrix_get (Bu, a, 0);
    	for ( b = 0; b < sp->k; b++) {
    	  gdouble cv = gsl_matrix_get (Bv, b, 0);
    	  A0[(sp->fs_index_x + ustart + a) + (vstart+b)*NUT][index1] = cu*cv;
    	}
      }
      p = spline2d_eval_point (sp, u, v);
      RHSX[index1] = p.x;
      RHSY[index1] = p.y;

      
    }

    // Continuity condition + normal derivative continuity between consecutive patches
    Spline2D * sp2 = sp->next == NULL ? splines: sp->next;    
    GrevillePoints * gr = sp->gr;
    for ( j = 0; j < sp->NXV; j++) {
      gint index1 = (sp->fs_index_x + sp->NXU - 1) + j*NUT;
      gint index2 = sp->fs_index_x + j*NUT;

      for ( i = 0; i < size; i++)
    	A0[i][index1] = A0[i][index2] = 0.;

      gdouble u = 1.;
      gdouble v = g_array_index (gr->vj, gdouble, j);

      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u);
      gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v);

      for ( m = 0; m < sp->k; m++) {
    	gdouble cu = gsl_matrix_get (Bu, m, 0);
	gdouble cdu = gsl_matrix_get (Bu, m, 1);
    	for ( n = 0; n < sp->k; n++) {
    	  A0[sp->fs_index_x + ustart + m + (vstart+n)*NUT][index1] =
    	    cu*gsl_matrix_get (Bv, n, 0);
	  A0[sp->fs_index_x + ustart + m + (vstart+n)*NUT][index2] =
    	    cdu*gsl_matrix_get (Bv, n, 0);
    	}
      }

      u = 0.;
      gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp2->wx_u, sp2->wxd_u);
      gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &vstart, &vend, sp2->w_v, sp2->wd_v);
      
      for ( m = 0; m < sp2->k; m++) {
    	gdouble cu = gsl_matrix_get (Bu, m, 0);
	gdouble cdu = gsl_matrix_get (Bu, m, 1);
    	for ( n = 0; n < sp2->k; n++) {
    	  A0[sp2->fs_index_x + ustart + m + (vstart+n)*NUT][index1] -=
    	    cu*gsl_matrix_get (Bv, n, 0);
	  A0[sp2->fs_index_x + ustart + m + (vstart+n)*NUT][index2] -=
    	    cdu*gsl_matrix_get (Bv, n, 0);
    	}
      }

      RHSX[index1] = RHSY[index1] = 0.;
      RHSX[index2] = RHSY[index2] = 0.;
    }
    sp = sp->next;
  }

  gsl_matrix_free (Bu);
  gsl_matrix_free (Bv);

  /* Storage in Compressed Column Storage (CCS) format for superlu */
  CCSProblem * css = ccs_problem_new ();
  gsl_vector * rhs_x = gsl_vector_alloc (size);
  gsl_vector * rhs_y = gsl_vector_alloc (size);
  gint count = 0;
  for ( i = 0; i < size; i++) {
    gsl_vector_set (rhs_x, i, RHSX[i]);
    gsl_vector_set (rhs_y, i, RHSY[i]);
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

  // Solve problem
  ccs_problem_lu_solve (css, rhs_x);
  ccs_problem_lu_solve (css, rhs_y);

  // Store solution
  sp = splines;
  while (sp) {

    for ( i = 0; i < sp->NXU; i++) {
      for ( j = 1; j < sp->NXV-1; j++) { // Don't copy the coefficients on the borders as we don't want them to change at all
  	gint indexi = (sp->fs_index_x + i) + j*NUT;
  	coeff_assign (sp, i, j, 0, gsl_vector_get (rhs_x, indexi));
  	coeff_assign (sp, i, j, 1, gsl_vector_get (rhs_y, indexi));
      }
    }

    sp = sp->next;
  }

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  ccs_problem_destroy (css);
}

/* void spline2d_elliptic_smoothing (Spline2D * sp) */
/* { */
/*   gint i, j, a, b, m, n; */
/*   g_assert (spl != NULL); */

/*   gint NU = sp->NU; */
/*   gint NV = sp->NV; */
/*   gint size = NU*NV; */
/*   gdouble A0[size][size]; */
/*   gdouble RHSX[size]; */
/*   gdouble RHSY[size]; */

/*   for ( i = 0; i < size; i++) { */
/*     RHSX[i] = 0.; */
/*     RHSY[i] = 0.; */
/*     for ( j = 0; j < size; j++) */
/*       A0[i][j] = 0.; */
/*   } */
/*   gint k = splines->k; */
/*   gsl_matrix * Bu = gsl_matrix_alloc (k, 3); */
/*   gsl_matrix * Bv = gsl_matrix_alloc (k, 3); */
/*   size_t ustart, uend, vstart, vend; */

  
/*     for ( i = 0; i < sp->NXU; i++) { */
/*       gdouble u = gsl_bspline_greville_abscissa (i, sp->wx_u); */
/*       gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 2, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u); */
/*       ustart += fs_index; */

/*       for ( j = 1; j < NV-1; j++) { // Inner domain */
/* 	gdouble v = gsl_bspline_greville_abscissa (j, sp->w_v); */
/* 	gsl_bspline_deriv_eval_nonzero (v, 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v); */
      
/* 	gdouble dxdu = 0., dydu = 0., dxdv = 0., dydv = 0.; */

/* 	for ( m = 0; m < sp->k; m++) { */
/* 	  gdouble cu = gsl_matrix_get (Bu, m, 0); */
/* 	  gdouble cdu = gsl_matrix_get (Bu, m, 1); */
/* 	  for ( n = 0; n < sp->k; n++) { */
/* 	    gdouble cvdu = cdu*gsl_matrix_get (Bv, n, 0); */
/* 	    gdouble cudv = cu*gsl_matrix_get (Bv, n, 1); */
/* 	    gdouble v0 = coeff (sp, ustart+m, vstart+n, 0); */
/* 	    gdouble v1 = coeff (sp, ustart+m, vstart+n, 1); */

/* 	    dxdu += v0*cvdu; */
/* 	    dxdv += v0*cudv; */
/* 	    dydu += v1*cvdu; */
/* 	    dydv += v1*cudv; */
/* 	  } */
/* 	} */


/* 	gdouble J = dxdu*dydv - dydu*dxdv; */
/* 	gdouble alpha = dxdv*dxdv + dydv*dydv; */
/* 	gdouble beta = dxdu*dxdv + dydu*dydv; // Typo in (Kho et al., 2011) beta and gamma are swapped */
/* 	gdouble gamma = dxdu*dxdu + dydu*dydu; */
	
/* 	gdouble lpp = /\* 70 *\/10; */
/* 	gdouble Q = lpp == 0. ? 0 : 100./(lpp*lpp)*exp(-0.5*(j-1)); // Grid control functions (Kho et al., 2011) page 856 */
/* 	gdouble P = 0.; */

/* 	gint indexi = fs_index + i + j*NUT; // no %NU as it is the index of the point and not that of the spline */
      
/* 	for ( a = 0; a < k; a++) { */
/* 	  gint atmp = (ustart + a); */
/* 	  for ( b = 0; b < k; b++) { */
/* 	    A0[atmp + (vstart+b)*NUT][indexi] += */
/* 	      alpha*gsl_matrix_get (Bu, a, 2)*gsl_matrix_get (Bv, b, 0) */
/* 	      -2.*beta*gsl_matrix_get (Bu, a, 1)*gsl_matrix_get (Bv, b, 1) */
/* 	      +gamma*gsl_matrix_get (Bu, a, 0)*gsl_matrix_get (Bv, b, 2) */
/* 	      + J*J*(P*gsl_matrix_get (Bu, a, 1)*gsl_matrix_get (Bv, b, 0) */
/* 		     + Q*gsl_matrix_get (Bu, a, 0)*gsl_matrix_get (Bv, b, 1)); */
/* 	  } */
/* 	} */
/* 	// NB: RHS is zero */
/*       } */
/*     } */
    
/*     // Keep boundary points at v = 0 and v = 1 where they are */
/*     for ( i = 0; i < sp->NXU; i++) { */
/*       gdouble u = gsl_bspline_greville_abscissa (i, sp->wx_u); */

/*       gint index1 = sp->fs_index_x + i; // Inner boundary */
/*       gdouble v = 0.; */
      
/*       for ( j = 0; j < size; j++) */
/*     	A0[j][index1] = 0.; */

/*       gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 2, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u); */
/*       gsl_bspline_deriv_eval_nonzero (v, 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v); */

/*       for ( a = 0; a < sp->k; a++) { */
/*     	gdouble cu = gsl_matrix_get (Bu, a, 0); */
/*     	for ( b = 0; b < sp->k; b++) { */
/*     	  gdouble cv = gsl_matrix_get (Bv, b, 0); */
/*     	  A0[(sp->fs_index_x + ustart + a) + (vstart+b)*NUT][index1] = cu*cv; */
/*     	} */
/*       } */

/*       Point p = spline2d_eval_point (sp, u, v); */
/*       RHSX[index1] = p.x; */
/*       RHSY[index1] = p.y; */

 
/*       index1 = sp->fs_index_x + i + (sp->NXV-1)*NUT; // Outer boundary */
/*       v = 1.; */
      
/*       for ( j = 0; j < size; j++) */
/*     	A0[j][index1] = 0.; */

/*       gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 2, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u); */
/*       gsl_bspline_deriv_eval_nonzero (v, 2, Bv, &vstart, &vend, sp->w_v, sp->wd_v); */

/*       for ( a = 0; a < sp->k; a++) { */
/*     	gdouble cu = gsl_matrix_get (Bu, a, 0); */
/*     	for ( b = 0; b < sp->k; b++) { */
/*     	  gdouble cv = gsl_matrix_get (Bv, b, 0); */
/*     	  A0[(sp->fs_index_x + ustart + a) + (vstart+b)*NUT][index1] = cu*cv; */
/*     	} */
/*       } */
/*       p = spline2d_eval_point (sp, u, v); */
/*       RHSX[index1] = p.x; */
/*       RHSY[index1] = p.y; */

      
/*     } */

/*     // Continuity condition + normal derivative continuity between consecutive patches */
/*     Spline2D * sp2 = sp->next == NULL ? splines: sp->next;     */
/*     GrevillePoints * gr = sp->gr; */
/*     for ( j = 0; j < sp->NXV; j++) { */
/*       gint index1 = (sp->fs_index_x + sp->NXU - 1) + j*NUT; */
/*       gint index2 = sp->fs_index_x + j*NUT; */

/*       for ( i = 0; i < size; i++) */
/*     	A0[i][index1] = A0[i][index2] = 0.; */

/*       gdouble u = 1.; */
/*       gdouble v = g_array_index (gr->vj, gdouble, j); */

/*       gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp->wx_u, sp->wxd_u); */
/*       gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &vstart, &vend, sp->w_v, sp->wd_v); */

/*       for ( m = 0; m < sp->k; m++) { */
/*     	gdouble cu = gsl_matrix_get (Bu, m, 0); */
/* 	gdouble cdu = gsl_matrix_get (Bu, m, 1); */
/*     	for ( n = 0; n < sp->k; n++) { */
/*     	  A0[sp->fs_index_x + ustart + m + (vstart+n)*NUT][index1] = */
/*     	    cu*gsl_matrix_get (Bv, n, 0); */
/* 	  A0[sp->fs_index_x + ustart + m + (vstart+n)*NUT][index2] = */
/*     	    cdu*gsl_matrix_get (Bv, n, 0); */
/*     	} */
/*       } */

/*       u = 0.; */
/*       gsl_bspline_deriv_eval_nonzero (MAX(0.+1e-12,MIN(1.-1e-12,u)), 1, Bu, &ustart, &uend, sp2->wx_u, sp2->wxd_u); */
/*       gsl_bspline_deriv_eval_nonzero (v, 1, Bv, &vstart, &vend, sp2->w_v, sp2->wd_v); */
      
/*       for ( m = 0; m < sp2->k; m++) { */
/*     	gdouble cu = gsl_matrix_get (Bu, m, 0); */
/* 	gdouble cdu = gsl_matrix_get (Bu, m, 1); */
/*     	for ( n = 0; n < sp2->k; n++) { */
/*     	  A0[sp2->fs_index_x + ustart + m + (vstart+n)*NUT][index1] -= */
/*     	    cu*gsl_matrix_get (Bv, n, 0); */
/* 	  A0[sp2->fs_index_x + ustart + m + (vstart+n)*NUT][index2] -= */
/*     	    cdu*gsl_matrix_get (Bv, n, 0); */
/*     	} */
/*       } */

/*       RHSX[index1] = RHSY[index1] = 0.; */
/*       RHSX[index2] = RHSY[index2] = 0.; */
/*     } */
  

/*   gsl_matrix_free (Bu); */
/*   gsl_matrix_free (Bv); */

/*   /\* Storage in Compressed Column Storage (CCS) format for superlu *\/ */
/*   CCSProblem * css = ccs_problem_new (); */
/*   gsl_vector * rhs_x = gsl_vector_alloc (size); */
/*   gsl_vector * rhs_y = gsl_vector_alloc (size); */
/*   gint count = 0; */
/*   for ( i = 0; i < size; i++) { */
/*     gsl_vector_set (rhs_x, i, RHSX[i]); */
/*     gsl_vector_set (rhs_y, i, RHSY[i]); */
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

/*   // Solve problem */
/*   ccs_problem_lu_solve (css, rhs_x); */
/*   ccs_problem_lu_solve (css, rhs_y); */

/*   // Store solution */
/*   sp = splines; */
/*   while (sp) { */

/*     for ( i = 0; i < sp->NXU; i++) { */
/*       for ( j = 1; j < sp->NXV-1; j++) { // Don't copy the coefficients on the borders as we don't want them to change at all */
/*   	gint indexi = (sp->fs_index_x + i) + j*NUT; */
/*   	coeff_assign (sp, i, j, 0, gsl_vector_get (rhs_x, indexi)); */
/*   	coeff_assign (sp, i, j, 1, gsl_vector_get (rhs_y, indexi)); */
/*       } */
/*     } */

/*     sp = sp->next; */
/*   } */

/*   gsl_vector_free (rhs_x); */
/*   gsl_vector_free (rhs_y); */
/*   ccs_problem_destroy (css); */
/* } */

gdouble regrid_x (SPPanel * spp, gint m, gint n, gpointer data)
{
  gdouble u = g_array_index (spp->outer->ui, gdouble, m);
  gdouble v = g_array_index (spp->outer->vj, gdouble, n);
  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval (sp, u, v, 0);
}

gdouble regrid_y (SPPanel * spp, gint m, gint n, gpointer data)
{
  gdouble u = g_array_index (spp->outer->ui, gdouble, m);
  gdouble v = g_array_index (spp->outer->vj, gdouble, n);
  
  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval (sp, u, v, 1);
}

gdouble regrid_z (SPPanel * spp, gint m, gint n, gpointer data)
{
  gdouble u = g_array_index (spp->outer->ui, gdouble, m);
  gdouble v = g_array_index (spp->outer->vj, gdouble, n);

  Spline2D * sp = (Spline2D *) data;

  return spline2d_eval (sp, u, v, 2);
}

Spline2D * spline2d_regrid (Spline2D * sp,  gint M, gint N)
{
  Spline2D * grid = spline2d_new (N, M, 3, 4, 3);

  coeff_set_var_to_zero (grid, 0);
  coeff_set_var_to_zero (grid, 1);
  coeff_set_var_to_zero (grid, 2);

  spline2d_init_panels (grid);

  gint size = grid->NU*grid->NV;

  CCSProblem * fit = spline2d_build_galerkin_fit_matrix_no_metric (grid);

  gsl_vector * rhs_x = build_galerkin_rhs_gauss_no_metric (grid, regrid_x, sp, NULL, NULL);
  gsl_vector * rhs_y = build_galerkin_rhs_gauss_no_metric (grid, regrid_y, sp, NULL, NULL);
  gsl_vector * rhs_z = build_galerkin_rhs_gauss_no_metric (grid, regrid_z, sp, NULL, NULL);
  
  // LU decomposition
  ccs_problem_lu_solve (fit, rhs_x);
  ccs_problem_lu_solve (fit, rhs_y);
  ccs_problem_lu_solve (fit, rhs_z);

  gint i, j;
  grid->copy_fit_solution (grid, rhs_x, 0);
  grid->copy_fit_solution (grid, rhs_y, 1);
  grid->copy_fit_solution (grid, rhs_z, 2);

  spline2d_reinit_panels_physical_quantities (grid);

  grid->fit = grid->build_fit_matrix (grid);

  gsl_vector_free (rhs_x);
  gsl_vector_free (rhs_y);
  gsl_vector_free (rhs_z);

  ccs_problem_destroy (fit);

  return grid;
}

void second_derivatives (Spline2D * sp, gdouble u, gdouble v)
{
  Vector grad;
  Vector xu, xv;
  
  xu.x = spline2d_derivative_eval (sp, u, v, 1, 0, 0);
  xu.y = spline2d_derivative_eval (sp, u, v, 1, 0, 1);
  xu.z = spline2d_derivative_eval (sp, u, v, 1, 0, 2);
  xv.x = spline2d_derivative_eval (sp, u, v, 0, 1, 0);
  xv.y = spline2d_derivative_eval (sp, u, v, 0, 1, 1);
  xv.z = spline2d_derivative_eval (sp, u, v, 0, 1, 2);

  /* gsl_matrix * S = gsl_matrix_alloc (6,5); */

  /* gsl_matrix_set (S, 0, 0, xu.x); */
  /* gsl_matrix_set (S, 0, 1, xu.y); */
  /* gsl_matrix_set (S, 0, 2, xu.z); */
  /* gsl_matrix_set (S, 0, 3, 0.); */
  /* gsl_matrix_set (S, 0, 4, 0.); */
  
  /* gsl_matrix_set (S, 1, 0, xv.x); */
  /* gsl_matrix_set (S, 1, 1, xv.y); */
  /* gsl_matrix_set (S, 1, 2, xv.z); */
  /* gsl_matrix_set (S, 1, 3, 0.); */
  /* gsl_matrix_set (S, 1, 4, 0.); */

  /* gsl_matrix_set (S, 2, 0, 0.); */
  /* gsl_matrix_set (S, 2, 1, xu.x); */
  /* gsl_matrix_set (S, 2, 2, 0.); */
  /* gsl_matrix_set (S, 2, 3, xu.y); */
  /* gsl_matrix_set (S, 2, 4, xu.z); */

  /* gsl_matrix_set (S, 3, 0, 0.); */
  /* gsl_matrix_set (S, 3, 1, xv.x); */
  /* gsl_matrix_set (S, 3, 2, 0.); */
  /* gsl_matrix_set (S, 3, 3, xv.y); */
  /* gsl_matrix_set (S, 3, 4, xv.z); */

  /* gsl_matrix_set (S, 4, 0, -xu.z); */
  /* gsl_matrix_set (S, 4, 1, 0.); */
  /* gsl_matrix_set (S, 4, 2, xu.x); */
  /* gsl_matrix_set (S, 4, 3, -xu.z); */
  /* gsl_matrix_set (S, 4, 4, xu.y); */

  /* gsl_matrix_set (S, 5, 0, -xv.z); */
  /* gsl_matrix_set (S, 5, 1, 0.); */
  /* gsl_matrix_set (S, 5, 2, xv.x); */
  /* gsl_matrix_set (S, 5, 3, -xv.z); */
  /* gsl_matrix_set (S, 5, 4, xv.y); */


  /* gsl_matrix_free (S); */
  /* gsl_matrix * ST = gsl_matrix_transpose (S); */

  gdouble S[6][5];

  S[0][0] = xu.x;
  S[0][1] = xu.y;
  S[0][2] = xu.z;
  S[0][3] = 0.;
  S[0][4] = 0.;

  S[1][0] = xv.x;
  S[1][1] = xv.y;
  S[1][2] = xv.z;
  S[1][3] = 0.;
  S[1][4] = 0.;

  S[2][0] = 0.;
  S[2][1] = xu.x;
  S[2][2] = 0.;
  S[2][3] = xu.y;
  S[2][4] = xu.z;

  S[3][0] = 0.;
  S[3][1] = xv.x;
  S[3][2] = 0.;
  S[3][3] = xv.y;
  S[3][4] = xv.z;

  S[4][0] = -xu.z;
  S[4][1] = 0;
  S[4][2] = xu.x;
  S[4][3] = -xu.z;
  S[4][4] = xu.y;

  S[5][0] = -xv.z;
  S[5][1] = 0;
  S[5][2] = xv.x;
  S[5][3] = -xv.z;
  S[5][4] = xv.y;

  /* double s_data[] = {  xu.x , xu.y, xu.z,  0.  , 0.  , */
  /* 		       xv.x , xv.y, xv.z,  0.  , 0.  , */
  /* 		       0.   , xu.x, 0.  ,  xu.y, xu.z, */
  /* 		       0.   , xv.x, 0.  ,  xv.y, xv.z, */
  /* 		       -xu.z, 0.  , xu.x, -xu.z, xu.y, */
  /* 		       -xv.z, 0.  , xv.x, -xv.z, xv.y}; */
     
  gint i, j, k, l;
  gdouble A[5][5];

  for ( i = 0; i < 5; i++ ) {
    for ( k = 0; k < 5; k++ ) {
      A[k][l] = 0.;
    }
  }

  for ( i = 0; i < 5; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      for ( k = 0; k < 6; k++ ) {
	A[i][j] += S[k][i]*S[k][j]; 
	}
      }
  }

  gsl_matrix * A0 = gsl_matrix_alloc (5,5);
  for ( i = 0; i < 5; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      gsl_matrix_set (A0, i, j, A[i][j]);
    }
  }

  gdouble B[6];
  gdouble RHS[5];
  
  for ( i = 0; i < 5; i++)
    RHS[i] = 0.;

  // RHS = S^T B
  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      RHS[j] += S[i][j]*B[i];
    }
  }

  
  gsl_vector * rhs0 = gsl_vector_alloc (5);
  for ( i = 0; i < 5; i++ )
    gsl_vector_set (rhs0, i, RHS[i]);




  gsl_vector * res = gsl_vector_alloc (5);

  int s;
     
  gsl_permutation * p = gsl_permutation_alloc (5);
     
  gsl_linalg_LU_decomp (A0, p, &s);
     
  gsl_linalg_LU_solve (A0, p, res, rhs0);
     
  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 6; j++) {
      A[i][j] = 0.;
    }
  }

  
  // Get dzz from Laplace equation
     
  gsl_permutation_free (p);
  gsl_vector_free (rhs0);
  gsl_vector_free (res);
  gsl_matrix_free (A0);
}

static gdouble gradPhix_rhs (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = gp->ui->len;

  Vector gradPhi = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  return gradPhi.x;
}

static gdouble gradPhiy_rhs (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = gp->ui->len;

  Vector gradPhi = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  return gradPhi.y;
}

static gdouble gradPhiz_rhs (SPPanel * spp, gint m, gint n, gpointer data)
{
  GaussPoints * gp = spp->outer;
  gint ng = gp->ui->len;

  Vector gradPhi = potential_gradient_on_surface_gauss_point (spp->sp, gp, m, n, 3);

  return gradPhi.z;
}

void evaluate_velocity_splines (GSList * patches)
{
  GSList * splines = patches;
  
  while (splines) {
    Spline2D * sp = splines->data;
    while (sp) {
      
      sp->rhs = sp->build_fit_rhs (sp, gradPhix_rhs, NULL, NULL, NULL, sp->rhs);
      ccs_problem_lu_solve (sp->fit, sp->rhs);
      sp->copy_fit_solution (sp, sp->rhs, 24);

      sp->rhs = sp->build_fit_rhs (sp, gradPhiy_rhs, NULL, NULL, NULL, sp->rhs);
      ccs_problem_lu_solve (sp->fit, sp->rhs);
      sp->copy_fit_solution (sp, sp->rhs, 25);

      sp->rhs = sp->build_fit_rhs (sp, gradPhiz_rhs, NULL, NULL, NULL, sp->rhs);
      ccs_problem_lu_solve (sp->fit, sp->rhs);
      sp->copy_fit_solution (sp, sp->rhs, 26);

      sp = sp->next;
    }
    splines = splines->next;
  }
}


gdouble phi_dzz (SPPanel * spp, gint m, gint n, gpointer data)
{
  Vector grad;
  Vector xu, xv;
  Spline2D * sp = spp->sp;

  GaussPoints * gp = spp->outer;
  gint ng = gp->ui->len; // Order of outer Gauss-Legendre rule
  gdouble u = g_array_index (gp->ui, gdouble, m);
  gdouble v = g_array_index (gp->vj, gdouble, n);

  xu.x = spline2d_derivative_eval (sp, u, v, 1, 0, 0);
  xu.y = spline2d_derivative_eval (sp, u, v, 1, 0, 1);
  xu.z = spline2d_derivative_eval (sp, u, v, 1, 0, 2);
  xv.x = spline2d_derivative_eval (sp, u, v, 0, 1, 0);
  xv.y = spline2d_derivative_eval (sp, u, v, 0, 1, 1);
  xv.z = spline2d_derivative_eval (sp, u, v, 0, 1, 2);

  gdouble phidxdu = spline2d_derivative_eval (sp, u, v, 1, 0, 24);
  gdouble phidxdv = spline2d_derivative_eval (sp, u, v, 0, 1, 24);
  gdouble phidydu = spline2d_derivative_eval (sp, u, v, 1, 0, 25);
  gdouble phidydv = spline2d_derivative_eval (sp, u, v, 0, 1, 25);
  gdouble phidzdu = spline2d_derivative_eval (sp, u, v, 1, 0, 26);
  gdouble phidzdv = spline2d_derivative_eval (sp, u, v, 0, 1, 26);

  gdouble S[6][5];

  S[0][0] = xu.x;
  S[0][1] = xu.y;
  S[0][2] = xu.z;
  S[0][3] = 0.;
  S[0][4] = 0.;

  S[1][0] = xv.x;
  S[1][1] = xv.y;
  S[1][2] = xv.z;
  S[1][3] = 0.;
  S[1][4] = 0.;

  S[2][0] = 0.;
  S[2][1] = xu.x;
  S[2][2] = 0.;
  S[2][3] = xu.y;
  S[2][4] = xu.z;

  S[3][0] = 0.;
  S[3][1] = xv.x;
  S[3][2] = 0.;
  S[3][3] = xv.y;
  S[3][4] = xv.z;

  S[4][0] = -xu.z;
  S[4][1] = 0;
  S[4][2] = xu.x;
  S[4][3] = -xu.z;
  S[4][4] = xu.y;

  S[5][0] = -xv.z;
  S[5][1] = 0;
  S[5][2] = xv.x;
  S[5][3] = -xv.z;
  S[5][4] = xv.y;

     
  gint i, j, k, l;
  gdouble A[5][5];

  for ( i = 0; i < 5; i++ ) {
    for ( k = 0; k < 5; k++ ) {
      A[k][i] = 0.;
    }
  }

  for ( i = 0; i < 5; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      for ( k = 0; k < 6; k++ ) {
	A[i][j] += S[k][i]*S[k][j]; 
	}
      }
  }

  /* for ( i = 0; i < 5; i++ ) { */
  /*   for ( j = 0; j < 5; j++ ) { */
  /*     fprintf (stdout, "%f ", A[i][j]); */
  /*   } */
  /*   fprintf (stdout, " \n"); */
  /* } */


  gsl_matrix * A0 = gsl_matrix_alloc (5,5);
  for ( i = 0; i < 5; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      gsl_matrix_set (A0, i, j, A[i][j]);
    }
  }

  /* for ( i = 0; i < 5; i++ ) { */
  /*   for ( j = 0; j < 5; j++ ) { */
  /*     fprintf (stdout, "%f ", gsl_matrix_get (A0,i,j)); */
  /*   } */
  /*   fprintf (stdout, " \n"); */
  /* } */

  gdouble B[6];
  B[0] = phidxdu;
  B[1] = phidxdv;
  B[2] = phidydu;
  B[3] = phidydv;
  B[4] = phidzdu;
  B[5] = phidzdv;

  /* fprintf (stdout, "\n B: "); */
  /* for ( i = 0; i < 6; i++ ) */
  /* fprintf (stdout, "%f ", B[i]); */
  /* fprintf (stdout, "\n"); */

  gdouble RHS[5];
  for ( i = 0; i < 5; i++)
    RHS[i] = 0.;

  // RHS = S^T B
  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      RHS[j] += S[i][j]*B[i];
    }
  }

  /* for ( i = 0; i < 5; i++ ) */
  /* fprintf (stdout, "%f ", RHS[i]); */
  /* fprintf (stdout, "\n"); */


  gsl_vector * rhs0 = gsl_vector_alloc (5);
  for ( i = 0; i < 5; i++ )
    gsl_vector_set (rhs0, i, RHS[i]);


  /* for ( i = 0; i < 5; i++ ) */
  /*   fprintf (stdout, "%f ", gsl_vector_get (rhs0, i)); */
  /* fprintf (stdout, "\n"); */

  gsl_vector * res = gsl_vector_alloc (5);

  int s;
     
  gsl_permutation * p = gsl_permutation_alloc (5);
     
  gsl_linalg_LU_decomp (A0, p, &s);

  gsl_matrix * Bi = gsl_matrix_alloc (5,5);
  gsl_linalg_LU_invert (A0, p, Bi);
     
  /* fprintf (stdout, "Invert \n"); */
  /* for ( i = 0; i < 5; i++ ) { */
  /*   for ( j = 0; j < 5; j++ ) { */
  /*     fprintf (stdout, "%f ", gsl_matrix_get (Bi,i,j)); */
  /*   } */
  /*   fprintf (stdout, " \n"); */
  /* } */
  /* fprintf (stdout, " \n"); */

  gdouble tt[5];
  for ( i = 0; i < 5; i++ )
    tt[i] = 0.;

  for ( i = 0; i < 5; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      tt[i] += gsl_matrix_get (Bi,i,j)*RHS[j];
    }
  }

  /* for ( i = 0; i < 5; i++ ) */
  /*   fprintf (stdout, "%e ", tt[i]); */
  /* fprintf (stdout, "\n"); */
  //  g_assert_not_reached ();

  gsl_linalg_LU_solve (A0, p, rhs0, res);
     
  /* for ( i = 0; i < 5; i++ ) */
  /*   fprintf (stdout, "%e ", gsl_vector_get (res,i)); */
  /* fprintf (stdout, "\n"); */

  
  // Get dzz from Laplace equation
  gdouble dzz = -gsl_vector_get (res, 0) - gsl_vector_get (res, 3);;
  /* Vector grad2 = potential_gradient_on_surface_gauss_point (sp, gp, m, n, 14); */
    
  /* fprintf (stdout, "\n %f %f %f \n", dzz, spline2d_eval (sp,u,v,15), grad2.z); */
  /* g_assert_not_reached (); */
  
  gsl_permutation_free (p);
  gsl_vector_free (rhs0);
  gsl_vector_free (res);
  gsl_matrix_free (A0);

  return dzz;
}

void calculate_dzz_direct (GSList * patches)
{
  GSList * splines = patches;

  evaluate_velocity_splines (patches);   
  
  while (splines) {
    Spline2D * sp = splines->data;
    while (sp) {
      
      sp->rhs = sp->build_fit_rhs (sp, phi_dzz, NULL, NULL, NULL, sp->rhs);
      ccs_problem_lu_solve (sp->fit, sp->rhs);
      sp->copy_fit_solution (sp, sp->rhs, 25);

      sp = sp->next;
    }
    splines = splines->next;
  }

}

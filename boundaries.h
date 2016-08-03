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

typedef Point (* BoundaryCurve)  (gdouble u, gdouble t, gpointer data);

/* Discrete and spline representation of the boundaries */
typedef struct {
  DCurve * dcb, * dct, * dcl, * dcr;
  BoundaryCurve curve_bottom;
  BoundaryCurve curve_top;
  BoundaryCurve curve_right;
  BoundaryCurve curve_left;
} Boundaries;

Boundaries *  boundaries_new                   ();
void          boundaries_init                  (Boundaries * b, gdouble time,
						gpointer data, gint N, gint M);
void          boundaries_print                 (Boundaries * b, FILE *fp);
void          boundaries_destroy               (Boundaries * b);

DCurve *      hybrid_curve_point_distribution  (gint m, BoundaryCurve curve,
						gdouble t, gpointer * data);

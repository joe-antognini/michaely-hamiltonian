/* -*- linux-c -*- */
/* fewbody_scat.c

   Copyright (C) 2002-2004 John M. Fregeau
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "fewbody.h"

/* JMA 10-14-13 -- Calculate the distance at which the pertpurber exceeds
 * the tidal tolerance. */
double fb_calc_rtid(fb_obj_t *interloper, fb_obj_t *bin, double tidaltol)
{
  double a, r, m1, m2, M;
  
  M = interloper->m;
  m1 = bin->obj[0]->m;
  m2 = bin->obj[1]->m;
  a = bin->a;
  r = sqrt(fb_sqr(bin->x[0] - interloper->x[0]) + \
           fb_sqr(bin->x[1] - interloper->x[1]) + \
           fb_sqr(bin->x[2] - interloper->x[2]));

  return tidaltol / (1. / (r - 2 * m1 * a / M) - 1. / \
          (r + 2 * m2 * a / M));
}

/* move objects in from infinity analytically along a hyperbolic orbit */
void fb_init_scattering(fb_obj_t *obj[2], double vinf, double b, double rtid)
{
  double m0, m1, M, r, rperi, x, y, vx, vy, qa, qb, qc, rad;
  
  /* some useful variables */
  m0 = obj[0]->m;
  m1 = obj[1]->m;
  M = m0 + m1;

  /* pericenter */
  rperi = (vinf==0.0?0.0:(M/fb_sqr(vinf)*(sqrt(1.0+fb_sqr(b*fb_sqr(vinf)/M))-1.0)));

  /* make sure r>=rperi, otherwise analytically moving the obj's below will give NANs */
  r = FB_MAX(rtid, rperi);

  /* solve for y */
  qa = M*M + fb_sqr(b*vinf*vinf);
  qb = 2.0*(M*r-fb_sqr(b*vinf))*b*vinf*vinf;
  qc = -2.0*fb_sqr(b*vinf)*M*r + fb_sqr(fb_sqr(b*vinf));
  rad = qb*qb-4.0*qa*qc;
  y = (-qb+(rad<=0.0?0.0:sqrt(rad)))/(2.0*qa);

  /* x can sometimes be zero, so need to worry about the square root here */
  x = (r*r-y*y<=0.0?0.0:sqrt(r*r-y*y));

  /* determine v_x from quadratic */
  if (x > 0.0) {
    qa = 1.0 + fb_sqr(y/x);
    qb = 2.0*y*b*vinf/(x*x);
    qc = fb_sqr(b*vinf/x) - 2.0*M/r - vinf*vinf;
    rad = qb*qb-4.0*qa*qc;
    vx = (-qb-(rad<=0.0?0.0:sqrt(rad)))/(2.0*qa);
    vy = (b*vinf+y*vx)/x;
  } else { /* x=0 */
    vx = -vinf;
    vy = 0.0;
  }

  /* set the positions and velocities */
  obj[0]->x[0] = - x / (1.0 + m0/m1);
  obj[0]->x[1] = - y / (1.0 + m0/m1);
  obj[0]->x[2] = 0.0;
  
  obj[0]->v[0] = - vx / (1.0 + m0/m1);
  obj[0]->v[1] = - vy / (1.0 + m0/m1);
  obj[0]->v[2] = 0.0;

  obj[1]->x[0] = x / (1.0 + m1/m0);
  obj[1]->x[1] = y / (1.0 + m1/m0);
  obj[1]->x[2] = 0.0;
  
  obj[1]->v[0] = vx / (1.0 + m1/m0);
  obj[1]->v[1] = vy / (1.0 + m1/m0);
  obj[1]->v[2] = 0.0;
}

/* JMA 10-8-2013 */
/* orient the interloper in a random direction */
void fb_randscat(fb_obj_t *obj, gsl_rng *rng)
{
  int i;
  double u, theta, phi;
  double k[3], kcrossx[3], kcrossv[3], xrot[3], vrot[3];

  /* pick a random vector -- see Wolfram MathWorld for eq. */
  u = 2. * gsl_rng_uniform(rng) - 1.;
  theta = 2 * FB_CONST_PI * gsl_rng_uniform(rng);
  k[0] = sqrt(1. - fb_sqr(u)) * cos(theta);
  k[1] = sqrt(1. - fb_sqr(u)) * sin(theta);
  k[2] = u;

  /* now pick a random angle through which to rotate. */
  phi = 2 * FB_CONST_PI * gsl_rng_uniform(rng);

  /* some preparations for the Rodrigues' rotation formula */
  fb_cross(k, obj->x, kcrossx);
  fb_cross(k, obj->v, kcrossv);

  /* use the Rodrigues' rotation formula */
  for (i=0; i<3; i++) {
    xrot[i] = obj->x[i] * cos(phi);
    xrot[i] += kcrossx[i] * sin(phi);
    xrot[i] += k[i] * fb_dot(k, obj->x) * (1 - cos(phi));
  }

  for (i=0; i<3; i++) {
    vrot[i] = obj->v[i] * cos(phi);
    vrot[i] += kcrossv[i] * sin(phi);
    vrot[i] += k[i] * fb_dot(k, obj->v) * (1 - cos(phi));
  }

  for (i=0; i<3; i++) {
    obj->x[i] = xrot[i];
    obj->v[i] = vrot[i];
  }
}

/* normalize all variables to our system of units */
void fb_normalize(fb_hier_t *hier, fb_units_t units)
{
  int i, k;

  /* just normalize everything, since it can't hurt */
  for (i=hier->hi[1]; i<=hier->hi[hier->nstarinit]; i++) {
    hier->hier[i].m /= units.m;
    hier->hier[i].R /= units.l;
    hier->hier[i].Eint /= units.E;
    hier->hier[i].a /= units.l;
    hier->hier[i].t /= units.t;
    for (k=0; k<3; k++) {
      hier->hier[i].x[k] /= units.l;
      hier->hier[i].v[k] /= units.v;
      hier->hier[i].Lint[k] /= units.m * units.l * units.v;
    }
  }
}

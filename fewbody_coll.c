/* -*- linux-c -*- */
/* fewbody_coll.c

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

/* the main collision criterion */
int fb_is_collision(double r, double R1, double R2) {
  /* JMA 8-2-2012 -- I am altering the merger criterion here.  Pray I do
   * not alter it further.
   *
   * Basically, the code spends a long time on the last part of the
   * inspiral.  This is computationally intensive but does not add much
   * time to the merger timescale.  Once the semi-major axis of the inner
   * binary has reached 1% of its original value, the rest of the merger is
   * essentially instantaneous. 
   */
  
  /*
   * lol jk
   */

  if (r < R1 + R2) {
    return(1);
  } else {
    return(0);
  }
}

int fb_collide(fb_hier_t *hier, double f_exp, fb_units_t units, gsl_rng *rng, double *t)
{
  int i, j=-1, k, retval=0, cont=1, sma_cont=1;
  double R[3], peinit;

  /* this is a non-recursive way to perform a recursive operation: keep going until there are no more
     mergers */
  while (cont) {
    cont = 0;
    sma_cont = 0;
    for (i=0; i<hier->nstar-1; i++) {
      for (j=i+1; j<hier->nstar; j++) {
        /* calculate relative separation */
        for (k=0; k<3; k++) {
          R[k] = hier->hier[hier->hi[1]+i].x[k] - hier->hier[hier->hi[1]+j].x[k];
        }

        /* test collision criterion */
        if (fb_is_collision(fb_mod(R), hier->hier[hier->hi[1]+i].R, hier->hier[hier->hi[1]+j].R)) {
          /* JMA 11-12-2012 -- Print the state of the system right before
           * merger.
           */
          fprintf(stdout, "%.12f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", *t,
            hier->hier[hier->hi[2]+0].a, hier->hier[hier->hi[2]+0].e,
            hier->hier[hier->hi[3]+0].a, hier->hier[hier->hi[3]+0].e,
            fb_dot(hier->hier[hier->hi[2]+0].Lhat, hier->hier[hier->hi[2]+1].Lhat),
            acos(fb_dot(hier->hier[hier->hi[2]+0].Ahat, hier->hier[hier->hi[2]+1].Ahat)) * 180 / FB_CONST_PI,
            hier->hier[hier->hi[2]+0].Lhat[0],
            hier->hier[hier->hi[2]+0].Lhat[1],
            hier->hier[hier->hi[2]+0].Lhat[2],
            hier->hier[hier->hi[1]+1].x[0] - hier->hier[hier->hi[1]+0].x[0], 
            hier->hier[hier->hi[1]+1].x[1] - hier->hier[hier->hi[1]+0].x[1],
            hier->hier[hier->hi[1]+1].x[2] - hier->hier[hier->hi[1]+0].x[2],
            hier->hier[hier->hi[1]+2].x[0],
            hier->hier[hier->hi[1]+2].x[1],
            hier->hier[hier->hi[1]+2].x[2],
            hier->hier[hier->hi[1]+1].v[0] - hier->hier[hier->hi[1]+0].v[0], 
            hier->hier[hier->hi[1]+1].v[1] - hier->hier[hier->hi[1]+0].v[1],
            hier->hier[hier->hi[1]+1].v[2] - hier->hier[hier->hi[1]+0].v[2]
            );

          cont = 1;
          /* break out of the double loop if there is a collision, so we can merge
             the stars immediately */
          break;
        }
      }
      /* break out of the double loop if there is a collision, so we can merge
         the stars immediately */
      if (cont) {
        break;
      }
    }

    /* merge two stars if necessary */
    if (cont) {
      /* return 1 if there is a collision */
      retval = 1;

      /* calculate the potential energy before the collision, since we're going to need
         to account for the change in potential energy, and put it in Eint */
      peinit = fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);      
      
      /* do the actual merger */
      fb_dprintf("fewbody: collide(): merging stars: i=%d j=%d\n", i, j);
      fb_merge(&(hier->hier[hier->hi[1]+i]), &(hier->hier[hier->hi[1]+j]), hier->nstarinit, f_exp, units, rng);
      fb_objcpy(&(hier->hier[hier->hi[1]+j]), &(hier->hier[hier->hi[1]+hier->nstar-1]));
      hier->nstar--;

      /* calculate the difference in potential energy before and after the collision, and put it
         in Eint for accounting */
      hier->hier[hier->hi[1]+i].Eint += peinit - fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);
    }

    /* JMA 8-2-12 -- I am adding a merging criterion.  If the semi-major
     * axis of a binary is less than 1% of the original semi-major axis of
     * the inner binary, the system merges.  
     * 
     * I am adding this criterion because the final steps of the inspiral
     * are very computationally expensive but do not add much time to the
     * total merger time.  Once the inner orbit reaches 1% of its original
     * value, the change in the merger time is about one part in one
     * million.
     *
     * WARNING: THIS WILL ONLY WORK ON TRIPLE SYSTEMS!
     */ 

    /* JMA 8-21-12 -- This criterion is being modified slightly.  Because
     * the energy calcuulations are classical, the semi-major axis
     * calculation is wrong when the orbit is extremely eccentric.  This
     * leads to large underestimates of the semi-major axis which causes
     * the code to terminate prematurely.  By requiring that the orbit be
     * somewhat circular when we stop, we can prevent this from happening.
     * This is something of a kludge, but it will work for now until the
     * energy calculations are done correctly.
     */
    
    //if (hier->hier[hier->hi[2]+0].a < .01 && hier->hier[hier->hi[2]+0].e < .9) {
    //  sma_cont = 1;
    //}

    /*
    if (hier->hier[hier->hi[2]+0].a < .03) {
      sma_cont = 1;
    }
    */
    
    /*
    for (i=0; i < hier->nstar - 1; i++) {
      if (hier->hier[hier->hi[2]+i].a < .01) {
        cont = 1;
        break;
      }
    }
    */
    
    if (sma_cont) {
      /* return 1 if there is a collision */
      retval = 1;

      /* calculate the potential energy before the collision, since we're going to need
         to account for the change in potential energy, and put it in Eint */
      peinit = fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);      
      
      /* do the actual merger */
      //fb_dprintf("fewbody: collide(): merging stars: i=%d j=%d\n", i, j);
      fb_merge(&(hier->hier[hier->hi[1]+0]), &(hier->hier[hier->hi[1]+1]), hier->nstarinit, f_exp, units, rng);
      fb_objcpy(&(hier->hier[hier->hi[1]+1]), &(hier->hier[hier->hi[1]+hier->nstar-1]));
      hier->nstar--;

      /* calculate the difference in potential energy before and after the collision, and put it
         in Eint for accounting */
      hier->hier[hier->hi[1]+0].Eint += peinit - fb_petot(&(hier->hier[hier->hi[1]]), hier->nstar);
    }

    /* END JMA */

  }

  return(retval);
}

void fb_merge(fb_obj_t *obj1, fb_obj_t *obj2, int nstarinit, double f_exp, fb_units_t units, gsl_rng *rng)
{
  int i;
  double x1[3], x2[3], v1[3], v2[3], l1[3], l2[3];
  double clight, vrel[3], x[3], y[3], z[3], theta, vkickmod;
  fb_obj_t tmpobj;
  
  clight = FB_CONST_C / units.v;

  /* sanity check */
  if (obj1->n != 1 || obj2->n != 1) {
    fprintf(stderr, "fb_merge: trying to merge an object that isn't single?!\n");
    exit(1);
  }

  /* need some temporary storage here */
  tmpobj.id = (long *) malloc(nstarinit * sizeof(long));

  /* merge id's */
  tmpobj.ncoll = obj1->ncoll + obj2->ncoll;
  for (i=0; i<obj1->ncoll; i++) {
    tmpobj.id[i] = obj1->id[i];
  }
  for (i=0; i<obj2->ncoll; i++) {
    tmpobj.id[obj1->ncoll + i] = obj2->id[i];
  }

  /* create idstring */
  snprintf(tmpobj.idstring, FB_MAX_STRING_LENGTH, "%s:%s", obj1->idstring, obj2->idstring);

  /* assume no mass loss */
  tmpobj.m = obj1->m + obj2->m;

  /* this is just a simple prescription */
  /* tmpobj.R = f_exp * (obj1->R + obj2->R); */
  tmpobj.R = FB_REFF_BH2 * tmpobj.m / fb_sqr(clight);

  /* set new position and velocity, calculate relative positions and velocities */
  for (i=0; i<3; i++) {
    tmpobj.x[i] = (obj1->m * obj1->x[i] + obj2->m * obj2->x[i]) / tmpobj.m;
    tmpobj.v[i] = (obj1->m * obj1->v[i] + obj2->m * obj2->v[i]) / tmpobj.m;
    x1[i] = obj1->x[i] - tmpobj.x[i];
    x2[i] = obj2->x[i] - tmpobj.x[i];
    v1[i] = obj1->v[i] - tmpobj.v[i];
    v2[i] = obj2->v[i] - tmpobj.v[i];
  }

  /* set internal energy, using the difference in kinetic energy; the difference in potential energy
     depends on the positions of the other stars, and will be calculated later and added to Eint */
  tmpobj.Eint = obj1->Eint + obj2->Eint +
    0.5 * (obj1->m * fb_dot(obj1->v, obj1->v) + obj2->m * fb_dot(obj2->v, obj2->v)) -
    0.5 * tmpobj.m * fb_dot(tmpobj.v, tmpobj.v);

  /* set internal angular momentum */
  fb_cross(x1, v1, l1);
  fb_cross(x2, v2, l2);
  for (i=0; i<3; i++) {
    tmpobj.Lint[i] = obj1->Lint[i] + obj2->Lint[i] + obj1->m * l1[i] + obj2->m * l2[i];
  }
  
  /* DEBUG: radiation rocket */
  /* determine random angle in orbital plane for radiation rocket kick */
  /* first calculate relative velocity, which we'll use for the y-vector */
  for (i=0; i<3; i++) {
    vrel[i] = obj1->v[i] - obj2->v[i];
  }
  /* the z-vector is the unit vector pointing in the direction of the angular momentum */
  for (i=0; i<3; i++) {
    y[i] = vrel[i]/fb_mod(vrel);
    z[i] = tmpobj.Lint[i]/fb_mod(tmpobj.Lint);
  }
  /* x equals y cross z */
  fb_cross(y, z, x);
  /* fprintf(stdout, "x=%g %g %g  y=%g %g %g  z=%g %g %g\n", 
     x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]); */
  /* we randomize the angle in the orbital plane */
  theta = 2.0 * FB_CONST_PI * gsl_rng_uniform(rng);
  /* the kick magnitude is determined from this subroutine */
  vkickmod = fb_vkick(obj1->m, obj2->m);
  /* and here we add the kick */
  for (i=0; i<3; i++) {
    tmpobj.v[i] += vkickmod/units.v * (cos(theta) * x[i] + sin(theta) * y[i]);
  }
  /* DEBUG */

  /* and better set these, too... */
  tmpobj.n = 1;
  
  tmpobj.obj[0] = NULL;
  tmpobj.obj[1] = NULL;

  /* finally, copy over the merger from temporary storage */
  fb_objcpy(obj1, &tmpobj);

  free(tmpobj.id);
}

/* radiation rocket kick speed; output is in CGS; input is in arbitrary units */
double fb_vkick(double m1, double m2)
{
  double q, fmax=0.01788854382;

  q =  FB_MIN(m1, m2) / FB_MAX(m1, m2);

  return(FB_VKICK * fb_sqr(q) * (1.0 - q) / pow(1.0 + q, 5) / fmax);
}

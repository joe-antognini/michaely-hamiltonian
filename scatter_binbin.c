/* -*- linux-c -*- */
/* scatter_binbin.c

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
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <fcntl.h>
#include <unistd.h>
//#include <assert.h>
#include "fewbody.h"
#include "scatter_binbin.h"

/* print the usage */
void print_usage(FILE *stream)
{
  fprintf(stream, "USAGE:\n");
  fprintf(stream, "  scatter_binbin [options...]\n");
  fprintf(stream, "\n");
  fprintf(stream, "OPTIONS:\n");
  fprintf(stream, "  -m --m00 <m00/MSUN>          : set mass of star 0 of of binary 0 [%.6g]\n", FB_M00/FB_CONST_MSUN);
  fprintf(stream, "  -n --m01 <m01/MSUN>          : set mass of star 1 of of binary 0 [%.6g]\n", FB_M01/FB_CONST_MSUN);
  fprintf(stream, "  -M --m10 <m10/MSUN>          : set mass of star 0 of of binary 1 [%.6g]\n", FB_M10/FB_CONST_MSUN);
  fprintf(stream, "  -o --m11 <m11/MSUN>          : set mass of star 1 of of binary 1 [%.6g]\n", FB_M10/FB_CONST_MSUN);
  fprintf(stream, "  -r --r00 <r00/R_SUN>         : set merge radius of star 1 of binary 0 [%.6g]\n", FB_R00/FB_CONST_RSUN);
  fprintf(stream, "  -H --r01 <r01/R_SUN>         : set merge radius of star 0 of binary 0 [%.6g]\n", FB_R01/FB_CONST_RSUN);
  fprintf(stream, "  -C --r10 <r10/R_SUN>         : set merge radius of star 0 of binary 1 [%.6g]\n", FB_R10/FB_CONST_RSUN);
  fprintf(stream, "  -W --r11 <r11/R_SUN>         : set merge radius of star 1 of binary 1  [%.6g]\n", FB_R11/FB_CONST_RSUN);
  fprintf(stream, "  -a --a0 <a0/AU>              : set semimajor axis of binary 0 [%.6g]\n", FB_A0/FB_CONST_AU);
  fprintf(stream, "  -q --a1 <a1/AU>              : set semimajor axis of binary 1 [%.6g]\n", FB_A1/FB_CONST_AU);
  fprintf(stream, "  -e --e0 <e0>                 : set eccentricity of binary 0 [%.6g]\n", FB_E0);
  fprintf(stream, "  -F --e1 <e1>                 : set eccentricity of binary 1 [%.6g]\n", FB_E1);
  fprintf(stream, "  -b --bmax <b/AU>             : set maximum impact parameter (-1 for Hut & Bahcall's choice) [%.6g]\n", FB_BMAX);
  fprintf(stream, "  -X --hutfactor <X>           : set factor times Hut & Bahcall's choice for bmax [%.6g]\n", FB_XBMAX);
  fprintf(stream, "  -v --vinf <vinf/km/s>        : set velocity at infinity [%.6g]\n", FB_VINF/FB_CONST_KMS);
  fprintf(stream, "  -t --tstop <tstop/t_cross>   : set stopping time in units of the overall crossing time [%d]\n", FB_TSTOP);
  fprintf(stream, "  -D --dt <dt/t_dyn>           : set approximate output dt [%.6g]\n", FB_DT);
  fprintf(stream, "  -c --tcpustop <tcpustop/sec> : set cpu stopping time [%.6g]\n", FB_TCPUSTOP);
  fprintf(stream, "  -A --absacc <absacc>         : set integrator's absolute accuracy [%.6g]\n", FB_ABSACC);
  fprintf(stream, "  -R --relacc <relacc>         : set integrator's relative accuracy [%.6g]\n", FB_RELACC);
  fprintf(stream, "  -N --ncount <ncount>         : set number of integration steps between calls\n");
  fprintf(stream, "                                 to fb_classify() [%d]\n", FB_NCOUNT);
  fprintf(stream, "  -O --outputfreq <outputfreq> : set the output frequency (-1 for no output) [%d]\n", FB_OUTFREQ);
  fprintf(stream, "  -z --tidaltol <tidaltol>     : set tidal tolerance [%.6g]\n", FB_TIDALTOL);
  fprintf(stream, "  -y --speedtol <speedtol>     : set speed tolerance [%.6g]\n", FB_SPEEDTOL);
  fprintf(stream, "  -P --PN1 <PN1>               : PN1 terms on? [%d]\n", FB_PN1);
  fprintf(stream, "  -Q --PN2 <PN2>               : PN2 terms on? [%d]\n", FB_PN2);
  fprintf(stream, "  -S --PN25 <PN25>             : PN2.5 terms on? [%d]\n", FB_PN25);
  fprintf(stream, "  -T --PN3 <PN3>               : PN3 terms on? [%d]\n", FB_PN3);
  fprintf(stream, "  -U --PN35 <PN35>             : PN3.5 terms on? [%d]\n", FB_PN35);
  fprintf(stream, "  -x --fexp <f_exp>            : set expansion factor of merger product [%.6g]\n", FB_FEXP);
  fprintf(stream, "  -k --ks                      : turn K-S regularization on or off [%d]\n", FB_KS);
  fprintf(stream, "  -s --seed                    : set random seed (0 to sample from /dev/urandom) [%ld]\n", FB_SEED);
  fprintf(stream, "  -d --debug                   : turn on debugging\n");
  fprintf(stream, "  -V --version                 : print version info\n");
  fprintf(stream, "  -h --help                    : display this help text\n");
}

/* calculate the units used */
int calc_units(fb_obj_t *obj[1], fb_units_t *units)
{
  double m0, m00, m01, a0;

  m0 = obj[0]->m;
  m00 = obj[0]->obj[0]->m;
  m01 = obj[0]->obj[1]->m;

  a0 = obj[0]->a;
  
  /* Unit of velocity is approximate relative orbital speed of inner binary,
     unit of length is semimajor axis of inner binary; 
     therefore, unit of time is approximately 1 inner orbital period. */
  units->v = sqrt(FB_CONST_G*(m00+m01)/a0);
  units->l = a0;
  units->t = units->l / units->v;
  units->m = units->l * fb_sqr(units->v) / FB_CONST_G;
  units->E = units->m * fb_sqr(units->v);
  
  return(0);
}

/* the main attraction */
int main(int argc, char *argv[])
{
  int i, j;
  unsigned long int seed, input_seed;
  double m00, m01, m10, m11, r00, r01, r10, r11, a0, a1, e0, e1;
  double input_bmax, bmax, input_vinf, vinf, vcrit, b, rtid, xbmax;
  int linecount, input_tstop;
  int random_data;
  ssize_t result;
  double Ei, Lint[3], Li[3], t;
  fb_hier_t hier;
  fb_input_t input;
  fb_ret_t retval;
  fb_units_t units;
  char string1[FB_MAX_STRING_LENGTH], string2[FB_MAX_STRING_LENGTH];
  gsl_rng *rng;
  const gsl_rng_type *rng_type=gsl_rng_mt19937;
  const char *short_opts = "m:n:M:o:r:H:W:C:a:q:e:F:b:X:v:t:D:c:A:R:N:O:z:y:P:Q:S:T:U:x:k:s:dVh";
  const struct option long_opts[] = {
    {"m00", required_argument, NULL, 'm'},
    {"m01", required_argument, NULL, 'n'},
    {"m10", required_argument, NULL, 'M'},
    {"m11", required_argument, NULL, 'o'},
    {"r00", required_argument, NULL, 'r'},
    {"r01", required_argument, NULL, 'H'},
    {"r10", required_argument, NULL, 'C'},
    {"r11", required_argument, NULL, 'W'},
    {"a0", required_argument, NULL, 'a'},
    {"a1", required_argument, NULL, 'q'},
    {"e0", required_argument, NULL, 'e'},
    {"e1", required_argument, NULL, 'F'},
    {"bmax", required_argument, NULL, 'b'},
    {"xbmax", required_argument, NULL, 'X'},
    {"vinf", required_argument, NULL, 'v'},
    {"tstop", required_argument, NULL, 't'},
    {"dt", required_argument, NULL, 'D'},
    {"tcpustop", required_argument, NULL, 'c'},
    {"absacc", required_argument, NULL, 'A'},
    {"relacc", required_argument, NULL, 'R'},
    {"ncount", required_argument, NULL, 'N'},
    {"outputfreq", required_argument, NULL, 'O'},
    {"tidaltol", required_argument, NULL, 'z'},
    {"ks", required_argument, NULL, 'k'},
    {"seed", required_argument, NULL, 's'},
    {"debug", no_argument, NULL, 'd'},
    {"version", no_argument, NULL, 'V'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}
  };

  /* set parameters to default values */
  m00 = FB_M00;
  m01 = FB_M01;
  m10 = FB_M10;
  m11 = FB_M11;
  r00 = FB_R00;
  r01 = FB_R01;
  r10 = FB_R10;
  r11 = FB_R11;
  a0 = FB_A0;
  a1 = FB_A1;
  e0 = FB_E0;
  e1 = FB_E1;
  input_bmax = FB_BMAX;
  xbmax = FB_XBMAX;
  input_vinf = FB_VINF;
  input.ks = FB_KS;
  input_tstop = FB_TSTOP;
  input.Dflag = 0;
  input.dt = FB_DT;
  input.tcpustop = FB_TCPUSTOP;
  input.absacc = FB_ABSACC;
  input.relacc = FB_RELACC;
  input.ncount = FB_NCOUNT;
  input.outfreq = FB_OUTFREQ;
  input.tidaltol = FB_TIDALTOL;
  input.fexp = FB_FEXP;
  input_seed = FB_SEED;
  input.speedtol = FB_SPEEDTOL;
  input.PN1 = FB_PN1;
  input.PN2 = FB_PN2;
  input.PN25 = FB_PN25;
  input.PN3 = FB_PN3;
  input.PN35 = FB_PN35;
  input.econs = FB_ECONS;
  input.lcons = FB_LCONS;
  fb_debug = FB_DEBUG;
  
  while ((i = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
    switch (i) {
    case 'm':
      m00 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'n':
      m01 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'M':
      m10 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'o':
      m11 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'r':
      r00 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'H':
      r01 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'C':
      r10 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'W':
      r11 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'a':
      a0 = atof(optarg) * FB_CONST_AU;
      break;
    case 'q':
      a1 = atof(optarg) * FB_CONST_AU;
      break;
    case 'e':
      e0 = atof(optarg);
      if (e0 >= 1.0) {
        fprintf(stderr, "e0 must be less than 1\n");
        return(1);
      }
      break;
    case 'F':
      e1 = atof(optarg);
      if (e1 >= 1.0) {
        fprintf(stderr, "e1 must be less than 1\n");
        return(1);
      }
      break;
    case 'b':
      input_bmax = atof(optarg) * FB_CONST_AU;
      if (input_bmax < 0 && input_bmax != -1) {
        fprintf(stderr, "bmax must be positive, or -1 for Hut & Bahcall's choice.\n");
        return(1);
      }
      break;
    case 'X':
      xbmax = atof(optarg);
      break;
    case 'v':
      input_vinf = atof(optarg) * FB_CONST_KMS;
      if (input_vinf < 0) {
        fprintf(stderr, "vinf must be positive.\n");
        return(1);
      }
      break;
    case 't':
      input_tstop = atof(optarg);
      break;
    case 'D':
      input.Dflag = 1;
      input.dt = atof(optarg);
      break;
    case 'c':
      input.tcpustop = atof(optarg);
      break;
    case 'A':
      input.absacc = atof(optarg);
      break;
    case 'R':
      input.relacc = atof(optarg);
      break;
    case 'N':
      input.ncount = atoi(optarg);
      break;
    case 'O':
      input.outfreq = atoi(optarg);
      break;
    case 'z':
      input.tidaltol = atof(optarg);
      break;
    case 'x':
      input.fexp = atof(optarg);
      break;
    case 'y':
      input.speedtol = atof(optarg);
      break;
    case 'P':
      input.PN1 = atoi(optarg);
      break;
    case 'Q':
      input.PN2 = atoi(optarg);
      break;
    case 'S':
      input.PN25 = atoi(optarg);
      break;
    case 'T':
      input.PN3 = atoi(optarg);
      break;
    case 'U':
      input.PN35 = atoi(optarg);
      break;
    case 'k':
      input.ks = atoi(optarg);
      break;
    case 's':
      input_seed = strtoul(optarg, NULL, 0);
      break;
    case 'd':
      fb_debug = 1;
      break;
    case 'V':
      fb_print_version(stdout);
      return(0);
    case 'h':
      fb_print_version(stdout);
      fprintf(stdout, "\n");
      print_usage(stdout);
      return(0);
    default:
      break;
    }
  }
  
  /* check to make sure there was nothing crazy on the command line */
  if (optind < argc) {
    print_usage(stdout);
    return(1);
  }

  // JMA 10-26-2013 -- If no seed given, draw random bits from
  // /dev/urandom.
  if (input_seed == FB_SEED) {
    random_data = open("/dev/urandom", O_RDONLY);
    result = read(random_data, &seed, sizeof seed);
    close(random_data);
  } else {
    seed = input_seed;
  }

  /* initialize a few things for integrator */
  t = 0.0;
  hier.nstarinit = 4;
  hier.nstar = 4;
  fb_malloc_hier(&hier);
  fb_init_hier(&hier);

  /* put stuff in log entry */
  snprintf(input.firstlogentry, FB_MAX_LOGENTRY_LENGTH, "  command line:");
  for (i=0; i<argc; i++) {
    snprintf(&(input.firstlogentry[strlen(input.firstlogentry)]), 
       FB_MAX_LOGENTRY_LENGTH-strlen(input.firstlogentry), " %s", argv[i]);
  }
  snprintf(&(input.firstlogentry[strlen(input.firstlogentry)]),
     FB_MAX_LOGENTRY_LENGTH-strlen(input.firstlogentry), "\n");

  /* initialize GSL rng */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(rng, seed);

  /* print out values of paramaters */
  fprintf(stderr, "PARAMETERS:\n");
  fprintf(stderr, "  ks=%d  seed=%lu\n", input.ks, seed);
  fprintf(stderr, "  a0=%.6g AU  e0=%.6g  m00=%.6g MSUN  m01=%.6g MSUN  r00=%.6g RSUN  r01=%.6g RSUN\n", \
    a0/FB_CONST_AU, e0, m00/FB_CONST_MSUN, m01/FB_CONST_MSUN, r00/FB_CONST_RSUN, r01/FB_CONST_RSUN);
  fprintf(stderr, "  a1=%.6g AU  e1=%.6g  m10=%.6g MSUN  m11=%.6g MSUN\n", \
    a1/FB_CONST_AU, e1, m01/FB_CONST_MSUN, m11/FB_CONST_MSUN);
  fprintf(stderr, "  tidaltol=%.6g  speedtol=%.6g  abs_acc=%.6g rel_acc=%.6g  ncount=%d  fexp=%.6g  outfreq=%d\n", \
    input.tidaltol, input.speedtol, input.absacc, input.relacc, input.ncount, input.fexp, input.outfreq);
  fprintf(stderr, "  PN1=%d  PN2=%d  PN25=%d  PN3=%d  PN35=%d\n", \
    input.PN1, input.PN2, input.PN25, input.PN3, input.PN35);


  /* create hierarchies */
  hier.narr[2] = 2;
  hier.narr[3] = 0;
  hier.narr[4] = 0;
  /* binary 0 */
  hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+0]);
  hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+1]);
  hier.hier[hier.hi[2]+0].t = t;
  // JMA 12-12-13 -- Interloping binary
  hier.hier[hier.hi[2]+1].obj[0] = &(hier.hier[hier.hi[1]+2]);
  hier.hier[hier.hi[2]+1].obj[1] = &(hier.hier[hier.hi[1]+3]);
  hier.hier[hier.hi[2]+1].t = t;

  /* give the objects some properties */
  for (j=0; j<hier.nstar; j++) {
    hier.hier[hier.hi[1]+j].ncoll = 1;
    hier.hier[hier.hi[1]+j].id[0] = j;
    snprintf(hier.hier[hier.hi[1]+j].idstring, FB_MAX_STRING_LENGTH, "%d", j);
    hier.hier[hier.hi[1]+j].n = 1;
    hier.hier[hier.hi[1]+j].obj[0] = NULL;
    hier.hier[hier.hi[1]+j].obj[1] = NULL;
    hier.hier[hier.hi[1]+j].Eint = 0.0;
    hier.hier[hier.hi[1]+j].Lint[0] = 0.0;
    hier.hier[hier.hi[1]+j].Lint[1] = 0.0;
    hier.hier[hier.hi[1]+j].Lint[2] = 0.0;
  }

  hier.hier[hier.hi[1]+0].R = r00;
  hier.hier[hier.hi[1]+1].R = r01;
  hier.hier[hier.hi[1]+2].R = r10;
  hier.hier[hier.hi[1]+3].R = r11;

  hier.hier[hier.hi[1]+0].m = m00;
  hier.hier[hier.hi[1]+1].m = m01;
  hier.hier[hier.hi[1]+2].m = m10;
  hier.hier[hier.hi[1]+3].m = m11;

  hier.hier[hier.hi[2]+0].m = m00 + m01;
  hier.hier[hier.hi[2]+1].m = m10 + m11;

  hier.hier[hier.hi[2]+0].a = a0;
  hier.hier[hier.hi[2]+1].a = a1;
  
  hier.hier[hier.hi[2]+0].e = e0;
  hier.hier[hier.hi[2]+1].e = e1;

  hier.nobj = 2;
  hier.obj[0] = &(hier.hier[hier.hi[2]+0]);
  hier.obj[1] = &(hier.hier[hier.hi[2]+1]);
  hier.obj[2] = NULL;
  hier.obj[3] = NULL;

  /* get the units and normalize */
  calc_units(hier.obj, &units);
  fb_normalize(&hier, units);

  vinf = input_vinf / units.v;
  vcrit = sqrt(FB_CONST_G / ((m00 + m01) * (m10 + m11) / (m00 + m01 + m10 \
    + m11)) * (m00 * m01 / a0 + m10 * m11 / a1)) / units.v;

  if (input_bmax == -1) {
    bmax = (4 / (vinf / vcrit) + 3) * (a0 + a1) * xbmax;
  } else {
    bmax = input_bmax;
  }

  bmax /= units.l;
  b = sqrt(gsl_rng_uniform(rng)) * bmax;

  rtid = pow(2.0*(hier.hier[hier.hi[2]+0].m + hier.hier[hier.hi[2]+1].m) / \
          ((hier.hier[hier.hi[2]+0].m) * input.tidaltol), 1. / 3) * \
          hier.hier[hier.hi[2]+0].a * (1 + hier.hier[hier.hi[2]+0].e);
  //fprintf(stdout, "%g\n", rtid);

  // JMA 3-12-14 -- Calculate the crossing time and set tstop.
  input.tstop = 2 * rtid / vinf * input_tstop;

  fprintf(stderr, "  vinf=%.6g km/s  vcrit=%.6g km/s  bmax=%.6g AU  b=%.6g AU\n", 
    vinf*units.v/FB_CONST_KMS, vcrit*units.v/FB_CONST_KMS, 
    bmax*units.l/FB_CONST_AU, b*units.l/FB_CONST_AU);

  fprintf(stderr, "  tstop=%.6g  tcpustop=%.6g\n\n", \
    input.tstop, input.tcpustop);

  fb_init_scattering(hier.obj, vinf, b, rtid);
  fb_randscat(hier.obj[1], rng);
 
  /* randomize binary orientations and downsync */
  fb_binaryorient(&(hier.hier[hier.hi[2]+0]), rng, -1, -1, FB_CONST_PI);
  fb_downsync(&(hier.hier[hier.hi[2]+0]), t);
  fb_binaryorient(&(hier.hier[hier.hi[2]+1]), rng, -1, -1, FB_CONST_PI);
  fb_dprintf("\nProperties of the interloping binary: a, e, L[0]: %.6f %.6f %.6f\n", hier.hier[hier.hi[2]+1].a, hier.hier[hier.hi[2]+1].e, hier.hier[hier.hi[2]+1].Lhat[2]);
  fb_downsync(&(hier.hier[hier.hi[2]+1]), t);
  
  fb_dprintf("triple x-coor: %g\n", hier.obj[0]->x[0]);
  fb_dprintf("triple x-coor: %g\n", hier.hier[hier.hi[3]].x[0]);
  fb_dprintf("binary x-coor: %g\n", hier.hier[hier.hi[2]+1].x[0]);
  fb_dprintf("star coors: %.16f %.16f %.16f\n", hier.hier[hier.hi[1]].x[0], hier.hier[hier.hi[1]+1].x[0], hier.hier[hier.hi[1]+2].x[0]);
  fb_dprintf("\n");


  /* JMA 7-10-12 -- Check to make sure the angle between the two binaries
   * is still the inclination.  This can be removed later if the code
   * proves to be working as expected.
   */
  /*
  fprintf(stderr, "inner: %g, %g, %g\n", hier.hier[hier.hi[2]+0].Lhat[0], hier.hier[hier.hi[2]+0].Lhat[1], hier.hier[hier.hi[2]+0].Lhat[2]);
  fprintf(stderr, "outer: %g, %g, %g\n", hier.hier[hier.hi[3]+0].Lhat[0], hier.hier[hier.hi[3]+0].Lhat[1], hier.hier[hier.hi[3]+0].Lhat[2]);
  fprintf(stderr, "dot: %g\n", fb_dot(hier.hier[hier.hi[2]+0].Lhat, hier.hier[hier.hi[3]+0].Lhat));
  fprintf(stderr, "cos i: %g\n", cos(inc));
  assert(fb_dot(hier.hier[hier.hi[2]+0].Lhat, hier.hier[hier.hi[3]+0].Lhat) == cos(inc));
  */

  fprintf(stderr, "UNITS:\n");
  fprintf(stderr, "  v=%.6g km/s  l=%.6g AU  t=t_dyn=%.6g yr\n", \
    units.v/1.0e5, units.l/FB_CONST_AU, units.t/FB_CONST_YR);
  fprintf(stderr, "  M=%.6g M_sun  E=%.6g erg\n\n", units.m/FB_CONST_MSUN, units.E);

  /* trickle down properties (not sure if this is actually needed here, but it doesn't harm anything) */
  fb_trickle(&hier, t);

  fb_dprintf("after first trickle...\n");
  fb_dprintf("triple x-coor: %g\n", hier.obj[0]->x[0]);
  fb_dprintf("triple x-coor: %g\n", hier.hier[hier.hi[3]].x[0]);
  fb_dprintf("binary x-coor: %g\n", hier.hier[hier.hi[2]].x[0]);
  fb_dprintf("\n");

  /* store the initial energy and angular momentum*/
  Ei = fb_petot(&(hier.hier[hier.hi[1]]), hier.nstar) + fb_ketot(&(hier.hier[hier.hi[1]]), hier.nstar) +
    fb_einttot(&(hier.hier[hier.hi[1]]), hier.nstar);
  fb_angmom(&(hier.hier[hier.hi[1]]), hier.nstar, Li);
  fb_angmomint(&(hier.hier[hier.hi[1]]), hier.nstar, Lint);
  for (j=0; j<3; j++) {
    Li[j] += Lint[j];
  }

  input.Ei = Ei;
  input.Li = fb_mod(Li);

  /* integrate along */
  fb_dprintf("calling fewbody()...\n");
  
  /* call fewbody! */
  retval = fewbody(input, units, &hier, &t, rng);

  /* print information to screen */
  fprintf(stderr, "OUTCOME:\n");
  if (retval.retval == 1) {
    fprintf(stderr, "  encounter complete:  t=%.6g (%.6g yr)  %s  (%s)\n\n",
      t, t * units.t/FB_CONST_YR,
      fb_sprint_hier(hier, string1),
      fb_sprint_hier_hr(hier, string2));
  } else {
    fprintf(stderr, "  encounter NOT complete:  t=%.6g (%.6g yr)  %s  (%s)\n\n",
      t, t * units.t/FB_CONST_YR,
      fb_sprint_hier(hier, string1),
      fb_sprint_hier_hr(hier, string2));
  }

  fb_dprintf("there were %ld integration steps\n", retval.count);
  fb_dprintf("fb_classify() was called %ld times\n", retval.iclassify);
  
  fprintf(stderr, "FINAL:\n");
  fprintf(stderr, "  t_final=%.6g (%.6g yr)  t_cpu=%.6g s\n", \
    t, t*units.t/FB_CONST_YR, retval.tcpu);

  fprintf(stderr, "  L0=%.6g  DeltaL/L0=%.6g  DeltaL=%.6g  MaxDeltaL=%.6g\n", 
    fb_mod(Li), retval.DeltaLfrac, retval.DeltaL, retval.maxDeltaLfrac);
  fprintf(stderr, "  E0=%.6g  DeltaE/E0=%.6g  DeltaE=%.6g  MaxDeltaE=%.6g\n", 
    Ei, retval.DeltaEfrac, retval.DeltaE, retval.maxDeltaEfrac);
  fprintf(stderr, "  Rmin=%.6g (%.6g RSUN)  Rmin_i=%d  Rmin_j=%d\n", \
    retval.Rmin, retval.Rmin*units.l/FB_CONST_RSUN, retval.Rmin_i, retval.Rmin_j);
  fprintf(stderr, "  Nosc=%d (%s)\n", retval.Nosc, (retval.Nosc>=1?"resonance":"non-resonance"));
  
  linecount = 0;
  for (i=0; i < hier.nobj; i++) {
    if (hier.obj[i]->n == 3) {
      fprintf(stderr, "  Triple endstate: mutual inclination (degrees): %.6g\n",
        180. / FB_CONST_PI * acos(fb_dot(hier.hier[hier.hi[3]].Lhat, hier.hier[hier.hi[2]].Lhat)));
      fprintf(stderr, "                   a0=%.6g AU  e0=%.6g  a1=%.6g AU  e1=%.6g\n",
        hier.hier[hier.hi[2]].a * units.l, hier.hier[hier.hi[2]].e, hier.hier[hier.hi[3]].a * units.l, hier.hier[hier.hi[3]].e);
      linecount += 2;
    }
  }

  for (i=0; i < hier.nobj; i++) {
    if (hier.obj[i]->n == 2) {
      fprintf(stderr, "  Binary endstate: a=%.6g AU  e=%.6g\n",
        hier.hier[hier.hi[2]].a * units.l, hier.hier[hier.hi[2]].e);
      linecount += 1;
    }
  }

  while (linecount < 2) {
    fprintf(stderr, "\n");
    linecount += 1;
  }
  
  fprintf(stderr, "######\n");
  
  /* free GSL stuff */
  gsl_rng_free(rng);

  /* free our own stuff */
  fb_free_hier(hier);

  /* done! */
  return(0);
}

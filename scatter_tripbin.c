/* -*- linux-c -*- */
/* triple.c

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
#include "scatter_tripbin.h"

/* print the usage */
void print_usage(FILE *stream)
{
  fprintf(stream, "USAGE:\n");
  fprintf(stream, "  scatter_tripsingle [options...]\n");
  fprintf(stream, "\n");
  fprintf(stream, "OPTIONS:\n");
  fprintf(stream, "  -M --m10 <m10/MSUN>          : set mass of star 0 of inner binary of interloping binary [%.6g]\n", FB_M10/FB_CONST_MSUN);
  fprintf(stream, "  -Z --m11 <m11/MSUN>          : set mass of star 1 of inner binary of interloping binary [%.6g]\n", FB_M10/FB_CONST_MSUN);
  fprintf(stream, "  -m --m000 <m000/MSUN>        : set mass of star 0 of inner binary of triple [%.6g]\n", FB_M000/FB_CONST_MSUN);
  fprintf(stream, "  -n --m001 <m001/MSUN>        : set mass of star 1 of inner binary of triple [%.6g]\n", FB_M001/FB_CONST_MSUN);
  fprintf(stream, "  -o --m01 <m01/MSUN>          : set mass of outer star of triple [%.6g]\n", FB_M01/FB_CONST_MSUN);
  fprintf(stream, "  -r --r10 <r10/R_SUN>         : set merge radius of interloper [%.6g]\n", FB_R10/FB_CONST_RSUN);
  fprintf(stream, "  -W --r11 <r11/R_SUN>         : set merge radius of interloper [%.6g]\n", FB_R11/FB_CONST_RSUN);
  fprintf(stream, "  -H --r01 <r01/R_SUN>         : set merge radius of tertiary [%.6g]\n", FB_R01/FB_CONST_RSUN);
  fprintf(stream, "  -C --r000 <r000/R_SUN>       : set merge radius of star 0 of inner binary [%.6g]\n", FB_R000/FB_CONST_RSUN);
  fprintf(stream, "  -G --r001 <r001/R_SUN>       : set merge radius of star 1 of inner binary [%.6g]\n", FB_R001/FB_CONST_RSUN);
  fprintf(stream, "  -a --a00 <a00/AU>            : set inner semimajor axis of triple [%.6g]\n", FB_A00/FB_CONST_AU);
  fprintf(stream, "  -q --a0 <a0/AU>              : set outer semimajor axis of triple [%.6g]\n", FB_A0/FB_CONST_AU);
  fprintf(stream, "  -L --a1 <a1/AU>              : set semimajor axis of interloping binary [%.6g]\n", FB_A1/FB_CONST_AU);
  fprintf(stream, "  -e --e00 <e00>               : set inner eccentricity of triple [%.6g]\n", FB_E00);
  fprintf(stream, "  -F --e0 <e0>                 : set outer eccentricity of triple [%.6g]\n", FB_E0);
  fprintf(stream, "  -E --e1 <e1>                 : set eccentricity of interloping binary [%.6g]\n", FB_E1);
  fprintf(stream, "  -I --inc <i/degree>          : set inclination of the tertiary\n");
  fprintf(stream, "                                 from angular momentum vector of inner (-1 for random) [%.6g]\n", FB_INC);
  fprintf(stream, "  -p --peri <p00/degree>       : set argument of periapsis for the inner binary (-1 for random) [%.6g]\n", FB_P00);
  fprintf(stream, "  -B --peri <p0/degree>        : set argument of periapsis for the outer binary (-1 for random) [%.6g]\n", FB_P0);
  fprintf(stream, "  -Y --peri <p1/degree>        : set argument of periapsis for the interloping binary (-1 for random) [%.6g]\n", FB_P1);
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
  double m0, m00, m01, m000, m001, a00, a0;

  m0 = obj[0]->m;
  m00 = obj[0]->obj[0]->m;
  m01 = obj[0]->obj[1]->m;
  m000 = obj[0]->obj[0]->obj[0]->m;
  m001 = obj[0]->obj[0]->obj[1]->m;

  a0 = obj[0]->a;
  a00 = obj[0]->obj[0]->a;
  
  /* Unit of velocity is approximate relative orbital speed of inner binary,
     unit of length is semimajor axis of inner binary; 
     therefore, unit of time is approximately 1 inner orbital period. */
  units->v = sqrt(FB_CONST_G*(m000+m001)/a00);
  units->l = a00;
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
  double m000, m001, m01, m10, m11, r000, r001, r01, r10, r11, a00, a0, a1, e00, e0, e1, inc;
  double p00, p0, p1, inc_in, inc_out;
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
  const char *short_opts = "m:n:o:M:Z:r:W:H:C:G:Q:g:i:a:q:L:e:F:E:p:B:Y:I:b:X:v:t:D:c:A:R:N:O:z:x:y:P:Q:S:T:U:k:s:dVh";
  const struct option long_opts[] = {
    {"m000", required_argument, NULL, 'm'},
    {"m001", required_argument, NULL, 'n'},
    {"m01", required_argument, NULL, 'o'},
    {"m10", required_argument, NULL, 'M'},
    {"m11", required_argument, NULL, 'Z'},
    {"r10", required_argument, NULL, 'r'},
    {"r11", required_argument, NULL, 'W'},
    {"r01", required_argument, NULL, 'H'},
    {"r000", required_argument, NULL, 'C'},
    {"r001", required_argument, NULL, 'G'},
    {"a00", required_argument, NULL, 'a'},
    {"a0", required_argument, NULL, 'q'},
    {"a1", required_argument, NULL, 'L'},
    {"e00", required_argument, NULL, 'e'},
    {"e0", required_argument, NULL, 'F'},
    {"e1", required_argument, NULL, 'E'},
    {"p00", required_argument, NULL, 'p'},
    {"p0", required_argument, NULL, 'B'},
    {"p1", required_argument, NULL, 'Y'},
    {"inc", required_argument, NULL, 'I'},
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
    {"fexp", required_argument, NULL, 'x'},
    {"ks", required_argument, NULL, 'k'},
    {"seed", required_argument, NULL, 's'},
    {"debug", no_argument, NULL, 'd'},
    {"version", no_argument, NULL, 'V'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}
  };

  /* set parameters to default values */
  m000 = FB_M000;
  m001 = FB_M001;
  m01 = FB_M01;
  m10 = FB_M10;
  m11 = FB_M11;
  r10 = FB_R10;
  r11 = FB_R11;
  r01 = FB_R01;
  r000 = FB_R000;
  r001 = FB_R001;
  a00 = FB_A00;
  a0 = FB_A0;
  a1 = FB_A1;
  e00 = FB_E00;
  e0 = FB_E0;
  e1 = FB_E1;
  p00 = FB_P00;
  p0 = FB_P0;
  p1 = FB_P1;
  inc = FB_INC;
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
      m000 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'n':
      m001 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'o':
      m01 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'M':
      m10 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'Z':
      m11 = atof(optarg) * FB_CONST_MSUN;
      break;
    case 'r':
      r10 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'W':
      r11 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'H':
      r01 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'C':
      r000 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'G':
      r001 = atof(optarg) * FB_CONST_RSUN;
      break;
    case 'a':
      a00 = atof(optarg) * FB_CONST_AU;
      break;
    case 'q':
      a0 = atof(optarg) * FB_CONST_AU;
      break;
    case 'L':
      a1 = atof(optarg) * FB_CONST_AU;
      break;
    case 'e':
      e00 = atof(optarg);
      if (e00 >= 1.0) {
        fprintf(stderr, "e00 must be less than 1\n");
        return(1);
      }
      break;
    case 'F':
      e0 = atof(optarg);
      if (e0 >= 1.0) {
        fprintf(stderr, "e0 must be less than 1\n");
        return(1);
      }
      break;
    case 'E':
      e1 = atof(optarg);
      if (e1 >= 1.0) {
        fprintf(stderr, "e1 must be less than 1\n");
        return(1);
      }
      break;
    case 'p':
      if (atof(optarg) > 0) {
        p00 = atof(optarg) * FB_CONST_PI / 180.;
      } else {
        p00 = atof(optarg);
      }
      if (p00 > 360 || (p00 < 0.0 && p00 != -1.0)) {
        fprintf(stderr, "argument of periapsis must be between 0 and 360 (or -1 for random argument).\n");
        return(1);
      }
      break;
    case 'B':
      if (atof(optarg) > 0) {
        p0 = atof(optarg) * FB_CONST_PI / 180.;
      } else {
        p0 = atof(optarg);
      }
      if (p0 > 360 || (p0 < 0.0 && p0 != -1.0)) {
        fprintf(stderr, "argument of periapsis must be between 0 and 360 (or -1 for random argument).\n");
        return(1);
      }
      break;
    case 'Y':
      if (atof(optarg) > 0) {
        p1 = atof(optarg) * FB_CONST_PI / 180.;
      } else {
        p1 = atof(optarg);
      }
      if (p1 > 360 || (p1 < 0.0 && p1 != -1.0)) {
        fprintf(stderr, "argument of periapsis must be between 0 and 360 (or -1 for random argument).\n");
        return(1);
      }
      break;
    case 'I':
      if (atof(optarg) > 0) {
        inc = atof(optarg) * FB_CONST_PI / 180.;
      } else {
        inc = atof(optarg);
      }
      if (inc > 360. || (inc < 0.0 && inc != -1.0)) {
        fprintf(stderr, "inclination must be between 0 and 180 (or -1 for random inclination).\n");
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
  hier.nstarinit = 5;
  hier.nstar = 5;
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

  if (inc == -1.) {
    inc = acos(2 * gsl_rng_uniform(rng) - 1);
  }
  if (p00 == -1.) {
    p00 = 2 * FB_CONST_PI * gsl_rng_uniform(rng);
  }
  if (p0 == -1.) {
    p0 = 2 * FB_CONST_PI * gsl_rng_uniform(rng);
  }
  if (p1 == -1.) {
    p1 = 2 * FB_CONST_PI * gsl_rng_uniform(rng);
  }

  fb_dprintf("a00=%.6g\n", a00);

  /* print out values of paramaters */
  fprintf(stderr, "PARAMETERS:\n");
  fprintf(stderr, "  ks=%d  seed=%lu\n", input.ks, seed);
  fprintf(stderr, "  a00=%.6g AU  e00=%.6g  m000=%.6g MSUN  m001=%.6g MSUN r000=%.6g RSUN  r001=%.6g RSUN\n", \
    a00/FB_CONST_AU, e00, m000/FB_CONST_MSUN, m001/FB_CONST_MSUN, r000/FB_CONST_RSUN, r001/FB_CONST_RSUN);
  fprintf(stderr, "  a0=%.6g AU  e0=%.6g  m01=%.6g MSUN\n", \
    a0/FB_CONST_AU, e0, m01/FB_CONST_MSUN);
  fprintf(stderr, "  a1=%.6g AU  e1=%.6g  m10=%.6g MSUN  m11=%.6g MSUN\n", \
    a1/FB_CONST_AU, e1, m01/FB_CONST_MSUN, m11/FB_CONST_MSUN);
  fprintf(stderr, "  inc=%.6g  p00=%.6g  p0=%.6g  p1=%.6g\n", \
    inc * 180 / FB_CONST_PI, p00 * 180 / FB_CONST_PI, p0 * 180 / FB_CONST_PI, p1 * 180 / FB_CONST_PI);
  fprintf(stderr, "  tidaltol=%.6g  speedtol=%.6g  abs_acc=%.6g rel_acc=%.6g  ncount=%d  fexp=%.6g  outfreq=%d\n", \
    input.tidaltol, input.speedtol, input.absacc, input.relacc, input.ncount, input.fexp, input.outfreq);
  fprintf(stderr, "  PN1=%d  PN2=%d  PN25=%d  PN3=%d  PN35=%d\n", \
    input.PN1, input.PN2, input.PN25, input.PN3, input.PN35);


  /* create hierarchies */
  hier.narr[2] = 2;
  hier.narr[3] = 1;
  hier.narr[4] = 0;
  /* inner binary of triple */
  hier.hier[hier.hi[2]+0].obj[0] = &(hier.hier[hier.hi[1]+0]);
  hier.hier[hier.hi[2]+0].obj[1] = &(hier.hier[hier.hi[1]+1]);
  hier.hier[hier.hi[2]+0].t = t;
  /* outer binary of triple */
  hier.hier[hier.hi[3]+0].obj[0] = &(hier.hier[hier.hi[2]+0]);
  hier.hier[hier.hi[3]+0].obj[1] = &(hier.hier[hier.hi[1]+2]);
  hier.hier[hier.hi[3]+0].t = t;
  // JMA 12-12-13 -- Interloping binary
  hier.hier[hier.hi[2]+1].obj[0] = &(hier.hier[hier.hi[1]+3]);
  hier.hier[hier.hi[2]+1].obj[1] = &(hier.hier[hier.hi[1]+4]);
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

  hier.hier[hier.hi[1]+0].R = r000;
  hier.hier[hier.hi[1]+1].R = r001;
  hier.hier[hier.hi[1]+2].R = r01;
  hier.hier[hier.hi[1]+3].R = r10;
  hier.hier[hier.hi[1]+4].R = r11;

  // JMA 6-1-12 
  /*
  hier.hier[hier.hi[1]+0].R = r000*2*FB_CONST_G*m000/(FB_CONST_C*FB_CONST_C);
  hier.hier[hier.hi[1]+1].R = r000*2*FB_CONST_G*m001/(FB_CONST_C*FB_CONST_C);
  hier.hier[hier.hi[1]+2].R = r000*2*FB_CONST_G*m01/(FB_CONST_C*FB_CONST_C);
  */

  hier.hier[hier.hi[1]+0].m = m000;
  hier.hier[hier.hi[1]+1].m = m001;
  hier.hier[hier.hi[1]+2].m = m01;
  hier.hier[hier.hi[1]+3].m = m10;
  hier.hier[hier.hi[1]+4].m = m11;

  hier.hier[hier.hi[2]+0].m = m000 + m001;
  hier.hier[hier.hi[3]+0].m = m000 + m001 + m01;

  hier.hier[hier.hi[2]+1].m = m10 + m11;

  hier.hier[hier.hi[2]+0].a = a00;
  hier.hier[hier.hi[3]+0].a = a0;
  hier.hier[hier.hi[2]+1].a = a1;
  
  hier.hier[hier.hi[2]+0].e = e00;
  hier.hier[hier.hi[3]+0].e = e0;
  hier.hier[hier.hi[2]+1].e = e1;

  hier.nobj = 2;
  hier.obj[0] = &(hier.hier[hier.hi[3]+0]);
  hier.obj[1] = &(hier.hier[hier.hi[2]+1]);
  hier.obj[2] = NULL;
  hier.obj[3] = NULL;

  /* get the units and normalize */
  calc_units(hier.obj, &units);
  fb_normalize(&hier, units);

  vinf = input_vinf / units.v;
  vcrit = sqrt((FB_CONST_G * (m000 + m001 + m01 + m10 + m11) * (a00 * a1 * 
    (m000 + m001) * m01 + a0 * (a1 * m000 * m001 + a00 * m10 * m11))) 
    / (a0 * a00 * a1 * (m000 + m001 + m01) * (m10 + m11))) / units.v;
  if (input_bmax == -1) {
    bmax = (4 / (vinf / vcrit) + 3) * (a0 + a1) * xbmax;
  } else {
    bmax = input_bmax;
  }

  bmax /= units.l;
  b = sqrt(gsl_rng_uniform(rng)) * bmax;

  rtid = pow(2.0*(hier.hier[hier.hi[2]+1].m + hier.hier[hier.hi[3]].m) / \
          ((hier.hier[hier.hi[3]].m) * input.tidaltol), 1. / 3) * \
          hier.hier[hier.hi[3]].a * (1 + hier.hier[hier.hi[3]].e);
  //fprintf(stdout, "%g\n", rtid);

  // JMA 3-12-14 -- Calculate the crossing time and set tstop.
  input.tstop = 2 * rtid / vinf * input_tstop;

  fprintf(stderr, "  vinf=%.6g km/s  vcrit=%.6g km/s  bmax=%.6g AU  b=%.6g AU\n", 
    vinf*units.v/FB_CONST_KMS, vcrit*units.v/FB_CONST_KMS, 
    bmax*units.l/FB_CONST_AU, b*units.l/FB_CONST_AU);

  fprintf(stderr, "  tstop=%.6g  tcpustop=%.6g\n\n", \
    input.tstop, input.tcpustop);

  fb_init_scattering(hier.obj, vinf, b, rtid);
  /*
  fprintf(stdout, "%.12f %.12f %g %g %g %g %g %g\n", t, sqrt(fb_sqr(hier.hier[hier.hi[1]+3].x[0]) + \
    fb_sqr(hier.hier[hier.hi[1]+3].x[1]) + fb_sqr(hier.hier[hier.hi[1]+3].x[2])),
    hier.hier[hier.hi[1]+3].x[0], hier.hier[hier.hi[1]+3].x[1], hier.hier[hier.hi[1]+3].x[2],
    hier.hier[hier.hi[1]+3].v[0], hier.hier[hier.hi[1]+3].v[1], hier.hier[hier.hi[1]+3].v[2]);
  */
  fb_randscat(hier.obj[1], rng);
  /*
  fprintf(stdout, "%.12f %.12f %g %g %g %g %g %g\n", t, sqrt(fb_sqr(hier.hier[hier.hi[1]+3].x[0]) + \
    fb_sqr(hier.hier[hier.hi[1]+3].x[1]) + fb_sqr(hier.hier[hier.hi[1]+3].x[2])),
    hier.hier[hier.hi[1]+3].x[0], hier.hier[hier.hi[1]+3].x[1], hier.hier[hier.hi[1]+3].x[2],
    hier.hier[hier.hi[1]+3].v[0], hier.hier[hier.hi[1]+3].v[1], hier.hier[hier.hi[1]+3].v[2]);
  */
 
  /* JMA 7-10-12 -- Partition the inclination about the invariant plane. */
  inc_out = fb_incpartition(hier.obj, inc);
  inc_in = inc - inc_out;

  /* JMA 7-11-12 -- Debugging code. */
  /*
  fprintf(stderr, "inc: %g\n", inc);
  fprintf(stderr, "inc_out: %g\n", inc_out);
  fprintf(stderr, "inc_in: %g\n", inc_in);
  */

  /* randomize binary orientations and downsync */
  fb_binaryorient(&(hier.hier[hier.hi[3]+0]), rng, inc_out, p0, 0.0);
  fb_downsync(&(hier.hier[hier.hi[3]+0]), t);
  //fprintf(stdout, "%g %g %g %g %g\n", m000, m1, input.tidaltol, a0, e0);
  fb_binaryorient(&(hier.hier[hier.hi[2]+0]), rng, inc_in, p00, FB_CONST_PI);
  fb_downsync(&(hier.hier[hier.hi[2]+0]), t);
  //fprintf(stdout, "%g %g %g %g %g\n", m000, m1, input.tidaltol, a0, e0);
  fb_binaryorient(&(hier.hier[hier.hi[2]+1]), rng, -1, p1, FB_CONST_PI);
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
      fprintf(stderr, "                   a0 = %.6g  e0 = %.6g  a1 = %.6g  e1 = %.6g\n",
        hier.hier[hier.hi[2]].a, hier.hier[hier.hi[2]].e, hier.hier[hier.hi[3]].a, hier.hier[hier.hi[3]].e);
      linecount += 2;
    }
  }

  for (i=0; i < hier.nobj; i++) {
    if (hier.obj[i]->n == 2) {
      fprintf(stderr, "  Binary endstate: a = %.6g  e = %.6g\n",
        hier.hier[hier.hi[2]].a, hier.hier[hier.hi[2]].e);
      linecount += 1;
    }
  }

  if (hier.obj[0]->n == 4) {
    fprintf(stderr, "  Quadruple endstate: inc0: %.6g  inc1: %.6g\n",
      180. / FB_CONST_PI * acos(fb_dot(hier.hier[hier.hi[4]].Lhat, hier.hier[hier.hi[2]].Lhat)),
      180. / FB_CONST_PI * acos(fb_dot(hier.hier[hier.hi[4]].Lhat, hier.hier[hier.hi[2]+1].Lhat)));
    fprintf(stderr, "  aout: %.6g eout: %.6g\n", hier.hier[hier.hi[4]].a, hier.hier[hier.hi[4]].e);
    fprintf(stderr, "  a0 = %.6g  e0 = %.6g  a1 = %.6g  e1 = %.6g\n",
      hier.hier[hier.hi[2]].a, hier.hier[hier.hi[2]].e, hier.hier[hier.hi[2]+1].a, hier.hier[hier.hi[2]+1].e);
    linecount += 3;
  }

  while (linecount < 3) {
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

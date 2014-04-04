/* -*- linux-c -*- */
/* triple.h

   Copyright (C) 2006 John M. Fregeau
   
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

#define FB_TIDALTOL 1.0e-5
#define FB_SPEEDTOL 1.0e-4

#define FB_PN1 1
#define FB_PN2 1
#define FB_PN25 1
#define FB_PN3 1
#define FB_PN35 1

/* numbering convention:

          triple
             0
	    / \
          00   01
         /  \
      000    001

*/

#define FB_M000 (1.0 * FB_CONST_MSUN)
#define FB_M001 (1.0 * FB_CONST_MSUN)
#define FB_M01 (1.0 * FB_CONST_MSUN)

#define FB_R000 (0.0 * FB_CONST_RSUN)
#define FB_R001 (0.0 * FB_CONST_RSUN)
#define FB_R01 (0.0 * FB_CONST_RSUN)

#define FB_A00 (1.0 * FB_CONST_AU)
#define FB_A0 (2.0 * FB_CONST_AU)

#define FB_E00 0.0
#define FB_E0 0.0

// JMA 6-8-12
#define FB_INC -1.0 /* randomized inclinations */
#define FB_PERIARG_IN -1.0 /* randomized argument of periapsis */
#define FB_PERIARG_OUT -1.0 /* randomized argument of periapsis */

#define FB_DT 1.0 /* approximate output dt */
#define FB_TSTOP 1.0e12 /* in units of t_dyn */
#define FB_TCPUSTOP 3600.0 /* in seconds */

#define FB_ABSACC 1.0e-14 /* absolute accuracy of integrator */
#define FB_RELACC 1.0e-14 /* relative accuracy of integrator */
#define FB_NCOUNT 1 /* number of timesteps between calls to classify() */
#define FB_OUTFREQ 1000 /* number of timesteps between printing orbital information */

#define FB_KS 0

#define FB_FEXP 3.0 /* expansion factor of merger product */

#define FB_SEED 0UL
#define FB_DEBUG 0

/* effective BH radius (in units of 2M); this has to be >~6.5, since otherwise 
   we can get v/c>~0.4 and the higher order PN terms will fail */
// JMA 6-8-12 -- We can actually get weird things happening at ~100 M.
// Since we are only concerned with the merger timescale, going from 6.5 to
// 100 decreases the merger time only negligibly.
#define FB_REFF_BH 10

void print_usage(FILE *stream);
int calc_units(fb_obj_t *obj[2], fb_units_t *units);

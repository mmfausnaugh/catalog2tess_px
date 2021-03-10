/* starsig1.c
 * Alan M. Levine
 * June 12, 2017
 * Heritage: starspx1.c
 *
 * For predicted star positions in terms of pixel coordinates,
 * add sigma(column) and sigma(row).
 * Perturb the predicted pixel coordinates by Gaussian random values based
 * on the sigma values.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "seed_tw_ran.h"
#include "twister2.h"

#define NSTAR 1000

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera

double row[NSTAR], col[NSTAR], tmag[NSTAR], ra[NSTAR], dec[NSTAR];
int starno[NSTAR], ccdno[NSTAR], nst;

/******************************************************************************/
// Read RA, Dec, TESS magnitude, and color (to be inserted later) for each star
//      fprintf(stdout,"%4d %11.5f %11.5f   %11.5f %2d %9.3f %9.3f\n",
//	      i,ra[i],dec[i],tmag,iccd+1,fitpx[0],fitpx[1]);

void read_spx_stars(FILE *fpin, FILE *fpd)
{
  int i;

  i = 0;
  while(fscanf(fpin,"%d %lf %lf %lf",&starno[i],&ra[i],&dec[i],&tmag[i]) == 4) {
    fscanf(fpin,"%d %lf %lf\n",&ccdno[i],&col[i],&row[i]);
    ++i;
  }
  nst = i;
}

/******************************************************************************/

void print_sig_star(int i, FILE *fpo, FILE *fpd, double sx, double sy)
{
  fprintf(fpo,"%4d %11.5f %11.5f  %9.3f  ",starno[i],ra[i],dec[i],tmag[i]);
  fprintf(fpo,"%d %9.3f %9.3f %9.3f %9.3f\n",ccdno[i],col[i],row[i],sx,sy);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  double sigxy, sx, sy;
  int i;
  unsigned long seed;

  seed = get_tw_seed();
  fprintf(stderr,"seed = %lu \n",seed);
  seedMT(seed);

  if(argc != 2) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }

  // sigma_x or sigma_y for a tmag=8.0 star near the center of the field
  sigxy = atof(argv[1]);

  read_spx_stars(stdin,stderr);
  for(i=0;i<nst;++i) {
    // if(i > 3) break;

    // This is an oversimplified treatment - I don't have data to allow
    // a better treatment at the time of coding this function.
    sx = sigxy;
    sy = sigxy;
    row[i] += gasdev()*sx;
    col[i] += gasdev()*sy;

    print_sig_star(i,stdout,stderr,sx,sy);
  }
}

/******************************************************************************/

/* starsig2.c
 * Alan M. Levine
 * June 29, 2017
 * Heritage: starsig1.c
 *
 * For predicted star positions in terms of pixel coordinates,
 * add sigma(column) and sigma(row).
 * Perturb the predicted pixel coordinates by Gaussian random values based
 * on the sigma values.
 *
 * June 29, 2017 - Base one sigma column and row uncertainties on the expected
 *                 statistical and systematic errors.
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
/* psf_sigma = equivalent 1 sigma 1-d width of the PSF (pixels)
 * syserr = equivalent 1 sigma 1-d systematic uncertainty (pixels)
 * bkmag = TESS magnitude of background light per pixel
 * npix = number of pixels used to determine the centroid
 * return value = equivalent 1 sigma 1-d uncertainty (pixels) from the
 *   combined effects of systematic and statistical errors.
 */

double centroid_error(double psf_sigma, double syserr, double bkmag, int npix,
		      double tmag)
{
  double sigt, vart, varst, xnphot, xnbkg;

  xnbkg = 2.0e8*pow(10.0,-0.4*bkmag)*npix;

  xnphot = 2.0e8*pow(10.0,-0.4*tmag);

  varst = (sqrt(xnphot + xnbkg)/xnphot)*psf_sigma;
  varst = varst*varst;

  vart = varst + (syserr*syserr);
  sigt = sqrt(vart);
  return(sigt);
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
  double psfsig, sigsys, sx, sy, bkmag;
  int i, npix;
  unsigned long seed;

  seed = get_tw_seed();
  fprintf(stderr,"seed = %lu \n",seed);
  seedMT(seed);

  if(argc != 5) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }

  psfsig = atof(argv[1]);   // ignores variation over the field
  // sigsys = 1 sigma uncertainty due to pixelization or other systematic effects
  sigsys = atof(argv[2]);
  bkmag = atof(argv[3]);
  npix = atoi(argv[4]);

  read_spx_stars(stdin,stderr);
  for(i=0;i<nst;++i) {
    // if(i > 3) break;

    // This is still a bit oversimplified
    sx = centroid_error(psfsig,sigsys,bkmag,npix,tmag[i]);
    sy = sx;
    row[i] += gasdev()*sx;
    col[i] += gasdev()*sy;

    print_sig_star(i,stdout,stderr,sx,sy);
  }
}

/******************************************************************************/

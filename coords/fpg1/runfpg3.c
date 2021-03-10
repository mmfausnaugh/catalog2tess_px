/* runfpg3.c
 * Alan M. Levine
 * August 3, 2017
 * Heritage: runfpg2.c
 *
 * Compute geometric parameters of a TESS camera or cameras.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"
#include "fpg3.h"
#include "mrq_fpg3.h"

/******************************************************************************/
/*  Outline

Get list of cameras to work on, S/C attitude parameters, S/C space velocity
Read in star coordinates and magnitudes.
Read in starting parameter values.
Perform other initial set up tasks.

Iterate until convergence is reached:
   Compute predicted positions, and get new chi-square
   Check for convergence
   If not converged:
      If chi-sq is lowered, get new parameter values, modify alamda
      If chi-sq is not lowered, modify alamda
   If converged:
      Recompute final chi-sq
      Print results, etc.

 (note) Need to change measured star pixel 1 sigma values in 'starsig1.c'
*/
/****************************************************************************/
/* Input data needed:
 *    Command line arguments
 *       1) camera number (1-4)
 *       2) Right ascension (J2000) corresponding to direction of spacecraft z axis 
 *       3) Declination (J2000) corresponding to direction of spacecraft z axis 
 *       4) Roll angle of spacecraft orientation 
 *       5) x component of S/C velocity in km/s relative to the Solar System barycenter
 *       6) y component of S/C velocity in km/s (ditto)
 *       7) z component of S/C velocity in km/s (ditto)
 *       8) small distance in units of mm used for computing partial derivatives
 *          (0.03 mm is good)
 *       9) delta-chi-square criterion used to determine whether the program has converged
 *       10) number of iterations over which the delta-chi-square criterion is applied
 *       11) number of pixels (type = double) of a border region around each CCD in
 *           stars will not be used for fitting parameters
 *
 * Input files needed:
 *    1) (stdin) stars - ra, dec, Tmag, column, row
 *       column and row are doubles giving the measured centroid locations in a FITS
 *       file-format image
 *    2) (fpg_pars.txt) file with instrument parameter values
 *
 * Output:
 *    1) (stdout) input values,
 *                the determined parameter values,
 *                number of stars used in the final fit,
 *                chi-square of the final fit
 *    2) (rf1.serr) diagnostic information
 *
 * Notes:
 *    A sufficient no. of stars must fall on the CCDs of the selected camera
 *    for the specified spacecraft attitude in order for this code to run
 *    successfully.
 *
 *    The star centroid locations should be computed on the basis of a set of
 *    fpg parameters that plausibly differ from the nominal parameters.
 *    The column and row centroids should be perturbed by plausible random
 *    errors.
 */
/****************************************************************************/

int main(int argc, char *argv[])
{
  double dxy_mm, delchi, pxbor;
  int iicam, ncon;
  FILE *fpf;

  if(argc != 12) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  iicam = atoi(argv[1]) - 1;
  ra_sc = atof(argv[2]);
  dec_sc = atof(argv[3]);
  roll_sc = atof(argv[4]);
  fprintf(stderr,"ra,dec,roll S/C = %f %f %f\n",ra_sc,dec_sc,roll_sc);
  radecroll[0] = ra_sc;
  radecroll[1] = dec_sc;
  radecroll[2] = roll_sc;

  // read spacecraft velocity in km/s
  sc_vel[0] = atof(argv[5]);
  sc_vel[1] = atof(argv[6]);
  sc_vel[2] = atof(argv[7]);
  fprintf(stderr,"sc_vel[3] = %f %f %f\n",sc_vel[0],sc_vel[1],sc_vel[2]);

  dtor = M_PI/180.0;

  // read dxy_mm from command line
  dxy_mm = atof(argv[8]);

  delchi = atof(argv[9]);
  ncon = atoi(argv[10]);
  pxbor = atof(argv[11]);

  read_fpg_stars(stdin,stderr,pxbor);

  fpf = fopen("fpg_pars.txt","r");
  read_fpg_pars(fpf,stderr,iicam);
  fclose(fpf);
  fprintf(stderr,"B\n");fflush(stderr);
  delta_for_deriv(dxy_mm,2048);
  fprintf(stderr,"C\n");fflush(stderr);
  fill_mrq_param_arrays();
  fprintf(stderr,"D\n");fflush(stderr);
  fill_mrq_stars();
  fprintf(stderr,"E\n");fflush(stderr);
  chisq_hist_init();
  fprintf(stderr,"F\n");fflush(stderr);
  print_fpg_pars_deriv(stderr,1);
  alamda = -1.0;
  mrqmin_init();
  print_fpg_pars_deriv(stderr,1);
  fprintf(stderr,"G\n");fflush(stderr);
  // convergence criterion needs to be evaluated in while statement
  while(mrq_test_converge(delchi,ncon) != 1) {     // Levenberg-Marquardt loop
    fprintf(stderr,"H\n");fflush(stderr);
    mrqmin_core();
    fprintf(stderr,"main():chisq,ochisq,alamda = %f %f %e\n",chisq,ochisq,alamda);
    print_fpg_pars_deriv(stderr,1);   // For testing diagnostics only
  }
  mrq_finish();
  // Print out results
  fprintf(stderr,"finished\n");
  print_fpg_pars_deriv(stdout,0);
  print_fpg_pars_deriv(stderr,1);
  fpf = fopen("runfpg_stars.txt","w");
  print_fpg_all_stars(fpf);
  fclose(fpf);
}
/******************************************************************************/

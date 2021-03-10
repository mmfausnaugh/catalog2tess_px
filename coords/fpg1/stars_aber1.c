/* stars_aber1.c
 * Alan M. Levine
 * August 5, 2017
 * Heritage: fpg3.c, aber_stars1.c
 *
 * Compute apparent (aberrated) star celestial coordinates, valid only
 * for highly nonrelativistic velocities.
 *
 * This code writes both the unaberrated and aberrated coordinates.
 *
 */

// Allow for a few different input and output formats
// Aug. 3, 2017 - only one format is coded

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"

#define C_LIGHT   (3.0e5)            // km/s

double dtor;

double ra, dec, ra_cor, dec_cor, tmag;

/****************************************************************************/
// Get the part of a vector v[3] that is perpendicular to another vector s[3].
// This is needed for stellar aberration calculations

void get_perp(double s[3], double v[3], double vp[3])
{
  double sdotv, sdots;
  int i;

  sdotv = dot(s,v);
  sdots = dot(s,s);
  if (sdots > 0.0) {
    for(i=0;i<3;++i) {
      vp[i] = v[i] - (sdotv*s[i]/sdots);
    }
  }
  else {
    for(i=0;i<3;++i) {
      vp[i] = v[i];
    }
  }
}

/****************************************************************************/
// Compute apparent ra and dec taking velocity aberration into account.

// ra, dec, raap, decap should have units of radians
// v[3] is the velocity vector; it should be given in km/s

void aber_star_pos(double v[3], double ra, double dec, double *raap,
                   double *decap)
{
  double s[3], sap[3], beta[3], bp[3];
  int i;

  for(i=0;i<3;++i) {
    beta[i] = v[i]/C_LIGHT;
  }
  dsphcr(ra,dec,s);
  get_perp(s,beta,bp);
  for(i=0;i<3;++i) {
    sap[i] = s[i] + bp[i];
  }
  dnorm(sap);
  dcrsph(sap,raap,decap);
}

/******************************************************************************/
// Read RA, Dec, TESS magnitude for each star

int read_star(FILE *fpin, FILE *fpd, int iin)
{
  int nrd;

  switch(iin) {
  case 1:
    nrd = fscanf(fpin,"%lf %lf %lf",&ra,&dec,&tmag);
    if(nrd==3)
      return(1);
    break;
  case 2:
    nrd = fscanf(fpin,"%lf %lf %lf",&ra,&dec,&tmag);
    if(nrd==3)
      return(1);
    break;
  case 3:
    nrd = fscanf(fpin,"%lf %lf %lf",&ra,&dec,&tmag);
    if(nrd==3)
      return(1);
    break;
  default:
    return(0);
    break;
  }  
}

/******************************************************************************/

void write_star(FILE *fpo, int iin)
{
  switch(iin) {
  case 1:
    fprintf(fpo,"%10.6f %10.6f %10.6f %10.6f %8.4f\n",ra,dec,ra_cor,dec_cor,tmag);
    break;
  case 2:
    fprintf(fpo,"%10.6f %10.6f %10.6f %10.6f %8.4f\n",ra,dec,ra_cor,dec_cor,tmag);
    break;
  case 3:
    fprintf(fpo,"%10.6f %10.6f %10.6f %10.6f %8.4f\n",ra,dec,ra_cor,dec_cor,tmag);
    break;
  default:
    break;
  }  

}

/******************************************************************************/

void usage(char *prog)
{
  fprintf(stderr,"Usage: %s <input_format_integer_code(1-3)> <v_sc_x(km/s)>",prog);
  fprintf(stderr,"  <v_sc_y(km/s)> <v_sc_z(km/s)> < input_file  >  output_file\n");
  exit(-1);
}

/******************************************************************************/
/* input format codes:
 * 1 =>
 * 2 =>
 * 3 =>
 */

int main(int argc, char *argv[])
{
  double v[3], ra_ab, dec_ab;
  int iin;

  if(argc != 5) {
    usage(argv[0]);
  }

  dtor = M_PI/180.0;

  iin = atoi(argv[1]);
  v[0] = atof(argv[2]);
  v[1] = atof(argv[3]);
  v[2] = atof(argv[4]);

  while(read_star(stdin,stderr,iin)==1) {
    aber_star_pos(v,ra*dtor,dec*dtor,&ra_ab,&dec_ab);
    ra_cor = ra_ab/dtor;
    dec_cor = dec_ab/dtor;
    write_star(stdout,iin);
  }
}

/******************************************************************************/

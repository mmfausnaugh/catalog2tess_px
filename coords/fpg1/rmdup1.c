/* rmdup1.c
 * Alan M. Levine
 * January 31, 2018
 * Heritage: anres1.c (loose connection)
 *
 * Compute statistics of star position fits that come from runfpg4.c.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NSTAR 10000

double tmag[NSTAR], ra[NSTAR], dec[NSTAR], ra_cor[NSTAR], dec_cor[NSTAR];
double col[NSTAR], row[NSTAR], sigx[NSTAR], sigy[NSTAR];
int starno[NSTAR], ccdno[NSTAR];
int nst;
double dtor;

/******************************************************************************/
// Read RA, Dec, TESS magnitude, and color (to be inserted later) for each star

void read_fpg_stars(FILE *fpin, FILE *fpd)
{
  int i;

  i = 0;
  while(fscanf(fpin,"%d %lf %lf %lf %lf %lf",
               &starno[i],&ra[i],&dec[i],&ra_cor[i],&dec_cor[i],&tmag[i]) == 6) {
    fscanf(fpin,"%d %lf %lf %lf %lf",
           &ccdno[i],&col[i],&row[i],&sigx[i],&sigy[i]);
    ++i;
  }
  nst = i;
}

/******************************************************************************/
void print_fpg_star(int i, FILE *fpo)
{
  fprintf(fpo,"%4d %11.5f %11.5f  %11.5f %11.5f  %9.3f  ",
          starno[i],ra[i],dec[i],ra_cor[i],dec_cor[i],tmag[i]);
  fprintf(fpo,"%d %10.4f %10.4f %9.4f %9.4f\n",
          ccdno[i],col[i],row[i],sigx[i],sigy[i]);
}

/******************************************************************************/
// check RA, Dec pairs for duplicates - print out those identified.

void find_duplicate_stars(FILE *fpdup)
{
  int i, j, idup;
  double r1, d1, r2, d2, rsq;

  for(i=0;i<nst;++i) {
    idup = 0;
    r1 = ra[i];
    d1 = dec[i];
    for(j=i+1;j<nst;++j) {
      r2 = ra[j];
      d2 = dec[j];
      rsq = (r1 - r2)*(r1 - r2) + (d1 - d2)*(d1 -d2);
      if(rsq < 1.0e-8) {
	fprintf(stderr,"no.,ra,dec, no.,ra,dec= %d %f %f       %d %f %f\n",
		starno[i],r1,d1,starno[j],r2,d2);
	idup = 1;
      }
    }
    if(idup==0)
      print_fpg_star(i,fpdup);
  }
}

/******************************************************************************/

int main(int argc, char *argv[])
{

  dtor = M_PI/180.0;

  read_fpg_stars(stdin,stderr);
  fprintf(stderr,"nst = %d\n",nst);

  find_duplicate_stars(stdout);
}

/******************************************************************************/

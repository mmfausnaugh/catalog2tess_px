/* crosschk1.c
 * Alan M. Levine
 * February 25, 2018
 * Heritage: rmdup1.c
 *
 * Cross check two star lists.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NSTAR 10000

double tmag[2][NSTAR], ra[2][NSTAR], dec[2][NSTAR], ra_cor[2][NSTAR], dec_cor[2][NSTAR];
double col[2][NSTAR], row[2][NSTAR], sigx[2][NSTAR], sigy[2][NSTAR];
int starno[2][NSTAR], ccdno[2][NSTAR];
int nst[2];
double dtor;

/******************************************************************************/
// Read RA, Dec, TESS magnitude, and color (to be inserted later) for each star

// inform = format code
// iseq = 0 or 1 for first or second list of stars to be checked

void read_stars_format(int inform, int iseq, FILE *fpin, FILE *fpd)
{
  int i;

  i = 0;

  if(inform == 1) {
    while(fscanf(fpin,"%d %lf %lf %lf %lf %lf",
		 &starno[iseq][i],&ra[iseq][i],&dec[iseq][i],&ra_cor[iseq][i],&dec_cor[iseq][i],&tmag[iseq][i]) == 6) {
      fscanf(fpin,"%d %lf %lf %lf %lf",
	     &ccdno[iseq][i],&col[iseq][i],&row[iseq][i],&sigx[iseq][i],&sigy[iseq][i]);
      ++i;
    }
  }
  else if(inform == 2) {
    while(fscanf(fpin,"%lf %lf %lf",
		 &ra[iseq][i],&dec[iseq][i],&tmag[iseq][i]) == 3) {
      ++i;
    }
  }

  nst[iseq] = i;
}

/******************************************************************************/
void print_star_format(int inform, int iseq, int i, FILE *fpo)
{

  if(inform == 1) {
    fprintf(fpo,"%4d %11.5f %11.5f  %11.5f %11.5f  %9.3f  ",
	    starno[iseq][i],ra[iseq][i],dec[iseq][i],ra_cor[iseq][i],dec_cor[iseq][i],tmag[iseq][i]);
    fprintf(fpo,"%d %10.4f %10.4f %9.4f %9.4f\n",
	    ccdno[iseq][i],col[iseq][i],row[iseq][i],sigx[iseq][i],sigy[iseq][i]);
  }
  else if(inform == 2) {
    fprintf(fpo,"%11.5f %11.5f  %9.3f\n",
	    ra[iseq][i],dec[iseq][i],tmag[iseq][i]);
  }
}

/******************************************************************************/
// check RA, Dec pairs for duplicates - print out those identified.

void cross_check_stars(int inform, FILE *fpnodup, FILE *fpdup)
{
  int i, j, idup;
  double r1, d1, r2, d2, rsq;

  for(i=0;i<nst[0];++i) {
    idup = 0;
    r1 = ra[0][i];
    d1 = dec[0][i];
    for(j=0;j<nst[1];++j) {
      r2 = ra[1][j];
      d2 = dec[1][j];
      rsq = (r1 - r2)*(r1 - r2) + (d1 - d2)*(d1 -d2);
      if(rsq < 1.0e-8) {
	idup = 1;
      }
    }
    if(idup==0)
      print_star_format(inform,0,i,fpnodup);
    else
      print_star_format(inform,0,i,fpdup);
  }
}

/******************************************************************************/
// check RA, Dec pairs for duplicates - print out those identified.

void find_duplicate_stars(int inform, int iseq, FILE *fpdup)
{
  int i, j, idup;
  double r1, d1, r2, d2, rsq;

  for(i=0;i<nst[iseq];++i) {
    idup = 0;
    r1 = ra[iseq][i];
    d1 = dec[iseq][i];
    for(j=i+1;j<nst[iseq];++j) {
      r2 = ra[iseq][j];
      d2 = dec[iseq][j];
      rsq = (r1 - r2)*(r1 - r2) + (d1 - d2)*(d1 -d2);
      if(rsq < 1.0e-8) {
	fprintf(stderr,"no.,ra,dec, no.,ra,dec= %d %f %f       %d %f %f\n",
		starno[iseq][i],r1,d1,starno[iseq][j],r2,d2);
	idup = 1;
      }
    }
    if(idup==0)
      print_star_format(inform,iseq,i,fpdup);
  }
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  int inf1, inf2;
  FILE *fpin1, *fpin2, *fpnod, *fpdup;

  dtor = M_PI/180.0;

  if(argc != 3) {
    fprintf(stderr,"Wrong number of argument.  Try again.\n");
    exit(-1);
  }
  inf1 = atoi(argv[1]);
  inf2 = atoi(argv[2]);

  fpin1 = fopen("stars1.txt","r");
  fpin2 = fopen("stars2.txt","r");
  fpnod = fopen("stars_nodup.txt","w");
  fpdup = fopen("stars_dup.txt","w");

  read_stars_format(inf1,0,fpin1,stderr);
  fprintf(stderr,"nst[0] = %d\n",nst[0]);
  read_stars_format(inf2,1,fpin2,stderr);
  fprintf(stderr,"nst[1] = %d\n",nst[1]);
  cross_check_stars(inf1,fpnod,fpdup);
}

/******************************************************************************/

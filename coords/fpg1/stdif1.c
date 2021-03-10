/* stdif1.c
 * Alan M. Levine
 * June 24, 2017
 * Heritage: fpg_a1.c
 *
 * Compare predicted star positions in the focal plane and pixel coordinates.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NSTAR 1000

double ra[2][NSTAR], dec[2][NSTAR], lngd[2][NSTAR], colatd[2][NSTAR];
double xfp[2][NSTAR], yfp[2][NSTAR], ccdx[2][NSTAR], ccdy[2][NSTAR];
double fitx[2][NSTAR], fity[2][NSTAR];
int iccd[2][NSTAR];
int nst[2];

double dtor;

/******************************************************************************/

void read_long_stars(FILE *fpin, FILE *fpd, int ifl)
{
  int i, j;

  i = 0;
  while(fscanf(fpin,"%d %lf %lf",&j,&ra[ifl][i],&dec[ifl][i]) == 3) {
    fscanf(fpin,"%lf %lf %lf %lf %lf %lf %d %lf %lf",
	   &lngd[ifl][i],&colatd[ifl][i],&xfp[ifl][i],&yfp[ifl][i],&ccdx[ifl][i],
	   &ccdy[ifl][i],&iccd[ifl][i],&fitx[ifl][i],&fity[ifl][i]);
    ++i;
    if(i == NSTAR) {
      fprintf(stderr,"i=%d; breaking read loop ...\n",i);
      break;
    }
  }
  nst[ifl] = i;

}

/******************************************************************************/

int main(int argc, char *argv[])
{
  char *filea, *fileb;
  double dx, dy, dcx, dcy, dfx, dfy, dissq, dis, sigd, vard, dchi, chit;
  int i;
  FILE *fpf, *fpa, *fpb, *fpc;

  dtor = M_PI/180.0;

  if(argc != 4) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  filea = argv[1];
  fileb = argv[2];
  sigd = atof(argv[3]);
  vard = sigd*sigd;
  fprintf(stderr,"%s %s %f %f\n",filea,fileb,sigd,vard);fflush(stderr);

  fpf = fopen(filea,"r");
  // fprintf(stderr,"A\n");fflush(stderr);
  read_long_stars(fpf,stderr,0);
  fclose(fpf);
  fpf = fopen(fileb,"r");
  read_long_stars(fpf,stderr,1);
  fclose(fpf);
  if(nst[0] != nst[1]) {
    fprintf(stderr,"nst[0], [1] = %d %d are unequal - exiting\n",nst[0],nst[1]);
    exit(-1);
  }

  fpa = fopen("std_xy.txt","w");
  fpb = fopen("std_ccd.txt","w");
  fpc = fopen("std_fits.txt","w");

  chit = 0.0;
  for(i=0;i<nst[0];++i) {
    // if(i > 3) break;
    dx = xfp[1][i] - xfp[0][i];
    dy = yfp[1][i] - yfp[0][i];
    dissq = dx*dx + dy*dy;
    dis = sqrt(dissq);
    dchi = dissq/vard;
    chit += dchi;

    fprintf(fpa,"%4d %11.5f %11.5f   %11.5f %11.5f  %9.5f %9.5f %9.3f %9.3f\n",
	    i,xfp[0][i],yfp[0][i],xfp[1][i],yfp[1][i],dx,dy,dis,dchi);

    dx = ccdx[1][i] - ccdx[0][i];
    dy = ccdy[1][i] - ccdy[0][i];
    dissq = dx*dx + dy*dy;
    dis = sqrt(dissq);
    fprintf(fpb,"%4d %11.5f %11.5f   %11.5f %11.5f  %9.5f %9.5f %9.3f %1d %1d\n",
	    i,ccdx[0][i],ccdy[0][i],ccdx[1][i],ccdy[1][i],dx,dy,dis,iccd[0][i],iccd[1][i]);

    dx = fitx[1][i] - fitx[0][i];
    dy = fity[1][i] - fity[0][i];
    dissq = dx*dx + dy*dy;
    dis = sqrt(dissq);
    // dchi = dissq/vard;
    // chit += dchi;
    fprintf(fpc,"%4d %11.5f %11.5f   %11.5f %11.5f  %9.5f %9.5f %9.3f %9.3f\n",
	    i,fitx[0][i],fity[0][i],fitx[1][i],fity[1][i],dx,dy,dis,0.0);
  }
  fprintf(stdout,"Total chi square = %f\n",chit);
}
/******************************************************************************/

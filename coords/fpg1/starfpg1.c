/* starfpg1.c
 * Alan M. Levine
 * June 11, 2018
 * Heritage: anres3.c
 *
 * Manipulate star_fpg files such as those made by starspx6.c.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "vec.h"
//#include "mat_ra3.h"

#define NSTAR 10000

#define CCDWD_T 2048
#define CCDHT_T 2058
#define CCDWD_A 2048
#define CCDHT_A 2048
#define ROWA 44
#define ROWB 44
#define COLDK_T 20
#define COLDK_A 30

#define NCCD 4  // no. of CCDs per camera

double tmag[NSTAR], ra[NSTAR], dec[NSTAR], ra_cor[NSTAR], dec_cor[NSTAR];
double colms[NSTAR], rowms[NSTAR], sigx[NSTAR], sigy[NSTAR];
double colmod[NSTAR], rowmod[NSTAR], cdif[NSTAR], rdif[NSTAR];
int starno[NSTAR], ccdno[NSTAR];
int nst;

double dtor;

/******************************************************************************/
// Read fpg residuals file

void read_fpg_input_star(FILE *fpin, FILE *fpd)
{
  int i, nrd;

  i = 0;
  while(i<NSTAR) {
    nrd = fscanf(fpin,"%d %lf %lf %lf %lf %lf",
		 &starno[i],&ra[i],&dec[i],&ra_cor[i],&dec_cor[i],&tmag[i]);
    nrd += fscanf(fpin,"%d %lf %lf %lf %lf",
	   &ccdno[i],&colms[i],&rowms[i],&sigx[i],&sigy[i]);
    if(nrd < 11)
      break;
    ++i;
  }
  nst = i;

}

/******************************************************************************/
void print_fpg_input_star(int i, FILE *fpo)
{
  fprintf(fpo,"%4d %11.5f %11.5f  %11.5f %11.5f  %9.3f  ",
          starno[i],ra[i],dec[i],ra_cor[i],dec_cor[i],tmag[i]);
  fprintf(fpo,"%d %9.3f %9.3f %9.3f %9.3f\n",
          ccdno[i],colms[i],rowms[i],sigx[i],sigy[i]);
}

/******************************************************************************/
void print_fpg_rvsd_star(int i, FILE *fpo)
{
  fprintf(fpo,"%4d %11.5f %11.5f  %11.5f %11.5f  %9.3f  ",
          starno[i],ra[i],dec[i],ra_cor[i],dec_cor[i],tmag[i]);
  fprintf(fpo,"%d %9.3f %9.3f %9.3f %9.3f",
          ccdno[i],colmod[i],rowmod[i],sigx[i],sigy[i]);
  fprintf(fpo,"\n");
}

/******************************************************************************/
void print_fpg_star_long(int i, FILE *fpo)
{
  fprintf(fpo,"%4d %11.5f %11.5f  %11.5f %11.5f  %9.3f  ",
          starno[i],ra[i],dec[i],ra_cor[i],dec_cor[i],tmag[i]);
  fprintf(fpo,"%d %9.3f %9.3f %9.3f %9.3f",
          ccdno[i],colms[i],rowms[i],sigx[i],sigy[i]);
  fprintf(fpo,"   %9.3f %9.3f %9.3f %9.3f\n",
          colmod[i],rowmod[i],cdif[i],rdif[i]);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  int i;
  double pxcor[2];

  dtor = M_PI/180.0;

  if(argc != 3) {
    fprintf(stderr,"argc = %d is not 3, exiting...\n",argc);
  }
  pxcor[0] = atof(argv[1]);
  pxcor[1] = atof(argv[2]);

  read_fpg_input_star(stdin,stderr);
  fprintf(stderr,"nst = %d\n",nst);

  for(i=0;i<nst;++i) {
    colmod[i] = colms[i] + pxcor[0];
    rowmod[i] = rowms[i] + pxcor[1];
  }

  // Write lines for stars with well-fit positions
  for(i=0;i<nst;++i) {
    print_fpg_rvsd_star(i,stdout);
  }

}

/******************************************************************************/

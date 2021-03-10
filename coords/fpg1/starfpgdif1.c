/* starfpgdif1.c
 * Alan M. Levine
 * June 13, 2018
 * Heritage: starfpg1.c
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
int starno[NSTAR], ccdno[NSTAR];
double tmagb[NSTAR], rab[NSTAR], decb[NSTAR], ra_corb[NSTAR], dec_corb[NSTAR];
double colmsb[NSTAR], rowmsb[NSTAR], sigxb[NSTAR], sigyb[NSTAR];
int starnob[NSTAR], ccdnob[NSTAR];
double colmod[NSTAR], rowmod[NSTAR], cdif[NSTAR], rdif[NSTAR];
int nst, nstb;

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
// Read fpg residuals file

void read_fpg_input_star_b(FILE *fpin, FILE *fpd)
{
  int i, nrd;

  i = 0;
  while(i<NSTAR) {
    nrd = fscanf(fpin,"%d %lf %lf %lf %lf %lf",
		 &starnob[i],&rab[i],&decb[i],&ra_corb[i],&dec_corb[i],&tmagb[i]);
    nrd += fscanf(fpin,"%d %lf %lf %lf %lf",
	   &ccdnob[i],&colmsb[i],&rowmsb[i],&sigxb[i],&sigyb[i]);
    if(nrd < 11)
      break;
    ++i;
  }
  nstb = i;
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
  char *infile1, *infile2;
  int i;
  double pxcor[2];
  FILE *fp;

  dtor = M_PI/180.0;

  if(argc != 3) {
    fprintf(stderr,"argc = %d is not 3, exiting...\n",argc);
  }
  infile1 = argv[1];
  infile2 = argv[2];

  fp = fopen(infile1,"r");
  read_fpg_input_star(fp,stderr);
  fprintf(stderr,"nst = %d\n",nst);
  fclose(fp);

  fp = fopen(infile2,"r");
  read_fpg_input_star_b(fp,stderr);
  fprintf(stderr,"nstb = %d\n",nstb);
  fclose(fp);

  if(nstb != nst) {
    fprintf(stderr,"nst,nstb = %d %d are not equal.\n",nst,nstb);
    exit(-1);
  }

  for(i=0;i<nst;++i) {
    colmod[i] = colmsb[i];
    rowmod[i] = rowmsb[i];
    cdif[i] = colmsb[i] - colms[i];
    rdif[i] = rowmsb[i] - rowms[i];
  }

  // Write lines for stars with well-fit positions
  for(i=0;i<nst;++i) {
    print_fpg_star_long(i,stdout);
  }

}

/******************************************************************************/

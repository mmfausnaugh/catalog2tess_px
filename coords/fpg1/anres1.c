/* anres1.c
 * Alan M. Levine
 * January 30, 2018
 * Heritage: chkc1.c (loose connection)
 *
 * Compute statistics of star position fits that come from runfpg4.c.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "vec.h"
//#include "mat_ra3.h"

#define PIXEL   0.015    // pixel side in mm
#define MAXPIX 5000.0

#define NSTAR 10000
#define MAXNQ 16

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
int kcmin[NCCD], kcmax[NCCD], krmin[NCCD], krmax[NCCD], kcb[2], krb[2];
double dtor;

/******************************************************************************/
// set pixel limits of the CCDs in a FITS full frame image.
void set_ccd_lims()
{
  kcmin[0] = CCDWD_T + ROWA + (2*ROWB);
  kcmax[0] = kcmin[0] + CCDWD_T;
  krmin[0] = CCDHT_T + (2*COLDK_T);
  krmax[0] = krmin[0] + CCDHT_T;

  kcmin[1] = ROWA;
  kcmax[1] = kcmin[1] + CCDWD_T;
  krmin[1] = CCDHT_T + (2*COLDK_T);
  krmax[1] = krmin[1] + CCDHT_T;

  kcmin[2] = ROWA;
  kcmax[2] = kcmin[2] + CCDWD_T;
  krmin[2] = 0;
  krmax[2] = krmin[2] + CCDHT_T;

  kcmin[3] = CCDWD_T + ROWA + (2*ROWB);
  kcmax[3] = kcmin[3] + CCDWD_T;
  krmin[3] = 0;
  krmax[3] = krmin[3] + CCDHT_T;

  kcb[0] = kcmin[1];
  kcb[1] = kcmax[0];
  krb[0] = krmin[2];
  krb[1] = krmax[1];
}

/******************************************************************************/
// Read fpg residuals file

void read_fpg_results_file(FILE *fpin, FILE *fpd)
{
  int i, nrd;

  i = 0;
  while(i<NSTAR) {
    nrd = fscanf(fpin,"%d %lf %lf %lf %lf %lf",
		 &starno[i],&ra[i],&dec[i],&ra_cor[i],&dec_cor[i],&tmag[i]);
    nrd += fscanf(fpin,"%d %lf %lf %lf %lf",
	   &ccdno[i],&colms[i],&rowms[i],&sigx[i],&sigy[i]);
    nrd += fscanf(fpin,"%lf %lf %lf %lf",
          &colmod[i],&rowmod[i],&cdif[i],&rdif[i]);
    if(nrd < 15)
      break;
    ++i;
  }
  nst = i;

}

/******************************************************************************/
void print_fpg_star(int i, FILE *fpo)
{
  fprintf(fpo,"%4d %11.5f %11.5f  %11.5f %11.5f  %9.3f  ",
          starno[i],ra[i],dec[i],ra_cor[i],dec_cor[i],tmag[i]);
  fprintf(fpo,"%d %9.3f %9.3f %9.3f %9.3f",
          ccdno[i],colms[i],rowms[i],sigx[i],sigy[i]);
  fprintf(fpo,"   %9.3f %9.3f %9.3f %9.3f\n",
          colmod[i],rowmod[i],cdif[i],rdif[i]);
}

/******************************************************************************/
// Divide kcb[] and krb[] area into nq x nq boxes

int get_box_section(double col, double row, int nq, int isec[2])
{
  double cdf, rdf;
  int jc, jcd, jr, jrd, iret;

  cdf = (kcb[1] - kcb[0])/( (double) nq);
  jc = (col - kcb[0])/cdf;
  rdf = (krb[1] - krb[0])/( (double) nq);
  jr = (row - krb[0])/rdf;
  isec[0] = jc;
  isec[1] = jr;
  if( (jc >= 0) && (jc < nq) && (jr >= 0) && (jr < nq) )
    return(1);
  else
    return(0);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  int i, j, n, iq, jq, nq, isec[2];
  double avgc, avgr, sc, sr, crms, rrms, adst, ad0, delad, rtdel, admax;
  FILE *fpd;

  dtor = M_PI/180.0;
  set_ccd_lims();

  if(argc != 2) {
    fprintf(stderr,"argc = %d is not 2, exiting...\n",argc);
  }
  nq = atoi(argv[1]);

  read_fpg_results_file(stdin,stderr);
  fprintf(stderr,"nst = %d\n",nst);
  j = 10;
  if(nst <10)
    j = nst;
  for(i=0;i<j;++i) {
    print_fpg_star(i,stderr);
  }

  // compute average and rms difference in each coordinate
  avgc = 0.0;
  avgr = 0.0;
  n = 0;
  sc = 0.0;
  sr = 0.0;
  adst = 0.0;

  for(i=0;i<nst;++i) {
    if( (colmod[i] > MAXPIX) || (rowmod[i] > MAXPIX) ) {
      continue;
    }
    ++n;
    avgc += cdif[i];
    avgr += rdif[i];
    sc += cdif[i]*cdif[i];
    sr += rdif[i]*rdif[i];
    adst += cdif[i]*cdif[i] + rdif[i]*rdif[i];
  }
  avgc /= n;
  avgr /= n;
  crms = (sc/n) - (avgc*avgc);
  if(crms >= 0.0)
    crms = sqrt(crms);
  rrms = (sr/n) - (avgr*avgr);
  if(rrms >= 0.0)
    rrms = sqrt(rrms);
  adst = sqrt(adst/n);
  printf("n,avgc,avgr,crms,rrms,adst= %d %f %f %f %f    %f\n",
	 n,avgc,avgr,crms,rrms,adst);
  ad0 = adst;

  // compute average and rms difference in each coordinate for each sector
  // nq x nq sectors in the focal plane
  for(iq=0;iq<nq;++iq) {
    for(jq=0;jq<nq;++jq) {
      avgc = 0.0;
      avgr = 0.0;
      n = 0;
      sc = 0.0;
      sr = 0.0;
      adst = 0.0;
      admax = 0.0;

      for(i=0;i<nst;++i) {
	if( (colmod[i] > MAXPIX) || (rowmod[i] > MAXPIX) ) {
	  continue;
	}
	if(get_box_section(colms[i],rowms[i],nq,isec) != 1)
	  continue;
	if( (isec[0] != iq) || (isec[1] != jq) )
	  continue;
	++n;
	avgc += cdif[i];
	avgr += rdif[i];
	sc += cdif[i]*cdif[i];
	sr += rdif[i]*rdif[i];
	delad = cdif[i]*cdif[i] + rdif[i]*rdif[i];
	rtdel = sqrt(delad);
	if(rtdel > admax)
	  admax = rtdel;
	if (rtdel >= (4.0*ad0)) {
	  fprintf(stdout,"starno,rtdel,rtde//ad0= %d %f %f\n",
		  starno[i],rtdel,(rtdel/ad0));
	}
	adst += delad;
      }
      avgc /= n;
      avgr /= n;
      crms = (sc/n) - (avgc*avgc);
      if(crms >= 0.0)
	crms = sqrt(crms);
      rrms = (sr/n) - (avgr*avgr);
      if(rrms >= 0.0)
	rrms = sqrt(rrms);
      adst = sqrt(adst/n);
      printf("iq,jq,n,avgc,avgr,crms,rrms,adst,admax= %d %d    %d %f %f %f %f    %f %f\n",
	     iq,jq,n,avgc,avgr,crms,rrms,adst,admax);
    }
  }
}

/******************************************************************************/

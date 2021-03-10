/* stgrid1.c
 * Alan M. Levine
 * June 23, 2017
 * Heritage: fpg_a1.c
 *
 * Compute (fake) star RAs and Decs based on a grid in x and y field angles.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"

#define PIXEL   0.015    // pixel side in mm
#define FLMM    146.986  // focal length in mm

#define CCDWD_T 2048
#define CCDHT_T 2058
#define CCDWD_A 2048
#define CCDHT_A 2048
#define ROWA 44
#define ROWB 44
#define COLDK_T 20
#define COLDK_A 30

#define FA_MAX 12.0
#define DEL_FA 1.0

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];

double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3], rmat5[3][3];
double ra_sc, dec_sc, roll_sc;

double dtor;
void prmat();

/******************************************************************************/
// All input angles are in degrees.

void read_fpg_pars(FILE *fpin, FILE *fpo, int icam)
{
  char pardesc[100];
  int iccd, j;

  for(j=0;j<3;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&keulcam[icam][j],&eulcam[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,keulcam[icam][j],eulcam[icam][j]);
  }
  for(j=0;j<6;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&koptcon[icam][j],&optcon[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,koptcon[icam][j],optcon[icam][j]);
  }
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdxy0[icam][iccd][j],&ccdxy0[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
    }
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kpixsz[icam][iccd][j],&pixsz[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
    }
    fscanf(fpin,"%s %d %lf",pardesc,&kccdang[icam][iccd],&ccdang[icam][iccd]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,kccdang[icam][iccd],ccdang[icam][iccd]);
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdtilt[icam][iccd][j],&ccdtilt[icam][iccd][j]);
      fprintf(fpo,"%s       %d   %f\n",pardesc,kccdtilt[icam][iccd][j],ccdtilt[icam][iccd][j]);
    }
  }
}

/******************************************************************************/
// Multiply a 3-vector in celestial coordinates by rmat1 to get a 3-vector
// in S/C coordinates.

void sky_to_sc_mat(double radecroll[3])
{
  double xeul[3];
  int i;

  xeul[0] = dtor*radecroll[0];
  xeul[1] = (M_PI/2.0) - (dtor*radecroll[1]);
  xeul[2] = dtor*radecroll[2];
  xeul[2] += M_PI;

  eulerm323(xeul,rmat1);
}

/******************************************************************************/
// Multiply a 3-vector in S/C coordinates by rmat2 to get a 3-vector
// in camera no. icam coordinates.

void sc_to_cam_mat(int icam, double euler[3])
{
  double xeul[3], angle, rm[3][3];
  int i;

  for(i=0;i<3;++i) {
    xeul[i] = dtor*euler[i];
  }

  eulerm323(xeul,rmat2);
  prmat(rmat2,stderr);
}

/******************************************************************************/
// Print a 3x3 matrix

void prmat(double rm[3][3], FILE *fp)
{
  int i, j;

  fprintf(fp,"\n");
  for(i=0;i<3;++i) {
    fprintf(fp,"%d",i);
    for(j=0;j<3;++j) {
      fprintf(fp,"  %f",rm[i][j]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
}

/******************************************************************************/
/* return 1 if star is in the field of view of camera no. icam
 * return 0 otherwise
 */

int star_in_fov(double lngdeg, double latdeg)
{
  int ifov;
  double latr, lngr, v[3];

  lngr = dtor*lngdeg;
  latr = dtor*latdeg;

  ifov = 0;
  if(latdeg > 73.0) {
    dsphcr(lngr,latr,v);
    dnorm(v);
    if( (fabs(atan(v[0]/v[2])) <= (12.0*dtor)) &&
	(fabs(atan(v[1]/v[2])) <= (12.0*dtor)) ) {
      ifov = 1;
    }
  }

  return(ifov);
}

/******************************************************************************/
/*  TSIG/AML optics model
 *
 * lng_deg, lat_deg are the spherical coordinates of the direction to the
 * star (degrees)
 *
 * xyfp[2] contains the x and y focal plane coordinates in mm
 */

void optics_fp(double lng_deg, double lat_deg, double xyfp[2], int icam)
{
  double thetar, tanth, cphi, sphi, ttsq, rfp0, rfp, pw;
  int i;

  thetar = (M_PI/2.0) - (lat_deg*dtor);
  // fprintf(stderr,"           %e %e    %24.14e\n",lng_deg,lat_deg,thetar);
  tanth = tan(thetar);
  cphi = cos(dtor*lng_deg);
  sphi = sin(dtor*lng_deg);
  /* ttsq = tanth*tanth;
     rfp = FLMM*tanth*(1.00000140 + (0.28174612*ttsq) + (-0.59667259*ttsq*ttsq)
     + (9.17151267*ttsq*ttsq*ttsq) + (-4.36928235*ttsq*ttsq*ttsq*ttsq));
  */
  rfp0 = optcon[icam][0]*tanth;
  rfp = 0.0;
  // fprintf(stderr,"con,tanth,rfp0 = %e %e %e\n",optcon[icam][0],tanth,rfp0);
  for(i=1;i<NOPTCON;++i) {
    pw = pow(tanth,2.0*(i-1));
    rfp += optcon[icam][i]*pw;
    // fprintf(stderr,"i,con,pw,rfp = %d %e %e %e\n",i,optcon[icam][i],pw,rfp);
  }
  xyfp[0] = -cphi*rfp0*rfp;
  xyfp[1] = -sphi*rfp0*rfp;
}

/******************************************************************************/

void get_ra_dec_roll(double rm[3][3], double angs[3])
{
  double eul323[3];

  mateuler323(rm,eul323);  // angles in eul323 will be in radians
  angs[0] = eul323[0]/dtor;
  if (angs[0] < 0.0)
    angs[0] += 360.0;
  angs[1] = ((M_PI/2.0) - eul323[1])/dtor;
  angs[2] = eul323[2]/dtor;
}

/******************************************************************************/
// xfa, yfa field angles should be in radians.

void fieldangles_to_vcam(double xfa, double yfa, double vcam[3])
{
  vcam[0] = tan(xfa);
  vcam[1] = tan(yfa);
  vcam[2] = 1.0;
  dnorm(vcam);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  double radecroll[3], eulc[3], vstar[3], vcam[3], lng, lat, lngd, latd;
  double radg, decdg, ras, decs, tmag;
  double rdr[3];
  double xfa, yfa;
  int icam, iccd, i, jccd, kccd, j;

  dtor = M_PI/180.0;

  icam = 0;
  read_fpg_pars(stdin,stderr,icam);

  if(argc != 5) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  icam = atoi(argv[1]) - 1;
  ra_sc = atof(argv[2]);
  radecroll[0] = ra_sc;
  dec_sc = atof(argv[3]);
  radecroll[1] = dec_sc;
  roll_sc = atof(argv[4]);
  radecroll[2] = roll_sc;
  fprintf(stderr,"ra,dec,roll S/C = %f %f %f\n",ra_sc,dec_sc,roll_sc);
  tmag = 8.0;

  sky_to_sc_mat(radecroll);
  prmat(rmat1,stderr);
  for(j=0;j<3;++j) {
    eulc[j] = eulcam[icam][j];
  }
  sc_to_cam_mat(icam,eulc);
  matmat(rmat2,rmat1,rmat4);
  prmat(rmat4,stderr);
  get_ra_dec_roll(rmat4,rdr);
  fprintf(stderr,"cam ra,dec,roll = %f %f %f\n",
	  rdr[0],rdr[1],rdr[2]);
  trans(rmat4,rmat5);

  for(xfa=-FA_MAX;xfa<=FA_MAX;xfa+=DEL_FA) {
    for(yfa=-FA_MAX;yfa<=FA_MAX;yfa+=DEL_FA) {
      fieldangles_to_vcam(xfa*dtor,yfa*dtor,vcam);
      matvec(rmat5,vcam,vstar);
      dcrsph(vstar,&ras,&decs);
      radg = ras/dtor;
      decdg = decs/dtor;
      fprintf(stdout,"%11.6f   %11.6f  %6.2f\n",radg,decdg,tmag);
    }
  }
}

/******************************************************************************/

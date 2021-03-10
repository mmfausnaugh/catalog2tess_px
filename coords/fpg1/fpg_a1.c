/* fpg_a1.c
 * Alan M. Levine
 * June 1, 2017
 * Heritage: guid1.c
 *
 * Compute predicted star positions in the focal plane and pixel coordinates.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"

#define PIXEL   0.015    // pixel side in mm
#define FLMM    146.986  // focal length in mm

#define NSTAR 1000

#define CCDWD_T 2048
#define CCDHT_T 2058
#define CCDWD_A 2048
#define CCDHT_A 2048
#define ROWA 44
#define ROWB 44
#define COLDK_T 20
#define COLDK_A 30

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];

double row[NSTAR], col[NSTAR], tmag[NSTAR], ra[NSTAR], dec[NSTAR];
int nst;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3];
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
    if( (fabs(atan(v[0]/v[2])) <= (13.0*dtor)) &&
	(fabs(atan(v[1]/v[2])) <= (13.0*dtor)) ) {
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
// Read RA, Dec, TESS magnitude, and color (to be inserted later) for each star

void read_stars(FILE *fpin, FILE *fpd)
{
  int i;

  i = 0;
  while(fscanf(fpin,"%lf %lf %lf",&ra[i],&dec[i],&tmag[i]) == 3) {
    // fprintf(fpd,"%f %f %f\n",ra[i],dec[i],tmag[i]);
    ++i;
  }
  nst = i;

}

/******************************************************************************/
// Find new coordinates after rotating coordinate system by angle_deg

void xyrotate(double angle_deg, double xyin[2], double xyout[2])
{
  double ca, sa;

  ca = cos(dtor*angle_deg);
  sa = sin(dtor*angle_deg);
  xyout[0] = (ca*xyin[0]) + (sa*xyin[1]);
  xyout[1] = (-sa*xyin[0]) + (ca*xyin[1]);
}

/******************************************************************************/
/* Convert focal plane coordinates in mm to pseudo-equivalent TESS
 * camera image CCD pixel numbers and FITS file pixel numbers.
 */

// CCD tilt is ignored
// No checking is done of whether the star is on a CCD active area

// Camera number = icam + 1
// CCD number = iccd + 1

int mm_to_pix(int icam, double xy[2], double ccdpx[2], double fitpx[2])
{
  int iccd;
  double xya[2], xyb[2], xyccd[2];

  // CHECK SIGNS IN TESTING
  xya[0] = xy[0];
  xya[1] = xy[1];

  if(xya[0] >= 0.0) {
    if(xya[1] >= 0.0) {
      iccd = 0;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = xyccd[0]/pixsz[icam][iccd][0];
      ccdpx[1] = xyccd[1]/pixsz[icam][iccd][1];
      fitpx[0] = (CCDWD_T - ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB;
      fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T;
    }
    else {
      iccd = 3;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = xyccd[0]/pixsz[icam][iccd][0];
      ccdpx[1] = xyccd[1]/pixsz[icam][iccd][1];
      fitpx[0] = (ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB;
      fitpx[1] = (ccdpx[1]);
    }
  }
  else {
    if(xya[1] >= 0.0) {
      iccd = 1;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = xyccd[0]/pixsz[icam][iccd][0];
      ccdpx[1] = xyccd[1]/pixsz[icam][iccd][1];
      fitpx[0] = (CCDWD_T - ccdpx[0]) + ROWA;
      fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T;
    }
    else {
      iccd = 2;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = xyccd[0]/pixsz[icam][iccd][0];
      ccdpx[1] = xyccd[1]/pixsz[icam][iccd][1];
      fitpx[0] = (ccdpx[0]) + ROWA;
      fitpx[1] = (ccdpx[1]);
    }
  }
  return(iccd);
}

/******************************************************************************/
/* Convert pseudo-equivalent TESS camera image FITS file pixel numbers to
 * focal plane coordinates in pixel units.
 */

// Cam num = icam + 1
// CCD num = iccd + 1

void from_fits_pix_t(int icam, double xy[2], double fitp[2])
{
  int iccd;
  double xyp[2], fpxlo[4], fpxhi[4], fpylo[4], fpyhi[4];

  fpxlo[0] = CCDWD_T + 2*ROWA + ROWB;
  fpxlo[1] = ROWA;
  fpxlo[2] = ROWA;
  fpxlo[3] = CCDWD_T + 2*ROWA + ROWB;
  fpxhi[0] = fpxlo[0] + CCDWD_T + 1;
  fpxhi[1] = fpxlo[1] + CCDWD_T + 1;
  fpxhi[2] = fpxlo[2] + CCDWD_T + 1;
  fpxhi[3] = fpxlo[3] + CCDWD_T + 1;
  fpylo[0] = CCDHT_T + 2*COLDK_T;
  fpylo[1] = CCDHT_T + 2*COLDK_T;
  fpylo[2] = 0;
  fpylo[3] = 0;
  fpyhi[0] = fpylo[0] + CCDWD_T + 1;
  fpyhi[1] = fpylo[1] + CCDWD_T + 1;
  fpyhi[2] = fpylo[2] + CCDWD_T + 1;
  fpyhi[3] = fpylo[3] + CCDWD_T + 1;
  /*
  for(iccd=0;iccd<4;++iccd) {
    fprintf(stderr,"%d %f %f %f %f\n",
	    iccd,fpxlo[iccd],fpxhi[iccd],fpylo[iccd],fpyhi[iccd]);
  }
  */

  if( (fitp[0] >= fpxlo[0]) ) { 
    // CCD 1 or 4
    if( (fitp[1] >= fpylo[0]) ) {
      iccd = 0;
      xy[0] = fitp[0] - fpxlo[iccd];
      xy[1] = fitp[1] - fpylo[iccd];
    }
    else {
      iccd = 3;
      xy[0] = fitp[0] - fpxlo[3];
      xy[1] = fitp[1] - fpyhi[3];
     }
  }
  else {
    if( (fitp[1] >= fpylo[0]) ) {
      iccd = 1;   
      xy[0] = fitp[0] - fpxhi[iccd];
      xy[1] = fitp[1] - fpylo[iccd];
    }
    else {
      iccd = 2;
      xy[0] = fitp[0] - fpxhi[iccd];
      xy[1] = fitp[1] - fpyhi[iccd];
    }
  }
  // CHECK PIXEL NUMBERING IN TESTING
  // xy[0] -= 1.0;
  // xy[1] -= 1.0;
  // CHECK SIGNS IN TESTING - REMOVE MINUS SIGNS UPON COMPLETION
  xy[0] = -xy[0];
  xy[1] = -xy[1];
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  double radecroll[3], eulc[3], vstar[3], vcam[3], lng, lat, lngd, latd;
  double rdr[3];
  double xyfp[2], ccdpx[2], fitpx[2];
  int icam, iccd, i, jccd, kccd, j;
  FILE *fpf;

  dtor = M_PI/180.0;

  fpf = fopen("fpg_pars.txt","r");
  icam = 0;
  read_fpg_pars(fpf,stderr,icam);

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

  read_stars(stdin,stderr);
  for(i=0;i<nst;++i) {
    // if(i > 3) break;
    dsphcr(dtor*ra[i],dtor*dec[i],vstar);
    matvec(rmat4,vstar,vcam);
    dcrsph(vcam,&lng,&lat);
    lngd = lng/dtor;
    latd = lat/dtor;
    if(star_in_fov(lngd,latd) == 1) {
      optics_fp(lngd,latd,xyfp,icam);
      // fprintf(stderr,"%d %f %f    %e %e    %e %e\n",i,ra[i],dec[i],lngd,latd,xyfp[0],xyfp[1]);
      iccd = mm_to_pix(icam,xyfp,ccdpx,fitpx);

      fprintf(stdout,"%4d %11.5f %11.5f   %11.5f %11.5f  %9.5f %9.5f %9.3f %9.3f"
	      " %2d %9.3f %9.3f\n",
	      i,ra[i],dec[i],lngd,90.0-latd,xyfp[0],xyfp[1],ccdpx[0],ccdpx[1],
	      iccd+1,fitpx[0],fitpx[1]);
    }
  }
}
/******************************************************************************/

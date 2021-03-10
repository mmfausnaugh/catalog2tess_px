/* guid1.c
 * Alan M. Levine
 * April 4, 2017
 *
 * Look at predicted guide star positions in the focal plane.
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

double row[NSTAR], col[NSTAR], bmin[NSTAR], bmax[NSTAR], ra[NSTAR], dec[NSTAR];
double latrm[NSTAR], lngrm[NSTAR], fxm[NSTAR], fym[NSTAR];
int nst;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3];
double ra_sc, dec_sc, roll_sc;

double dtor;
double camoff[4] = { -36.0, -12.0, 12.0, 36.0 };
double camoff2[4] = { -37.5, -12.5, 12.5, 37.5 };
void prmat();

/******************************************************************************/
double rmat1[3][3], rmat2[3][3], rmat3[3][3];

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

void sc_to_cam_mat(int icam, double euler[3])
{
  double xeul[3], angle, rm[3][3];
  int i;

  for(i=0;i<3;++i) {
    xeul[i] = dtor*euler[i];
  }

  eulerm323(xeul,rmat2);
  if (icam < 2) {
    angle = M_PI/2.0;
  }
  else {
    angle = -M_PI/2.0;
  }
  rotm1(2,angle,rm);
  matmat(rm,rmat2,rmat3);
  prmat(rmat2,stderr);
  prmat(rm,stderr);
  prmat(rmat3,stderr);
}

/******************************************************************************/

void prmat(double rm[3][3], FILE *fp)
{
  int i, j;

  for(i=0;i<3;++i) {
    fprintf(fp,"%d",i);
    for(j=0;j<3;++j) {
      fprintf(fp,"  %f",rm[i][j]);
    }
    fprintf(fp,"\n");
  }
}

/******************************************************************************/
//  TSIG/AML optics model

void optics_fp_t(double lng_deg, double lat_deg, double xyfp[2])
{
  double thetar, tanth, cphi, sphi, ttsq, rfp;

  thetar = (M_PI/2.0) - (lat_deg*dtor);
  tanth = tan(thetar);
  cphi = cos(dtor*lng_deg);
  sphi = sin(dtor*lng_deg);
  ttsq = tanth*tanth;
  rfp = FLMM*tanth*(1.00000140 + (0.28174612*ttsq) + (-0.59667259*ttsq*ttsq) + (9.17151267*ttsq*ttsq*ttsq) +
		  (-4.36928235*ttsq*ttsq*ttsq*ttsq));
  xyfp[0] = -cphi*rfp;
  xyfp[1] = -sphi*rfp;
}

/******************************************************************************/
// NASA Ames optics model
// This routine will only be accurate when a corner of the CCD active area is at the
// very center of the focal plane (projection of camera boresight).

double na_coef[7] = { 0.0049791024, 2851.175, 10.547364, 56.481487, 33.51133, 
		      -69.38692, 51.90312 };
double naorig = 0.0;
double naoff = 0.0;
double nascale = 1.654982e-5;

void optics_fp_a(double lng_deg, double lat_deg, double xyfp[2])
{
  double thetar, tanth, cphi, sphi, ttsq, rfp;
  double sv[3], eta, xi, rho, poly, ppwr, g;
  int ior, nor;

  dsphcr(lng_deg*dtor,lat_deg*dtor,sv);
  eta = -sv[0]/sv[2];
  xi = -sv[1]/sv[2];
  rho = sqrt(eta*eta + xi*xi);
  g = (rho*3600.0/dtor - naorig)*nascale + naoff;

  nor = 6;
  poly = 0.0;
  ppwr = 1.0;
  for(ior=0;ior<=nor;++ior) {
    poly += na_coef[ior]*ppwr;
    ppwr *= g;
  }

  cphi = cos(dtor*lng_deg);
  sphi = sin(dtor*lng_deg);

  xyfp[0] = eta*poly/rho;
  xyfp[1] = xi*poly/rho;
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

void read_guides(FILE *fpin, FILE *fpd)
{
  int i;

  i = 0;
  while(fscanf(fpin,"%lf %lf",&row[i],&col[i]) == 2) {
    fscanf(fpin,"%lf %lf %*lf %*lf %*d %lf %lf",&bmin[i],&bmax[i],&ra[i],&dec[i]);
    fprintf(fpd,"%f %f %f %f\n",row[i],col[i],ra[i],dec[i]);
    ++i;
  }
  nst = i;

}

/******************************************************************************/

void read_guides2(FILE *fpin, FILE *fpd)
{
  int i;

  i = 0;
  while(fscanf(fpin,"%lf %lf",&row[i],&col[i]) == 2) {
    fscanf(fpin,"%lf %lf %*lf %*lf %*d %lf %lf",&bmin[i],&bmax[i],&ra[i],&dec[i]);
    fprintf(fpd,"%f %f %f %f   ",row[i],col[i],ra[i],dec[i]);
    fscanf(fpin,"%*d %lf %lf %lf %lf",&latrm[i],&lngrm[i],&fxm[i],&fym[i]);
    fprintf(fpd,"%f %f %f %f\n",latrm[i]/dtor,lngrm[i]/dtor,fxm[i]/PIXEL,fym[i]/PIXEL);
    ++i;
  }
  nst = i;

}

/******************************************************************************/
/* Convert focal plane coordinates in pixel units to pseudo-equivalent TESS
 * camera image FITS file pixel numbers.
 */

// Cam num = icam + 1
// CCD num = iccd + 1

int fits_pix_t(int icam, double xy[2], double fitp[2])
{
  int iccd;
  double xyp[2];

  // CHECK SIGNS IN TESTING - REMOVE MINUS SIGNS UPON COMPLETION
  xyp[0] = -xy[0];
  xyp[1] = -xy[1];

  if(xyp[0] >= 0.0) {
    if(xyp[1] >= 0.0) {
      iccd = 0;
      fitp[0] = xyp[0] + CCDWD_T + 2*ROWA + ROWB;
      fitp[1] = xyp[1] + CCDHT_T + 2*COLDK_T;
    }
    else {
      iccd = 3;
      fitp[0] = xyp[0] + CCDWD_T + 2*ROWA + ROWB;
      fitp[1] = xyp[1] + CCDHT_T;
    }
  }
  else {
    if(xyp[1] >= 0.0) {
      iccd = 1;   
      fitp[0] = xyp[0] + CCDWD_T + ROWA;
      fitp[1] = xyp[1] + CCDHT_T + 2*COLDK_T;
    }
    else {
      iccd = 2;
      fitp[0] = xyp[0] + CCDWD_T + ROWA;
      fitp[1] = xyp[1] + CCDHT_T;
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
/* Convert focal plane coordinates in pixel units to pseudo-equivalent TESS
 * camera image FITS file pixel numbers.
 */

// Cam num = icam + 1
// CCD num = iccd + 1

int fits_pix_a(int icam, double xy[2], double fitp[2])
{
  int iccd;
  double xyp[2];

  // CHECK SIGNS IN TESTING - REMOVE MINUS SIGNS UPON COMPLETION
  xyp[0] = -xy[0];
  xyp[1] = -xy[1];

  if(xyp[0] >= 0.0) {
    if(xyp[1] >= 0.0) {
      iccd = 0;
      fitp[0] = xyp[0] + CCDWD_A + 2*ROWA + ROWB;
      fitp[1] = xyp[1] + CCDHT_A + 2*COLDK_A;
    }
    else {
      iccd = 3;
      fitp[0] = xyp[0] + CCDWD_A + 2*ROWA + ROWB;
      fitp[1] = xyp[1] + CCDHT_A;
    }
  }
  else {
    if(xyp[1] >= 0.0) {
      iccd = 1;   
      fitp[0] = xyp[0] + CCDWD_A + ROWA;
      fitp[1] = xyp[1] + CCDHT_A + 2*COLDK_A;
    }
    else {
      iccd = 2;
      fitp[0] = xyp[0] + CCDWD_A + ROWA;
      fitp[1] = xyp[1] + CCDHT_A;
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

void from_fits_pix_a(int icam, double xy[2], double fitp[2])
{
  int iccd;
  double xyp[2], fpxlo[4], fpxhi[4], fpylo[4], fpyhi[4];

  fpxlo[0] = CCDWD_A + 2*ROWA + ROWB;
  fpxlo[1] = ROWA;
  fpxlo[2] = ROWA;
  fpxlo[3] = CCDWD_A + 2*ROWA + ROWB;
  fpxhi[0] = fpxlo[0] + CCDWD_A + 1;
  fpxhi[1] = fpxlo[1] + CCDWD_A + 1;
  fpxhi[2] = fpxlo[2] + CCDWD_A + 1;
  fpxhi[3] = fpxlo[3] + CCDWD_A + 1;
  fpylo[0] = CCDHT_A + 2*COLDK_A;
  fpylo[1] = CCDHT_A + 2*COLDK_A;
  fpylo[2] = 0;
  fpylo[3] = 0;
  fpyhi[0] = fpylo[0] + CCDWD_A + 1;
  fpyhi[1] = fpylo[1] + CCDWD_A + 1;
  fpyhi[2] = fpylo[2] + CCDWD_A + 1;
  fpyhi[3] = fpylo[3] + CCDWD_A + 1;
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
  double xyfp[2], xypix[2], fitpna[2], r1, r2, phi1, phi2, xyfpa[2];
  double radecroll[3], eulcam[3], vstar[3], vcam[3], lng, lat, lngd, latd;
  double rdr[3], fitp[2], fitpb[2];
  int icam, iccd, i, jccd, kccd;

  dtor = M_PI/180.0;

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

  eulcam[0] = 0.0;
  eulcam[1] = camoff[icam];
  eulcam[2] = 0.0;
  sc_to_cam_mat(icam,eulcam);
  matmat(rmat3,rmat1,rmat4);
  prmat(rmat4,stderr);
  get_ra_dec_roll(rmat4,rdr);
  fprintf(stderr,"cam ra,dec,roll = %f %f %f\n",
	  rdr[0],rdr[1],rdr[2]);


  read_guides2(stdin,stderr);
  for(i=0;i<nst;++i) {
    // if(i > 3) break;
    dsphcr(dtor*ra[i],dtor*dec[i],vstar);
    matvec(rmat4,vstar,vcam);
    dcrsph(vcam,&lng,&lat);
    lngd = lng/dtor;
    latd = lat/dtor;
    optics_fp_t(lngd,latd,xyfp);
    xypix[0] = xyfp[0]/PIXEL;
    xypix[1] = xyfp[1]/PIXEL;
    optics_fp_a(lngd,latd,xyfpa);
    jccd = fits_pix_t(icam,xypix,fitp);
    kccd = fits_pix_t(icam,xyfpa,fitpb);
    fitpna[0] = col[i];
    fitpna[1] = row[i];
    from_fits_pix_t(icam,xyfp,fitpna);
    r1 = sqrt(xypix[0]*xypix[0] + xypix[1]*xypix[1]);
    phi1 = atan2(xypix[1],xypix[0])/dtor;
    r2 = sqrt(xyfp[0]*xyfp[0] + xyfp[1]*xyfp[1]);
    phi2 = atan2(xyfp[1],xyfp[0])/dtor;

    fprintf(stdout,"%4d %12.5f %12.5f     %12.5f %12.5f    %12.5f %12.5f   %12.5f %12.5f"
	    "%2d %12.5f %12.5f   %12.5f %12.5f\n",
	    i,ra[i],dec[i],col[i],row[i],xypix[0],xypix[1],90.0-latd,lngd,
	    jccd+1,fitp[0],fitp[1],xyfp[0],xyfp[1]);
    fprintf(stdout,"   %14.5f %14.5f    %14.5f %14.5f   %14.5f %14.5f    %14.5f %14.5f"
	    "    %14.5f %14.5f\n",
	    r1,phi1,r2,phi2,r1-r2,phi1-phi2,xyfpa[0],xyfpa[1],fitpb[0],fitpb[1]);
  }
}
/******************************************************************************/

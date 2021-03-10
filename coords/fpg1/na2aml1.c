/* na2aml1.c
 * Alan M. Levine
 * November 3, 2017
 * Heritage: aml2na1.c
 *
 * perform various computations related to deriving AML-type geometric parameters
 * from a set of NASA Ames type geometric parameters.
 *
 * Aug. 1, 2017:  Adopt flight software pixel numbering for fitpx[].
 *                Use zero-base for fitpx[] indices,
 *                and denote the center of a pixel by the pixel number exactly
 *                for floating point pixel locations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"

#define CCDWD_T 2048
#define CCDHT_T 2058
#define ROWA 44
#define ROWB 44
#define COLDK_T 20

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
double asymang[NCAM], asymfac[NCAM];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];

double col, row, ra, dec, ra_ab, dec_ab;
int icam;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3], rmat5[3][3];
double ra_sc, dec_sc, roll_sc;

double fpxlo[NCCD], fpxhi[NCCD], fpylo[NCCD], fpyhi[NCCD];

double dtor;
void prmat();

/******************************************************************************/
// All input angles are in degrees.

void read_fpg_pars(FILE *fpin, FILE *fpo)
{
  char pardesc[100];
  int nrd, iccd, j;

  for(j=0;j<3;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&keulcam[icam][j],&eulcam[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,keulcam[icam][j],eulcam[icam][j]);
  }
  for(j=0;j<6;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&koptcon[icam][j],&optcon[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,koptcon[icam][j],optcon[icam][j]);
  }
  fscanf(fpin,"%s %d %lf",pardesc,&kasymang[icam],&asymang[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymang[icam],asymang[icam]);
  fscanf(fpin,"%s %d %lf",pardesc,&kasymfac[icam],&asymfac[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymfac[icam],asymfac[icam]);
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdxy0[icam][iccd][j],&ccdxy0[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
    }
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kpixsz[icam][iccd][j],&pixsz[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %.8f\n",pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
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

void sc_to_cam_mat(double euler[3])
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
// stretch predicted focal plane position by 'asymfac' parallel to
// azimuthal angle 'asymang'

void make_az_asym(double xyin[2], double xyout[2])
{
  double xyp[2], xypa[2];

  xyrotate(asymang[icam],xyin,xyp);
  xypa[0] = asymfac[icam]*xyp[0];
  xypa[1] = xyp[1];
  xyrotate(-asymang[icam],xypa,xyout);
}

/******************************************************************************/
// correct for stretch predicted focal plane position by 'asymfac' parallel to
// azimuthal angle 'asymang'

void rev_az_asym(double xyin[2], double xyout[2])
{
  double xyp[2], xypa[2];

  xyrotate(asymang[icam],xyin,xyp);
  xypa[0] = xyp[0]/asymfac[icam];
  xypa[1] = xyp[1];
  xyrotate(-asymang[icam],xypa,xyout);
}

/******************************************************************************/
/*  AML optics model
 *
 * lng_deg, lat_deg are the spherical coordinates of the direction to the
 * star (degrees)
 *
 * xyfp[2] contains the x and y focal plane coordinates in mm
 */

void optics_fp(double lng_deg, double lat_deg, double xyfp[2])
{
  double thetar, tanth, cphi, sphi, ttsq, rfp0, rfp, xytmp[2];
  int i;

  thetar = (M_PI/2.0) - (lat_deg*dtor);
  tanth = tan(thetar);
  cphi = cos(dtor*lng_deg);
  sphi = sin(dtor*lng_deg);
  rfp0 = optcon[icam][0]*tanth;
  rfp = 0.0;
  for(i=1;i<NOPTCON;++i) {
    rfp += optcon[icam][i]*pow(tanth,2.0*(i-1));
  }
  xytmp[0] = -cphi*rfp0*rfp;
  xytmp[1] = -sphi*rfp0*rfp;
  make_az_asym(xytmp,xyfp);
}

/******************************************************************************/

double r_of_tanth(double tanth)
{
  double rfp0, rfp, pw, rp;
  int i;

  rfp0 = optcon[icam][0]*tanth;
  rfp = 0.0;
  // fprintf(stderr,"con,tanth,rfp0 = %e %e %e\n",optcon[icam][0],tanth,rfp0);
  for(i=1;i<NOPTCON;++i) {
    pw = pow(tanth,2.0*(i-1));
    rfp += optcon[icam][i]*pw;
    // fprintf(stderr,"i,con,pw,rfp = %d %e %e %e\n",i,optcon[icam][i],pw,rfp);
  }
  rp = rfp0*rfp;
  return(rp);
}

/******************************************************************************/

double dr_dtanth(double tanth)
{
  double z, drdz, c0, pw;
  int i;

  z = tanth;
  c0 = optcon[icam][0];
  drdz = 0.0;
  for(i=1;i<NOPTCON;++i) {
    pw = pow(z,2.0*(i - 1));
    drdz += (2.0*i - 1.0)*c0*optcon[icam][i]*pw;
    // fprintf(stderr,"i,con,pw,rfp = %d %e %e %e\n",i,optcon[icam][i],pw,drdz);
  }
  return(drdz);
}

/******************************************************************************/
// Get tanth(r) using Newton's method

double tanth_of_r(double r)
{
  double zi, ri, c0, drdz, delr, delz;

  if(fabs(r) < 1.0e-10)
    return(0.0);

  c0 = optcon[icam][0];
  zi = r/c0;
  ri = r_of_tanth(zi);
  delr = r - ri;
  while(fabs(delr) > 1.0e-10) {
    delz = delr/dr_dtanth(zi);
    fprintf(stderr,"r,ri,delr,delz,zi,zi+delz= %f %f %e %e %f %f\n",r,ri,delr,delz,zi,zi+delz);
    zi += delz;
    ri = r_of_tanth(zi);
    delr = r - ri;
  }
  fprintf(stderr,"r,ri,delr,zi= %f %f %e %20.10f\n",r,ri,delr,zi);

  return(zi);
}

/******************************************************************************/
/*
 * Input: xyfp[2] contains the x and y focal plane coordinates in mm
 *
 * Output: lnglat_deg[2] will have the spherical coordinates of the direction to the
 * star (degrees) in the camera frame.
 *
 */

void fp_optics(double lnglat_deg[2], double xyfp[2])
{
  double thetar, phir, thetad, phid, rfp, tanth;
  double xy[2];

  rev_az_asym(xyfp,xy);
  rfp = sqrt(xy[0]*xy[0] + xy[1]*xy[1]);
  phir = atan2(-xy[1],-xy[0]);
  phid = phir/dtor;
  if(phid < 0.0)
    phid += 360.0;

  tanth = tanth_of_r(rfp);
  thetar = atan(tanth);
  thetad = thetar/dtor;

  lnglat_deg[0] = phid;
  lnglat_deg[1] = 90.0 - thetad;
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
/* Convert focal plane coordinates in mm to pseudo-equivalent TESS
 * camera image CCD pixel numbers and FITS file pixel numbers.
 */

// CCD tilt is ignored
// No checking is done of whether the star is on a CCD active area

// Camera number = icam + 1
// CCD number = iccd + 1

// pixel numbers in fitpx[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention

int mm_to_pix(double xy[2], double ccdpx[2], double fitpx[2])
{
  int iccd;
  double xya[2], xyb[2], xyccd[2];

  xya[0] = xy[0];
  xya[1] = xy[1];

  if(xya[0] >= 0.0) {
    if(xya[1] >= 0.0) {
      iccd = 0;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (CCDWD_T - ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB - 1.0;
      fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0;
    }
    else {
      iccd = 3;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
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
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (CCDWD_T - ccdpx[0]) + ROWA - 1.0;
      fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0;
    }
    else {
      iccd = 2;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (ccdpx[0]) + ROWA;
      fitpx[1] = (ccdpx[1]);
    }
  }
  /* fprintf(stderr,"mm_to_pix(): %d %f %f      %f %f      %f %f      %f %f\n",
	  iccd,xy[0],xy[1],xyccd[0],xyccd[1],ccdpx[0],ccdpx[1],fitpx[0],fitpx[1]);
  */
  return(iccd);
}

/******************************************************************************/
// pixel numbers in fpxlo[] and fpxhi[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention

/*
#define CCDWD_T 2048
#define CCDHT_T 2058
#define ROWA 44
#define ROWB 44
#define COLDK_T 20
*/

void set_fpx_lims()
{
  int iccd;

  fpxlo[1] = ROWA - 0.5;
  fpxhi[1] = fpxlo[1] + CCDWD_T;
  fpxlo[2] = fpxlo[1];
  fpxhi[2] = fpxhi[1];
  fpxlo[0] = fpxhi[1] + ROWA + ROWB;  
  fpxhi[0] = fpxlo[0] + CCDWD_T;
  fpxlo[3] = fpxlo[0];
  fpxhi[3] = fpxhi[0];

  fpylo[2] = -0.5;
  fpyhi[2] = fpylo[2] + CCDHT_T;
  fpylo[3] = fpylo[2];
  fpyhi[3] = fpyhi[2];
  fpylo[0] = fpyhi[2] + 2*COLDK_T;
  fpyhi[0] = fpylo[0] + CCDHT_T;
  fpylo[1] = fpylo[0];
  fpyhi[1] = fpyhi[0];

  for(iccd=0;iccd<4;++iccd) {
    fprintf(stderr,"set_fpx_lims(): %d %f %f %f %f\n",
	    iccd,fpxlo[iccd],fpxhi[iccd],fpylo[iccd],fpyhi[iccd]);
  }
}

/******************************************************************************/
/* Convert pseudo-equivalent TESS camera image FITS file pixel numbers to
 * focal plane coordinates in mm.
 */

// Cam num = icam + 1
// CCD num = iccd + 1

// pixel numbers in fitpx[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention

int pix_to_mm(double xy[2], double fitpx[2])
{
  int iccd;
  double ccdpx[2];
  double xyb[2], xyccd[2];

  if( (fitpx[0] >= fpxlo[0]) ) { 
    // CCD 1 or 4
    if( (fitpx[1] >= fpylo[0]) ) {
      iccd = 0;
    }
    else {
      iccd = 3;
     }
  }
  else {
    if( (fitpx[1] >= fpylo[0]) ) {
      iccd = 1;   
    }
    else {
      iccd = 2;
    }
  }

  if(iccd < 2) {  // CCD no. = 1 or 2
    ccdpx[0] = fpxhi[iccd] - fitpx[0];
    ccdpx[1] = fpyhi[iccd] - fitpx[1];
  }
  else {          // CCD no. = 3 or 4
    ccdpx[0] = fitpx[0] - fpxlo[iccd];
    ccdpx[1] = fitpx[1] - fpylo[iccd];
  }
  xyccd[0] = ccdpx[0]*pixsz[icam][iccd][0];
  xyccd[1] = ccdpx[1]*pixsz[icam][iccd][1];
  /*fprintf(stderr,"AA: %f %f      %.8f %.8f      %f %f \n",
	  ccdpx[0],ccdpx[1],pixsz[icam][iccd][0],pixsz[icam][iccd][1],xyccd[0],xyccd[1]);
  */
  xyrotate(-ccdang[icam][iccd],xyccd,xyb);
  xy[0] = xyb[0] + ccdxy0[icam][iccd][0];
  xy[1] = xyb[1] + ccdxy0[icam][iccd][1];
  /*fprintf(stderr,"pix_to_mm(): %d %f %f      %f %f      %f %f      %f %f\n",
	  iccd,ccdpx[0],ccdpx[1],xyccd[0],xyccd[1],xyb[0],xyb[1],xy[0],xy[1]);
  */
  return(iccd);
}

/******************************************************************************/
// Read RA, Dec, TESS magnitude, and color (to be inserted later) for each star

int read_pixel_coords(FILE *fpin, FILE *fpd)
{
  int nrd;

  nrd = fscanf(fpin,"%lf %lf",&col,&row);
  return(nrd);
}

/******************************************************************************/
/* Coordinate convention:
 * The center of a pixel has pixel coordinates equal to the pixel indices.
 * lo = coordinate value at low edge of low pixel
 * hi = coordinate value at high edge of high pixel
 */

void mm_ccd_corners(int iccd, double xycorn[4][2])
{
  int ic, jccd;
  double fpx[2], xycc[2];

  for(ic=0;ic<4;++ic) {
    if( (ic==0) || (ic==1) )
      fpx[0] = fpxlo[iccd];
    else {
      fpx[0] = fpxhi[iccd];
    }
    if( (ic==0) || (ic==3) )
      fpx[1] = fpylo[iccd];
    else {
      fpx[1] = fpyhi[iccd];
    }
    jccd = pix_to_mm(xycc,fpx);
    fprintf(stderr,"mm_ccd_corners(): iccd,jccd,fpx,xycc = ");
    fprintf(stderr,"   %d %d %f %f   %f %f\n",
	    iccd,jccd,fpx[0],fpx[1],xycc[0],xycc[1]);
    xycorn[ic][0] = xycc[0];
    xycorn[ic][1] = xycc[1];
  }
}

/******************************************************************************/
// All input angles are in degrees.

#define NPMAX 7
double xang[NCCD], yang[NCCD], zang[NCCD];
double offset[NCCD], scale[NCCD], origin;[NCCD]
int fnum[NCCD], rnum[NCCD];
double fcoeff[NCCD][NPMAX], rcoeff[NCCD][NPMAX];

void read_na_xml(FILE *fpin, FILE *fpo)
{
  char pardesc[100];
  int nrd, iccd, j;

  for(j=0;j<3;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&keulcam[icam][j],&eulcam[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,keulcam[icam][j],eulcam[icam][j]);
  }
  for(j=0;j<6;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&koptcon[icam][j],&optcon[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,koptcon[icam][j],optcon[icam][j]);
  }
  fscanf(fpin,"%s %d %lf",pardesc,&kasymang[icam],&asymang[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymang[icam],asymang[icam]);
  fscanf(fpin,"%s %d %lf",pardesc,&kasymfac[icam],&asymfac[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymfac[icam],asymfac[icam]);
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdxy0[icam][iccd][j],&ccdxy0[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
    }
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kpixsz[icam][iccd][j],&pixsz[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %.8f\n",pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
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

void get_ames_model(int npoly, double **polyco, double sck, double eta,
		    double xi, double xyr[3])
{
  int ipol;
  double rho, rhok, poly, x, y, xymag, z;

  rho = sqrt(eta*eta + xi*xi);
  rhok = sck*rho;
  // fprintf(stderr,"g_a_m(): eta,xi,rho,rhok = %f %f %f %f\n",eta,xi,rho,rhok);
  poly = 0.0;
  for(ipol=0;ipol<npoly;++ipol) {
    poly += polyco[ipol][0]*pow(rhok,ipol);
    /* fprintf(stderr,"ipol,coeff,pow,poly = %d %f %f %f\n",
	    ipol,polyco[ipol][0],pow(rhok,ipol),poly);
    */
  }
  xyr[2] = poly;
  z = 1.0/sqrt(1.0 + eta*eta + xi*xi);
  x = eta*z;
  y = xi*z;
  xymag = sqrt(x*x + y*y);
  xyr[0] = -(x/xymag)*poly;
  xyr[1] = -(y/xymag)*poly;
}

/******************************************************************************/

int old_main(int argc, char *argv[])
{
  double eulc[3], eulna[3], eulnad[3];
  double xyfp[2], ccdpx[2], fitpx[2];
  double xyc[4][2], xycc[2], xyin[4][2];
  double xcam[3], ycam[3], zcam[3], xccd[3], yccd[3], zccd[3], czang, zangr;
  double pole[3], poly;
  double lng, lngd, lat, latd;
  double xfad, yfad, txfa, tyfa, v[3], vccd[3];
  double eta, xi, rho, rhok, psir, px, py, ral, sck, sckin, dmt, xfaoff, yfaoff;
  double scalex, inscalex;
  double **cmat, **cvec, **dmat, **dvec;
  double **pxaml, **pyaml, **pxna, **pyna, **etana, **xina, xyr[3], **rpxaml, **rhona;
  double rmssq, dpx, dpy, dsq, dpmax;
  int ifa, jfa, nmat, ka, ib, jg, nfa;
  int iccd, i, jccd, kccd, j, ic, ipol;
  int nval;
  FILE *fpf;

  dtor = M_PI/180.0;
  set_fpx_lims();

                           // move below reading of command line parameters
  fpf = fopen("fpg_pars.txt","r");   
  icam = 0;
  read_fpg_pars(fpf,stderr);
  fclose(fpf);

  if(argc != 2) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  icam = atoi(argv[1]) - 1;   //   RENAME from icam
  fprintf(stdout,"Camera no. %d\n\n",icam+1);
  fprintf(stderr,"Camera no. %d\n\n",icam+1);
  icam = 0;

  for(j=0;j<3;++j) {
    eulc[j] = eulcam[icam][j];
  }
  sc_to_cam_mat(eulc);   // puts rotation matrix in rmat2

  scalex = 1.654982e-5;
  sck = (3600.0/dtor)*scalex;    // = 3.4136...  
  nmat = 7;
  cmat = matrix(nmat,nmat);
  cvec = matrix(nmat,1);
  dmat = matrix(nmat,nmat);
  dvec = matrix(nmat,1);
  xfaoff = 0.1;
  yfaoff = 0.1;

  nfa = 13;
  pxaml = matrix(nfa,nfa);
  pyaml = matrix(nfa,nfa);
  pxna = matrix(nfa,nfa);
  pyna = matrix(nfa,nfa);
  etana = matrix(nfa,nfa);
  xina = matrix(nfa,nfa);
  rhona = matrix(nfa,nfa);
  rpxaml = matrix(nfa,nfa);

  for(iccd=0;iccd<4;++iccd) {
    fprintf(stdout,"\n");
    fprintf(stdout,"CCD no. %d\n\n",iccd+1);
    fprintf(stderr,"\n");
    fprintf(stderr,"CCD no. %d\n\n",iccd+1);

    // find x, y focal plane coordinates of the corners of the active area
    // of CCD no. 'iccd'
    mm_ccd_corners(iccd,xyc);
    for(ic=0;ic<4;++ic) {
      xycc[0] = xyc[ic][0];
      xycc[1] = xyc[ic][1];
      jccd = mm_to_pix(xycc,ccdpx,fitpx);
    }
    if(iccd==0) {
      xyin[iccd][0] = xyc[0][0];
      xyin[iccd][1] = xyc[0][1];
    }
    else if(iccd==1) {
      xyin[iccd][0] = xyc[3][0];
      xyin[iccd][1] = xyc[3][1];
    }
    else if(iccd==2) {
      xyin[iccd][0] = xyc[2][0];
      xyin[iccd][1] = xyc[2][1];
    }
    else if(iccd==3) {
      xyin[iccd][0] = xyc[1][0];
      xyin[iccd][1] = xyc[1][1];
    }
    fprintf(stdout,"iccd = CCD no. - 1 = %d\n",iccd);
    fprintf(stdout,"Focal plane coordinates of centermost corner of active area of this CCD:\n");
    fprintf(stdout,"x, y (xyin[]) = %f %f\n",xyin[iccd][0],xyin[iccd][1]);

    /* construct 'z_ccd' unit vectors
       compute x_ccd and y_ccd unit vectors:
       pole = norm(z_ccd cross zcam)
       cos(zangr) = dot(z_ccd,zcam)
       rotate x, y by zangr around pole in correct sense
    */

    xcam[0] = 1.0;
    xcam[1] = 0.0;
    xcam[2] = 0.0;
    ycam[0] = 0.0;
    ycam[1] = 1.0;
    ycam[2] = 0.0;
    zcam[0] = 0.0;
    zcam[1] = 0.0;
    zcam[2] = 1.0;

    zccd[0] = -xyin[iccd][0];
    zccd[1] = -xyin[iccd][1];
    zccd[2] = optcon[icam][0];
    dnorm(zccd);

    czang = dot(zcam,zccd);
    zangr = acos(czang);
    dcross(zcam,zccd,pole);
    dnorm(pole);
    drotate(xcam,pole,zangr,xccd);
    dnorm(xccd);
    drotate(ycam,pole,zangr,yccd);
    dnorm(yccd);

    /* Compute xAngle, yAngle, zAngle ala Ames using my geometrical model.
     *
     * First compute camera to CCD rotation matrix
     */

    for(i=0;i<3;++i) {
      rmat3[0][i] = xccd[i];
      rmat3[1][i] = yccd[i];
      rmat3[2][i] = zccd[i];
    }
    prmat(rmat3,stderr);
    matmat(rmat3,rmat2,rmat4);
    prmat(rmat4,stderr);
    mateuler123(rmat4,eulna);
    for(i=0;i<3;++i) {
      eulnad[i] = eulna[i]/dtor;
    }
    fprintf(stdout,"Ames Euler angles (x-y-z; degrees) = %f %f %f\n",
	    eulnad[0],eulnad[1],eulnad[2]);
    fprintf(stderr,"Ames Euler angles (x-y-z; degrees) = %f %f %f\n",
	    eulnad[0],eulnad[1],eulnad[2]);


    /* Now get NASA-Ames polynomial by fitting my model results (fitpx[2]) 
     * as a function of the tangent of the angle of the input rays from zccd[].
     */

    for(ka=0;ka<nmat;++ka) {
      for(ib=0;ib<nmat;++ib) {
	cmat[ka][ib] = 0.0;
      }
      cvec[ka][0] = 0.0;
    }

    for(ifa=0;ifa<13;++ifa) {
      xfad = ifa + xfaoff;
      txfa = tan(dtor*xfad);
      if( (iccd==0) || (iccd==3) )
	txfa *= -1.0;
      for(jfa=0;jfa<13;++jfa) {
	yfad = jfa + yfaoff;
	tyfa = tan(dtor*yfad);
	if( (iccd==0) || (iccd==1) )
	  tyfa *= -1.0;
	v[0] = txfa;
	v[1] = tyfa;
	v[2] = 1.0;
	dnorm(v);

	vccd[0] = dot(v,xccd);
	vccd[1] = dot(v,yccd);
	vccd[2] = dot(v,zccd);
	dnorm(vccd);
	eta = vccd[0]/vccd[2];
	xi = vccd[1]/vccd[2];
	etana[ifa][jfa] = eta;
	xina[ifa][jfa] = xi;
	rho = sqrt(eta*eta + xi*xi);
	rhona[ifa][jfa] = rho;
	rhok = sck*rho;

	dcrsph(v,&lng,&lat);
	lngd = lng/dtor;
	latd = lat/dtor;
	/* fprintf(stderr,"A %d %d %d %f %f       %f %f\n",
		iccd,ifa,jfa,xfad,yfad,lngd,latd);
	*/

	pxaml[ifa][jfa] = 9999.0;  // flag that value was not set
	pyaml[ifa][jfa] = 9999.0;
	rpxaml[ifa][jfa] = 9999.0;
	pxna[ifa][jfa] = 9999.0;
	pyna[ifa][jfa] = 9999.0;
	if(star_in_fov(lngd,latd) == 1) {
	  optics_fp(lngd,latd,xyfp);
	  kccd = mm_to_pix(xyfp,ccdpx,fitpx);
	  /* fprintf(stderr,"AA %d    %f %f    %f %f       %f %f       %f %f\n",
		  kccd,lngd,latd,xyfp[0],xyfp[1],ccdpx[0],ccdpx[1],fitpx[0],fitpx[1]);
	  */
	  if( (kccd == iccd) && (fitpx[0] >= fpxlo[iccd]) && (fitpx[0] <= fpxhi[iccd]) &&
	      (fitpx[1] >= fpylo[iccd]) && (fitpx[1] <= fpyhi[iccd]) ) {
	    // compute properly conditioned pixel numbers px, py
	    if( (iccd == 0) || (iccd == 2) ) {
	      px = 2048.0 - ccdpx[0];
	      py = 2058.0 - ccdpx[1];
	    }
	    else {
	      px = ccdpx[0];
	      py = 2058.0 - ccdpx[1];
	    }
	    pxaml[ifa][jfa] = px;
	    pyaml[ifa][jfa] = py;
	    ral = sqrt(px*px + py*py);
	    rpxaml[ifa][jfa] = ral;
	    /* fprintf(stderr,"B %d    %f %f    %f %f %f        %f\n",
		    kccd,fitpx[0],fitpx[1],px,py,ral,rhok);
	    */

	    for(ka=0;ka<nmat;++ka) {
	      for(ib=0;ib<nmat;++ib) {
		dmt = pow(rhok,ka+ib);
		// fprintf(stderr,"D %d %d    %f %f %f\n",ka,ib,rhok,dmt,ral);
		cmat[ka][ib] += pow(rhok,ka+ib);
	      }
	      cvec[ka][0] += ral*pow(rhok,ka);
	    }
	  } // if(kccd==iccd
	} // if(star_in_fov...
	else {
	}
      } // for(jfa...
    } // for(ifa...
    for(ka=0;ka<nmat;++ka) {
      fprintf(stderr,"%f           ",cvec[ka][0]);
      for(ib=0;ib<nmat;++ib) {
	fprintf(stderr,"%f ",cmat[ka][ib]);
      }
      fprintf(stderr,"\n");
    }
    fflush(stderr);
    jg = gaussj(cmat,nmat,cvec,1);
    if(jg == 1) {
      // print message
      fprintf(stderr,"iccd=%d; singular matrix\n",iccd);
    }
    else {
      // print cmat, cvec
      printf("\n");
      printf("plateScalePoly\n");
      printf("originx = %e\n",0.0);
      printf("scalex = %e\n",scalex);
      printf("offsetx = %e\n",0.0);
      printf("power    poly_coeff   matrix_row\n");
      for(ka=0;ka<nmat;++ka) {
	printf("%d      %f         ",ka,cvec[ka][0]);
	fprintf(stderr,"%f           ",cvec[ka][0]);
	for(ib=0;ib<nmat;++ib) {
	  printf("%f ",cmat[ka][ib]);
	  fprintf(stderr,"%f ",cmat[ka][ib]);
	}
	printf("\n");
	fprintf(stderr,"\n");
      }
      printf("\n");
      fprintf(stderr,"\n");
    }

    /* Now get NASA-Ames inverse polynomial by fitting rhona[][] 
     * as a function of rpxaml[][].
     */

    for(ka=0;ka<nmat;++ka) {
      for(ib=0;ib<nmat;++ib) {
	dmat[ka][ib] = 0.0;
      }
      dvec[ka][0] = 0.0;
    }

    /* Compute scale factor for use with inverse polynomial calculations
     */
    poly = 0.0;
    for(ipol=0;ipol<nmat;++ipol) {
      poly += cvec[ipol][0];
    }
    sckin = 1.0/poly;
    fprintf(stderr,"inverse sum poly = %f\n",poly);
    fprintf(stderr,"inverse scale factor = %f\n",sckin);

    for(ifa=0;ifa<13;++ifa) {
      for(jfa=0;jfa<13;++jfa) {
	rho = rhona[ifa][jfa];
	ral = rpxaml[ifa][jfa];
	if(fabs(ral) > 9000.0)
	  continue;
	for(ka=0;ka<nmat;++ka) {
	  for(ib=0;ib<nmat;++ib) {
	    dmt = pow(rhok,ka+ib);
	    // fprintf(stderr,"D %d %d    %f %f %f\n",ka,ib,rhok,dmt,ral);
	    dmat[ka][ib] += pow(ral*sckin,ka+ib);
	  }
	  dvec[ka][0] += (3600.0/dtor)*rho*pow(ral*sckin,ka);
	}
      } // for(jfa...
    } // for(ifa...
    for(ka=0;ka<nmat;++ka) {
      fprintf(stderr,"%f           ",dvec[ka][0]);
      for(ib=0;ib<nmat;++ib) {
	fprintf(stderr,"%f ",dmat[ka][ib]);
      }
      fprintf(stderr,"\n");
    }
    fflush(stderr);
    jg = gaussj(dmat,nmat,dvec,1);
    if(jg == 1) {
      // print message
      fprintf(stderr,"inverse polynomial code - iccd=%d; singular matrix\n",iccd);
    }
    else {
      // print dmat, dvec
      printf("inversePlateScalePoly\n");
      printf("originx = %e\n",0.0);
      printf("scalex = %e\n",sckin);
      printf("offsetx = %e\n",0.0);
      printf("power    poly_coeff   matrix_row\n");
      for(ka=0;ka<nmat;++ka) {
	printf("%d     %f        ",ka,dvec[ka][0]);
	fprintf(stderr,"%f           ",dvec[ka][0]);
	for(ib=0;ib<nmat;++ib) {
	  printf("%f ",dmat[ka][ib]);
	  fprintf(stderr,"%f ",dmat[ka][ib]);
	}
	printf("\n");
	fprintf(stderr,"\n");
      }
      printf("\n");
      fprintf(stderr,"\n");
    }

    // compute Ames pixel numbers and then Ames - AML pixel number differences
    nval = 0;
    rmssq = 0.0;
    dpmax = 0.0;
    for(ifa=0;ifa<nfa;++ifa) {
       for(jfa=0;jfa<nfa;++jfa) {
	 get_ames_model(nmat,cvec,sck,etana[ifa][jfa],xina[ifa][jfa],xyr);
	 if( (iccd==0) || (iccd==3) ) 
	   pxna[ifa][jfa] = xyr[0];
	 else
	   pxna[ifa][jfa] = -xyr[0];
	 if( (iccd==0) || (iccd==1) ) 
	   pyna[ifa][jfa] = xyr[1];
	 else
	   pyna[ifa][jfa] = -xyr[1];
	 if(pxaml[ifa][jfa] != 9999.0) {
	   /* fprintf(stderr,"XX %d %d %f %f      %f %f\n",
		   ifa,jfa,pxaml[ifa][jfa],pyaml[ifa][jfa],pxna[ifa][jfa],pyna[ifa][jfa]);
	   */
	   dpx = pxaml[ifa][jfa] - pxna[ifa][jfa];
	   dpy = pyaml[ifa][jfa] - pyna[ifa][jfa];
	   dsq = dpx*dpx + dpy*dpy;
	   if(dsq > dpmax)
	     dpmax = dsq;
	   ++nval;
	   rmssq += dsq;
	 }
       }
    }
    rmssq = sqrt(rmssq/nval);
    dpmax = sqrt(dpmax);
    fprintf(stderr,"nval,rmssq,dpmax = %d %f %f\n",nval,rmssq,dpmax);

  }  // for(iccd...

}

/******************************************************************************/

int main(int argc, char *argv[])
{
  double eulc[3], eulna[3], eulnad[3];
  double xyfp[2], ccdpx[2], fitpx[2];
  double xyc[4][2], xycc[2], xyin[4][2];
  double xcam[3], ycam[3], zcam[3], xccd[3], yccd[3], zccd[3], czang, zangr;
  double pole[3], poly;
  double lng, lngd, lat, latd;
  double xfad, yfad, txfa, tyfa, v[3], vccd[3];
  double eta, xi, rho, rhok, psir, px, py, ral, sck, sckin, dmt, xfaoff, yfaoff;
  double scalex, inscalex;
  double **cmat, **cvec, **dmat, **dvec;
  double **pxaml, **pyaml, **pxna, **pyna, **etana, **xina, xyr[3], **rpxaml, **rhona;
  double rmssq, dpx, dpy, dsq, dpmax;
  int ifa, jfa, nmat, ka, ib, jg, nfa;
  int iccd, i, jccd, kccd, j, ic, ipol;
  int nval;
  FILE *fpf;

  dtor = M_PI/180.0;
  set_fpx_lims();

  if(argc != 2) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  icam = atoi(argv[1]) - 1;
  fprintf(stdout,"Camera no. %d\n\n",icam+1);
  fprintf(stderr,"Camera no. %d\n\n",icam+1);

  fpf = fopen("na_pars.txt","r");
  read_na_xml(fpf,stderr);
  fclose(fpf);

  // Average xang, yang, xang over 4 CCDs to get angles that approximately describe the
  // camera boresight.
  // Convert these angles to AML camera Euler angles,
  // i.e., convert 123 angles to 323 angles.

  double naeul[3], rmat1;
  int ieul;

  for(ieul=0;ieul<3;++ieul) {
    naeul[ieul] = 0.0;
  }
  for(iccd=0;iccd<4;++iccd) {
    naeul[0] += xang[iccd];
    naeul[1] += yang[iccd];
    naeul[2] += zang[iccd];
  }
  for(ieul=0;ieul<3;++ieul) {
    naeul[ieul] /= 4.0;
  }
  eulerm123(naeul,rmat1);
  mateuler323(rmat1,&eulcam[icam][0]);
  for(ieul=0;ieul<3;++ieul) {
    eulcam[icam][ieul] /= dtor;
  }
  // correct for 90 degree offsets due to differences in coordinate frame 
  // definitions
  if(icam < 2)
    eulcam[icam][2] += 90.0;
  else
    eulcam[icam][2] += -90.0;

  // Print camera boresight Euler angles


  // Determine focal plane locations of CCD innermost corners.
  // Set by focal length, boresight direction, and zccd vectors
  // Get focal length
  // Get boresight direction from rmat1

  // CCD rotation angles?
  // Find angle between, e.g., xcam[3] and the part of xccd[iccd][3] that is normal to
  // the boresight.

  // Determine locations of the other CCD corners.


  // Get optical radial distortion coefficients

  // pixel size - may depend on focal length?

  // set tilts = 0.0
  // set asymang = 0.0, asymfac = 1.0
















  for(j=0;j<3;++j) {
    eulc[j] = eulcam[icam][j];
  }
  sc_to_cam_mat(eulc);   // puts rotation matrix in rmat2

  scalex = 1.654982e-5;
  sck = (3600.0/dtor)*scalex;
  nmat = 7;
  cmat = matrix(nmat,nmat);
  cvec = matrix(nmat,1);
  dmat = matrix(nmat,nmat);
  dvec = matrix(nmat,1);
  xfaoff = 0.1;
  yfaoff = 0.1;

  nfa = 13;
  pxaml = matrix(nfa,nfa);
  pyaml = matrix(nfa,nfa);
  pxna = matrix(nfa,nfa);
  pyna = matrix(nfa,nfa);
  etana = matrix(nfa,nfa);
  xina = matrix(nfa,nfa);
  rhona = matrix(nfa,nfa);
  rpxaml = matrix(nfa,nfa);

  for(iccd=0;iccd<4;++iccd) {
    fprintf(stdout,"\n");
    fprintf(stdout,"CCD no. %d\n\n",iccd+1);
    fprintf(stderr,"\n");
    fprintf(stderr,"CCD no. %d\n\n",iccd+1);

    // find x, y focal plane coordinates of the corners of the active area
    // of CCD no. 'iccd'
    mm_ccd_corners(iccd,xyc);
    for(ic=0;ic<4;++ic) {
      xycc[0] = xyc[ic][0];
      xycc[1] = xyc[ic][1];
      jccd = mm_to_pix(xycc,ccdpx,fitpx);
    }
    if(iccd==0) {
      xyin[iccd][0] = xyc[0][0];
      xyin[iccd][1] = xyc[0][1];
    }
    else if(iccd==1) {
      xyin[iccd][0] = xyc[3][0];
      xyin[iccd][1] = xyc[3][1];
    }
    else if(iccd==2) {
      xyin[iccd][0] = xyc[2][0];
      xyin[iccd][1] = xyc[2][1];
    }
    else if(iccd==3) {
      xyin[iccd][0] = xyc[1][0];
      xyin[iccd][1] = xyc[1][1];
    }
    fprintf(stdout,"iccd = CCD no. - 1 = %d\n",iccd);
    fprintf(stdout,"Focal plane coordinates of centermost corner of active area of this CCD:\n");
    fprintf(stdout,"x, y (xyin[]) = %f %f\n",xyin[iccd][0],xyin[iccd][1]);

    /* construct 'z_ccd' unit vectors
       compute x_ccd and y_ccd unit vectors:
       pole = norm(z_ccd cross zcam)
       cos(zangr) = dot(z_ccd,zcam)
       rotate x, y by zangr around pole in correct sense
    */

    xcam[0] = 1.0;
    xcam[1] = 0.0;
    xcam[2] = 0.0;
    ycam[0] = 0.0;
    ycam[1] = 1.0;
    ycam[2] = 0.0;
    zcam[0] = 0.0;
    zcam[1] = 0.0;
    zcam[2] = 1.0;

    zccd[0] = -xyin[iccd][0];
    zccd[1] = -xyin[iccd][1];
    zccd[2] = optcon[icam][0];
    dnorm(zccd);

    czang = dot(zcam,zccd);
    zangr = acos(czang);
    dcross(zcam,zccd,pole);
    dnorm(pole);
    drotate(xcam,pole,zangr,xccd);
    dnorm(xccd);
    drotate(ycam,pole,zangr,yccd);
    dnorm(yccd);

    /* Compute xAngle, yAngle, zAngle ala Ames using my geometrical model.
     *
     * First compute camera to CCD rotation matrix
     */

    for(i=0;i<3;++i) {
      rmat3[0][i] = xccd[i];
      rmat3[1][i] = yccd[i];
      rmat3[2][i] = zccd[i];
    }
    prmat(rmat3,stderr);
    matmat(rmat3,rmat2,rmat4);
    prmat(rmat4,stderr);
    mateuler123(rmat4,eulna);
    for(i=0;i<3;++i) {
      eulnad[i] = eulna[i]/dtor;
    }
    fprintf(stdout,"Ames Euler angles (x-y-z; degrees) = %f %f %f\n",
	    eulnad[0],eulnad[1],eulnad[2]);
    fprintf(stderr,"Ames Euler angles (x-y-z; degrees) = %f %f %f\n",
	    eulnad[0],eulnad[1],eulnad[2]);


    /* Now get NASA-Ames polynomial by fitting my model results (fitpx[2]) 
     * as a function of the tangent of the angle of the input rays from zccd[].
     */

    for(ka=0;ka<nmat;++ka) {
      for(ib=0;ib<nmat;++ib) {
	cmat[ka][ib] = 0.0;
      }
      cvec[ka][0] = 0.0;
    }

    for(ifa=0;ifa<13;++ifa) {
      xfad = ifa + xfaoff;
      txfa = tan(dtor*xfad);
      if( (iccd==0) || (iccd==3) )
	txfa *= -1.0;
      for(jfa=0;jfa<13;++jfa) {
	yfad = jfa + yfaoff;
	tyfa = tan(dtor*yfad);
	if( (iccd==0) || (iccd==1) )
	  tyfa *= -1.0;
	v[0] = txfa;
	v[1] = tyfa;
	v[2] = 1.0;
	dnorm(v);

	vccd[0] = dot(v,xccd);
	vccd[1] = dot(v,yccd);
	vccd[2] = dot(v,zccd);
	dnorm(vccd);
	eta = vccd[0]/vccd[2];
	xi = vccd[1]/vccd[2];
	etana[ifa][jfa] = eta;
	xina[ifa][jfa] = xi;
	rho = sqrt(eta*eta + xi*xi);
	rhona[ifa][jfa] = rho;
	rhok = sck*rho;

	dcrsph(v,&lng,&lat);
	lngd = lng/dtor;
	latd = lat/dtor;
	/* fprintf(stderr,"A %d %d %d %f %f       %f %f\n",
		iccd,ifa,jfa,xfad,yfad,lngd,latd);
	*/

	pxaml[ifa][jfa] = 9999.0;  // flag that value was not set
	pyaml[ifa][jfa] = 9999.0;
	rpxaml[ifa][jfa] = 9999.0;
	pxna[ifa][jfa] = 9999.0;
	pyna[ifa][jfa] = 9999.0;
	if(star_in_fov(lngd,latd) == 1) {
	  optics_fp(lngd,latd,xyfp);
	  kccd = mm_to_pix(xyfp,ccdpx,fitpx);
	  /* fprintf(stderr,"AA %d    %f %f    %f %f       %f %f       %f %f\n",
		  kccd,lngd,latd,xyfp[0],xyfp[1],ccdpx[0],ccdpx[1],fitpx[0],fitpx[1]);
	  */
	  if( (kccd == iccd) && (fitpx[0] >= fpxlo[iccd]) && (fitpx[0] <= fpxhi[iccd]) &&
	      (fitpx[1] >= fpylo[iccd]) && (fitpx[1] <= fpyhi[iccd]) ) {
	    // compute properly conditioned pixel numbers px, py
	    if( (iccd == 0) || (iccd == 2) ) {
	      px = 2048.0 - ccdpx[0];
	      py = 2058.0 - ccdpx[1];
	    }
	    else {
	      px = ccdpx[0];
	      py = 2058.0 - ccdpx[1];
	    }
	    pxaml[ifa][jfa] = px;
	    pyaml[ifa][jfa] = py;
	    ral = sqrt(px*px + py*py);
	    rpxaml[ifa][jfa] = ral;
	    /* fprintf(stderr,"B %d    %f %f    %f %f %f        %f\n",
		    kccd,fitpx[0],fitpx[1],px,py,ral,rhok);
	    */

	    for(ka=0;ka<nmat;++ka) {
	      for(ib=0;ib<nmat;++ib) {
		dmt = pow(rhok,ka+ib);
		// fprintf(stderr,"D %d %d    %f %f %f\n",ka,ib,rhok,dmt,ral);
		cmat[ka][ib] += pow(rhok,ka+ib);
	      }
	      cvec[ka][0] += ral*pow(rhok,ka);
	    }
	  } // if(kccd==iccd
	} // if(star_in_fov...
	else {
	}
      } // for(jfa...
    } // for(ifa...
    for(ka=0;ka<nmat;++ka) {
      fprintf(stderr,"%f           ",cvec[ka][0]);
      for(ib=0;ib<nmat;++ib) {
	fprintf(stderr,"%f ",cmat[ka][ib]);
      }
      fprintf(stderr,"\n");
    }
    fflush(stderr);
    jg = gaussj(cmat,nmat,cvec,1);
    if(jg == 1) {
      // print message
      fprintf(stderr,"iccd=%d; singular matrix\n",iccd);
    }
    else {
      // print cmat, cvec
      printf("\n");
      printf("plateScalePoly\n");
      printf("originx = %e\n",0.0);
      printf("scalex = %e\n",scalex);
      printf("offsetx = %e\n",0.0);
      printf("power    poly_coeff   matrix_row\n");
      for(ka=0;ka<nmat;++ka) {
	printf("%d      %f         ",ka,cvec[ka][0]);
	fprintf(stderr,"%f           ",cvec[ka][0]);
	for(ib=0;ib<nmat;++ib) {
	  printf("%f ",cmat[ka][ib]);
	  fprintf(stderr,"%f ",cmat[ka][ib]);
	}
	printf("\n");
	fprintf(stderr,"\n");
      }
      printf("\n");
      fprintf(stderr,"\n");
    }

    /* Now get NASA-Ames inverse polynomial by fitting rhona[][] 
     * as a function of rpxaml[][].
     */

    for(ka=0;ka<nmat;++ka) {
      for(ib=0;ib<nmat;++ib) {
	dmat[ka][ib] = 0.0;
      }
      dvec[ka][0] = 0.0;
    }

    /* Compute scale factor for use with inverse polynomial calculations
     */
    poly = 0.0;
    for(ipol=0;ipol<nmat;++ipol) {
      poly += cvec[ipol][0];
    }
    sckin = 1.0/poly;
    fprintf(stderr,"inverse sum poly = %f\n",poly);
    fprintf(stderr,"inverse scale factor = %f\n",sckin);

    for(ifa=0;ifa<13;++ifa) {
      for(jfa=0;jfa<13;++jfa) {
	rho = rhona[ifa][jfa];
	ral = rpxaml[ifa][jfa];
	if(fabs(ral) > 9000.0)
	  continue;
	for(ka=0;ka<nmat;++ka) {
	  for(ib=0;ib<nmat;++ib) {
	    dmt = pow(rhok,ka+ib);
	    // fprintf(stderr,"D %d %d    %f %f %f\n",ka,ib,rhok,dmt,ral);
	    dmat[ka][ib] += pow(ral*sckin,ka+ib);
	  }
	  dvec[ka][0] += (3600.0/dtor)*rho*pow(ral*sckin,ka);
	}
      } // for(jfa...
    } // for(ifa...
    for(ka=0;ka<nmat;++ka) {
      fprintf(stderr,"%f           ",dvec[ka][0]);
      for(ib=0;ib<nmat;++ib) {
	fprintf(stderr,"%f ",dmat[ka][ib]);
      }
      fprintf(stderr,"\n");
    }
    fflush(stderr);
    jg = gaussj(dmat,nmat,dvec,1);
    if(jg == 1) {
      // print message
      fprintf(stderr,"inverse polynomial code - iccd=%d; singular matrix\n",iccd);
    }
    else {
      // print dmat, dvec
      printf("inversePlateScalePoly\n");
      printf("originx = %e\n",0.0);
      printf("scalex = %e\n",sckin);
      printf("offsetx = %e\n",0.0);
      printf("power    poly_coeff   matrix_row\n");
      for(ka=0;ka<nmat;++ka) {
	printf("%d     %f        ",ka,dvec[ka][0]);
	fprintf(stderr,"%f           ",dvec[ka][0]);
	for(ib=0;ib<nmat;++ib) {
	  printf("%f ",dmat[ka][ib]);
	  fprintf(stderr,"%f ",dmat[ka][ib]);
	}
	printf("\n");
	fprintf(stderr,"\n");
      }
      printf("\n");
      fprintf(stderr,"\n");
    }

    // compute Ames pixel numbers and then Ames - AML pixel number differences
    nval = 0;
    rmssq = 0.0;
    dpmax = 0.0;
    for(ifa=0;ifa<nfa;++ifa) {
       for(jfa=0;jfa<nfa;++jfa) {
	 get_ames_model(nmat,cvec,sck,etana[ifa][jfa],xina[ifa][jfa],xyr);
	 if( (iccd==0) || (iccd==3) ) 
	   pxna[ifa][jfa] = xyr[0];
	 else
	   pxna[ifa][jfa] = -xyr[0];
	 if( (iccd==0) || (iccd==1) ) 
	   pyna[ifa][jfa] = xyr[1];
	 else
	   pyna[ifa][jfa] = -xyr[1];
	 if(pxaml[ifa][jfa] != 9999.0) {
	   /* fprintf(stderr,"XX %d %d %f %f      %f %f\n",
		   ifa,jfa,pxaml[ifa][jfa],pyaml[ifa][jfa],pxna[ifa][jfa],pyna[ifa][jfa]);
	   */
	   dpx = pxaml[ifa][jfa] - pxna[ifa][jfa];
	   dpy = pyaml[ifa][jfa] - pyna[ifa][jfa];
	   dsq = dpx*dpx + dpy*dpy;
	   if(dsq > dpmax)
	     dpmax = dsq;
	   ++nval;
	   rmssq += dsq;
	 }
       }
    }
    rmssq = sqrt(rmssq/nval);
    dpmax = sqrt(dpmax);
    fprintf(stderr,"nval,rmssq,dpmax = %d %f %f\n",nval,rmssq,dpmax);

  }  // for(iccd...

}

/******************************************************************************/

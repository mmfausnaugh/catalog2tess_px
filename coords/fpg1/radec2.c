/* radec2.c
 * Alan M. Levine
 * August 9, 2018
 * Heritage: radec1.c
 *
 * Compute RA, Dec based on star pixel coordinates.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"

#define C_LIGHT   (3.0e5)            // km/s

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
double asymang[NCAM], asymfac[NCAM];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];

double fpxlo[4], fpxhi[4], fpylo[4], fpyhi[4];
double cpxlo[4], cpxhi[4], cpylo[4], cpyhi[4];
double col, row, ra, dec, ra_ab, dec_ab;
int icam;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3], rmat5[3][3];
double ra_sc, dec_sc, roll_sc;

double dtor;

/****************************************************************************/
// Get the part of a vector v[3] that is perpendicular to another vector s[3].
// This is needed for stellar aberration calculations

void get_perp(double s[3], double v[3], double vp[3])
{
  double sdotv, sdots;
  int i;

  sdotv = dot(s,v);
  sdots = dot(s,s);
  if (sdots > 0.0) {
    for(i=0;i<3;++i) {
      vp[i] = v[i] - (sdotv*s[i]/sdots);
    }
  }
  else {
    for(i=0;i<3;++i) {
      vp[i] = v[i];
    }
  }
}

/****************************************************************************/
// Compute apparent ra and dec taking velocity aberration into account.

// ra, dec, raap, decap should have units of radians
// v[3] is the velocity vector; it should be given in km/s

void aber_star_pos(double v[3], double ra, double dec, double *raap,
                   double *decap)
{
  double s[3], sap[3], beta[3], bp[3];
  int i;

  for(i=0;i<3;++i) {
    beta[i] = v[i]/C_LIGHT;
  }
  dsphcr(ra,dec,s);
  get_perp(s,beta,bp);
  for(i=0;i<3;++i) {
    sap[i] = s[i] + bp[i];
  }
  dnorm(sap);
  dcrsph(sap,raap,decap);
}

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
/*  TSIG/AML optics model
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
// cpx, cpy limits are expressed in the ccdpx system
// fpx, fpy limits are expressed in the fitpx system

void ccd_pix_limits()
{
  int iccd;

  fpxlo[0] = CCDWD_T + 2*ROWA + ROWB - 0.5;
  fpxlo[1] = ROWA - 0.5;
  fpxlo[2] = ROWA - 0.5;
  fpxlo[3] = CCDWD_T + 2*ROWA + ROWB - 0.5;
  fpxhi[0] = fpxlo[0] + CCDWD_T;
  fpxhi[1] = fpxlo[1] + CCDWD_T;
  fpxhi[2] = fpxlo[2] + CCDWD_T;
  fpxhi[3] = fpxlo[3] + CCDWD_T;
  fpylo[0] = CCDHT_T + 2*COLDK_T - 0.5;
  fpylo[1] = CCDHT_T + 2*COLDK_T - 0.5;
  fpylo[2] = -0.5;
  fpylo[3] = -0.5;
  fpyhi[0] = fpylo[0] + CCDHT_T;
  fpyhi[1] = fpylo[1] + CCDHT_T;
  fpyhi[2] = fpylo[2] + CCDHT_T;
  fpyhi[3] = fpylo[3] + CCDHT_T;
  
  for(iccd=0;iccd<4;++iccd) {
    fprintf(stderr,"%d %f %f %f %f\n",
            iccd,fpxlo[iccd],fpxhi[iccd],fpylo[iccd],fpyhi[iccd]);
  }
  fprintf(stderr,"\n");

  for(iccd=0;iccd<4;++iccd) {
    cpxlo[iccd] = -0.5;
    cpxhi[iccd] = cpxlo[iccd] + CCDWD_T;
    cpylo[iccd] = -0.5;
    cpyhi[iccd] = cpylo[iccd] + CCDHT_T;
  }

  for(iccd=0;iccd<4;++iccd) {
    fprintf(stderr,"%d %f %f %f %f\n",
            iccd,cpxlo[iccd],cpxhi[iccd],cpylo[iccd],cpyhi[iccd]);
  }
  fprintf(stderr,"\n");
}

/******************************************************************************/
// return value = 0 => not on a ccd
//              = 1 => on a ccd

int ccdpx_on_ccd(double ccdpx[2])
{
  int iccd, jccd, jon;

  jon = 0;
  for(jccd=0;jccd<4;++jccd) {
    if( (ccdpx[0] >= cpxlo[jccd]) && (ccdpx[0] <= cpxhi[jccd])
        && (ccdpx[1] >= cpylo[jccd]) && (ccdpx[1] <= cpyhi[jccd])  ) { 
      jon = 1;
      break;
    }
  }
  return(jon);
}

/******************************************************************************/

int fitpx_ccd_num(double fitpx[2])
{
  int iccd, jccd;

  iccd = -1;
  for(jccd=0;jccd<4;++jccd) {
    if( (fitpx[0] >= fpxlo[jccd]) && (fitpx[0] <= fpxhi[jccd])
        && (fitpx[1] >= fpylo[jccd]) && (fitpx[1] <= fpyhi[jccd])  ) { 
      iccd = jccd;
      break;
    }
  }
  return(iccd);
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
  if(ccdpx_on_ccd(ccdpx) != 1)
    iccd = -1;

  return(iccd);
}

/******************************************************************************/
/* Convert pseudo-equivalent TESS camera image FITS file pixel numbers to
 * focal plane coordinates in pixel .
 */

// Cam num = icam + 1
// CCD num = iccd + 1

// July 1, 2018 - converted to flight software pixel numbering convention

int pix_to_mm(double xy[2], double fitpx[2])
{
  int iccd, jccd;
  double ccdzpx[2];
  double xyb[2], xyccd[2];

  iccd = -1;
  xy[0] = -999.0;
  xy[1] = -999.0;

  iccd = fitpx_ccd_num(fitpx);

  if( (iccd == 0) || (iccd == 1) ) {
    ccdzpx[0] = CCDWD_T + fpxlo[iccd] - fitpx[0];
    ccdzpx[1] = CCDHT_T + fpylo[iccd] - fitpx[1];
  }
  else if( (iccd == 2) || (iccd == 3) ){
    ccdzpx[0] = fitpx[0] - fpxlo[iccd];
    ccdzpx[1] = fitpx[1] - fpylo[iccd];
  }

  if(iccd >= 0) {
    xyccd[0] = ccdzpx[0]*pixsz[icam][iccd][0];
    xyccd[1] = ccdzpx[1]*pixsz[icam][iccd][1];
    xyrotate(-ccdang[icam][iccd],xyccd,xyb);
    xy[0] = xyb[0] + ccdxy0[icam][iccd][0];
    xy[1] = xyb[1] + ccdxy0[icam][iccd][1];
  }

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
// DON'T FORGET - Add effects due to stellar aberration later

int main(int argc, char *argv[])
{
  char *pixbase;
  double radecroll[3], eulc[3], vstar[3], vcam[3], ras, decs, radg, decdg;
  double rdr[3], racor, deccor, racdg, deccdg;
  double xyfp[2], ccdpx[2], fitpx[2], lnglat_deg[2], vel[3];
  int icam, iccd, i, jccd, kccd, j;
  FILE *fpf;

  dtor = M_PI/180.0;
  ccd_pix_limits();

  fpf = fopen("fpg_pars.txt","r");
  icam = 0;
  read_fpg_pars(fpf,stderr);
  fclose(fpf);

  if(argc != 8) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  ra_sc = atof(argv[1]);
  radecroll[0] = ra_sc;
  dec_sc = atof(argv[2]);
  radecroll[1] = dec_sc;
  roll_sc = atof(argv[3]);
  radecroll[2] = roll_sc;
  fprintf(stderr,"ra,dec,roll S/C = %f %f %f\n",ra_sc,dec_sc,roll_sc);
  pixbase = argv[4];
  vel[0] = -1.0*atof(argv[5]);
  vel[1] = -1.0*atof(argv[6]);
  vel[2] = -1.0*atof(argv[7]);

  sky_to_sc_mat(radecroll);
  prmat(rmat1,stderr);
  for(j=0;j<3;++j) {
    eulc[j] = eulcam[icam][j];
  }
  sc_to_cam_mat(eulc);
  matmat(rmat2,rmat1,rmat4);
  prmat(rmat4,stderr);
  /* get_ra_dec_roll(rmat4,rdr);
  fprintf(stderr,"cam ra,dec,roll = %f %f %f\n",
	  rdr[0],rdr[1],rdr[2]);
  */
  trans(rmat4,rmat5);

  while(read_pixel_coords(stdin,stderr) == 2) {
    if(strcmp(pixbase,"ZERO_BASE") == 0) {
      fitpx[0] = col;
      fitpx[1] = row;
    }
    else if(strcmp(pixbase,"ONE_BASE") == 0) {
      fitpx[0] = col - 1.0;
      fitpx[1] = row - 1.0;
    }
    else {
      fprintf(stderr,"ERROR: pixbase = %s; exiting...\n",pixbase);
      exit(-1);
    } 
    iccd = pix_to_mm(xyfp,fitpx);
    fp_optics(lnglat_deg,xyfp);
    dsphcr(dtor*lnglat_deg[0],dtor*lnglat_deg[1],vcam);
    matvec(rmat5,vcam,vstar);
    dcrsph(vstar,&ras,&decs);
    radg = ras/dtor;
    decdg = decs/dtor;
    aber_star_pos(vel,ras,decs,&racor,&deccor);
    racdg = racor/dtor;
    deccdg = deccor/dtor;
    fprintf(stdout,"col,row= %11.4f %11.4f aber_pos= %11.6f  %11.6f unaber_pos= %11.6f %11.6f\n",
	    col,row,radg,decdg,racdg,deccdg);
  }
}

/******************************************************************************/

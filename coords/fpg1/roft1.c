/* roft1.c
 * Alan M. Levine
 * July 13, 2018
 * Heritage: radldst1.c
 *
 * Compute radial distances in the focal plane that correspond to
 * specified angles from the boresight.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "vec.h"
//#include "mat_ra3.h"

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

#define NTHETA 4

double tharr[NTHETA] = { 2.0, 6.0, 9.0, 12.0 };  // degrees

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
double asymang[NCAM], asymfac[NCAM];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];

double fpxlo[4], fpxhi[4], fpylo[4], fpyhi[4];
double cpxlo[4], cpxhi[4], cpylo[4], cpyhi[4];
int icam;
double dtor;

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
  /* ttsq = tanth*tanth;
  rfp = FLMM*tanth*(1.00000140 + (0.28174612*ttsq) + (-0.59667259*ttsq*ttsq)
        + (9.17151267*ttsq*ttsq*ttsq) + (-4.36928235*ttsq*ttsq*ttsq*ttsq));
  */
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

void usage(char *prog)
{
  fprintf(stderr,"usage: %s <fpg_pars.txt file> > output_file\n",
          prog);
  exit(1);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  char *filea;
  int i;
  double thetadeg, tanth, r;
  FILE *fp;

  dtor = M_PI/180.0;

  if(argc != 2) {
    usage(argv[0]);
    exit(-1);
  }
  filea = argv[1];  // fpg_pars.txt file
  fprintf(stderr,"fpg_pars.txt file = %s\n",filea);

  icam = 0;
  fp = fopen(filea,"r");
  read_fpg_pars(fp,stderr);
  fclose(fp);

  for(i=0;i<NTHETA;++i) {
    thetadeg = tharr[i];
    tanth = tan(dtor*thetadeg);
    r = r_of_tanth(tanth);
    printf("%10.5f ",r);
  }
  printf("\n");

}

/******************************************************************************/

/* starspx8.c
 * Alan M. Levine
 * June 23, 2018
 * Heritage: starspx7.c
 *
 * Compute predicted star positions in the focal plane and pixel coordinates.
 */

// August 3, 2017
// pixel numbers in fitpx[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention
// This addition of one is done below for consistency with ground FITS files.

// March 6, 2018 - Add columns in the output for sigx and sigy.

// March 19, 2018 - Add command line parameter to choose zero or one base for
//                  pixel numbering.
//                  zero base - numbering I use for pixels in a FITS image
//                  one base - numbering used for pixels in a ds9 display

// June 16, 2018 - Widen field definition in star_in_fov().

// June 23, 2018 - Read and process one star at a time, so there is no limit
//                 on the number of stars to be processed.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vec.h"
#include "mat_ra3.h"

// #define NSTAR 10000

#define CCDWD_T 2048
#define CCDHT_T 2058
#define ROWA 44
#define ROWB 44
#define COLDK_T 20

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

#define SIGX 1.0    // placeholder
#define SIGY 1.0    // ditto

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
double asymang[NCAM], asymfac[NCAM];

int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];

double row, col, tmag, ra, dec;
double ra_ab, dec_ab;
int nst, icam;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3];
double ra_sc, dec_sc, roll_sc;

double dtor;
void prmat();

/******************************************************************************/
// All input angles are in degrees.

void read_fpg_pars(FILE *fpin, FILE *fpo, int iicam)
{
  char pardesc[100];
  int nrd, iccd, j;

  icam = iicam;   // save for later use
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
    // 12.0 -> 12.2 in the following 2 lines (June 16, 2018)
    if( (fabs(atan(v[0]/v[2])) <= (12.2*dtor)) &&
	(fabs(atan(v[1]/v[2])) <= (12.2*dtor)) ) {
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
// return value = 1 if a star is read successfully, = 0 otherwise
// ipr = 1 => print results

int read_star(FILE *fpin, FILE *fpd, int ipr)
{
  int nrd;

  nrd = fscanf(fpin,"%lf %lf %lf %lf %lf",
	       &ra,&dec,&ra_ab,&dec_ab,&tmag);
  if(nrd == 5) {
    if(ipr == 1)
      fprintf(fpd,"%f %f %f %f %f\n",ra,dec,ra_ab,dec_ab,tmag);
    return(1);
  }
  else {
    return(0);
  }
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

void usage(char *prog, int argc)
{
  fprintf(stderr,"Usage: %s camera_num(1-4) RA_sc Dec_sc Roll_sc base_string\n",prog);
  fprintf(stderr,"  < star_list_5_columns > output_file\n");
  fprintf(stderr,"  where base_string = ZERO_BASE or ONE_BASE for FFI image pixel numbers\n");
  exit(-1);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  char *pixbase;
  double radecroll[3], eulc[3], vstar[3], vcam[3], lng, lat, lngd, latd;
  double rdr[3];
  double xyfp[2], ccdpx[2], fitpx[2];
  int iccd, i, jccd, kccd, j;
  FILE *fpf;

  dtor = M_PI/180.0;

  if(argc != 6) {
    fprintf(stderr,"ERROR: argc = %d;\n",argc);
    usage(argv[0],argc);
  }
  icam = atoi(argv[1]) - 1;
  ra_sc = atof(argv[2]);
  radecroll[0] = ra_sc;
  dec_sc = atof(argv[3]);
  radecroll[1] = dec_sc;
  roll_sc = atof(argv[4]);
  radecroll[2] = roll_sc;
  fprintf(stderr,"ra,dec,roll S/C = %f %f %f\n",ra_sc,dec_sc,roll_sc);
  pixbase = argv[5];

  fpf = fopen("fpg_pars.txt","r");
  read_fpg_pars(fpf,stderr,icam);

  sky_to_sc_mat(radecroll);
  prmat(rmat1,stderr);
  for(j=0;j<3;++j) {
    eulc[j] = eulcam[icam][j];
  }
  sc_to_cam_mat(eulc);
  matmat(rmat2,rmat1,rmat4);
  prmat(rmat4,stderr);
  get_ra_dec_roll(rmat4,rdr);
  fprintf(stderr,"cam ra,dec,roll = %f %f %f\n",
	  rdr[0],rdr[1],rdr[2]);

  i = 0;
  while(read_star(stdin,stderr,0)==1) {
    dsphcr(dtor*ra_ab,dtor*dec_ab,vstar);
    matvec(rmat4,vstar,vcam);
    dcrsph(vcam,&lng,&lat);
    lngd = lng/dtor;
    latd = lat/dtor;
    if(i<10) {
      fprintf(stderr,"i,vstar = %d  %f %f %f\n",i,vstar[0],vstar[1],vstar[2]);
      fprintf(stderr,"i,vcam = %d  %f %f %f\n",i,vcam[0],vcam[1],vcam[2]);
      fprintf(stderr,"i,lngd,latd = %d %f %f\n",i,lngd,latd);
    }
    if(star_in_fov(lngd,latd) == 1) {
      optics_fp(lngd,latd,xyfp);
      iccd = mm_to_pix(xyfp,ccdpx,fitpx);

      if(strcmp(pixbase,"ONE_BASE") == 0) {
	fprintf(stdout,"%4d %11.5f %11.5f %11.5f %11.5f    %11.5f %2d %11.4f %11.4f   %11.4f %11.4f\n",
		i,ra,dec,ra_ab,dec_ab,tmag,iccd+1,fitpx[0]+1.0,fitpx[1]+1.0,
		SIGX,SIGY);
      }
      else if(strcmp(pixbase,"ZERO_BASE") == 0) {
	fprintf(stdout,"%4d %11.5f %11.5f %11.5f %11.5f    %11.5f %2d %11.4f %11.4f   %11.4f %11.4f\n",
		i,ra,dec,ra_ab,dec_ab,tmag,iccd+1,fitpx[0],fitpx[1],
		SIGX,SIGY);
      }
      else {
	usage(argv[0],argc);
      }
    }
    ++i;
  }
  fprintf(stderr,"%d stars were processed.\n",i);
}
/******************************************************************************/

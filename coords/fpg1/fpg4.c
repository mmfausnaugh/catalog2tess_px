/* fpg4.c
 * Alan M. Levine
 * August 5, 2017
 * Heritage: fpg3.c
 *
 * Compute predicted star positions in the focal plane and pixel coordinates.
 *
 * June 30, 2017 - Add capability to work with azimuthal asymmetry parameters.
 *
 * Aug. 3, 2017 - Start to check pixel numbering and centroid definitions.
 *
 * Aug. 5, 2017 - Do velocity aberration corrections prior to running this code.
 */

// WILL ONLY DO ONE CAMERA AT A TIME

// pixel numbers in fitpx[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"
#include "fpg4.h"
#include "mrq_fpg4.h"

#define NSTAR 10000

#define CCDWD_T 2048
#define CCDHT_T 2058
#define ROWA 44
#define ROWB 44
#define COLDK_T 20

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

double fpxlo[4], fpxhi[4], fpylo[4], fpyhi[4];

double eulcam0[NCAM][3], optcon0[NCAM][6], ccdxy00[NCAM][NCCD][2];
double pixsz0[NCAM][NCCD][2], ccdang0[NCAM][NCCD], ccdtilt0[NCAM][NCCD][2];
double asymang0[NCAM], asymfac0[NCAM];
double eulcam[NCAM][3], optcon[NCAM][6], ccdxy0[NCAM][NCCD][2];
double pixsz[NCAM][NCCD][2], ccdang[NCAM][NCCD], ccdtilt[NCAM][NCCD][2];
double asymang[NCAM], asymfac[NCAM];
double deleul[NCAM][3], delopt[NCAM][6], delxy0[NCAM][NCCD][2];
double delsz[NCAM][NCCD][2], delang[NCAM][NCCD], deltilt[NCAM][NCCD][2];
double dasymang[NCAM], dasymfac[NCAM];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];
int npar, nparfit, icam;
int **iparcode;

double row[NSTAR], col[NSTAR], tmag[NSTAR], ra[NSTAR], dec[NSTAR];
double ra_cor[NSTAR], dec_cor[NSTAR], sigx[NSTAR], sigy[NSTAR];
int starno[NSTAR], ccdno[NSTAR], istfl[NSTAR];
int nst;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3];
double ra_sc, dec_sc, roll_sc, radecroll[3];

double dtor;

/******************************************************************************/
// Multiply a 3-vector in celestial coordinates by rmat1 to get a 3-vector
// in S/C coordinates.

void sky_to_sc_mat()
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
  // prmat(rmat2,stderr);
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

int star_in_border(int iccd, double fcol, double frow, double pxbor)
{
  int infl;

  /* fprintf(stderr,"star_in_border(): iccd,fcol,frow,pxbor = %d %f %f %f\n",
     iccd,fcol,frow,pxbor); */
  infl = 0;
  if( (fcol >= (fpxlo[iccd] + pxbor)) && (fcol <= (fpxhi[iccd] - pxbor))  &&
      (frow >= (fpylo[iccd] + pxbor)) && (frow <= (fpyhi[iccd] - pxbor)) ) {
    infl = 1;
  }
  return(infl);
}

/******************************************************************************/

int get_ccdno(int istar)
{
  return(ccdno[istar] - 1);
}

/******************************************************************************/

void setup_mats(int ipr)
{
  int j;
  double eulc[3], rdr[3];

  sky_to_sc_mat();
  // prmat(rmat1,stderr);
  for(j=0;j<3;++j) {
    eulc[j] = eulcam[icam][j];
  }
  sc_to_cam_mat(eulc);
  matmat(rmat2,rmat1,rmat4);
  // prmat(rmat4,stderr);
  get_ra_dec_roll(rmat4,rdr);
  if(ipr == 1)
    fprintf(stderr,"cam ra,dec,roll = %f %f %f\n",rdr[0],rdr[1],rdr[2]);
}

/******************************************************************************/

int star_pixels(int i, double ccdpx[2], double fitpx[2], int icode, int ideriv)
{
  double vstar[3], vcam[3], lng, lat, lngd, latd, xyfp[2];
  int iccd;

  // PRE-MAKE VELOCITY ABERRATION CORRECTION - done
  dsphcr(dtor*ra_cor[i],dtor*dec_cor[i],vstar);
  setup_mats(0);
  matvec(rmat4,vstar,vcam);
  dcrsph(vcam,&lng,&lat);
  lngd = lng/dtor;
  latd = lat/dtor;
  optics_fp(lngd,latd,xyfp);
  iccd = mm_to_pix(xyfp,ccdpx,fitpx);

  /* fprintf(stdout,"%4d %1d %1d %11.5f %11.5f   %11.5f %11.5f      %11.5f %11.5f  %9.5f %9.5f %9.3f %9.3f"
	  " %2d %9.3f %9.3f\n",
	  i,icode,ideriv,ra[i],dec[i],ra_cor[i],dec_cor[i],lngd,90.0-latd,xyfp[0],xyfp[1],
	  ccdpx[0],ccdpx[1],iccd+1,fitpx[0],fitpx[1]); 
   */
  // fprintf(stderr,"star_pixels(): iccd = %d\n",iccd);
  return(iccd);
}

/******************************************************************************/
// For testing and other uses
// Not for use in Levenberg-Marquardt fitting

void all_star_pixels()
{
  double ccdpx[2], fitpx[2];
  int i, iccd;

  setup_mats(0);
  for(i=0;i<nst;++i) {
    // if(i > 3) break;
    iccd = star_pixels(i,ccdpx,fitpx,-1,0);
  }
}

/******************************************************************************/
/* INTERFACE FROM mrqcof()
 *
 * Translate parameters from a[][] array to variable names used in this file.
 *
 * For a given star,
 *  compute predicted pixel coordinates
 *  compute partial derivatives
 */

// This will be done numerically.
/* For each parameter, a small change in the parameter needs to be given.
 * Each 'delta' should produce approximately the same change in the pixel
 * coordinate, in most cases.
 *
 * Both d(x_pixel)/da and d(y_pixel)/da are needed.
 *
 * Use dpixel/da =  (pixel(a+da) - pixel(a-da))/(2*da)
 */

// PRECORRECT ra's and dec's for aberration

int one_star_fit_inputs(int istar, double fitpx[2])
{
  double ccdpx[2], cpxp[2], fpxp[2], cpxm[2], fpxm[2], delta, temp;
  int icc[2], iccd, jccd, kccd, mccd;
  int ipf, jp, idec[4];

  // fprintf(stderr,"one_star_fit_inputs(): beginning...\n");

  if(istfl[istar]==0)
    return(-2);

  // fprintf(stderr,"one_star_fit_inputs(): istar= %d\n",istar);
  iccd = star_pixels(istar,ccdpx,fitpx,-1,0);

  ipf = 0;
  for(jp=0;jp<npar;++jp) {
    if(lista[jp][0] == 1) {
      decode_param(iparcode[jp][0],idec);
      /* fprintf(stderr,"jp= %d, ,iparcode,idec[0],[1],[2],[3] = %d   %d %d %d %d\n",
	 jp,iparcode[jp][0],idec[0],idec[1],idec[2],idec[3]); */
      switch(idec[0]) {
      case 0:
	temp = eulcam[idec[1]][idec[3]];
	delta = deleul[idec[1]][idec[3]];
	eulcam[idec[1]][idec[3]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,0,1);
	eulcam[idec[1]][idec[3]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,0,2);
	eulcam[idec[1]][idec[3]] = temp;
	break;
      case 1:
	temp = optcon[idec[1]][idec[3]];
	delta = delopt[idec[1]][idec[3]];
	optcon[idec[1]][idec[3]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,1,1);
	optcon[idec[1]][idec[3]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,1,2);
	optcon[idec[1]][idec[3]] = temp;
	break;
      case 2:
	temp = asymang[idec[1]];
	delta = dasymang[idec[1]];
	asymang[idec[1]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,2,1);
	asymang[idec[1]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,2,2);
	asymang[idec[1]] = temp;
	break;
      case 3:
	temp = asymfac[idec[1]];
	delta = dasymfac[idec[1]];
	asymfac[idec[1]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,3,1);
	asymfac[idec[1]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,3,2);
	asymfac[idec[1]] = temp;
	break;
      case 4:
	temp = ccdxy0[idec[1]][idec[2]][idec[3]];
	delta = delxy0[idec[1]][idec[2]][idec[3]];
	ccdxy0[idec[1]][idec[2]][idec[3]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,4,1);
	ccdxy0[idec[1]][idec[2]][idec[3]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,4,2);
	ccdxy0[idec[1]][idec[2]][idec[3]] = temp;
	break;
      case 5:
	temp = pixsz[idec[1]][idec[2]][idec[3]];
	delta = delsz[idec[1]][idec[2]][idec[3]];
	pixsz[idec[1]][idec[2]][idec[3]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,5,1);
	pixsz[idec[1]][idec[2]][idec[3]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,5,2);
	pixsz[idec[1]][idec[2]][idec[3]] = temp;
	break;
      case 6:
	temp = ccdang[idec[1]][idec[2]];
	delta = delang[idec[1]][idec[2]];
	ccdang[idec[1]][idec[2]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,6,1);
	ccdang[idec[1]][idec[2]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,6,2);
	ccdang[idec[1]][idec[2]] = temp;
	break;
      case 7:
        temp = ccdtilt[idec[1]][idec[2]][idec[3]];
	delta = deltilt[idec[1]][idec[2]][idec[3]];
	ccdtilt[idec[1]][idec[2]][idec[3]] = temp + delta;
	jccd = star_pixels(istar,cpxp,fpxp,7,1);
	ccdtilt[idec[1]][idec[2]][idec[3]] = temp - delta;
	kccd = star_pixels(istar,cpxm,fpxm,7,2);
	ccdtilt[idec[1]][idec[2]][idec[3]] = temp;
	break;
      default:
	// ERROR - print diagnostic and exit
	exit(-2);
      }
      if( (jccd != iccd) || (kccd != iccd) ) {
	mccd = -1;
	break;  // don't use this star
      }
      else {
	mccd = iccd;
      }
      if(fabs(delta) > 0.0) {
	dxda[jp][0] = (fpxp[0] - fpxm[0])/(2.0*delta);
	dyda[jp][0] = (fpxp[1] - fpxm[1])/(2.0*delta);
      }
      else {
	dxda[jp][0] = 1.0e-12;
	dyda[jp][0] = 1.0e-12;
      }  
      /* fprintf(stderr,"jp,delta,dxda,dyda= %d %f %f %f\n",
	 jp,delta,dxda[jp][0],dyda[jp][0]);  */
    } // if(lista == 1)
    ++ipf;
  }  // for(jp=0 ...)
  return(mccd);
}

/******************************************************************************/
/* I need a way to go from index number j of a parameter in a[][0] to the type
 * of parameter, camera number, CCD number, and component number.
 *
 * Type   0-5
 * Camera 0-3
 * CCD    0-3
 * num    0-5
 */

int code_param(int itype, int iicam, int iccd, int inum)
{
  int icode;

  icode = itype +(10*iicam) + (100*iccd) + (1000*inum);
  return(icode);
}
/******************************************************************************/

void decode_param(int icode, int idec[4])
{
  idec[0] = icode % 10;         // itype
  idec[1] = (icode/10) % 10;    // iicam
  idec[2] = (icode/100) % 10;   // iccd
  idec[3] = (icode/1000) % 10;  // inum
}
/******************************************************************************/

void fill_fpg_param_arrays(int iarr)
{
  int ipar, idec[4];

  fprintf(stderr,"fill_fpg_param_arrays(): beginning...\n");

  for(ipar=0;ipar<npar;++ipar) {
    decode_param(iparcode[ipar][0],idec);
    switch(idec[0]) {
    case 0:
      if(iarr == 0)
	eulcam[idec[1]][idec[3]] = a[ipar][0];
      else
	eulcam[idec[1]][idec[3]] = atry[ipar][0];
      break;
    case 1:
      if(iarr == 0)
	optcon[idec[1]][idec[3]] = a[ipar][0];
      else
	optcon[idec[1]][idec[3]] = atry[ipar][0];
      break;
    case 2:
      if(iarr == 0)
	asymang[idec[1]] = a[ipar][0];
      else
	asymang[idec[1]] = atry[ipar][0];
      break;
    case 3:
      if(iarr == 0)
	asymfac[idec[1]] = a[ipar][0];
      else
	asymfac[idec[1]] = atry[ipar][0];
      break;
    case 4:
      if(iarr == 0)
	ccdxy0[idec[1]][idec[2]][idec[3]] = a[ipar][0];
      else
	ccdxy0[idec[1]][idec[2]][idec[3]] = atry[ipar][0];
      break;
    case 5:
      if(iarr == 0)
	pixsz[idec[1]][idec[2]][idec[3]] = a[ipar][0];
      else
	pixsz[idec[1]][idec[2]][idec[3]] = atry[ipar][0];
      break;
    case 6:
      if(iarr == 0)
	ccdang[idec[1]][idec[2]] = a[ipar][0];
      else
	ccdang[idec[1]][idec[2]] = atry[ipar][0];
      break;
    case 7:
      if(iarr == 0)
	ccdtilt[idec[1]][idec[2]][idec[3]] = a[ipar][0];
      else
	ccdtilt[idec[1]][idec[2]][idec[3]] = atry[ipar][0];
      break;
    default:
      // ERROR - print diagnostic and exit
      fprintf(stderr,"fill_fpg_param_arrays(): decode error\n");
      fprintf(stderr,"    type = %d;  exiting ...\n",idec[0]);
      exit(-3);
    }
  }
}

/******************************************************************************/
/* This function should be called after read_fpg_pars() and read_stars().
 *
 * Put all parameters in a single array and fill lista[][].
 */

// double deleul[NCAM][3], delopt[NCAM][6], delxy0[NCAM][NCCD][2];
// double delsz[NCAM][NCCD][2], delang[NCAM][NCCD], deltilt[NCAM][NCCD][2];
// double dasymang[NCAM], dasymfac[NCAM];

void fill_mrq_param_arrays()
{
  int j, k, iccd;

  // fprintf(stderr,"fmpa A, npar=%d\n",npar);fflush(stderr);
  fprintf(stderr,"fill_mrq_param_arrays(): beginning...\n");

  mrq_define(nst,npar,nparfit);
  iparcode = imatrix(npar,1);

  k = 0;
  for(j=0;j<3;++j) {
    a[k][0] = eulcam[icam][j];
    damax[k][0] = fabs(deleul[icam][j]);
    lista[k][0] = keulcam[icam][j];
    iparcode[k][0] = code_param(0,icam,0,j);
   ++k;
  }
  for(j=0;j<6;++j) {
    a[k][0] = optcon[icam][j];
    damax[k][0] = fabs(delopt[icam][j]);
    lista[k][0] = koptcon[icam][j];
    iparcode[k][0] = code_param(1,icam,0,j);
    ++k;
  }
  a[k][0] = asymang[icam];
  damax[k][0] = fabs(dasymang[icam]);
  lista[k][0] = kasymang[icam];
  iparcode[k][0] = code_param(2,icam,0,j);
  ++k;
  a[k][0] = asymfac[icam];
  damax[k][0] = fabs(dasymfac[icam]);
  lista[k][0] = kasymfac[icam];
  iparcode[k][0] = code_param(3,icam,0,j);
  ++k;
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      a[k][0] = ccdxy0[icam][iccd][j];
      damax[k][0] = fabs(delxy0[icam][iccd][j]);
      lista[k][0] = kccdxy0[icam][iccd][j];
      iparcode[k][0] = code_param(4,icam,iccd,j);
      ++k;
    }
    for(j=0;j<2;++j) {
      a[k][0] = pixsz[icam][iccd][j];
      damax[k][0] = fabs(delsz[icam][iccd][j]);
      lista[k][0] = kpixsz[icam][iccd][j];
      iparcode[k][0] = code_param(5,icam,iccd,j);
      ++k;
    }
    a[k][0] = ccdang[icam][iccd];
    damax[k][0] = fabs(delang[icam][iccd]);
    lista[k][0] = kccdang[icam][iccd];
    iparcode[k][0] = code_param(6,icam,iccd,0);
    ++k;
    for(j=0;j<2;++j) {
      a[k][0] = ccdtilt[icam][iccd][j];
      damax[k][0] = fabs(deltilt[icam][iccd][j]);
      lista[k][0] = kccdtilt[icam][iccd][j];
      iparcode[k][0] = code_param(7,icam,iccd,j);
      ++k;
    }
  }
  set_indexa();
  print_index_lists();
}

/******************************************************************************/
// All input angles are in degrees.

void read_fpg_pars(FILE *fpin, FILE *fpo, int iicam)
{
  char pardesc[100];
  int nrd, iccd, j;

  fprintf(stderr,"read_fpg_pars(): beginning...\n");
  npar = 0;
  nparfit = 0;
  icam = iicam;   // save for later use
  for(j=0;j<3;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&keulcam[icam][j],&eulcam[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,keulcam[icam][j],eulcam[icam][j]);
    eulcam0[icam][j] = eulcam[icam][j];
    ++npar;
    if(keulcam[icam][j] == 1) ++nparfit;
  }
  for(j=0;j<6;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&koptcon[icam][j],&optcon[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,koptcon[icam][j],optcon[icam][j]);
    optcon0[icam][j] = optcon[icam][j];
    ++npar;
    if(koptcon[icam][j] == 1) ++nparfit;
  }
  // double asymang[NCAM][NCCD], asymfac[NCAM][NCCD];
  fscanf(fpin,"%s %d %lf",pardesc,&kasymang[icam],&asymang[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymang[icam],asymang[icam]);
  asymang0[icam] = asymang[icam];
  ++npar;
  if(kasymang[icam] == 1) ++nparfit;
  fscanf(fpin,"%s %d %lf",pardesc,&kasymfac[icam],&asymfac[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymfac[icam],asymfac[icam]);
  asymfac0[icam] = asymfac[icam];
  ++npar;
  if(kasymfac[icam] == 1) ++nparfit;
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdxy0[icam][iccd][j],&ccdxy0[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
      ccdxy00[icam][iccd][j] = ccdxy0[icam][iccd][j];
      ++npar;
      if(kccdxy0[icam][iccd][j] == 1) ++nparfit;
    }
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kpixsz[icam][iccd][j],&pixsz[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
      pixsz0[icam][iccd][j] = pixsz[icam][iccd][j];
      ++npar;
      if(kpixsz[icam][iccd][j] == 1) ++nparfit;
    }
    fscanf(fpin,"%s %d %lf",pardesc,&kccdang[icam][iccd],&ccdang[icam][iccd]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,kccdang[icam][iccd],ccdang[icam][iccd]);
    ccdang0[icam][iccd] = ccdang[icam][iccd];
    ++npar;
    if(kccdang[icam][iccd] == 1) ++nparfit;
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdtilt[icam][iccd][j],&ccdtilt[icam][iccd][j]);
      fprintf(fpo,"%s       %d   %f\n",pardesc,kccdtilt[icam][iccd][j],ccdtilt[icam][iccd][j]);
      ccdtilt0[icam][iccd][j] = ccdtilt[icam][iccd][j];
      ++npar;
      if(kccdtilt[icam][iccd][j] == 1) ++nparfit;
    }
  }
}

/******************************************************************************/
// All angles are in degrees.

void print_fpg_pars_deriv(FILE *fpo, int idoder)
{
  char pardesc[100];
  int iccd, j, jcam;

  jcam = icam + 1;
  fprintf(stderr,"print_fpg_pars_deriv(): beginning...\n");
  for(j=0;j<3;++j) {
    sprintf(pardesc,"ang%1d_cam%1d",j+1,jcam);
    fprintf(fpo,"%-18s        %d  %f",
	      pardesc,keulcam[icam][j],eulcam[icam][j]);
    if(idoder == 1)
      fprintf(fpo," %f\n",deleul[icam][j]);
    else
      fprintf(fpo,"\n");
 }
  sprintf(pardesc,"fl_cam%1d",jcam);
  fprintf(fpo,"%-18s        %d  %f",
	  pardesc,koptcon[icam][0],optcon[icam][0]);
  if(idoder == 1)
    fprintf(fpo," %f\n",delopt[icam][0]);
  else
    fprintf(fpo,"\n");
  for(j=1;j<6;++j) {
    sprintf(pardesc,"opt_coef%1d_cam%1d",j,jcam);
    fprintf(fpo,"%-18s        %d  %11.8f",
	    pardesc,koptcon[icam][j],optcon[icam][j]);
    if(idoder == 1)
      fprintf(fpo," %f\n",delopt[icam][j]);
    else
      fprintf(fpo,"\n");
  }
  sprintf(pardesc,"asymang_cam%1d",jcam);
  fprintf(fpo,"%-18s        %d  %11.8f",
	    pardesc,kasymang[icam],asymang[icam]);
  if(idoder == 1)
    fprintf(fpo," %f\n",dasymang[icam]);
  else
    fprintf(fpo,"\n");
  sprintf(pardesc,"asymfac_cam%1d",jcam);
  fprintf(fpo,"%-18s        %d  %11.8f",
	    pardesc,kasymfac[icam],asymfac[icam]);
  if(idoder == 1)
    fprintf(fpo," %f\n",dasymfac[icam]);
  else
    fprintf(fpo,"\n");
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      if(j==0)
	sprintf(pardesc,"x0_ccd%1d_cam%1d",iccd+1,jcam);
      else
	sprintf(pardesc,"y0_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s        %d  %f",
	      pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
      if(idoder == 1)
	fprintf(fpo," %f\n",delxy0[icam][iccd][j]);
      else
	fprintf(fpo,"\n");
    }
    for(j=0;j<2;++j) {
      if(j==0)
	sprintf(pardesc,"pix_x_ccd%1d_cam%1d",iccd+1,jcam);
      else
	sprintf(pardesc,"pix_y_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s        %d  %f",
	      pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
      if(idoder == 1)
	fprintf(fpo," %f\n",delsz[icam][iccd][j]);
      else
	fprintf(fpo,"\n");
    }
    sprintf(pardesc,"ang_ccd%1d_cam%1d",iccd+1,jcam);
    fprintf(fpo,"%-18s        %d  %f",
	    pardesc,kccdang[icam][iccd],ccdang[icam][iccd]);
    if(idoder == 1)
      fprintf(fpo," %f\n",delang[icam][iccd]);
    else
      fprintf(fpo,"\n");
    for(j=0;j<2;++j) {
      if(j==0)
	sprintf(pardesc,"tilt_x_ccd%1d_cam%1d",iccd+1,jcam);
      else
	sprintf(pardesc,"tilt_y_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s        %d  %f",
	      pardesc,kccdtilt[icam][iccd][j],ccdtilt[icam][iccd][j]);
      if(idoder == 1)
	fprintf(fpo," %f\n",deltilt[icam][iccd][j]);
      else
	fprintf(fpo,"\n");
   }
  }
}

/******************************************************************************/
// Set deltas for use in calculating dx/dpar and dy/dpar
// Do this after parameter values have been read in.

// 'dxy_mm' is a small distance of order the size of a pixel (in mm).
// 'npixccd' = 2048 - approximately the no. of pixels in a row or column

void delta_for_deriv(double dxy_mm, int npixccd)
{
  int j, k, iccd;
  double th_fov;

  // fprintf(stderr,"delta_for_deriv(): beginning...\n");
  th_fov = 0.25;   // rough value in radians (center to corner of fov)
  deleul[icam][0] = (dxy_mm/optcon[icam][0])/dtor;
  deleul[icam][1] = (dxy_mm/optcon[icam][0])/dtor;
  deleul[icam][2] = (1.0/th_fov)*(dxy_mm/optcon[icam][0])/dtor;
  delopt[icam][0] = (1.0/th_fov)*dxy_mm;
  for(j=1;j<6;++j) {
    delopt[icam][j] = (4.0*j*j)*(dxy_mm/optcon[icam][0])*pow(th_fov,-(1.0*j-1.0));
  }

  // CHECK NEXT TWO LINES
  // arbitrary choice - based on expectation that asymfac will be very close to one
  dasymang[icam] =  10.0;
  dasymfac[icam] =  1.0*dxy_mm/(1.0*npixccd*pixsz[icam][0][0]);
  /* fprintf(stderr,"delta_for_deriv(): dxy_mm,dasymang,dasymfac,np,sz = %f %f %f %d %f\n",
     dxy_mm,dasymang[icam],dasymfac[icam],npixccd,pixsz[icam][0][0]); */
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      delxy0[icam][iccd][j] = dxy_mm;
      delsz[icam][iccd][j] = 1.0*dxy_mm/npixccd;
      deltilt[icam][iccd][j] = (dxy_mm/(npixccd*pixsz[icam][iccd][j]))/dtor;
    }
    delang[icam][iccd] = (dxy_mm/(2.0*npixccd*pixsz[icam][iccd][0]))/dtor;
  }
}

/******************************************************************************/
// Read RA, Dec, TESS magnitude, and color (to be inserted later) for each star

void read_fpg_stars(FILE *fpin, FILE *fpd, double pxbor)
{
  int i;
  double raapr, decapr;

  fprintf(stderr,"read_fpg_stars(): beginning...\n");
  set_fpx_lims();
  i = 0;
  while(fscanf(fpin,"%d %lf %lf %lf %lf %lf",
	       &starno[i],&ra[i],&dec[i],&ra_cor[i],&dec_cor[i],&tmag[i]) == 6) {
    fscanf(fpin,"%d %lf %lf %lf %lf",
	   &ccdno[i],&col[i],&row[i],&sigx[i],&sigy[i]);
    istfl[i] = star_in_border(ccdno[i]-1,col[i],row[i],pxbor);
    fprintf(fpd,"%d %f %f    %f %f      %f     %f %f    %f %f      %d\n",
	    i,ra[i],dec[i],ra_cor[i],dec_cor[i],tmag[i],col[i],row[i],sigx[i],
	    sigy[i],istfl[i]);
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
	  ccdno[i],col[i],row[i],sigx[i],sigy[i]);
  fprintf(fpo,"   %9.3f %9.3f %9.3f %9.3f\n",
	  xmod[i][0],ymod[i][0],col[i]-xmod[i][0],row[i]-ymod[i][0]);
}

/******************************************************************************/

void print_fpg_all_stars(FILE *fpo)
{
  int i;

  for(i=0;i<nst;++i) {
    print_fpg_star(i,fpo);
  }
}

/******************************************************************************/

void fill_mrq_stars()
{
  int i;

  fprintf(stderr,"fill_mrq_stars(): beginning...\n");
  for(i=0;i<nst;++i) {
    x[i][0] = col[i];
    y[i][0] = row[i];
    xsig[i][0] = sigx[i];
    ysig[i][0] = sigy[i];
  }
}

/******************************************************************************/

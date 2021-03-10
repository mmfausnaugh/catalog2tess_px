/* fpgtr1.c
 * Alan M. Levine
 * February 9, 2018
 * Heritage: runfpg5.c
 *
 * Perform a coordinate transform of geometric parameters of a TESS camera.
 * The result may be exact or approximate according to the desired transformation.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"
#include "fpg4.h"

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

// input set
extern double eulcam[NCAM][3], optcon[NCAM][6], ccdxy0[NCAM][NCCD][2];
extern double pixsz[NCAM][NCCD][2], ccdang[NCAM][NCCD], ccdtilt[NCAM][NCCD][2];
extern double asymang[NCAM], asymfac[NCAM];

// A set
double aeulcam[NCAM][3], aoptcon[NCAM][6], accdxy0[NCAM][NCCD][2];
double apixsz[NCAM][NCCD][2], accdang[NCAM][NCCD], accdtilt[NCAM][NCCD][2];
double aasymang[NCAM], aasymfac[NCAM];

// B set
double beulcam[NCAM][3], boptcon[NCAM][6], bccdxy0[NCAM][NCCD][2];
double bpixsz[NCAM][NCCD][2], bccdang[NCAM][NCCD], bccdtilt[NCAM][NCCD][2];
double basymang[NCAM], basymfac[NCAM];

int icam;
double rm;

/******************************************************************************/
// Input angles are in degrees.
/* Rotate camera frame according to a set of 1-2-3 Euler angles.
 *
 * 1) Properly adjust ang1, ang2, and ang3.
 * 2) For each of the four CCDs, adjust x0, y0, and angle.
 */

void rotate_cam_frame(double eultr[3])
{
  double eulcr[3], eultrr[3], eult[3];
  double rm1[3][3], rm2[3][3], rm3[3][3];
  int i;

  for(i=0;i<3;++i) {
    eulcr[i] = dtor*eulcam[icam][i];
    eultrr[i] = dtor*eultr[i];
  }

  eulerm323(eulcr,rm1);
  prmat(rm1,stderr);
  eulerm123(eultrr,rm2);
  prmat(rm2,stderr);
  matmat(rm2,rm1,rm3);
  prmat(rm3,stderr);
  mateuler323(rm3,eult);
  for(i=0;i<3;++i) {
    eulcam[icam][i] = eult[i]/dtor;
  }
  fprintf(stderr,"New Euler angles(3-2-3)= %f %f %f\n",
          eulcam[icam][0],eulcam[icam][1],eulcam[icam][2]);

}

/****************************************************************************/
/* ifrom or ito = 1 => input set
 * ifrom or ito = 2 => A set
 * ifrom or ito = 3 => B set
 */

void copy_fpg_pars(int ifrom, int ito)
{
  double tmp;
  int i, iccd;

  for(i=0;i<3;++i) {
    if(ifrom == 1)
      tmp = eulcam[icam][i];
    else if(ifrom == 2)
      tmp = aeulcam[icam][i];
    else if(ifrom == 3)
      tmp = beulcam[icam][i];
    if(ito == 1)
      eulcam[icam][i] = tmp;
    else if(ito == 2)
      aeulcam[icam][i] = tmp;
    else if(ito == 3)
      beulcam[icam][i] = tmp;
  }

  for(i=0;i<6;++i) {
    if(ifrom == 1)
      tmp = optcon[icam][i];
    else if(ifrom == 2)
      tmp = aoptcon[icam][i];
    else if(ifrom == 3)
      tmp = boptcon[icam][i];
    if(ito == 1)
      optcon[icam][i] = tmp;
    else if(ito == 2)
      aoptcon[icam][i] = tmp;
    else if(ito == 3)
      boptcon[icam][i] = tmp;
  }

 for(iccd=0;iccd<NCCD;++iccd) {
   for(i=0;i<2;++i) {
     if(ifrom == 1)
       tmp = ccdxy0[icam][iccd][i];
     else if(ifrom == 2)
       tmp = accdxy0[icam][iccd][i];
     else if(ifrom == 3)
       tmp = bccdxy0[icam][iccd][i];
     if(ito == 1)
       ccdxy0[icam][iccd][i] = tmp;
     else if(ito == 2)
       accdxy0[icam][iccd][i] = tmp;
     else if(ito == 3)
       bccdxy0[icam][iccd][i] = tmp;
   }
 }

 for(iccd=0;iccd<NCCD;++iccd) {
   for(i=0;i<2;++i) {
     if(ifrom == 1)
       tmp = pixsz[icam][iccd][i];
     else if(ifrom == 2)
       tmp = apixsz[icam][iccd][i];
     else if(ifrom == 3)
       tmp = bpixsz[icam][iccd][i];
     if(ito == 1)
       pixsz[icam][iccd][i] = tmp;
     else if(ito == 2)
       apixsz[icam][iccd][i] = tmp;
     else if(ito == 3)
       bpixsz[icam][iccd][i] = tmp;
   }
 }

 for(iccd=0;iccd<NCCD;++iccd) {
   if(ifrom == 1)
     tmp = ccdang[icam][iccd];
   else if(ifrom == 2)
     tmp = accdang[icam][iccd];
   else if(ifrom == 3)
     tmp = bccdang[icam][iccd];
   if(ito == 1)
     ccdang[icam][iccd] = tmp;
   else if(ito == 2)
     accdang[icam][iccd] = tmp;
   else if(ito == 3)
     bccdang[icam][iccd] = tmp;
 }

 for(iccd=0;iccd<NCCD;++iccd) {
   for(i=0;i<2;++i) {
     if(ifrom == 1)
       tmp = ccdtilt[icam][iccd][i];
     else if(ifrom == 2)
       tmp = accdtilt[icam][iccd][i];
     else if(ifrom == 3)
       tmp = bccdtilt[icam][iccd][i];
     if(ito == 1)
       ccdtilt[icam][iccd][i] = tmp;
     else if(ito == 2)
       accdtilt[icam][iccd][i] = tmp;
     else if(ito == 3)
       bccdtilt[icam][iccd][i] = tmp;
   }
 }

 if(ifrom == 1)
   tmp = asymang[icam];
 else if(ifrom == 2)
   tmp = aasymang[icam];
 else if(ifrom == 3)
   tmp = basymang[icam];
 if(ito == 1)
   asymang[icam] = tmp;
 else if(ito == 2)
   aasymang[icam] = tmp;
 else if(ito == 3)
   basymang[icam] = tmp;

 if(ifrom == 1)
   tmp = asymfac[icam];
 else if(ifrom == 2)
   tmp = aasymfac[icam];
 else if(ifrom == 3)
   tmp = basymfac[icam];
 if(ito == 1)
   asymfac[icam] = tmp;
 else if(ito == 2)
   aasymfac[icam] = tmp;
 else if(ito == 3)
   basymfac[icam] = tmp;
}

/****************************************************************************/
/* xtr - translate in +x direction by xtr mm
 * Take A set as input, B set as output
 */

void xtr_focal_plane(double xtr)
{
  int iccd;

  copy_fpg_pars(2,3);
  for(iccd=0;iccd<NCCD;++iccd) {
    ccdxy0[icam][iccd][0] += xtr;
  }
  copy_fpg_pars(3,2);
}

/****************************************************************************/
/* ytr - translate in +y direction by ytr mm
 * Take A set as input, B set as output
 */

void ytr_focal_plane(double ytr)
{
  int iccd;

  copy_fpg_pars(2,3);
  for(iccd=0;iccd<NCCD;++iccd) {
    bccdxy0[icam][iccd][1] += ytr;
  }
  copy_fpg_pars(3,2);
}

/****************************************************************************/




/****************************************************************************/
/* zang - rotate around z-axis by zang degrees
 * Take A set as input, B set as output
 */

//  double aeulcam[NCAM][3], aoptcon[NCAM][6], accdxy0[NCAM][NCCD][2];
//  double apixsz[NCAM][NCCD][2], accdang[NCAM][NCCD], accdtilt[NCAM][NCCD][2];
//  double aasymang[NCAM], aasymfac[NCAM];

void zrot_focal_plane(double zang)
{
  int iccd;
  double xya[2], xyb[2];

  copy_fpg_pars(2,3);
  for(iccd=0;iccd<NCCD;++iccd) {
    xyrotate(-zang,&accdxy0[icam][iccd][0],&bccdxy0[icam][iccd][0]);
  }
  for(iccd=0;iccd<NCCD;++iccd) {
    bccdang[icam][iccd] = accdang[icam][iccd] + zang;
  }
  copy_fpg_pars(3,2);
}

/****************************************************************************/

int main(int argc, char *argv[])
{
  double fpdx, fpdy, fprot, xang, yang, zang, eul123[3];
  FILE *fpf;

  icam = 0;
  dtor = M_PI/180.0;

  if(argc != 4) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }
  fpdx = atof(argv[1]);   // mm
  fpdy = atof(argv[2]);   // mm
  fprot = atof(argv[3]);  // degrees

  fpf = fopen("fpg_pars.txt","r");
  read_fpg_pars(fpf,stderr,0);
  fclose(fpf);

  yang = (fpdx/optcon[icam][0])/dtor;
  xang = (-fpdy/optcon[icam][0])/dtor;
  zang = -fprot;
  eul123[0] = xang;
  eul123[1] = yang;
  eul123[2] = zang;

  copy_fpg_pars(1,2);
  xtr_focal_plane(fpdx);
  ytr_focal_plane(fpdy);
  zrot_focal_plane(fprot);
  copy_fpg_pars(2,1);
  rotate_cam_frame(eul123);

  print_fpg_pars_deriv(stdout,0);
}

/******************************************************************************/

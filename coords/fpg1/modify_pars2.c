/* modify_pars2.c
 * Alan M. Levine
 * June 11, 2017
 * Heritage: modify_pars1.c, fpg2.c
 *
 * Write a randomly modified set of FPG parameters for use in testing.
 */

// WILL ONLY DO ONE CAMERA AT A TIME

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "seed_tw_ran.h"
#include "twister2.h"

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

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
int icam;

/******************************************************************************/
// All input angles are in degrees.

void read_fpg_pars(FILE *fpin, FILE *fpo)
{
  char pardesc[100];
  int nrd, iccd, j;

  for(j=0;j<3;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&keulcam[icam][j],&eulcam[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,keulcam[icam][j],eulcam[icam][j]);
    eulcam0[icam][j] = eulcam[icam][j];
  }
  for(j=0;j<6;++j) {
    fscanf(fpin,"%s %d %lf",pardesc,&koptcon[icam][j],&optcon[icam][j]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,koptcon[icam][j],optcon[icam][j]);
    optcon0[icam][j] = optcon[icam][j];
  }
  // double asymang[NCAM][NCCD], asymfac[NCAM][NCCD];
  fscanf(fpin,"%s %d %lf",pardesc,&kasymang[icam],&asymang[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymang[icam],asymang[icam]);
  asymang0[icam] = asymang[icam];
  fscanf(fpin,"%s %d %lf",pardesc,&kasymfac[icam],&asymfac[icam]);
  fprintf(fpo,"%s        %d  %f\n",pardesc,kasymfac[icam],asymfac[icam]);
  asymfac0[icam] = asymfac[icam];
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdxy0[icam][iccd][j],&ccdxy0[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
      ccdxy00[icam][iccd][j] = ccdxy0[icam][iccd][j];
    }
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kpixsz[icam][iccd][j],&pixsz[icam][iccd][j]);
      fprintf(fpo,"%s        %d  %f\n",pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
      pixsz0[icam][iccd][j] = pixsz[icam][iccd][j];
    }
    fscanf(fpin,"%s %d %lf",pardesc,&kccdang[icam][iccd],&ccdang[icam][iccd]);
    fprintf(fpo,"%s        %d  %f\n",pardesc,kccdang[icam][iccd],ccdang[icam][iccd]);
    ccdang0[icam][iccd] = ccdang[icam][iccd];
    for(j=0;j<2;++j) {
      fscanf(fpin,"%s %d %lf",pardesc,&kccdtilt[icam][iccd][j],&ccdtilt[icam][iccd][j]);
      fprintf(fpo,"%s       %d   %f\n",pardesc,kccdtilt[icam][iccd][j],ccdtilt[icam][iccd][j]);
      ccdtilt0[icam][iccd][j] = ccdtilt[icam][iccd][j];
    }
  }
}

/******************************************************************************/
// All angles are in degrees.

void print_fpg_pars(FILE *fpo)
{
  char pardesc[100];
  int iccd, j, jcam;


  jcam = icam + 1;
  for(j=0;j<3;++j) {
    sprintf(pardesc,"ang%1d_cam%1d",j+1,jcam);
    fprintf(fpo,"%-18s      %d     %f\n",pardesc,keulcam[icam][j],eulcam[icam][j]);
  }
  sprintf(pardesc,"fl_cam%1d",jcam);
  fprintf(fpo,"%-18s      %d     %f\n",pardesc,koptcon[icam][0],optcon[icam][0]);
  for(j=1;j<6;++j) {
    sprintf(pardesc,"opt_coef%1d_cam%1d",j,jcam);
    fprintf(fpo,"%-18s      %d     %.8f\n",pardesc,koptcon[icam][j],optcon[icam][j]);
  }
  sprintf(pardesc,"asymang_cam%1d",jcam);
  fprintf(fpo,"%-18s      %d     %11.8f\n",
	    pardesc,kasymang[icam],asymang[icam]);
  sprintf(pardesc,"asymfac_cam%1d",jcam);
  fprintf(fpo,"%-18s      %d     %11.8f\n",
	    pardesc,kasymfac[icam],asymfac[icam]);
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      if(j==0)
	sprintf(pardesc,"x0_ccd%1d_cam%1d",iccd+1,jcam);
      else
	sprintf(pardesc,"y0_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s      %d     %f\n",pardesc,kccdxy0[icam][iccd][j],ccdxy0[icam][iccd][j]);
    }
    for(j=0;j<2;++j) {
      if(j==0)
	sprintf(pardesc,"pix_x_ccd%1d_cam%1d",iccd+1,jcam);
      else
	sprintf(pardesc,"pix_y_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s      %d     %f\n",pardesc,kpixsz[icam][iccd][j],pixsz[icam][iccd][j]);
    }
    sprintf(pardesc,"ang_ccd%1d_cam%1d",iccd+1,jcam);
    fprintf(fpo,"%-18s      %d     %f\n",pardesc,kccdang[icam][iccd],ccdang[icam][iccd]);
    for(j=0;j<2;++j) {
      if(j==0)
	sprintf(pardesc,"tilt_x_ccd%1d_cam%1d",iccd+1,jcam);
      else
	sprintf(pardesc,"tilt_y_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s      %d     %f\n",pardesc,kccdtilt[icam][iccd][j],ccdtilt[icam][iccd][j]);
    }
  }
}

/******************************************************************************/

void usage(char *prog)
{
  fprintf(stderr,"%s has incorrect command line parameters. Exiting...\n",prog);
  exit(-1);
}

/****************************************************************************/
/*  Typical input values are:
 *    dcam_mn        0.15 degrees
 *    dfl_mn         0.3 mm
 *    dopt_mn        0.002 radians
 *    dasf_mn        0.00001 
 *    dxy0_mn        0.03 mm
 *    dpsz_mn        0.000015 mm
 *    dang_mn        1.0 degrees
 *    dtilt_mn       0.02 degrees
 */

int main(int argc, char *argv[])
{
  double dcam_mn, dfl_mn, dopt_mn, dxy0_mn, dpsz_mn, dang_mn, dtilt_mn;
  double dasf_mn;
  double th_fov;
  int iccd, j;

  unsigned long seed;

  if(argc != 9) {
    usage(argv[0]);
  }

  seed = get_tw_seed();
  fprintf(stderr,"seed = %lu \n",seed);
  seedMT(seed);

  th_fov = 0.2;   // rough value in radians 

  icam = 0;                 // Camera no. is immaterial in this program
  // The follow parameters are standard deviations. 
  dcam_mn = atof(argv[1]);      // Camera euler angle
  dfl_mn = atof(argv[2]);       // focal length at enter of field
  dopt_mn = atof(argv[3]);      // optical distortion coefficients
  dasf_mn = atof(argv[4]);      // asymmetry factor, the asymmetry angle will be completely random
  dxy0_mn = atof(argv[5]);      // CCD origin positions
  dpsz_mn = atof(argv[6]);      // Pixel sizes
  dang_mn = atof(argv[7]);      // CCD in-plane rotation
  dtilt_mn = atof(argv[8]);     // CCD tilts

  read_fpg_pars(stdin,stderr);

  // int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
  // int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
  // int kasymang[NCAM], kasymfac[NCAM];
  if(keulcam[icam][0] == 1)
    eulcam[icam][0] = eulcam0[icam][0] + gasdev()*dcam_mn;
  if(keulcam[icam][1] == 1)
    eulcam[icam][1] = eulcam0[icam][1] + gasdev()*dcam_mn;
  if(keulcam[icam][2] == 1)
    eulcam[icam][2] = eulcam0[icam][2] + gasdev()*(2.0/th_fov)*dcam_mn;

  if(koptcon[icam][0] == 1)
    optcon[icam][0] = optcon0[icam][0] + gasdev()*dfl_mn;  
  for(j=1;j<NOPTCON;++j) {
    if(koptcon[icam][j] == 1) {
      optcon[icam][j] = optcon0[icam][j] + 
	(gasdev()*(j*j)*dopt_mn*pow(th_fov,-(1.0*j - 1.0)));
    }
  }
  if(kasymang[icam] == 1)
    asymang[icam] = 180.0*(ranMT() - 0.5);
  if(kasymfac[icam] == 1)
    asymfac[icam] = 1.0 + gasdev()*dasf_mn;
  for(iccd=0;iccd<NCCD;++iccd) {
    for(j=0;j<2;++j) {
      if(kccdxy0[icam][iccd][j] == 1)
	ccdxy0[icam][iccd][j] = ccdxy00[icam][iccd][j] + gasdev()*dxy0_mn;
      if(kpixsz[icam][iccd][j] == 1)
	pixsz[icam][iccd][j] = pixsz0[icam][iccd][j] + gasdev()*dpsz_mn;
      if(kccdtilt[icam][iccd][j] == 1)
	ccdtilt[icam][iccd][j] = ccdtilt0[icam][iccd][j] + gasdev()*dtilt_mn;
    }
    if(kccdang[icam][iccd] == 1)
      ccdang[icam][iccd] = ccdang0[icam][iccd] + gasdev()*dang_mn;
  }

  print_fpg_pars(stdout);
  exit(1);
}

/****************************************************************************/

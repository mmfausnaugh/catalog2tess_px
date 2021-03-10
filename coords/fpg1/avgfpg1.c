/* avgfpg1.c
 * Alan M. Levine
 * March 26, 2018
 * Heritage: chkc1.c, fpg4.c
 *
 * Read a number of fpg_pars.txt files and write such a file with the averages
 * of the various parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6
#define MAX_FILES 100   

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
double asymang[NCAM], asymfac[NCAM];
int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];
double deleul[NCAM][3], delopt[NCAM][6], delxy0[NCAM][NCCD][2];
double delsz[NCAM][NCCD][2], delang[NCAM][NCCD], deltilt[NCAM][NCCD][2];
double dasymang[NCAM], dasymfac[NCAM];
int icam;

double eulcam2[4][3], optcon2[4][6], ccdxy02[4][4][2];
double pixsz2[4][4][2], ccdang2[4][4], ccdtilt2[4][4][2];
double asymang2[NCAM], asymfac2[NCAM];

double eulcam3[4][3], optcon3[4][6], ccdxy03[4][4][2];
double pixsz3[4][4][2], ccdang3[4][4], ccdtilt3[4][4][2];
double asymang3[NCAM], asymfac3[NCAM];
double wttot;

double dtor;

char fpgfile[MAX_FILES][128];
int numin;
double wtin[MAX_FILES];

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

int read_fpg_name(FILE *fpin, int jin)
{
  int nrd;

  nrd = fscanf(fpin,"%s %lf",&fpgfile[jin][0],&wtin[jin]);
  return(nrd);
}

/******************************************************************************/

/*
double eulcam3[4][3], optcon3[4][6], ccdxy03[4][4][2];
double pixsz3[4][4][2], ccdang3[4][4], ccdtilt3[4][4][2];
double asymang3[NCAM], asymfac3[NCAM];
*/

void init_fpg_avgs(int icam)
{
  int j, k, iccd;

  for(j=0;j<3;++j) {
    eulcam3[icam][j] = 0.0;
  }
  for(j=0;j<6;++j) {
    optcon3[icam][j] = 0.0;
  }
  for(j=0;j<3;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdxy03[icam][iccd][j] = 0.0;
    }
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      pixsz3[icam][iccd][j] = 0.0;
    }
  }
  for(iccd=0;iccd<4;++iccd) {
    ccdang3[icam][iccd] = 0.0;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdtilt3[icam][iccd][j] = 0.0;
    }
  }
  asymang3[icam] = 0.0;
  asymfac3[icam] = 0.0;
}

/******************************************************************************/

void incr_fpg_avgs(int icam, double wt)
{
  int j, k, iccd;

  wttot += wt;
  for(j=0;j<3;++j) {
    eulcam3[icam][j] += wt*eulcam[icam][j];
  }
  for(j=0;j<6;++j) {
    optcon3[icam][j] += wt*optcon[icam][j];
  }
  for(j=0;j<3;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdxy03[icam][iccd][j] += wt*ccdxy0[icam][iccd][j];
    }
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      pixsz3[icam][iccd][j] += wt*pixsz[icam][iccd][j];
    }
  }
  for(iccd=0;iccd<4;++iccd) {
    ccdang3[icam][iccd] += wt*ccdang[icam][iccd];
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdtilt3[icam][iccd][j] += wt*ccdtilt[icam][iccd][j];
    }
  }
  asymang3[icam] += wt*asymang[icam];
  asymfac3[icam] += wt*asymfac[icam];
}

/******************************************************************************/

void compute_avg_fpg(int icam, double wt)
{
  int j, k, iccd;
  double xwt;

  xwt = 1.0/wt;
  for(j=0;j<3;++j) {
    eulcam[icam][j] = xwt*eulcam3[icam][j];
  }
  for(j=0;j<6;++j) {
    optcon[icam][j] = xwt*optcon3[icam][j];
  }
  for(j=0;j<3;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdxy0[icam][iccd][j] = xwt*ccdxy03[icam][iccd][j];
    }
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      pixsz[icam][iccd][j] = xwt*pixsz3[icam][iccd][j];
    }
  }
  for(iccd=0;iccd<4;++iccd) {
    ccdang[icam][iccd] = xwt*ccdang3[icam][iccd];
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdtilt[icam][iccd][j] = xwt*ccdtilt3[icam][iccd][j];
    }
  }
  asymang3[icam] = xwt*asymang3[icam];
  asymfac3[icam] = xwt*asymfac3[icam];
}

/******************************************************************************/
// maximum number of fpg_pars.txt files that may be read and averaged

int main(int argc, char *argv[])
{
  int icam, jin, nrd;
  FILE *fpf;

  if(argc != 1) {
    fprintf(stderr,"ERROR: argc = %d; exiting...\n",argc);
    exit(-1);
  }

  dtor = M_PI/180.0;
  icam = 0;
  init_fpg_avgs(icam);
  wttot = 0.0;
  numin = 0;

  // stdin - read a list of "fpg_pars.txt"-type files to be averaged with weights
  for(jin=0;jin<MAX_FILES;++jin) {
    nrd = read_fpg_name(stdin,jin);
    if(nrd < 2) break;
    fpf = fopen(&fpgfile[jin][0],"r");
    read_fpg_pars(fpf,stderr);
    fclose(fpf);
    ++numin;
    incr_fpg_avgs(icam,wtin[jin]);
  }
  compute_avg_fpg(icam,wttot);
  print_fpg_pars_deriv(stdout,0);
}

/******************************************************************************/

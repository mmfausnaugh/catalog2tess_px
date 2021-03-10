/* statfpg1.c
 * Alan M. Levine
 * March 27, 2018
 * Heritage: avgfpg1.c
 *
 * Read a number of fpg_pars.txt files and write such a file with the averages
 * and rms variations of the various parameters.
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
double wttot, swtsq;

double eulcam4[4][3], optcon4[4][6], ccdxy04[4][4][2];
double pixsz4[4][4][2], ccdang4[4][4], ccdtilt4[4][4][2];
double asymang4[NCAM], asymfac4[NCAM];

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

void print_fpg_pars_stats(FILE *fpo, int ido)
{
  char pardesc[100];
  int iccd, j, jcam;

  jcam = icam + 1;
  fprintf(stderr,"print_fpg_pars_deriv(): beginning...\n");
  for(j=0;j<3;++j) {
    sprintf(pardesc,"ang%1d_cam%1d",j+1,jcam);
    fprintf(fpo,"%-18s        %d  %f",
              pardesc,keulcam[icam][j],eulcam[icam][j]);
    if(ido == 1)
      fprintf(fpo," %f %f\n",eulcam3[icam][j],eulcam4[icam][j]);
    else
      fprintf(fpo,"\n");
 }
  sprintf(pardesc,"fl_cam%1d",jcam);
  fprintf(fpo,"%-18s        %d  %f",
          pardesc,koptcon[icam][0],optcon[icam][0]);
  if(ido == 1)
    fprintf(fpo," %f %f\n",optcon3[icam][0],optcon4[icam][0]);
  else
    fprintf(fpo,"\n");
  for(j=1;j<6;++j) {
    sprintf(pardesc,"opt_coef%1d_cam%1d",j,jcam);
    fprintf(fpo,"%-18s        %d  %11.8f",
            pardesc,koptcon[icam][j],optcon[icam][j]);
    if(ido == 1)
      fprintf(fpo," %f %f\n",optcon3[icam][j],optcon4[icam][j]);
    else
      fprintf(fpo,"\n");
  }
  sprintf(pardesc,"asymang_cam%1d",jcam);
  fprintf(fpo,"%-18s        %d  %11.8f",
            pardesc,kasymang[icam],asymang[icam]);
  if(ido == 1)
    fprintf(fpo," %f %f\n",asymang3[icam],asymang4[icam]);
  else
    fprintf(fpo,"\n");
  sprintf(pardesc,"asymfac_cam%1d",jcam);
  fprintf(fpo,"%-18s        %d  %11.8f",
            pardesc,kasymfac[icam],asymfac[icam]);
  if(ido == 1)
    fprintf(fpo," %f %f\n",asymfac3[icam],asymfac4[icam]);
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
      if(ido == 1)
        fprintf(fpo," %f %f\n",ccdxy03[icam][iccd][j],ccdxy04[icam][iccd][j]);
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
      if(ido == 1)
        fprintf(fpo," %f %f\n",pixsz3[icam][iccd][j],pixsz4[icam][iccd][j]);
      else
        fprintf(fpo,"\n");
    }
    sprintf(pardesc,"ang_ccd%1d_cam%1d",iccd+1,jcam);
    fprintf(fpo,"%-18s        %d  %f",
            pardesc,kccdang[icam][iccd],ccdang[icam][iccd]);
    if(ido == 1)
      fprintf(fpo," %f %f\n",ccdang3[icam][iccd],ccdang4[icam][iccd]);
    else
      fprintf(fpo,"\n");
    for(j=0;j<2;++j) {
      if(j==0)
        sprintf(pardesc,"tilt_x_ccd%1d_cam%1d",iccd+1,jcam);
      else
        sprintf(pardesc,"tilt_y_ccd%1d_cam%1d",iccd+1,jcam);
      fprintf(fpo,"%-18s        %d  %f",
              pardesc,kccdtilt[icam][iccd][j],ccdtilt[icam][iccd][j]);
      if(ido == 1)
        fprintf(fpo," %f %f\n",ccdtilt3[icam][iccd][j],ccdtilt4[icam][iccd][j]);
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
    eulcam2[icam][j] = 0.0;
    eulcam3[icam][j] = 0.0;
    eulcam4[icam][j] = 0.0;
  }
  for(j=0;j<6;++j) {
    optcon2[icam][j] = 0.0;
    optcon3[icam][j] = 0.0;
    optcon4[icam][j] = 0.0;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdxy02[icam][iccd][j] = 0.0;
      ccdxy03[icam][iccd][j] = 0.0;
      ccdxy04[icam][iccd][j] = 0.0;
    }
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      pixsz2[icam][iccd][j] = 0.0;
      pixsz3[icam][iccd][j] = 0.0;
      pixsz4[icam][iccd][j] = 0.0;
    }
  }
  for(iccd=0;iccd<4;++iccd) {
    ccdang2[icam][iccd] = 0.0;
    ccdang3[icam][iccd] = 0.0;
    ccdang4[icam][iccd] = 0.0;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      ccdtilt2[icam][iccd][j] = 0.0;
      ccdtilt3[icam][iccd][j] = 0.0;
      ccdtilt4[icam][iccd][j] = 0.0;
    }
  }
  asymang2[icam] = 0.0;
  asymang3[icam] = 0.0;
  asymang4[icam] = 0.0;
  asymfac2[icam] = 0.0;
  asymfac3[icam] = 0.0;
  asymfac4[icam] = 0.0;
}

/******************************************************************************/
/* wtavg = sum(wt*val)/sum(wt)
 * wtvar = sum(wt*val*val)/sum(wt) - wtavg*wtavg
*/

void incr_fpg_avgs(int icam, double wt)
{
  int j, k, iccd;
  double v;

  wttot += wt;
  swtsq += wt*wt;
  for(j=0;j<3;++j) {
    v = eulcam[icam][j];
    eulcam3[icam][j] += wt*v;
    eulcam2[icam][j] += wt*v*v;
  }
  for(j=0;j<6;++j) {
    v = optcon[icam][j];
    optcon3[icam][j] += wt*v;
    optcon2[icam][j] += wt*v*v;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      v = ccdxy0[icam][iccd][j];
      ccdxy03[icam][iccd][j] += wt*v;
      ccdxy02[icam][iccd][j] += wt*v*v;
    }
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      v = pixsz[icam][iccd][j];
      pixsz3[icam][iccd][j] += wt*v;
      pixsz2[icam][iccd][j] += wt*v*v;
    }
  }
  for(iccd=0;iccd<4;++iccd) {
    v = ccdang[icam][iccd];
    ccdang3[icam][iccd] += wt*v;
    ccdang2[icam][iccd] += wt*v*v;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      v = ccdtilt[icam][iccd][j];
      ccdtilt3[icam][iccd][j] += wt*v;
      ccdtilt2[icam][iccd][j] += wt*v*v;
    }
  }
  v = asymang[icam];
  asymang3[icam] += wt*v;
  asymang2[icam] += wt*v*v;
  v = asymfac[icam];
  asymfac3[icam] += wt*v;
  asymfac2[icam] += wt*v*v;
}

/******************************************************************************/

void compute_avg_fpg(int icam, double wt)
{
  int j, k, iccd;
  double av, var, xwt;

  xwt = 1.0/wt;
  for(j=0;j<3;++j) {
    av = xwt*eulcam3[icam][j];
    eulcam[icam][j] = av;
    var = xwt*eulcam2[icam][j] - av*av;
    if(var > 0.0) { 
      eulcam3[icam][j] = sqrt(var);
      eulcam4[icam][j] = xwt*sqrt(swtsq*var);
    }
    else
      eulcam3[icam][j] = 0.0;
  }
  for(j=0;j<6;++j) {
    av = xwt*optcon3[icam][j];
    optcon[icam][j] = av;
    var = xwt*optcon2[icam][j] - av*av;
    if(var > 0.0) {
      optcon3[icam][j] = sqrt(var);
      optcon4[icam][j] = xwt*sqrt(swtsq*var);
    }
    else
      optcon3[icam][j] = 0.0;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      av = xwt*ccdxy03[icam][iccd][j];
      ccdxy0[icam][iccd][j] = av;
      var = xwt*ccdxy02[icam][iccd][j] - av*av;
      if(var > 0.0) {
	ccdxy03[icam][iccd][j] = sqrt(var);
	ccdxy04[icam][iccd][j] = xwt*sqrt(swtsq*var);
      }
      else
	ccdxy03[icam][iccd][j] = 0.0;
    }
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      av = xwt*pixsz3[icam][iccd][j];
      pixsz[icam][iccd][j] = av;
      var = xwt*pixsz2[icam][iccd][j] - av*av;
      if(var > 0.0) {
	pixsz3[icam][iccd][j] = sqrt(var);
	pixsz4[icam][iccd][j] = xwt*sqrt(swtsq*var);
      }
      else
	pixsz3[icam][iccd][j] = 0.0;
    }
  }
  for(iccd=0;iccd<4;++iccd) {
    av = xwt*ccdang3[icam][iccd];
    ccdang[icam][iccd] = av;
    var = xwt*ccdang2[icam][iccd] - av*av;
    if(var > 0.0) {
      ccdang3[icam][iccd] = sqrt(var);
      ccdang4[icam][iccd] = xwt*sqrt(swtsq*var);
    }
    else
      ccdang3[icam][iccd] = 0.0;
  }
  for(j=0;j<2;++j) {
    for(iccd=0;iccd<4;++iccd) {
      av = xwt*ccdtilt3[icam][iccd][j];
      ccdtilt[icam][iccd][j] = av;
      var = xwt*ccdtilt2[icam][iccd][j] - av*av;
      if(var > 0.0) {
	ccdtilt3[icam][iccd][j] = sqrt(var);
	ccdtilt4[icam][iccd][j] = xwt*sqrt(swtsq*var);
      }
      else
	ccdtilt3[icam][iccd][j] = 0.0;
    }
  }
  av = xwt*asymang3[icam];
  asymang[icam] = av;
  var = xwt*asymang2[icam] - av*av;
  if(var > 0.0) {
    asymang3[icam] = sqrt(var);
    asymang4[icam] = xwt*sqrt(swtsq*var);
  }
  else
    asymang3[icam] = 0.0;
  av = xwt*asymfac3[icam];
  asymfac[icam] = av;
  var = xwt*asymfac2[icam] - av*av;
  if(var > 0.0) {
    asymfac3[icam] = sqrt(var);
    asymfac4[icam] = xwt*sqrt(swtsq*var);
  }
  else
    asymfac3[icam] = 0.0;
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
  swtsq = 0.0;
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
  print_fpg_pars_stats(stdout,1);
}

/******************************************************************************/

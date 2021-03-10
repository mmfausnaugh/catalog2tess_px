/* orient1.c
 * Alan M. Levine
 * May 9, 2018
 * Heritage: starspx6.c
 *
 * Compute relative orientations of cameras and the spacecraft frame.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6
#define MAX_N_EUL 500

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
double asymang[NCAM], asymfac[NCAM];

int keulcam[NCAM][3], koptcon[NCAM][6], kccdxy0[NCAM][NCCD][2];
int kpixsz[NCAM][NCCD][2], kccdang[NCAM][NCCD], kccdtilt[NCAM][NCCD][2];
int kasymang[NCAM], kasymfac[NCAM];

double ecam[MAX_N_EUL][4][3];
int nset, icam;
double rmat1[3][3];
double rma1[3][3], rma2[3][3], rma3[3][3], rma4[3][3];
double rmb1[3][3], rmb2[3][3], rmb3[3][3], rmb4[3][3];
double rmc1[3][3], rmc2[3][3], rmc3[3][3], rmc4[3][3];
double rmd1[3][3], rmd2[3][3], rmd3[3][3], rmd4[3][3];
double rm0a1[3][3], rm0a2[3][3], rm0a3[3][3], rm0a4[3][3];
double rm0b1[3][3], rm0b2[3][3], rm0b3[3][3], rm0b4[3][3];
double rm0c1[3][3], rm0c2[3][3], rm0c3[3][3], rm0c4[3][3];
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

void sc_to_cam_mat_arg(double euler[3], double rmt[3][3])
{
  double xeul[3], angle;
  int i;

  for(i=0;i<3;++i) {
    xeul[i] = dtor*euler[i];
  }

  eulerm323(xeul,rmt);
  // prmat(rmt,stderr);
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

void read_euler_set(int imno, int jcam, char *file, FILE *fdiag)
{
  int j;
  FILE *fin;


  fin = fopen(file,"r");
  read_fpg_pars(fin,fdiag,0);
  fclose(fin);
  fprintf(fdiag," iset,jcam,file = %d %d %s\n",imno,jcam,file);

  for(j=0;j<3;++j) {
    ecam[imno][jcam-1][j] = eulcam[0][j];
    fprintf(fdiag,"imno,jcam,j,ecam = %d %d %d %f\n",
	    imno,jcam,j,ecam[imno][jcam-1][j]);
  }

}

/******************************************************************************/

int check_cam_set(int kc[4])
{
  int i;

  for(i=0;i<4;++i) {
    if(kc[i] == 0) {
      return(0);
    }
  }

  return(1);
}

/******************************************************************************/

void read_fpg_pars_euler(FILE *fpf, FILE *fdiag)
{
  char filein[256];
  int i, k, imno, jcam, kcam[4];

  imno = 0;
  for(k=0;k<4;++k) kcam[k] = 0;

  while(fscanf(fpf,"%d %s",&jcam,filein) != EOF) {
    if( (jcam >= 1) && (jcam <= 4) ) {
      kcam[jcam-1] = 1;
    }
    else {
      fprintf(fdiag,"read_fpg_pars_euler(): Error - jcam = %d, exiting.\n",jcam);
      exit(-1);
    }
    read_euler_set(imno,jcam,filein,fdiag);

    for(i=0;i<3;++i) {
      fscanf(fpf,"%d %s",&jcam,filein);
      if( (jcam >= 1) && (jcam <= 4) ) {
	kcam[jcam-1] = 1;
      }
      else {
	fprintf(fdiag,"read_fpg_pars_euler(): Error - jcam = %d, exiting.\n",jcam);
	exit(-1);
      }
      read_euler_set(imno,jcam,filein,fdiag);
    }

    if(check_cam_set(kcam) == 0) {
      // error return
      fprintf(fdiag,"check_cam_set() returned 0 - exiting\n");
      exit(-1);
    }
    else {
      ++imno;
    }
  }
  nset = imno;
}

/******************************************************************************/

void copy_rmat(int kcam, double rmat[4][3][3], double rm[3][3])
{
  int j, k;

  for(j=0;j<3;++j) {
    for(k=0;k<3;++k) {
      rm[j][k] = rmat[kcam][j][k];
    }
  }

}

/******************************************************************************/

void rmat_copy(int kcam, double rm[3][3], double rmat[4][3][3])
{
  int j, k;

  for(j=0;j<3;++j) {
    for(k=0;k<3;++k) {
      rmat[kcam][j][k] = rm[j][k];
    }
  }

}

/******************************************************************************/

void rm_copy(double rma[3][3], double rmb[3][3])
{
  int j, k;

  for(j=0;j<3;++j) {
    for(k=0;k<3;++k) {
      rmb[j][k] = rma[j][k];
    }
  }

}

/******************************************************************************/

void usage(char *prog, int argc)
{
  fprintf(stderr,"Usage: %s camera_num(1-4) RA_sc Dec_sc Roll_sc\n",prog);
  fprintf(stderr,"  < fpg_pars_list > output_file\n");
  exit(-1);
}

/******************************************************************************/
/*

1) Read S/C RA, Dec, roll from the command line
2) Read a list with
cam_no.,  fpg_pars.txt file (path) names
There should be files for each of 4 cameras for each image number.
The first four should have the nominal Euler angles for each camera.
3) Read each file saving 3 Euler angles
4) Keep track of the number of files read

*/
/******************************************************************************/

int main(int argc, char *argv[])
{
  double radecroll[3], eulc[3];
  double rm[3][3], rmt[3][3];
  int iset, j, jcam;
  FILE *fpf;

  dtor = M_PI/180.0;

  if(argc != 4) {
    fprintf(stderr,"ERROR: argc = %d;\n",argc);
    usage(argv[0],argc);
  }
  ra_sc = atof(argv[1]);
  radecroll[0] = ra_sc;
  dec_sc = atof(argv[2]);
  radecroll[1] = dec_sc;
  roll_sc = atof(argv[3]);
  radecroll[2] = roll_sc;
  fprintf(stderr,"ra,dec,roll S/C = %f %f %f\n",ra_sc,dec_sc,roll_sc);

  fpf = fopen("fpg_list.txt","r");
  read_fpg_pars_euler(fpf,stderr);
  fclose(fpf);
  fprintf(stdout,"nset = %d\n",nset);
  if(nset < 3) {
    fprintf(stderr,"nset = %d not enough to work with; exiting\n",nset);
    exit(-1);
  }

  sky_to_sc_mat(radecroll);
  fprintf(stderr,"sky to S/C rotation matrix:   \n");
  prmat(rmat1,stderr);
  fprintf(stderr,"\n");

  // Compute the nominal camera orientation matrices.
  iset = 0;
  for(jcam=0;jcam<4;++jcam) {
    for(j=0;j<3;++j) {
      eulc[j] = ecam[iset][jcam][j];
    }
    sc_to_cam_mat_arg(eulc,rm);
    trans(rm,rmt);
    if(jcam==0) {
      rm_copy(rm,rm0a1);
      rm_copy(rmt,rm0b1);
    }
    else if(jcam==1) {
      rm_copy(rm,rm0a2);
      rm_copy(rmt,rm0b2);
    }
    else if(jcam==2) {
      rm_copy(rm,rm0a3);
      rm_copy(rmt,rm0b3);
    }
    else if(jcam==3) {
      rm_copy(rm,rm0a4);
      rm_copy(rmt,rm0b4);
    }
  }

  /* Compute inverse matrix for nominal camera 1-3 orientation relative
   * to nominal camera 4 orientation.
   */
  matmat(rm0a4,rm0b1,rm0c1);
  matmat(rm0a4,rm0b2,rm0c2);
  matmat(rm0a4,rm0b3,rm0c3);

  for(iset=1;iset<nset;++iset) {

    for(jcam=0;jcam<4;++jcam) {
      for(j=0;j<3;++j) {
	eulc[j] = ecam[iset][jcam][j];
      }
      sc_to_cam_mat_arg(eulc,rm);
      trans(rm,rmt);
      if(jcam==0) {
	rm_copy(rm,rma1);
	rm_copy(rmt,rmb1);
      }
      else if(jcam==1) {
	rm_copy(rm,rma2);
	rm_copy(rmt,rmb2);
      }
      else if(jcam==2) {
	rm_copy(rm,rma3);
	rm_copy(rmt,rmb3);
      }
      else if(jcam==3) {
	rm_copy(rm,rma4);
	rm_copy(rmt,rmb4);
      }
    }

    // Use camera 4 as a reference for the other 3 cameras.
    /* Compute matrix for inferred camera 1-3 orientation relative
     * to inferred camera 4 orientation.
     * Compute product of matrix products.
     */
    matmat(rma1,rmb4,rmc1);
    matmat(rm0c1,rmc1,rmd1);
    printf("Cam 1 rel to 4 (iset = %d): %f %f %f\n",iset,rmd1[0][1],rmd1[0][2],rmd1[1][2]);
    fprintf(stderr,"Cam 1 rel to Cam 4 rotation difference (iset = %d):\n",iset);
    prmat(rmd1,stderr);
    fprintf(stderr,"\n");

    matmat(rma2,rmb4,rmc2);
    matmat(rm0c2,rmc2,rmd2);
    printf("Cam 2 rel to 4 (iset = %d): %f %f %f\n",iset,rmd2[0][1],rmd2[0][2],rmd2[1][2]);
    fprintf(stderr,"Cam 2 rel to Cam 4 rotation difference (iset = %d):\n",iset);
    prmat(rmd2,stderr);
    fprintf(stderr,"\n");

    matmat(rma3,rmb4,rmc3);
    matmat(rm0c3,rmc3,rmd3);
    printf("Cam 3 rel to 4 (iset = %d): %f %f %f\n",iset,rmd3[0][1],rmd3[0][2],rmd3[1][2]);
    fprintf(stderr,"Cam 3 rel to Cam 4 rotation difference (iset = %d):\n",iset);
    prmat(rmd3,stderr);
    fprintf(stderr,"\n");
  }
}
/******************************************************************************/

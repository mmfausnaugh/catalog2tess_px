/* combeul1.c
 * Alan M. Levine
 * March 16, 2018
 *
 * Combine (RA,Dec,roll) with a set of Euler angles or
 * two sets of Euler angles to yield a single set of angles.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"

#define DG_RAD (M_PI/180.0)

double dtor;
double rmat1[3][3], rmatc1[3][3], rmat2[3][3], rmat3[3][3];

/******************************************************************************/

void deg_to_rad(int n, double ain[], double aout[])
{
  int i;

  for(i=0;i<n;++i) {
    aout[i] = DG_RAD*ain[i];
  }
}
/******************************************************************************/

void rad_to_deg(int n, double ain[], double aout[])
{
  int i;

  for(i=0;i<n;++i) {
    aout[i] = ain[i]/DG_RAD;
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
// Multiply a 3-vector in celestial equatorial coordinates by rmat1 to get a 3-vector
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

void sc_to_cam_mat(int icam, double euler[3])
{
  double xeul[3], angle, rm[3][3], rm2[3][3];
  int i;

  for(i=0;i<3;++i) {
    xeul[i] = dtor*euler[i];
  }

  eulerm323(xeul,rm2);
  if (icam < 2) {
    angle = M_PI/2.0;
  }
  else {
    angle = -M_PI/2.0;
  }
  rotm1(2,angle,rm);
  matmat(rm,rm2,rmatc1);
  // prmat(rm,stderr);
  // prmat(rm2,stderr);
  // prmat(rmatc1,stderr);
}

/******************************************************************************/

void get_sc_ra_dec_roll(double rm[3][3], double angs[3])
{
  double eul323[3];

  mateuler323(rm,eul323);  // angles in eul323 will be in radians
  angs[0] = eul323[0]/dtor;
  if (angs[0] < 0.0)
    angs[0] += 360.0;
  angs[1] = ((M_PI/2.0) - eul323[1])/dtor;
  angs[2] = (eul323[2] - M_PI)/dtor;
  if (angs[2] < -180.0)
    angs[2] += 360.0;
}

/******************************************************************************/

void usage(char *prog)
{
  fprintf(stderr,"Usage: %s iin1(1-2), iin2(2), iout(1-2) ang1_0 ang1_1, ang1_2 ang2_0 ang2_1 ang2_2\n",prog);
  fprintf(stderr," > output_file \n");
  exit(-1);
}

/******************************************************************************/
/* input format codes:
 * 1 =>  Ra, Dec, Roll as used for the spacecraft
 * 2 =>  Euler 3-2-3 
 * 3 =>  none
 */

int main(int argc, char *argv[])
{
  int iin1, iin2, iout, i;
  double ang1[3], ang2[3], ang3[3];
  double ang1r[3], ang2r[3], ang3r[3];

  dtor = M_PI/180.0;

  if(argc != 10) {
    usage(argv[0]);
  }

  iin1 = atoi(argv[1]);
  iin2 = atoi(argv[2]);
  iout = atoi(argv[3]);
  ang1[0] = atof(argv[4]);
  ang1[1] = atof(argv[5]);
  ang1[2] = atof(argv[6]);
  ang2[0] = atof(argv[7]);
  ang2[1] = atof(argv[8]);
  ang2[2] = atof(argv[9]);

  if(iin1 == 1) {
    sky_to_sc_mat(ang1);   // => rmat1
  }
  else if(iin1 == 2) {
    deg_to_rad(3,ang1,ang1r);
    eulerm323(ang1r,rmat1);
  }
  else {
    // error - print message and exit
  }

  if(iin2 == 2) {
    deg_to_rad(3,ang2,ang2r);
    eulerm323(ang2r,rmat2);
  }
  else {
    // error - print message and exit
  }

  matmat(rmat2,rmat1,rmat3);
  prmat(rmat1,stdout);
  prmat(rmat2,stdout);
  prmat(rmat3,stdout);

  if(iout == 1) {
    get_sc_ra_dec_roll(rmat3,ang3);
  }
  else if(iout == 2) {
    mateuler323(rmat3,ang3r);
    rad_to_deg(3,ang3r,ang3);
  }
  else {
    // error - print message and exit
  }

  printf("ang1 = %f %f %f\n",ang1[0],ang1[1],ang1[2]);
  printf("ang2 = %f %f %f\n",ang2[0],ang2[1],ang2[2]);
  printf("ang3 = %f %f %f\n",ang3[0],ang3[1],ang3[2]);
}

/******************************************************************************/

/* runsunpos1.c
 * Alan M. Levine
 * February 8, 2018
 *
 * Get the position and velocity of the Earth in its orbit.
 *
 */ 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "astron_rb.h"

/******************************************************************************/

int main(int argc, char *argv[])
{
  double mjd, possun[2], rsun[3], vearth[3];

  if(argc != 2) {
    fprintf(stderr,"Error - program call is incorrect.\n");
    fprintf(stderr,"  Try: %s <MJD>\n",argv[0]);
    exit(-1);
  }
  mjd = atof(argv[1]);

  sunpos(mjd,possun,rsun,vearth);
  printf("mjd,possun= %f    %f %f\n",mjd,possun[0],possun[1]);
  printf(" rsun= %f %f %f\n",rsun[0],rsun[1],rsun[2]);
  printf(" vearth= %f %f %f\n",vearth[0],vearth[1],vearth[2]);
  exit(1);
}

/******************************************************************************/

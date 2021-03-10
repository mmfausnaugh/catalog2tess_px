/* rhist1.c
 * Alan M. Levine
 * July 12, 2018
 * Heritage: radldst1.c
 *
 * Compute histogram of radial coordinates of stars in the focal plane.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "vec.h"
//#include "mat_ra3.h"

#define NSTAR 10000
#define MAXN 200
#define MAXRAD 50.0

double col[NSTAR], row[NSTAR], xfp[NSTAR], yfp[NSTAR], r[NSTAR];
int hist[MAXN];

/******************************************************************************/

void usage(char *prog)
{
  fprintf(stderr,"usage: %s < input file > output_file\n",
          prog);
  exit(1);
}

/******************************************************************************/

int main(int argc, char *argv[])
{
  int i, j, nmax, ib;
  double col, row, x, y, r;

  if(argc != 2) {
    usage(argv[0]);
    exit(-1);
  }
  nmax = atoi(argv[1]);  // max number of bins (less than or equal to MAXN)
  fprintf(stderr,"nmax = %d\n",nmax);

  for(i=0;i<nmax;++i) {
    hist[i] = 0;
  }

  while(scanf("%d %lf %lf %lf %lf %lf",&j,&col,&row,&x,&y,&r)==6) {
    ib = nmax*r/MAXRAD;
    fprintf(stderr,"%f %d\n",r,ib);
    if((ib >= 0) && (ib <= nmax)) {
      ++hist[ib];
    }
  }

  for(i=0;i<nmax;++i) {
    printf("%d %d\n",i,hist[i]);
  }
}

/******************************************************************************/

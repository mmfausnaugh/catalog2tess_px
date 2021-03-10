#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "seed_tw_ran.h"
#include "twister2.h"

#define NMX 100
#define NTEST1 50000

int h1[NMX];

void test1_ranMT()
{
  double x, xmean, xsig, del, chsq, rchsq;
  int i, j, nh, ix;

  // Test uniform distribution

  nh = NMX;
  for(j=0;j<nh;++j) {
    h1[j] = 0;
  }

  for(i=0;i<NTEST1;++i) {
    x = ranMT();
    ix = NMX*x;
    if (  (ix >= 0) && (ix < NMX) )
      ++h1[ix];
  }

  xmean = ((double) NTEST1)/NMX;
  xsig = sqrt(xmean);
  chsq = 0.0;
  fprintf(stderr,"xmean, xsig = %f %f\n",xmean,xsig);
  for(j=0;j<nh;++j) {
    del = h1[j] - xmean;
    chsq += del*del/(xsig*xsig);
    fprintf(stderr,"j,h1[j],del,nsig = %d %d %f %f\n",j,h1[j],del,del/xsig);
    if ( (del < (-10.0*xsig)) || (del > (10.0*xsig)) ) {
      fprintf(stderr,"***j,h1[j],del,nsig = %d %d %f %f\n",j,h1[j],del,del/xsig);
    }
  }
  rchsq = chsq/nh;
  fprintf(stderr,"chsq, rchsq = %f %f\n",chsq,rchsq);
}

int main(int argc, char *argv[])
{
  unsigned long seed;

  seed = get_tw_seed();
  fprintf(stderr,"seed = %lu \n",seed);
  seedMT(seed);

  test1_ranMT();
}

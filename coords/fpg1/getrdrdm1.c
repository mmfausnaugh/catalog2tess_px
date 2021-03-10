/* getrdrdm1.c
 * Alan M. Levine
 * May 16, 2018
 *
 * Read an 11-column "runfpg_fmt" file and write 5 columns
 * (Ra, Dec, Aber Ra, Aber Dec, mag).
 */

#include <stdio.h>
#include <stdlib.h>

main()
{
  double v[5];

  while(scanf("%*d %lf %lf %lf %lf %lf %*d %*lf %*lf %*lf %*lf",
	      &v[0],&v[1],&v[2],&v[3],&v[4]) == 5) {
    printf("%13.6f %13.6f %13.6f %13.6f %.3f\n",
	   v[0],v[1],v[2],v[3],v[4]);
  }

}

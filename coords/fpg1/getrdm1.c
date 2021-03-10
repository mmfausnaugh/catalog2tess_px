/* getrdm1.c
 * Alan M. Levine
 * May 16, 2018
 *
 * Read a 3-column TICid, RA, Dec file and write 3 columns
 * RA, Dec, mag, where mag is just a place holder.
 */

#include <stdio.h>
#include <stdlib.h>

main()
{
  double v[3], tmag;

  tmag = 8.0;   // place holder magnitude

  while(scanf("%*d %lf %lf",
	      &v[0],&v[1]) == 2) {
    printf("%13.6f %13.6f %.3f\n",
	   v[0],v[1],tmag);
  }

}

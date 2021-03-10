/* get3colrdm1.c
 * Alan M. Levine
 * August 30, 2018
 *
 * Read an 11-column guide star file and write 3 columns
 * (Ra, Dec, mag).
 */

#include <stdio.h>
#include <stdlib.h>

main()
{
  char arr[256];
  double v[5];
  int i;

  // Discard first 5 lines
  for(i=0;i<5;++i)
    fgets(arr,250,stdin);

  while(scanf("%*lf %*lf %*lf %*lf %*lf %*lf %*d %lf %lf %lf %*d",
	      &v[0],&v[1],&v[2]) == 3) {
    printf("%13.6f %13.6f %.3f\n",
	   v[0],v[1],v[2]);
  }

}

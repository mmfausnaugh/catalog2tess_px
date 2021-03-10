/* getpix1.c
 * Alan M. Levine
 * July 8, 2017
 *
 */

#include <stdio.h>
#include <stdlib.h>

main()
{
  double col, row;

  while(scanf("%*d %*lf %*lf %*lf %*d %lf %lf",&col,&row) == 2) {
    printf("%13.6f %13.6f\n",col,row);
  }

}

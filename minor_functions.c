#include "minor_functions.h"

double gcd(double x, double y)
{
	printf("Calculating the gcd between %lf and %lf.",x,y);
        int i, x_, y_;
        x_ = floor(fabs(x));
        y_ = floor(fabs(y));
        for(i=x_;i>=1;i--)
        {
        	if(x_%i == 0 && y_%i == 0) break;
        }
      printf(" The result is %lf.\n",(double) i);
        return (double) i;
}

/*********************************************************************************************************************** 

	Used for qsort.
	
***********************************************************************************************************************/

int comparison_function2(const void *a, const void *b)
{
	return ( *(int*)b - *(int*)a );
}

int comparison_function4(const void *a, const void *b)
{
	if( ((const struct store_it *)a)->ratio < ((const struct store_it *)b)->ratio ) return 1;
	return 0;
}


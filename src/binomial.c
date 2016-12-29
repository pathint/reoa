/**************************************************************************
 *  stat.c: Basic Statistics Functions
 *
 *  Xianlong Wang, Ph.D. 
 *  University of Electronic Science and Technology of China. 
 *  Email: Wang.Xianlong@139.com
 *
 *  Initialization. Nov. 11, 2016.
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "stat.h"


int main(int argc, char* argv[]) {
	if(argc<3){
		printf("ERROR: Two few agruments!\n");
		printf("Usage: %s n m\n", argv[0]);
		printf("where m and n are non-negative integer numbers,\n");
		printf("n is the total number of trials,\n");
		printf("m is the number of one of two possible outcomes which have the same proportions (p = 50%%) in the underlying population.\n");
		return -1;
	} 
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	if(n<2||m<0||m>n/2||n>255){
		printf("ERROR: Out-of-range Numbers!!\n");
	}
	else{
		printf("The one-tailed binomial test P value is %.15g\n", bino_p(n, m));
	}
	return 0;
}

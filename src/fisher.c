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
	if(argc<5){
		printf("ERROR: Two few agruments!\n");
		printf("Usage: %s a b c d\n", argv[0]);
		printf("where a, b, c and d are four non-negative integer numbers,\n");
		printf("a and b are observed frequencies for two outcomes in one group, e.g. the control group,\n");
		printf("c and d are observed frequencies for two outcomes in the other group, e.g., the treatment group.\n");
		return -1;
	} 
	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	int c = atoi(argv[3]);
	int d = atoi(argv[4]);
	if(a<0||b<0||c<0||d<0){
		printf("ERROR: Non-negative integeters only!\n");
	}
	else{
		printf("The two-tailed Fisher exact test P value is %.15g\n", fisher_test(a, b, c, d));
		//printf("Right tail value is %d\n", right_tail(a, b, c, d));
	}
	return 0;
}

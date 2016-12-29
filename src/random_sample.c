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
		printf("Usage: %s population_size sample_size\n", argv[0]);
		printf("       sample_size should be less than population_size and both should be non-negative integers.\n");
		exit(-1);
	} 
	unsigned int n = atoi(argv[1]);
	unsigned int m = atoi(argv[2]);
	if(m>n||m==0||n==0){
		printf("ERROR: population size is less than sample size or one of them is 0!!\n");
		exit(-1);
	}
	unsigned int *sample;
	sample = malloc(sizeof(unsigned int)*m);
	int status;
	status = random_sample(n, m, sample); 
	if(status==0){
	  unsigned int i;
	  for(i=0;i<m;i++) printf("%d\n",sample[i]);
	}
	free(sample);
	return 0;
}

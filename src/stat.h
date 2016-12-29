#include <stdbool.h>
#include "pcg_basic.h"
bool entropy_get(void* dest, size_t n);
bool member_of(unsigned int value, unsigned int n, unsigned int data[n]);
bool random_init();
unsigned char random_bin();
uint32_t random_uint();
uint32_t random_bounded(uint32_t bound);
float random_float();
int random_sample(unsigned int population_size, unsigned int sample_size, 
                   unsigned int sample[sample_size]);

int comp_int(const void *a, const void *b);
int comp_uint(const void *a, const void *b);

double bino_p(int n, int m);
int bino_ctrl(int n, double p);

int right_tail(int a, int b, int c, int d);
double hypergeo_test(int ng, int nl, int cg, int cl);
double   fisher_test(int ng, int nl, int cg, int cl);

struct FDR{
  int index;
  double value;
  double p;
};

double  bh_threshold(int n, double pvalues[n], double alpha);

struct FDR bh_ctrl(int sample_size, 
	   int upper_limit, int counts[upper_limit], int total, 
	   double alpha);

struct FDR bh_eval(int sample_size, 
	   int upper_limit, int counts[upper_limit], int total);

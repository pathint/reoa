/**************************************************************************
 *  stat.c: Basic Statistics Functions
 *
 *  Xianlong Wang, Ph.D. 
 *  University of Electronic Science and Technology of China. 
 *  Email: Wang.Xianlong@139.com
 *
 *  Initialization. Oct. 25, 2016.
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
//#include <math.h>
#include "bino.h" // bino_table, Binomial Cumulated Probability Table (Fair Coin)
#include "stat.h"
//PCG Random Number Generator. See http://www.pcg-random.org/
//#include <fcntl.h>
//#include <stdbool.h>
//#include "pcg_basic.h"

//WARNING: OS dependent?
bool entropy_get(void* dest, size_t n)
{
    FILE* fd = fopen("/dev/urandom", "r");
    if (fd == NULL)
        return false;
    size_t sz=fread(dest, sizeof(uint64_t), n, fd);
    if(sz<n)
	return false;
    return fclose(fd) == 0;
}

bool member_of(unsigned int value, unsigned int n, unsigned int data[n])
{
    unsigned int i;
    bool result=false;
    for(i=0;i<n;i++){
	if(data[i]==value){
          result=true;
	  break;
	}
    }
   return result; 
}

bool random_init() {
    uint64_t seeds[2];
    bool status=entropy_get((void*)seeds, 2);
    pcg32_srandom(seeds[0], seeds[1]);
    return status;
}
//return either 0, 1 randomly
unsigned char random_bin(){
   return pcg32_boundedrand(2);
   //rand()%2;
}

uint32_t random_uint(){
   return pcg32_random();
}

uint32_t random_bounded(uint32_t bound){
   return pcg32_boundedrand(bound);
}

float random_float() {
   return (float) pcg32_random()*(1.0/4294967296.0);   
}

// //Random Generator using OpenSSL implementation
// /* http://stackoverflow.com/questions/822323/how-to-generate-a-random-number-in-c */
// /* random integer in [0, limit) */
// #include <openssl/rand.h>
// unsigned int random_uint(unsigned int limit) {
//     union {
//        unsigned int i;
//        unsigned char c[sizeof(unsigned int)];
//     } u;
// 
//     do {
//       if (!RAND_bytes(u.c, sizeof(u.c))) {
// 	fprintf(stderr, "Can't get random bytes!\n");
// 	exit(1);
// 	}
//     } while (u.i < (-limit % limit)); /* u.i < (2**size % limit) */
//     return u.i % limit;
// }
// 
// /* random double in [0.0, 1.0) */
// double random_double() {
//    union {
//       uint64_t i;
//       unsigned char c[sizeof(uint64_t)];
//    } u;
//     if (!RAND_bytes(u.c, sizeof(u.c))) {
// 	fprintf(stderr, "Error: can't get random bytes!\n");
// 	exit(1);
//    }
//     /* 53 bits / 2**53 */
//     return (u.i >> 11) * (1.0/9007199254740992.0);
// }

// direct translation of random.py in python
// return sample_size unique random elements from a population, 
// sampling without replacement
int random_sample(unsigned int population_size, unsigned int sample_size, 
		  unsigned int sample[sample_size])
{
   if(sample_size>population_size||sample_size<=0){
      fprintf(stderr, "ERROR: sample_size is out of range.\n");
      exit(1);
   }
   random_init();
   unsigned int set_size = 21;
   if(sample_size>5) set_size += (unsigned int) pow(4.0, ceil(log(3*sample_size)/log(4.0)));
   unsigned int i, j;
   if(population_size<=set_size)
   {//pool tracking method
     // pool = list(population)
     unsigned int *pool;
     pool = malloc(sizeof(unsigned int)*population_size);
     for(i=0;i<population_size;i++) pool[i] = i;
     for(i=0;i<sample_size;i++)
     {
	j = (unsigned int) (random_float()*(population_size-i));
	sample[i] = pool[j];
	pool[j] = pool[population_size-i-1];//mov non-selected item into vacancy
     }
     free(pool);
   }
   else
   {  
      unsigned int *selected;
      int selected_pointer=-1; //use array to mimick list
      bool exist;
      selected = malloc(sizeof(int)*set_size);
      for(i=0;i<set_size;i++) selected[i] = -1;//an impossible value
      for(i=0;i<sample_size;i++)
      {
	do{
	   j = (unsigned int) (random_float()*population_size);
	   // check if j in selected
	   exist = false; //not exist 
	   int k = 0;
	   while(k<=selected_pointer)
           {
	      if(selected[k]==j){
	        exist = true;//exist already
	        break; //exit earlier
	      }//if
	      k += 1;
	   }
	}while(exist);
	//selected_pointer += 1;
	//selected = realloc(selected, sizeof(unsigned int)*(selected_pointer+1));
	//selected[selected_pointer]=j;
	selected[++selected_pointer]=j;
	sample[i] = j;
      }//end for
      free(selected);
   }//end else
   return 0;
}

/*******************************
 * bino_table look-up wrapper
 *   n, sample size, greater than or equal to 2
 *   m, number of one outcom, greater than  or equal to 0, 
 *      less than or equal to Ceiling[(n-1)/2]
 ********************************/
//WARNING: upper limit for n is 255, which depends on the table size
//           not checked here. This may lead to Segment Fault error.
double bino_p(int n, int m) {
  if (n < 2 || m <0 || m > (n /2)) return -1;
  int i;
  if ((n-2) % 2) { //odd
    i = (n-3)/2;
    return bino_table[(i+2)*(i+1) + i + m];
  }
  else { // even
    i = (n-2)/2;
    return bino_table[i*(i+3) + m];
  }
}

//given the sample size and the p value
//return the minimal value for one output to produce p value larger than the given p
int bino_ctrl(int n, double p)
{
   int i;
   for (i=0; i<(n+1)/2; i++)
      if (bino_p(n, i)>p) return i;
   return i;
}//end of bio_ctrl

//comparison function for sorting
int comp_double(const void *a, const void *b) {
   double *x = (double *) a;
   double *y = (double *) b;
   if (*x < *y) return -1;
   else if (*x > *y) return 1;
   return 0;
}

int comp_int(const void *a, const void *b) {
   return (*(int *)a - *(int *)b);
}
int comp_uint(const void *a, const void *b) {
   //return (*(int *)v1 - *(int *)v2);
   return (*(unsigned int *)a - *(unsigned int *)b);
}

// pvalues
double bh_threshold(int n, double pvalues[n], double alpha){
   qsort(pvalues, n, sizeof(*pvalues), comp_double);  
   int i;
   for(i=0;i<n;i++)
     if(pvalues[i]>((double)i*alpha/(double)n)) 
		     return pvalues[i>0?i-1:0];
   return pvalues[n-1];
}

struct FDR bh_ctrl(int sample_size, 
	   int upper_limit, int counts[upper_limit], int total,
	   double alpha) {
   int i, acc = 0;
   double p;
   struct FDR result;
   for(i=0;i<upper_limit;i++) {
     acc += counts[i];
     p = bino_p(sample_size, i); 
     if(p> ((double)acc*alpha/(double)total)){
       if(i==0) {
         result.index = i;
	 result.p     = p;
         result.value = p*(double)total/(double)acc;
       }
       else {
         result.index = i-1;
	 p = bino_p(sample_size, i-1);
	 result.p = p;
         result.value = p*(double)total/(double)(acc-counts[i]);
       }
       return result; 
     }
   }
   result.index = upper_limit -1 ;
   p = bino_p(sample_size, upper_limit-1);
   result.p = p;
   result.value = p*(double)total/(double)acc;
   return result;
}


struct FDR bh_eval(int sample_size, 
	   int upper_limit, int counts[upper_limit], int total) {
   int i, acc = 0, index = upper_limit;
   struct FDR result;
   for(i=0;i<=upper_limit;i++)
     if (counts[i]>0) {
	acc += counts[i];
	index = i;
     }
   result.index = index;
   result.value = bino_p(sample_size, index)*(double)total/(double)acc;
   return result;
}


//HyperQuick algorithm
// Aleš Berkopec, HyperQuick algorithm for discrete hypergeometric distribution
// Journal of Discrete Algorithms
// Volume 5, Issue 2, June 2007, Pages 341–347
//         http://dx.doi.org/10.1016/j.jda.2006.01.001
#define ACCURACY DBL_EPSILON
// Eq. (6) 
long double inv_jm(int n, int x, int N, int m) {
   //return (1.0- (long double)x/((long double)m+1.0))/(1.0-((long double)n-1.0-(long double)x)/((long double)N-1.0-(long double)m));
   return (long double)(1.0- x/(m+1.0))/(1.0-(n-1.0-x)/(N-1.0-m));
}

long double hypergeo_p(int n, int x, int N, int M, double eps) {
   int k;
   long double s = 1.0, ak, bk, ck, epsk, jjm, result;
   //printf("Input Values:%d, %d, %d, %d\n", n, x, N, M);
   if ( (n==N && x==M)|| M==0||N==0 ) return 1.0;
   for (k=x;k<(M-1);k++)
     s = s*inv_jm(n, x, N, k) + 1.0;
   ak = s;
   bk = s;
   k = M -1;
   epsk = 2.0*eps;
   while ((k<(N-(n-x))) && (epsk>eps)) {
     ck = ak/bk;
     jjm = inv_jm(n, x, N, k);
     ak = ak*jjm;
     bk = bk*jjm + 1.0;
     epsk = (N-(n-x)-1-k)*(ck-ak/bk);
     k += 1;
   }
   result = 1.0 - (ak/bk - epsk/2.0);
   return result;
}

//Hypergeometric test, one tail Fisher exact test
//if both values in a row or column are zero, the p value is 1
double hypergeo_test(int ng, int nl, int cg, int cl){
   if(ng<0||ng<0||cg<0||cl<0)
   {
	  printf("ERROR: hypergeometric test -  all the values must be nonnegative integers!");
	  return -1;
   };
   if(ng+nl==0||cg+cl==0||ng+cg==0||nl+cl==0) return 1.0;

   int N = ng + nl + cg + cl;
   int nm = (ng<nl)?ng:nl;
   int cm = (cg<cl)?cg:cl;
   int n;
  // printf("N.C:G.L %d, %d, %d, %d\n", ng, nl, cg, cl);
   if(nm<cm){
     n = ng + nl ;
     if(ng<nl)//ng is the smallest cell
       return hypergeo_p(n,ng,N,ng+cg, ACCURACY);
     else
       return hypergeo_p(n,nl,N,nl+cl, ACCURACY);
   }
   else{
     n = cg + cl ;
     if(cg<cl)//ng is the smallest cell
       return hypergeo_p(n,cg,N,ng+cg, ACCURACY);
     else
       return hypergeo_p(n,cl,N,nl+cl, ACCURACY);
   }
}

//find the other tail
//p(a+1)/p(a) = b*c/(a+1)(d+1)
// increase a till the cumulated ratio is less than 1 for the first time.
// a should be the smallest cell
int right_tail(int a, int b, int c, int d)
{
	long double da= a, db = b, dc = c, dd = d;
	long double ratio;
	ratio = db*dc/((da+1.0)*(dd+1.0));
	//WARNING, machine precison
	while(ratio-1.0>DBL_EPSILON && db >=0 && dc>=0)
	{
	  da += 1.0;
	  db -= 1.0;
	  dc -= 1.0;
	  dd += 1.0;
	  ratio = ratio*db*dc/((da+1.0)*(dd+1.0));
	}
	return (int)da+1;
}

//Fisher exatc test, two-tailed hypergergeometric test
double fisher_test(int ng, int nl, int cg, int cl){
   if(ng<0||ng<0||cg<0||cl<0)
   {
	  printf("ERROR: Fisher exact test -  all the values must be nonnegative integers!");
	  return -1;
   };
   if(ng+nl==0||cg+cl==0||ng+cg==0||nl+cl==0) return 1.0;
   int N = ng + nl + cg + cl;
   int nm = (ng<nl)?ng:nl;
   int cm = (cg<cl)?cg:cl;
   int n, o;
  // printf("N.C:G.L %d, %d, %d, %d\n", ng, nl, cg, cl);
   if(nm<cm){
     n = ng + nl ;
     if(ng<nl) 
     {//ng is the smallest cell
       o = right_tail(ng, nl, cg, cl);
       if(o<=n)
         return hypergeo_p(n,ng,N,ng+cg, ACCURACY) + hypergeo_p(n,n-o,N,nl+cl, ACCURACY);
       else
         return hypergeo_p(n,ng,N,ng+cg, ACCURACY); 
     }
     else
     {
       o = right_tail(nl, ng, cl, cg);
       if(o<=n)
         return hypergeo_p(n,nl,N,nl+cl, ACCURACY) +  hypergeo_p(n,n-o,N,ng+cg, ACCURACY);
       else
         return  hypergeo_p(n,nl,N,nl+cl, ACCURACY); 
     }
   }
   else
   {
     n = cg + cl ;
     if(cg<cl)//ng is the smallest cell
     {
       o = right_tail(cg, cl, ng, nl);
       if(o<=n)
         return  hypergeo_p(n,cg,N,ng+cg, ACCURACY) + hypergeo_p(n,n-o,N,nl+cl, ACCURACY);
       else
         return  hypergeo_p(n,cg,N,ng+cg, ACCURACY); 
     }
     else 
     {
       o = right_tail(cl, cg, nl, ng);
       if(o<=n)
         return hypergeo_p(n,cl,N,nl+cl, ACCURACY) + hypergeo_p(n,n-o,N,ng+cg, ACCURACY);
       else
         return hypergeo_p(n,cl,N,nl+cl, ACCURACY); 
     }
   }
}

//int main() {
//  printf("%d\t%d\t%d\t%d\t%d\n", 1/2, 2/2, 3/2,4/2, 5/2);
//  printf("%f\t%f\n", bino_p(2,0), bino_p(2, 1));
//  printf("%f\t%f\t%f\n", bino_p(3,0), bino_p(3, 1), bino_p(3, 2));
//  printf("%f\t%f\t%f\n", bino_p(4,0), bino_p(4, 1), bino_p(4, 2));
//  int i;
//  for(i=0;i<20;i++)
//     printf("%d\t%.15g\n", i, bino_p(45,i));
//  printf("\n");
//  printf("HP:%.17g\n", hypergeo_p(10,   5,    100,    50, ACCURACY));
//  printf("HP:%.17g\n", hypergeo_p(100  ,50,   1000,   500, ACCURACY));
//  printf("HP:%.17g\n", hypergeo_p(1000 ,500,  10000,  5000, ACCURACY));
//  printf("HP:%.17g\n", hypergeo_p(10000,5000, 100000, 50000, ACCURACY));
//  printf("HP:%.17g\n", hypergeo_p(9, 1,  23, 11, ACCURACY));
//  printf("HP:%.17g\n", hypergeo_p(14,10, 23, 11, ACCURACY));
//  printf("Test:%.17g\n", hypergeo_test(1,8,10,4));
//  printf("Test:%.17g\n", hypergeo_test(0,8,0,8));
//  printf("Test:%.17g\n", hypergeo_test(9830,8436,10043,8223));
//  printf("Test:%.17g\n", hypergeo_test(18074,1486,17362,2198));
//  printf("Test:%.17g\n", hypergeo_test(15217,3673,15117,3773));
//  printf("Test:%.17g\n", hypergeo_test(3789,14757,3558,14988));
//
//  short a = 0x0;
//  short b = 0x1;
//  short c = 0x2;
//  printf("%d\n", a | b);
//  printf("%d\n", a | c);
//  printf("%d\n", b | c);
//
//  //unsigned int seed = random_init();
//  //printf("seed, %d\n", seed);
//  //for(i=0;i<20;i++)
//  //   printf("%d\t%d\n", i, random_bin());
//}

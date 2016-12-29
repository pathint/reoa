#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <omp.h>
#include <math.h>
#include "stat.h"

#define DATATYPE_VALUE float
#define DATATYPE_SAMPLESIZE unsigned int
#define DATATYPE_GENESIZE unsigned int
#define CONV_THRESHOLD 50
#define MEMBLOCK 1024
#define BIGNINT -1000000
#define MASK_FDR    0b00000001
#define MASK_FILTER 0b00000010
#define MAX_FILES  255
#define MAX_SAMPLESIZE 255
#define DEFAULT_FDR 0.05

//for one sample only
struct pair {
  DATATYPE_GENESIZE h, l;
  DATATYPE_SAMPLESIZE count;
};

//for two samples, without count information, after FDR control or exceptions thresholding
struct pair2 {
  DATATYPE_GENESIZE i, j;
  unsigned char state; //1 for i<j, 0 for i>j; lower bit for sample 1, higer bit for sample 2
  //0b00 = 0, i>j in both samples, 0b11=3, i<j in both samples
  //0b01 = 1, i<j in sample 1, but i>j in sample 2
  //0b10 = 2, i>j in sample 1, but i<j in sample 2
};

/* for RankComp V2.0 */
struct gene_state {
   unsigned char state;//rightmost two bits for in the same or reversed pairs 
   DATATYPE_GENESIZE ng; // control group, greater than the current gene 
   DATATYPE_GENESIZE nl; // control group, less than the current gene
   DATATYPE_GENESIZE cg; // case group, greater than the current gene 
   DATATYPE_GENESIZE cl; // case group, less than the current gene
   double p; //hyergeometric test value
}; 

/* for RankComp (original) */
struct gene_state2 {
   unsigned char state;//rightmost two bits for in the same or reversed pairs
   DATATYPE_GENESIZE ng; // control group, greater than the current gene
   DATATYPE_GENESIZE nl; // control group, less than the current gene
   DATATYPE_GENESIZE g2l; // reversed pair, the number of genes which are greater than the current gene but lower in the case group
   DATATYPE_GENESIZE l2g; // similar as above, but in the reversed order
   double p; //hyergeometric test value
};


/* for Bayesian  */
struct gene_state3 {
   DATATYPE_GENESIZE cg; // concordant pairs, greater than the current gene 
   DATATYPE_GENESIZE cl; // concordant pairs, less than the current gene
   DATATYPE_GENESIZE g2l; // reversed pair, the number of genes which are greater than the current gene but lower in the case group
   DATATYPE_GENESIZE l2g; // similar as above, but in the reversed order
   unsigned char state;
   double p; //hyergeometric test value
};

//new solution: arrays of pointers to struct instead of lists
struct gene_pair
{
  DATATYPE_GENESIZE i, j;
  DATATYPE_SAMPLESIZE nx, cx;
  unsigned char state;
};
//revised on Nov. 8, 2016
//Use a 64-bit int to represent i, j, nx, cx and state
// Bits 0-6 for state, 0-3 for '<' and '>' test, 2-3 for '=' test
// 0 - if   count < threshold in normal group (sample 1)
// 1 - if m-count < threshold in normal group (sample 1)
// 2 - if   count < threshold in case group (sample 2)
// 3 - if m-count < threshold in case group (sample 2)
// 4 - if count0 < threshold in normal group (sample 1)
// 5 - if count0 < threshold in case group (sample 2)
//  6-13, nx, count of <|> (smaller one) in normal group, capacity 2^8-1 = 255
// 14-21, cx, ibid in case group, sample size 255*2=510
// 22-42, gi, gene index, capacity, 2097151 (2097k) gene number
// 43-63, gj, gene index
typedef long unsigned int UB8;

bool   VERBOSE;
double EPSILON; 

struct CHANGE
{
   DATATYPE_GENESIZE n; // number of genes
   float    level;      // change level, Fold-Change (FC)
};

int  less(DATATYPE_VALUE a, DATATYPE_VALUE b);
int equal(DATATYPE_VALUE a, DATATYPE_VALUE b);

//Stable gene-pair ordering for one sample 
int stable_pairs_one(DATATYPE_GENESIZE n, //number of genes
	             DATATYPE_SAMPLESIZE m, //sample size
	             DATATYPE_VALUE data[n*m], //data matrix
	             DATATYPE_SAMPLESIZE max,//threshold, max number of exceptions
	             DATATYPE_GENESIZE max0, //max allowed number of equal pairs
		     int nthreads,
		     struct pair *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     int    *exceptions[nthreads]//for FDR control, number of pairs with the specified number of exceptions
	            );

//Stable gene-pair ordering for two samples;
//Concordant and reversed pairs 
int stable_pairs_two(DATATYPE_GENESIZE n, //number of genes
	             DATATYPE_SAMPLESIZE m1, //sample size
	             DATATYPE_VALUE data1[n*m1], //data matrix
	             DATATYPE_SAMPLESIZE max1,//threshold, max number of exceptions
	             DATATYPE_SAMPLESIZE m2, //sample size
	             DATATYPE_VALUE data2[n*m2], //data matrix
	             DATATYPE_SAMPLESIZE max2,//threshold, max number of exceptions
	             DATATYPE_GENESIZE max0, //max allowed number of equal pairs
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     int    *exceptions1[nthreads],//for FDR control, number of pairs with the specified number of exceptions
		     int    *exceptions2[nthreads]//for FDR control, number of pairs with the specified number of exceptions
	            );

/* Stable gene-pair ordering for one  sample, using UB8 */
int stable_pairs_one2(DATATYPE_GENESIZE n, //number of genes
	             DATATYPE_SAMPLESIZE m, //sample size
	             DATATYPE_VALUE data[n*m], //data matrix
	             DATATYPE_SAMPLESIZE max,//threshold, max number of exceptions
	             DATATYPE_GENESIZE max0, //max allowed number of equal pairs
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     int    *exceptions[nthreads]//for FDR control, number of pairs with the specified number of exceptions
	            );

/* Stable gene-pair ordering for two samples,
 * the second sample is one column, individual sample,
 * stable pairs for the first sample is already given*/
int stable_pairs_ind(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
	             DATATYPE_VALUE column[n] //data matrix
	            );

int filter_gene_dirs(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     struct gene_state *states[n],
		     double alpha, //FDR alpha level for regulation direction
		     int max_cycles,
		     int conv_threshold
	            );

int filter_gene_orig(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     struct gene_state2 *states[n],
		     double alpha, //FDR alpha level for regulation direction
		     int max_cycles,
		     int conv_threshold
	            );

int filter_gene_bayesian(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     struct gene_state3 *states[n]/*,
		     double alpha, //FDR alpha level for regulation direction
		     int max_cycles,
		     int conv_threshold*/
	            );

int  gen_random_sample( DATATYPE_GENESIZE n, //number of genes
	      		DATATYPE_SAMPLESIZE m, //sample size 
	      		DATATYPE_VALUE normal[n*m], //data matrix
	      		DATATYPE_VALUE cases[n*m],  //to be generated
			DATATYPE_GENESIZE n_changes,
			struct CHANGE changes[n_changes], 
			FILE     *out_file
	              );

DATATYPE_GENESIZE similarity_between_two_cols( DATATYPE_GENESIZE n, //number of genes
	      		            DATATYPE_VALUE col1[n], // first column data 
	      		            DATATYPE_VALUE col2[n] // second column
				   );

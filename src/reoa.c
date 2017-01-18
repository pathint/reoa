/**************************************************************************
 *  FSRP: Fast Stably Ranked Pair
 *        identify the stably ranked pairs given a matrix of size n*m,
 *        where n is the gene number and m is the sample size.
 *
 *  Xianlong Wang, Ph.D.
 *  University of Electronic Science and Technology of China.
 *  Email: Wang.Xianlong@139.com
 *
 *  Initialization. Oct. 25, 2016.
 *  Revised on Nov.  8, 2016.
 *  Revised on Nov. 28, 2016.
 **************************************************************************/

#include "reoa.h"

#pragma omp declare simd
int less(DATATYPE_VALUE a, DATATYPE_VALUE b)
{
   if (a<b) return 1;
   else 
   {
	if(fabs(a-b)<EPSILON) return random_bin();
   	else return 0;
   }
}

#pragma omp declare simd
int equal(DATATYPE_VALUE a, DATATYPE_VALUE b)
{
   if (fabs(a-b)<EPSILON) return 1;
   else return 0;
}

#define MEMCHECK(var)								\
	if(var==NULL) {								\
		fprintf(stderr, "# Memory Allocation Error! More Memory May Be Needed!\n");\
		exit(-1);								   \
	}

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
	            )
{
  random_init();

  #pragma omp parallel
  {

    int ithread  = omp_get_thread_num();
    struct pair *new_array;

    int count_memblock = 0;
    count_pairs[ithread] = 0;
    int count, count0;
    DATATYPE_SAMPLESIZE k;
    exceptions[ithread] = malloc(sizeof(int)*(max+1)) ;
    pairs[ithread] = malloc(sizeof(struct pair)*MEMBLOCK); 
    for(k=0;k<max+1;k++)
    	exceptions[ithread][k]=0;

    unsigned char state;
    DATATYPE_GENESIZE i;
    #pragma omp for
    for(i=0;i<n;i++)
    {
      DATATYPE_GENESIZE j;
      for(j=i+1;j<n;j++)
      {
           //fprintf(stderr, "thread=%d, counts=%d\n",ithread, count_pairs[ithread]);
           if(count_pairs[ithread]>=MEMBLOCK*count_memblock)
           {
             count_memblock += 1;
             new_array = realloc(pairs[ithread], sizeof(struct pair)*MEMBLOCK*count_memblock);
             if(new_array) //avoid memory leak
               pairs[ithread] = new_array;
             else
             {
               fprintf(stderr, "# Memory is not allocated on thread %d.\n", ithread);
               exit(-1);
             }
           } //end if

    	   count=0; count0=0;

           #pragma omp simd
           for(k=0;k<m;k++) // count1 [i]<[j]
           	count += less(data[i*m+k], data[j*m+k]);
           #pragma omp simd
           for(k=0;k<m;k++) // count1 [i]<[j]
           	count0 += equal(data[i*m+k], data[j*m+k]);

           state =          ((  count<=max && count0<=max0)?1:0);
           state = state | (((m-count<=max && count0<=max0)?1:0)<<1);
           if(count<=max && count0<=max0) // [i]>[j]
           {
             //insert to the pairs
             pairs[ithread][count_pairs[ithread]].h = i;
             pairs[ithread][count_pairs[ithread]].l = j;
             pairs[ithread][count_pairs[ithread]].count = count;
             count_pairs[ithread] += 1;
             exceptions[ithread][count] += 1;
           } //end if
           else if(m-count<=max && count0<=max0) // [i] <[j]
           {
             //insert to the pairs
             pairs[ithread][count_pairs[ithread]].h = j;
             pairs[ithread][count_pairs[ithread]].l = i;
             pairs[ithread][count_pairs[ithread]].count = m-count;
             count_pairs[ithread] += 1;
    	   //could be improved
             exceptions[ithread][m-count] += 1;
           } //end if
      }//end for j
    }//end for i
  }//end parallel
  //sum all the counts into exceptions[0]

  DATATYPE_SAMPLESIZE k;
  int ithread;
  for(k=0;k<max+1;k++)
     for(ithread=1;ithread<nthreads;ithread++)
        exceptions[0][k] += exceptions[ithread][k];

  return EXIT_SUCCESS;
}


/* Stable gene-pair ordering for two samples */
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
	            )
{
  random_init();

  #pragma omp parallel
  {

    int ithread  = omp_get_thread_num();
    UB8 *new_array, current;

    int count_memblock = 0;
    count_pairs[ithread] = 0;
    unsigned char count1, count2, count01, count02;
    DATATYPE_SAMPLESIZE k;
    exceptions1[ithread] = malloc(sizeof(int)*(max1+1)) ;
    exceptions2[ithread] = malloc(sizeof(int)*(max2+1)) ;
    pairs[ithread] = malloc(sizeof(UB8)*MEMBLOCK); 
    for(k=0;k<max1+1;k++)
    	exceptions1[ithread][k]=0;
    for(k=0;k<max2+1;k++)
    	exceptions2[ithread][k]=0;

    unsigned char state;
    DATATYPE_GENESIZE i;
    #pragma omp for
    for(i=0;i<n;i++)
    {
      DATATYPE_GENESIZE j;
      for(j=i+1;j<n;j++)
      {
           //fprintf(stderr, "thread=%d, counts=%d\n",ithread, count_pairs[ithread]);
           if(count_pairs[ithread]>=MEMBLOCK*count_memblock)
           {
             count_memblock += 1;
             new_array = realloc(pairs[ithread], sizeof(UB8)*MEMBLOCK*count_memblock);
             if(new_array) //avoid memory leak
               pairs[ithread] = new_array;
             else
             {
               fprintf(stderr, "# Memory is not allocated on thread %d.\n", ithread);
               exit(-1);
             }
           } //end if

    	   count1 = 0; count01 = 0;
    	   count2 = 0; count02 = 0;

           #pragma omp simd
           for(k=0;k<m1;k++) // count1 [i]<[j]
           	 count1 += less(data1[i*m1+k], data1[j*m1+k]);
           #pragma omp simd
           for(k=0;k<m1;k++) // count01 [i]==[j]
           	count01 += equal(data1[i*m1+k], data1[j*m1+k]);

           #pragma omp simd
           for(k=0;k<m2;k++) // count2 [i]<[j]
           	 count2 += less(data2[i*m2+k], data2[j*m2+k]);
           #pragma omp simd
           for(k=0;k<m2;k++) // count02 [i]==[j]
           	count02 += equal(data2[i*m2+k], data2[j*m2+k]);

           state =          ((   count1<=max1 && count01<=max0)?1:0);
           state = state | ((((m1-count1)<=max1 && count01<=max0)?1:0)<<1);
           state = state | (((   count2<=max2 && count02<=max0)?1:0)<<2);
           state = state | ((((m2-count2)<=max2 && count02<=max0)?1:0)<<3);

	   if(state!=0) 
	   {
             if (state==1||state==5||state== 9) 
		exceptions1[ithread][   count1] += 1;
             if (state==2||state==6||state==10) 
		exceptions1[ithread][m1-count1] += 1;
             if (state==4||state==5||state== 6) 
		exceptions2[ithread][   count2] += 1;
             if (state==8||state==9||state==10) 
		exceptions2[ithread][m2-count2] += 1;

             state = state | (((count01<=max0)?1:0)<<4);
             state = state | (((count02<=max0)?1:0)<<5);
             //insert to the pairs 
	     current = (UB8) state;
	     current = current | ((UB8) count1 <<6 );
	     current = current | ((UB8) count2 <<14);
	     current = current | ((UB8)      i <<22);
	     current = current | ((UB8)      j <<43);
	     pairs[ithread][count_pairs[ithread]] = current;
	     count_pairs[ithread] += 1;
	   }
      }//end for j
    }//end for i
  }//end parallel

  /* sum all the counts into exceptions[0] */
  DATATYPE_SAMPLESIZE k;
  int ithread;
  for(k=0;k<max1+1;k++)
     for(ithread=1;ithread<nthreads;ithread++)
        exceptions1[0][k] += exceptions1[ithread][k];
  for(k=0;k<max2+1;k++)
     for(ithread=1;ithread<nthreads;ithread++)
        exceptions2[0][k] += exceptions2[ithread][k];

  return EXIT_SUCCESS;
}

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
	            )
{
  random_init();

  #pragma omp parallel
  {

    int ithread  = omp_get_thread_num();
    UB8 *new_array, current;

    int count_memblock = 0;
    count_pairs[ithread] = 0;
    unsigned char count, count0;
    DATATYPE_SAMPLESIZE k;
    exceptions[ithread] = malloc(sizeof(int)*(max+1)) ;
    pairs[ithread] = malloc(sizeof(UB8)*MEMBLOCK); 
    for(k=0;k<max+1;k++)
    	exceptions[ithread][k]=0;

    unsigned char state;
    DATATYPE_GENESIZE i;
    #pragma omp for
    for(i=0;i<n;i++)
    {
      DATATYPE_GENESIZE j;
      for(j=i+1;j<n;j++)
      {
           if(count_pairs[ithread]>=MEMBLOCK*count_memblock)
           {
             count_memblock += 1;
             new_array = realloc(pairs[ithread], sizeof(UB8)*MEMBLOCK*count_memblock);
             if(new_array) //avoid memory leak
               pairs[ithread] = new_array;
             else
             {
               fprintf(stderr, "# Memory is not allocated on thread %d.\n", ithread);
               exit(-1);
             }
           } //end if

    	   count = 0; count0 = 0;

           #pragma omp simd
           for(k=0;k<m;k++) // count [i]<[j]
           	 count += less(data[i*m+k], data[j*m+k]);
           #pragma omp simd
           for(k=0;k<m;k++) // count0 [i]==[j]
           	count0 += equal(data[i*m+k], data[j*m+k]);


           state =          ((    count<=max && count0<=max0)?1:0);
           state = state | ((((m-count)<=max && count0<=max0)?1:0)<<1);

	   if(state!=0) 
	   {
             if (state==1) 
		exceptions[ithread][  count] += 1;
	     else 
		exceptions[ithread][m-count] += 1;

             state = state | (((count0<=max0)?1:0)<<4);
             //insert to the pairs 
	     current = (UB8) state;
	     current = current | ((UB8) count <<6 );
	     current = current | ((UB8)     i <<22);
	     current = current | ((UB8)     j <<43);
	     pairs[ithread][count_pairs[ithread]] = current;
	     count_pairs[ithread] += 1;
	   }
	   else 
	     pairs[ithread][count_pairs[ithread]] = 0;

      }//end for j
    }//end for i
  }//end parallel

  /* sum all the counts into exceptions[0] */
  DATATYPE_SAMPLESIZE k;
  int ithread;
  for(k=0;k<max+1;k++)
     for(ithread=1;ithread<nthreads;ithread++)
        exceptions[0][k] += exceptions[ithread][k];

  return EXIT_SUCCESS;
}

/* Stable gene-pair ordering for two samples,
 * the second sample is one column, individual sample,
 * stable pairs for the first sample is already given*/
int stable_pairs_ind(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
	             DATATYPE_VALUE column[n] //data matrix
	            )
{
  #pragma omp parallel 
  {
    int ithread  = omp_get_thread_num();	      
    unsigned int index, i, j;
    unsigned char state;
    UB8 current;
    for(index=0;index<count_pairs[ithread];index++)
    {
      current = pairs[ithread][index];
      state = ((unsigned char)current) & 0b11;
      i = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
      j = ((unsigned int)(current >>43)) & 0x1FFFFF;
      switch(state)
      {
        case 1: // [i] > [j] in normal
	  if (less(column[i], column[j]) == 1) // [i] < [j] in case
	     state = 9; //1 + 8;
	  else 
	     if (equal(column[i], column[j]) != 1)
	        state = 5; // 1 + 4;
          break;
        case 2:// [i] < [j] in normal
	  if (less(column[i], column[j]) == 1)// [i] < [j] in case
	     state = 10; //2 + 8;
	  else 
	     if (equal(column[i], column[j]) != 1)
	        state = 6; // 2 + 4;
          break;
      } // end switch
      current = (UB8) state;
      current = current | ((UB8) i <<22);
      current = current | ((UB8) j <<43);
      pairs[ithread][index] = current;
     }// end for
  } //end parallel

  return EXIT_SUCCESS;
}

/****************************************************************************
 * Filter dys-regulated genes, judge directions, RanComp V2.0 
 * Input: consistent pairs
 * Output: up and down-regulated genes 
 *****************************************************************************/
#define STATE_INC(i, ngv, nlv, cgv, clv)     		                     	\
	    if(states[i] != NULL){						\
	      states[i]->ng += ngv;   					\
	      states[i]->nl += nlv;   					\
	      states[i]->cg += cgv;   					\
	      states[i]->cl += clv;   					\
	    }									\
	    else {								\
	      states[i] = malloc(sizeof(struct gene_state));		\
	      states[i]->ng = ngv;						\
	      states[i]->nl = nlv;						\
	      states[i]->cg = cgv;						\
	      states[i]->cl = clv;						\
	      states[i]->state = 0;					\
	    }
#define STATE_UPDATE(iv, jv) {							\
	      if (states[j]->state==0)					\
		  states[i]->iv += 1;					\
	      if (states[i]->state==0)					\
		  states[j]->jv += 1;					\
	   }

int filter_gene_dirs(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     struct gene_state *states[n],
		     double alpha, //FDR alpha level for regulation direction
		     int max_cycles,
		     int conv_threshold
	            )
{
  unsigned int i, j;

  //initialization of gene states
  for(i=0; i<nthreads; i++)
    for(j=0; j<count_pairs[i]; j++)
    {
      if (pairs[i][j]!=0)
      {
        UB8 current = pairs[i][j];
        unsigned char state = ((unsigned char)current) & 0b1111;
        unsigned int     c1 = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
        unsigned int     c2 = ((unsigned int)(current >>43)) & 0x1FFFFF;
  
        switch(state){
          case 5: // 1 + 4, [i]>[j] in both 
            STATE_INC(c1, 0, 1, 0, 0)
            STATE_INC(c2, 1, 0, 0, 0)
            STATE_INC(c1, 0, 0, 0, 1)
            STATE_INC(c2, 0, 0, 1, 0)
            break;
          case 10:// 2 + 8, [i]<[j] in both
            STATE_INC(c1, 1, 0, 0, 0)
            STATE_INC(c2, 0, 1, 0, 0)
            STATE_INC(c1, 0, 0, 1, 0)
            STATE_INC(c2, 0, 0, 0, 1)
            break;
          // Stable in the two groups, different direction
          case 6: // 2 + 4
            STATE_INC(c1, 1, 0, 0, 0)
            STATE_INC(c2, 0, 1, 0, 0)
            STATE_INC(c1, 0, 0, 0, 1)
            STATE_INC(c2, 0, 0, 1, 0)
            break;
          case 9: // 1 + 8
            STATE_INC(c1, 0, 1, 0, 0)
            STATE_INC(c2, 1, 0, 0, 0)
            STATE_INC(c1, 0, 0, 1, 0)
            STATE_INC(c2, 0, 0, 0, 1)
            break;
          case 0: //other, remove 
            break;
        } //end switch
      }//end if
    }//end double-for

  //screening for up or down direction using the two-tailed hypergeometric test
  int cycles = 0;
  double *pvalues;
  double p_upper;
  int  gene_up=BIGNINT, gene_down=BIGNINT, gene_flat=BIGNINT;
  int  gene_up_pre, gene_down_pre, gene_flat_pre;

  do {
    int count = 0;
    //two-tailed hypergeometric test
    #pragma omp parallel for reduction(+:count)
    for(i=0;i<n;i++)
      if(states[i]!=NULL) 
      {
         states[i]->p = fisher_test( (int) states[i]->ng, (int) states[i]->nl,
        			     (int) states[i]->cg, (int) states[i]->cl);
        count += 1;
      } //end if

    // copy p values to pvalues
    fprintf(stdout, "#  Cycle %d, %d genes have been tested using the Fisher exact test (aka. two-tail hypergeometric test).\n", cycles, count);
 
    pvalues = malloc(sizeof(double)*count);
    count = 0;
    for(i=0;i<n;i++)
      if(states[i]!=NULL)
          pvalues[count++] = states[i]->p;
 
    //FDR control
    p_upper = bh_threshold(count, pvalues, alpha);
    fprintf(stdout, "#  FDR = %f controlled P value threshold, %.17g. \n", alpha, p_upper);
    free(pvalues);

    //assign states
    gene_up_pre   = gene_up;
    gene_down_pre = gene_down;
    gene_flat_pre = gene_flat;
    gene_up   = 0;
    gene_down = 0;
    gene_flat = 0;
    #pragma omp parallel for reduction(+:gene_up,gene_down,gene_flat)
    for(i=0;i<n;i++)
      if(states[i]!=NULL) {
        if((states[i]->p)<p_upper){ //in case group, g/l ratio is higher, i is down-regulated.
           if(((states[i]->ng+1.0e-5)/(states[i]->nl+1.0e-5)) <
              ((states[i]->cg+1.0e-5)/(states[i]->cl+1.0e-5))) 
	   {
             states[i]->state = 1;
             gene_down += 1;
           }
           else {
             states[i]->state = 2;
             gene_up += 1;
           }
        } 
        else {
          states[i]->state = 0;
          gene_flat += 1;
        }
      } //end if for
    fprintf(stdout, "#  %d genes have been tested as up-regulated.\n", gene_up);
    fprintf(stdout, "#  %d genes have been tested as down-regulated.\n", gene_down);
    fprintf(stdout, "#  %d genes have no particular direction.\n", gene_flat);
 
    //count ng, nl, cg and cl, if gene_state is 0 (flat).
    //clear first
    #pragma omp parallel for
    for(i=0;i<n;i++)
      if(states[i]!=NULL) {
        states[i]->ng = 0;
        states[i]->nl = 0;
        states[i]->cg = 0;
        states[i]->cl = 0;
      }
 
    #pragma omp parallel 
    {
      int ithread  = omp_get_thread_num();	      
      unsigned int index, i, j;
      unsigned char state;
      UB8 current;
      for(index=0;index<count_pairs[ithread];index++)
      {
        current = pairs[ithread][index];
        state = ((unsigned char)current) & 0b1111;
        i = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
        j = ((unsigned int)(current >>43)) & 0x1FFFFF;
        #pragma omp critical
        switch(state){
          case 5: // 1 + 4, [i]>[j]
            STATE_UPDATE(nl, ng)
            STATE_UPDATE(cl, cg)
            break;
          case 10:// 2 + 8, [i]<[j]
            STATE_UPDATE(ng, nl)
            STATE_UPDATE(cg, cl)
            break;
          // Stable in the two groups, different direction
          case 6: // 2 + 4
            STATE_UPDATE(ng, nl)
            STATE_UPDATE(cl, cg)
            break;
          case 9: // 1 + 8
            STATE_UPDATE(nl, ng)
            STATE_UPDATE(cg, cl)
            break;
        } // end switch
       }// end for
    } //end parallel
    cycles +=1;
   } while(cycles<max_cycles && (abs(gene_flat-gene_flat_pre)>conv_threshold||  \
			         abs(gene_up-gene_up_pre)    >conv_threshold||  \
			         abs(gene_down-gene_down_pre)>conv_threshold)); 


  if(abs(gene_flat-gene_flat_pre)<=conv_threshold &&  
     abs(gene_up-gene_up_pre)    <=conv_threshold &&    
     abs(gene_down-gene_down_pre)<=conv_threshold)
    fprintf(stdout, "#  Convergence has been reached after %d cycles.\n", cycles);
  else 
    fprintf(stdout, "#  Max cycles  %d have been exceeded before the convergence.\n", max_cycles);

  return EXIT_SUCCESS;
}//end filter_gene_dirs


/****************************************************************************
 * Filter dys-regulated genes, judge directions, RanComp original 
 * Input: consistent pairs
 * Output: up and down-regulated genes 
 *****************************************************************************/
//Direction up or down?
#define STATE_IN2(i,vi,j,vj) 							\
	  if(states[i] != NULL){						\
	     states[i]->vi += 1;   						\
	  }									\
	  else {							        \
	    states[i] = malloc(sizeof(struct gene_state2));			\
	    states[i]->ng = 0;						\
	    states[i]->nl = 0;						\
	    states[i]->g2l = 0;						\
	    states[i]->l2g = 0;						\
	    states[i]->state = 0;						\
	    states[i]->vi += 1;						\
	  }									\
	  if(states[j] != NULL){						\
	     states[j]->vj += 1;   						\
	  }									\
	  else {								\
	    states[j] = malloc(sizeof(struct gene_state2));			\
	    states[j]->ng = 0;						\
	    states[j]->nl = 0;						\
	    states[j]->g2l = 0;						\
	    states[j]->l2g = 0;						\
	    states[j]->state = 0;						\
	    states[j]->vj += 1;						\
	  }

int filter_gene_orig(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     struct gene_state2 *states[n],
		     double alpha, //FDR alpha level for regulation direction
		     int max_cycles,
		     int conv_threshold
	            )
{

  unsigned int i, j;

  //initialization of gene states
  for(i=0; i<nthreads; i++)
    for(j=0; j<count_pairs[i]; j++)
      if (pairs[i][j]!=0)
      {
        UB8 current = pairs[i][j];
        unsigned char state = ((unsigned char)current) & 0b1111;
        unsigned int     c1 = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
        unsigned int     c2 = ((unsigned int)(current >>43)) & 0x1FFFFF;
  
        switch(state){
          case 5: // 1 + 4, [i] > [j]
	    STATE_IN2(c1, nl, c2, ng)
            break;
          case 10:// 2 + 8, [i] < [j]
	    STATE_IN2(c1, ng, c2, nl)
            break;
          // Stable in the two groups, different direction
          case 6: // 2 + 4
	    STATE_IN2(c1, ng,  c2, nl)
	    STATE_IN2(c1, g2l, c2, l2g)
            break;
          case 9: // 1 + 8
	    STATE_IN2(c1, nl,  c2, ng)
	    STATE_IN2(c1, l2g, c2, g2l)
            break;
          case 0: //other, remove 
            break;
        } //end switch
       } //end if

  int cycles = 0;
  double *pvalues;
  double p_upper;
  int gene_up=BIGNINT, gene_down=BIGNINT, gene_flat=BIGNINT;
  int gene_up_pre, gene_down_pre, gene_flat_pre;
  do {
    int count = 0;
    //two-tailed hypergeometric test
    #pragma omp parallel for reduction(+:count)
    for(i=0;i<n;i++)
      if(states[i]!=NULL)
      {
         states[i]->p =fisher_test((int) states[i]->ng, (int) states[i]->nl,
        			   (int) states[i]->ng + states[i]->l2g - states[i]->g2l, 
        			   (int) states[i]->nl - states[i]->l2g + states[i]->g2l
        			   );
        count += 1;
      } //end if

    // copy p values to pvalues
    printf("#  Cycle %d, %d genes have been tested using the Fisher exact test (aka. two-tail hypergeometric test). \n", cycles, count);
 
    pvalues = malloc(sizeof(double)*count);
    count = 0;
    for(i=0;i<n;i++)
      if(states[i]!=NULL)
        pvalues[count++] = states[i]->p;
 
    //FDR control
    p_upper = bh_threshold(count, pvalues, alpha);
    printf("#  FDR = %f controlled P value threshold, %.17g. \n", alpha, p_upper);
    free(pvalues);
    //assign states
    gene_up_pre   = gene_up;
    gene_down_pre = gene_down;
    gene_flat_pre = gene_flat;
    gene_up   = 0;
    gene_down = 0;
    gene_flat = 0;
    #pragma omp parallel for reduction(+:gene_up,gene_down,gene_flat)
    for(i=0;i<n;i++)
      if(states[i]!=NULL) {
        if((states[i]->p)<p_upper)
	{ //in case group, g/l ratio is higher, i is down-regulated.
	   double a, b, c, d;
	   a = (double) states[i]->ng;
	   b = (double) states[i]->nl;
	   c = (double) states[i]->ng + states[i]->l2g - states[i]->g2l;
	   d = (double) states[i]->nl - states[i]->l2g + states[i]->g2l;

           if(((a+1.0e-5)/(b+1.0e-5)) < ((c+1.0e-5)/(d +1.0e-5))) {
             states[i]->state = 1;
             gene_down += 1;
           }
           else {
             states[i]->state = 2;
             gene_up += 1;
           }
        } 
        else {
          states[i]->state = 0;
          gene_flat += 1;
        }
      } //end if for
    printf("#  %d genes have been tested as up-regulated.\n", gene_up);
    printf("#  %d genes have been tested as down-regulated.\n", gene_down);
    printf("#  %d genes have no particular direction.\n", gene_flat);
 
    //count ng, nl, cg and cl, if gene_state is 0 (flat).
    //clear first
    #pragma omp parallel for
    for(i=0;i<n;i++)
      if(states[i]!=NULL) {
         states[i]->g2l = 0;
         states[i]->l2g = 0;
      }
 
    #pragma omp parallel 
    {
      int ithread  = omp_get_thread_num();	      
      unsigned int index, i, j;
      unsigned char state;
      UB8 current;
      for(index=0;index<count_pairs[ithread];index++)
      {
        current = pairs[ithread][index];
        state = ((unsigned char)current) & 0b1111;
        i = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
        j = ((unsigned int)(current >>43)) & 0x1FFFFF;
        #pragma omp critical
        switch(state)
	{
          // Stable in the two groups, different direction
          case 6: // 2 + 4, reverse, [i]<[j] -> [i]>[j]
	    if (states[i]->state !=2) //if i is up-regulated, i is not counted as g2l for j
		states[j]->l2g += 1;

	    if (states[j]->state !=1) //if j is down-regulated, j is not counted as l2g for i 
	        states[i]->g2l += 1;
            break;
          case 9: // 1 + 8, reverse [i]>[j] -> [i]<[j]
	    if (states[i]->state !=1) //if i is down-regulated, i is not counted as l2g for j
	        states[j]->g2l += 1;

	    if (states[j]->state !=2) //if i is up-regulated, i is not counted as g2l for j
	        states[i]->l2g += 1;
            break;
        } // end switch
      }// end for
     }//end parallel
    cycles +=1;
   } while(cycles<max_cycles && (abs(gene_flat-gene_flat_pre)>conv_threshold||  \
			         abs(gene_up-gene_up_pre)>conv_threshold    ||  \
			         abs(gene_down-gene_down_pre)>conv_threshold)); 

  if(abs(gene_flat-gene_flat_pre)<=conv_threshold &&  
     abs(gene_up-gene_up_pre)    <=conv_threshold &&    
     abs(gene_down-gene_down_pre)<=conv_threshold)
    printf("#  Convergence has been reached after %d cycles.\n", cycles);
  else 
    printf("#  Max cycles  %d have been exceeded before the convergence.\n", max_cycles);

  return EXIT_SUCCESS;
}// end filter_gene_orig


#define STATE_IN3(i, cgv, clv, g2lv, l2gv)     		                     	\
	    if(states[i] != NULL){						\
	      states[i]->cg  += cgv;   					\
	      states[i]->cl  += clv;   					\
	      states[i]->g2l += g2lv;   					\
	      states[i]->l2g += l2gv;   					\
	    }									\
	    else {								\
	      states[i] = malloc(sizeof(struct gene_state3));		\
	      states[i]->cg  = cgv;						\
	      states[i]->cl  = clv;						\
	      states[i]->g2l = g2lv;						\
	      states[i]->l2g = l2gv;						\
	      states[i]->state = 0;					\
	    }

int filter_gene_bayesian(DATATYPE_GENESIZE n, //number of genes
		     int nthreads,
		     UB8 *pairs[nthreads],
		     DATATYPE_GENESIZE count_pairs[nthreads],
		     struct gene_state3 *states[n]/*,
		     double alpha, //FDR alpha level for regulation direction
		     int max_cycles,
		     int conv_threshold*/
	            )
{
  DATATYPE_GENESIZE i, j;

  //initialization of gene states
  for(i=0;i<nthreads;i++)
    for(j=0;j<count_pairs[i];j++)
    {
      if (pairs[i][j]!=0)
      {
        UB8 current = pairs[i][j];
        unsigned char state = ((unsigned char)current) & 0b1111;
        unsigned int     c1 = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
        unsigned int     c2 = ((unsigned int)(current >>43)) & 0x1FFFFF;
  
        switch(state){
          case 5: // 1 + 4, [i]>[j], cl for i, cg for j 
            STATE_IN3(c1, 0, 1, 0, 0)
            STATE_IN3(c2, 1, 0, 0, 0)
            break;
          case 10:// 2 + 8, [i]<[j], cg for i and cl for j
            STATE_IN3(c1, 1, 0, 0, 0)
            STATE_IN3(c2, 0, 1, 0, 0)
            break;
          // Stable in the two groups, different direction
          case 6: // 2 + 4, [i]<[j] to [i]>[j], g2l for i, l2g for j
            STATE_IN3(c1, 0, 0, 1, 0)
            STATE_IN3(c2, 0, 0, 0, 1)
            break;
          case 9: // 1 + 8, [i]>[j] to [i]<[j], l2g for i, g2l for j
            STATE_IN3(c1, 0, 0, 0, 1)
            STATE_IN3(c2, 0, 0, 1, 0)
            break;
          case 0: //other, remove 
            break;
        } //end switch
      }//end if
    }//end double-for

  return EXIT_SUCCESS;
}//end filter_gene_dirs




int  gen_random_sample( DATATYPE_GENESIZE n, //number of genes
	      		DATATYPE_SAMPLESIZE m, //sample size 
	      		DATATYPE_VALUE normal[n*m], //data matrix
	      		DATATYPE_VALUE cases[n*m],  //to be generated
			DATATYPE_GENESIZE n_changes,
			struct CHANGE changes[n_changes], 
			FILE    *out_file
	              )
{  
   DATATYPE_GENESIZE i, total_changes=0, current_group=0, current_sum=0;
   
   for(i=0;i<n_changes;i++) 
      total_changes += changes[i].n;
   
   if(total_changes>n) {
	   fprintf(stderr, "ERROR: Number of changes exceed the number of gene genes.\n");
	   exit(EXIT_FAILURE);
   }
   DATATYPE_GENESIZE *sample;
   sample = malloc(sizeof(DATATYPE_GENESIZE)*total_changes);
   random_sample(n, total_changes, sample);//initialization is done inside the function

   //assign the sample to each category
   DATATYPE_SAMPLESIZE j;
   for(i=0;i<total_changes;i++)
   { 
      if(i==current_sum+changes[current_group].n)
	 current_sum += changes[current_group++].n;
      fprintf(out_file, "%d\t%d\n", current_group, sample[i]);;
      for(j=0;j<m;j++)
	 cases[sample[i]*m+j]=normal[sample[i]*m+j]*changes[current_group].level;
   }//end for "i"

   free(sample);
   return EXIT_SUCCESS;
} // end of gen_random_sample


DATATYPE_GENESIZE similarity_between_two_cols( DATATYPE_GENESIZE n, //number of genes
	      		            DATATYPE_VALUE col1[n], // first column data 
	      		            DATATYPE_VALUE col2[n] // second column
				   )
{
   DATATYPE_GENESIZE i, j, count = 0;
   #pragma omp parallel for reduction(+:count)
   for (i=0; i<n; i++)
     #pragma omp simd
     for (j=i+1; j<n; j++)
	 if (less(col1[i], col1[j]) == less(col2[i], col2[j])) count +=1;

   //return (double) count * 2.0 / (double) (n*(n-1)); 
   return count; 
} //end of similarity_between_two_cols

DATATYPE_GENESIZE similarity_between_three_cols( DATATYPE_GENESIZE n, //number of genes
	      		            DATATYPE_VALUE col1[n], // first column data 
	      		            DATATYPE_VALUE col2[n], // second column
	      		            DATATYPE_VALUE col3[n] // second column
				   )
{
   DATATYPE_GENESIZE i, j, count = 0;
   #pragma omp parallel for reduction(+:count)
   for (i=0; i<n; i++)
     #pragma omp simd
     for (j=i+1; j<n; j++)
	 if (less(col1[i], col1[j]) == less(col2[i], col2[j]) &&
	     less(col2[i], col2[j]) == less(col3[i], col3[j]) ) count +=1;

   //return (double) count * 2.0 / (double) (n*(n-1)); 
   return count; 
} //end of similarity_between_two_cols


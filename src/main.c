
#include "rgpa.h"
#include "main.h"

// could be better
static void display_usage(FILE *f, const char *program_name)
{
     fprintf(f, 
	  "Usage: %s [options] n_files n_genes {data_file...} {sample_size...} [{fdr_level or exception_number...}]\n\n"
	  "Where,\n"
	  "          n_files: the number of input data files\n"
	  "          n_genes: the number of genes in each data file, aka. total row number.\n"
	  "   {data_file...}: list of data file names. All files should have the same number of rows.\n"
	  " {sample_size...}: list of sample sizes (aka. total column number) of the corresponding data files\n"
	  " {fdr_level or\n"
	  "exception_number}: list of FDR level or exception number, which are used to select significantly stable pairs\n"
	  "                   Float numbers between 0. and 1. are taken as FDR levels; otherwise, they are taken as exception numbers\n"
	  "                   If no fdr_level or exception_number is given, the value of --fdr option is also used to select stable pairs. If sample size is less than 3, this argument should be set as 0 only. \n\n"
	  "  e.g. %s [options] 1 20380 normal.dat 32  \n"
	  "       %s [options] 1 20380 normal.dat 32 0.05\n"
	  "       %s [options] 1 20380 normal.dat 32 2\n"
	  "       %s [options] 2 20380 normal.dat cancer.dat 32 26 0.05 0.01\n"
	  "       %s [options] 2 20380 normal.dat cancer.dat 32 26 2 2\n"
	  "       %s [options] 4 20380 normal.dat cancer1.dat cancer2.dat cancer3.dat 32 25 29 17\n\n"
"Function: Select (significantly) stable pairs among one sample set.\n"
"          e.g. to select the stable pairs, respectively, in three data sets, run the program: rgpa -j 0 3 20380 data1.dat data2.dat data3.dat 45 43 48\n"
"          Select concordant and reversed stable pairs between two sample sets.\n"
"          e.g to select the concordant and reversed stable pairs between data sets, run the program: rgpa -j 1 2 20380 data1.dat data2.dat 45 43\n"
"          to use data1.dat as the control group and  select the concordant and reversed stable pairs in data2.dat and data3.dat, respectively, versus data1.dat, run the program:rgpa -j 1 3 20380 data1.dat data2.dat data3.dat 45 43 48. See also --pair option.\n"
"          Identify the dysregulated (up or down-regulated) genes using the RankComp (original and/or V2.0) algorithm.\n"
"          e.g. rgpa -j 2 2 20380 data1.dat data2.dat 45 43. See above example for the application in  multiple groups. \n"
"          Generate simulation random data sets based on a template data set.\n"
"          e.g. rgpa -j 4 -l 10 -c changes.conf -t 5 1 20380 data1.dat 45\n\n"
"          And other miscellaneous functions.\n"
"Limitation: Gene size, number of rows, should be no more than 2097151.\n "
"            Sample size, number or columns, should be no more than 255.\n\n"
	  "\nOptions:\n"
 "                   -h, --help: print this message\n"
 "                -V, --version: print version and copyright information\n"
 "                -v, --verbose: verbose output. Additional information is printed in the result file.\n"
 "     -e, --epsilon = EPSILON : to determine if two values are equal or not. If |a-b|<EPSILON, a and b are equal.\n"
 "                               default: FLT_EPSILON, i.e. 1.192092896e–07\n"
 "                               A reasonable setting should be based on the resolution of the experimental measurements.\n"
 "     -q, --equals = MAXEQUALS: upper limit for equal relations in judging the stable gene pairs, default: 256\n"
 "                               See also --epsilon\n"
 "                               If the equal relation holds (a==b) in more than MAXEQUALS samples (columns), the pair ({a, b}) will NOT be considered as a stable pair. For example if MAXEQUALS is set to 3 and a<=b holds in all 45 samples, but in 4 samples the relation, a==b, holds, the pair {a, b} will NOT be considered as a stable pair if exception_number is 0 (though this condition is  met).\n"
 "        -f, --fdr = FDR_LEVEL: FDR (False Discovery Rate) level, default: 0.05\n"
 "                               This value is set for identifying dysregulated genes. If no fdr_level or exception_number (the 5th group of argument lists, see above) is given, this value is also used for selecting stable gene pairs.\n"
 "          -s, --sample = 0, 1: sample type, default: 0\n"
 "                               0 - cohort\n"
 "                               1 - individual\n"
 "                                   For the individual type samples, the first data set should be the control. Only the first FDR level or exception number (the 5th group of argument list) is applicable to select stable gene pairs.\n"
 "           -j, --job = 0 - 6: job type, default: 2\n"
 "                               0 - select stable gene pairs only. This job requires at least one data set.\n"
 "                               1 - select concordant and reversed gene pairs. This job requires at least two data sets.\n"
 "                               2 - identify dysregulated genes (with directions). job=1 is a pre-step for this job. \n"
 "                               3 - generate simulation disease data set using a normal data set as template\n"
 "                                   this requires --changes to be set. This job requires one data set and one fold-change (FC) setting file. Additional data sets will be omitted.\n"
 "                               4 - generate simulation disease data set and filter dysregulated genes. This job combines job=3 and job=2\n"
 "                               5 - calculate concordance scores of gene pairs between two samples (columns) within one data set. This job requires at least one data set.\n"
 "                                   The concordance score is defined as the number of the same-ordered pairs divided by the total number of pairs.\n"
 "                               6 - calculate concordance scores among three samples (columns) within one data set. This job requires at least one data set.\n"
 "         -p, --pair = 0, 1, 2: comparison looping mode for calculations applied to multiple groups of data sets, default: 0\n"
 "                               0 - one vs many, the first data set is used as the control, i.e. 0v1, 0v2, 0v3...\n"
 "                               1 - many vs many,  all possible pair-wise combinations, i.e. 0v1, 0v2, 0v3..., 1v2, 1v3...\n" 
 "                               2 - consecutive pair comparison, i.e. 0v1, 2v3, 4v5...\n" 
 "    -c, --changes = file_name: configuration file for simulation data set generation, see --job = 3 or 4\n"
 "                               Each line contains two numbers: an integer followed by a real number \n"
 "                               the integer is the number of genes to be modified and the real number is\n"
 "                               the fold-change (FC) level\n. The number of lines should be given by --lines."
 "            -l, --lines = NUM: number of record lines in the change configuration file, see --changes\n"
 "          -t, --times = TIMES: number of simulation data set to be generated, default: 1\n"
 "                               TIMES should be a positive integer. The simulation data set will be generated for TIMES times. See --job = 3 or 4\n"
 "    -a, --algorithm = 0, 1, 2: choice of algorithm, default: 0\n"
 "                               0 - RankComp V2.0 only\n"
 "                               1 - RankComp (original) only\n"
 "                               2 - Both the original and the new V2.0 algorithms\n"
 "    -m, --cycles = MAX_CYCLES: max cycles for filtering dysregulated genes, default: 128\n"
 "-o, --cycles_orig= MAX_CYCLES: max cycles used in the original RankComp algorithm, default: 2 (no filtering iterations)\n"
 "-n, --convergence = THRESHOLD: convergence threshold (difference in number of dysregulated genes between two consecutive cycles), default: 50\n"
 "    -M, --memory = PERCENTAGE: a real number between 0 and 1 to control memory usage (expert only). A small value reduces the  memory usage and speeds up the calculation. However, if the value is too small, the actual FDR level may be unnecessarily lower than the control FDR level. Default: automatically determined value.\n\n"
	  "Please cite the following paper when you publish results from this program.\n"
	  "Hongwei Wang, Qiang Sun, Wenyuan Zhao, Lishuang Qi, Yunyan Gu, Pengfei Li, Mengmeng Zhang, "
	  "Yang Li, Shu-Lin Liu, and Zheng Guo. (2015). Individual-level analysis of differential "
	  "expression of genes and pathways for personalized medicine. Bioinformatics. 31(1):62-8. DOI:10.1093/bioinformatics/btu522\n\n"
	  , program_name, program_name, program_name, program_name, program_name, program_name, program_name);
} //end display_usage

static void display_version(FILE *f, const char *program_name)
{
     fprintf(f, 
	  "Program: %s Version 0.1\n"
	  "\n"
	  "Copyright (C) 2016  by Xianlong Wang.\n\n"
	  "This program is free software; you can redistribute it and/or modify it "
	  "under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as published by the Free Software Foundation . "
	  "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; "
	  "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "
	  "See the GNU General Public License for more details.\n\n"
	  "Please cite the following paper when you publish results from this program.\n"
	  "Hongwei Wang, Qiang Sun, Wenyuan Zhao, Lishuang Qi, Yunyan Gu, Pengfei Li, Mengmeng Zhang, "
	  "Yang Li, Shu-Lin Liu, and Zheng Guo. (2015). Individual-level analysis of differential "
	  "expression of genes and pathways for personalized medicine. Bioinformatics. 31(1):62-8. DOI:10.1093/bioinformatics/btu522\n\n"
	  , program_name);
} //end display_version

int read_changes(const char *file_name, int n_changes, struct CHANGE changes[n_changes])
{
   FILE *file;
   file = fopen(file_name, "r");
   int  status;
   int i;
   for (i=0; i<n_changes; i++)
   {
      status = fscanf(file, "%d", &changes[i].n);
      status = fscanf(file, "%g", &changes[i].level);
   }
   return EXIT_SUCCESS;
}

int read_data(char *file_name, int n_genes, int sample_size, DATATYPE_VALUE *data)
{
   int j, k, status;
   FILE *file = fopen(file_name, "r");

   if (file == NULL) 
   {
      fprintf(stderr, "ERROR: file %s cannot be open to read.\n", file_name);
      exit(EXIT_FAILURE);
   }

   for (j=0; j<n_genes; j++)
     for (k=0; k<sample_size; k++)
	 status = fscanf(file, "%g", &data[j*sample_size+k]);
   fclose(file);
   return EXIT_SUCCESS;
}

int extract_column(unsigned int n, unsigned int m, DATATYPE_VALUE data[n*m], 
		   unsigned int c, DATATYPE_VALUE column[n])
{
   unsigned int i;
   for (i=0; i<n; i++)
      column[i] = data[i*m+c];
   return EXIT_SUCCESS;
}
char *my_itoa(int value, char *str)
{
   int len = snprintf( NULL, 0, "%d", value);  
   str = realloc(str, sizeof(char)*(len + 1));
   snprintf(str, len + 1, "%d", value);
   return str;
}

char *make_file_name(char *stem, int index, char *ext, char *fname)
{
   int len = snprintf( NULL, 0, "%s%d%s", stem, index, ext);  
   fname = realloc(fname, sizeof(char)*(len + 1));
   snprintf(fname, len + 1, "%s%d%s", stem, index, ext);
   return fname;
}

char *make_file_name2(char *stem, int i1, char *in, int i2, char *ext, char *fn)
{
   int len = snprintf( NULL, 0, "%s%d%s%d%s", stem, i1, in, i2, ext);  
   fn = realloc(fn, sizeof(char)*(len + 1));
   snprintf(fn, len + 1, "%s%d%s%d%s", stem, i1, in, i2, ext);
   return fn;
}

char *make_file_name3(char *stem, int i1, char *in1, int i2, char *in2, int i3, char *ext, char *fn)
{
   int len = snprintf( NULL, 0, "%s%d%s%d%s%d%s", stem, i1, in1, i2, in2, i3, ext);  
   fn = realloc(fn, sizeof(char)*(len + 1));
   snprintf(fn, len + 1, "%s%d%s%d%s%d%s", stem, i1, in1, i2, in2, i3, ext);
   return fn;
}

/* Job type 0: select significantly stable gene pairs and save into files */
/* Input: Data Files, FDR Levels (or Max Exception Number), max_equals
 * Output: stable gene pairs 
 */
int select_stable_pairs(int               n_files, 
		        char      *files[n_files],
			int sample_sizes[n_files],
			float        fdr[n_files],
			int              n_genes ,
			int            max_equals
		       ) 
{
  int i;

  /* load data into memory */
  for (i=0; i<n_files; i++)
  {
     if (VERBOSE) fprintf(stdout, "#Processing data file %d...\n", i);
     /* load data into memory */
     int             size_current = sample_sizes[i];
     char           *file_current = files[i];
     DATATYPE_VALUE *data_current = malloc( sizeof(DATATYPE_VALUE) * n_genes * size_current);
     if (data_current == NULL) 
     {
	fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i, file_current);
	break;
	//exit(EXIT_FAILURE);
     }
     read_data(file_current, n_genes, size_current, data_current);

     /* max exceptions */
     float fdr_current = fdr[i];
     int mode = 0, max_exceptions;
     if (fdr_current >=1.|| fdr_current == 0.) 
     {  //selection based on max exceptions
	max_exceptions = (int) fdr_current;
	mode = 1;
     }
     else
     {
	max_exceptions = (size_current - 1)/2;//WARNING: +1 is more reasonable here to cover all possibities, 
	                                      // but not okay in the next function, contagious problem
        if (MEMORY_USAGE < 0)
	   max_exceptions = bino_ctrl(size_current, fdr_current);
	else
	   max_exceptions = MEMORY_USAGE * max_exceptions;
	mode = 0;
     } //end if...else...

     int n_threads;
     #pragma omp parallel shared(n_threads)
     {
       #pragma omp single
       n_threads =  omp_get_num_threads(); 
     }

     struct pair **pairs;
     int    **exceptions;
     DATATYPE_GENESIZE *count_pairs;
     pairs      = malloc( sizeof(struct pair *) * n_threads); 
     exceptions = malloc( sizeof(int *) * n_threads); 
     count_pairs= malloc( sizeof(int *) * n_threads); 
     stable_pairs_one(n_genes, size_current, data_current, max_exceptions, max_equals, 
		      n_threads, pairs, count_pairs, exceptions);

     struct FDR result;
     /* FDR control mode*/
     if(mode==0)
     {
       result = bh_ctrl(size_current, max_exceptions, exceptions[0], n_genes*(n_genes-1)/2 ,fdr_current);
       max_exceptions = result.index;
     }
     else
     {
       result = bh_eval(size_current, max_exceptions, exceptions[0], n_genes*(n_genes-1)/2);
       result.p = bino_p(size_current, result.index);
     }

     /* print ouf statistics */
     fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i);
     fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr_current);
     fprintf(stdout, "#  FDR actual level: %g\n", result.value);
     fprintf(stdout, "#  Max. exception numbers: %d\n", max_exceptions);
     fprintf(stdout, "#  Max. binomial test P value: %g\n", result.p);
     int j, k;
     for(j=0; j<=max_exceptions; j++)
        fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions[0][j]);

     char *out_name;
     FILE *out_file;
     out_name = malloc(sizeof(char *)*10);
     out_name = make_file_name(STABLE, i, ".dat", out_name);
     out_file = fopen(out_name, "w");
     if (out_file == NULL) 
     {
        out_name = "stderr";
	out_file = stderr;
     }
     fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file:%s\n", out_name);
     free(out_name);
     for (j=0; j<n_threads; j++)
	for (k=0; k<count_pairs[j]; k++)
	{ 
	   struct pair pair_current = pairs[j][k];
          if (pair_current.count <= max_exceptions)
	  {
	    fprintf(out_file, "%d\t%d", pair_current.h, pair_current.l);
	    if (VERBOSE)
	       fprintf(out_file, "\t%d\t%.6g", pair_current.count, bino_p(size_current, pair_current.count));
	    fprintf(out_file, "\n");
	  }
	}
     /* clean up */
     if (out_file != stderr) fclose(out_file);
     for (j=0; j<n_threads; j++) 
     {
        free(pairs[j]);
        free(exceptions[j]);
     }
     free(pairs);
     free(count_pairs);
     free(exceptions);
     free(data_current);

  }//end for

  return EXIT_SUCCESS;

} // end of select_stable_pairs

/******************************************************************************
 * Job types 1 and 2: select concordant and reversed genes pairs
 *  Input: Data Files, FDR Levels (or Max Exception Number), max_equals, mode
 * Output: concordant and reversed gene pairs 
 ******************************************************************************/
int select_consistent_pairs(int           n_files, 
		        char      *files[n_files],
			int sample_sizes[n_files],
			float        fdr[n_files],
			int              n_genes ,
			int            max_equals,
			char		     job,
			char	       pair_mode,
			float              alpha, 
			int           max_cycles, 
			int           ori_cycles, 
			int        max_threshold,
			char           algorithm
		       ) 
{
  int i1, i2, end1, step1;;
  switch (pair_mode)
  {
    case 0:
       end1 = 1;
      step1 = 1;
      break;

    case 1:
       end1 = n_files -1;
      step1 = 1;
      break;

    case 2:
       end1 = n_files -1;
      step1 = 2;
      break;

  }//end switch 

  for (i1 =0; i1<end1; i1 += step1)
  {
     /* load data1 into memory */
     int             size1 = sample_sizes[i1];
     char           *file1 =        files[i1];
     DATATYPE_VALUE *data1 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size1);
     if (data1 == NULL) 
     {
        fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i1, file1);
        break;
        //exit(EXIT_FAILURE);
     }
     read_data(file1, n_genes, size1, data1);

     /* max exceptions */
     float fdr1 = fdr[i1];
     int  mode1 = 0, max1;
     if (fdr1 >=1.|| fdr1 == 0.) 
     {  //selection based on max exceptions
         max1 = (int) fdr1;
	 max1 = max1<0?0:max1;
        mode1 = 1;
     }
     else
     {
         max1 = (size1 - 1)/2;//WARNING: boundary issue
	 if (MEMORY_USAGE < 0)
	    max1 = (int) bino_ctrl(size1, fdr1 );
	 else
	    max1 = MEMORY_USAGE * max1;  
        mode1 = 0;
     } //end if...else...

     /*loop2 control*/
     int begin2, step2;
     switch (pair_mode)
     {
       case 0:
       case 1:
        begin2 = i1 + 1;
         step2 = 1;
         break;

       case 2:
        begin2 = 1;
         step2 = 2;
         break;

     }//end switch 
     /* loop over the rest of the data files */
     for (i2=begin2; i2<n_files; i2 += step2)
     {
        if (VERBOSE) fprintf(stdout, "#Processing data file %d...\n", i2);
        /* load data2 into memory */
        int             size2 = sample_sizes[i2];
        char           *file2 = files[i2];
        DATATYPE_VALUE *data2 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size2);
        if (data2 == NULL) 
        {
           fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i2, file2);
           break;
           //exit(EXIT_FAILURE);
        }
        read_data(file2, n_genes, size2, data2);

        /* max exceptions */
        float fdr2 = fdr[i2];
        int  mode2 = 0, max2;
        if (fdr2 >=1.|| fdr2 == 0.) 
        {  //selection based on max exceptions
            max2 = (int) fdr2;
	    max2 = max2<0?0:max2;
           mode2 = 1;
        }
        else
        {
            max2 = (size2 - 1)/2;//WARNING: boundary issue
	    if (MEMORY_USAGE < 0)
	       max2 = bino_ctrl(size2, fdr2 );
	    else
	       max2 = MEMORY_USAGE * max2;  
           mode2 = 0;
        } //end if...else...

        int n_threads;
        #pragma omp parallel shared(n_threads)
        {
          #pragma omp single
          n_threads =  omp_get_num_threads(); 
        }

        UB8 **pairs;
        int    **exceptions1, **exceptions2;
        DATATYPE_GENESIZE *count_pairs;
        pairs       = malloc( sizeof(UB8 *) * n_threads); 
        exceptions1 = malloc( sizeof(int *) * n_threads); 
        exceptions2 = malloc( sizeof(int *) * n_threads); 
        count_pairs = malloc( sizeof(int *) * n_threads); 

	stable_pairs_two(n_genes, size1, data1, max1, size2, data2, max2, max_equals,
		         n_threads, pairs, count_pairs, exceptions1, exceptions2);

        struct FDR result1, result2;
        int j, k;
        char *out_name1, *out_name2, *out_name3, *out_name4;
        FILE *out_file1, *out_file2, *out_file3, *out_file4;
	if (i1==0)
	{
          /* FDR control mode*/
          if (mode1==0)
          {
            result1 = bh_ctrl(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2 ,fdr1);
               max1 = result1.index;
          }
          else
          {
            result1  = bh_eval(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2);
           result1.p =  bino_p(size1, result1.index);
          }
          /* print ouf statistics */
          fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i1);
          fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr1);
          fprintf(stdout, "#  FDR actual level: %g\n", result1.value);
          fprintf(stdout, "#  Max. exception numbers: %d\n", max1);
          fprintf(stdout, "#  Max. binomial test P value: %g\n", result1.p);
          for(j=0; j<=max1; j++)
             fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions1[0][j]);

	  /* open file to write stable pairs */
	  if (job == 1 || VERBOSE)
	  {
             out_name1 = malloc(sizeof(char *)*10);
             out_name1 = make_file_name(STABLE, i1, EXT, out_name1);
             out_file1 = fopen(out_name1, "w");
             if (out_file1 == NULL) 
             {
                out_name1 = "stderr";
                out_file1 = stderr;
             }
             fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file: %s\n", out_name1);
             free(out_name1);
	  }
	}

        /* FDR control mode*/
        if (mode2==0)
        {
          result2 = bh_ctrl(size2, max2, exceptions2[0], n_genes*(n_genes-1)/2 ,fdr2);
             max2 = result2.index;
        }
        else
        {
          result2  = bh_eval(size2, max2, exceptions2[0], n_genes*(n_genes-1)/2);
         result2.p =  bino_p(size2, result2.index);
        }

        /* print ouf statistics */
        fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i2);
        fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr2);
        fprintf(stdout, "#  FDR actual level: %g\n", result2.value);
        fprintf(stdout, "#  Max. exception numbers: %d\n", max2);
        fprintf(stdout, "#  Max. binomial test P value: %g\n", result2.p);
        for(j=0; j<=max2; j++)
           fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions2[0][j]);

	/* open file to write stable pairs */
	if (job == 1 ||VERBOSE)
	{
           out_name2 = malloc(sizeof(char *)*10);
           out_name2 = make_file_name(STABLE, i2, EXT, out_name2);
           out_file2 = fopen(out_name2, "w");
           if (out_file2 == NULL) 
           {
              out_name2 = "stderr";
              out_file2 = stderr;
           }
           fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file: %s\n", out_name2);
           free(out_name2);

	   /* open file to write consistent pairs */
           out_name3 = malloc(sizeof(char *)*10);
           out_name3 = make_file_name2(CONCORDANT, i2, "_", i1, EXT, out_name3);
           out_file3 = fopen(out_name3, "w");
           if (out_file3 == NULL) 
           {
              out_name3 = "stderr";
              out_file3 = stderr;
           }
           fprintf(stdout, "#  Writing concordant pairs (high, low) to file: %s\n", out_name3);
           free(out_name3);

           out_name4 = malloc(sizeof(char *)*10);
           out_name4 = make_file_name2(REVERSED, i2, "_", i1, EXT, out_name4);
           out_file4 = fopen(out_name4, "w");
           if (out_file4 == NULL) 
           {
              out_name4 = "stderr";
              out_file4 = stderr;
           }
           fprintf(stdout, "#  Writing reversed pairs (high, low) [in the first data set]  to file: %s\n", out_name4);
           free(out_name4);
	}

	int count_same = 0, count_reverse = 0;

	/* write stable pairs to files */
        for (j=0; j<n_threads; j++)
           for (k=0; k<count_pairs[j]; k++)
           { 
             UB8 current = pairs[j][k];
	     unsigned char state;
	     unsigned char ne = ((unsigned char)(current >>4)) & 0b00000001;
	     unsigned char ce = ((unsigned char)(current >>5)) & 0b00000001;;
	     unsigned char c1 = (unsigned char)(current >> 6); 
	     unsigned char c2 = (unsigned char)(current >> 14); 
             state =          ((      c1<=max1 && ne == 1)?1:0);
             state = state | (((size1-c1<=max1 && ne == 1)?1:0)<<1);
             state = state | (((      c2<=max2 && ce == 1)?1:0)<<2);
             state = state | (((size2-c2<=max2 && ce == 1)?1:0)<<3);

             unsigned int  vi = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
             unsigned int  vj = ((unsigned int)(current >>43)) & 0x1FFFFF;

	     if (job == 1 ||VERBOSE)
	     {
               /* the first */
	       if (i1==0)
	       {
                  if (state==1||state==5||state==9)  // 1+{0, 4, 8}
	          {
	             fprintf(out_file1, "%d\t%d", vj, vi);
                     if (VERBOSE)
                        fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                     fprintf(out_file1, "\n");
	          }
	          else if (state==2||state==6||state==10) // 2+{0, 4, 8}
	          {
	             c1 = size1 - c1;
	             fprintf(out_file1, "%d\t%d", vj, vi);
                     if (VERBOSE)
                        fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                     fprintf(out_file1, "\n");
	          }
	       }//end if i1==0

               /* the rest */
               if (state==4||state==5||state==6)  // 4 +{0,1,2}
	       {
	           fprintf(out_file2, "%d\t%d", vj, vi);
                   if (VERBOSE)
                      fprintf(out_file2, "\t%d\t%.6g", c2, bino_p(size2, c2));
                   fprintf(out_file2, "\n");
	        }
	        else if (state==8||state==9||state==10) // 8+{0, 1, 2} 
	        {
	           c2 = size2 - c2;
	           fprintf(out_file2, "%d\t%d", vj, vi);
                   if (VERBOSE)
                      fprintf(out_file2, "\t%d\t%.6g", c2, bino_p(size2, c2));
                   fprintf(out_file2, "\n");
	        }
	     }

	     /* keep only the consistent pairs */
	     switch(state)
	     {
	        case  5: // 1+4
                  if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vi, vj);
		  count_same += 1;
	          break;
	        case 10: // 2+8
                  if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vj, vi);
		  count_same += 1;
	          break;
	        case  6: // 2+4
                  if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vj, vi);
		  count_reverse += 1;
	          break;
	        case  9: // 1+8
                  if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vi, vj);
		  count_reverse += 1;
	          break;
		default:
		  pairs[j][k] = 0;//clear the other data
	     }//end switch

           }//end for
	fprintf(stdout, "#  Number of concordant pairs between data set %d and %d: %d\n", i2, i1, count_same);
	fprintf(stdout, "#  Number of reversed pairs between data set %d and %d: %d\n", i2, i1, count_reverse);
        /* clean up */
	if (job==1||VERBOSE) 
	{
           if (i1==0 && out_file1 != stderr) fclose(out_file1);
           if (out_file2 != stderr) fclose(out_file2);
           if (out_file3 != stderr) fclose(out_file3);
           if (out_file4 != stderr) fclose(out_file4);
	}

        for (j=0; j<n_threads; j++) 
        {
           //free(pairs[j]);
           free(exceptions1[j]);
           free(exceptions2[j]);
        }
        //free(pairs);
        //free(count_pairs);
        free(exceptions1);
        free(exceptions2);
        free(data2);

	if (job==1) {
           for (j=0; j<n_threads; j++) 
              free(pairs[j]);
           free(pairs);
           free(count_pairs);
	   continue;
	}
	/* filter dysregulated genes */
	if (algorithm == 0||algorithm ==2)
	{
          //state for each gene, not pairs
          struct gene_state **genes;//an array of pointers to struct 
          //allocate an array of n pointers to struct gene_states
          genes = (struct gene_state **) malloc(sizeof(struct gene_state *)*n_genes);
          //initialization
          for (j=0; j<n_genes; j++) genes[j] = (struct gene_state *) NULL;
	  //call
	  fprintf(stdout, "#  Filter dysregulated genes...\n");
	  filter_gene_dirs(n_genes, n_threads, pairs, count_pairs, 
	  		 genes, alpha, max_cycles, max_threshold);

	  /* write dysregulated genes into files*/
	  /* open file to write stable pairs */
          out_name3 = malloc(sizeof(char *)*10);
          out_name3 = make_file_name2(UP, i2, "_", i1, EXT, out_name3);
          out_file3 = fopen(out_name3, "w");
          if (out_file3 == NULL) 
          {
             out_name3 = "stderr";
             out_file3 = stderr;
          }
          fprintf(stdout, "#  Writing up-regulated genes into file: %s\n", out_name3);
          free(out_name3);

          out_name4 = malloc(sizeof(char *)*10);
          out_name4 = make_file_name2(DOWN, i2, "_", i1, EXT, out_name4);
          out_file4 = fopen(out_name4, "w");
          if (out_file4 == NULL) 
          {
             out_name4 = "stderr";
             out_file4 = stderr;
          }
          fprintf(stdout, "#  Writing down-regulated genes into file: %s\n", out_name4);
          free(out_name4);
          for(j=0;j<n_genes;j++)
             if(genes[j]!=NULL) 
             {
	         switch(genes[j]->state)
	         {
	  	  case 2:
                       fprintf(out_file3, "%d", j);    
		       if (VERBOSE) fprintf(out_file3, "\t%.4g", genes[j]->p);
		       //if (genes[j]->p <1) fprintf(out_file3, "\t%d\t%d\t%d\t%d\t%g\t%g", 
			//	       genes[j]->ng, genes[j]->nl, genes[j]->cg, genes[j]->cl, genes[j]->p,
			//	       fisher_test(genes[j]->ng, genes[j]->nl, genes[j]->cg, genes[j]->cl));
		       fprintf(out_file3, "\n");
	  	     break;
	  	  case 1:
                       fprintf(out_file4, "%d", j);    
		       if (VERBOSE) fprintf(out_file4, "\t%.4g", genes[j]->p);
		       //if (genes[j]->p <1) fprintf(out_file4, "\t%d\t%d\t%d\t%d\t%g\t%g", 
			//	       genes[j]->ng, genes[j]->nl, genes[j]->cg, genes[j]->cl, genes[j]->p,
			//	       fisher_test(genes[j]->ng, genes[j]->nl, genes[j]->cg, genes[j]->cl));
		       fprintf(out_file4, "\n");
	  	     break;
	         }//end switch
                 free(genes[j]);
             }
          if (out_file3 != stderr) fclose(out_file3);
          if (out_file4 != stderr) fclose(out_file4);
	  free(genes);
	}//end if algo.

	/* filter dysregulated genes */
	if (algorithm == 1||algorithm ==2)
	{
          //state for each gene, not pairs
          struct gene_state2 **genes;//an array of pointers to struct 
          //allocate an array of n pointers to struct gene_states
          genes = (struct gene_state2 **) malloc(sizeof(struct gene_state2 *)*n_genes);
          //initialization
          for (j=0; j<n_genes; j++) genes[j] = (struct gene_state2 *) NULL;
	  //call
	  fprintf(stdout, "#  Filter dysregulated genes with the old RankComp algo...\n");
	  filter_gene_orig(n_genes, n_threads, pairs, count_pairs, 
	  		   genes, alpha, ori_cycles, max_threshold);

	  /* write dysregulated genes into files*/
	  /* open file to write stable pairs */
          out_name3 = malloc(sizeof(char *)*10);
          out_name3 = make_file_name2(UP, i2, "v", i1, EXT, out_name3);
          out_file3 = fopen(out_name3, "w");
          if (out_file3 == NULL) 
          {
             out_name3 = "stderr";
             out_file3 = stderr;
          }
          fprintf(stdout, "#  Writing up-regulated genes into file: %s\n", out_name3);
          free(out_name3);

          out_name4 = malloc(sizeof(char *)*10);
          out_name4 = make_file_name2(DOWN, i2, "v", i1, EXT, out_name4);
          out_file4 = fopen(out_name4, "w");
          if (out_file4 == NULL) 
          {
             out_name4 = "stderr";
             out_file4 = stderr;
          }
          fprintf(stdout, "#  Writing down-regulated genes into file: %s\n", out_name4);
          free(out_name4);
          for(j=0;j<n_genes;j++)
             if(genes[j]!=NULL) 
             {
	         switch(genes[j]->state)
	         {
	  	  case 2:
                       fprintf(out_file3, "%d", j);    
		       if (VERBOSE) fprintf(out_file3, "\t%.4g", genes[j]->p);
		       fprintf(out_file3, "\n");
	  	     break;
	  	  case 1:
                       fprintf(out_file4, "%d", j);    
		       if (VERBOSE) fprintf(out_file4, "\t%.4g", genes[j]->p);
		       fprintf(out_file4, "\n");
	  	     break;
	         }//end switch
                 free(genes[j]);
             }
          if (out_file3 != stderr) fclose(out_file3);
          if (out_file4 != stderr) fclose(out_file4);
	  free(genes);
	}//end if algo 

	/* clean up*/
        for (j=0; j<n_threads; j++) 
           free(pairs[j]);
        free(pairs);
        free(count_pairs);

     }//end inner for
     /*clean up*/
     free(data1);

  }//end outer for
  return EXIT_SUCCESS;

} // end of select_stable_pairs


/******************************************************************************
 * Sample Type 1: individual samples, each column is a sample
 * Job types 1 and 2: select concordant and reversed genes pairs
 *  Input: Data Files, FDR Levels (or Max Exception Number), max_equals, mode
 * Output: concordant and reversed gene pairs 
 ******************************************************************************/
int select_consistent_pairs_ind(int           n_files, 
		        char      *files[n_files],
			int sample_sizes[n_files],
			float        fdr[n_files],
			int              n_genes ,
			int            max_equals,
			char		     job,
			char	       pair_mode,
			float              alpha, 
			int           max_cycles, 
			int           ori_cycles, 
			int        max_threshold,
			char           algorithm
		       ) 
{
  int i1, i2, end1, step1;
  switch (pair_mode)
  {
    case 0:
       end1 = 1;
      step1 = 1;
      break;

    case 1:
       end1 = n_files -1;
      step1 = 1;
      break;

    case 2:
       end1 = n_files -1;
      step1 = 2;
      break;

  }//end switch 

  for (i1 =0; i1<end1; i1 += step1)
  {
     /* load data1 into memory */
     int             size1 = sample_sizes[i1];
     char           *file1 =        files[i1];
     DATATYPE_VALUE *data1 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size1);
     if (data1 == NULL) 
     {
        fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i1, file1);
        break;
        //exit(EXIT_FAILURE);
     }
     read_data(file1, n_genes, size1, data1);

     /* max exceptions */
     float fdr1 = fdr[i1];
     int  mode1 = 0, max1;
     if (fdr1 >=1.|| fdr1 == 0.) 
     {  //selection based on max exceptions
         max1 = (int) fdr1;
	 max1 = max1<0?0:max1;
        mode1 = 1;
     }
     else
     {
         max1 = (size1 - 1)/2;//WARNING: boundary issue
        mode1 = 0;
     } //end if...else...

     /* select stable pairs in sample i1 */
     int n_threads;
     #pragma omp parallel shared(n_threads)
     {
       #pragma omp single
       n_threads =  omp_get_num_threads(); 
     }

     UB8 **pairs;
     int **exceptions1;
     DATATYPE_GENESIZE *count_pairs;
     pairs       = malloc( sizeof(UB8 *) * n_threads); 
     exceptions1 = malloc( sizeof(int *) * n_threads); 
     count_pairs = malloc( sizeof(int *) * n_threads); 

     stable_pairs_one2(n_genes, size1, data1, max1, max_equals, 
		       n_threads, pairs, count_pairs, exceptions1);

     struct FDR result1;
     int j, k;

     /* FDR control and print out stable pairs */
     /* FDR control mode*/
     if (mode1==0)
     {
       result1 = bh_ctrl(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2 ,fdr1);
          max1 = result1.index;
     }
     else
     {
       result1  = bh_eval(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2);
      result1.p =  bino_p(size1, result1.index);
     }
     /* print ouf statistics */
     fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i1);
     fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr1);
     fprintf(stdout, "#  FDR actual level: %g\n", result1.value);
     fprintf(stdout, "#  Max. exception numbers: %d\n", max1);
     fprintf(stdout, "#  Max. binomial test P value: %g\n", result1.p);
     for(j=0; j<=max1; j++)
        fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions1[0][j]);

     char *out_name1, *out_name3, *out_name4;
     FILE *out_file1, *out_file3, *out_file4;
     if (job==1||VERBOSE)
     {
        /* open file to write stable pairs */
        out_name1 = malloc(sizeof(char *)*10);
        out_name1 = make_file_name(STABLE, i1, EXT, out_name1);
        out_file1 = fopen(out_name1, "w");
        if (out_file1 == NULL) 
        {
           out_name1 = "stderr";
           out_file1 = stderr;
        }
        fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file: %s\n", out_name1);
        free(out_name1);
     }

     /* write stable pairs to files */
     for (j=0; j<n_threads; j++)
        for (k=0; k<count_pairs[j]; k++)
         { 
            UB8 current = pairs[j][k];
            if (current == 0) continue;

	    unsigned char state;
	    unsigned char e1 = ((unsigned char)(current >>4)) & 0b00000001;
	    unsigned char c1 = (unsigned char)(current >> 6); 

            state =          ((      c1<=max1 && e1 == 1)?1:0);
            state = state | (((size1-c1<=max1 && e1 == 1)?1:0)<<1);

            unsigned int  vi = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
            unsigned int  vj = ((unsigned int)(current >>43)) & 0x1FFFFF;

            if (job==1||VERBOSE)
	    {
               if (state==1)  // i<j
	       {
	          fprintf(out_file1, "%d\t%d", vj, vi);
                  if (VERBOSE)
                     fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                  fprintf(out_file1, "\n");
	       }
	       else if (state==2) // i> j
	       {
	          c1 = size1 - c1;
	          fprintf(out_file1, "%d\t%d", vi, vj);
                  if (VERBOSE)
                     fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                  fprintf(out_file1, "\n");
	       }
	    }

	    if (state==0)
	       current = 0;
	    else
	    {
	       current = (UB8) state;
	       current = current | ((UB8) vi <<22);
	       current = current | ((UB8) vj <<43);
	    }
	    pairs[j][k] = current;

           }//end for

     if (job==1||VERBOSE)
        if (out_file1 != stderr) fclose(out_file1);


     /*loop2 control*/
     int begin2, step2;
     switch (pair_mode)
     {
       case 0:
       case 1:
        begin2 = i1 + 1;
         step2 = 1;
         break;

       case 2:
        begin2 = 1;
         step2 = 2;
         break;

     }//end switch 
     /* loop over the rest of the data files */
     for (i2=begin2; i2<n_files; i2 += step2)
     {
        if (VERBOSE) fprintf(stdout, "#Processing data file %d...\n", i2);
        /* load data2 into memory */
        int             size2 = sample_sizes[i2];
        char           *file2 = files[i2];
        DATATYPE_VALUE *data2 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size2);
        if (data2 == NULL) 
        {
           fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i2, file2);
           break;
           //exit(EXIT_FAILURE);
        }
        read_data(file2, n_genes, size2, data2);

	int i3;
	DATATYPE_VALUE *column;
	column = malloc( sizeof(DATATYPE_VALUE) * n_genes);

	for (i3=0; i3<size2; i3++)
	{
           fprintf(stdout, "# Processing column %d...\n", i3);

	   /* column i3*/
	   extract_column(n_genes, size2, data2, i3, column);

           /* consistent pairs */
	   stable_pairs_ind(n_genes, n_threads, pairs, count_pairs, column);

	   if (job==1||VERBOSE)
	   {
	      /* write out consistent pairs */
	      /* open file to write stable pairs */
              out_name3 = malloc(sizeof(char *)*10);
              out_name3 = make_file_name3(CONCORDANT, i2, "c", i3, "_", i1, EXT, out_name3);
              out_file3 = fopen(out_name3, "w");
              if (out_file3 == NULL) 
              {
                 out_name3 = "stderr";
                 out_file3 = stderr;
              }
              fprintf(stdout, "#  Writing concordant pairs (high, low) to file: %s\n", out_name3);
              free(out_name3);

              out_name4 = malloc(sizeof(char *)*10);
              out_name4 = make_file_name3(REVERSED, i2, "c", i3, "_", i1, EXT, out_name4);
              out_file4 = fopen(out_name4, "w");
              if (out_file4 == NULL) 
              {
                 out_name4 = "stderr";
                 out_file4 = stderr;
              }
              fprintf(stdout, "#  Writing reversed pairs (high, low) [in the first data set]  to file: %s\n", out_name4);
              free(out_name4);
	   }

	   int count_same = 0, count_reverse = 0;

	   /* write stable pairs to files */
           for (j=0; j<n_threads; j++)
              for (k=0; k<count_pairs[j]; k++)
              { 
                UB8 current = pairs[j][k];
	        unsigned char state = ((unsigned char)current) & 0b1111;
                unsigned int  vi = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
                unsigned int  vj = ((unsigned int)(current >>43)) & 0x1FFFFF;

	        /* keep only the consistent pairs */
	        switch(state)
	        {
	           case  5: // 1+4
                     if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vi, vj);
	   	     count_same += 1;
	             break;
	           case 10: // 2+8
                     if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vj, vi);
	   	     count_same += 1;
	             break;
	           case  6: // 2+4
                     if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vj, vi);
	   	     count_reverse += 1;
	             break;
	           case  9: // 1+8
                     if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vi, vj);
	   	     count_reverse += 1;
	             break;
	        }//end switch
              }//end for
	   fprintf(stdout, "#  Number of concordant pairs between data set %d column %d and %d: %d\n", i2, i3,  i1, count_same);
	   fprintf(stdout, "#  Number of reversed pairs between data set %d column %d and %d: %d\n", i2, i3, i1, count_reverse);
           /* clean up */

	   if (job==1 || VERBOSE)
	   {
             if (out_file3 != stderr) fclose(out_file3);
             if (out_file4 != stderr) fclose(out_file4);
	   }

	   /* filter dysregulated genes */
	   if (algorithm == 0||algorithm ==2)
	   {
             //state for each gene, not pairs
             struct gene_state **genes;//an array of pointers to struct 
             //allocate an array of n pointers to struct gene_states
             genes = (struct gene_state **) malloc(sizeof(struct gene_state *)*n_genes);
             //initialization
             for (j=0; j<n_genes; j++) genes[j] = (struct gene_state *) NULL;
	     //call
	     fprintf(stdout, "#  Filter dysregulated genes...\n");
	     filter_gene_dirs(n_genes, n_threads, pairs, count_pairs, 
	     		 genes, alpha, max_cycles, max_threshold);

	     /* write dysregulated genes into files*/
	     /* open file to write stable pairs */
             out_name3 = malloc(sizeof(char *)*10);
             out_name3 = make_file_name3(UP, i2, "c", i3, "_", i1, EXT, out_name3);
             out_file3 = fopen(out_name3, "w");
             if (out_file3 == NULL) 
             {
                out_name3 = "stderr";
                out_file3 = stderr;
             }
             fprintf(stdout, "#  Writing up-regulated genes into file: %s\n", out_name3);
             free(out_name3);

             out_name4 = malloc(sizeof(char *)*10);
             out_name4 = make_file_name3(DOWN, i2, "c", i3, "_", i1, EXT, out_name4);
             out_file4 = fopen(out_name4, "w");
             if (out_file4 == NULL) 
             {
                out_name4 = "stderr";
                out_file4 = stderr;
             }
             fprintf(stdout, "#  Writing down-regulated genes into file: %s\n", out_name4);
             free(out_name4);
             for(j=0;j<n_genes;j++)
                if(genes[j]!=NULL) 
                {
	            switch(genes[j]->state)
	            {
	     	      case 2:
                         fprintf(out_file3, "%d", j);    
		         if (VERBOSE) fprintf(out_file3, "\t%.4g", genes[j]->p);
		         fprintf(out_file3, "\n");
	     	         break;
	     	       case 1:
                          fprintf(out_file4, "%d", j);    
		          if (VERBOSE) fprintf(out_file4, "\t%.4g", genes[j]->p);
		           fprintf(out_file4, "\n");
	     	          break;
	            }//end switch
                    free(genes[j]);
                }
             if (out_file3 != stderr) fclose(out_file3);
             if (out_file4 != stderr) fclose(out_file4);
	     free(genes);
	   }//end if algo.

	   /* filter dysregulated genes */
	   if (algorithm == 1||algorithm ==2)
	   {
             //state for each gene, not pairs
             struct gene_state2 **genes;//an array of pointers to struct 
             //allocate an array of n pointers to struct gene_states
             genes = (struct gene_state2 **) malloc(sizeof(struct gene_state2 *)*n_genes);
             //initialization
             for (j=0; j<n_genes; j++) genes[j] = (struct gene_state2 *) NULL;
	     //call
	     fprintf(stdout, "#  Filter dysregulated genes with the old RankComp algo...\n");
	     filter_gene_orig(n_genes, n_threads, pairs, count_pairs, 
	     		   genes, alpha, ori_cycles, max_threshold);

	     /* write dysregulated genes into files*/
	     /* open file to write stable pairs */
             out_name3 = malloc(sizeof(char *)*10);
             out_name3 = make_file_name3(UP, i2, "c", i3, "v", i1, EXT, out_name3);
             out_file3 = fopen(out_name3, "w");
             if (out_file3 == NULL) 
             {
                out_name3 = "stderr";
                out_file3 = stderr;
             }
             fprintf(stdout, "#  Writing up-regulated genes into file: %s\n", out_name3);
             free(out_name3);

             out_name4 = malloc(sizeof(char *)*10);
             out_name4 = make_file_name3(DOWN, i2, "c",i3, "v", i1, EXT, out_name4);
             out_file4 = fopen(out_name4, "w");
             if (out_file4 == NULL) 
             {
                out_name4 = "stderr";
                out_file4 = stderr;
             }
             fprintf(stdout, "#  Writing down-regulated genes into file: %s\n", out_name4);
             free(out_name4);
             for(j=0;j<n_genes;j++)
                if(genes[j]!=NULL) 
                {
	            switch(genes[j]->state)
	            {
	     	     case 2:
                       fprintf(out_file3, "%d", j);    
		       if (VERBOSE) fprintf(out_file3, "\t%.4g", genes[j]->p);
		       fprintf(out_file3, "\n");
	     	       break;
	     	      case 1:
                       fprintf(out_file4, "%d", j);    
		       if (VERBOSE) fprintf(out_file4, "\t%.4g", genes[j]->p);
		       fprintf(out_file4, "\n");
	     	       break;
	            }//end switch
                    free(genes[j]);
                }
             if (out_file3 != stderr) fclose(out_file3);
             if (out_file4 != stderr) fclose(out_file4);
	     free(genes);
	   }//end if algo 

	}//end for i3

     }//end inner for
     /*clean up*/
     /* clean up*/
     for (j=0; j<n_threads; j++) 
        free(pairs[j]);
     free(pairs);
     free(count_pairs);
     free(data1);

  }//end outer for
  return EXIT_SUCCESS;

} // end of select_stable_pairs_ind

/* New Output Format */
int select_consistent_pairs_ind2(int           n_files, 
		        char      *files[n_files],
			int sample_sizes[n_files],
			float        fdr[n_files],
			int              n_genes ,
			int            max_equals,
			char		     job,
			char	       pair_mode,
			float              alpha, 
			int           max_cycles, 
			int           ori_cycles, 
			int        max_threshold,
			char           algorithm
		       ) 
{
  int i1, i2, end1, step1;
  switch (pair_mode)
  {
    case 0:
       end1 = 1;
      step1 = 1;
      break;

    case 1:
       end1 = n_files -1;
      step1 = 1;
      break;

    case 2:
       end1 = n_files -1;
      step1 = 2;
      break;

  }//end switch 

  for (i1 =0; i1<end1; i1 += step1)
  {
     /* load data1 into memory */
     int             size1 = sample_sizes[i1];
     char           *file1 =        files[i1];
     DATATYPE_VALUE *data1 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size1);
     if (data1 == NULL) 
     {
        fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i1, file1);
        break;
        //exit(EXIT_FAILURE);
     }
     read_data(file1, n_genes, size1, data1);

     /* max exceptions */
     float fdr1 = fdr[i1];
     int  mode1 = 0, max1;
     if (fdr1 >=1.|| fdr1 == 0.) 
     {  //selection based on max exceptions
         max1 = (int) fdr1;
	 max1 = max1<0?0:max1;
        mode1 = 1;
     }
     else
     {
         max1 = (size1 - 1)/2;//WARNING: boundary issue
        mode1 = 0;
     } //end if...else...

     /* select stable pairs in sample i1 */
     int n_threads;
     #pragma omp parallel shared(n_threads)
     {
       #pragma omp single
       n_threads =  omp_get_num_threads(); 
     }

     UB8 **pairs;
     int **exceptions1;
     DATATYPE_GENESIZE *count_pairs;
     pairs       = malloc( sizeof(UB8 *) * n_threads); 
     exceptions1 = malloc( sizeof(int *) * n_threads); 
     count_pairs = malloc( sizeof(int *) * n_threads); 

     stable_pairs_one2(n_genes, size1, data1, max1, max_equals, 
		       n_threads, pairs, count_pairs, exceptions1);

     struct FDR result1;
     int j, k;

     /* FDR control and print out stable pairs */
     /* FDR control mode*/
     if (mode1==0)
     {
       result1 = bh_ctrl(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2 ,fdr1);
          max1 = result1.index;
     }
     else
     {
       result1  = bh_eval(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2);
      result1.p =  bino_p(size1, result1.index);
     }
     /* print ouf statistics */
     fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i1);
     fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr1);
     fprintf(stdout, "#  FDR actual level: %g\n", result1.value);
     fprintf(stdout, "#  Max. exception numbers: %d\n", max1);
     fprintf(stdout, "#  Max. binomial test P value: %g\n", result1.p);
     for(j=0; j<=max1; j++)
        fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions1[0][j]);

     char *out_name1, *out_name3, *out_name4;
     FILE *out_file1, *out_file3, *out_file4;
     if (job==1||VERBOSE)
     {
        /* open file to write stable pairs */
        out_name1 = malloc(sizeof(char *)*10);
        out_name1 = make_file_name(STABLE, i1, EXT, out_name1);
        out_file1 = fopen(out_name1, "w");
        if (out_file1 == NULL) 
        {
           out_name1 = "stderr";
           out_file1 = stderr;
        }
        fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file: %s\n", out_name1);
        free(out_name1);
     }

     /* write stable pairs to files */
     for (j=0; j<n_threads; j++)
        for (k=0; k<count_pairs[j]; k++)
         { 
            UB8 current = pairs[j][k];
            if (current == 0) continue;

	    unsigned char state;
	    unsigned char e1 = ((unsigned char)(current >>4)) & 0b00000001;
	    unsigned char c1 = (unsigned char)(current >> 6); 

            state =          ((      c1<=max1 && e1 == 1)?1:0);
            state = state | (((size1-c1<=max1 && e1 == 1)?1:0)<<1);

            unsigned int  vi = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
            unsigned int  vj = ((unsigned int)(current >>43)) & 0x1FFFFF;

            if (job==1||VERBOSE)
	    {
               if (state==1)  // i<j
	       {
	          fprintf(out_file1, "%d\t%d", vj, vi);
                  if (VERBOSE)
                     fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                  fprintf(out_file1, "\n");
	       }
	       else if (state==2) // i> j
	       {
	          c1 = size1 - c1;
	          fprintf(out_file1, "%d\t%d", vi, vj);
                  if (VERBOSE)
                     fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                  fprintf(out_file1, "\n");
	       }
	    }

	    if (state==0)
	       current = 0;
	    else
	    {
	       current = (UB8) state;
	       current = current | ((UB8) vi <<22);
	       current = current | ((UB8) vj <<43);
	    }
	    pairs[j][k] = current;

           }//end for

     if (job==1||VERBOSE)
        if (out_file1 != stderr) fclose(out_file1);

     /*loop2 control*/
     int begin2, step2;
     switch (pair_mode)
     {
       case 0:
       case 1:
        begin2 = i1 + 1;
         step2 = 1;
         break;

       case 2:
        begin2 = 1;
         step2 = 2;
         break;

     }//end switch 
     /* loop over the rest of the data files */
     for (i2=begin2; i2<n_files; i2 += step2)
     {
        if (VERBOSE) fprintf(stdout, "#Processing data file %d...\n", i2);
        /* load data2 into memory */
        int             size2 = sample_sizes[i2];
        char           *file2 = files[i2];
        DATATYPE_VALUE *data2 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size2);
        if (data2 == NULL) 
        {
           fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i2, file2);
           break;
           //exit(EXIT_FAILURE);
        }
        read_data(file2, n_genes, size2, data2);

	int i3;
	DATATYPE_VALUE *column;
	column = malloc( sizeof(DATATYPE_VALUE) * n_genes);

	/* table to store all gene directions */
	// bits 0,1 for state in the RandComp V2.0 algo. 00=0, flat; 10=2, up; 01=1, down 
	// bits 2,3 for state in the original algo. 00xx, flat; 10xx=4+, up; 01xx=4+, down
	unsigned char *genes_result;
	genes_result = malloc( sizeof(unsigned char) * n_genes * size2);
	for (j=0; j<n_genes; j++)
	   for(k=0; k<size2; k++) genes_result[j*size2+k] = 0; 

	for (i3=0; i3<size2; i3++)
	{
           fprintf(stdout, "# Processing column %d...\n", i3);

	   /* column i3*/
	   extract_column(n_genes, size2, data2, i3, column);

           /* consistent pairs */
	   stable_pairs_ind(n_genes, n_threads, pairs, count_pairs, column);

	   if (job==1||VERBOSE)
	   {
	      /* write out consistent pairs */
	      /* open file to write stable pairs */
              out_name3 = malloc(sizeof(char *)*10);
              out_name3 = make_file_name3(CONCORDANT, i2, "c", i3, "_", i1, EXT, out_name3);
              out_file3 = fopen(out_name3, "w");
              if (out_file3 == NULL) 
              {
                 out_name3 = "stderr";
                 out_file3 = stderr;
              }
              fprintf(stdout, "#  Writing concordant pairs (high, low) to file: %s\n", out_name3);
              free(out_name3);

              out_name4 = malloc(sizeof(char *)*10);
              out_name4 = make_file_name3(REVERSED, i2, "c", i3, "_", i1, EXT, out_name4);
              out_file4 = fopen(out_name4, "w");
              if (out_file4 == NULL) 
              {
                 out_name4 = "stderr";
                 out_file4 = stderr;
              }
              fprintf(stdout, "#  Writing reversed pairs (high, low) [in the first data set]  to file: %s\n", out_name4);
              free(out_name4);
	   }

	   int count_same = 0, count_reverse = 0;

	   /* write stable pairs to files */
           for (j=0; j<n_threads; j++)
              for (k=0; k<count_pairs[j]; k++)
              { 
                UB8 current = pairs[j][k];
	        unsigned char state = ((unsigned char)current) & 0b1111;
                unsigned int  vi = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
                unsigned int  vj = ((unsigned int)(current >>43)) & 0x1FFFFF;

	        /* keep only the consistent pairs */
	        switch(state)
	        {
	           case  5: // 1+4
                     if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vi, vj);
	   	     count_same += 1;
	             break;
	           case 10: // 2+8
                     if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vj, vi);
	   	     count_same += 1;
	             break;
	           case  6: // 2+4
                     if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vj, vi);
	   	     count_reverse += 1;
	             break;
	           case  9: // 1+8
                     if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vi, vj);
	   	     count_reverse += 1;
	             break;
	        }//end switch
              }//end for
	   fprintf(stdout, "#  Number of concordant pairs between data set %d column %d and %d: %d\n", i2, i3,  i1, count_same);
	   fprintf(stdout, "#  Number of reversed pairs between data set %d column %d and %d: %d\n", i2, i3, i1, count_reverse);
           /* clean up */

	   if (job==1 || VERBOSE)
	   {
             if (out_file3 != stderr) fclose(out_file3);
             if (out_file4 != stderr) fclose(out_file4);
	   }

	   /* filter dysregulated genes */
	   if (algorithm == 0||algorithm ==2)
	   {
             //state for each gene, not pairs
             struct gene_state **genes;//an array of pointers to struct 
             //allocate an array of n pointers to struct gene_states
             genes = (struct gene_state **) malloc(sizeof(struct gene_state *)*n_genes);
             //initialization
             for (j=0; j<n_genes; j++) genes[j] = (struct gene_state *) NULL;
	     //call
	     fprintf(stdout, "#  Filter dysregulated genes...\n");
	     filter_gene_dirs(n_genes, n_threads, pairs, count_pairs, 
	     		 genes, alpha, max_cycles, max_threshold);

	     /* write dysregulated genes into files*/
             for(j=0;j<n_genes;j++)
                if(genes[j]!=NULL) 
                {
                  genes_result[j*size2+i3] = genes[j]->state;
                  free(genes[j]);
                }
	     free(genes);
	   }//end if algo.

	   /* filter dysregulated genes */
	   if (algorithm == 1||algorithm ==2)
	   {
             //state for each gene, not pairs
             struct gene_state2 **genes;//an array of pointers to struct 
             //allocate an array of n pointers to struct gene_states
             genes = (struct gene_state2 **) malloc(sizeof(struct gene_state2 *)*n_genes);
             //initialization
             for (j=0; j<n_genes; j++) genes[j] = (struct gene_state2 *) NULL;
	     //call
	     fprintf(stdout, "#  Filter dysregulated genes with the old RankComp algo...\n");
	     filter_gene_orig(n_genes, n_threads, pairs, count_pairs, 
	     		   genes, alpha, ori_cycles, max_threshold);

	     /* write dysregulated genes into files*/
             for(j=0;j<n_genes;j++)
                if(genes[j]!=NULL) 
                {
	            genes_result[j*size2+i3] = genes_result[j*size2+i3] | ((genes[j]->state) << 2);
                    free(genes[j]);
                }
	     free(genes);
	   }//end if algo 

	}//end for i3

	/*Print out genes_result */
        /* open file to write dysregulated genes */
        out_name1 = malloc(sizeof(char *)*10);
        out_name1 = make_file_name(GENE_STATE, i2, EXT, out_name1);
        out_file1 = fopen(out_name1, "w");
        if (out_file1 == NULL) 
        {
           out_name1 = "stderr";
           out_file1 = stderr;
        }
        fprintf(stdout, "#  Writing gene states to file: %s\n", out_name1);
        fprintf(stdout, "#  Gene state code: rightmost 1st and 2nd  two bits for results given by the default RankComp 2.0 algorithm.\n");
        fprintf(stdout, "#                   the 3rd and 4th bits for results given by the original RankComp algorithm.\n");
        fprintf(stdout, "#                   00 for flat genes, 01 for down-regulated genes, 10 for up-regulated genes.\n");
        fprintf(stdout, "#  Example: 0, flat by both algo. 1, down given by the RankComp 2.0 algo but flat by the original algo.\n");
        fprintf(stdout, "#           5, down given by both algo. 10, up given by both algo.\n");
        free(out_name1);
	for (j=0; j<n_genes; j++)
	{
	  for(k=0; k<size2; k++) 
		  fprintf(out_file1, "%d\t", genes_result[j*size2+k]); 
	  fprintf(out_file1, "\n");
	}
	if (out_file1 != stderr) fclose(out_file1);


     }//end inner for
     /*clean up*/
     /* clean up*/
     for (j=0; j<n_threads; j++) 
        free(pairs[j]);
     free(pairs);
     free(count_pairs);
     free(data1);

  }//end outer for
  return EXIT_SUCCESS;

} // end of select_stable_pairs_ind

/******************************************************************************
 * Job type 7: new algo. based on the Bayesian idea
 *  Input: Data Files, FDR Levels (or Max Exception Number), max_equals, mode
 * Output: concordant and reversed gene pairs 
 ******************************************************************************/
int select_genes__bayesian(int           n_files, 
		        char      *files[n_files],
			int sample_sizes[n_files],
			float        fdr[n_files],
			int              n_genes ,
			int            max_equals,
			char		     job,
			char	       pair_mode/*,
			float              alpha, 
			int           max_cycles, 
			int           ori_cycles, 
			int        max_threshold,
			char           algorithm*/
		       ) 
{
  int i1, i2, end1, step1;
  switch (pair_mode)
  {
    case 0:
       end1 = 1;
      step1 = 1;
      break;

    case 1:
       end1 = n_files -1;
      step1 = 1;
      break;

    case 2:
       end1 = n_files -1;
      step1 = 2;
      break;

  }//end switch 

  for (i1 =0; i1<end1; i1 += step1)
  {
     /* load data1 into memory */
     int             size1 = sample_sizes[i1];
     char           *file1 =        files[i1];
     DATATYPE_VALUE *data1 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size1);
     if (data1 == NULL) 
     {
        fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i1, file1);
        break;
        //exit(EXIT_FAILURE);
     }
     read_data(file1, n_genes, size1, data1);

     /* max exceptions */
     float fdr1 = fdr[i1];
     int  mode1 = 0, max1;
     if (fdr1 >=1.|| fdr1 == 0.) 
     {  //selection based on max exceptions
         max1 = (int) fdr1;
	 max1 = max1<0?0:max1;
        mode1 = 1;
     }
     else
     {
         max1 = (size1 - 1)/2;//WARNING: boundary issue
	 if (MEMORY_USAGE < 0)
	    max1 = (int) bino_ctrl(size1, fdr1 );
	 else
	    max1 = MEMORY_USAGE * max1;  
        mode1 = 0;
     } //end if...else...

     /*loop2 control*/
     int begin2, step2;
     switch (pair_mode)
     {
       case 0:
       case 1:
        begin2 = i1 + 1;
         step2 = 1;
         break;

       case 2:
        begin2 = 1;
         step2 = 2;
         break;

     }//end switch 
     /* loop over the rest of the data files */
     for (i2=begin2; i2<n_files; i2 += step2)
     {
        if (VERBOSE) fprintf(stdout, "#Processing data file %d...\n", i2);
        /* load data2 into memory */
        int             size2 = sample_sizes[i2];
        char           *file2 = files[i2];
        DATATYPE_VALUE *data2 = malloc( sizeof(DATATYPE_VALUE) * n_genes * size2);
        if (data2 == NULL) 
        {
           fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i2, file2);
           break;
           //exit(EXIT_FAILURE);
        }
        read_data(file2, n_genes, size2, data2);

        /* max exceptions */
        float fdr2 = fdr[i2];
        int  mode2 = 0, max2;
        if (fdr2 >=1.|| fdr2 == 0.) 
        {  //selection based on max exceptions
            max2 = (int) fdr2;
	    max2 = max2<0?0:max2;
           mode2 = 1;
        }
        else
        {
            max2 = (size2 - 1)/2;//WARNING: boundary issue
	    if (MEMORY_USAGE < 0)
	       max2 = bino_ctrl(size2, fdr2 );
	    else
	       max2 = MEMORY_USAGE * max2;  
           mode2 = 0;
        } //end if...else...

        int n_threads;
        #pragma omp parallel shared(n_threads)
        {
          #pragma omp single
          n_threads =  omp_get_num_threads(); 
        }

        UB8 **pairs;
        int    **exceptions1, **exceptions2;
        DATATYPE_GENESIZE *count_pairs;
        pairs       = malloc( sizeof(UB8 *) * n_threads); 
        exceptions1 = malloc( sizeof(int *) * n_threads); 
        exceptions2 = malloc( sizeof(int *) * n_threads); 
        count_pairs = malloc( sizeof(int *) * n_threads); 

	stable_pairs_two(n_genes, size1, data1, max1, size2, data2, max2, max_equals,
		         n_threads, pairs, count_pairs, exceptions1, exceptions2);

        struct FDR result1, result2;
        unsigned int j, k;
        char *out_name1, *out_name2, *out_name3, *out_name4;
        FILE *out_file1, *out_file2, *out_file3, *out_file4;
	if (i1==0)
	{
          /* FDR control mode*/
          if (mode1==0)
          {
            result1 = bh_ctrl(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2 ,fdr1);
               max1 = result1.index;
          }
          else
          {
            result1  = bh_eval(size1, max1, exceptions1[0], n_genes*(n_genes-1)/2);
           result1.p =  bino_p(size1, result1.index);
          }
          /* print ouf statistics */
          fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i1);
          fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr1);
          fprintf(stdout, "#  FDR actual level: %g\n", result1.value);
          fprintf(stdout, "#  Max. exception numbers: %d\n", max1);
          fprintf(stdout, "#  Max. binomial test P value: %g\n", result1.p);
          for(j=0; j<=max1; j++)
             fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions1[0][j]);

	  /* open file to write stable pairs */
	  if (job == 1 || VERBOSE)
	  {
             out_name1 = malloc(sizeof(char *)*10);
             out_name1 = make_file_name(STABLE, i1, EXT, out_name1);
             out_file1 = fopen(out_name1, "w");
             if (out_file1 == NULL) 
             {
                out_name1 = "stderr";
                out_file1 = stderr;
             }
             fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file: %s\n", out_name1);
             free(out_name1);
	  }
	}//end if i1==0

        /* FDR control mode*/
        if (mode2==0)
        {
          result2 = bh_ctrl(size2, max2, exceptions2[0], n_genes*(n_genes-1)/2 ,fdr2);
             max2 = result2.index;
        }
        else
        {
          result2  = bh_eval(size2, max2, exceptions2[0], n_genes*(n_genes-1)/2);
         result2.p =  bino_p(size2, result2.index);
        }

        /* print ouf statistics */
        fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i2);
        fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr2);
        fprintf(stdout, "#  FDR actual level: %g\n", result2.value);
        fprintf(stdout, "#  Max. exception numbers: %d\n", max2);
        fprintf(stdout, "#  Max. binomial test P value: %g\n", result2.p);
        for(j=0; j<=max2; j++)
           fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions2[0][j]);

	/* open file to write stable pairs */
	if (job == 1 ||VERBOSE)
	{
           out_name2 = malloc(sizeof(char *)*10);
           out_name2 = make_file_name(STABLE, i2, EXT, out_name2);
           out_file2 = fopen(out_name2, "w");
           if (out_file2 == NULL) 
           {
              out_name2 = "stderr";
              out_file2 = stderr;
           }
           fprintf(stdout, "#  Writing stable pairs (high, low, [ exceptions, P]) to file: %s\n", out_name2);
           free(out_name2);

	   /* open file to write consistent pairs */
           out_name3 = malloc(sizeof(char *)*10);
           out_name3 = make_file_name2(CONCORDANT, i2, "_", i1, EXT, out_name3);
           out_file3 = fopen(out_name3, "w");
           if (out_file3 == NULL) 
           {
              out_name3 = "stderr";
              out_file3 = stderr;
           }
           fprintf(stdout, "#  Writing concordant pairs (high, low) to file: %s\n", out_name3);
           free(out_name3);

           out_name4 = malloc(sizeof(char *)*10);
           out_name4 = make_file_name2(REVERSED, i2, "_", i1, EXT, out_name4);
           out_file4 = fopen(out_name4, "w");
           if (out_file4 == NULL) 
           {
              out_name4 = "stderr";
              out_file4 = stderr;
           }
           fprintf(stdout, "#  Writing reversed pairs (high, low) [in the first data set]  to file: %s\n", out_name4);
           free(out_name4);
	}//end if 

	int count_same = 0, count_reverse = 0;

	/* write stable pairs to files */
        for (j=0; j<n_threads; j++)
           for (k=0; k<count_pairs[j]; k++)
           { 
             UB8 current = pairs[j][k];
	     unsigned char state;
	     unsigned char ne = ((unsigned char)(current >>4)) & 0b00000001;
	     unsigned char ce = ((unsigned char)(current >>5)) & 0b00000001;;
	     unsigned char c1 = (unsigned char)(current >> 6); 
	     unsigned char c2 = (unsigned char)(current >> 14); 
             state =          ((      c1<=max1 && ne == 1)?1:0);
             state = state | (((size1-c1<=max1 && ne == 1)?1:0)<<1);
             state = state | (((      c2<=max2 && ce == 1)?1:0)<<2);
             state = state | (((size2-c2<=max2 && ce == 1)?1:0)<<3);

             unsigned int  vi = ((unsigned int)(current >>22)) & 0x1FFFFF;//clear left bits
             unsigned int  vj = ((unsigned int)(current >>43)) & 0x1FFFFF;

	     if (job == 1 ||VERBOSE)
	     {
               /* the first */
	       if (i1==0)
	       {
                  if (state==1||state==5||state==9)  // 1+{0, 4, 8}
	          {
	             fprintf(out_file1, "%d\t%d", vj, vi);
                     if (VERBOSE)
                        fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                     fprintf(out_file1, "\n");
	          }
	          else if (state==2||state==6||state==10) // 2+{0, 4, 8}
	          {
	             c1 = size1 - c1;
	             fprintf(out_file1, "%d\t%d", vj, vi);
                     if (VERBOSE)
                        fprintf(out_file1, "\t%d\t%.6g", c1, bino_p(size1, c1));
                     fprintf(out_file1, "\n");
	          }
	       }//end if i1==0

               /* the rest */
               if (state==4||state==5||state==6)  // 4 +{0,1,2}
	       {
	           fprintf(out_file2, "%d\t%d", vj, vi);
                   if (VERBOSE)
                      fprintf(out_file2, "\t%d\t%.6g", c2, bino_p(size2, c2));
                   fprintf(out_file2, "\n");
	        }
	        else if (state==8||state==9||state==10) // 8+{0, 1, 2} 
	        {
	           c2 = size2 - c2;
	           fprintf(out_file2, "%d\t%d", vj, vi);
                   if (VERBOSE)
                      fprintf(out_file2, "\t%d\t%.6g", c2, bino_p(size2, c2));
                   fprintf(out_file2, "\n");
	        }
	     }//end if job==1 || VERBOSE

	     /* keep only the consistent pairs */
	     switch(state)
	     {
	        case  5: // 1+4
                  if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vi, vj);
		  count_same += 1;
	          break;
	        case 10: // 2+8
                  if (job==1||VERBOSE) fprintf(out_file3, "%d\t%d\n", vj, vi);
		  count_same += 1;
	          break;
	        case  6: // 2+4
                  if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vj, vi);
		  count_reverse += 1;
	          break;
	        case  9: // 1+8
                  if (job==1||VERBOSE) fprintf(out_file4, "%d\t%d\n", vi, vj);
		  count_reverse += 1;
	          break;
		default:
		  pairs[j][k] = 0;//clear the other data
	     }//end switch

           }//end for
	fprintf(stdout, "#  Number of concordant pairs between data set %d and %d: %d\n", i2, i1, count_same);
	fprintf(stdout, "#  Number of reversed pairs between data set %d and %d: %d\n", i2, i1, count_reverse);
        /* clean up */
	if (job==1||VERBOSE) 
	{
           if (i1==0 && out_file1 != stderr) fclose(out_file1);
           if (out_file2 != stderr) fclose(out_file2);
           if (out_file3 != stderr) fclose(out_file3);
           if (out_file4 != stderr) fclose(out_file4);
	}

        for (j=0; j<n_threads; j++) 
        {
           //free(pairs[j]);
           free(exceptions1[j]);
           free(exceptions2[j]);
        }
        //free(pairs);
        //free(count_pairs);
        free(exceptions1);
        free(exceptions2);
        free(data2);

	/* filter dysregulated genes */
        //state for each gene, not pairs
        struct gene_state3 **genes;//an array of pointers to struct 
        //allocate an array of n pointers to struct gene_states
        genes = (struct gene_state3 **) malloc(sizeof(struct gene_state3 *)*n_genes);
        //initialization
        for (j=0; j<n_genes; j++) genes[j] = (struct gene_state3 *) NULL;
	//call
	fprintf(stdout, "#  Filter genes based Bayesian algo...\n");
	filter_gene_bayesian(n_genes, n_threads, pairs, count_pairs, genes);

	/* write statistics on genes into files*/
	/* open file to write stable pairs */
        out_name3 = malloc(sizeof(char *)*10);
        out_name3 = make_file_name2(BAYESIAN, i2, "_", i1, EXT, out_name3);
        out_file3 = fopen(out_name3, "w");
        if (out_file3 == NULL) 
        {
           out_name3 = "stderr";
           out_file3 = stderr;
        }
        fprintf(stdout, "#  Writing statistics on  genes into file: %s\n", out_name3);
        fprintf(stdout, "#  FORMAT: # greater and # less than in concordant pairs, # greater to less, # less to greater in reversed pairs\n");
        free(out_name3);

        for(j=0;j<n_genes;j++)
           if(genes[j]!=NULL) 
	   {
              fprintf(out_file3, "%d\t%d\t%d\t%d\t%d\n", j, genes[j]->cg, genes[j]->cl, genes[j]->g2l, genes[j]->l2g);    
              free(genes[j]);
           }

        if (out_file3 != stderr) fclose(out_file3);
	free(genes);

	/* clean up*/
        for (j=0; j<n_threads; j++) 
           free(pairs[j]);
        free(pairs);
        free(count_pairs);

     }//end inner for
     /*clean up*/
     free(data1);

  }//end outer for
  return EXIT_SUCCESS;
} // end of select_genes_bayesian

/* generate simulation data set */
int  gen_simulation_data(int  times,
		        char *fname,
			int   size,
			int   n_genes,
			int   n_changes,
			char *cname,
			char *simu_names[times+1],
			int   max_equals,
			char  job,
			float alpha, 
			int   max_cycles, 
			int   ori_cycles, 
			int   max_threshold,
			char  algorithm
	              )
{

  /* load data1 into memory */
  DATATYPE_VALUE *template = malloc(sizeof(DATATYPE_VALUE) * n_genes * size);
  if (template == NULL) 
  {
     fprintf(stderr, "ERROR: memory is not allocated for template file  %s.\n", fname);
     exit(EXIT_FAILURE);
  }
  read_data(fname, n_genes, size, template);

  /* read fold change into memory */
  struct CHANGE *changes;
  changes = malloc(sizeof(struct CHANGE) * n_changes);
  read_changes(cname, n_changes, changes);

  /* generate simulation data */
  //char *simu_names[times+1];

  int len = snprintf(NULL, 0, "%s", fname);
  simu_names[0] = malloc(sizeof(char*)*len);
  snprintf(simu_names[0], len+1, "%s", fname);

  int i;
  for (i=0; i<times; i++)
  {
    DATATYPE_VALUE *cases;
    cases = malloc(sizeof(DATATYPE_VALUE) * n_genes * size); 
    int j, k;
    for(j=0;j<n_genes;j++) 
       for(k=0;k<size;k++)
          cases[j*size+k]=template[j*size+k];

    if(cases==NULL) {
      printf("ERROR: array cases is not allocated. NOT ENOUGH MEMORY!\n");
      exit(EXIT_FAILURE);
    }


    /* Write sample into file */
    /* open file to write consistent pairs */
    char *out_name = malloc(sizeof(char *)*10);
          out_name = make_file_name(SAMPLED_GENES, i, EXT, out_name);
    FILE *out_file = fopen(out_name, "w");
    if (out_file == NULL) 
    {
       out_name = "stderr";
       out_file = stderr;
    }
    fprintf(stdout, "#  Writing sampled genes index to file: %s\n", out_name);
    fprintf(stdout, "#  FORMAT: group #, gene index # \n");

    gen_random_sample(n_genes, size, template, cases, n_changes, changes, out_file);
    
    if (out_file != stderr) fclose(out_file);
     

    char *current_name;
    FILE *current_file;
    current_name = malloc(sizeof(char *)*10);
    current_name = make_file_name(SIMULATION, i, EXT, current_name);
    current_file = fopen(current_name, "w");
    if (current_file == NULL) 
    {
       current_name = "stderr";
       current_file = stderr;
    }
    fprintf(stdout, "#  Writing simulation data into file: %s\n", current_name);

    len = snprintf(NULL, 0, "%s", current_name);
    simu_names[i+1] = malloc(sizeof(char*)*len);
    snprintf(simu_names[i+1], len+1, "%s", current_name);
    free(current_name);

    for(j=0;j<n_genes;j++) 
    {
       for(k=0;k<size;k++)
          fprintf(current_file, "%.6g\t", cases[j*size + k]);
       fprintf(current_file, "\n");
    }//end for

    if (current_file != stderr) 
       fclose(current_file);
    free(cases);
  }//end for

  /* filtering genes */
  /* clean up */
  //for (i=0; i<times+1; i++)
  //   free(simu_names[i]);
  free(changes);
  free(template);
  return EXIT_SUCCESS;

}// end of gen_simulation_data

/*Job Type 5 */
int calculate_similarity(DATATYPE_GENESIZE   n,
		         DATATYPE_SAMPLESIZE m,
			 char *file_name, 
			 char *out_name
		        )
{

  DATATYPE_VALUE *data = malloc( sizeof(DATATYPE_VALUE) * n * m);
  if (data == NULL) 
  {
     fprintf(stderr, "ERROR: memory is not allocated for data file  %s.\n", file_name);
     exit(EXIT_FAILURE);
  }
  read_data(file_name, n, m, data);

  FILE *out_file = fopen(out_name, "w");
  if (out_file == NULL) 
  {
     out_name = "stderr";
     out_file = stderr;
  }
  fprintf(stdout, "#  Writing gene-pair similarity into file: %s\n", out_name);
  fprintf(stdout, "#  FORMAT: col1 #, col2 #, count of the same pair order, percentage \n");

  DATATYPE_SAMPLESIZE i, j;
  DATATYPE_VALUE *col1, *col2;
  col1 = malloc( sizeof(DATATYPE_VALUE) * n);
  col2 = malloc( sizeof(DATATYPE_VALUE) * n);
  double total = (double) n*(n-1)/2;
  for (i=0; i<m-1; i++)
  {
    extract_column(n, m, data, i, col1);
    for (j=i+1; j<m;j++)
      {
       extract_column(n, m, data, j, col2);
       DATATYPE_GENESIZE count = similarity_between_two_cols(n, col1, col2);
       fprintf(out_file, "%d\t%d\t%d\t%.4g\n", i, j, count, (double) count / total);
      }
  }
  fclose(out_file);
  free(col2);
  free(col1);
  free(data);
  return EXIT_SUCCESS;
}//end function calculate_similarity


/*Job Type 6 */
int calculate_simi_three(DATATYPE_GENESIZE   n,
		         DATATYPE_SAMPLESIZE m,
			 char *file_name, 
			 char *out_name
		        )
{

  DATATYPE_VALUE *data = malloc( sizeof(DATATYPE_VALUE) * n * m);
  if (data == NULL) 
  {
     fprintf(stderr, "ERROR: memory is not allocated for data file  %s.\n", file_name);
     exit(EXIT_FAILURE);
  }
  read_data(file_name, n, m, data);

  FILE *out_file = fopen(out_name, "w");
  if (out_file == NULL) 
  {
     out_name = "stderr";
     out_file = stderr;
  }
  fprintf(stdout, "#  Writing gene-pair similarity into file: %s\n", out_name);
  fprintf(stdout, "#  FORMAT: col1 #, col2 #, col3 #, count of the same pair order, percentage \n");

  DATATYPE_SAMPLESIZE i, j, k;
  DATATYPE_VALUE *col1, *col2, *col3;
  col1 = malloc( sizeof(DATATYPE_VALUE) * n);
  col2 = malloc( sizeof(DATATYPE_VALUE) * n);
  col3 = malloc( sizeof(DATATYPE_VALUE) * n);
  double total = (double) n*(n-1)/2;
  for (i=0; i<m-1; i++)
  {
    extract_column(n, m, data, i, col1);
    for (j=i+1; j<m; j++)
      {
       extract_column(n, m, data, j, col2);
       for (k=j+1; k<m; k++)
	  //if ( !(i==j && j==k) )
	  {
            extract_column(n, m, data, k, col3);
            DATATYPE_GENESIZE count = similarity_between_three_cols(n, col1, col2, col3);
            fprintf(out_file, "%d\t%d\t%d\t%d\t%.4g\n", i, j, k, count, (double) count / total);
	  }
      }
  }
  fclose(out_file);
  free(col3);
  free(col2);
  free(col1);
  free(data);
  return EXIT_SUCCESS;
}//end function calculate_similarity

/* Job type: 5 or 6 */
int calc_sample_similarity(int               n_files, 
		           char      *files[n_files],
			   int sample_sizes[n_files],
			   int              n_genes,
			   char                 job
		          ) 
{
  int i;

  for (i=0; i<n_files; i++)
  {
     fprintf(stdout, "#Processing data file %d...\n", i);
     int  size_current = sample_sizes[i];
     char *file_current = files[i];
     char *out_name;
     out_name = malloc(sizeof(char *)*10);
     out_name = make_file_name(SIMILARITY, i, ".dat", out_name);
     if (job == 5)
        calculate_similarity(n_genes, size_current, file_current, out_name);
     else
        calculate_simi_three(n_genes, size_current, file_current, out_name);
     free(out_name);
  } //end if 
  return EXIT_SUCCESS;
}//end function calc_sample_similarity 

int main (int argc, char **argv)
{
  /* default values */
        VERBOSE           = false;
        EPSILON           = FLT_EPSILON;
  int   max_equals        = 256;
  float fdr_level         = 0.05;
  char  sample_type       = 0;
  char  job_type          = 2;
  char  pair_mode         = 0;
  int   times             = 1;
  char  algorithm         = 0; 
  int   max_cycles        = 128;
  int   ori_cycles        = 2;
  int   conv_threshold    = 50;
  char *cname             = "";
  int   n_changes         = 0;

  int   int_current;
  float float_current;

  /* check options */
  int opt;
  int option_index = 0;
  while((opt = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1)
    {
      switch (opt)
        {
        case 'h':
	  display_usage(stdout, argv[0]);
	  exit(EXIT_SUCCESS);
          break;

        case 'V':
	  display_version(stdout, argv[0]);
	  exit(EXIT_SUCCESS);
          break;

        case 'v':
	  VERBOSE = true;
          break;

        case 'e':
	  float_current = atof(optarg);
	  if (float_current < 0.)
	  {
	    fprintf(stderr, "ERROR: --epsilon should be given an non-negative argument.\n");
	    exit(EXIT_FAILURE);
	  }
	  else EPSILON = (double) float_current;
          break;

        case 'q':
	  int_current = atoi(optarg);
	  if (int_current < 0)
	  {
	    fprintf(stderr, "ERROR: --equals should be given an non-negative integer.\n");
	    exit(EXIT_FAILURE);
	  }
	  else max_equals = int_current;
          break;

        case 'f':
	  float_current = atof(optarg);
	  if (float_current < 0.|| float_current > 1.)
	  {
	    fprintf(stderr, "ERROR: --fdr should be given a float between 0 and 1.\n");
	    exit(EXIT_FAILURE);
	  }
	  else fdr_level = float_current;
          break;

        case 's':
	  int_current = atoi(optarg);
	  if (int_current < 0 || int_current > 1)
	  {
	    fprintf(stderr, "ERROR: --sample should be given an integer between 0 and 1.\n");
	    exit(EXIT_FAILURE);
	  }
	  else sample_type = int_current;
          break;

        case 'j':
	  int_current = atoi(optarg);
	  if (int_current < 0 || int_current > 7)
	  {
	    fprintf(stderr, "ERROR: --job should be given an integer between 0 and 7.\n");
	    exit(EXIT_FAILURE);
	  }
	  else job_type = int_current;
          break;

        case 'p':
	  int_current = atoi(optarg);
	  if (int_current < 0 || int_current > 2)
	  {
	    fprintf(stderr, "ERROR: --pair should be given an integer between 0 and 2.\n");
	    fprintf(stderr, "ERROR: default mode, %d, will be used.\n", pair_mode);
	    //exit(EXIT_FAILURE);
	  }
	  else pair_mode = int_current;
          break;

        case 'c':
          if ( access(optarg, R_OK) == -1)
          {
              fprintf(stderr, "ERROR: change configuration file, --changes %s does not exist or not readable.\n", optarg);
              exit(EXIT_FAILURE);
          }
	  else cname = optarg;
          break;

        case 'l':
	  int_current = atoi(optarg);
	  if (int_current < 1)
	  {
	    fprintf(stderr, "ERROR: --lines should be given a positive  integer.\n");
	    exit(EXIT_FAILURE);
	  }
	  else n_changes  = int_current;
          break;

        case 't':
	  int_current = atoi(optarg);
	  if (int_current < 0)
	  {
	    fprintf(stderr, "ERROR: --times should be given an non-negative integer.\n");
	    exit(EXIT_FAILURE);
	  }
	  else times = int_current;
          break;

        case 'a':
          algorithm = atoi(optarg);
	  if (algorithm < 0 || algorithm > 2) algorithm = 0; 
          break;

        case 'm':
          max_cycles = atoi(optarg);
	  if (max_cycles < 0) max_cycles = 0; 
          break;

        case 'o':
          ori_cycles = atoi(optarg);
	  if (ori_cycles < 0) ori_cycles = 0; 
          break;

        case 'n':
          conv_threshold = atoi(optarg);
	  if (conv_threshold < 0) conv_threshold = 0; 
          break;

        case 'M':
	  float_current = atof(optarg);
	  if (float_current < 0.|| float_current > 1.)
	  {
	    fprintf(stderr, "ERROR: --memory should be given a float between 0 and 1.\n");
	    fprintf(stderr, "ERROR: automatically determined value will be used .\n");
	    //exit(EXIT_FAILURE);
	  }
	  else MEMORY_USAGE = float_current;
          break;

        case '?':
          /* getopt_long already printed an error message. */
	  display_usage(stderr, argv[0]);
	  exit(EXIT_FAILURE);

        default:
          abort ();
        } //end switch
    } //end while

  int n_files, n_genes, sample_sizes[MAX_FILES];
  float fdr_levels[MAX_FILES];
  char *files[MAX_FILES];

  /* check number of files */
  if (optind < argc)
    {
      n_files = atoi(argv[optind++]);
      if (n_files<=0)
      {
	  fprintf(stderr, "ERROR: n_files should be an non-negative integer.\n");
	  exit(EXIT_FAILURE);
      }
    }
  else
    {
      fprintf(stderr, "ERROR: too few arguments.\n\n"); 	  
      display_usage(stderr, argv[0]); 	  
      exit(EXIT_FAILURE);
    }

  /* check number of genes */
  if (optind < argc)
    {
      n_genes = atoi(argv[optind++]);
      if (n_genes<=0)
      {
	  fprintf(stderr, "ERROR: n_genes should be an non-negative integer.\n\n");
	  exit(EXIT_FAILURE);
      }
    }
  else
    {
      fprintf(stderr, "ERROR: too few arguments.\n\n"); 	  
      display_usage(stderr, argv[0]); 	  
      exit(EXIT_FAILURE);
    }

  int i;
  char *file_name;

  /* check file names  */
  for(i=0;i<n_files;i++)
    {
      if (optind < argc)
        {
          file_name = argv[optind++];
          if ( access(file_name, R_OK) == -1)
          {
              fprintf(stderr, "ERROR: data file %s does not exist or not readable.\n", file_name);
              exit(EXIT_FAILURE);
          }
	  else files[i] = file_name;
        }
      else
        {
          fprintf(stderr, "ERROR: number of file names does not match with n_files.\n\n"); 	  
          display_usage(stderr, argv[0]); 	  
          exit(EXIT_FAILURE);
        }//end if...else...

    }//end for

  /* check sample size */
  int size_current;
  for(i=0;i<n_files;i++)
    {
      if (optind < argc)
        {
          size_current = atoi(argv[optind++]);
          if (size_current<=0 || size_current > MAX_SAMPLESIZE)
          {
              fprintf(stderr, "ERROR: sample size %d is too small(<=0) or too big(>255).\n", size_current);
              exit(EXIT_FAILURE);
          }
	  else sample_sizes[i]=size_current;
        }
      else
        {
          fprintf(stderr, "ERROR: number of sample sizes does not match with n_files.\n\n"); 	  
          display_usage(stderr, argv[0]); 	  
          exit(EXIT_FAILURE);
        }//end if...else...
    }//end for

  /* check fdr_level or exception number */
  float fdr_current;
  i = 0;
  while (optind < argc && i<n_files)
    {
      fdr_current = atof(argv[optind++]);
      if (fdr_current<0.)
         {
           fprintf(stderr, "ERROR: FDR level or exception threshold  %f should be non-negative.\n", fdr_current);
           exit(EXIT_FAILURE);
         }
      else fdr_levels[i++]=fdr_current;
    } //end while
  while (i<n_files)
     fdr_levels[i++]=DEFAULT_FDR;

  char *simu_names[times+1];
  /* running the jobs*/
  switch(job_type)
  {
    case 0:
      fprintf(stdout, "#Job type 0\n");
      select_stable_pairs(n_files, files, sample_sizes, fdr_levels, n_genes, max_equals);
      break;

    case 5:
      fprintf(stdout, "#Job type 5\n");
    case 6:
      fprintf(stdout, "#Job type 6\n");
      calc_sample_similarity(n_files, files, sample_sizes, n_genes, job_type);
      break;

    case 7:
      fprintf(stdout, "#Job type 6\n");
      select_genes__bayesian(n_files, files, sample_sizes, fdr_levels, n_genes, max_equals, \
		             job_type, pair_mode);
      break;

    case 1:
      fprintf(stdout, "#Job type 1\n");
    case 2:
      fprintf(stdout, "#Job type 2\n");
      if (sample_type==0)
      {
        fprintf(stdout, "#Sample Type: Cohort Sample\n");
        select_consistent_pairs(n_files, files, sample_sizes, fdr_levels, n_genes, max_equals, \
		              job_type, pair_mode, fdr_level, max_cycles, ori_cycles, conv_threshold, algorithm);
      }
      else
      {
        fprintf(stdout, "#Sample Type: Individual Sample\n");
	if (VERBOSE)
           select_consistent_pairs_ind(n_files, files, sample_sizes, fdr_levels, n_genes, max_equals, \
		              job_type, pair_mode, fdr_level, max_cycles, ori_cycles, conv_threshold, algorithm);
	else
           select_consistent_pairs_ind2(n_files, files, sample_sizes, fdr_levels, n_genes, max_equals, \
		              job_type, pair_mode, fdr_level, max_cycles, ori_cycles, conv_threshold, algorithm);
      }
      break;

    case 3:
      fprintf(stdout, "#Job type 3\n");
      gen_simulation_data(times, files[0], sample_sizes[0], n_genes, n_changes, cname, simu_names,
			max_equals, job_type, fdr_level, max_cycles, ori_cycles,conv_threshold, algorithm);
      break;

    case 4:
      fprintf(stdout, "#Job type 4\n");
      gen_simulation_data(times, files[0], sample_sizes[0], n_genes, n_changes, cname, simu_names,
			max_equals, job_type, fdr_level, max_cycles, ori_cycles,conv_threshold, algorithm);

      float local_fdr_levels[times+1];
      int   local_sizes[times+1];
      for(i=0;i<times+1;i++)
      {
	 local_fdr_levels[i] = fdr_levels[0];
	 local_sizes[i] = sample_sizes[0];
      }
      select_consistent_pairs(times+1, simu_names, local_sizes, local_fdr_levels, n_genes, max_equals, \
		            2, 0, fdr_level, max_cycles, ori_cycles, conv_threshold, algorithm);
      for(i=0;i<times+1;i++)
	 free(simu_names[i]);
      break;

  } //end switch

  return EXIT_SUCCESS;
}

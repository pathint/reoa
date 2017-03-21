#include <getopt.h>
#include <unistd.h>

#define STABLE      "stable_pairs_"
#define CONCORDANT  "concordant_pairs_"
#define REVERSED    "reversed_pairs_"
#define UP	    "up_regulated_"
#define DOWN	    "down_regulated_"
#define UP_ORIG	    "up_orig_algo_"
#define DOWN_ORIG   "down_orig_algo_"
#define GENE_STATE  "gene_state_"
#define SIMULATION  "simulation_"
#define SIMILARITY  "similarity_"
#define BAYESIAN    "bayesian_"
#define SAMPLED_GENES "sampled_genes_"
#define EXT	    ".dat"

static struct option long_options[] =
  {
    {"verbose", no_argument,       NULL, 'v'},
    {"Version", no_argument,       NULL, 'V'},
    {"help",  	no_argument,       NULL, 'h'},
    {"epsilon", required_argument, NULL, 'e'},
    {"equals",  required_argument, NULL, 'q'},
    {"fdr",     required_argument, NULL, 'f'},
    {"sample",  required_argument, NULL, 's'},
    {"job",     required_argument, NULL, 'j'},
    {"pair",    required_argument, NULL, 'p'},
    {"platform",required_argument, NULL, 'P'},
    {"changes", required_argument, NULL, 'c'},
    {"lines",   required_argument, NULL, 'l'},
    {"plines",  required_argument, NULL, 'L'},
    {"time",    required_argument, NULL, 't'},
    {"algorithm",required_argument,NULL, 'a'},
    {"cycles",  required_argument, NULL, 'm'},
    {"cycles_orig",required_argument, NULL, 'o'},
    {"convergence",required_argument, NULL, 'n'},
    {"memory",  required_argument, NULL, 'M'},
    {NULL, 	no_argument, 	   NULL, 0}
  };

static const char *opt_string = "e:q:f:s:j:p:P:c:l:L:t:a:m:o:n:M:vVh?";

float MEMORY_USAGE = -1.;

// could be better
static void display_usage(FILE *f, const char *program_name);
static void display_version(FILE *f, const char *program_name);

char *my_itoa(int value, char *str);
char *make_file_name(char *stem, int index, char *ext, char *fname);
char *make_file_name2(char *stem, int i1, char *in, int i2, char *ext, char *fn);
char *make_file_name3(char *stem, int i1, char *in1, int i2, char *in2, int i3, char *ext, char *fn);

int read_changes(const char *file_name, int n_changes, struct CHANGE changes[n_changes]);
int read_data(char *file, int n_genes, int sample_size, DATATYPE_VALUE *data);
int extract_column(unsigned int n, unsigned int m, DATATYPE_VALUE data[n*m], 
		   unsigned int c, DATATYPE_VALUE column[n]);

int select_stable_pairs(int               n_files, 
		        char      *files[n_files],
			int sample_sizes[n_files],
			float        fdr[n_files],
			int              n_genes ,
			int            max_equals
		       ); 

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
		       );

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
		       );

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
	              );

int calculate_similarity(DATATYPE_GENESIZE   n,
		         DATATYPE_SAMPLESIZE m,
			 char *file_name, 
			 char *out_name
		        );

/*Job Type 6 */
int calculate_simi_three(DATATYPE_GENESIZE   n,
		         DATATYPE_SAMPLESIZE m,
			 char *file_name, 
			 char *out_name
		        );

int calc_sample_similarity(int               n_files, 
		           char      *files[n_files],
			   int sample_sizes[n_files],
			   int              n_genes,
			   char                 job
		          ); 

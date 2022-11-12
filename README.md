## Update on Nov 12, 2022
Users are advised to use the updated RankCompV3 algorithm (https://github.com/pathint/RankCompV3.jl) which is implemented with Julia.

# Relative Expression Ordering Analysis (REOA) package

This program implements the original RankComp and the RankCompV2 algorithms to detect the dysregulated genes. It is applicable to gene expression, methylation and other gene microarray data sets. Both large-scale cohort and individual samples or small-scale cell-line samples with two or three technical replicates can be analyzed. 



Besides the main program, `reoa`, this package also includes several derivative programs which apply the method to specific scenarios. 

1. *CellComp*.
The `cellcomp` program is applicable to small-scale cell-line data which include only a few (e.g. two or three) technical replicates. See [below](https://github.com/pathint/reoa#small-scale-cell-line-data-with-two-or-three-technical-replicates) for usage.

2. *OneComp*.
The `onecomp` program is applicable to the datasets with only one case sample and possibly one control sample. See [below](https://github.com/pathint/reoa/blob/master/README.md#datasets-with-only-one-control-sample-and-one-case-sample) for usage.


If you publish results from using this program, please cite the following publications.

Xiangyu Li, Hao Cai, Xianlong Wang, Lu Ao, You Guo, Jun He, Yunyan Gu, Lishuang Qi, Qingzhou Guan, Xu Lin, Zheng Guo.
A rank-based algorithm of differential expression analysis for small cell line data with statistical control.
Briefings in Bioinformatics, Volume 20, Issue 2, March 2019, Pages 482–491, https://doi.org/10.1093/bib/bbx135.

Hongwei Wang, Qiang Sun, Wenyuan Zhao, Lishuang Qi, Yunyan Gu, Pengfei Li, Mengmeng Zhang, Yang Li, Shu-Lin Liu, and Zheng Guo. (**2015**). Individual-level analysis of differential expression of genes and pathways for personalized medicine. *Bioinformatics*. 31(1):62-8. [DOI:10.1093/bioinformatics/btu522](http://dx.doi.org/10.1093/bioinformatics/btu522).


## Algorithm

### Significantly Stable Gene Pairs
A pair of genes, {*a*, *b*}, is considered as statistically stable if they hold the same order relationship (*a* < *b* or *a* >*b*) in most of the samples. Binomial distribution is used to calculate the *P* value under the null hypothesis (a and b does not have a stable order relation) for the large-scale samples. The [Benjamini–Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate) is used to control the false discovery rate (FDR) at level alpha, which is 0.05 by default. For small-scale samples, such as technical replicates in the cell-line studies, the exception number is used to screen the stable pairs. For example, if the exception number is set as 0, the pair holds the same order relation in all the samples. 

*Note on Notation Abuse*: *a*, *b* are used to represent both genes and the expression values (or methylation level etc) of the genes. Readers should be able to discern what the notations represent with no trouble. 

### Concordant and Reversal Gene Pairs
If a pair of genes, {*a*, *b*}, is significantly stable in both the normal group and the disease group and the orders are the same (either *a* < *b* or *a* > *b*), the pair is a concordant pair. Likewise, if a pair is stable in both groups, but the orders are different (*a* < *b* in one group, but *a* > *b* in the other group), the pair is called reversal pair.  Concordant and reversal gene pairs are the basic data to detect dysregulated genes.   


### Dysregulated Genes, the RankCompV2 Algorithm
Whether a gene is dysregulated or not in a disease cohort compared with normal cohort is identified through the original RankComp algorithm or through the improved RankCompV2 algorithm. Here is a brief description on the algorithms. Interested users are advised to read the related papers.

All the tested genes are going to be classified into three possible regulation directions, up-regulated, down-regulated and non-dysregulated. 

For one gene, *a*, we count the numbers of genes whose expression (or methylation etc) levels are higher and lower than this gene in the normal cohort and in the disease cohort. The counting is carried out for the concordant and reversal gene pairs only. Thus, we obtain the following contingency table.

Group | Number of genes whose expression levels are higher than *a* | Number of genes whose expression levels are less than *a* 
---- | ------------ | -------------
Normal | *n*<sub>g</sub> | *n*<sub>l</sub>
Disease |*d*<sub>g</sub> | *d*<sub>l</sub>

[Fisher exact test](https://en.wikipedia.org/wiki/Fisher's_exact_test) is then used to calculate the *P* value under the null hypothesis which states that the number of genes whose expression levels are higher or lower than the tested gene has no association with disease state (absence or presence). The [Benjamini–Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate) is used to control the false discovery rate (FDR) at level alpha, which is 0.05 by default. If the null hypothesis is rejected, the gene *a* is identified as a potential dysregulated gene (DEG). The dysregulation direction is judged by comparing *n*<sub>g</sub> / *n*<sub>l</sub> with *d*<sub>g</sub> / *d*<sub>l</sub> (equivalently, *n*<sub>g</sub> / (*n*<sub>g</sub> + *n*<sub>l</sub>) vs. *d*<sub>g</sub> / ( *d*<sub>g</sub> + *d*<sub>l</sub>)). 

The above describes the first step in both the algorithms. In this step, all the genes are considered as non-DEGs when the contingency table is constructed. After the gene states are assigned, we repeat the above steps to further filter potential DEGs. However, the counting to build the contingency tables does not include all stable pairs. Whether a stable pair, either concordant or reversal, is included depends on the dysregulation direction of the partner gene in the pair.  For a stable pair {*a*, *b*}, only if gene *b* is assigned as non-dysregulated gene, the pair is counted for the construction of the contingency table of gene *a*. (Initially, all the genes are assumed to be non-dysregulated, that is to say all the stable pairs will be used to construct the table in the initial cycle.)  After the contingency tables are renewed, Fisher exact test with FDR control are called for the null hypothesis test again. This iteration continues until the number of DEGs does not change significantly anymore.              

### Dysregulated Genes, the original RankComp Algorithm
In the original RankComp algorithm, the later iterations use a different strategy to renew the contingency tables. We redo the counting based on the reversal pairs only. For gene *a*, we count the number of genes whose expression level is greater than that of *a* in the normal group but becomes less than that of *a* in the disease group and the number of genes whose expression level is less than that of *a* in the normal group but becomes greater than that of *a* in the disease group. The numbers are labeled as *g2l* and *l2g*, respectively. For a reversal pair, *a* < *b* in the normal group and *a* > *b* in the disease group, if *a* is up-regulated, *a* will not be counted as *g2l* for *b* and if *b* is down-regulated, *b* is not counted as *l2g* for *a*. The recounted values *g2l* and *l2g*, and the values of *n*<sub>g</sub> and *n*<sub>l</sub> from the first step, are used to construct the above contingency tables. Fisher exact test is again used to identify the dysregulated genes.

## Usage
The online help can be invoked by running the progrom with `-h` or `--help` option. 
```
reoa -h
```
or
```
reoa --help
```
The input data sets should be given as text-based data matrix data files. One file contains the microarray value matrix of one group of samples. The number of rows corresponds to the total number of gene probes and the number of columns corresponds to the sample size.  The values should be tab or space delimited. If multiple data sets are given as input, the number of genes should be the same. 

The output gene pairs or dysregulated genes are given by the indices (starting from 0) of the genes. If the verbose option (`-v` or `--verbose`) is turn on, extensive amount of output will be printed out.  Please read the standard (screen) output for the names of the output files and the content format. 

The most important option is `-j` (`--job`) which sets up the job type. Job types 0-2 identify stable gene pairs and dysregulated genes. Job types 3-4 generate simulation data sets and identify dysregulated genes in the simulated data sets. Job types 5-6 calculate the concordance scores between two samples in one data set. 

For most of job types, multiple data sets can be given. Three comparison modes are avaiable for selecting stable pairs and identifying dysregulated genes and related jobs. The choice can be made through the `-p` (`--pair`) option.  

The choice of algorithm is specified by the `-a` (`--algorithm`) option. 

For the small-scale samples, such as cell line technical replicates, exception numbers should be used to select stable pairs instead of FDR levels. And the exception number is usually set as 0. Furthermore, the option of `-q` (`--equals`) should be set as 0 as well. 

The default sample type is 'cohort'. For individual-type samples, the option `-s` (`--sample`) should be set as 1. Each sample (column) in a data set is analyzed separately.    

See the folder `test` for test cases. 

## Application Scenarios

### Small-Scale Cell Line Data with Two or Three Technical Replicates
For the small-scale cell line expriments with only a few (e.g. two or three) technical replicates, please use the program `cellcomp` to detect differentially expressed genes (DEGs). 

Example. (The sample expression files are available under the `test` folder.) 
```
cellcomp cellcomp_control.dat cellcomp_case.dat
```
Detailed usage information for `cellcomp` is as follows. 

#### NAME
     cellcomp, a bash script driver for application of reoa
     to small-scale data set with only a few technical replicates
#### SYNOPSIS
     cellcomp [-h|-V]
     cellcomp [OPTIONS] control_file case_file
#### DESCRIPTION
     CellComp (cellcomp) is a program to apply the RanComp and 
     RankCompV2 algorithms to detect differentially expressed 
     genes (DEGs) in small-scale cell line experiments. In such 
     experiments, the gene expression profiles contain only two 
     or three technical replicates and most conventional methods
     lack statistical power when applied to such small-scale data
     sets. The RankComp and RanCompV2 are based on the analysis of
     within-sample relative expression orderings (REOs) of gene pairs.
     See the references more details.
#### OPTIONS
     -h, --help, --usage
       Show this message
     -V, --version
       Show program version
     -a, --algorithm = ALGORITHM
       Choice of algorithm. Default: 0
       0 - RankCompV2 only
       1 - Original RankComp only
       2 - Both RankComp and RankCompV2
     -f, --fdr = FDR
       False discovery rate (FDR) level. Default: 0.05
     -c, --cycles = CYCLES
       Maximum iteration cycles. Default: 128
     -t, --threshold = THRESHOLD
       Convergence criterion. Default: 50 
       (Maximum fluctuation in number of DEGs)
     -v, --verbose
       Output extra information

### Datasets with Only One Control Sample and One Case Sample
Under this application scenario, background gene pairs, which show stable REOs in many normal samples, must be obtained before-hand. This can be achieved by applying `reoa` to the merged normal samples of a particular tissue or a cell line from the same or different data sources. The program `onecomp` first customizes the background gene pairs by filtering out those pairs with reversed REOs in the control sample(s) of the current dataset, then detects the DEGs in a case sample with the filtered background gene pairs.   

#### NAME
     onecomp, a bash script driver for application of reoa
     to detect dysregulated genes in one treated sample
     given a list of stable gene pairs and one paired control
     sample.
     
#### SYNOPSIS
     onecomp [-h|-V]
     onecomp [OPTIONS] pair_file control_file treated_file
     
#### DESCRIPTION
     OneComp (onecomp) is a program to apply the original RankComp 
     and RankCompV2 algorithms to detect differentially expressed 
     genes (DEGs) in one treated sample given the paired control 
     sample and a list of predetermined gene pairs with stable REOs.
     The control sample is used to customize the gene pair list
     and 'paired' here does not necessarily mean the expermental
     design but a compatiable corresponding relation between the two.
     The 'control_file' and 'treated_file contains the expression
     profiles of a group of control samples and a group of treated
     samples. They should have the same number of rows (gene probes) 
     and the same number of columns (number of samples). The i-th
     column of 'control_file' and the i-th column of 'treated_file' 
     constitue one paired sample and they are processed by the algo.
     independent of the other columns.
     The pair_file contains the list of predetermined gene pairs which
     are given by the gene probe indices (0-based) which have the same
     order as the expression files (control_file or treated_file).
     The RankComp and RanCompV2 are based on the analysis of
     within-sample relative expression orderings (REOs) of gene pairs.
     See the references more details.
     
#### OPTIONS
     -h, --help, --usage
       Show this message
     -V, --version
       Show program version
     -a, --algorithm = VALUE
       Choice of algorithm. Default: 0
       0 - RankCompV2 only
       1 - Original RankComp only
       2 - Both RankComp and RankCompV2
     -m, --mode = VALUE
       Choice of filter mode. Default: 0
       0 - Use one control sample to filter the pairs to detect DEGs
           in the corresponding treated sample
       1 - Use all the control samples to filter the pairs to detect
           DEGs in each treated samples
     -f, --fdr = FDR
       False discovery rate (FDR) level. Default: 0.05
     -c, --cycles = CYCLES
       Maximum iteration cycles. Default: 128
     -t, --threshold = THRESHOLD
       Convergence criterion. Default: 50 
       (Maximum fluctuation in number of DEGs)
     -v, --verbose
       Output extra information



## Install 
There are two possible ways to install the program.

1. *Install from the Precompiled Execulables*.
Just copy the binary file to the folder where the executables are located, such as `/usr/local/bin`, and rename the binary file as `reoa`. It is recommended to use `reoa_linux64_gomp4` which has the OpenMP 4.0 support.  See `README.md` under `bin` for further details.  

2. *Compile from the Source Files*. 
Under the `src` folder, run `make` command to compile the sources and run `make install` to install the executables. The latter operation may require adminstation privilege (aka. `root`). (If `make install` does not work, just copy the compiled `reoa` file and other executables to `/usr/local/bin` or other executable folder). If the compilation does not finish successfully, you may need to change the settings in `Makefile` manually. See `README.md` under `src` for further details. 

## License 
This program is free software; it is released under the GNU GENERAL PUBLIC LICENSE Version 3. 

## Contact Us

If you would like to receive updates from us regarding bug fixes, patches, feature updates or if you would like to contact us, please write to [1353023@qq.com](1353023@qq.com) or contact us via QQ or WeChat: 1353023.


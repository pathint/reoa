# Rank-based Gene Pair Analysis (rGPA)
This program implements the original RankComp and the RankCompV2 algorithms to detect the dysregulated genes. It is applicable to gene expression, methylation and other gene microarray data sets. Both large-scale cohort and individual samples or small-scale cell-line samples with two or three technical replicates can be analyzed. 

If you publish results from using this progam, please cite the following paper.

Hongwei Wang, Qiang Sun, Wenyuan Zhao, Lishuang Qi, Yunyan Gu, Pengfei Li, Mengmeng Zhang, Yang Li, Shu-Lin Liu, and Zheng Guo. (**2015**). Individual-level analysis of differential expression of genes and pathways for personalized medicine. *Bioinformatics*. 31(1):62-8. [DOI:10.1093/bioinformatics/btu522](http://dx.doi.org/10.1093/bioinformatics/btu522).


## Algorithm

### Significantly Stable Gene Pairs
A pair of genes, {*a*, *b*}, is considered as statistically stable if they hold the same order relationship (*a* < *b* or *a* >*b*) in most of the samples. Binomial distribution is used to calculate the *P* value under the null hypothesis (a and b does not have a stable order relation) for the large-scale samples. The [Benjamini–Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate) is used to control the false discovery rate (FDR) at level alpha, which is 0.05 by default. For small-scale samples, such as technical replicates in the cell-line studies, the exception number is used to screen the stable pairs. For example, if the exception number is set as 0, the pair holds the same order relation in all the samples. 

*Note on Notation Abuse*: *a*, *b* are used to represent both genes and the expression values (or methylation level etc) of the genes. Readers should be able to discern what the notations represent with no trouble. 

### Concordant and Reversal Gene Pairs
If a pair of genes, {*a*, *b*}, is significantly stable in both the normal group and the disease group and the orders are the same (either *a* < *b* or *a* > *b*), the pair is a concordant pair. Likewise, if a pair is stable in both groups, but the orders are different (*a* < *b* in one group, but *a* > *b* in the other group), the pair is called reversal pair.  Concordant and reversal gene pairs are the basic data to detect dysregulated genes.   


### Dysregulated Genes, the Original RankComp Algorithm
Whether a gene is dysregulated or not in a disease cohort compared with normal cohort is identified through the original RankComp algorithm or through the improved RankCompV2 algorithm. Here is a breif description on the algorithms. Interested users are advised to read the related papers.

For one gene, *a*, we count the numbers of genes whose expression (or methylation etc) levels are higher and lower than this gene in the normal cohort and in the disease cohort. The counting is carried out for the concordant and reversal gene pairs only. Thus, we obtain the following contingency table.

Group | the number of genes whose expression levels are higher than *a* | the number of genes whose expression levels are less than *a* 
---- | ------------ | -------------
Normal | *n*<sub>g</sub> | *n*<sub>l</sub>
Disease |*d*<sub>g</sub> | *d*<sub>l</sub>

[Fisher exact test](https://en.wikipedia.org/wiki/Fisher's_exact_test) is then used to calculate the *P* value under the null hypothesis which states that the number of genes whose expression levels are higher or lower has no association with disease state (absence or presence). Then the [Benjamini–Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate) is used to control the false discovery rate (FDR) at level alpha, which is 0.05 by default. If the null hypothesis is rejected, the gene *a* is identified as a dysregulated gene. The dysregulation direction is judged by comparing *n*<sub>g</sub> / *n*<sub>l</sub> with *d*<sub>g</sub> / *d*<sub>l</sub>. 

The above describes the first step in both the algorithms. In the original RankComp algorithm, one further step is carried out. We redo the counting based on the reversal pairs. For gene *a*, we count the number of genes whose expression level is greater than that of *a* in the normal group but becomes less than that of *a* in the disease group and the number of genes whose expression level is less than that of *a* in the normal group but becomes greater than that of *a* in the disease group. The numbers are labeled as *g2l* and *l2g*, respectively. After the initial step, each gene has been assigned to one of three possible states, up-regulated, down-regulated and not dysregulated. In the second step, for a reversal pair, *a* < *b* in the normal group and *a* > *b* in the disease group, if *a* is up-regulated, *a* will not be counted as *g2l* for *b* and if *b* is not down-regulated, *b* is not counted as *l2g* for *a*. This is to say, the reversed order may be due to the up-regulation of gene *a* or the down-regulation of gene *b* alone.  

The recounted values *g2l* and *l2g*, and the values of *n*<sub>g</sub> and *n*<sub>l</sub> from the first step, are used to construct the above contingency table. Fisher exact test is again used to identify the dysregulated genes.         


### Dysregulated Genes, the RankCompV2 Algorithm
Because of the low statistical power of the original RankComp algorithm, we further developed the RankCompV2 algorithm.  In the RankCompV2 algorithm, the first step as described above is repeated until the number of up or down-regulated genes no longer changes. In the first step, all the genes are assumed to be not dysregulated and the contingency table is constructed based on all the concordant and reversal pairs. In the later cycles, if one gene of a pair is determined to be dysregulated in the previous cycle, the gene will not be counted for the other gene of the pair.        


## Usage
how to use.

### Scenario 1
description

### Scenario 2
description

### Scenario 3
description

## Install 
There are two possible ways to install the program.

1. *Install from the Precompiled Execulables*.
Just copy the binary file to the folder where the executables are located, such as `/usr/local/bin`. It is recommended to use `rgpa_linux64_gomp4` which has the OpenMP 4.0 support.  See `README.md` under `bin` for further details.  

2. *Compile from the Source Files*. 
Under the `src` folder, run `make` command. Then copy the compiled `rgpa` to `/usr/local/bin` or other executable folder. If the compile does not finish successfully, you may need to change the settings in `Makefile` manually. See `README.md` under `src` for further details. 

## License 
This program is free software; it is released under the GNU GENERAL PUBLIC LICENSE Version 3. 

## Contact Us

If you would like to receive updates from us regarding bug fixes, patches, feature updates or if you would like to contact us, please write to [1353023@qq.com](1353023@qq.com) or contact us via QQ or WeChat: 1353023.



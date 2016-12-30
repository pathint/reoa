# Rank-based Gene Pair Analysis (rGPA)
This program implements the original RankComp and the RankCompV2 algorithms to detect the dysregulated genes. It is applicable to gene expression, methylation and other gene microarray data sets. Both large-scale cohort and individual samples or small-scale cell-line samples with two or three technical replicates can be analyzed. 

If you publish results from using this progam, please cite the following paper.

Hongwei Wang, Qiang Sun, Wenyuan Zhao, Lishuang Qi, Yunyan Gu, Pengfei Li, Mengmeng Zhang, Yang Li, Shu-Lin Liu, and Zheng Guo. (**2015**). Individual-level analysis of differential expression of genes and pathways for personalized medicine. *Bioinformatics*. 31(1):62-8. [DOI:10.1093/bioinformatics/btu522](http://dx.doi.org/10.1093/bioinformatics/btu522).


## Algorithm

### Significally Stable Gene Pairs
A pair of genes, {*a*, *b*}, is considered as statistically stable if they hold the same order relationship (*a* < *b* or *a* >*b*) in most of the samples. The binomial distribution is used to calculate the *P* value under the null hypothesis (a and b does not have a stable order relation) for the large-scale samples. The [Benjamini–Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate) is used to control the false discovery rate (FDR) at level alpha, which is 0.05 by default. For small-scale samples, such as technial replicates in the cell-line studies, the exception number is used to screen the stable pairs. For example, if the exception number is set as 0, the pair holds the same order relation in all the samples. 

### Concordant and Reversal Gene Pairs
If a pair of genes, {*a*, *b*}, is significanlly stable in both the normal group and the disease group and the orders are the same (either *a* < *b* or *a* > *b*), the pair is a concordant pair. Likewise, if a pair is stable in both groups, but the orders are different (*a* < *b* in one group, but *a* > *b* in the other group), the pair is called reversal pair.  Concordant and reveral gene pairs are the basic data to detect dysregulated genes.   


### Dysregulated Genes, the Original RankComp Algorithm
description

### Dysregulated Genes, the RankCompV2 Algorithm
description

## Usage
how to use.

### Scenario 1
description

### Scenario 2
description

### Scenario 3
description

## Install 
description

```
module XOR(output out1,  input in1, in2);
  wire w1, w2, w3, w4;
  not (w1, in1);
  not (w2, in2);
  not (w3, in1, w2);
  not (w4, in2, w1);
  or (out1, w3, w4);
endmodule
```


## Contact Us
If you would like to receive updates from us regarding bug fixes, patches, feature updates or if you would like to contact us, please write to [1353023@qq.com](1353023@qq.com) or contact us via QQ or WeChat: 1353023.



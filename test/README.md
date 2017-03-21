## Test Cases

### CellComp Test Case
Small-scale cell line data sets (three technical replicates). To test, run the following command.
```
cellcomp cellcomp_control.dat cellcomp_case.dat
```
### OneComp Test Case
Besides the two input files `onecomp_control.dat` and `onecomp_treated.dat` provideded here, there is an additional file `HCT116.txt` which is also required but is too large to be put here. Interested readers may contact the authors for the file. The file `HCT116.txt` contains the background stable gene pairs which were constructed from normal samples with a sufficent large sample size.  

Once the background file is obtained, run the following command to obtain dysregulated genes.
```
onecomp HCT116.txt onecomp_control.dat onecomp_treated.dat
```
The above running mode uses the default background gene pair filtering mode: one control sample is used for filtering to detect DEGs in the corresponding case sample.

The following command, with setting `-m 1`, uses all the control samples to filter the background gene pairs.
```
onecomp -m 1  HCT116.txt onecomp_control.dat onecomp_treated.dat
```

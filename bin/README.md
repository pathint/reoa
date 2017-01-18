## Install from the Binary Executables

1. Copy the executable (`reoa*`) file to the folder where your PATH can find, such as `/usr/local/bin/` and rename it as `reoa`. Copy the scripts (`cellcomp` and others) to the folder where your PATH can find.
2.  Make sure the files have the executable permission.  If it does not, use `chmod 755 reoa` to make the modification.
3. Try it! For example, running the following command will return the usage of the programs. 

```
cp reoa_linux64_gomp4 /usr/local/bin/reoa
cp ../src/cellcomp /usr/local/bin/cellcomp
reoa -h
cellcomp -h
```
 
## About the Binary Files
1. `reoa_linux64_gomp3` is compiled with gcc-4.7.3 with OpenMP Application Program Interface v3.0 enabled.
2. `reoa_linux64_gomp4` is compiled with gcc-4.9.3 with OpenMP Application Program Interface v4.0 enabled. It requires the installation of OpenMP (`libgomp`) v4.0.



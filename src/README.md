## Install from the Source File
1. Run the command `make` under this folder (`src`).
2. Run the command `make install` to install the package. The default installation folder is `/usr/local/bin`. This requires `root` privilege. Alternatively, you may manually copy the executables (`reoa`, `cellcomp` and `onecomp`) to the folder where your PATH can find, such as `/usr/local/bin/`.
3. Make sure the files have the correct executable permission (e.g. 755).  If it does not, use `chmod 755 reoa` to make the modification.
4. Besides the main program `reoa` and specific scenario application scripts (`cellcomp` and `onecomp`), four additional utility programs, `binomial`, `hypergeom`, `fisher` and `random_sample`, should be also generated. They are used to test the statistical functions implemented in the package.  
5. Try it! For example, running the following command will return the usage of the program. 

```
reoa -h
cellcomp -h
onecomp -h
```


## Trouble Shooting
1. Make sure you have the basic compiling tools installed on your system, such as `gcc`, `as` and `make`.
2. Try to replace `gcc`, `CCFLAGS`, and `INCLUDES` in the `Makefile` to fit your system settings if you fail with the default settings. Then run the `make` command.
3. On macOS, the compiler in `Makefile` is set as `gcc-mp-5`, and `INCLUDES` is set as `-I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.2.0/include`. You may need to edit them to fit your system. The default `clang` compiler on macOS seems does not support OpenMP v4.0. Therefore, you need to install `gcc` first or `clang-openmp` using `brew`. 
4. Alternatively, you may install the binary executable directly which are available under the `bin` folder. Please choose the appropriate version that fits your system.  
5. Windows users should be able to compile the source code using **Cygwin** or other similar cross-platform tools.        


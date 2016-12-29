## Install from the Source File
1. Run the command `make` under this folder (`src`).
2. Copy the executable (`rgpa`) file to the folder where your PATH can find, such as `/usr/local/bin/`.
2. Make sure the file have the executable permission.  If it does not, use `chmod +x rgpa` to make the modifications.
3. Try it! For example, running the following command will return the usage of the program. 

```
rgpa -h 
```

## Trouble Shooting
1. Make sure you have the basic tools installed on your system, such as `gcc`, `as` and `make`.
2. Try to replace `gcc`, `CCFLAGS`, and `INCLUDES` in the `Makefile` to fit your system settings. Then run the `make` command.
3. On macOS, the compiler in `Makefile` is set as `gcc-mp-5`, and `INCLUDES` is set as `-I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.2.0/include`. You may need to edit them to fit your system. The default `clang` compiler on macOS seems does not support OpenMP v4.0. Therefore, you may need install `gcc` first.
4. Alternatively, you may install the binary executable directly which are avaiable under the `bin` folder. Please choose the appropriate version that fits your system.  



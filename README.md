# MATLAB experiments with integer-bassed matrix multiply
This repository contains the data plotted in [Sec. 4, 1] and the MATLAB code used to generate it.

The scripts rely on the [gemmi](https://github.com/north-numerical-computing/gemmi) library for simulating integer-based GEMM.

The errors in [data](./data) were produced with MATLAB version R2024a Update 2 and the gemmi version under the SHA [ddc562ea9404b6ced325ef66a6a804ef383b73e5](https://github.com/north-numerical-computing/gemmi/tree/ddc562ea9404b6ced325ef66a6a804ef383b73e5).
Run the script [run_tests.m](./run_test.m) to regenerate the data files in [data](./data) or run the shell script [plot.sh](./plot.sh) to regenerate the data files and plot the figures in [diagrams](./diagrams).

### References

 [1] A. Abdelfattah, J. Dongarra, M. Fasi, M. Mikaitis, and F. Tisseur [*Error Analysis of Floating-Point Matrix Multiplication Computed via Low-Precision Integer Arithmetic*](https://arxiv.org/pdf/2506.11277). arXiv:2506.11277 [math.NA]. June, 2025.

### License

This software is distributed under the terms of the 2-clause BSD software license (see [LICENSE](./LICENSE)).

# Arbitrary Precision Integer BLAS (APIB)
This repository contains the code associated with the bachelor's thesis "Arbitrary Precision Integer BLAS" written by Jonas Kühne at the Advanced Computing Laboratory, Department of Computer Science, ETH Zurich in 2025. The thesis was supervised by Prof. Markus Püschel and co-supervised by Tommaso Pegolotti.

The repository contains the source code of the tool developed during the thesis and instructions on how to use it. All scripts were tested on an Intel Xeon Silver 4410Y (Sapphire Rapids) CPU running Ubuntu 22.04.

## The Library
APIB is a library for the multiplication of unsigned integer matrices of the form $D\leftarrow \alpha\cdot C + \beta\cdot A\times B$, where $\alpha$ and $\beta$ are also unsigned integers. The number of bits that the elements in $A$, $B$, and $C$ take in memory must be given at the time of compilation. The same holds for the sizes of the scalars $\alpha$ and $\beta$. We support sizes of up to 64 bits. Our library supports splitting the matrix $A$ into horizontal slices, each slice having its own number of bits per element. With this additional information, we represent the matrices in a compressed data format, reducing memory traffic during the computation of the matrix-matrix product.

## Prerequisites
### Hardware
Our code requires the following extensions to the x86 ISA: AVX-512, AVX-512IFMA52, AVX-512_VNNI, AMX-TILE, and AMX-INT8.

### OS
Our scripts expect a unix OS and that the g++ compiler and python3 (to create plots) are installed on the system.

### Additional Libraries
| Library | Purpose |
| --- | --- |
| [PAPI](https://github.com/icl-utk-edu/papi/wiki/Downloading-and-Installing-PAPI) | Conduct Measurements for Benchmarks |
| [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&linux-install=apt), [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) | Comparison Libraries |
| [NumPy](https://numpy.org/), [Matplotlib](https://matplotlib.org/) | Create Plots |

To use APIB without comparing it to other libraries, only [PAPI](https://github.com/icl-utk-edu/papi/wiki/Downloading-and-Installing-PAPI) is required.

### Environment Variables
Our scripts expect the environment variable PAPI_DIR to be set as described in the [instructions](https://github.com/icl-utk-edu/papi/wiki/Downloading-and-Installing-PAPI).

We also use an environment variable EIGEN_DIR that should be set to the directory where Eigen is installed. 

To set the environment variables needed by Intel MKL, the included setup script has to be executed, it is usually located at /opt/intel/oneapi/setvars.sh. If any error mentioning MKL occurs, it often helps to run this script. 

If APIB is meant to be used in your own code, you also need to set the environment variable APIB_DIR to where this directory is located.

## Installation
Due to how we set paths, it is important to run all of the following commands from the directory in which APIB is installed.
To setup APIB, simply run the followoing command:
```
make init
```
If you are using another CPU than the Intel Xeon Silver 4410Y (Sapphire Rapids), we recommend running the following command, which auto-tunes some parameters. It will take a while to complete.
```
make bench_blocks
```
If you are using another CPU than the Intel Xeon Silver 4410Y (Sapphire Rapids) and want to create plots from the benchmark measurements, also make sure to run the following command, which estimated the memory bandwidth of your CPU.
```
make bench_bandwidth
```
To make sure that everything worked, run the following command.
```
make test
```
## Reproduce Results
To reproduce our results, run the following command. It will take some time to complete.
```
make reproduce_results_all
```
If you only want to reproduce the results for some libraries, you can run a subset of the following commands, each benchmarks the library mentioned in the command. The benchmarks of Eigen will take a while to complete.
```
make reproduce_results_apib
make reproduce_results_apib_slice
make reproduce_results_eigen
make reproduce_results_mkl
make reproduce_results_lut
```
The results of the benchmarks will be stored here [here](results/measurements). Before running your own benchmarks, this directory contains our measurements.
To create plots that illustrate these results, use the following command.
```
python3 util/roofline.py
```
The plots will be stored here [here](results/rooflines). Before creating your own plots, this directory contains our plots.

## Using APIB
If you want to use APIB in your code, you can look at the examples [here](examples).
To compile, we recommend the following compiler options. If you do not want to use multiple threads, remove the flag -fopenmp.
```
g++ your_code.cpp -o your_binary -march=native -O3 -fno-strict-aliasing -fno-tree-vectorize -I/${APIB_DIR}/include -L/${APIB_DIR}/lib -lAPIB -fopenmp
```


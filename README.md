# CUDA-Based Parallelization of the Temporal PC-algorithm with Multiple Imputation
The implementation is an extention of the [cuPC](https://github.com/LIS-Laboratory/cupc/tree/master) implementation for multiple imputations of Gaussian data, utilizing the algorithmic procedure from the [gaussMItest](https://github.com/bips-hb/micd/blob/master/R/gaussMItest.R) from the [micd](https://github.com/bips-hb/micd/) library.
It is possible to include tiered background information as seen in [tpc](https://github.com/bips-hb/tpc). 

For installation of dependencies, we refer to the previously mentioned libraries.

## Compilation 
Execute "nvcc -O3 --shared -Xcompiler -fPIC -o SkeletonMI.so cuPC-S-MI.cu" to compile .cu files

## Example
An example of usage of the temporal PC-algorithm with multiple imputations utilizing CUDA-based GPU parallelization can be seen in
```
use_CUTPC-MI.R
```
It is also possible to utilize the version without tiers, an example of this can be seen in 
```
use_CUPC-MI.R
```

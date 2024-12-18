# CUDA-Based Parallelization of the Temporal PC-algorithm with Multiple Imputation
This implementation is made for the Master Thesis _Scalable Constraint-Based Causal Discovery with Multiple Imputation for Incomplete Data_ by Frederik Fabricius-Bjerre, December 2024. 

The implementation is an extension of the [cuPC](https://github.com/LIS-Laboratory/cupc/tree/master) implementation for multiple imputations of Gaussian data, utilizing the algorithmic procedure from the [gaussMItest](https://github.com/bips-hb/micd/blob/master/R/gaussMItest.R) from the [micd](https://github.com/bips-hb/micd/) library.
It is possible to include tiered background information as seen in [tpc](https://github.com/bips-hb/tpc). 

For installation of dependencies, we refer to the previously mentioned libraries.

## Compilation 
Execute ```nvcc -O3 --shared -Xcompiler -fPIC -o SkeletonMI.so cuPC-S-MI.cu``` to compile .cu files

## Example
An example of usage of the temporal PC-algorithm with multiple imputations utilizing CUDA-based GPU parallelization can be seen in
```
use_CUTPC-MI.R
```
It is also possible to utilize the version without tiers, an example of this can be seen in 
```
use_CUPC-MI.R
```

## References 
### tpc, micd
Foraita R, Friemel J, Günther K, Behrens T, Bullerdiek J, Nimzyk R, Ahrens W, Didelez V (2020). Causal discovery of gene regulation with incomplete data. Journal of the Royal Statistical Society: Series A (Statistics in Society), 183(4), 1747-1775. URL https://rss.onlinelibrary.wiley.com/doi/10.1111/rssa.12565.

Foraita R, Witte J, Börnhorst C, Gwozdz W, Pala V, Lissner L, Lauria F, Reisch L, Molnár D, De Henauw S, Moreno L, Veidebaum T, Tornaritis M, Pigeot I, Didelez V. A longitudinal causal graph analysis investigating modifiable risk factors and obesity in a European cohort of children and adolescents. 2022; medRxiv [2022.05.18.22275036](https://www.medrxiv.org/content/10.1101/2022.05.18.22275036v1).

Witte J, Foraita R, Didelez V (2022). Multiple imputation and test-wise deletion for causal discovery with incomplete cohort data. Statistics in Medicine, <https://doi.org/10.1002/sim.9535>.


### pcalg

Markus Kalisch, Martin Mächler, Diego Colombo, Marloes H. Maathuis, Peter Bühlmann (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47(11), 1-26. URL [www.jstatsoft.org/v47/i11/](https://www.jstatsoft.org/article/view/v047i11).

Alain Hauser, Peter Bühlmann (2012). Characterization and greedy learning of interventional Markov equivalence classes of directed acyclic graphs. Journal of Machine Learning Research, 13, 2409-2464. URL <https://www.jmlr.org/papers/v13/hauser12a.html>.

### cuPC
Behrooz Zarebavani, Foad Jafarinejad, Matin Hashemi, Saber Salehkaleybar (2020). {cuPC}: CUDA-based Parallel PC Algorithm for Causal Structure Learning on GPU, IEEE Transactions on Parallel and Distributed Systems (TPDS), Vol. 31, No. 3, March 2020. Url: <https://ieeexplore.ieee.org/document/8823064>, 
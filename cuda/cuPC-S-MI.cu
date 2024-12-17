#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "gpuerrors.h"
#include "cuPC-S-MI.h"
#include "MI-p-val.cuh"

//========================> Main Function Parameter <========================
//Description : this function just calculate one Stage of PC stable algorithm
//@param C          = Correlation matrix
//@param VarSize    = Number of Nodes in Dataset
//@param Stage      = Number of Neighbor in each dimension of Neighbor Matrix
//@param G          = Is the Graph array
//@param Alpha      = The alpha significance level for deleting each edge
//@param Nbr        = Neighbor Matrix with format of:
//[i , j , |Neighbor idx 1|,|Neighbor idx 2| , ...]
//@param Nrow       = Number Of row in Nbr matrix
//@param Ncol       = Number of Col in Nbr matrix
//============================================================================

// global df method constant
__constant__ int df_method;

void SkeletonMI(double* C, int *P, int *Nrows, int *m, int *G, double *Alpha, int *l, int *maxlevel, double *pMax, int* SepSet, int* tiers, int* DF_method)
{
    double *C_cuda;         // copy of C array in GPU
    double *pMax_cuda;
    int    *G_cuda;         // copy of G Array in GPU
    int    *nprime_cuda;
    int    *SepSet_cuda;
    int    *GPrime_cuda;
    int    *mutex_cuda;
    int    *tiers_cuda;

    int    n = *P;
    int    nrows = *Nrows;
    int    M = *m;
    double alpha = *Alpha;
	int    nprime = 0;

    // copy df_method as global variable
    cudaMemcpyToSymbol(df_method, DF_method, sizeof(int));

    dim3   BLOCKS_PER_GRID;
    dim3   THREADS_PER_BLOCK;
    
    bool    FinishFlag = false;

    *l = 0;
    HANDLE_ERROR( cudaMalloc((void**)&mutex_cuda,  n * n * sizeof(int)) );
    HANDLE_ERROR( cudaMalloc((void**)&nprime_cuda,     1 * sizeof(int)) );
    HANDLE_ERROR( cudaMalloc((void**)&SepSet_cuda,  n * n * ML * sizeof(int)) );
    HANDLE_ERROR( cudaMalloc((void**)&GPrime_cuda,     n * n * sizeof(int)) );
    HANDLE_ERROR( cudaMalloc((void**)&C_cuda,     n * n * M * sizeof(double)) );
    HANDLE_ERROR( cudaMalloc((void**)&G_cuda,     n * n * sizeof(int)) );
    HANDLE_ERROR( cudaMalloc((void**)&tiers_cuda,     n * sizeof(int)) );
    HANDLE_ERROR( cudaMalloc((void**)&pMax_cuda,  n * n * sizeof(double)) );
    
    // copy correlation matrix from CPU to GPU
    HANDLE_ERROR( cudaMemcpy(C_cuda, C,       n * n * M * sizeof(double), cudaMemcpyHostToDevice) );
    HANDLE_ERROR( cudaMemcpy(tiers_cuda, tiers,       n * sizeof(int),    cudaMemcpyHostToDevice) );
    // initialize a 0 matrix 
    HANDLE_ERROR( cudaMemset(mutex_cuda, 0, n * n * sizeof(int)) );
    
    CudaCheckError();
    //----------------------------------------------------------
    for (*l = 0; *l <= ML && !FinishFlag && *l <= *maxlevel; *l = *l + 1){
        if (*l == 0){
            BLOCKS_PER_GRID = dim3(n * n, 1, 1);
            THREADS_PER_BLOCK = dim3(ML, 1, 1);
            SepSet_initialize<<< BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(SepSet_cuda, n);
            CudaCheckError();
            BLOCKS_PER_GRID = dim3(n, n, 1);
            THREADS_PER_BLOCK = dim3(1, 1, 1);
            pMax_initialize<<< BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(pMax_cuda, n);
            CudaCheckError();
            if ( (n * n) < 1024) {
                // BLOCKS_PER_GRID   = dim3( 1, 1 ,1);
                // THREADS_PER_BLOCK = dim3(32, 32, 1);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL1, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL1, 1, 1);
                cal_Indepl0 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>  (C_cuda, G_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else {
                // BLOCKS_PER_GRID   = dim3(ceil( ( (double) (n)) / 32.0), ceil( ( (double) (n)) / 32.0), 1);
                // THREADS_PER_BLOCK = dim3(32, 32, 1);
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL1, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL1, 1, 1);
                cal_Indepl0 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK >>> (C_cuda, G_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
        } else {
            //================================> Start Scan Process <===============================
            HANDLE_ERROR( cudaMemset(nprime_cuda, 0, 1 * sizeof(int)) );
            BLOCKS_PER_GRID = dim3(1, n, 1);
            THREADS_PER_BLOCK = dim3(1024, 1, 1);
            scan_compact <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, n * sizeof(int) >>> (GPrime_cuda, G_cuda, n, nprime_cuda);
            CudaCheckError();
            HANDLE_ERROR( cudaMemcpy(&nprime, nprime_cuda, 1 * sizeof(int), cudaMemcpyDeviceToHost) );

            //================================> Begin The Gaussian CI Test  <==============================
            // check whether a CI test is possible
            if (nprime - 1 < *l){ // if not:
                *l = *l - 1;
                FinishFlag = true;
                break;
            }

            if (*l == 1){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL1, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL1, 1, 1);
                cal_Indepl1 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda, G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
                HANDLE_ERROR( cudaDeviceSynchronize() ) ;
                CudaCheckError();
            }
            else if (*l == 2){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL2, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL2, 1, 1);
                cal_Indepl2 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda, G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 3){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL3, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL3, 1, 1);
                cal_Indepl3 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 4){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL4, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL4, 1, 1);
                cal_Indepl4 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 5){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL5, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL5, 1, 1);
                cal_Indepl5 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 6){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL6, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL6, 1, 1);
                cal_Indepl6 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 7){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL7, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL7, 1, 1);
                cal_Indepl7 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 8){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL8, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL8, 1, 1);
                cal_Indepl8 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 9){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL9, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL9, 1, 1);
                cal_Indepl9 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 10){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL10, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL10, 1, 1);
                cal_Indepl10 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 11){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL11, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL11, 1, 1);
                cal_Indepl11 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 12){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL12, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL12, 1, 1);
                cal_Indepl12 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 13){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL13, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL13, 1, 1);
                cal_Indepl13 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            }
            else if(*l == 14){
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeL14, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenL14, 1, 1);
                cal_Indepl14 <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M);
                CudaCheckError();
            } else{
                // if l > 14, we call something that takes up more memory, but essentially is the same. 
                // this works up to l = ML = 32. This can be changed if the machine can handle it. 
                BLOCKS_PER_GRID = dim3(NumOfBlockForEachNodeLAbove14, n, 1);
                THREADS_PER_BLOCK = dim3(ParGivenLAbove14, 1, 1);
                cal_Indep <<< BLOCKS_PER_GRID, THREADS_PER_BLOCK, nprime * sizeof(int) >>> (C_cuda,  G_cuda, GPrime_cuda, mutex_cuda, SepSet_cuda, pMax_cuda, tiers_cuda, alpha, n, nrows, M, *l);
                CudaCheckError();
            }
        }
    } // if l > 0

    // copy Graph G from GPU to CPU
    HANDLE_ERROR( cudaMemcpy(G, G_cuda, n * n * sizeof(int), cudaMemcpyDeviceToHost) );
    // copy separation set from GPU to CPU
    HANDLE_ERROR( cudaMemcpy(SepSet, SepSet_cuda,   n * n * ML * sizeof(int), cudaMemcpyDeviceToHost) );  
    // copy  Pmax from GPU to CPU
    HANDLE_ERROR( cudaMemcpy(pMax, pMax_cuda, n * n * sizeof(double), cudaMemcpyDeviceToHost) );
    
    // preprocess pMax
    double temp = 0;
    for (int i = 0; i < n; i++){
        pMax[i * n + i] = 1;
        for(int j = (i + 1); j < n; j++){
            if(G[i * n + j] == 0){
                temp = fmax(pMax[j * n + i], pMax[i * n + j]);
                pMax[j * n + i] = temp;
                pMax[i * n + j] = temp;
            }
            else{
                pMax[j * n + i] = -100000;
                pMax[i * n + j] = -100000;
            }
            
        }
    }
    // free allocated space
    HANDLE_ERROR( cudaFree(SepSet_cuda) );
    HANDLE_ERROR( cudaFree(C_cuda) );
    HANDLE_ERROR( cudaFree(GPrime_cuda) );
    HANDLE_ERROR( cudaFree(G_cuda) );
    HANDLE_ERROR( cudaFree(mutex_cuda) );
    HANDLE_ERROR( cudaFree(pMax_cuda) );
    HANDLE_ERROR( cudaFree(tiers_cuda) );
}// SkeletonMI


__global__ void SepSet_initialize(int *SepSet, int size){
    int row = bx;
    SepSet[row * ML + tx] = -1;
}
__global__ void pMax_initialize(double *pMax, int size){
    int row = bx;
    pMax[row * size + tx] = -1;
}


__global__ void cal_Indepl0(
    double *C,       // correlation matrices  
    int *G,          // adjacency matrix        
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values     
    int *tiers,      // tiers vector      
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations               
)
{
    // initialize variables
    int row = blockDim.x * bx + tx;
    int col = blockDim.y * by + ty;
    double z_m[MAX_M];
    double p_val;
    int ord = 0;

    if(row < col && col < n){
        // loop over all M imputations
        for (int m = 0; m < M; m++) {
            // compute the index into the 1D C array
            int C_index = m * n * n + row * n + col;
            // compute the correlation coefficient for this imputation
            double rho_m = C[C_index];
            // compute Fisher's Z-transformation
            rho_m = 0.5 * log((1 + rho_m) / (1 - rho_m));
            
            z_m[m] = rho_m;
        }
        // compute MI p-value
        p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
        if (p_val >= alpha) {
            // asign values to pMax and remove edges
            pMax[row * n + col] = p_val;
            G[row * n + col] = 0;
            G[col * n + row] = 0;
            // add placeholder value, so we now sepset is the empty set
            Sepset[(row * n + col) * ML] = 0; 
        } else {
            G[row * n + col] = 1;
            G[col * n + row] = 1;
        }
    }
    if (row == col && col < n){
        // remove self-loops
        G[row * n + col] = 0;
        G[col * n + row] = 0; 
    }
}


__global__ void cal_Indepl1(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{
    // initialize variables
    int YIdx;
    int XIdx = by;              
    int NbrIdxPointer;
    int NbrIdx;
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];          
    double M0_m;
    double M1[2];               
    double H[2][2];             
    double rho_m;
    double p_val;
    int ord = 1;

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];  // number of neighbours for node X
    if ((SizeOfArr % ParGivenL1) == 0) {
        NumberOfJump = SizeOfArr / ParGivenL1;
    } else {
        NumberOfJump = SizeOfArr / ParGivenL1 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++) {
        if ((tx + cnt * ParGivenL1) < SizeOfArr) {
            G_Chunk[tx + cnt * ParGivenL1] = GPrime[XIdx * n + tx + cnt * ParGivenL1];
        }
        __syncthreads();
    }

    // calculate the number of iterations over neighbor chunks
    if ((SizeOfArr % (ParGivenL1 * NumOfBlockForEachNodeL1)) == 0) {
        NumOfGivenJump = SizeOfArr / (ParGivenL1 * NumOfBlockForEachNodeL1);
    } else {
        NumOfGivenJump = SizeOfArr / (ParGivenL1 * NumOfBlockForEachNodeL1) + 1;
    }
    __syncthreads();
    // main loop
    for (int d1 = 0; d1 < NumOfGivenJump; d1++) {
        __syncthreads();
        if (NoEdgeFlag == 1) {
            return;
        }
        __syncthreads();
        NbrIdxPointer = tx + bx * ParGivenL1 + d1 * ParGivenL1 * NumOfBlockForEachNodeL1;
        NoEdgeFlag = 1;
        __syncthreads();

        if (NbrIdxPointer < SizeOfArr) {
            NbrIdx = G_Chunk[NbrIdxPointer];
            
            // loop over neighbors
            for (int d2 = 0; d2 < SizeOfArr; d2++) {
                if (d2 == NbrIdxPointer) {
                    continue;
                }
                YIdx = G_Chunk[d2];
                
                // tiers constraint
                if(tiers[NbrIdx] > tiers[XIdx]){
                    // skip, since the tier constraint is violated
                    continue;
                }
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // correlation coefs
                        M1[0] = C[m * n * n + XIdx * n + NbrIdx];  
                        M0_m  = C[m * n * n + XIdx * n + YIdx];    
                        M1[1] = C[m * n * n + YIdx * n + NbrIdx];  

                        // compute H
                        H[0][0] = 1.0  - (M1[0] * M1[0]);
                        H[0][1] = M0_m - (M1[0] * M1[1]);
                        H[1][1] = 1.0  - (M1[1] * M1[1]);

                        // compute partial correlation
                        rho_m = H[0][1] / (sqrt(fabs(H[0][0] * H[1][1])));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // compute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);

                    if (p_val >= alpha) {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0) {
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            Sepset[(XIdx * n + YIdx) * ML] = NbrIdx + 1;
                        }
                    }
                }
            }
        }
    }
}



__global__ void cal_Indepl2(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{
    // initialize variables
    int YIdx;
    int XIdx = by;              
    int NbrIdxPointer[2];
    int NbrIdx[2];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];         
    double M0_m;
    double M1[2][2];           
    double M2[2][2];           
    double M2Inv[2][2];        
    double M1MulM2Inv[2][2];   
    double H[2][2];            
    double rho_m;
    double p_val;
    int ord = 2;

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1];  // number of neighbours for node X

    // if number of neighbours are smaller than potential sepset size
    if (SizeOfArr <= 2){
        return;
    }

    // calculate the number of chunks to process neighbors
    if ((SizeOfArr % ParGivenL2) == 0) {
        NumberOfJump = SizeOfArr / ParGivenL2;
    } else {
        NumberOfJump = SizeOfArr / ParGivenL2 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++) {
        if ((tx + cnt * ParGivenL2) < SizeOfArr) {
            G_Chunk[tx + cnt * ParGivenL2] = GPrime[XIdx * n + tx + cnt * ParGivenL2];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose 2)
    BINOM(SizeOfArr, 2, &NumOfComb);

    // calculate the number of iterations over neighbor combinations
    if ((NumOfComb % (ParGivenL2 * NumOfBlockForEachNodeL2)) == 0) {
        NumOfGivenJump = NumOfComb / (ParGivenL2 * NumOfBlockForEachNodeL2);
    } else {
        NumOfGivenJump = NumOfComb / (ParGivenL2 * NumOfBlockForEachNodeL2) + 1;
    }

    __syncthreads();
    // main loop
    for (int d1 = 0; d1 < NumOfGivenJump; d1++) {
        __syncthreads();
        if (NoEdgeFlag == 1) {
            return;
        }
        int combIdx = tx + bx * ParGivenL2 + d1 * ParGivenL2 * NumOfBlockForEachNodeL2;
        if (combIdx < NumOfComb) {
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();

            // get the indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 2, combIdx + 1);
            NbrIdx[0] = G_Chunk[NbrIdxPointer[0] - 1];
            NbrIdx[1] = G_Chunk[NbrIdxPointer[1] - 1];

            // loop over neighbors
            for (int d2 = 0; d2 < SizeOfArr; d2++) {
                if ((d2 == (NbrIdxPointer[0] - 1)) || (d2 == (NbrIdxPointer[1] - 1))) {
                    continue;
                }
                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {

                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx]; 
                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // update M2 matrix with current imputation
                        M2[0][0] = 1.0;
                        M2[0][1] = C[m * n * n + NbrIdx[0] * n + NbrIdx[1]]; 
                        M2[1][0] = M2[0][1];
                        M2[1][1] = 1.0;

                        // compute the inverse of M2 for current imputation
                        pseudoinversel2(M2, M2Inv);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++) {
                            for (int c2 = 0; c2 < 2; c2++) {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < 2; c3++) {
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                                }
                            }
                        }

                        // compute H
                        H[0][0] = 1.0 - (M1MulM2Inv[0][0] * M1[0][0] + M1MulM2Inv[0][1] * M1[0][1]);
                        H[0][1] = M0_m - (M1MulM2Inv[0][0] * M1[1][0] + M1MulM2Inv[0][1] * M1[1][1]);
                        H[1][1] = 1.0 - (M1MulM2Inv[1][0] * M1[1][0] + M1MulM2Inv[1][1] * M1[1][1]);

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // compute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha) {
                        if (atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0) {
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            Sepset[(XIdx * n + YIdx) * ML] = NbrIdx[0] + 1;
                            Sepset[(XIdx * n + YIdx) * ML + 1] = NbrIdx[1] + 1;
                        }
                    }
                }
            }
        }
    }
}


__global__ void cal_Indepl3(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = by;              
    int NbrIdxPointer[3];
    int NbrIdx[3];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][3];
    double M2[3][3];
    double M2Inv[3][3];
    double M1MulM2Inv[2][3];
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 3;

    NoEdgeFlag = 0;
    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbours for node X
    if (SizeOfArr <= 3){
        return;
    }
    // calculate the number of chunks to process neighbors
    if ((SizeOfArr % ParGivenL3) == 0) {
        NumberOfJump = SizeOfArr / ParGivenL3;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL3 + 1;
    }
    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( tx + cnt * ParGivenL3 ) < SizeOfArr){
            G_Chunk[ tx + cnt * ParGivenL3 ] =  GPrime[ XIdx * n + tx + cnt * ParGivenL3];
        }
        __syncthreads();
    }
    // calculate the number of combinations (n choose 2)
    BINOM(SizeOfArr, 3, &NumOfComb);
    if( (NumOfComb % (ParGivenL3 * NumOfBlockForEachNodeL3)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL3 * NumOfBlockForEachNodeL3);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL3 * NumOfBlockForEachNodeL3) + 1;
    }

    __syncthreads();
    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = tx + bx * ParGivenL3 + d1 * ParGivenL3 * NumOfBlockForEachNodeL3;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();

            // get the indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 3, tx + bx * ParGivenL3 + d1 * ParGivenL3 * NumOfBlockForEachNodeL3 + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {    
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx];
                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }
                        
                        // update M2 matrix with current imputation
                        M2[0][0] = 1.0;
                        M2[0][1] = C[m * n * n + NbrIdx[0] * n + NbrIdx[1]];
                        M2[0][2] = C[m * n * n + NbrIdx[0] * n + NbrIdx[2]];
                        M2[1][0] = M2[0][1];
                        M2[1][1] = 1.0;
                        M2[1][2] = C[m * n * n + NbrIdx[1] * n + NbrIdx[2]];
                        M2[2][0] = M2[0][2];
                        M2[2][1] = M2[1][2];
                        M2[2][2] = 1.0;
                        
                        // compute the pseudoinverse of M2 for current imputation
                        pseudoinversel3(M2, M2Inv);
                        
                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        
                        // compute H matrix
                        H[0][0]   = 1  - H[0][0];
                        H[0][1]   = M0_m - H[0][1];
                        H[1][1]   = 1  - H[1][1];

                        // compute partial correlation
                        rho_m     =  H[0][1] / ( sqrt( fabs(H[0][0] * H[1][1]) ) );
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // compute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha) {
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock        
                            // update G and pMax                
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl4(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[4];
    int NbrIdx[4];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][4];
    double M2[4][4];
    double M2Inv[4][4];
    double M1MulM2Inv[2][4];
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 4;               

    // initialize variables for SVD pseudoinverse
    double v[4][4];
    double w[4], rv1[4];
    double res1[4][4];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 4){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL4) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL4;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL4 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL4 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL4 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL4];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 4, &NumOfComb);
    if( (NumOfComb % (ParGivenL4 * NumOfBlockForEachNodeL4)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL4 * NumOfBlockForEachNodeL4);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL4 * NumOfBlockForEachNodeL4) + 1;
    }
    __syncthreads();

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL4 + d1 * ParGivenL4 * NumOfBlockForEachNodeL4;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();

            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 4, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {    
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {

                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx]; 
                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel4(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }

                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


__global__ void cal_Indepl5(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{
    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[5];
    int NbrIdx[5];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][5];
    double M2[5][5];
    double M2Inv[5][5];
    double M1MulM2Inv[2][5];
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 5;              

    // initialize variables for SVD pseudoinverse
    double v[5][5];
    double w[5], rv1[5];
    double res1[5][5];
    
    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 5){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL5) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL5;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL5 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL5 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL5 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL5];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 5, &NumOfComb);
    if( (NumOfComb % (ParGivenL5 * NumOfBlockForEachNodeL5)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL5 * NumOfBlockForEachNodeL5);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL5 * NumOfBlockForEachNodeL5) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL5 + d1 * ParGivenL5 * NumOfBlockForEachNodeL5;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 5, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;
                
                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;

                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {    
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        
                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx]; 
                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel5(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }

                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation 
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


__global__ void cal_Indepl6(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{
    // initialize variables
    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[6];
    int NbrIdx[6];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];          
    double M0_m;
    double M1[2][6];           
    double M2[6][6];           
    double M2Inv[6][6];        
    double M1MulM2Inv[2][6];   
    double H[2][2];            
    double rho_m;
    double p_val;
    int ord = 6;               

    // initialize pseudoinverse variables
    double v[6][6];
    double w[6], rv1[6];
    double res1[6][6];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 6){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL6) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL6;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL6 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL6 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL6 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL6];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 6, &NumOfComb);
    if( (NumOfComb % (ParGivenL6 * NumOfBlockForEachNodeL6)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL6 * NumOfBlockForEachNodeL6);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL6 * NumOfBlockForEachNodeL6) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL6 + d1 * ParGivenL6 * NumOfBlockForEachNodeL6;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();

            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 6, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;
                
                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {    
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                       
                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx]; 
                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel6(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }

                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            } 
                        }
                    }
                }
            }         
        }
    }
}


__global__ void cal_Indepl7(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{
    // initialize variables
    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[7];
    int NbrIdx[7];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][7];
    double M2[7][7];
    double M2Inv[7][7];
    double M1MulM2Inv[2][7];
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 7;               

    // initialize variables for SVD pseudoinverse
    double v[7][7];
    double w[7], rv1[7];
    double res1[7][7];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 7){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL7) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL7;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL7 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL7 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL7 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL7];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 7, &NumOfComb);
    if( (NumOfComb % (ParGivenL7 * NumOfBlockForEachNodeL7)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL7 * NumOfBlockForEachNodeL7);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL7 * NumOfBlockForEachNodeL7) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL7 + d1 * ParGivenL7 * NumOfBlockForEachNodeL7;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 7, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;
                
                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {    
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        
                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx]; 
                        
                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel7(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }

                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // compute MI p-value                   
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl8(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[8];
    int NbrIdx[8];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][8];
    double M2[8][8];
    double M2Inv[8][8];
    double M1MulM2Inv[2][8];
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 8;               

    // initialize variables for SVD pseudoinverse
    double v[8][8];
    double w[8], rv1[8];
    double res1[8][8];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 8){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL8) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL8;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL8 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL8 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL8 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL8];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 8, &NumOfComb);
    if( (NumOfComb % (ParGivenL8 * NumOfBlockForEachNodeL8)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL8 * NumOfBlockForEachNodeL8);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL8 * NumOfBlockForEachNodeL8) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL8 + d1 * ParGivenL8 * NumOfBlockForEachNodeL8;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 8, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;
                
                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {    
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx]; 

                        // compute M1 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel8(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }   
                        // compute H matrix
                        H[0][0]   = 1.0  - H[0][0];
                        H[0][1]   = M0_m - H[0][1];
                        H[1][1]   = 1.0  - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));  

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }       
        }
    }
}


__global__ void cal_Indepl9(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[9];
    int NbrIdx[9];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][9];
    double M2[9][9];
    double M2Inv[9][9];
    double M1MulM2Inv[2][9];
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 9;               

    // initialize variables for SVD pseudoinverse
    double v[9][9];
    double w[9], rv1[9];
    double res1[9][9];
    
    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 9){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL9) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL9;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL9 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL9 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL9 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL9];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 9, &NumOfComb);
    if( (NumOfComb % (ParGivenL9 * NumOfBlockForEachNodeL9)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL9 * NumOfBlockForEachNodeL9);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL9 * NumOfBlockForEachNodeL9) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL9 + d1 * ParGivenL9 * NumOfBlockForEachNodeL9;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 9, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {

                        // correlation coefs
                        M0_m = C[m * n * n + XIdx * n + YIdx];
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel9(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation 
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


__global__ void cal_Indepl10(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[10];
    int NbrIdx[10];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][10];
    double M2[10][10];          
    double M2Inv[10][10];       
    double M1MulM2Inv[2][10];
    double H[2][2];             // H matrix
    double rho_m;
    double p_val;
    int ord = 10;               

    // initialize variables for SVD pseudoinverse
    double v[10][10];
    double w[10], rv1[10];
    double res1[10][10];
    
    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 10){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL10) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL10;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL10 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL10 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL10 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL10];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 10, &NumOfComb);
    if( (NumOfComb % (ParGivenL10 * NumOfBlockForEachNodeL10)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL10 * NumOfBlockForEachNodeL10);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL10 * NumOfBlockForEachNodeL10) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL10 + d1 * ParGivenL10 * NumOfBlockForEachNodeL10;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 10, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx];

                        // compute M1 matrices
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel10(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


__global__ void cal_Indepl11(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[11];
    int NbrIdx[11];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][11];          
    double M2[11][11];         
    double M2Inv[11][11];      
    double M1MulM2Inv[2][11];  
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 11;              

    // initialize variables for SVD pseudoinverse
    double v[11][11];
    double w[11], rv1[11];
    double res1[11][11];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 11){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL11) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL11;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL11 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL11 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL11 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL11];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 11, &NumOfComb);
    if( (NumOfComb % (ParGivenL11 * NumOfBlockForEachNodeL11)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL11 * NumOfBlockForEachNodeL11);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL11 * NumOfBlockForEachNodeL11) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL11 + d1 * ParGivenL11 * NumOfBlockForEachNodeL11;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 11, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx];

                        // compute M1 matrices
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel11(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


__global__ void cal_Indepl12(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[12];
    int NbrIdx[12];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][12];          
    double M2[12][12];         
    double M2Inv[12][12];      
    double M1MulM2Inv[2][12];  
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 12;              

    // initialize variables for SVD pseudoinverse
    double v[12][12];
    double w[12], rv1[12];
    double res1[12][12];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 12){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL12) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL12;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL12 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL12 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL12 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL12];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 12, &NumOfComb);
    if( (NumOfComb % (ParGivenL12 * NumOfBlockForEachNodeL12)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL12 * NumOfBlockForEachNodeL12);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL12 * NumOfBlockForEachNodeL12) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL12 + d1 * ParGivenL12 * NumOfBlockForEachNodeL12;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 12, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx];

                        // compute M1 matrices
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel12(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl13(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[13];
    int NbrIdx[13];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][13];          
    double M2[13][13];         
    double M2Inv[13][13];      
    double M1MulM2Inv[2][13];  
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 13;              

    // initialize variables for SVD pseudoinverse
    double v[13][13];
    double w[13], rv1[13];
    double res1[13][13];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 13){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL13) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL13;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL13 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL13 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL13 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL13];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 13, &NumOfComb);
    if( (NumOfComb % (ParGivenL13 * NumOfBlockForEachNodeL13)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL13 * NumOfBlockForEachNodeL13);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL13 * NumOfBlockForEachNodeL13) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL13 + d1 * ParGivenL13 * NumOfBlockForEachNodeL13;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 13, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx];

                        // compute M1 matrices
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel13(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indepl14(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M            // number of imputations   
)
{

    int YIdx;
    int XIdx = blockIdx.y;      
    int NbrIdxPointer[14];
    int NbrIdx[14];
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double M1[2][14];          
    double M2[14][14];         
    double M2Inv[14][14];      
    double M1MulM2Inv[2][14];  
    double H[2][2];
    double rho_m;
    double p_val;
    int ord = 14;              

    // initialize variables for SVD pseudoinverse
    double v[14][14];
    double w[14], rv1[14];
    double res1[14][14];

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= 14){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenL14) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenL14;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenL14 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenL14 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenL14 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenL14];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, 14, &NumOfComb);
    if( (NumOfComb % (ParGivenL14 * NumOfBlockForEachNodeL14)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenL14 * NumOfBlockForEachNodeL14);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenL14 * NumOfBlockForEachNodeL14) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenL14 + d1 * ParGivenL14 * NumOfBlockForEachNodeL14;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, 14, combIdx + 1);
            for(int idx = 0; idx < ord; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < ord; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < ord; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx];

                        // compute M1 matrices
                        for (int c1 = 0; c1 < ord; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < ord; c1++){
                            for(int c2 = 0; c2 < ord; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinversel14(M2, M2Inv, v, rv1, w, res1);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < ord; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c2 < ord; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < ord; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, ord, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < ord; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void cal_Indep(
    double *C,       // correlation matrices   
    int *G,          // adjacency matrix      
    int *GPrime,     // number of neighbors for nodes
    int *mutex,      // mutex array for atomic operations 
    int *Sepset,     // separation set matrix 
    double *pMax,    // maximum p-values    
    int *tiers,      // tiers vector
    double alpha,    // significance level                  
    int n,           // number of columns/nodes/variables   
    int nrows,       // number of rows/samples/observations 
    int M,           // number of imputations   
    int order        // order of cond set
)
{
    // check if order exceeds ML
    if (order > ML) {
        // ensure only one thread prints the error message
        if (threadIdx.x == 0 && blockIdx.x == 0 && blockIdx.y == 0) {
            printf("Error: 'order' (%d) exceeds the maximum allowed value of %d.\n", order, ML);
        }
        // terminate the kernel early
        return;
    }


    int YIdx;
    int XIdx = blockIdx.y;      
    int SizeOfArr;
    int NumberOfJump;
    int NumOfGivenJump;
    int NumOfComb;
    __shared__ int NoEdgeFlag;
    extern __shared__ int G_Chunk[];  // shared memory for neighbor indices

    // initialize MI variables
    double z_m[MAX_M];
    double M0_m;
    double rho_m;
    double p_val;
    int NbrIdxPointer[ML];        
    int NbrIdx[ML];                
    double M1[2][ML];              
    double M2[ML][ML];             
    double M2Inv[ML][ML];          
    double M1MulM2Inv[2][ML];      
    double H[2][2];                

    // initialize variables for SVD pseudoinverse
    double v[ML][ML];              
    double w[ML];                  
    double rv1[ML];                
    double res1[ML][ML];                       

    NoEdgeFlag = 0;

    SizeOfArr = GPrime[XIdx * n + n - 1]; // number of neighbors for node X
    if (SizeOfArr <= order){
        return;
    }

    // calculate the number of iterations needed to copy neighbor indices to shared memory
    if( (SizeOfArr % ParGivenLAbove14) == 0 ){
        NumberOfJump = SizeOfArr / ParGivenLAbove14;
    }
    else{
        NumberOfJump = SizeOfArr / ParGivenLAbove14 + 1;
    }

    // copy neighbor indices from global memory to shared memory
    for (int cnt = 0; cnt < NumberOfJump; cnt++){
        if( ( threadIdx.x + cnt * ParGivenLAbove14 ) < SizeOfArr){
            G_Chunk[ threadIdx.x + cnt * ParGivenLAbove14 ] =  GPrime[ XIdx * n + threadIdx.x + cnt * ParGivenLAbove14];
        }
        __syncthreads();
    }

    // calculate the number of combinations (n choose l)
    BINOM(SizeOfArr, order, &NumOfComb);
    if( (NumOfComb % (ParGivenLAbove14 * NumOfBlockForEachNodeLAbove14)) == 0 ){
        NumOfGivenJump = NumOfComb / (ParGivenLAbove14 * NumOfBlockForEachNodeLAbove14);
    }
    else{
        NumOfGivenJump = NumOfComb / (ParGivenLAbove14 * NumOfBlockForEachNodeLAbove14) + 1;
    }

    // main loop
    for(int d1 = 0; d1 < NumOfGivenJump; d1++){
        __syncthreads();
        if(NoEdgeFlag == 1){
            return;
        }
        int combIdx = threadIdx.x + blockIdx.x * ParGivenLAbove14 + d1 * ParGivenLAbove14 * NumOfBlockForEachNodeLAbove14;
        if(combIdx < NumOfComb){
            __syncthreads();
            NoEdgeFlag = 1;
            __syncthreads();
            // get indices of the combination
            IthCombination(NbrIdxPointer, SizeOfArr, order, combIdx + 1);
            for(int idx = 0; idx < order; idx++){
                NbrIdx[idx] = G_Chunk[NbrIdxPointer[idx] - 1];
            }

            // loop over neighbors
            for(int d2 = 0; d2 < SizeOfArr; d2++){
                bool skip = false;
                for (int idx = 0; idx < order; idx++) {
                    if (d2 == (NbrIdxPointer[idx] - 1)) {
                        skip = true;
                        break;
                    }
                }
                if (skip) continue;

                YIdx = G_Chunk[d2];

                // tiers constraint
                bool tiersConstraintViolation = false;
                
                for (int idx = 0; idx < order; idx++){
                    if(tiers[NbrIdx[idx]] > tiers[XIdx]){
                        // skip, since the tier constraint is violated
                        tiersConstraintViolation = true;
                        break;
                    }
                }
                if (tiersConstraintViolation){ 
                    continue;
                }
                
                if (G[XIdx * n + YIdx] == 1) {
                    NoEdgeFlag = 0;

                    // loop over all M imputations
                    for (int m = 0; m < M; m++) {
                        // compute M0_m
                        M0_m = C[m * n * n + XIdx * n + YIdx];

                        // compute M1 matrices
                        for (int c1 = 0; c1 < order; c1++){
                            M1[0][c1] = C[m * n * n + XIdx * n + NbrIdx[c1]];
                            M1[1][c1] = C[m * n * n + YIdx * n + NbrIdx[c1]];
                        }

                        // compute M2 matrix
                        for (int c1 = 0; c1 < order; c1++){
                            for(int c2 = 0; c2 < order; c2++){
                                if(c1 > c2){
                                    M2[c1][c2] = M2[c2][c1];
                                }
                                else if(c1 == c2){
                                    M2[c1][c1] = 1.0;
                                }
                                else{
                                    M2[c1][c2] = C[m * n * n + NbrIdx[c1] * n + NbrIdx[c2]];
                                }
                            }
                        }

                        // compute pseudoinverse of M2
                        pseudoinverse(M2, M2Inv, v, rv1, w, res1, order);

                        // compute M1 * M2Inv
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < order; c2++)
                            {
                                M1MulM2Inv[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < order; c3++)
                                    M1MulM2Inv[c1][c2] += M1[c1][c3] * M2Inv[c3][c2];
                            }
                        }

                        // compute H matrix
                        for (int c1 = 0; c1 < 2; c1++)
                        {
                            for (int c2 = 0; c2 < 2; c2++)
                            {
                                H[c1][c2] = 0.0;
                                for (int c3 = 0; c3 < order; c3++)
                                    H[c1][c2] += M1MulM2Inv[c1][c3] * M1[c2][c3];
                            }
                        }
                        // compute H matrix
                        H[0][0] = 1.0 - H[0][0];
                        H[0][1] = M0_m - H[0][1];
                        H[1][1] = 1.0 - H[1][1];

                        // compute partial correlation
                        rho_m = H[0][1] / sqrt(fabs(H[0][0] * H[1][1]));

                        // fisher's z-transformation 
                        double Z_m = 0.5 * log((1.0 + rho_m) / (1.0 - rho_m));
                        z_m[m] = Z_m;
                    }
                    // comptute MI p-value
                    p_val = compute_MI_p_value(z_m, M, nrows, order, df_method);
                    if (p_val >= alpha){
                        if(atomicCAS(&mutex[XIdx * n + YIdx], 0, 1) == 0){ // lock
                            // update G and pMax
                            G[XIdx * n + YIdx] = 0;
                            G[YIdx * n + XIdx] = 0;
                            pMax[XIdx * n + YIdx] = p_val;
                            // assign sepset (+ 1 since R is one-indexed)
                            for (int idx = 0; idx < order; idx++) {
                                Sepset[(XIdx * n + YIdx) * ML + idx] = NbrIdx[idx] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void Scan (int* G_ScanOut, int* G_ScanIn, int Step, int GSize){
    int index = tx + blockDim.x * bx;
    int row   = by;
    if ( (index < Step) && (index < GSize) ) {
        G_ScanOut[row * GSize + index] = G_ScanIn[row * GSize + index];
    }
    if ( (index >= Step) && (index < GSize)){
        G_ScanOut[row * GSize + index] = G_ScanIn[row * GSize + index] + G_ScanIn[row * GSize + index - Step];
    }
}

__global__ void Compact (int* G_Compact, const int* G, const int* G_ScanRes, int GSize){
    int index = tx + blockDim.x * bx;
    int row   = by;
    int CompactIdx;
    if(index < GSize){
        if( (G[row * GSize + index] == 1) ){
            CompactIdx = G_ScanRes[row * GSize + index] - 1;
            G_Compact[row * GSize + CompactIdx] = index;
        }
        if(index >= G_ScanRes[row * GSize + GSize - 1]){
            if( index != (GSize - 1) ){
                G_Compact[row * GSize + index] = 0;
            }else{
                G_Compact[row * GSize + GSize - 1] = G_ScanRes[row * GSize + GSize - 1];
            }
        }
    }

}

__device__ double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt){
         ct = bt / at;
          result = at * sqrt(1.0 + ct * ct);
    }else if (bt > 0.0) {
        ct = at / bt;
        result = bt * sqrt(1.0 + ct * ct);
    }
    else{
         result = 0.0;
    }
    return(result);
}

__device__ void pseudoinversel2(double M2[][2], double M2Inv[][2])
{
    double A[2][2];
    double M[2][2];
    double L[2][2];
    double newL[2][2];
    double temp[2][2];
    double temp1[2][2];
    double temp2[2][2];
    double temp3[2][2];    

    double tol = 999.99;
    double aux = 0.0;
    double det = 0.0;

    int r = 0;
    int size = 2;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j]     = 0.0;
            M[i][j]     = 0.0;
            L[i][j]     = 0.0;
            newL[i][j]  = 0.0;
            temp[i][j]  = 0.0;
            temp1[i][j] = 0.0;
            temp2[i][j] = 0.0;
            temp3[i][j] = 0.0;
            M2Inv[i][j] = 0.0;
        }
    }
        
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                A[i][j] += M2[i][k] * M2[k][j];
            }
        }
    }

    
    for (int i = 0; i < size; i++) {
        if (tol > A[i][i] && A[i][i] > 0) {
            tol = A[i][i];
        }
    }

    //tol = tol * 1e-9 accroding to paper
    tol = tol * (1e-20);

    for (int k = 0; k < size; k++) {
        
        if (r == 0) {
            for (int i = k; i < size; i++) {
                L[i][r] = A[i][k];
            }
        } else {
            
            for (int i = k; i < size; ++i) {
                for (int l = 0; l < r; l++) {
                    temp[i][k] += L[i][l] * L[k][l];
                }
            }
            for (int i = k; i < size; i++) {
                L[i][r] = A[i][k] - temp[i][k];
            }
        }
        //check with threshold
        if (L[k][r] > tol) {
            L[k][r] = sqrt(L[k][r]);
            if (k < size) {
                for (int i = k + 1; i < size; i++) {
                    L[i][r] = L[i][r] / L[k][r];
                }
            }
        } else {
            r--;
        }
        r++;
    }
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < r; j++) {
            newL[i][j] = L[i][j];
        }
    }

    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < r; j++) {
            for (int k = 0; k < size; k++) {
                M[i][j] += newL[k][i] * newL[k][j];
            }
        }
    }

    /*
    * it's time to compute inv(M) in this stage M is 2*2 so
    * I use close form of 2*2 
    */
   if(r == 1){
        M[0][0] = 1/M[0][0];
   }else if( r == 2){
        aux = 0.0;
        det = 1 / (M[0][0] * M[1][1] - M[0][1] * M[1][0]);
        aux = M[0][0];
        M[0][0] = det * M[1][1];
        M[1][1] = det * aux;
        M[0][1] = (-1 * det) * M[0][1];
        M[1][0] = (-1 * det) * M[1][0];
   }

    
    /*At the final step we must compute L * M * M * L' * G'
     * at first I compute   temp1 = L  * M
     * after that I compute temp2 = L' * G'
     * after that I compute temp3 = M  * temp2
     * finally I compute   output = temp1 * temp3
    */

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < r; j++) {
            for (int k = 0; k < r; k++) {
                temp1[i][j] += newL[i][k] * M[k][j];
            }
        }
    }
   
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                temp2[i][j] += newL[k][i] * M2[k][j];
            }
        }
    }

    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                temp3[i][j] += M[i][k] * temp2[k][j];
            }
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                M2Inv[i][j] += temp1[i][k] * temp3[k][j];
            }
        }
    }
 
}

__device__ void pseudoinversel3(double M2[][3], double M2Inv[][3])
{
    double A[3][3];
    double M[3][3];
    double tempM[3][3];
    double L[3][3];
    double newL[3][3];
    double temp[3][3];
    double temp1[3][3];
    double temp2[3][3];
    double temp3[3][3];

    double tol = 999.99;
    double det = 0.0;

    int r = 0;
    int size = 3;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j]     = 0.0;
            M[i][j]     = 0.0;
            L[i][j]     = 0.0;
            newL[i][j]  = 0.0;
            temp[i][j]  = 0.0;
            temp1[i][j] = 0.0;
            temp2[i][j] = 0.0;
            temp3[i][j] = 0.0;
            M2Inv[i][j] = 0.0;
            tempM[i][j] = 0.0;
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                A[i][j] += M2[i][k] * M2[k][j];
            }
        }
    }


    for (int i = 0; i < size; i++) {
        if (tol > A[i][i] && A[i][i] > 0) {
            tol = A[i][i];
        }
    }

    //tol = tol * 1e-9 accroding to paper
    tol = tol * (1e-20);

    for (int k = 0; k < size; k++) {

        if (r == 0) {
            for (int i = k; i < size; i++) {
                L[i][r] = A[i][k];
            }
        } else {

            for (int i = k; i < size; ++i) {
                for (int l = 0; l < r; l++) {
                    temp[i][k] += L[i][l] * L[k][l];
                }
            }
            for (int i = k; i < size; i++) {
                L[i][r] = A[i][k] - temp[i][k];
            }
        }
        //check with threshold
        if (L[k][r] > tol) {
            L[k][r] = sqrt(L[k][r]);
            if (k < size) {
                for (int i = k + 1; i < size; i++) {
                    L[i][r] = L[i][r] / L[k][r];
                }
            }
        } else {
            r--;
        }
        r++;
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < r; j++) {
            newL[i][j] = L[i][j];
        }
    }


    for (int i = 0; i < r; i++) {
        for (int j = 0; j < r; j++) {
            for (int k = 0; k < size; k++) {
                tempM[i][j] += newL[k][i] * newL[k][j];
            }
        }
    }

    /*
    * it's time to compute inv(M) in this stage M is 2*2 so
    * I use close form of 2*2
    */
   if(r == 1){
        M[0][0] = 1/tempM[0][0];
   }else if( r == 2){
        det = 1 / (tempM[0][0] * tempM[1][1] - tempM[0][1] * tempM[1][0]);
        M[0][0] = det * tempM[1][1];
        M[1][1] = det * tempM[0][0];
        M[0][1] = (-1 * det) * tempM[0][1];
        M[1][0] = (-1 * det) * tempM[1][0];
   }else{
        inverse(tempM, M);
   }


    /*At the final step we must compute L * M * M * L' * G'
     * at first I compute   temp1 = L  * M
     * after that I compute temp2 = L' * G'
     * after that I compute temp3 = M  * temp2
     * finally I compute   output = temp1 * temp3
    */

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < r; j++) {
            for (int k = 0; k < r; k++) {
                temp1[i][j] += newL[i][k] * M[k][j];
            }
        }
    }


    for (int i = 0; i < r; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                temp2[i][j] += newL[k][i] * M2[k][j];
            }
        }
    }


    for (int i = 0; i < r; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                temp3[i][j] += M[i][k] * temp2[k][j];
            }
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                M2Inv[i][j] += temp1[i][k] * temp3[k][j];
            }
        }
    }
}

__device__ void pseudoinversel4(double M2[][4], double M2Inv[][4], double v[][4], double *rv1, double *w, double res1[][4] )
{
    int m = 4;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;

    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel5(double M2[][5], double M2Inv[][5], double v[][5], double *rv1, double *w, double res1[][5] )
{
    int m = 5;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}
__device__ void pseudoinversel6(double M2[][6], double M2Inv[][6], double v[][6], double *rv1, double *w, double res1[][6] )
{
    int m = 6;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel7(double M2[][7], double M2Inv[][7], double v[][7], double *rv1, double *w, double res1[][7] )
{
    int m = 7;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel8(double M2[][8], double M2Inv[][8], double v[][8], double *rv1, double *w, double res1[][8] )
{
    int m = 8;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel9(double M2[][9], double M2Inv[][9], double v[][9], double *rv1, double *w, double res1[][9] )
{
    int m = 9;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel10(double M2[][10], double M2Inv[][10], double v[][10], double *rv1, double *w, double res1[][10] )
{
    int m = 10;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel11(double M2[][11], double M2Inv[][11], double v[][11], double *rv1, double *w, double res1[][11] )
{
    int m = 11;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel12(double M2[][12], double M2Inv[][12], double v[][12], double *rv1, double *w, double res1[][12] )
{
    int m = 12;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel13(double M2[][13], double M2Inv[][13], double v[][13], double *rv1, double *w, double res1[][13] )
{
    int m = 13;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

__device__ void pseudoinversel14(double M2[][14], double M2Inv[][14], double v[][14], double *rv1, double *w, double res1[][14])
{
    int m = 14;
    int flag, its,i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    /* Householder reduction to bidiagonal form */
    for (i = 0; i < m; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs(M2[k][i]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    M2[k][i] = (M2[k][i]/scale);
                    s += (M2[k][i] * M2[k][i]);
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += (M2[k][i] * M2[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            M2[k][j] += (f * M2[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    M2[k][i] = (M2[k][i]*scale);
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != m - 1)
        {
            for (k = l; k < m; k++)
                scale += fabs(M2[i][k]);
            if (scale)
            {
                for (k = l; k < m; k++)
                {
                    M2[i][k] = (M2[i][k]/scale);
                    s += (M2[i][k] * M2[i][k]);
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < m; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                            s += (M2[j][k] * M2[i][k]);
                        for (k = l; k < m; k++)
                            M2[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < m; k++)
                    M2[i][k] = M2[i][k] * scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        if (i < m - 1)
        {
            if (g)
            {
                for (j = l; j < m; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[i][k] * v[k][j]);
                    for (k = l; k < m; k++)
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < m; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = m - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < m - 1)
            for (j = l; j < m; j++)
                M2[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != m - 1)
            {
                for (j = l; j < m; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += (M2[k][i] * M2[k][j]);
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < m; k++)
                        M2[k][j] += (f * M2[k][i]);
                }
            }
            for (j = i; j < m; j++)
                M2[j][i] = (M2[j][i] * g);
        }
        else
        {
            for (j = i; j < m; j++)
                M2[j][i] = 0.0;
        }
        ++M2[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = m - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = (y * c + z * s);
                            M2[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < m; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                printf("Not converged\n");
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < m; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = (y * c + z * s);
                    M2[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    //start compute inverse matrix

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
        }
    }

    for(int rowNumber = 0; rowNumber < m; rowNumber++)
    {
        for(int colNumber = 0; colNumber < m; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0;
            for(int thirdIndex = 0; thirdIndex < m; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] = M2Inv[rowNumber][colNumber]+ res1[rowNumber][thirdIndex]*M2[colNumber][thirdIndex];
            }
        }
    }
    
}

// general pseudoinverse
__device__ void pseudoinverse(double M2[][ML], double M2Inv[][ML], double v[][ML], double *rv1, double *w, double res1[][ML], int order)
{
    int flag, its, i, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;

    /* Householder reduction to bidiagonal form */
    for (i = 0; i < order; i++)
    {
        /* Left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < order)
        {
            for (k = i; k < order; k++)
                scale += fabs(M2[k][i]);
            if (scale != 0.0)
            {
                for (k = i; k < order; k++)
                {
                    M2[k][i] /= scale;
                    s += M2[k][i] * M2[k][i];
                }
                f = M2[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][i] = f - g;
                if (i != order - 1)
                {
                    for (j = l; j < order; j++)
                    {
                        s = 0.0;
                        for (k = i; k < order; k++)
                            s += M2[k][i] * M2[k][j];
                        f = s / h;
                        for (k = i; k < order; k++)
                            M2[k][j] += f * M2[k][i];
                    }
                }
                for (k = i; k < order; k++)
                    M2[k][i] *= scale;
            }
        }
        w[i] = scale * g;

        /* Right-hand reduction */
        g = s = scale = 0.0;
        if (i < order && i != order - 1)
        {
            for (k = l; k < order; k++)
                scale += fabs(M2[i][k]);
            if (scale != 0.0)
            {
                for (k = l; k < order; k++)
                {
                    M2[i][k] /= scale;
                    s += M2[i][k] * M2[i][k];
                }
                f = M2[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                M2[i][l] = f - g;
                for (k = l; k < order; k++)
                    rv1[k] = M2[i][k] / h;
                if (i != order - 1)
                {
                    for (j = l; j < order; j++)
                    {
                        s = 0.0;
                        for (k = l; k < order; k++)
                            s += M2[j][k] * M2[i][k];
                        for (k = l; k < order; k++)
                            M2[j][k] += s * rv1[k];
                    }
                }
                for (k = l; k < order; k++)
                    M2[i][k] *= scale;
            }
        }
        anorm = MAX(anorm, fabs(w[i]) + fabs(rv1[i]));
    }

    /* Accumulate the right-hand transformation */
    for (i = order - 1; i >= 0; i--)
    {
        l = i + 1;
        if (i < order - 1)
        {
            if (g != 0.0)
            {
                for (j = l; j < order; j++)
                    v[j][i] = (M2[i][j] / M2[i][l]) / g;
                for (j = l; j < order; j++)
                {
                    s = 0.0;
                    for (k = l; k < order; k++)
                        s += M2[i][k] * v[k][j];
                    for (k = l; k < order; k++)
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j < order; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
    }

    /* Accumulate the left-hand transformation */
    for (i = order - 1; i >= 0; i--)
    {
        l = i + 1;
        g = w[i];
        if (i < order - 1)
            for (j = l; j < order; j++)
                M2[i][j] = 0.0;
        if (g != 0.0)
        {
            g = 1.0 / g;
            if (i != order - 1)
            {
                for (j = l; j < order; j++)
                {
                    s = 0.0;
                    for (k = l; k < order; k++)
                        s += M2[k][i] * M2[k][j];
                    f = (s / M2[i][i]) * g;
                    for (k = i; k < order; k++)
                        M2[k][j] += f * M2[k][i];
                }
            }
            for (j = i; j < order; j++)
                M2[j][i] *= g;
        }
        else
        {
            for (j = i; j < order; j++)
                M2[j][i] = 0.0;
        }
        M2[i][i] += 1.0;
    }

    /* Diagonalize the bidiagonal form */
    for (k = order - 1; k >= 0; k--)
    {   /* Loop over singular values */
        for (its = 0; its < 30; its++)
        {   /* Loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {   /* Test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (nm >= 0 && (fabs(w[nm]) + anorm == anorm))
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = -f * h;
                        for (j = 0; j < order; j++)
                        {
                            y = M2[j][nm];
                            z = M2[j][i];
                            M2[j][nm] = y * c + z * s;
                            M2[j][i] = z * c - y * s;
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {   /* Convergence */
                if (z < 0.0)
                {   /* Make singular value nonnegative */
                    w[k] = -z;
                    for (j = 0; j < order; j++)
                        v[j][k] = -v[j][k];
                }
                break;
            }
            if (its >= 30)
            {
                printf("Not converged\n");
                break;
            }

            /* Shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z)*(x + z) + h * (y / (f + SIGN(g, f)) - h)) / x;

            /* Next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < order; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z != 0.0)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 0; jj < order; jj++)
                {
                    y = M2[jj][j];
                    z = M2[jj][i];
                    M2[jj][j] = y * c + z * s;
                    M2[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }

    /* Compute inverse matrix */
    // compute res1 = v * (1 / w)
    for (int rowNumber = 0; rowNumber < order; rowNumber++)
    {
        for (int colNumber = 0; colNumber < order; colNumber++)
        {
            if (w[colNumber] != 0.0)
                res1[rowNumber][colNumber] = v[rowNumber][colNumber] / w[colNumber];
            else
                res1[rowNumber][colNumber] = 0.0;
        }
    }

    // compute M2Inv = res1 * M2^T
    for (int rowNumber = 0; rowNumber < order; rowNumber++)
    {
        for (int colNumber = 0; colNumber < order; colNumber++)
        {
            M2Inv[rowNumber][colNumber] = 0.0;
            for (int thirdIndex = 0; thirdIndex < order; thirdIndex++)
            {
                M2Inv[rowNumber][colNumber] += res1[rowNumber][thirdIndex] * M2[colNumber][thirdIndex];
            }
        }
    }
}

__global__ void scan_compact(int* G_Compact, const int* G, const int n, int *nprime){
	const int row = by;
	const int section = (n + blockDim.x - 1) / blockDim.x;
	int thid = 0;
	int tmp = 0;
	int stepSize = 0; 
	extern __shared__ int G_shared[];
	// copy a row of data into shared memory
    for (int cnt = 0; cnt < section; cnt++){
		thid = tx + blockDim.x * cnt;
        if( thid  < n){
			G_shared[thid] = G[row * n + thid];
        }
	}

	__syncthreads();
	for (int sec = 0; sec < section; sec++){
		thid = tx + blockDim.x * sec;
		stepSize = ( (n - sec * blockDim.x) / blockDim.x) > 0 ? blockDim.x : (n - sec * blockDim.x);
		for (int step = 1; step < stepSize; step = step * 2){
			if(thid < n){
				if ( tx < step ) {
					tmp = G_shared[thid];
				} else if (tx >= step){
					tmp = G_shared[thid] + G_shared[thid - step];
				}
			}
			__syncthreads();
			if(thid < n){
				G_shared[thid] = tmp;
			}
			__syncthreads();
		}
		if ( thid == (blockDim.x * (sec + 1) - 1) && sec != (section - 1) ){	
			G_shared[thid + 1] = G_shared[thid + 1] + G_shared[thid];
		}
		__syncthreads();
	}
	// ===============> Compact <===============
	const int row_size = G_shared[n - 1];
	
	for (int sec = 0; sec < section; sec++){
		thid = tx + blockDim.x * sec;
		if( thid  < n && thid > 0){
			if (G_shared[thid] != G_shared[thid - 1]){
				G_Compact[row * n + G_shared[thid] - 1] = thid;
			}
			if (thid >= row_size && thid != n - 1){
				G_Compact[row * n + thid] = 0;
			}
			if (thid == n - 1){
				atomicMax(nprime, G_shared[n - 1]);
				G_Compact[row * n + n - 1] = G_shared[n - 1];
			}
		}
	}

	if (tx == 0 && G[row * n] == 1){
		G_Compact[row * n] = 0;
	}
}


__device__ void inverse(double M2[][3], double M2Inv[][3])
{
    double det =  M2[0][0] * (M2[2][2] * M2[1][1]) - M2[0][0] * (M2[2][1] * M2[1][2])
                - M2[1][0] * (M2[2][2] * M2[0][1]) + M2[1][0] * (M2[2][1] * M2[0][2])
                + M2[2][0] * (M2[1][2] * M2[0][1]) - M2[2][0] * (M2[1][1] * M2[0][2]);
    double tmp = 1.0 / det;
    M2Inv[0][0] = tmp * (M2[1][1] * M2[2][2] - M2[1][2] * M2[2][1]);
    M2Inv[0][1] = tmp * (M2[0][2] * M2[2][1] - M2[0][1] * M2[2][2]);
    M2Inv[0][2] = tmp * (M2[0][1] * M2[1][2] - M2[0][2] * M2[1][1]);

    M2Inv[1][0] = tmp * (M2[1][2] * M2[2][0] - M2[1][0] * M2[2][2]);
    M2Inv[1][1] = tmp * (M2[0][0] * M2[2][2] - M2[0][2] * M2[2][0]);
    M2Inv[1][2] = tmp * (M2[0][2] * M2[1][0] - M2[0][0] * M2[1][2]);

    M2Inv[2][0] = tmp * (M2[1][0] * M2[2][1] - M2[1][1] * M2[2][0]);
    M2Inv[2][1] = tmp * (M2[0][1] * M2[2][0] - M2[0][0] * M2[2][1]);
    M2Inv[2][2] = tmp * (M2[0][0] * M2[1][1] - M2[0][1] * M2[1][0]);
}


__device__ void BINOM(int n, int k, int *out)
{
    int P, N1, R;
    // between n - k and k, N1 should be Max(n-k, k) and P should be Min(n-k, k);
    N1 = k;
    P = n - k;
    if (N1 <= P){
        N1 = P;
        P = k;
    }
    if(P == 0){
        R = 1;
    }
    else if( P == 1){
        R = N1 + 1;
    }
    else{
        R = N1 + 1;
        for (int i = 2; i < (P + 1); i++){
            R = ( R * (N1 + i) ) / i;
        }
    }
    *out = R; 
}

__device__  void IthCombination(int out[], int N, int P, int L)
{
    //The out[p] can be calculated  using formula out[p] = out[p - 1] + L - K. note that out[p] is in 1-base indexing
    int P1 = P - 1;
    int R;
    int k = 0;
    for (int i = 0; i < P1; i++){
        out[i] = 0;
        if(i > 0){
            out[i] = out[i - 1];
        }
        while(k < L){
            out[i] = out[i] + 1;
            BINOM(N - out[i], P - (i + 1), &R);
            k = k + R;
        }
        k = k - R;
    }
    out[P1] = out[P1 - 1] + L - k;
}
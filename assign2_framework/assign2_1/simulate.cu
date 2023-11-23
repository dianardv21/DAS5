#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>

#include "simulate.hh"

using namespace std;


/* Utility function, use to do error checking for CUDA calls
 *
 * Use this function like this:
 *     checkCudaCall(<cuda_call>);
 *
 * For example:
 *     checkCudaCall(cudaMalloc((void **) &deviceRGB, imgS * sizeof(color_t)));
 * 
 * Special case to check the result of the last kernel invocation:
 *     kernel<<<...>>>(...);
 *     checkCudaCall(cudaGetLastError());
**/

static void check(cudaError_t result) {
    if (result != cudaSuccess) {
        cerr << "cuda error: " << cudaGetErrorString(result) << endl;
        exit(EXIT_FAILURE);
    }
}

// kernel computing each data point
__global__ void waveKernel(const long i_max, double *old, double *curr, double *next) {

    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < i_max) // if data is unevenly distributed, skip non-existing data
    next[i] = 2*curr[i] - old[i] + 0.15 * (curr[i-1] - (2*curr[i] - curr[i+1]));

    __syncthreads();

}

__constant__ double c = 0.15;

double *simulate(const long i_max, const long t_max, const long block_size,
                 double *old_array, double *current_array, double *next_array) {

    // init memory size of vectors
    int memSize = i_max * sizeof(double);
    // init devices
    double *deviceOld = NULL;
    double *deviceCurr = NULL;
    double *deviceNext = NULL;
    // allocate mem for devices
    check( cudaMalloc((void **) &deviceOld,  memSize) );
    check( cudaMalloc((void **) &deviceCurr, memSize) );
    check( cudaMalloc((void **) &deviceNext, memSize) );
    // copy host arrays to GPU devices
    check( cudaMemcpy(deviceOld,  old_array,     memSize, cudaMemcpyHostToDevice) );
    check( cudaMemcpy(deviceCurr, current_array, memSize, cudaMemcpyHostToDevice) );
    check( cudaMemcpy(deviceNext, next_array,    memSize, cudaMemcpyHostToDevice) );

    // account for uneven distribution of threads
    int mod = i_max % block_size;
    if (mod != 0) mod = 1;
    int grid_size = i_max/block_size + mod;

    // calculate wave function
    for (int t = 0; t < t_max; t++) {
        
        // calc wave function
        waveKernel<<<grid_size, block_size>>>(i_max, deviceOld, deviceCurr, deviceNext);
        // swap buffers
        

    }
    
    // retrieve result from device to CPU
    check( cudaMemcpy(next_array, deviceNext, memSize, cudaMemcpyDeviceToHost) );

    cudaFree(deviceOld);
    cudaFree(deviceCurr);
    cudaFree(deviceNext);

    
    return current_array;
}

// cudaMemcpyToSymbol() for initializing c constant
// what's __shared__ used for? _syncthreads()
// why use __constant__ if you can just plug it in
// nvcc -Xptxas="-v" to see how many data per thread fits in fast registers
// prevent race conditions with atomics: order of operation still undefined though, good for counting
// copying data costs time, divide in chunks and let copying and processing go simultaneously
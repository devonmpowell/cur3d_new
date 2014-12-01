#include <stdio.h>
#include <cuda_runtime.h>

void curas_fill(float* data, int nx, int ny, int nz); 


void curas_checkerr(cudaError_t err, char* msg);

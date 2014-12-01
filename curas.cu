#include "curas.ch"

__global__ void k_fill(float* data, int nx, int ny, int nz) {
	for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < nx*ny*nz; 
         i += blockDim.x * gridDim.x) 
      {

          data[i] = sinf(i/1024.);
      
	}
}

void curas_fill(float* data_h, int nx, int ny, int nz) {

	cudaError_t e = cudaSuccess;

	float* data_d;
	e = cudaMalloc((void**) &data_d, nx*ny*nz*sizeof(float));
	curas_checkerr(e, "allocate device memory");


    // Launch the Vector Add CUDA Kernel
	int threadsPerBlock = 256;
	int stridePerBlock = 256;

	/*int numSMs;*/
	/*cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);*/

	int blocksPerGrid = (nx*ny*nz + threadsPerBlock - 1) / (threadsPerBlock * stridePerBlock);
	/*int blocksPerGrid = 256*numSMs;*/
	/*int threadsPerBlock = (nx*ny*nz)/blocksPerGrid;*/



    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    k_fill<<<blocksPerGrid, threadsPerBlock>>>(data_d, nx, ny, nz);
    e = cudaGetLastError();
	curas_checkerr(e, "kernel call");

	e = cudaMemcpy(data_h, data_d, nx*ny*nz*sizeof(float), cudaMemcpyDeviceToHost);
	curas_checkerr(e, "copy to host");

	e = cudaFree(data_d);
	curas_checkerr(e, "free device memory");


	return;
}
  
void curas_checkerr(cudaError_t err, char* msg) {
	if (err != cudaSuccess) {
		printf("CUDA Error: %s at %s.\n", cudaGetErrorString(err), msg);
		exit(0);
	}
}

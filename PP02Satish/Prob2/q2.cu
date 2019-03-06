#include <stdio.h>
#include <assert.h>

///////////////////////////////////////////////////////////////////////////////
//////Coalesced Kernal Function ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
 __global__ void Ax1(float *d_A, float *d_x, float *d_b, int *d_rs, int *d_NumBlocks,int *d_NumThreads){

	int ThreadId = threadIdx.x;
	int BlockId = blockIdx.x;
	int i, cs;

	cs = ThreadId + BlockId*(*d_NumBlocks);
	for(i=0; i<(*d_rs); i++)
		d_b[i*(*d_rs)+cs] = d_A[i*(*d_rs)+cs]*d_x[cs];
}

///////////////////////////////////////////////////////////////////////////////
//////Non-Coalesced Kernal Function ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
 __global__ void Ax2(float *d_A, float *d_x, float *d_b, int *d_rs, int *d_NumBlocks,int *d_NumThreads){

	int ThreadId = threadIdx.x;
	int BlockId = blockIdx.x;
	int i, cs; 
	
	cs = ThreadId + BlockId*(*d_NumBlocks);
    	for(i=0; i<(*d_rs); i++)
		d_b[i+cs*(*d_rs)] = d_A[cs*(*d_rs)+i]*d_x[i]
}


int main( int argc, char** argv)
{	
	int NumRow = 12800, NumCol = 12800;
        float A[NumRow][NumCol], x[NumCol], b[NumRow];
    
	int NumThreads = 256, NumBlocks = 50;
	int rs,cs;
	
	int i,j; 	
	struct timespec ti, tf;
	float dt;	
	
	for(i=0;i<NumRow;i++)				// Assigning values for A and x
		for(j=0;j<NumCol;j++)
			A[i][j] = 0.1/i+j;
	for(j=0;j<NumCol;j++)
		x[j] = 0.1/j;


   	float *d_A, *d_x,*d_b;				// Device Variables
	int *d_NumThreads,*d_NumBlocks;
	int *d_rs,*d_cs;

	cudaMalloc((void**)&d_A, sizeof(float)*NumRow*NumCol);
	cudaMalloc((void**)&d_x, sizeof(float)*NumCol);
	cudaMalloc((void**)&d_b, sizeof(float)*NumRow*NumCol);
        cudaMalloc((void**)&d_cs, sizeof(int));
        cudaMalloc((void**)&d_rs, sizeof(int));
	cudaMalloc((void**)&d_NumThreads, sizeof(int));
        cudaMalloc((void**)&d_NumBlocks, sizeof(int));

	cudaMemcpy(d_A, A, sizeof(float)*NumRow*NumCol , cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, sizeof(float)*NumCol, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(float)*NumRow*NumCol, cudaMemcpyHostToDevice);
        cudaMemcpy(d_cs, &rs, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rs, &cs, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_NumThreads, &NumThreads, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_NumBlocks, &NumBlocks, sizeof(int), cudaMemcpyHostToDevice);
	cudaThreadSynchronize();

        clock_gettime(CLOCK_MONOTONIC,&ti);
	Ax1<<<50,256>>>( d_A,d_x,d_b, d_rs,d_NumBlocks, d_NumThreads);	// Coalesced + Wraped
	cudaThreadSynchronize();
        clock_gettime(CLOCK_MONOTONIC, &tf);
        dt = (float)((tf.tv_sec-ti.tv_sec)*1e9 + tf.tv_nsec - ti.tv_nsec);
	printf("Coaleasced + Wraped: %lf \n", dt);
	
        clock_gettime(CLOCK_MONOTONIC,&ti);
	Ax1<<<128,100>>>( d_A,d_x,d_b, d_rs,d_NumBlocks, d_NumThreads);	// Coalesced + Non-Wraped
	cudaThreadSynchronize();
        clock_gettime(CLOCK_MONOTONIC, &tf);
        dt = (float)((tf.tv_sec-ti.tv_sec)*1e9 + tf.tv_nsec - ti.tv_nsec);
	printf("Coa + Non-Wraped: %lf \n", dt);

	clock_gettime(CLOCK_MONOTONIC,&ti);
	Ax2<<<50,256>>>( d_A,d_x,d_b, d_rs,d_NumBlocks, d_NumThreads);	// Non-Coalesced + Wraped
	cudaThreadSynchronize();
        clock_gettime(CLOCK_MONOTONIC, &tf);
        dt = (float)((tf.tv_sec-ti.tv_sec)*1e9 + tf.tv_nsec - ti.tv_nsec);
	printf("Non-Coa + Wraped: %lf \n", dt);

        clock_gettime(CLOCK_MONOTONIC,&ti);
	Ax2<<<128,100>>>( d_A,d_x,d_b, d_rs,d_NumBlocks, d_NumThreads);	// Non-Coalesced + Non-Wraped
	cudaThreadSynchronize();
        clock_gettime(CLOCK_MONOTONIC, &tf);
        dt = (float)((tf.tv_sec-ti.tv_sec)*1e9 + tf.tv_nsec - ti.tv_nsec);
	printf("Non-Coa + Non-Wraped: %lf \n", dt);

	cudaMemcpy(b, d_b, sizeof(float)*NumRow*NumCol, cudaMemcpyDeviceToHost);    
        cudaThreadSynchronize();
    
	return 0;
}

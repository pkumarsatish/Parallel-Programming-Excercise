#include <stdio.h>

__global__ void csr_mv(int m, double *AVal, int *ARowP, int *AColI, double *x, double *b)
{
	int i = ((blockDim.x * blockIdx.x) + threadIdx.x), j;
	if (i < m) {
		for (j = ARowP[i]; j < ARowP[i+1]; j++)
		{
			b[i] += AVal[j]*x[AColI[j]];
		}}
}

int main()
{
	// Define the Sparse Matrix
	int n0 = 80000;   // No. of no. zero enteries
	int m = 8000;     // m rows n collumns for A
	int n = m;
	int ARowP[m+1];
	int AColI[n0];
	double AVal[n0];

	double x[n];	// b=Ax
	double b[m];
	int thrd = 1000;
	int bloc = m/1000;

	int i,j;
//	Value Assignment
	
	for (i=0;i<n0;i++){
		AVal[i] = 1+i;
		AColI[i] = rand() % m;
        }

	for (i=0;i<(m+1);i++)
		ARowP[i] = (i-1)*10;

	for (i=0;i<n;i++)
		x[i] = 1;

	for (i=0;i<m;i++)
		b[i] = 0;


//  Data comminication

	double *d_AVal, *d_x, *d_b;
	int *d_AColI, *d_ARowP;

	cudaMalloc(&d_AVal, n0*sizeof(double));
	cudaMalloc(&d_x, n*sizeof(double));
	cudaMalloc(&d_b, m*sizeof(double));
	cudaMalloc(&d_ARowP, (m+1)*sizeof(int));
	cudaMalloc(&d_AColI, n0*sizeof(int));

	cudaMemcpy(d_AVal, AVal, n0*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x, n*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, m*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ARowP, ARowP, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_AColI, AColI, n0*sizeof(int), cudaMemcpyHostToDevice);

	csr_mv <<< bloc, thrd >>> ( m, d_AVal, d_ARowP, d_AColI, d_x, d_b);

	cudaMemcpy(b, d_b, m*sizeof(double), cudaMemcpyDeviceToHost);

	for (i = 0; i < 100; i++)
		printf("b[%d] = %lf\n", i, b[i]);

	return 0;
}

#include <stdio.h>
#include "mpi.h"

void UpdateRechable (int *RowP, int n, int *ColI, int nnz, int *Deg, int *r, int MinIdx, int *visited, int *Perm){

	int i,j,k;
	//printf("%d %d",visited[0],visited[1]);
	// Find rechable set r (Neighbour) for elemenated node MinIdx
	for (j=RowP[MinIdx]; j<RowP[MinIdx+1]; j++){
	if (visited[Perm[ColI[j]]] == 0){
		if(Deg[Perm[ColI[j]]] == n+1){
			visited[Perm[ColI[j]]] = 1;
			//printf("Visited: %d %d %d\n",visited[ColI[j]],ColI[j], itr);	
			UpdateRechable(RowP, n, ColI, nnz, Deg, r, ColI[j], visited, Perm);		 // enode
		}
		else{
			//printf("NotVisited: %d %d %d\n",visited[ColI[j]],ColI[j], itr);
			r[Perm[ColI[j]]] = 1;
		}
	}
	}
}

void UpdateDegree (int *RowP, int n, int *ColI, int nnz, int *Deg, int *ds, int Idx, int *visited, int *Perm){

	int i,j,k;

	// Find each Neighbour od Idx node
	for (j=RowP[Idx]; j<RowP[Idx+1]; j++){
	if (visited[Perm[ColI[j]]] == 0){
		if (Deg[Perm[ColI[j]]] == n+1){
			visited[Perm[ColI[j]]] = 1;
			UpdateDegree(RowP, n, ColI, nnz, Deg, ds, ColI[j], visited, Perm);		 // enode
		}
		else{
			ds[Perm[ColI[j]]] = 1;
			//printf("%d \n", ColI[j]);
		}
	}}
}

int No_fills_csr(int *RowP, int n, int *ColI, int nnz, int *Perm){

	int Deg[n], PreDeg[n]; // Deg = n+1: enode
	int r[n], ds[n]; // Rechable set (NB), r[i] = 1 => ith node is in r, else 0
	int visited[n];
	int nF;

	int i, j, k;

// Initialize Degree
	for (i=0; i<n; i++)	Deg[i] = RowP[i+1]-(RowP[i]+1);
	for (i=0; i<n; i++)     PreDeg[i] = Deg[i];
	nF = 0;

	for (i=0; i<n; i++){
		// Make Perm[i] as eNode
		Deg[Perm[i]] = n+1;
		
		// Rechable Nodes for node i
		for (j=0; j<n; j++)	r[j]=0;
		for (j=0; j<n; j++)	visited[j]=0;
		visited[Perm[i]] = 1;
		UpdateRechable (RowP, n, ColI, nnz, Deg, r, Perm[i], visited, Perm);
	    
		// Update Degree for i as a minidx
		for (j=0; j<n; j++){
			if(r[j] == 1){		// For all reachable of node MinIdx(i) as j
				Deg[j] = 0;
				for (k=0; k<n; k++)	ds[k]=0;
				for (k=0; k<n; k++)	visited[k]=0;
				visited[j] = 1;	
				UpdateDegree(RowP, n, ColI, nnz, Deg, ds ,j, visited, Perm);
				for (k=0; k<n; k++)	Deg[j] = Deg[j] + ds[k];
				nF = nF+(Deg[j]-PreDeg[j])+1;
				//printf("%d %d: %d %d\n", Deg[j], PreDeg[j], nF, j);
				PreDeg[j] = Deg[j];
			}
		}
	}
	nF = nF/2;
	return(nF);
}

void min_deg(int *RowP, int n, int *ColI, int nnz, int *Perm){

	int Deg[n], PreDeg[n]; // Deg = n+1: enode
	int MinIdx, MinDeg;
	int r[n], ds[n]; // Rechable set (NB), r[i] = 1 => ith node is in r, else 0
	int visited[n];
	//int enode_count;
	int IniPerm[n];

	int i, j, k;
// Initialize Degree
	for (i=0; i<n; i++)	Deg[i] = RowP[i+1]-(RowP[i]+1);
	for (i=0; i<n; i++)     PreDeg[i] = Deg[i];
	for (i=0; i<n; i++)	IniPerm[i] = Perm[i];	
// Find min Degree and Min Idx
	MinDeg = n;	
	for (i=0; i<n; i++){
		if (Deg[i] < MinDeg){
			MinDeg = Deg[i];
			MinIdx = i;
		}
	}

// Elemenating Min Deg Node
	Deg[MinIdx] = n+1; // enode
	Perm[0] = MinIdx;
	printf("Swap Nodes: %d %d\n", MinIdx, Perm[0]);	

// Do the same for each node except enode
	for (i=1; i<n; i++){
		// Update Rechable
		for (j=0; j<n; j++)	r[j]=0;
		for (j=0; j<n; j++)	visited[j]=0;
		visited[MinIdx] = 1;
		//printf("i = %d :",i);
		UpdateRechable (RowP, n, ColI, nnz, Deg, r, MinIdx, visited, IniPerm);
		//printf("Ho Gya");		
		//for (j=0; j<n; j++)	printf("%d ", r[j]);
		//Update Ranks for reachables
		for (j=0; j<n; j++){
			if(r[j] == 1){
				Deg[j] = 0;
				for (k=0; k<n; k++)	ds[k]=0;
				for (k=0; k<n; k++)	visited[k]=0;
				visited[j] = 1;	
				UpdateDegree(RowP, n, ColI, nnz, Deg, ds ,j, visited, IniPerm);
				for (k=0; k<n; k++)	Deg[j] = Deg[j] + ds[k];
				//printf("%d: %d\n", j, Deg[j]);
			}
		}
		// Find min Degree and Min Idx
		MinDeg = n;	
		for (j=0; j<n; j++){
			if (Deg[j] < MinDeg){
				MinDeg = Deg[j];
				MinIdx = j;
			}
		}
		// Elemenating Min Deg Node
	Deg[MinIdx] = n+1; // enode
	Perm[i] = MinIdx;
	//printf("Swap Node: %d %d\n", i, Perm[i]);
	}
	
// 
}


int main()
{

	int m, n = 4000, nnz = 15560;
	int ColI[nnz], RowP[n+1];
	int k = n/2;
	int i;
	char skip1[100];
	int temp[100];
	int Perm[n];	
	int nF;
	
	double startTime, endTime,totalTime;
  	startTime = MPI_Wtime();

	FILE *fp1, *fp2;
	fp1 = fopen("row.dat","r");
	if(fp1 == NULL){
		printf("Couldn't open file Row");
		return 1;
	}

	fp2 = fopen("colsr.dat","r");
        if(fp2 == NULL){
                printf("Couldn't open file Col");
                return 1;
        }
	
	for (i=0; i<=n; i++)	fscanf(fp1, " %d\n",&RowP[i]);
	for (i=0; i<nnz; i++){	
		fscanf(fp2, "%d\n",&ColI[i]);
		ColI[i] = ColI[i] -0;
	}
	for (i=0; i<n; i++) Perm[i]=i;
	int V[k], A[k], B[k];
	nF = No_fills_csr(RowP, n, ColI, nnz, Perm);
	printf ("No. of fills in Original Graph: %d\n", nF);

	min_deg(RowP, n, ColI, nnz, Perm);
	printf("MIN DEGREE");

	nF = No_fills_csr(RowP, n, ColI, nnz, Perm);
	printf ("No. of fills After Min-Degree Reordering: %d\n", nF);

	endTime = MPI_Wtime();
  	printf("total time %1.2f\n", endTime-startTime);
	return 0;
}

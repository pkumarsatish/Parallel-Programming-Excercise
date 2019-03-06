#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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


void KLND (int *RowP, int n, int *ColI, int nnz, int *Perm, int *IniPartition){
	
	int Partition[n];	// Partion - FOR NODE i: if a - (-1) or if b (+1)
	int Belonging[nnz];	// bELONGING - +1 if internal, -1 if external
	int Ecount, Icount, D[n];
	int k = n/2;
	int gain, max_gain, com_gain[k][3], max_com_gain_idx, max_com_gain = 1;
	int MarkCount = 0, Flag[n];
	int No_Sep=0;
	int i, j, count, i1, i2, i3=0, i4;


while (max_com_gain > 0){
	i3 = i3 + 1;
	for (i=0; i<n; i++)	Partition[i] = IniPartition[i];
	for (i=0; i<n; i++) Flag[i] = 0; // Not flaged, 1 - flaged
	max_gain = -(n+10);
	max_com_gain_idx = -1;
	max_com_gain = 0;
/////////////////////////////////////////////////////////////////////////////
	for (i=0; i<n; i++){
		Ecount = 0;
		Icount = 0;
		for (j=RowP[i]; j<RowP[i+1]; j++){
			if (Partition[ColI[j]] == Partition[i]){
				Belonging[j] = 1;
				Icount = Icount + 1;
			}
			else{
				Belonging[j] = -1;
				Ecount = Ecount + 1;
			}
		}
		D[i] = Ecount - Icount;
		//printf("Count: %d %d\n", Ecount, Icount);
	}
	//printf("Gain: %d %d %d\n", D[0], D[1], D[2]);

// Gain Maximization
	for (i=0; i<n; i++){
		for (j=i+1; j<n; j++){			
			if (Partition[i] != Partition[j]){
				gain = D[i] + D[j];
				for (i1=RowP[i]; i1<RowP[i+1]; i1++){
					if (ColI[i1] == j){		// (i,j) is non zero
						gain = gain - 2;
					}
				}
				if (gain > max_gain){
					max_gain = gain;
					com_gain[0][1] = i;
					com_gain[0][2] = j;
				}
			}
		}
	}
	com_gain[0][0] = max_gain;
	MarkCount = 1;
	//printf("Gain: %d %d %d\n", com_gain[0][0], com_gain[0][1], com_gain[0][2]);

///////////////////////////////////////////////////////////////////////////////////
// Do until all pairs (a,b) are marked

	for (i2=1; i2<k; i2++){
		// Update the partition (vertual swaping) and Belongings and D(n)
		Partition[com_gain[i2-1][1]] = -1*Partition[com_gain[i2-1][1]];
		Partition[com_gain[i2-1][2]] = -1*Partition[com_gain[i2-1][2]];
		Flag[com_gain[i2-1][1]] = 1;
		Flag[com_gain[i2-1][2]] = 1;
		
		for (i=0; i<n; i++){
			if(Flag[i] == 0){
			Ecount = 0;
			Icount = 0;
			for (j=RowP[i]; j<RowP[i+1]; j++){
				if (Partition[ColI[j]] == Partition[i]){
					Belonging[j] = 1;
					Icount = Icount + 1;
				}
				else{
					Belonging[j] = -1;
					Ecount = Ecount + 1;
				}
			}
			D[i] = Ecount - Icount;
			}
		}

		// Gain Maximization
		max_gain = -(n+10);
		for (i=0; i<n; i++){
			for (j=i+1; j<n; j++){
				if ((Flag[i] == 0) && (Flag[j] == 0)){			
				if (Partition[i] != Partition[j]){
					gain = D[i] + D[j];
					for (i1=RowP[i]; i1<RowP[i+1]; i1++){
						if (ColI[i1] == j){		// (i,j) is non zero
							gain = gain - 2;
						}
					}
					if (gain > max_gain){
						max_gain = gain;
						com_gain[i2][1] = i;
						com_gain[i2][2] = j;
					}
				}
				}
			}
		}

	com_gain[i2][0] = max_gain;
	MarkCount = MarkCount+1;
	//printf("Gain: %d %d %d\n", com_gain[i2][0], com_gain[i2][1], com_gain[i2][2]);
	}
////////////////////////////////////////////////////////////////////////////////////////////
	max_com_gain = com_gain[0][0];	
	for (i=1; i<k; i++){
		if(com_gain[i][0] > max_com_gain){
			max_com_gain = com_gain[i][0];
			max_com_gain_idx = i;
		}
	}
	//printf("%d ", max_com_gain_idx);

	if (max_com_gain > 0){
	for (i=0; i<=max_com_gain_idx; i++){
		IniPartition[com_gain[i][1]] = -1*IniPartition[com_gain[i][1]];
		IniPartition[com_gain[i][2]] = -1*IniPartition[com_gain[i][2]];
	}
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i=0; i<4; i++){
	//	printf("%d %d %d\n", IniPartition[i], IniPartition[i+4], IniPartition[i+8]);
	//}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//printf("%d: %d\n", i3, max_com_gain);
	if (i3 == 10) break;
	}
///////////////////////////////////////////
	i3 = 1;
	for (i=0; i<n; i++){
		if (IniPartition[i] == -1){
			for (j=RowP[i]; j<RowP[i+1]; j++){
				if (IniPartition[ColI[j]] == 1){
					if (i3 == -1)	IniPartition[i] = 0;
					if (i3 == 1)	IniPartition[ColI[j]] = 0;
					i3 = i3*-1;
					No_Sep = No_Sep + 1;
				}
			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i4=0; i4<4; i4++){
	//printf("%d %d %d\n", IniPartition[i4], IniPartition[i4+4], IniPartition[i4+8]);
	//}	
/////////////////////////////////////////////////////////////////////////////
// Re-Ordering
	i=0;
	j=k-(No_Sep/2);
	i3=n-No_Sep;

	for(i4=0; i4<n; i4++){
		if (IniPartition[i4] == -1){
			Perm[i4] = i;
			i=i+1;
		}
		else if (IniPartition[i4] == 1){
			Perm[i4] = j;
			j=j+1;
		}
		else if (IniPartition[i4] == 0){
			Perm[i4] = i3;
			i3=i3+1;
		}
	}
	printf("%d %d %d", i, j, i3);
////////////////////////////////////////////////////	
}
		

int main(int argc,char ** argv){
	int myid, nproc;
	MPI_Status status;
	int n=12, nnz=46;
	int ColI[nnz], RowP[n+1];
	int k = n/2;
	int level;
	int klfirstid, klnproc, klnodes;
	int Perm[n], myabv[n];
	int i, j;
	int ABV_count = 0, nF=0;

	double startTime, endTime,totalTime;
  	startTime = MPI_Wtime();

  	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

/////////////////////////////////////////////////////////
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
	
	for (i=0; i<=n; i++)	fscanf(fp1, "%d\n",&RowP[i]);
	for (i=0; i<nnz; i++){	
		fscanf(fp2, "%d\n",&ColI[i]);
		ColI[i] = ColI[i] -0;
	}

/////////////////////////////////////////////////////////////
	for (j=0; j<n; j++)	Perm[j] = j;
	for (j=0; j<k; j++)	myabv[j] = -1;
	for (j=k; j<n; j++)	myabv[j] = 1;	// -1 for rest	


/*	if(i==0)
		for (j=0; j<n; j++)	Perm[j] = j;
		for (j=0; j<k; j++)	myabv[j] = 1;
		for (j=k; j<n; j++)	myabv[j] = 2;	// -1 for rest	
*/

/*	for (i=0; i<6; i++){	// For Each Level Till Each process does thier respective parallel work
		level = pow(2.0,i);
		klnproc = nproc/(i+1);
		klnodes = n/klnproc;
		ABV_count = 0;
		//for (j=0; j<level; j++){	// @ ith level - 2^i nodes
		klfirstid = myid/klnproc;
		klfirstid = klfirstid*klnproc;
		//if (myid >= klfirstid && myid < klfirstid+klnproc){

		//}

		// Send by klfirstid, Perm ABV to 0 proc
		if (myid == klfirstid && myid != 0){
			//Perm_ABV_Add[]--> Recieve From klfirstid to 0;
		}
		if (myid == 0){
			for (j=1; j<level; j++) {
				//Recieve ABV From j*klnproc; @ Perm_ABV_Add[3*j]
			}
			//MPI_Barrier(MPI_COMM_WORLD);

			// Send G in Perm Mat to each of the proc 
		}
		
		//ABV_count = ABV_count + 3;
	}
*/
	//nF = No_fills_csr(RowP, n, ColI, nnz, Perm);
	//printf ("No. of fills in Original Graph: %d\n", nF);

	KLND (RowP, n, ColI, nnz, Perm, myabv);	
	nF = No_fills_csr(RowP, n, ColI, nnz, Perm);
	printf ("No. of fills After Min-Degree Reordering: %d\n", nF);

	endTime = MPI_Wtime();
  	printf("%d total time %1.2f\n",myid, endTime-startTime);
	MPI_Finalize();
	return 0;
}

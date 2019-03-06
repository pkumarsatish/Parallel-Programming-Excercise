#include <stdio.h>
#include "mpi.h"

int BinarySearch(int array[], int n, int search)
{
	int a,b,c, first, last, middle;
 	first = 0;
	last = n - 1;
	middle = (first+last)/2;
 
	while( first <= last ){
		if ( array[middle] < search )	first = middle + 1;    
		else if ( array[middle] == search ) {
			c = middle;
			break;
		}
		else	last = middle - 1;
 	
		middle = (first + last)/2;
	}

	if ( first > last ){
		a = abs(array[last]-search);
		b = abs(array[first]-search);		
		if ( a > b) c = first;
		else c = last;
   }
 
   return c;   
}


int main(int argc, char **argv)
{
	int myid, nproc;

	const int NumVertices = 50000, NumEdg = 200000;
	int RowP[NumVertices+1], ColI[NumEdg];
	double PageRank[50000] = { [ 0 ... 49999 ] = (1.0/NumVertices) };
	double Local_Sum[NumVertices],Global_Sum[NumVertices]; // = { [ 0 ... 49999 ] = 0.0 };	

	char skip1[100], skip2;
	int skip3;
	int i,j=1,k, temp1, temp2=1, NumNb, MaxItr = 10;
	double c1 = 0.15/NumVertices;

	double startTime, endTime,totalTime;
	startTime = MPI_Wtime();

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int rowi[NumVertices], rowf[NumVertices];
//////////////////////////////////////////////////////////////////////////////////
////Reading data from graph.dat and storing the adjecency matrix in CSR Formate///
/////////////////////////////////////////////////////////////////////////////////
	FILE *fp;
	fp = fopen("graph.dat", "r");
	if(fp == NULL){
		printf("Couldn't open file");
		return 1;
	}

	for(i=0;i<8;i++)	fgets(skip1, 100, fp);
	
	RowP[0] = 0;
	for(i=0;i<NumEdg;i++){
		temp1 = temp2;
		fscanf(fp, "%c %d %d %d\n",&skip2,&temp2,&ColI[i],&skip3);
		ColI[i]-=1;
		if (temp1 != temp2){
			for(j=temp1;j<temp2;j++)	RowP[j] = i;
		}
	}
	for(j=temp2;j<(NumVertices+1);j++)	RowP[j] = i;

	rewind(fp);
	fclose(fp);
//////////////////////////////////////////////////////////////////////////////////
////  Finding the first and last Row for each Process ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////
	//rowf = BinarySearch (RowP, (NumVertices+1), (myid+1)*(NumEdg/nproc));
	//rowi = BinarySearch (RowP, (NumVertices+1), (myid)*(NumEdg/nproc));
	for (i=0;i<nproc;i++){	
		rowf[i] = BinarySearch (RowP, (NumVertices+1), (i+1)*(NumEdg/nproc));
		rowi[i] = BinarySearch (RowP, (NumVertices+1), (i)*(NumEdg/nproc));
		rowi[i] = rowi[i] + 1;
	}
	
	//if (myid == 0)	rowi = 0;
	//if (myid == (nproc-1))	rowf = NumVertices-1;
	rowi[0] = 0;
	rowf[nproc-1] = 49999;
	//printf("%d : %d %d :%d\n", myid, rowi[myid], rowf[myid], (RowP[rowf[myid]]-RowP[rowi[myid]]));	


///////////////////////////////////////////////////////////////////////////////
/////////////// Sum Evaluation  //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
	for (j=0; j<NumVertices; j++){
		Local_Sum[j] = 0;
		Global_Sum[j] = 0;
	}

	for(k=0;k<30;k++){	
		for (j=rowi[myid]; j<=rowf[myid]; j++){
			NumNb = RowP[j+1] - RowP[j];
			for (i=RowP[j];i<RowP[j+1];i++){
				Local_Sum[ColI[i]] += PageRank[j]/NumNb;
			}
		}

	//MPI_Barrier( MPI_COMM_WORLD);
	for (i=0; i<nproc; i++){
		MPI_Reduce(&Local_Sum[rowi[i]], &Global_Sum[rowi[i]],(rowf[i]-rowi[i]+1),MPI_DOUBLE,MPI_SUM, i ,MPI_COMM_WORLD );
	}

	//MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
     //          MPI_Op op, int root, MPI_Comm comm)

	//MPI_Reduce(Local_Sum, Global_Sum, NumVertices, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//////////////////////////////////////////////////////////////////////////////////
///////////////// Page Rank Updation ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
		for (j=rowi[myid]; j<=rowf[myid]; j++){	
			PageRank[j] = (0.15/NumVertices + (0.85*Global_Sum[j]));
		}
		for (j=0; j<NumVertices; j++){	
			Local_Sum[j] = 0;
			Global_Sum[j] = 0;
		}	
	} // End of iteration

	//printf("%d : %lf %lf %lf %lf \n", myid, PageRank[0],PageRank[1], PageRank[49999], PageRank[2222]);
	endTime = MPI_Wtime();
  	printf("%d total time %1.2f\n",myid, endTime-startTime);
	MPI_Finalize();
	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char **argv)
{
  int myid, m;
  int n= 2000;                     // total grid points
  double u[8][8];
  memset(u, 0, sizeof(u));      // Solution matrix
  int nloc;                     // local no. of grid points 
  MPI_Status status;

  double startTime, endTime,totalTime;
  startTime = MPI_Wtime();

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &m);

  // nloc = n/m;              // no. of data per process

  int i,j,k;
  int maxitr = 2;

// /*********************************************************************/
// // /**Creating Cartesian 2D topology & Finding neibhours process rank****/
// // /*********************************************************************/


  int ndims = 2;
  int n1,n2,m1,m2;
  m1 = 4;
  m2 = 8;
  n1 = n/m1;
  n2 = n/m2;
  // int nloc1 = sqrt(nloc); // nloc1 x nloc1 grids per process 
  int dim[2]={m1,m2};
  int periods[2]={0,0};
  int reorder =0;
  MPI_Comm comm_cart;
  int mycartrank;
  int mycoords[2];
  int rightrank, leftrank, downrank, uprank;

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dim, periods, reorder, &comm_cart);
  MPI_Comm_rank(comm_cart, &mycartrank);
  MPI_Cart_coords(comm_cart, mycartrank, ndims, mycoords);

  MPI_Cart_shift(comm_cart, 0, -1, &downrank, &uprank);
  MPI_Cart_shift(comm_cart, 1, -1, &rightrank, &leftrank);

  double uloc[n1][n2];         // local sol matrix for each process
  double uold[n1][n2];         // temp for ghost rows
  double leftg[n1], rightg[n1],downg[n2],upg[n2];
  double leftb[n1], rightb[n1];

//  printf ("(%d,%d)\n",m1,nloc1);

/*********************************************************************/
/************** Initialization (0 at bounries)***********************/
/*********************************************************************/

  for(i=0;i<n1;i++){
    for(j=0;j<n2;j++){
      uloc[i][j] = (double)rand() / (double)RAND_MAX;
    }
  }

// printf("(%d,%d): %lf %lf %lf %lf\n",mycoords[0],mycoords[1],uloc[0][0],uloc[0][1],uloc[0][2],uloc[0][3]);

/****************************************************************************/
/***************************Communications & Computations********************/
/****************************************************************************/
  for (k=1;k<maxitr;k++){

    for(i=0;i<n1;i++){
      for(j=0;j<n2;j++){
        uold[i][j] = uloc[i][j];
      }
    }

    MPI_Sendrecv(&uloc[0][0], n2, MPI_DOUBLE, uprank, 1, upg, n2, MPI_DOUBLE, uprank, 1, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&uloc[n1-1][0], n2, MPI_DOUBLE, downrank, 1, downg, n2, MPI_DOUBLE, downrank, 1, MPI_COMM_WORLD, &status);

    for (i=0;i<n1;i++){
      leftb[i] = uloc[i][0];
      rightb[i] = uloc[i][n2-1];
    }

      MPI_Sendrecv(rightb, n1, MPI_DOUBLE, rightrank, 1, &rightg, n1, MPI_DOUBLE, rightrank, 1, MPI_COMM_WORLD, &status);
      MPI_Sendrecv(&leftb, n1, MPI_DOUBLE, leftrank, 1, &leftg, n1, MPI_DOUBLE, leftrank, 1, MPI_COMM_WORLD, &status);

  if(leftrank==MPI_PROC_NULL)
    for(i=0;i<n1;i++)
      leftg[i]=0;
  if(rightrank==MPI_PROC_NULL)
    for(i=0;i<n1;i++)
      rightg[i]=0;
  if(uprank==MPI_PROC_NULL)
	for(i=0;i<n2;i++)
	  upg[i]=0;
  if(downrank==MPI_PROC_NULL)
    for(i=0;i<n2;i++)
	  downg[i]=0;


    for(i=1;i<(n1-1);i++){
      for(j=1;j<(n2-1);j++){
        uloc[i][j] = (uold[i-1][j] + uold[i+1][j] + uold[i][j] + uold[i][j-1] + uold[i][j+1])/5;
      }
    }

    uloc[0][0] = (upg[0] + uold[1][0] + uold[0][0] + leftg[0] + uold[0][1])/5;
    uloc[0][n2-1] = (upg[n2-1] + uold[1][n2-1] + uold[0][n2-1] + uold[0][n2-2] + rightg[0])/5;
    uloc[n1-1][0] = (uold[n1-2][0] + downg[0] + uold[n1-1][0] + leftg[n1-1] + uold[n1-1][1])/5;
    uloc[n1-1][n2-1] = (uold[n1-2][n2-1] + downg[n2-1] + uold[n1-1][n2-1] + uold[n1-1][n2-2] + rightg[n1-1])/5;
// ...............

    for (i = 1; i<(n1-1); i++){
      uloc[i][0] = (uold[i-1][0] + uold[i+1][0] + uold[i][0] + leftg[i] + uold[i][1])/5;
      uloc[i][n1-1] = (uold[i-1][n1-1] + uold[i+1][n1-1] + uold[i][n1-1] + uold[i][n1-2] + rightg[i])/5;
      }
    for (i = 1; i<(n2-1); i++){
      uloc[0][i] = (upg[i] + uold[1][i] + uold[0][i] + uold[0][i-1] + uold[0][i+1])/5;
      uloc[n2-1][i] = (uold[n2-2][i] + downg[i] + uold[n2-1][i] + uold[n2-1][i-1] + uold[n2-1][i+1])/5;
    }

    MPI_Barrier( MPI_COMM_WORLD);
  }


/****************************************************************************/
/***************************Storing/Displaying Solution ********************/
/****************************************************************************/
//  printf("buufer: (%d,%d): %lf %lf %lf %lf\n",mycoords[0],mycoords[1],upg[0],upg[1],upg[2],upg[3]);
//  MPI_Barrier( MPI_COMM_WORLD);
//  printf("left buufer: (%d,%d): %lf %lf %lf %lf\n",mycoords[0],mycoords[1],leftg[0],leftg[1],leftg[2],leftg[3]);
//  MPI_Barrier( MPI_COMM_WORLD);
//  printf("After Itr:(%d,%d): %lf %lf %lf %lf\n",mycoords[0],mycoords[1],uloc[0][0],uloc[0][1],uloc[0][12],uloc[0][15]);

  endTime = MPI_Wtime();
  printf("(%d,%d) took total time %1.2f\n",n,m,endTime-startTime);
  MPI_Finalize();
  return 0;
}

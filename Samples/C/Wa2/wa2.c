#include <stdio.h>
#include "mpi.h"

int main(int argc, char* argv[]) {
  MPI_Status istatus;
  int  ierr;
  int  myid, nprocs;
  int  inlogp; 
  int  i, k;
  double  dsendbuf, drecvbuf;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  /* === 各プロセッサにおける総和演算したい値。ここではmyidとする。*/
  dsendbuf = myid;
  printf ("myid:%d, dsendbuf=%4.2lf \n", myid, dsendbuf);
  drecvbuf = 0.0;

  /* プロセッサー台数nprocsの2を底とする対数log_2(nprocs) */
  inlogp = 7;
  k = 1;
  for(i=0; i<=inlogp-1; i++) {
     if ( (myid&k)  == k ) {
       ierr = MPI_Recv(&drecvbuf, 1, MPI_DOUBLE, myid-k, i, MPI_COMM_WORLD, &istatus);
       dsendbuf = dsendbuf + drecvbuf;
       k *= 2;
     } else {
       ierr = MPI_Send(&dsendbuf, 1, MPI_DOUBLE, myid+k, i, MPI_COMM_WORLD);
       break;
     }
  }
  if (myid == nprocs-1) printf ("Total = %4.2lf \n", dsendbuf);
  ierr = MPI_Finalize();
}



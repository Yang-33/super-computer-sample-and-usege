#include <stdio.h>
#include "mpi.h"

int main(int argc, char* argv[]) {
  MPI_Status istatus;
  int  ierr;
  int  myid, nprocs; 
  double  dsendbuf, drecvbuf;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  /* === $B3F%W%m%;%C%5$K$*$1$kAmOB1i;;$7$?$$CM!#$3$3$G$O(Bmyid$B$H$9$k!#(B*/
  dsendbuf = myid;
  printf ("myid:%d, dsendbuf=%4.2lf \n", myid, dsendbuf);
  drecvbuf = 0.0;
  if (myid != 0) {
     ierr = MPI_Recv(&drecvbuf, 1, MPI_DOUBLE, myid-1, 0, MPI_COMM_WORLD, &istatus);
  }  
  dsendbuf = dsendbuf + drecvbuf;
  if (myid != nprocs-1) {
    ierr = MPI_Send(&dsendbuf, 1, MPI_DOUBLE, myid+1, 0, MPI_COMM_WORLD);
  }
  if (myid == nprocs-1) printf ("Total = %4.2lf \n", dsendbuf);

  ierr = MPI_Finalize();
}



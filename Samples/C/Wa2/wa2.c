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

  /* === $B3F%W%m%;%C%5$K$*$1$kAmOB1i;;$7$?$$CM!#$3$3$G$O(Bmyid$B$H$9$k!#(B*/
  dsendbuf = myid;
  printf ("myid:%d, dsendbuf=%4.2lf \n", myid, dsendbuf);
  drecvbuf = 0.0;

  /* $B%W%m%;%C%5!<Bf?t(Bnprocs$B$N(B2$B$rDl$H$9$kBP?t(Blog_2(nprocs) */
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



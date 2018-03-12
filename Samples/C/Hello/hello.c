#include <stdio.h>
#include "mpi.h"

int  main(int argc, char* argv[]) {

     int    myid, numprocs;
     int    ierr, rc;
      
     ierr = MPI_Init(&argc, &argv);
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

     printf("Hello parallel world!  Myid:%d \n", myid);

     rc = MPI_Finalize();

}



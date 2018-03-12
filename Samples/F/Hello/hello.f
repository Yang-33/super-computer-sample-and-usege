       program main
       include "mpif.h"
       common /mpienv/myid,numprocs

       integer  myid, numprocs
       integer  ierr

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

       print *, "Hello parallel world!  Myid:", myid

       call MPI_FINALIZE(ierr)

       stop
       end


       program main
       include 'mpif.h'
       common /mpienv/myid,numprocs

       integer  istatus(MPI_STATUS_SIZE)
       integer  myid, numprocs
       integer  ierr
       double precision dsendbuf, drecvbuf

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

c       === 各プロセッサにおける総和演算したい値。ここではmyidとする。
       dsendbuf = myid


       print *, "myid:", myid, "dsendbuf=", dsendbuf

       drecvbuf = 0.0
       if (myid .ne. 0) then
         call MPI_RECV(drecvbuf, 1, MPI_DOUBLE_PRECISION, myid-1, 0, 
     &          MPI_COMM_WORLD, istatus, ierr)
       endif
       dsendbuf = dsendbuf + drecvbuf
       if (myid .ne. numprocs-1) then
         call MPI_SEND(dsendbuf, 1, MPI_DOUBLE_PRECISION, myid+1, 0, 
     &         MPI_COMM_WORLD, ierr)
       endif
       if (myid .eq. numprocs-1) then
         print *, "Total = ", dsendbuf
       endif

       call MPI_FINALIZE(ierr)

       stop
       end

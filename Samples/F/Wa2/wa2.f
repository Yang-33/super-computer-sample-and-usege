       program main
       include 'mpif.h'
       common /mpyienv/myid,numprocs
       integer  istatus(MPI_STATUS_SIZE)
       integer  myid, numprocs
       integer  ierr
       integer  ilogp 
       integer  i, k
       double precision dsendbuf, drecvbuf;

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

c      === $B3F%W%m%;%C%5$K$*$1$kAmOB1i;;$7$?$$CM!#$3$3$G$O(Bmyid$B$H$9$k!#(B
       dsendbuf = myid
       print *, "myid:", myid, "dsendbuf=", dsendbuf
       drecvbuf = 0.0

c      ===  $B%W%m%;%C%5!<Bf?t(Bnprocs$B$N(B2$B$rDl$H$9$kBP?t(Blog_2(nprocs) 
       ilogp = 7
       k = 1
       do i=0, ilogp-1
        if ( iand(myid, k) .eq. k ) then
           call MPI_RECV(drecvbuf, 1, MPI_DOUBLE_PRECISION, myid-k, i, 
     &            MPI_COMM_WORLD, istatus, ierr)
           dsendbuf = dsendbuf + drecvbuf
           k = k * 2
        else           
           call MPI_SEND(dsendbuf, 1, MPI_DOUBLE_PRECISION, myid+k, i, 
     &            MPI_COMM_WORLD, ierr)
           goto 10
        endif     
       enddo

 10    CONTINUE

       if (myid .eq. numprocs-1) then
         print *, "Total = ", dsendbuf
       endif

       call MPI_FINALIZE(ierr)

       stop
       end


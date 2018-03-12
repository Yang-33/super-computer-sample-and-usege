c**********************************************************************
c   pi.c - compute pi by integrating f(x) = 4/(1 + x**2)     
c     
c   Each node: 
c    1) receives the number of rectangles used in the approximation.
c    2) calculates the areas of it's rectangles.
c    3) Synchronizes for a global summation.
c   Node 0 prints the result.
c
c  Variables:
c
c    pi  the calculated result
c    n   number of points of integration.  
c    x           midpoint of each rectangle's interval
c    f           function to integrate
c    sum,pi      area of rectangles
c    tmp         temporary scratch space for global summation
c    i           do loop index
c****************************************************************************

       program main
       include 'mpif.h'
       common /mpienv/myid,numprocs

       integer  myid, numprocs
       integer  ierr

       double precision  PI25DT 
       parameter (PI25DT= 3.141592653589793238462643)
       double precision  mypi, pi, h, sum, x;
       integer n, i;

       double precision t0, t1, t2, t_w

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
   
       if ( myid .eq. 0 ) then
         print *, "Enter the number of intervals"
         read (*,*) n
         print *, "n=",n
       endif

       call MPI_BCAST(n, 1, MPI_INTEGER,0,MPI_COMM_WORLD, ierr);

c      === check for quit signal 
       if ( n .le. 0 )  then
         stop
       endif


       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t1 = MPI_WTIME(ierr)

c      === calculate the interval size 
       h = 1.0 / n

       sum  = 0.0;

       do i=myid+1, n, numprocs
         x = h * (i - 0.5)
         sum = sum + 4.0 / (1.0 + x*x)
       enddo
       mypi = h * sum

c      === collect all the partial sums
       call MPI_REDUCE(mypi, pi, 1, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, 0, MPI_COMM_WORLD, ierr)


       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t2 = MPI_WTIME(ierr)

       t0 =  t2 - t1
       call MPI_REDUCE(t0, t_w, 1, MPI_DOUBLE_PRECISION, 
     &         MPI_MAX, 0, MPI_COMM_WORLD, ierr)


c      ===  myid 0 prints the answer
       if (myid .eq. 0) then
         print *, "  pi is approximately:", pi, 
     &         "Error is:", dabs(pi-PI25DT)
         print *, "Time(sec.) = ", t_w
       endif

       call MPI_FINALIZE(ierr)

       stop
       end




      




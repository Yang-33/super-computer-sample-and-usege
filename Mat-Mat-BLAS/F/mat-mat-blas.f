       program main

       include 'mpif.h'
       include 'mat-mat-blas.inc'

       integer DEBUG
       parameter (DEBUG=1)
       double precision EPS
       parameter (EPS=1.0e-18)


       double precision  A(NN, NN)
       double precision  B(NN, NN)
       double precision  C(NN, NN)

       integer N
       double precision  t0, t1, t2
       double precision  d_mflops, dtemp
    
       integer  i, j
       integer  iflag

       N = NN

c      === matrix generation --------------------------
       if (DEBUG .eq. 1) then
          do j=1, N
            do i=1, N
              A(i, j) = 1.0;
              B(i, j) = 1.0;
              C(i, j) = 0.0;
            enddo
          enddo
       else 
          call RANDOM_SEED
          do j=1, N
            do i=1, N
              call RANDOM_NUMBER(dtemp)   
              A(i, j) = dtemp
              call RANDOM_NUMBER(dtemp)   
              B(i, j) = dtemp
              C(i, j) = 0.0
            enddo
          enddo
       endif
c      === end of matrix generation ------------------------

c      === Start of mat-mat routine ----------------------------
       t1 = MPI_Wtime()

       call MyMatMat(C, A, B, N)

       t2 = MPI_Wtime()

       t0 =  t2 - t1 
c      === End of mat-mat routine ---------------------------

       print *, "N  = ", N
       print *, "Mat-Mat time[sec.] = ", t0

       d_mflops = 2.0*dble(N)*dble(N)*dble(N)/t0
       d_mflops = d_mflops * 1.0e-6
       print *, "MFLOPS = ", d_mflops

       if (DEBUG .eq. 1) then
c      === Verification routine ----------------- 
         iflag = 0
         do j=1, N
           do i=1, N 
             if (dabs(C(i,j) - dble(N)) > EPS) then
               print *, " Error! in (", i, ",", j, ") th argument."
               iflag = 1
               goto 10
             endif 
           enddo
         enddo
 10      continue
c        -----------------------------------------

         if (iflag .eq. 0) then
           print *, " OK!"
         endif
      
       endif
c      -----------------------------------------

       stop
       end



       subroutine MyMatMat(C ,A, B, n)
       include 'mat-mat-blas.inc'
       double precision C(NN, NN)
       double precision A(NN, NN)
       double precision B(NN, NN)
       integer n

       integer  i, j, k
       double precision ALPHA,BETA

       do i=1, n
         do j=1, n
           do k=1, n 
             C(i, j) = C(i, j) + A(i, k) * B(k, j) 
           enddo
         enddo
       enddo

c      Please write BLAS call in here.
c      CALL DGEMM(....)


      return
      end

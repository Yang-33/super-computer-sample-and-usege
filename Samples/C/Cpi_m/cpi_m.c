/*
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
*/

#include <stdio.h>
#include <math.h>
#include <mpi.h>


int main(int argc, char* argv[]) {

     double PI25DT = 3.141592653589793238462643;
     double mypi, pi, h, sum, x;
     double t1, t2, t0, t_w;
     int    n, myid, numprocs, i, rc;
     int    ierr;
      
     ierr = MPI_Init(&argc, &argv);
     ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
     ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

     /*
      sizetype   = 1;
      sumtype    = 2; 
     */
      
    if ( myid == 0 ) {
        printf("Enter the number of intervals \n");
        scanf("%d",&n);
        printf("n=%d \n",n);
     }
      
      ierr = MPI_Bcast(&n, 1, MPI_INT,0,MPI_COMM_WORLD);


      ierr = MPI_Barrier(MPI_COMM_WORLD);
      t1 = MPI_Wtime();

/*                                 check for quit signal */

     if ( n <= 0 ) {
        exit(0);
     }

/*
                                  calculate the interval size */
      h = 1.0 / n;

      sum  = 0.0;

      for (i = myid+1; i<=n; i+= numprocs) {
         x = h * (i - 0.5);
         sum = sum + 4.0 / (1.0 + x*x);
      }
      mypi = h * sum;

/*
                                  collect all the partial sums */
       MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE,
                  MPI_SUM, 0, MPI_COMM_WORLD);

       ierr = MPI_Barrier(MPI_COMM_WORLD);
       t2 = MPI_Wtime();
       t0 =  t2 - t1; 
       ierr = MPI_Reduce(&t0, &t_w, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


/*                                 node 0 prints the answer. */
       if (myid == 0) {
         printf("  pi is approximately: %18.16lf Error is: %18.16lf \n", 
                 pi, fabs(pi-PI25DT));
         printf("  execution time = : %8.4lf  [sec.] \n", t_w);         
       }

       rc = MPI_Finalize();

     }



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <mpi.h>

#define  N      1000
#define  DEBUG  1
#define  EPS    1.0e-18

/* Please define here for the matrices. */
static double  A[N][N];
static double  B[N][N];
static double  C[N][N];

void MyMatMat(double [][], double [][], double [][], int); 


void main(int argc, char* argv[]) {

     double  t1, t2, t_w;
     double  dc_inv, d_mflops;

     int     i, j;      
     int     iflag, iflag_t;


     /* matrix generation --------------------------*/
     if (DEBUG == 1) {
       for(j=0; j<N; j++) {
         for(i=0; i<N; i++) {
           A[j][i] = 1.0;
           B[j][i] = 1.0;
           C[j][i] = 0.0;
         }
       }
     } else {
       srand(1);
       dc_inv = 1.0/(double)RAND_MAX;
       for(j=0; j<N; j++) {
         for(i=0; i<N; i++) {
           A[j][i] = rand()*dc_inv;
           B[j][i] = rand()*dc_inv;
           C[j][i] = 0.0;
         }
       }
     } /* end of matrix generation --------------------------*/

     /* Start of mat-mat routine ----------------------------*/
     t1 = MPI_Wtime();

     MyMatMat(C, A, B, N);

     t2 = MPI_Wtime();
     t_w =  t2 - t1; 
     /* End of mat-mat routine --------------------------- */

     printf("N  = %d \n",N);
     printf("Mat-Mat time  = %lf [sec.] \n",t_w);

     d_mflops = 2.0*(double)N*(double)N*(double)N/t_w;
     d_mflops = d_mflops * 1.0e-6;
     printf(" %lf [MFLOPS] \n", d_mflops);

     /* Verification routine ----------------- */
     if (DEBUG == 1) {
       iflag = 0;
       for(j=0; j<N; j++) { 
         for(i=0; i<N; i++) { 
           if (fabs(C[j][i] - (double)N) > EPS) {
             printf(" Error! in ( %d , %d ) th argument. \n",j, i);
             iflag = 1;
             break;
           } 
         }
       }
       if (iflag == 0) printf("OK! \n");
     }
     /* ------------------------------------- */
     printf("nenett\n");
     exit(0);

}

void MyMatMat(double C[N][N], double A[N][N], double B[N][N], int n) 
{
  int  i, j, k;
  double ALPHA, BETA;
  char TEX[1] = {'N'};

  //  for(i=0; i<n; i++) {
  // for(j=0; j<n; j++) {
  //   for(k=0; k<n; k++) {
  //     C[i][j] += A[i][k] * B[k][j]; 
  //   }
  //  }
  //}

  //  Please write BLAS call in here.
  //  dgemm_(....);
  ALPHA = BETA = 1;
  dgemm_(&TEX,&TEX ,&n ,&n ,&n ,&ALPHA ,A ,&n ,B ,&n ,&BETA ,C ,&n );
   
}


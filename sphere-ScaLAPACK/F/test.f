!C
!C****
!C****  tightly connected sphere elements
!C****
!C
      implicit REAL*8 (A-H,O-Z)
      include 'mpif.h'


      real(kind=8) :: VOL, AREA, QVOL, COND0, COND, RADImax
      real(kind=8) :: HCONV, T0, SURF, DEL, coef1, coef2

      real(kind=8), dimension(:)  , allocatable :: XC, YC, ZC
      real(kind=8), dimension(:,:), allocatable :: AMAT
      real(kind=8), dimension(:)  , allocatable :: RHS

      real(kind=8), dimension(:), allocatable :: XX

      real(kind=8) :: DFLOPS
      real(kind=8) :: tw, t1, t2

*     ---------------------- ScaLAPACK variables 
*     .. Parameters ..
      INTEGER          DLEN_, IA, JA, IB, JB, MM, NN, MB, NB, RSRC,
     $                 CSRC, MXLLDA, MXLLDB, NRHS, NBRHS, 
     $                 MXLOCR, MXLOCC, MXRHSC
      PARAMETER      ( DLEN_ = 1000, IA = 1, JA = 1, IB = 1, JB = 1,
     $                 MM = 1000, NN = 1000, MB = 32, NB = 32, RSRC = 0,
     $                 CSRC = 0, MXLLDA = 128, MXLLDB = 128, NRHS = 1,
     $                 NBRHS = 1, NOUT = 6, MXLOCR = 128, MXLOCC = 128,
     $                 MXRHSC = 1 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )

*     .. Local Scalars ..
      INTEGER            ICTXT, INFO, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   EPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ),
     $                   IPIV( MXLOCR+NB )
      DOUBLE PRECISION   A( MXLLDA, MXLOCC ),
     $                   B( MXLLDB, MXRHSC )

*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLAMCH, PDLANGE, PDDOT
      EXTERNAL           PDLAMCH, PDLANGE, PDDOT

*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   DESCINIT, PDGEMM, PDGESV, PDLACPY,
     $                   SL_INIT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Data statements ..
      DATA               NPROW / 8 / , NPCOL / 8 /
*     ---------------------- ScaLAPACK variables 


      integer, dimension(:), allocatable :: PIV
      integer :: NX, NY, NZ, N
      integer :: myid, nprocs

!C
!C +------+
!C | INIT |
!C +------+
!C===

      call MPI_Init( ierr )
      call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
      call MPI_Comm_rank( MPI_COMM_WORLD, myid, ierr )


      open (11, file= 'inp1', status='unknown')
        read (11,*) NX, NY, NZ
        read (11,*) DX, DY, DZ
        read (11,*) VOL, AREA, QVOL, COND0
        read (11,*) HCONV, T0, SURF
        read (11,*) RADImax
        read (11,*) coef1, coef2
      close (11)

      N= NX*NY*NZ
      
      allocate (XC(N), YC(N), ZC(N))

      icou= 0
      do k= 1, NZ
        do j= 1, NY
          do i= 1, NX
            icou= icou + 1
            XC(icou)= dfloat(i-1)*DX
            YC(icou)= dfloat(j-1)*DY
            ZC(icou)= dfloat(k-1)*DZ
          enddo
        enddo
      enddo
!C===      

!C
!C +--------+
!C | MATRIX |
!C +--------+
!C===
      allocate (AMAT(N,N), RHS(N))

      AMAT= 0.d0
      RHS = 0.d0
      do i= 1, N
        do j= 1, N
          if (j.ne.i) then
            DEL= dsqrt((XC(i)-XC(j))**2 + (YC(i)-YC(j))**2              &
     &                                  + (ZC(i)-ZC(j))**2)

            if (DEL.le.RADImax) then
              COND= COND0/(10.d0**dmin1(DEL,20.d0))
              coef= COND*AREA / DEL
              AMAT(i,j)= coef
              AMAT(i,i)= AMAT(i,i) - coef
            endif
          endif
        enddo
      enddo

      do i= 1, N
        DELQ= dsqrt(XC(i)**2 + YC(i)**2 + ZC(i)**2)
        RHS(i)= -QVOL*(coef1*VOL + coef2*DELQ*VOL)
      enddo
      
      i= 1
      do k= 1, NZ
      do j= 1, NY
        ic= (k-1)*NX*NY + (j-1)*NX + i
        AMAT(ic,ic)= -HCONV*SURF    + AMAT(ic,ic)
        RHS (ic   )= -HCONV*SURF*T0 + RHS (ic)
      enddo
      enddo
!C===

!C
!C +---------------+
!C | ScaLAPACK Call|
!C +---------------+
!C===

*     .. Executable Statements ..
*
*     INITIALIZE THE PROCESS GRID
*
      CALL SL_INIT( ICTXT, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     If I'm not in the process grid, go to the end of the program
*
      IF( MYROW.EQ.-1 )
     $   GO TO 10
*
*     DISTRIBUTE THE MATRIX ON THE PROCESS GRID
*     Initialize the array descriptors for the matrices A and B
*
      CALL DESCINIT( DESCA, MM, NN,   MB, NB,    RSRC, CSRC, ICTXT, 
     $               MXLLDA, INFO )
      CALL DESCINIT( DESCB, NN, NRHS, NB, NBRHS, RSRC, CSRC, ICTXT,
     $               MXLLDB, INFO) 

*
*     Generate matrices A and B and distribute to the process grid
*
      DO J = 1, N
        DO I = 1, N
          CALL PDELSET( A, I, J, DESCA, AMAT(I,J))
        ENDDO
        CALL PDELSET( B, J, 1,  DESCB, RHS(J))
      ENDDO

*
       
      t1 = MPI_WTIME()

*     CALL THE SCALAPACK ROUTINE
*     Solve the linear system A * X = B
*
c     ------------------------------------------------
c     Please write ScaLAPACK call in here.
c      CALL PDGESV( .......... )
*
      CALL PDGESV(NN,NRHS,A,IA,JA,DESCA,IPIV,B,IB,JB,DESCB,INFO)
      t2 = MPI_WTIME()
      tw = t2 - t1


!C********************************************************************


  900 continue
!C===


      call MPI_ALLREDUCE(tw, t1, 1, MPI_DOUBLE_PRECISION, 
     &             MPI_MAX, MPI_COMM_WORLD, ierr)   

      if ( (MYCOL .eq. 0).and.(MYROW .eq. 0) ) then

        PRINT *, "TIME[sec] =",t1
 
        DFLOPS = 2.0d0/3.0d0*DBLE(N)*DBLE(N)*DBLE(N)
        DFLOPS = DFLOPS + 7.0d0/2.0d0*DBLE(N)*DBLE(N)
        DFLOPS = DFLOPS + 4.0d0/3.0d0*DBLE(N)
        DFLOPS = DFLOPS / DBLE(TW) / 1.0E6
        PRINT *, "MFLOPS =",DFLOPS


        IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
           WRITE( NOUT, FMT = 9999 )
           WRITE( NOUT, FMT = 9998 )MM, NN, NB
           WRITE( NOUT, FMT = 9997 )NPROW*NPCOL, NPROW, NPCOL
           WRITE( NOUT, FMT = 9996 )INFO
        END IF

      endif


!C
!C +--------+
!C | OUTPUT |
!C +--------+
!C===
      N1= 1
      N3= 3

!C
!C-- MESH
      open (22, file='sphere.fld', status='unknown')
      write (22,'(a)')    '# AVS field file'
      write (22,'(a,i5)') 'ndim=', N3
      write (22,'(a,i5)') 'dim1=', NX
      write (22,'(a,i5)') 'dim2=', NY
      write (22,'(a,i5)') 'dim3=', NZ
      write (22,'(a,i5)') 'nspace=', N3
      write (22,'(a,i5)') 'veclen=', N1
      write (22,'(a,i5)') 'data= float'
      write (22,'(a,i5)') 'field= uniform'
      write (22,'(a,i5)') 'label= temperature'

      write (22,'(a,i5)') 'variable 1 file=./spheredata filetype=ascii'
      close (22)

!C
!C-- RESULTS
c      open (22, file='spheredata', status='unknown')
c      do i= 1, N
c        write (22,'(1pe16.6)') RHS(i)
c      enddo

c      if ( (MYCOL .eq. 0).and.(MYROW .eq.1) ) then
c
c        do i= 1, MXLLDB
c          write (*,'(i8,1pe16.6)') i,B(i,1)
c        enddo
c
c      endif
!C===

      allocate (XX(N))

      open (22, file='spheredata_ans', status='unknown')
      do i= 1, N
        read (22,'(1pe16.6)') XX(i)
      enddo


c     --- copy answer to XX
      do i=1, MXLLDB
        A(i,1) = 0.0d0
      enddo
      DO J = 1, N
        CALL PDELSET(A, J, 1, DESCB, XX(J))
      ENDDO
    

c     --- EPS = dot(B,XX)
      if (MYCOL .eq. 0) then
        EPS = 0.0d0
        do i=1, MXLLDB
          EPS = EPS + (A(i,1)-B(i,1))*(A(i,1)-B(i,1))
        enddo 
      endif

c     ---- Reduction of EPS
      CALL DGSUM2D(ICTXT,'C',' ' ,1,1,EPS,1,0,0)
 
      if ((MYCOL .eq. 0).and.(MYROW .eq. 0)) then
        EPS = DSQRT(EPS)
        PRINT *, " EPS =",EPS
      endif



      deallocate(XX)


 910  continue





 9999   FORMAT( / 'ScaLAPACK Example Program -- Sep, 2010' )
 9998   FORMAT( / 'Solving Ax=b where A is a ', I5, ' by ', I5,
     $      ' matrix with a block size of ', I3 )
 9997   FORMAT( 'Running on ', I3, ' processes, where the process grid',
     $      ' is ', I3, ' by ', I3 )
 9996   FORMAT( / 'INFO code returned by PDGESV = ', I3 )

*

*
*     RELEASE THE PROCESS GRID
*     Free the BLACS context
*
      CALL BLACS_GRIDEXIT( ICTXT )
 10    CONTINUE



*     Exit the BLACS
*
      CALL BLACS_EXIT( 0 )
*
c      call MPI_FINALIZE(MPI_COMM_WORLD,ierr)

*
      stop
      end   




      SUBROUTINE SL_INIT( ICTXT, NPROW, NPCOL )
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, NPCOL, NPROW
*     ..
*
*  Purpose
*  =======
*
*  SL_INIT initializes an NPROW x NPCOL process grid using a row-major
*  ordering  of  the  processes. This routine retrieves a default system
*  context  which  will  include all available processes. In addition it
*  spawns the processes if needed.
*
*  Arguments
*  =========
*
*  ICTXT   (global output) INTEGER
*          ICTXT specifies the BLACS context handle identifying the
*          created process grid.  The context itself is global.
*
*  NPROW   (global input) INTEGER
*          NPROW specifies the number of process rows in the grid
*          to be created.
*
*  NPCOL   (global input) INTEGER
*          NPCOL specifies the number of process columns in the grid
*          to be created.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IAM, NPROCS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
*     If machine needs additional set up, do it now
*
      IF( NPROCS.LT.1 ) THEN
         IF( IAM.EQ.0 )
     $      NPROCS = NPROW*NPCOL
         CALL BLACS_SETUP( IAM, NPROCS )
      END IF
*
*     Define process grid
*
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
*
      RETURN
*
*     End of SL_INIT
*
      END



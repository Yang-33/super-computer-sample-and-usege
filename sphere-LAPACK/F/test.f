!C
!C****
!C****  tightly connected sphere elements
!C****
!C
      include "mpif.h"

      implicit REAL*8 (A-H,O-Z)

      real(kind=8) :: VOL, AREA, QVOL, COND0, COND, RADImax
      real(kind=8) :: HCONV, T0, SURF, DEL, DEL0, coef1, coef2
      real(kind=8) :: EPS

      real(kind=8), dimension(:)  , allocatable :: XC, YC, ZC
      real(kind=8), dimension(:,:), allocatable :: AMAT
      real(kind=8), dimension(:)  , allocatable :: RHS

      real(kind=8), dimension(:), allocatable :: XX

      real(kind=8) :: DFLOPS
      real(kind=8) :: tw, t1, t2 

      integer, dimension(:), allocatable :: PIV
      integer :: INC, INFO

      integer :: NX, NY, NZ, N
      integer :: R, Z, P, Q, DD

!C
!C +------+
!C | INIT |
!C +------+
!C===
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
!C | LAPACK Call   |
!C +---------------+
!C===

      INC = 1
      ALLOCATE(PIV(N))

      t1 = MPI_Wtime()

c      -----------------------------------------------------
c      Please write LAPACK Call in here.
c      call DGESV(.....)

      t2 = MPI_Wtime()

      tw =  t2 - t1

      PRINT *, "TIME[sec] =",tw

      DFLOPS = 2.0d0/3.0d0*DBLE(N)*DBLE(N)*DBLE(N)
      DFLOPS = DFLOPS + 7.0d0/2.0d0*DBLE(N)*DBLE(N)
      DFLOPS = DFLOPS + 4.0d0/3.0d0*DBLE(N)
      DFLOPS = DFLOPS / DBLE(tw) / 1.0E6
      PRINT *, "MFLOPS =",DFLOPS


      DEALLOCATE(PIV)
      DEALLOCATE(AMAT)

!C********************************************************************


  900 continue
!C===

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
      open (22, file='spheredata', status='unknown')
      do i= 1, N
        write (22,'(1pe16.6)') RHS(i)
      enddo

      do i= N-NX+1, N
        write (*,'(i8,1pe16.6)') i,RHS(i)
      enddo
!C===

      allocate (XX(N))

      open (22, file='spheredata_ans', status='unknown')
      do i= 1, N
        read (22,'(1pe16.6)') XX(i)
      enddo

      EPS = 0.0d0
      do i=1, N
        EPS = EPS + (XX(i)-RHS(i))*(XX(i)-RHS(i))
      enddo
      EPS = DSQRT(EPS)
      PRINT *, "EPS =",EPS
 
      deallocate(XX)


 910  continue
      stop
      end   

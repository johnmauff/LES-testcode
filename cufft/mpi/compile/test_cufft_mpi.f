! ======================================================================
!   code to test xderiv, yderiv, and 2d FFT's using cuFFT -- no mpi
! ======================================================================
!   different from the LES code, this code is currently configured to 
!   transform horizontal slabs of data at once.  did this to retain 
!   consistency with the LES code where xderiv is called within a k-loop. 
!   can easily shift to transforming 3d tubes of data by adding an 
!   izs:ize dimension on x_in/x_out (and equivalent) and increasing 
!   batch to match.
! ======================================================================

      module pars

         integer, parameter :: i_fft = 2   ! == 2 -> cuFFT

         integer, parameter :: nnx = 32
         integer, parameter :: nny = 32
         integer, parameter :: nnz = 1

         integer, parameter :: ncpu_s = 8

         real, parameter :: pi2 = 8.*atan(1.0)
         real, parameter :: xl = pi2
         real, parameter :: yl = pi2

         ! ---- set some variables to retain calling structure

         integer :: izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,
     +              mxs, mxe, iss, ise, iys, iye, jys, jye

         integer :: myid, numprocs, ncpu_z, maxp
         logical :: l_root

         integer, allocatable, dimension(:) ::
     +              ix_s, ix_e, jx_s, jx_e,
     +              kx_s, kx_e, mx_s, mx_e,
     +              iy_s, iy_e, jy_s, jy_e,
     +              is_s, is_e, iz_s, iz_e

         integer :: ncx,ncy
         real, dimension(nnx) :: xkn, xk
         real, dimension(nny) :: ykn, yk

      end module pars

! ======================================================================

      module fields
         real, allocatable, dimension(:,:,:) :: a,b,ay
         real, allocatable, dimension(:,:)   :: ax
      end module fields

! ======================================================================

      module cufft_wrk

        ! ---- setup plans and arrays for using cuFFT
        !      to use legacy fortran calls

        use cufft

        integer :: pln_xf, pln_xb, 
     +             pln_yf, pln_yb, 
     +             pln_cf, pln_cb

        integer :: nx_c, ny_c

        integer, device :: d_is,d_ir

        real, device, allocatable :: xk_d(:), yk_d(:)

        real, device, allocatable :: x_in(:,:), x_out(:,:), 
     +                               y_in(:,:), y_out(:,:),
     +                               c_in(:,:,:), c_out(:,:,:)

        contains

          subroutine cufft_config(xk,yk)

          use pars, only : iys,iye,ixs,ixe,jxs,jxe
         
          real, intent(in) :: xk(nx_c), yk(ny_c)
          integer :: jj, batch

          ! ---- allocate and copy wavenumbers to device

          allocate(xk_d(nx_c), yk_d(ny_c))

          do i = 1,nx_c
             xk_d(i) = xk(i)
          enddo
          do j = 1,ny_c
             yk_d(j) = yk(j)
          enddo

          ! ---- build plans 
          !
          !  https://docs.nvidia.com/cuda/pdf/CUFFT_Library.pdf (pp 7-8)
          !
          !  for out-of-place R2C or C2R transforms, input 
          !  and output sizes match the logical size N and 
          !  the non-redundant size N/2+1, respectively.
          !
          !  for in-place R2C and C2R transforms, the input
          !  size must be padded to N/2+1 complex elements.

          allocate( x_in(nx_c,  iys:iye)) 
          allocate(x_out(nx_c+2,iys:iye))

          allocate( y_in(ny_c,  ixs:ixe))
          allocate(y_out(ny_c+2,ixs:ixe))

          jj = (jxe-jxs+1)/2
          allocate(c_in(2,ny_c,jj), c_out(2,ny_c,jj))

          ! ---- plans for x

          batch = iye-iys+1       ! --- number of 1D ffts to execute

          ierr = cufftPlan1D(pln_xf,nx_c,CUFFT_D2Z,batch)  ! R2C => single prec, D2Z => double
          ierr = cufftPlan1D(pln_xb,nx_c,CUFFT_Z2D,batch)  ! C2R => single prec, Z2D => double

          ! ---- plans for y

          batch = ixe-ixs+1       ! --- number of 1D ffts to execute

          ierr = cufftPlan1D(pln_yf,ny_c,CUFFT_D2Z,batch)
          ierr = cufftPlan1D(pln_yb,ny_c,CUFFT_Z2D,batch)

          ! ---- plans for complex in y

          batch = (jxe-jxs+1)/2   ! --- number of 1D ffts to execute

          ierr = cufftPlan1D(pln_cf,ny_c,CUFFT_Z2Z,batch)
          ierr = cufftPlan1D(pln_cb,ny_c,CUFFT_Z2Z,batch)

          return
          end subroutine cufft_config

          subroutine cufft_finalize

          ! ---- destroy the plans

          ierr = cufftDestroy(pln_xf)
          ierr = cufftDestroy(pln_xb)
          ierr = cufftDestroy(pln_yf)
          ierr = cufftDestroy(pln_yb)
          ierr = cufftDestroy(pln_cf)
          ierr = cufftDestroy(pln_cb)
    
          ! ---- deallocate vars
    
          deallocate(xk_d, yk_d)
          deallocate(x_in, x_out)
          deallocate(y_in, y_out)
          deallocate(c_in, c_out)
    
          return
          end subroutine cufft_finalize

      end module cufft_wrk

! ======================================================================

      module fftwk
         real, allocatable :: trigx(:,:), trigc(:)
      end module fftwk

! ======================================================================

      program test_cufft

      use pars
      use fields
      use cufft
      use cufft_wrk
      use fftwk
      include 'mpif.h'

      ! --- initialize mpi

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)

      call gridd
      call build_wavenumbers
      call mpi_range

      ! ---- allocate vars

      allocate( a(nnx+2,iys:iye,izs-1:ize+1) )
      allocate( b(nny,  jxs:jxe,izs-1:ize+1) )

      allocate( ax(nnx,iys:iye) )
      allocate( ay(nnx,iys:iye,izs:ize) )

      ! ---- initialize the ffts

      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2),
     +         trigc(4*nq_trig+15))


      if (i_fft == 2) then  ! setup cufft plans
         nx_c = nnx
         ny_c = nny
         call cufft_config(xk,yk)
      endif

      ! ---- perform tests

      jj = 13

      call test_xderiv(jj)
      call test_yderiv(jj)
      call test_fft2d(jj)

      ! ---- clean up

      call cufft_finalize
    
      deallocate(a)
      deallocate(b)
      deallocate(ax)
      deallocate(ay)
      deallocate(trigx)
      deallocate(trigc)

      call mpi_finalize(ierr)

      stop
      end program test_cufft

! ======================================================================

      subroutine test_fft2d(jj)

      ! ---- routine to test fft2d

      use pars
      use fields
      use cufft
      use cufft_wrk
      use fftwk

      integer, intent(in) :: jj

      real, dimension(nnx) :: ai

      ! ---- initialize array on host

      dx = xl / dble(nnx)

      do k = izs,ize
         do j = iys,iye
         do i = 1,nnx
            a(i,j,k) = sin(dble(i-1)*dx)
         end do
         end do
      end do

      ! ---- grab desired line from the initial array for printout

      if (jj >= iys .and. jj <= iye) then
         do i = 1,nnx
            ai(i)  = a(i,jj,1)
         enddo
      endif

      ! ---- forward transform

      call fft2d_mpi(a(1,iys,izs),b(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)

      ! ---- backward transform

      call fft2d_mpi(a(1,iys,izs),b(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)

      ! ---- check ouput array

      if (jj >= iys .and. jj <= iye) then
         write(*,*)
         write(*,*) 'fft2d:'

         do i = 1,nnx
            write(*,100) dble(i-1)*dx, ai(i), a(i,jj,1)
         enddo
 100     format(' x = ',f,' , send = ',f,' , recv = ',f)
      endif

      return
      end subroutine test_fft2d

! ======================================================================

      subroutine test_xderiv(jj)

      ! ---- routine to test xderiv

      use pars
      use fields
      use cufft
      use cufft_wrk
      use fftwk

      integer, intent(in) :: jj

      ! ---- initialize array on host

      dx = xl / dble(nnx)

      do k = izs,ize
         do j = iys,iye
         do i = 1,nnx
            a(i,j,k) = sin(dble(i-1)*dx)
            ax(i,j) = a(i,j,k)
         end do
         end do
         call xderivp(ax(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      end do

      if (jj >= iys .and. jj <= iye) then
         write(*,*)
         write(*,*) 'xderiv:'

         do i = 1,nnx
            write(*,100) dble(i-1)*dx,a(i,jj,1),ax(i,jj)
         enddo
 100     format(' x = ',f,' , a = ',f,' , ax = ',f)
      endif

      return
      end subroutine test_xderiv

! ======================================================================

      subroutine test_yderiv(jj)

      ! ---- routine to test yderiv

      use pars
      use fields
      use cufft
      use cufft_wrk
      use fftwk

      real ::  at(nny,ixs:ixe,izs:ize)
      real :: ayt(nny,ixs:ixe,izs:ize)

      integer, intent(in) :: jj

      dy = yl / dble(nny)

      do k = izs,ize
         do j = iys,iye
         do i = 1,nnx
             a(i,j,k) = sin(dble(j-1)*dy)
            ay(i,j,k) = a(i,j,k)
         end do
         end do
      enddo
      call yd_mpi(ay(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

      ! ---- for printout

      call xtoy_trans(a(1,iys,izs),at,nnx,nny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      call xtoy_trans(ay,ayt,nnx,nny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

      if (jj >= ixs .and. jj <= ixe) then
         write(*,*)
         write(*,*) 'yderiv:'

         i = jj
         do j = 1,nny
            write(*,100) dble(j-1)*dy,at(j,i,1),ayt(j,i,1)
         enddo
 100     format(' y = ',f,' , a = ',f,' , ay = ',f)
      endif

      return
      end subroutine test_yderiv

! ======================================================================

      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)

      ! ---- choose which fft to use for 2D fft

      use pars, only : i_fft

      if (i_fft == 2) then
         call fft2d_cuda(ax,at,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)
      endif

      return
      end subroutine fft2d_mpi

! ======================================================================

      subroutine fft2d_cuda(ax,at,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)

      ! ---- get 2d fft using cuFFT routines
      !
      !     isgn = -1 do forward transform, get coefficients
      !               incoming array is ax(nx+2,iys:iye,iz1:iz2)
      !               outgoing array is ax(nx+2,iys:iye,iz1:iz2)
      !
      !     isgn = -2 do forward transform, get coefficients
      !               incoming array is ax(nx+2,iys:iye,iz1:iz2)
      !               outgoing array is at(ny,jxs:jxe,iz1:iz2)
      !
      !     isgn =  1 do inverse transform, move to physical space
      !               incoming array is ax(nx+2,iys:iye,iz1:iz2)
      !               outgoing array is ax(nx+2,iys:iye,iz1:iz2)
      !
      !     isgn =  2 do inverse transform, move to physical space
      !               incoming array is at(ny,jxs:jxe,iz1:iz2)
      !               outgoing array is ax(nx+2,iys:iye,iz1:iz2)

      use cufft
      use cufft_wrk

      real :: ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2)

      integer, dimension(0:np-1) :: jx_s,jx_e,iy_s,iy_e
      integer :: ij

      nxp2 = nx + 2

      if (isgn < 0) then

         fn = 1./dble(nx*ny)

         ! ---- 1d fft in x over [iys,iye] for all z

         do k = iz1,iz2
            do j = iys,iye
               do i = 1,nx
                  x_in(i,j) = ax(i,j,k)*fn                ! note: currently shifting one x-y slab to device
               enddo
            enddo
            ierr = cufftExecD2Z(pln_xf,x_in,x_out)        ! perform forward R2C 1D ffts in x-direction for [iys,iye]
            do j = iys,iye
               do i = 1,nxp2
                  ax(i,j,k) = x_out(i,j)                  ! note: shifting data back to host 
               enddo
            enddo
         enddo
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

         ! ---- 1d fft in y over [jxs,jxe] for all z

         do k = iz1,iz2
            ij = 0
            do i = jxs,jxe,2
               ij = ij + 1
               do j = 1,ny
                  c_in(1,j,ij) = at(j,i,k)                 ! note: shifting data to device here
                  c_in(2,j,ij) = at(j,i+1,k)
               enddo
            enddo
            ierr = cufftExecZ2Z(pln_cf, c_in, c_out, CUFFT_FORWARD) ! perform forward C2C 1D ffts in y for [jxs,jxe]/2
            ij = 0
            do i = jxs,jxe,2
               ij = ij + 1
               do j = 1,ny
                  at(j,i,k)   = c_out(1,j,ij)              ! note: shifting data back to host here
                  at(j,i+1,k) = c_out(2,j,ij)
               enddo
            enddo
         enddo

         ! ---- decide whether to transpose back or leave as is

         if (isgn == -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif

      else

         ! ---- decide whether to first transpose or leave as is

         if (isgn == 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif

         ! ---- 1d fft in y over [jxs,jxe] for all z

         do k = iz1,iz2
            ij = 0
            do i = jxs,jxe,2
               ij = ij + 1
               do j = 1,ny
                  c_in(1,j,ij) = at(j,i,k)
                  c_in(2,j,ij) = at(j,i+1,k)
               enddo
            enddo
            ierr = cufftExecZ2Z(pln_cb, c_in, c_out, CUFFT_INVERSE)
            ij = 0
            do i = jxs,jxe,2
               ij = ij + 1
               do j = 1,ny
                  at(j,i,k)   = c_out(1,j,ij)
                  at(j,i+1,k) = c_out(2,j,ij)
               enddo
            enddo
         enddo
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

         ! ----  1d fft in x over [iys,iye] for all z

         do k = iz1,iz2
            do j = iys,iye
               do i = 1,nxp2
                  x_out(i,j) = ax(i,j,k)               ! note: shifting data to device
               enddo
            enddo
            ierr = cufftExecZ2D(pln_xb, x_out, x_out)  ! perform backward 1D C2R ffts in x-direction for entire yz-tube
            do j = iys,iye
               do i = 1,nx
                  ax(i,j,k) = x_out(i,j)               ! note: shifting data back to host
               enddo
            enddo
         enddo
      endif

      return
      end subroutine fft2d_cuda

! ====================================================================

      subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
 
      ! ---- transpose array  f(nx,iys:iye,iz1:iz2)
      !                  ---> g(ny,ixs:ixe,iz1:iz2)
      !      load balanced version, march 2017
 
      include 'mpif.h'
      integer istatus(mpi_status_size)
 
      real :: f(nx,iys:iye,iz1:iz2),
     +        g(ny,ixs:ixe,iz1:iz2)
      real :: ft(nx*(iye+1-iys)*(iz2-iz1+1)),
     +        gt(ny*(ixe+1-ixs)*(iz2-iz1+1))

      real, device :: d_ft(nx*(iye+1-iys)*(iz2-iz1+1)),
     +                d_gt(ny*(ixe+1-ixs)*(iz2-iz1+1))
      integer :: nsend,nrecv,is,ir,ireqs,ireqr
      integer, device :: d_nsend,d_nrecv,d_is,d_ir,d_ireqs,d_ireqr

      integer, intent(in), dimension(0:np-1) :: ix_s,ix_e,iy_s,iy_e

      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
 
      ! ---- loop over cpus on a slab for given myid
 
      iss = (myid/ncpu_s)*ncpu_s
 
      do i = 1,ncpu_s
         is    = mod(myid - iss + i,ncpu_s) + iss
         ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
         nsend = (ix_e(is) - ix_s(is) + 1)*jk
         nrecv = (iy_e(ir) - iy_s(ir) + 1)*ik

         d_is = is
         d_ir = ir
         d_nsend = nsend
         d_nrecv = nrecv

         if(is == myid) then
            call send_xtoy(f,gt(1),nx,ix_s(is),ix_e(is),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
         else
            call send_xtoy(f,ft(1),nx,ix_s(is),ix_e(is),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)

            d_ft = ft

            call mpi_irecv(d_gt(1),nrecv,mpi_real8,d_ir,1,
     +                     mpi_comm_world,ireqr,ierr)
            call mpi_isend(d_ft(1),nsend,mpi_real8,d_is,1,
     +                     mpi_comm_world,ireqs,ierr)
            call mpi_wait(ireqs,istatus,ierr)
            call mpi_wait(ireqr,istatus,ierr)

            gt = d_gt

         endif
         call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(ir),iy_e(ir),iz1,iz2)
      enddo
 
      return
      end

! ====================================================================

      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
 
      ! ---- transpose array g(ny,ixs:ixe,iz1:iz2)
      !                 ---> f(nx,iys:iye,iz1:iz2)
      !      load balanced version, march 2017
 
      include 'mpif.h'
      integer istatus(mpi_status_size)
 
      real :: f(nx,iys:iye,iz1:iz2),
     +        g(ny,ixs:ixe,iz1:iz2)
      real :: ft(nx*(iye+1-iys)*(iz2-iz1+1)),
     +        gt(ny*(ixe+1-ixs)*(iz2-iz1+1))
 
      integer, intent(in), dimension(0:np-1) :: ix_s,ix_e,iy_s,iy_e
 
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
 
      ! ---- loop thru cpus on a slab
 
      iss   = (myid/ncpu_s)*ncpu_s
      do i = 1,ncpu_s
         is    = mod(myid - iss + i,ncpu_s) + iss
         ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
         nsend = (iy_e(is) - iy_s(is) + 1)*ik
         nrecv = (ix_e(ir) - ix_s(ir) + 1)*jk
         if(is == myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(is),iy_e(is),iz1,iz2)
         else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(is),iy_e(is),iz1,iz2)
 
            call mpi_irecv(ft(1),nrecv,mpi_real8,ir,1,
     +                     mpi_comm_world,ireqr,ierr)
            call mpi_isend(gt(1),nsend,mpi_real8,is,1,
     +                     mpi_comm_world,ireqs,ierr)
            call mpi_wait(ireqs,istatus,ierr)
            call mpi_wait(ireqr,istatus,ierr)
 
         endif
         call recv_ytox(f,ft(1),nx,ix_s(ir),ix_e(ir),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
      enddo
 
      return
      end

! ====================================================================

      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
 
      ! ---- grab correct chunk of array to be sent
 
      real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
 
      do k = izs,ize
         do j = iys,iye
         do i = ixs,ixe
            ft(i,j,k) = f(i,j,k)
         enddo
         enddo
      enddo
 
      return
      end

! ====================================================================

      subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
      real g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
 
      do k = izs,ize
         do j = iys,iye
         do i = ixs,ixe
            g(j,i,k) = gt(i,j,k)
         enddo
         enddo
      enddo

      return
      end

! ====================================================================

      subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
 
      ! ---- grab correct chunk of array to be sent
 
      real g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
 
      do k = izs,ize
         do i = ixs,ixe
         do j = iys,iye
            gt(j,i,k) = g(j,i,k)
         enddo
         enddo
      enddo
 
      return
      end

! ====================================================================

      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
 
      do k = izs,ize
         do i = ixs,ixe
         do j = iys,iye
            f(i,j,k) = ft(j,i,k)
         enddo
         enddo
      enddo

      return
      end

! ======================================================================

      subroutine build_wavenumbers
      use pars

      ncx = nnx/2 + 1
      ncy = nny/2 + 1

      do i = 1,nnx
         xkn(i) = dble(i-1)*pi2/xl
         if (i > ncx) xkn(i) = -dble(nnx-i+1)*pi2/xl
      enddo
      fn = 1.0/dble(nnx)
      do i = 1,nnx
         xk(i) = xkn(i)*fn
      enddo

      do i = 1,nny
         ykn(i) = dble(i-1)*pi2/yl
         if (i > ncy) ykn(i) = -dble(nny-i+1)*pi2/yl
      enddo
      fn = 1.0/dble(nny)
      do i = 1,nny
         yk(i) = ykn(i)*fn
      enddo

      return
      end subroutine build_wavenumbers

! ======================================================================

      subroutine xderivp(ax,trigx,xk,nnx,iys,iye)

      ! ---- choose fft routine to use for x-derivatives

      use pars, only : i_fft

      if (i_fft == 2) then
         call xd_cuda(ax,nnx,iys,iye)
      endif

      return
      end subroutine xderivp

! ======================================================================

      subroutine xd_cuda(ax,nx,iys,iye)

      ! ----- get multiple x derivatives using cuda routines

      use cufft
      use cufft_wrk

      integer, intent(in) :: nx,iys,iye
      real, intent(inout), dimension(nx,iys:iye) :: ax

      ! ---- move input array to device

!     fn = 1.0/dble(nx)
      fn = 1.0
      do j = iys,iye
         do i = 1,nx
            x_in(i,j) = ax(i,j)*fn
         enddo
      enddo

      ! ---- forward 1D fft in x for all [iys,iye]

      ierr = cufftExecD2Z(pln_xf, x_in, x_out)

      ! ---- spectral derivative

!$cuf kernel do <<< *,* >>>                    ! execute these loops on the device
      do j = iys,iye
         ii           = 1
         x_out(1,j)   = 0.0
         x_out(2,j)   = 0.0
         do i = 3,nx-1,2
            ii           = ii + 1
            temp         = x_out(i,j)
            x_out(i,j)   = -xk_d(ii)*x_out(i+1,j)
            x_out(i+1,j) = xk_d(ii)*temp
         enddo
         x_out(nx+1,j) = 0.0
         x_out(nx+2,j) = 0.0
      enddo

      ! ---- backward fft

      ierr = cufftExecZ2D(pln_xb, x_out, x_out)

      ! ---- copy data back to host

      do j = iys,iye
         do i = 1,nx
            ax(i,j) = x_out(i,j)
         enddo
      enddo

      return
      end subroutine xd_cuda

! ======================================================================

      subroutine yd_mpi(ay,trigy,yk,
     +           nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

      ! ---- choose fft routine to use for y-derivatives

      use pars, only : i_fft

      if (i_fft == 2) then
         call yd_cuda(ay,yk,nx,ny,ixs,ixe,ix_s,ix_e,
     +                iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
      endif

      return
      end subroutine yd_mpi

! ======================================================================

      subroutine yd_cuda(ay,yk,
     +                   nx,ny,ixs,ixe,ix_s,ix_e,
     +                   iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

      ! ---- get multiple y derivatives using cuda routines

      use cufft
      use cufft_wrk

      integer, intent(in) :: nx,ny,ixs,ixe,iys,iye,iz1,iz2,myid,ncpu,np
      integer, intent(in), dimension(0:np-1) :: ix_s, ix_e, iy_s, iy_e

      real, intent(in), dimension(ny) :: yk
      real, intent(inout), dimension(nx,iys:iye,iz1:iz2) :: ay

      real :: ayt(ny,ixs:ixe,iz1:iz2)

      ! ---- transpose so all y resides locally on each task

      call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

      ! ---- loop over z

      do k = iz1,iz2

         ! ---- copy y-x slab to device

!        fn = 1.0/dble(ny)
         fn = 1.0
         do i = ixs,ixe
            do j = 1,ny
               y_in(j,i) = ayt(j,i,k)*fn
            enddo
         enddo

         ! ---- forward 1D fft in y for all [ixs:ixe]

         ierr = cufftExecD2Z(pln_yf, y_in, y_out)

         ! ---- spectral derivative

!$cuf kernel do <<< *,* >>>                    ! execute these loops on the device
         do i = ixs,ixe
            ii         = 1
            y_out(1,i) = 0.0
            y_out(2,i) = 0.0
            do j = 3,ny-1,2
               ii           = ii + 1
               temp         = y_out(j,i)
               y_out(j,i)   = -yk_d(ii)*y_out(j+1,i)
               y_out(j+1,i) = yk_d(ii)*temp
            enddo
            y_out(ny+1,i) = 0.0
            y_out(ny+2,i) = 0.0
         enddo

         ! ---- backward fft

         ierr = cufftExecZ2D(pln_yb, y_out, y_out)

         ! ---- copy data back to host

         do i = ixs,ixe
            do j = 1,ny
               ayt(j,i,k) = y_out(j,i)
            enddo
         enddo

      enddo

      ! ---- transpose back so all x resides locally on each task

      call ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

      return
      end subroutine yd_cuda

! ===============================================================================

      subroutine gridd
      use pars

      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s


      allocate(ix_s(0:maxp), ix_e(0:maxp),
     +         jx_s(0:maxp), jx_e(0:maxp),
     +         kx_s(0:maxp), kx_e(0:maxp),
     +         mx_s(0:maxp), mx_e(0:maxp),
     +         iy_s(0:maxp), iy_e(0:maxp),
     +         jy_s(0:maxp), jy_e(0:maxp),
     +         is_s(0:maxp), is_e(0:maxp),
     +         iz_s(0:maxp), iz_e(0:maxp))

      return
      end subroutine gridd

! =================================================================

      subroutine mpi_range

      ! ---- get ranges for x and z variables
      !      note x range is based on nnx+2 fourier modes

      use pars

      ii = -1
      do nn = 0,ncpu_z-1
         call range(1,nnx+2,ncpu_z,nn,lx_s,lx_e)
         call range(1,nnx,ncpu_z,nn,nx_s,nx_e)
         call range(1,nny,ncpu_z,nn,ly_s,ly_e)
         call range(1,nnz,ncpu_z,nn,mz_s,mz_e)
         do mm = 0,ncpu_s-1
            call range(1,nny,ncpu_s,mm,ny_s,ny_e)
            call range(1,nnx,ncpu_s,mm,nxy_s,nxy_e)
            call range(1,ncx,ncpu_s,mm,l2x_s,l2x_e)
            ii       = ii + 1

            ix_s(ii) = nxy_s
            ix_e(ii) = nxy_e
            jx_s(ii) = (l2x_s - 1)*2 + 1
            jx_e(ii) = l2x_e*2
            kx_s(ii) = lx_s
            kx_e(ii) = lx_e
            mx_s(ii) = nx_s
            mx_e(ii) = nx_e
 
            iy_s(ii) = ny_s
            iy_e(ii) = ny_e
            jy_s(ii) = ly_s
            jy_e(ii) = ly_e
 
            iz_s(ii) = mz_s
            iz_e(ii) = mz_e
 
            is_s(ii) = (ii/ncpu_s)*ncpu_s
            is_e(ii) = is_s(ii) + ncpu_s - 1
         enddo
      enddo

      iys = iy_s(myid)
      iye = iy_e(myid)
      jys = jy_s(myid)
      jye = jy_e(myid)
      ixs = ix_s(myid)
      ixe = ix_e(myid)
      jxs = jx_s(myid)
      jxe = jx_e(myid)
      kxs = kx_s(myid)
      kxe = kx_e(myid)
      mxs = mx_s(myid)
      mxe = mx_e(myid)
      izs = iz_s(myid)
      ize = iz_e(myid)

!     write(*,123)myid,izs,ize,iys,iye,iss,ise
!123  format('myid = ',i6,' izs,ize = ',2i6,
!    +       ' iys,iye = ',2i6,' iss,ise = ',2i6)

!     write(*,124)myid,ixs,ixe,jxs,jxe
!124  format('myid = ',i6,' ixs,ixe = ',2i6,
!    +       ' jxs,jxe = ',2i6)

      return
      end

! =================================================================

      subroutine range(n1,n2,nprocs,irank,ista,iend)
 
      ! ---- the ibm range finder to balance load
 
      iwork1 = (n2 - n1 + 1)/nprocs
      iwork2 = mod(n2 - n1 +1, nprocs)
      ista = irank*iwork1 + n1 + min(irank,iwork2)
      iend = ista + iwork1 - 1
      if (iwork2 > irank) iend = iend + 1
 
      return
      end

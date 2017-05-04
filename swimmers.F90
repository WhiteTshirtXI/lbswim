!************************************************************************
!*                                                                      *
!*     InitSwimmers                                                     *
!*                                                                      *
!************************************************************************

! ... Initializes swimmer positions and orientations 

subroutine InitSwimmers
 
   use LBModule

   implicit none

   integer(4) :: i, ierr, buflen
   real(8)    :: ori(3), nsq
   real(8)    :: Random
 
#if defined (MPI)
   call mpi_bcast(nSwim,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(l,1,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(fswim,1,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(vswim,1,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(lambda,1,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(ltumbles,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(lswims,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(ladvects,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(lrotates,1,mpi_logical,rootid,comm,ierr)
#endif

!  fswim = 6*Pi*eta*a*vswim
   tumbleProb = vswim/lambda 

   allocate(r(1:3,1:nSwim)) 
   allocate(n(1:3,1:nSwim)) 
   buflen = max(nSwim,nTrac) 
   allocate(rbuf(1:3,buflen)) 

   if(master) then                           ! ... Let master do everything
      do i = 1, nSwim
         r(1,i) = Random(iseed)*nx           ! ... Swimmer domain is bounded by (0,nx)
         r(2,i) = Random(iseed)*ny 
         r(3,i) = Random(iseed)*nz

         ori(1) = Random(iseed)-0.5d0
         ori(2) = Random(iseed)-0.5d0
         ori(3) = Random(iseed)-0.5d0
         nsq = ori(1)**2 + ori(2)**2 + ori(3)**2
         n(1:3,i) = ori(1:3)/sqrt(nsq)
      end do
   end if

#if defined (MPI)
   call mpi_bcast(r(1:3,1:nSwim),3*nSwim,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(n(1:3,1:nSwim),3*nSwim,mpi_real8,rootid,comm,ierr)
#endif

end subroutine

!************************************************************************
!*                                                                      *
!*     UpdateSwimmers                                                   *
!*                                                                      *
!************************************************************************

! ... Integrates swimmer positions and orientations 

subroutine UpdateSwimmers
 
   use LBModule

   implicit none
   
   integer(4) :: i

   real(8) :: rdot(3), ndot(3), vplus(3), vminus(3), rminus(3), norm
   real(8) :: tumble, ori(3), nsq
   real(8) :: Random 
 
   do i = iplow(myid), ipupp(myid)

      rdot(1:3) = 0.0d0

!     if(lswims) rdot(1:3) = rdot(1:3) + fswim*n(1:3,i)*(1.0d0/a - 1.0d0/h)/(6.0d0*Pi*eta)
      if(lswims) rdot(1:3) = rdot(1:3) + vswim*n(1:3,i)

      if(ladvects) then
         call velocity(r(1:3,i),vplus)
         rdot(1:3) = rdot(1:3) + vplus(1:3)
      end if

      r(1:3,i) = r(1:3,i) + rdot(1:3)  
      r(1,i) = r(1,i) - int(r(1,i)/nx)*nx         ! ... PBCs - box going from 0 to nx
      r(2,i) = r(2,i) - int(r(2,i)/ny)*ny
      r(3,i) = r(3,i) - int(r(3,i)/nz)*nz

      tumble = Random(iseed)
      if(ltumbles .and. (tumble < tumbleProb)) then     ! ... Tumble with probability tumbleProb
         ori(1) = Random(iseed)-0.5d0
         ori(2) = Random(iseed)-0.5d0
         ori(3) = Random(iseed)-0.5d0
         nsq = ori(1)**2 + ori(2)**2 + ori(3)**2
         n(1:3,i) = ori(1:3)/sqrt(nsq) 
      else if(lrotates) then                            ! ... Otherwise just rotate with the flow
         rminus(1:3) = r(1:3,i)-l*n(1:3,i)
         call velocity(rminus,vminus)
         ndot(1:3) = (vplus(1:3)-vminus(1:3))/l
         n(1:3,i) = n(1:3,i) + ndot(1:3)
         norm = dsqrt(n(1,i)**2 + n(2,i)**2 + n(3,i)**2)
         n(1:3,i) = n(1:3,i)/norm  
      end if

   end do

end subroutine

!************************************************************************
!*                                                                      *
!*     Velocity                                                         *
!*                                                                      *
!************************************************************************
 
! ... Returns the interpolated fluid velocity in a point
 
subroutine velocity(r0,v0)
   
   use LBModule

   implicit none

   real(8), intent(in)  :: r0(3)   
   real(8), intent(out) :: v0(3)

   integer(4) :: lat(3), x0, x, ind(1:3,0:3)
   integer(4) :: i,j,k,d
   real(8)    :: deltas(1:3,0:3), delta3d
   real(8)    :: Peskin

   lat(1) = nx
   lat(2) = ny
   lat(3) = nz

   do d = 1, 3                                                     ! ... Loop over the full support (2*2*2) of the Peskin delta
      x0 = ceiling(r0(d)-2.0d0)
      do i = 0, 3                                                 
         x = x0+i                                               
         ind(d,i) = x                                              ! ... Indices of the lattice sites in the support
         if(ind(d,i) < 0) ind(d,i) = ind(d,i) + lat(d)             ! ... Periodic boundary conditions
         if(ind(d,i) > (lat(d)-1)) ind(d,i) = ind(d,i) - lat(d)
         deltas(d,i) = Peskin(r0(d)-x)                              
      end do
   end do

   v0(1:3) = 0.0d0
   do k = 0, 3                                                     ! ... delta3d(x,y,z) = delta(x)*delta(y)*delta(z)
      do j = 0, 3
         do i = 0, 3
            delta3d = deltas(1,i)*deltas(2,j)*deltas(3,k)
            do d = 1, 3
               v0(d) = v0(d) + delta3d*u(d,ind(1,i),ind(2,j),ind(3,k))
            end do
         end do
      end do
   end do

end subroutine


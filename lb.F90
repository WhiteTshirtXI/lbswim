!************************************************************************
!*                                                                      *
!*     InitLattice                                                      *
!*                                                                      *
!************************************************************************

! ... Sets weights and lattice vectors and initializes DFs 

subroutine InitLattice
 
   use LBModule

   implicit none

   integer(4) :: i, j, k, p, ierr
   

#if defined (MPI)
   call mpi_bcast(nx,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(ny,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(nz,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(eta,1,mpi_real8,rootid,comm,ierr)
#endif

   do i = 0, nproc-1
      zlow(i) = i*nz/nproc
      zupp(i) = (i+1)*nz/nproc-1
      zlen(i) = zupp(i)-zlow(i)+1
   end do

   allocate(f(0:14,0:nx-1,0:ny-1,0:nz-1))
   allocate(force(1:3,0:nx-1,0:ny-1,0:nz-1)) 
   allocate(u(1:3,0:nx-1,0:ny-1,0:nz-1)) 
   allocate(rho(0:nx-1,0:ny-1,0:nz-1)) 
   allocate(help(0:14,0:nx-1,0:ny-1,0:nz-1)) 

   tau = eta*3.0d0
   omega = 1.0d0/(0.5d0 + tau)

   force(1:3,0:nx-1,0:ny-1,0:nz-1) = 0.0d0
   u(1:3,0:nx-1,0:ny-1,0:nz-1) = 0.0d0

   ! ... Set up lattice weights

   w(0) = 2.0d0/9.0d0
   w(1:6) = 1.0d0/9.0d0
   w(7:14) = 1.0d0/72.0d0 

   ! ... Velocity vectors

   ci(1:3,0) = (/ 0, 0, 0 /)

   ci(1:3,1) = (/ 1, 0, 0 /)
   ci(1:3,2) = (/-1, 0, 0 /)
   ci(1:3,3) = (/ 0, 1, 0 /)
   ci(1:3,4) = (/ 0,-1, 0 /)
   ci(1:3,5) = (/ 0, 0, 1 /)
   ci(1:3,6) = (/ 0, 0,-1 /)

   ci(1:3,7)  = (/ 1, 1, 1 /)
   ci(1:3,8)  = (/ 1, 1,-1 /)
   ci(1:3,9)  = (/ 1,-1, 1 /)
   ci(1:3,10) = (/ 1,-1,-1 /)
   ci(1:3,11) = (/-1, 1, 1 /)
   ci(1:3,12) = (/-1, 1,-1 /)
   ci(1:3,13) = (/-1,-1, 1 /)
   ci(1:3,14) = (/-1,-1,-1 /)

   ! Initialize DFs

   do k = 0, nz-1
      do j = 0, ny-1
         do i = 0, nx-1
            f(0:14,i,j,k) = w(0:14)
         end do
      end do
   end do

end subroutine

!************************************************************************
!*                                                                      *
!*     Collide                                                          *
!*                                                                      *
!************************************************************************

! ... Relaxes distribution functions

subroutine Collide
 
   use LBModule

   implicit none

   integer(4) :: i, j, k, p
   real(8) :: feq(0:14), phi(0:14)
#if defined (MPI)
   integer(4) :: ierr, sreq1, sreq2, rreq1, rreq2, istatus(mpi_status_size)
#endif

   do k = zlow(myid), zupp(myid)
      do j = 0, ny-1
         do i = 0, nx-1
            call CalcHydroSite(i,j,k)
            call CalcEquil(i,j,k,feq)
            call CalcPhi(i,j,k,phi)
            do p = 0, 14
               f(p,i,j,k) = f(p,i,j,k) - omega*(f(p,i,j,k)-feq(p)) + omega*tau*phi(p)
            end do
         end do
      end do
   end do
#if defined (MPI)
   call mpi_isend(f(0:14,0:nx-1,0:ny-1,zlow(myid)), 15*nx*ny,mpi_real8,cpudown,0,comm,sreq1,ierr)
   call mpi_irecv(f(0:14,0:nx-1,0:ny-1,zlow(cpuup)),15*nx*ny,mpi_real8,cpuup,mpi_any_tag,comm,rreq1,ierr)
   call mpi_isend(f(0:14,0:nx-1,0:ny-1,zupp(myid)),   15*nx*ny,mpi_real8,cpuup,0,comm,sreq2,ierr)
   call mpi_irecv(f(0:14,0:nx-1,0:ny-1,zupp(cpudown)),15*nx*ny,mpi_real8,cpudown,mpi_any_tag,comm,rreq2,ierr)
   call mpi_wait(sreq1,istatus,ierr)
   call mpi_wait(rreq1,istatus,ierr)
   call mpi_wait(sreq2,istatus,ierr)
   call mpi_wait(rreq2,istatus,ierr)
#endif

end subroutine

!************************************************************************
!*                                                                      *
!*     Stream                                                           *
!*                                                                      *
!************************************************************************

! ... Propagates distribution functions (with periodic boundary conditions)

subroutine Stream
 
   use LBModule

   implicit none

   integer(4) :: i, j, k, p, imod, jmod, kmod

   help(0:14,0:nx-1,0:ny-1,0:nz-1) = f(0:14,0:nx-1,0:ny-1,0:nz-1) 
   do k = zlow(myid), zupp(myid)
      do j = 0, ny-1
         do i = 0, nx-1
            do p = 0, 14
               imod = mod((i-ci(1,p)+nx),nx)
               jmod = mod((j-ci(2,p)+ny),ny)
               kmod = mod((k-ci(3,p)+nz),nz)
               f(p,i,j,k) = help(p,imod,jmod,kmod)
            end do 
         end do
      end do
   end do

end subroutine

!************************************************************************
!*                                                                      *
!*     UpdateHydroVars                                                  *
!*                                                                      *
!************************************************************************

! ... Calculates hydrodynamic variables for the whole lattice 

subroutine UpdateHydroVars
 
   use LBModule

   implicit none

   integer(4) :: i, j, k, ierr
   
   do k = zlow(myid), zupp(myid)
      do j = 0, ny-1
         do i = 0, nx-1
            call CalcHydroSite(i,j,k)
         end do
      end do
   end do

#if defined (MPI)
   call mpi_allgatherv(u(1:3,0:nx-1,0:ny-1,zlow(myid):zupp(myid)),3*nx*ny*zlen(myid),mpi_real8,help(1:3,0:nx-1,0:ny-1,0:nz-1),3*nx*ny*zlen(0:nproc-1),3*nx*ny*zlow(0:nproc-1),mpi_real8,comm,ierr)
   u(1:3,0:nx-1,0:ny-1,0:nz-1) = help(1:3,0:nx-1,0:ny-1,0:nz-1)
#endif

end subroutine


!************************************************************************
!*                                                                      *
!*     UpdateForces                                                     *
!*                                                                      *
!************************************************************************

! ... Interpolates forces from swimmers to the lattice 

subroutine UpdateForces
 
   use LBModule

   implicit none
   
   integer(4) :: i, ierr
   real(8)    :: rplus(3), rminus(3), fminus(3), swimforce(3)
#if defined (MPI)
   integer(4) :: sreq, rreq, istatus(mpi_status_size)
#endif

   force(1:3,0:nx-1,0:ny-1,0:nz-1) = 0.0d0
   do i = iplow(myid), ipupp(myid)
      swimforce(1:3) = fswim*n(1:3,i) 
      call AddForce(r(1:3,i),swimforce(1:3))   ! ... Force from swimmer body
      rminus(1:3) = r(1:3,i)-l*n(1:3,i)
      fminus(1:3) = -swimforce(1:3)            ! ... Force from tail
      call AddForce(rminus,fminus)
   end do

#if defined (MPI)
   call mpi_allreduce(force(1:3,0:nx-1,0:ny-1,0:nz-1),help(1:3,0:nx-1,0:ny-1,0:nz-1),3*nx*ny*nz,mpi_real8,mpi_sum,comm,ierr)
   force(1:3,0:nx-1,0:ny-1,zlow(myid):zupp(myid)) = help(1:3,0:nx-1,0:ny-1,zlow(myid):zupp(myid))
#endif

!  call mpi_reduce(force(1:3,0:nx-1,0:ny-1,0:nz-1),help(1:3,0:nx-1,0:ny-1,0:nz-1),3*nx*ny*nz,mpi_real8,mpi_sum,rootid,comm,ierr)
!  if(master) force(1:3,0:nx-1,0:ny-1,zlow(myid):zupp(myid)) = help(1:3,0:nx-1,0:ny-1,zlow(myid):zupp(myid))

!  do i = 0, nproc-1
!     if(i == rootid) cycle
!     if(master)    call mpi_isend(help(1:3,0:nx-1,0:ny-1,zlow(i):zupp(i)),3*nx*ny*zlen(i),mpi_real8,i,i,comm,sreq,ierr)
!     if(i == myid) call mpi_irecv(force(1:3,0:nx-1,0:ny-1,zlow(myid):zupp(myid)),3*nx*ny*zlen(myid),mpi_real8,rootid,myid,comm,rreq,ierr)
!     if(master) call mpi_wait(sreq,istatus,ierr)
!     if(i == myid) call mpi_wait(rreq,istatus,ierr)
!     print*, i, myid, sreq, rreq
!  end do
!  if(master) then
!     do i = 0, nproc-1
!        if(i == rootid) cycle
!        call mpi_wait(sreq(i),istatus,ierr)
!     end do
!  else
!     call mpi_wait(rreq,istatus,ierr)
!  end if
      
!  print*, myid, sreq(1:3), rreq(1:3)  
  
!  do i = 0, nproc-1
!     if(i == rootid) cycle
!     if(master) call mpi_wait(sreq(i),istatus,ierr)
!     if(i == myid) call mpi_wait(rreq(myid),istatus,ierr)
!  end do

end subroutine

!************************************************************************
!*                                                                      *
!*     AddForce                                                         *
!*                                                                      *
!************************************************************************

! ... Interpolates a single off-lattice force to the lattice using the Peskin delta

subroutine AddForce(r0,f0)
 
   use LBModule

   implicit none
   
   real(8), intent(in) :: r0(3), f0(3)

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

   do k = 0, 3                                                     ! ... delta3d(x,y,z) = delta(x)*delta(y)*delta(z)
      do j = 0, 3
         do i = 0, 3
            delta3d = deltas(1,i)*deltas(2,j)*deltas(3,k)
            do d = 1, 3
               force(d,ind(1,i),ind(2,j),ind(3,k)) = force(d,ind(1,i),ind(2,j),ind(3,k)) + delta3d*f0(d)
            end do
         end do
      end do
   end do

end subroutine


!************************************************************************
!*                                                                      *
!*     CalcHydroSite                                                    *
!*                                                                      *
!************************************************************************

! ... Calculates hydrodynamic variables for a single lattice site

subroutine CalcHydroSite(i,j,k)
 
   use LBModule

   implicit none

   integer(4), intent(in) :: i, j, k
   integer(4) :: p, d

   real(8) :: mom(1:3), rholoc 

   rholoc = 0.0d0
   mom(1:3) = 0.5d0*force(1:3,i,j,k)
   do p = 0, 14
      rholoc = rholoc + f(p,i,j,k)
      do d = 1, 3
         mom(d) = mom(d) + f(p,i,j,k)*ci(d,p)
      end do
   end do 
          
   rho(i,j,k) = rholoc
   u(1:3,i,j,k) = mom(1:3)/rholoc

end subroutine

!************************************************************************
!*                                                                      *
!*     CalcEquil                                                        *
!*                                                                      *
!************************************************************************

! ... Calculates equilibrium distribution functions

subroutine CalcEquil(i,j,k,feq)
 
   use LBModule

   implicit none

   integer(4), intent(in) :: i, j, k
   real(8), intent(out)   :: feq(0:14)

   integer(4) :: p
   real(8) :: term1, term2, udotc, usqhalf, ueq(3)

   ueq(1:3) = u(1:3,i,j,k) 
   usqhalf = 0.5d0*(ueq(1)**2 + ueq(2)**2 + ueq(3)**2)

   do p = 0, 14
      udotc = ueq(1)*ci(1,p) + ueq(2)*ci(2,p) + ueq(3)*ci(3,p)
      term1 = (udotc - usqhalf)/cs2     
      term2 = (udotc*udotc)/(2.0d0*cs2*cs2)     
      feq(p) = w(p)*rho(i,j,k)*(1.0d0 + term1 + term2)
   end do

end subroutine

!************************************************************************
!*                                                                      *
!*     CalcPhi                                                          *
!*                                                                      *
!************************************************************************

! ... Calculates forcing term

subroutine CalcPhi(i,j,k,phi)
 
   use LBModule

   implicit none

   integer(4), intent(in) :: i, j, k
   real(8), intent(out)   :: phi(0:14)

   integer(4) :: p
   real(8) :: udotc, fdotc, udotf

   udotf = u(1,i,j,k)*force(1,i,j,k) + u(2,i,j,k)*force(2,i,j,k) + u(3,i,j,k)*force(3,i,j,k)
   do p = 0, 14
      udotc = u(1,i,j,k)*ci(1,p) + u(2,i,j,k)*ci(2,p) + u(3,i,j,k)*ci(3,p)
      fdotc = force(1,i,j,k)*ci(1,p) + force(2,i,j,k)*ci(2,p) + force(3,i,j,k)*ci(3,p)
      phi(p) = w(p)*(fdotc - udotf + udotc*fdotc/cs2)/cs2
   end do

end subroutine

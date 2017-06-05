!************************************************************************
!*                                                                      *
!*     InitTracers                                                      *
!*                                                                      *
!************************************************************************

! ... Initializes swimmer positions and orientations 

subroutine InitTracers 
 
   use LBModule

   implicit none

   integer(4) :: i, ierr
   real(8)    :: fran(3)

#if defined (MPI)
   call mpi_bcast(nTrac,1,mpi_integer4,rootid,comm,ierr)
#endif

   if(nTrac <= 0) return

   allocate(rtr(1:3,1:nTrac)) 
   allocate(str(1:3,1:nTrac)) 

   do i = 0, nproc-1
      trlow(i) = i*nTrac/nproc + 1
      trupp(i) = (i+1)*nTrac/nproc
      trlen(i) = trupp(i)-trlow(i)+1
   end do

   if(master) then                           ! ... Let master do everything
      do i = 1, nTrac
         call random_number(fran)
         rtr(1,i) = fran(1)*nx           ! ... Swimmer domain is bounded by (0,nx)
         rtr(2,i) = fran(2)*ny 
         rtr(3,i) = fran(3)*nz
      end do
   end if

#if defined (MPI)
   call mpi_bcast(rtr(1:3,1:nSwim),3*nTrac,mpi_real8,rootid,comm,ierr)
#endif

   str(1:3,1:nTrac) = rtr(1:3,1:nTrac)

end subroutine

!************************************************************************
!*                                                                      *
!*     UpdateTracers                                                    *
!*                                                                      *
!************************************************************************

! ... Integrates tracer positions and orientations 

subroutine UpdateTracers
 
   use LBModule

   implicit none
   
   integer(4) :: i

   real(8) :: rdot(3) 
   
   if(nTrac <= 0) return

   do i = trlow(myid), trupp(myid)
      call velocity(rtr(1:3,i),rdot)
      rtr(1:3,i) = rtr(1:3,i) + rdot(1:3)  
      rtr(1,i) = rtr(1,i) - floor(rtr(1,i)/nx)*nx         ! ... PBCs - box going from 0 to nx
      rtr(2,i) = rtr(2,i) - floor(rtr(2,i)/ny)*ny
      rtr(3,i) = rtr(3,i) - floor(rtr(3,i)/nz)*nz
      str(1:3,i) = str(1:3,i) + rdot(1:3)  
   end do

end subroutine

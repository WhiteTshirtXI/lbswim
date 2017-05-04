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
   real(8)    :: Random

#if defined (MPI)
   call mpi_bcast(nTrac,1,mpi_integer4,rootid,comm,ierr)
#endif

   if(nTrac <= 0) return

   allocate(rtr(1:3,1:nTrac)) 
   allocate(str(1:3,1:nTrac)) 

   if(master) then                           ! ... Let master do everything
      do i = 1, nTrac
         rtr(1,i) = Random(iseed)*nx           ! ... Swimmer domain is bounded by (0,nx)
         rtr(2,i) = Random(iseed)*ny 
         rtr(3,i) = Random(iseed)*nz
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
      rtr(1,i) = rtr(1,i) - int(r(1,i)/nx)*nx         ! ... PBCs - box going from 0 to nx
      rtr(2,i) = rtr(2,i) - int(r(2,i)/ny)*ny
      rtr(3,i) = rtr(3,i) - int(r(3,i)/nz)*nz
      str(1:3,i) = str(1:3,i) + rdot(1:3)  
   end do

end subroutine

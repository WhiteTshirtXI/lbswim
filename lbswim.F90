!************************************************************************
!************************************************************************
!**                                                                    **
!**  Joakim Stenhammar                                                 **
!**  Physical Chemistry                                                **
!**  Center for Chemistry and Chemical Engineering                     **
!**  Lund University                                                   **
!**  Sweden                                                            **
!**                                                                    **
!**  All rights reserved. The code may not be modified or              **
!**  redistributed without the written conscent of the copyright       **
!**  owner. The copyright owner does not take any responsibility       **
!**  for any error in the code or its documentation.                   **
!**                                                                    **
!************************************************************************
!************************************************************************

! ... D3Q15 Lattice Boltzmann code with forcing

!************************************************************************
!*                                                                      *
!*     LBSwim                                                           *
!*                                                                      *
!************************************************************************

! ... Main routine

program LBSwim 
   use LBModule

   implicit none

   real(8) :: usq, starttime, endtime, interval, tdump
   integer(4) :: uin, uout, ulog, istep, ierr, i 
   character(30) :: fin, fout, flog
   logical :: lexists

#if defined (MPI)
   call mpi_init(ierr)
   call mpi_comm_size(comm,nproc,ierr)
   if(nproc > mnproc) then
      write(*,*) 'Too many threads'
      stop
   end if
   call mpi_comm_rank(comm,myid,ierr)
   master = .false.
   if(myid == rootid) master = .true.
   if(master) starttime = mpi_wtime()
#else
   nproc = 1
   myid = 0   
   master = .true.
   call cpu_time(starttime) 
#endif  


   namelist /nmlRun/ nstep, iseed

   namelist /nmlLB/ nx, ny, nz, eta

   namelist /nmlSwimmers/ nSwim, l, fswim, vswim, lambda, ltumbles, lswims, ladvects, lrotates

   namelist /nmlTracers/ nTrac 

   namelist /nmlIO/ idump, ldumpswim, ldumplat, lformatted, lrestore
 
   uin   = 1
   uout  = 2
   ulog  = 3
   fin   = "lb.in"
   flog  = "lb.log"
   interval = 3600.0d0
   startstep = 1

! .. open units
   if(master) open(uin, file = fin)
   if(master) open(ulog, file = flog)

   if(master) read(uin,nmlRun)
#if defined (MPI)
   call mpi_bcast(nstep,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(iseed,1,mpi_integer4,rootid,comm,ierr)
#endif
   iseed = iseed + myid          ! ... Ensure each process has its own random seed

   if(master) read(uin,nmlLB)
   call InitLattice

   if(master) read(uin,nmlSwimmers)
   call InitSwimmers

   nTrac = 0
   if(master) read(uin,nmlTracers)
   call InitTracers

   if(master) read(uin,nmlIO)
#if defined (MPI)
   call mpi_bcast(idump,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(ldumpswim,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(ldumplat,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(lformatted,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(lrestore,1,mpi_logical,rootid,comm,ierr)
#endif

   if(master) close(uin)

   inquire(file='checkpoint.out',exist=lexists) 
   if(lrestore .and. lexists) call RestoreFromCP

   call InitParallel

   tdump = starttime
   do istep = startstep, nstep
      call UpdateForces
      call Collide
      call Stream
      call UpdateHydroVars
      call UpdateSwimmers
      call UpdateTracers
      if(mod(istep,idump) == 0) then
         call DoDump(uout,istep)
         usq = sum(u(1:3,0:nx-1,0:ny-1,0:nz-1)**2)/(nx*ny*nz)
         if(master) write(ulog,'(i8,es20.10)') istep, dsqrt(usq)/vswim
         flush(ulog)
      end if
      call Checkpoint(tdump,interval,istep+1)
   end do


#if defined (MPI)
   if(master) endtime = mpi_wtime()
#else
   call cpu_time(endtime) 
#endif 

   if(master) write(ulog,'(a16,i4)'), trim('Number of tasks:'), nproc
   if(master) write(ulog,'(a19,f15.3)'), trim('Execution time (s):'), endtime-starttime
   if(master) close(ulog)

#if defined (MPI)
   call mpi_finalize(ierr)
#endif

end program LBSwim 


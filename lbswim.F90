!************************************************************************
!************************************************************************
!**                                                                    **
!**  Joakim Stenhammar                                                 **
!**  Physical Chemistry                                                **
!**  Center for Chemistry and Chemical Engineering                     **
!**  Lund University                                                   **
!**  Sweden                                                            **
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

   real(8) :: usq, starttime, endtime, interval, tdump, elapsed, rate
   integer(4) :: uin, uout, ulog, istep, ierr, i, ran_seed(2), istart, iend, buflen 
   character(40) :: fin, fout, flog, fchk, basename, txpath
   logical :: lexists, lmakedir

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
   call system_clock(istart) 
#endif  
   cpudown = myid-1
   if(myid == 0) cpudown = nproc-1
   cpuup = myid+1
   if(myid == nproc-1) cpuup = 0

   namelist /nmlRun/ nstep, iseed

   namelist /nmlLB/ nx, ny, nz, eta

   namelist /nmlSwimmers/ nSwim, l, fswim, vswim, lambda, ltumbles, lswims, ladvects, lrotates

   namelist /nmlTracers/ nTrac 

   namelist /nmlIO/ idump, ldumpswim, ldumplat, lformatted, lrestore
 
   if(master) then
      call get_command_argument(1,basename)
      if(len_trim(basename) == 0) then
         write(*,*) "ERROR! Syntax: program name + project name"
         stop
      end if
   end if
#if defined (MPI)
   call mpi_bcast(basename,40,mpi_character,rootid,comm,ierr)
#endif

   uin   = 1
   uout  = 2
   ulog  = 3
   fin   = trim(basename)//'.in'
   flog  = trim(basename)//'.log'
   fchk  = trim(basename)//'.chk'
   txpath = './'//trim(basename)
   interval = 3600.0d0  ! ... Checkpointing interval (seconds)
   startstep = 1

   lmakedir = makedirqq(txpath)

! .. open units
   if(master) open(uin, file = fin)
   if(master) open(ulog, file = flog)

   if(master) read(uin,nmlRun)
#if defined (MPI)
   call mpi_bcast(nstep,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(iseed,1,mpi_integer4,rootid,comm,ierr)
#endif
   iseed = iseed + myid          ! ... Ensure each process has its own random seed
   ran_seed(1) = iseed
   ran_seed(2) = myid
   call random_seed(PUT=ran_seed)

   if(master) read(uin,nmlLB)
   call InitLattice

   if(master) read(uin,nmlSwimmers)
   call InitSwimmers

   nTrac = 0
   if(master) read(uin,nmlTracers)
   call InitTracers

   buflen = max(nSwim,nTrac) 
   allocate(rbuf(1:3,buflen)) 

   if(master) read(uin,nmlIO)
#if defined (MPI)
   call mpi_bcast(idump,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(ldumpswim,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(ldumplat,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(lformatted,1,mpi_logical,rootid,comm,ierr)
   call mpi_bcast(lrestore,1,mpi_logical,rootid,comm,ierr)
#endif

   if(master) close(uin)

   inquire(file=fchk,exist=lexists) 
   if(lrestore .and. lexists) call RestoreFromCP(fchk)

   if(master) tdump = starttime
   do istep = startstep, nstep
      call UpdateForces
      call Collide
      call Stream
      call UpdateHydroVars
      call UpdateSwimmers
      call UpdateTracers
      if(mod(istep,idump) == 0) then
         call DoDump(uout,txpath,istep)
         usq = sum(u(1:3,0:nx-1,0:ny-1,0:nz-1)**2)/(nx*ny*nz)
         if(master) write(ulog,'(i8,es20.10)') istep, dsqrt(usq)/vswim
         flush(ulog)
      end if
      call Checkpoint(tdump,interval,istep+1,fchk)
   end do


#if defined (MPI)
   if(master) then
      endtime = mpi_wtime()
      elapsed = endtime - starttime
   end if
#else
   call system_clock(iend,rate)
   elapsed = (iend-istart)/rate
#endif 

   if(master) write(ulog,'(a16,i4)'), trim('Number of tasks:'), nproc
   if(master) write(ulog,'(a19,f15.3)'), trim('Execution time (s):'), elapsed
   if(master) close(ulog)

#if defined (MPI)
   call mpi_finalize(ierr)
#endif

end program LBSwim 


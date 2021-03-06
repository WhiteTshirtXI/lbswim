!************************************************************************
!*                                                                      *
!*     DoDump                                                           *
!*                                                                      *
!************************************************************************

! ... Dump data to disk

subroutine DoDump(unit,path,step)

   use LBModule
  
   implicit none

   integer(4), intent(in) :: unit, step
   character(40), intent(in) :: path
   character(50) :: root, outfile, txformat, txifile
   integer(4) :: i, j, k, ierr

#if defined (MPI)
   if(ldumpswim) then
      call mpi_gatherv(r(1:3,iplow(myid):ipupp(myid)),3*iplen(myid),mpi_real8,rbuf(1:3,1:nSwim),3*iplen(0:nproc-1),3*(iplow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
      if(master) r(1:3,1:nSwim) = rbuf(1:3,1:nSwim)
      call mpi_gatherv(n(1:3,iplow(myid):ipupp(myid)),3*iplen(myid),mpi_real8,rbuf(1:3,1:nSwim),3*iplen(0:nproc-1),3*(iplow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
      if(master) n(1:3,1:nSwim) = rbuf(1:3,1:nSwim)
   endif

   if(nTrac > 0) then
      call mpi_gatherv(str(1:3,trlow(myid):trupp(myid)),3*trlen(myid),mpi_real8,rbuf(1:3,1:nTrac),3*trlen(0:nproc-1),3*(trlow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
      if(master) str(1:3,1:nTrac) = rbuf(1:3,1:nTrac)
   endif
#endif

   if(.not. master) return

   if(step < 10) then
      txformat = '(i1)'
      write(txifile,txformat) step
      root  = trim(path)//'/00000'//trim(txifile)
   else if(step < 100) then
      txformat = '(i2)'
      write(txifile,txformat) step
      root  = trim(path)//'/0000'//trim(txifile)
   else if(step < 1000) then
      txformat = '(i3)'
      write(txifile,txformat) step
      root = trim(path)//'/000'//trim(txifile)
   else if(step < 10000) then
      txformat = '(i4)'
      write(txifile,txformat) step
      root = trim(path)//'/00'//trim(txifile)
   else if(step < 100000) then
      txformat = '(i5)'
      write(txifile,txformat) step
      root = trim(path)//'/0'//trim(txifile)
   else if(step < 1000000) then
      txformat = '(i6)'
      write(txifile,txformat) step
      root = trim(path)//'/'//trim(txifile)
   else 
      write(*,*) "nstep too big"
      stop
   end if

   if(ldumpswim) then
      outfile = trim(root)//'.swimmers' 
      if(lformatted) then
         open(unit, file = outfile, status = 'unknown', form = 'formatted', iostat = ierr)
         write(unit,'(3i8)') nx, ny, nz
         write(unit,'(i8)') nSwim
         do i = 1, nSwim
            write(unit,'(6f12.7)') r(1:3,i), n(1:3,i) 
         end do 
         close(unit) 
      else
         open(unit, file = outfile, status = 'unknown', form = 'unformatted', iostat = ierr)
         write(unit) nx, ny, nz
         write(unit) nSwim
         write(unit) r(1:3,1:nSwim) 
         write(unit) n(1:3,1:nSwim)
         close(unit) 
      end if
   end if

   if(ldumplat) then
      outfile = trim(root)//'.lattice'
      if(lformatted) then
         open(unit, file = outfile, status = 'unknown', form = 'formatted', iostat = ierr)
         write(unit,'(3i8)') nx, ny, nz
         do i = 0, nx-1
            do j = 0, ny-1
               do k = 0, nz-1
                  write(unit,'(3es17.7)') u(1:3,i,j,k)
               end do
            end do
         end do
         close(unit) 
      else
         open(unit, file = outfile, status = 'unknown', form = 'unformatted', iostat = ierr)
         write(unit) nx, ny, nz
         write(unit) u(1:3,0:nx-1,0:ny-1,0:nz-1)
         close(unit)
      end if
   end if

   if(nTrac > 0) then
      outfile = trim(root)//'.tracers' 
      if(lformatted) then
         open(unit, file = outfile, status = 'unknown', form = 'formatted', iostat = ierr)
         write(unit,'(i8)') nTrac
         do i = 1, nTrac
            write(unit,'(3f12.7)') str(1:3,i)
         end do 
         close(unit) 
      else
         open(unit, file = outfile, status = 'unknown', form = 'unformatted', iostat = ierr)
         write(unit) nTrac
         write(unit) str(1:3,1:nTrac) 
         close(unit) 
      end if
   end if
end subroutine

!************************************************************************
!*                                                                      *
!*     Checkpoint                                                       *
!*                                                                      *
!************************************************************************

! ... Checks if checkpointing is required

subroutine Checkpoint(tdump,interval,step,filename)

   use LBModule

   implicit none

   character(40), intent(in) :: filename
   real(8), intent(inout) :: tdump
   real(8), intent(in)    :: interval, step
   real(8) :: t0
   integer(4) :: ierr
   logical :: lCP

#if defined (MPI)
   lCP = .false.
   if(master .and. (mpi_wtime()-tdump) > interval) lCP = .true.
   call mpi_bcast(lCP,1,mpi_logical,rootid,comm,ierr)
   if(lCP) then
      call DoCP(step,filename)
      if(master) tdump = mpi_wtime()
   end if
#else 
   call cpu_time(t0)
   if(t0-tdump > interval) then
      call DoCP(step,filename)
      call cpu_time(tdump)
   end if
#endif      
  
end subroutine

!************************************************************************
!*                                                                      *
!*     DoCP                                                             *
!*                                                                      *
!************************************************************************

! ... Dumps checkpoint file to disk

subroutine DoCP(step,filename)

   use LBModule

   implicit none

   integer(4), intent(in) :: step

   character(40), intent(in) :: filename
   integer(4) :: unit, ierr

   unit = 33

   if(master) then
      deallocate(fbuf)
      allocate(fbuf(0:14,0:nx-1,0:ny-1,0:nz-1))   ! ... master needs to collect global f array for output
   end if
#if defined (MPI)
   call mpi_gatherv(f(0:14,0:nx-1,0:ny-1,1:zlen(myid)),15*nx*ny*zlen(myid),mpi_real8,fbuf(0:14,0:nx-1,0:ny-1,0:nz-1),15*nx*ny*zlen(0:nproc-1),15*nx*ny*zlow(0:nproc-1),mpi_real8,rootid,comm,ierr)
   call mpi_gatherv(r(1:3,iplow(myid):ipupp(myid)),3*iplen(myid),mpi_real8,rbuf(1:3,1:nSwim),3*iplen(0:nproc-1),3*(iplow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
   if(master) r(1:3,1:nSwim) = rbuf(1:3,1:nSwim)
   call mpi_gatherv(n(1:3,iplow(myid):ipupp(myid)),3*iplen(myid),mpi_real8,rbuf(1:3,1:nSwim),3*iplen(0:nproc-1),3*(iplow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
   if(master) n(1:3,1:nSwim) = rbuf(1:3,1:nSwim)
   if(nTrac > 0) then
      call mpi_gatherv(rtr(1:3,trlow(myid):trupp(myid)),3*trlen(myid),mpi_real8,rbuf(1:3,1:nTrac),3*trlen(0:nproc-1),3*(trlow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
      if(master) rtr(1:3,1:nTrac) = rbuf(1:3,1:nTrac)
      call mpi_gatherv(str(1:3,trlow(myid):trupp(myid)),3*trlen(myid),mpi_real8,rbuf(1:3,1:nTrac),3*trlen(0:nproc-1),3*(trlow(0:nproc-1)-1),mpi_real8,rootid,comm,ierr)
      if(master) str(1:3,1:nTrac) = rbuf(1:3,1:nTrac)
   endif
#else
   fbuf(0:14,0:nx-1,0:ny-1,0:nz-1) = f(0:14,0:nx-1,0:ny-1,1:nz)
#endif

   if(.not. master) return 

   open(unit, file = trim(filename), status = 'unknown', form = 'unformatted', iostat = ierr)
   write(unit) step
   write(unit) nx, ny, nz
   write(unit) fbuf(0:14,0:nx-1,0:ny-1,0:nz-1)
   write(unit) nSwim
   write(unit) r(1:3,1:nSwim) 
   write(unit) n(1:3,1:nSwim)
   write(unit) nTrac
   if(nTrac > 0) write(unit) rtr(1:3,1:nTrac) 
   if(nTrac > 0) write(unit) str(1:3,1:nTrac) 
   close(unit) 
   if(master) then
      deallocate(fbuf)
      allocate(fbuf(0:14,0:nx-1,0:ny-1,0:zlen(myid)+1)) 
   end if

end subroutine

!************************************************************************
!*                                                                      *
!*     RestoreFromCP                                                    *
!*                                                                      *
!************************************************************************

! ... Restores checkpointed configuration from disk

subroutine RestoreFromCP(filename)

   use LBModule

   implicit none

   character(40), intent(in) :: filename
   integer(4) :: unit, ierr

   unit = 33
 
   if(master) then
      open(unit, file = trim(filename), status = 'unknown', form = 'unformatted', iostat = ierr)
      read(unit) startstep
      read(unit) nx, ny, nz
      deallocate(fbuf)
      allocate(fbuf(0:14,0:nx-1,0:ny-1,0:nz-1))   ! ... master needs to collect global f array for output
      read(unit) fbuf(0:14,0:nx-1,0:ny-1,0:nz-1)
      read(unit) nSwim
      read(unit) r(1:3,1:nSwim) 
      read(unit) n(1:3,1:nSwim)
      read(unit) nTrac
      if(nTrac > 0) read(unit) rtr(1:3,1:nTrac) 
      if(nTrac > 0) read(unit) str(1:3,1:nTrac) 
      close(unit) 
   end if
#if defined (MPI)
   call mpi_scatterv(fbuf(0:14,0:nx-1,0:ny-1,0:nz-1),15*nx*ny*zlen(0:nproc-1),15*nx*ny*zlow(0:nproc-1),mpi_real8,f(0:14,0:nx-1,0:ny-1,1:zlen(myid)),15*nx*ny*zlen(myid),mpi_real8,rootid,comm,ierr)
   call mpi_bcast(startstep,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(nx,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(ny,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(nz,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(nSwim,1,mpi_integer4,rootid,comm,ierr)
   call mpi_bcast(r(1:3,1:nSwim),3*nSwim,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(n(1:3,1:nSwim),3*nSwim,mpi_real8,rootid,comm,ierr)
   call mpi_bcast(nTrac,1,mpi_integer4,rootid,comm,ierr)
   if(nTrac > 0) call mpi_bcast(rtr(1:3,1:nTrac),3*nTrac,mpi_real8,rootid,comm,ierr)
   if(nTrac > 0) call mpi_bcast(str(1:3,1:nTrac),3*nTrac,mpi_real8,rootid,comm,ierr)
#else
   f(0:14,0:nx-1,0:ny-1,1:nz) = fbuf(0:14,0:nx-1,0:ny-1,0:nz-1)
#endif
   if(master) then
      deallocate(fbuf)
      allocate(fbuf(0:14,0:nx-1,0:ny-1,0:zlen(myid)+1)) 
   end if
   if(master) write(*,*) "Restored configuration from checkpoint file"
 
end subroutine

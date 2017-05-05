!************************************************************************
!*                                                                      *
!*     LBModule                                                         *
!*                                                                      *
!************************************************************************

! ... module containing global variables

module LBModule

#if defined (MPI)
   use mpi
#endif

   real(8), parameter :: cs2 = 1.0d0/3.0d0
   real(8), parameter :: Pi = 4.0d0*datan(1.0d0) 
   integer(4), parameter :: mnproc = 20             ! Maximum number of processes
   integer(4), parameter :: rootid = 0
#if defined (MPI)
   integer(4), parameter :: comm = mpi_comm_world
#endif
   
   integer(4) :: nproc                                ! Number of MPI processes
   integer(4) :: myid                                 ! Process ID
   integer(4) :: iplow(0:mnproc-1), ipupp(0:mnproc-1), iplen(0:mnproc-1)  ! Local upper and lower swimmer index
   integer(4) :: trlow(0:mnproc-1), trupp(0:mnproc-1), trlen(0:mnproc-1)  ! Local upper and lower tracer index
   integer(4) :: zlow(0:mnproc-1), zupp(0:mnproc-1), zlen(0:mnproc-1)     ! Local upper and lower lattice limits
   integer(4) :: cpudown, cpuup                       ! Neighbours in the processor array

   integer(4) :: iseed = 243453                   ! Seed for random number generator
   logical    :: master 

   real(8) :: w(0:14)                             ! Lattice weights
   integer(4) :: ci(3,0:14)                       ! Velocity vectors
   integer(4) :: nx,ny,nz                         ! Side lengths
   integer(4) :: nstep                            ! Number of LB steps
   integer(4) :: idump                            ! Interval between dumping configs 
   integer(4) :: startstep                        ! Initial timestep (if restored from checkpoint)
 
   real(8) :: eta, tau, omega                     ! Viscosity etc

   integer(4) :: nSwim                            ! Number of swimmers
   integer(4) :: nTrac                            ! Number of swimmers
   real(8)    :: vswim                            ! Swim speed
   real(8)    :: fswim                            ! Swimmer force
   real(8)    :: l                                ! Swimmer length
!  real(8)    :: a                                ! Swimmer body radius
!  real(8)    :: h                                ! Correction to hydrodynamic radius
   real(8)    :: lambda                           ! Run length
   real(8)    :: tumbleProb                       ! Tumble rate = vswim/lambda

   logical    :: ltumbles = .true.                ! Flag to control tumbling
   logical    :: lswims = .true.                  ! Self-propulsion on or off
   logical    :: ladvects = .true.                ! Advection by the fluid on or off
   logical    :: lrotates = .true.                ! Rotation by the fluid on or off
   logical    :: ldumpswim = .false. 
   logical    :: ldumplat = .false.               ! Flags to control data dumping
   logical    :: lformatted = .false.             ! Formatted or unformatted output?
   logical    :: lrestore = .false.

   real(8), allocatable :: f(:,:,:,:)         ! Distribution functions
   real(8), allocatable :: force(:,:,:,:)     ! Forces
   real(8), allocatable :: u(:,:,:,:)         ! Fluid velocity
   real(8), allocatable :: rho(:,:,:)         ! Density
   real(8), allocatable :: help(:,:,:,:)      ! Temporary array (should hold f)
   real(8), allocatable :: rbuf(:,:)          ! Temporary array (should hold r)

   real(8), allocatable :: r(:,:)         ! swimmer positions
   real(8), allocatable :: n(:,:)         ! swimmer orientations
   real(8), allocatable :: rtr(:,:)       ! tracer positions
   real(8), allocatable :: str(:,:)       ! tracer positions (unwrapped)

end module LBModule


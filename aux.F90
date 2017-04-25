
!************************************************************************
!*                                                                      *
!*     Peskin                                                           *
!*                                                                      *
!************************************************************************
 
! ... Peskin's regularisation of the Dirac delta function
 
function Peskin(x)

   implicit none

   real(8), intent(in)  :: x   
   real(8)  :: Peskin
   real(8) :: abs_x, root, phi

   abs_x = abs(x)
   
   if(abs_x >= 2.0d0) then
      Peskin = 0.0d0
   else if(abs_x >= 1.0d0) then
      root = -4.0d0*x*x + 12.0d0*abs_x - 7.0d0  
      phi = -2.0d0*abs_x + 5.0d0 - dsqrt(root)
      Peskin = 0.125d0*phi
   else 
      root = -4.0d0*x*x + 4.0d0*abs_x + 1.0d0
      phi = -2.0d0*abs_x + 3.0d0 + dsqrt(root)
      Peskin = 0.125d0*phi
   end if 

end function
 

!************************************************************************
!*                                                                      *
!*     Random                                                           *
!*                                                                      *
!************************************************************************
 
! ... return a random number in the range of 0 to 1
!     by Learmonth and Lewis 1973 (32 bits integer)
 
function Random(iseed)

   implicit none
   integer(4), parameter :: k = 16087, l = 0, nb = 31
   real(8), parameter :: f = 2.0d0**(-nb)
   integer(4), intent(inout) :: iseed    ! iseed should not be equal to Zero
   real(8) :: Random

   iseed = iseed*k+l
   iseed = ibclr(iseed,nb)
   Random = f*iseed

end function Random

!> Module containing all globally defined parameters
module param_mod
implicit none
complex, parameter :: i=(0.0,1.0)
real, parameter :: pi= 4*atan(1.0)
integer :: Nspecies
integer :: sign
real :: theta
real, parameter :: zero=  0.0, half=  0.5, one=  1.0, &
                       r3o2=  1.5, two=   2.0, r5o2= 2.5, &
                       three= 3.0, r7o2=  3.5, four= 4.0, &
                       r9o2=  4.5, five=  5.0, six=  6.0, &
                       seven= 7.0, eight= 8.0, nine= 9.0, &
                       ten=   1.0

complex :: z0= (zero, zero), z1= (one, zero), zi= (zero, one)
integer, allocatable, dimension (:) :: mode
real, allocatable, dimension (:) :: mu
real, allocatable, dimension (:) :: q
real, allocatable, dimension (:) :: dens
real, allocatable, dimension (:) :: drift
real, allocatable, dimension (:) :: beta_para
real, allocatable, dimension (:):: beta_perp
real, allocatable, dimension (:) :: beta_ratio
real, allocatable, dimension (:,:,:) :: distribution
real, allocatable, dimension (:,:) :: vpara,vperp
integer, allocatable, dimension (:) :: npara, nperp
integer :: npara_max, nperp_max
integer :: narb
real :: delta 
real :: rf_error
real :: eps_error
! Blocos de interfaces
abstract interface
   pure function fabst(x)
   !import :: dp
   real :: fabst
   real, intent(in) :: x
   end function fabst
!***
   function fzabst(z)
   !import :: dp
   complex :: fzabst
   complex, intent(in) :: z
   end function fzabst
end interface
!
end module param_mod

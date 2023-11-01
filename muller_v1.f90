subroutine muller_rg(logu, omega_start, k, sol, splcoeff1, splcoeff2)
use param_mod ! Importe os parâmetros do módulo param_mod
implicit none

integer, intent(in) :: logu
real(8), intent(in) :: k
complex(16), intent(in) :: omega_start
complex(16), intent(out) :: sol
complex(16) :: disp_det
complex(16), dimension(1:4) :: fx
complex(16), dimension(1:4) :: omega
complex(16) :: a, b, c, d1, d2
integer :: n, j

complex(16) :: om, ga, om_old, ga_old
real(8), parameter :: rf_error = 1.0e-12 ! Defina rf_error conforme necessário

! Muller's method requires three starting points
! One point is chosen by the user - the other two points are taken to be slightly left and right of it 
omega(1) = 0.999D0 * omega_start
omega(2) = omega_start
omega(3) = 1.001D0 * omega_start

fx(1) = disp_det(omega(1), k, splcoeff1, splcoeff2)
fx(2) = disp_det(omega(2), k, splcoeff1, splcoeff2)
fx(3) = disp_det(omega(3), k, splcoeff1, splcoeff2)

! Perform Muller iteration
n = 0

do while (.true.)

   ! Determine coefficients from preceding three points

   a = ((omega(2) - omega(3)) * (fx(1) - fx(3)) - (omega(1) - omega(3)) * (fx(2) - fx(3))) / &
       & ((omega(1) - omega(3)) * (omega(2) - omega(3)) * (omega(1) - omega(2)))
   b = ((omega(1) - omega(3))**2 * (fx(2) - fx(3)) - (omega(2) - omega(3))**2 * (fx(1) - fx(3))) / &
       & ((omega(1) - omega(3)) * (omega(2) - omega(3)) * (omega(1) - omega(2)))
   c = fx(3)

   d1 = b + sqrt(b**2 - 4.0D0 * a * c)
   d2 = b - sqrt(b**2 - 4.0D0 * a * c)

   ! Compute new root from coefficients

   if (abs(d1) >= abs(d2)) then
      omega(4) = omega(3) - 2.0D0 * c / d1
   else
      omega(4) = omega(3) - 2.0D0 * c / d2
   endif

   fx(4) = disp_det(omega(4), k, splcoeff1, splcoeff2)

   ! Measure the accuracy of iterated root and check exit-condition

   om = real(omega(4))
   ga = aimag(omega(4))

   om_old = real(omega(3))
   ga_old = aimag(omega(3))

   if (((((om >= om_old) .and. (abs(1.0D0 - abs(om_old / om)) < rf_error)) .or. &
        & ((om < om_old) .and. (abs(1.0D0 - abs(om / om_old)) < rf_error))) .and. &
       & (((ga >= ga_old) .and. (abs(1.0D0 - abs(ga_old / ga)) < rf_error)) .or. &
       & ((ga < ga_old) .and. (abs(1.0D0 - abs(ga / ga_old)) < rf_error)))) .or. &
      & (((((om >= om_old) .and. (abs(1.0D0 - abs(om_old / om)) < rf_error)) .or. &
        & ((om < om_old) .and. (abs(1.0D0 - abs(om / om_old)) < rf_error))) .and. &
       & ((abs(ga) < 1.0D-10) .and. (abs(ga_old) < 1.0D-10)))) .or. &
      & (((((ga >= ga_old) .and. (abs(1.0D0 - abs(ga_old / ga)) < rf_error)) .or. &
        & ((ga < ga_old) .and. (abs(1.0D0 - abs(ga / ga_old)) < rf_error))) .and. &
       & ((abs(om) < 1.0D-10) .and. (abs(om_old) < 1.0D-10)))) .or. &
      & ((abs(om) < 1.0D-10) .and. (abs(ga) < 1.0D-10))) then
      exit ! Exit the loop if convergence is achieved
   endif

   ! Stop iteration if last step was ineffective
   if ((abs((real(fx(4)) - real(fx(3))) / real(fx(4))) < 1.0D-12) .and. &
      & (abs((aimag(fx(4)) - aimag(fx(3))) / aimag(fx(4))) < 1.0D-12)) then
      write(logu, *) 'Last step in Muller iteration was ineffective'
      exit
   endif

   do j = 1, 3
      omega(j) = omega(j+1)
      fx(j) = fx(j+1)
   end do

   if (n > 40) then
      write(logu, *) 'Error: Muller method did not converge'
      exit
   endif

   n = n + 1

end do

! Solution of root finding procedure
sol = omega(4)

end subroutine muller_rg

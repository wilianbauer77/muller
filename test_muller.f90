program test_muller
  implicit none
  integer, parameter :: nnov = 1
  integer, parameter :: nprev = 2
  integer, parameter :: maxit = 100
  real, parameter :: errabs = 1.0e-6
  real, parameter :: errrel = 1.0e-6
  real, parameter :: fator = 1.0
  complex, dimension(nnov + nprev) :: czeros
  real, dimension(nprev + 1) :: iterac
  integer :: i

  ! Função cúbica: f(x) = x^3 - 2x^2 + x - 1
  interface
    function cubic_fn(x)
      complex :: cubic_fn
      real, intent(in) :: x
    end function cubic_fn
  end interface

  ! Definição da função cúbica
  function cubic_fn(x)
    real, intent(in) :: x
    complex :: cubic_fn
    cubic_fn = cmplx(x**3 - 2.0*x**2 + x - 1.0, 0.0)
  end function cubic_fn

  ! Inicialização dos valores iniciais das raízes
  czeros = [(0.5, 0.0), (1.0, 0.0), (2.0, 0.0)]

  ! Chamada da rotina Muller
  call muller(cubic_fn, nnov, nprev, maxit, errabs, errrel, czeros, fator, iterac)

  ! Impressão dos resultados
  do i = nprev + 1, nnov + nprev
    print "('Root ', I0, ': (', F8.6, ', ', F8.6, ')')", i - nprev, real(czeros(i)), aimag(czeros(i))
    if (present(iterac)) then
      print "('Number of iterations for root ', I0, ': ', I0)", i - nprev, iterac(i - nprev)
    end if
  end do

end program test_muller

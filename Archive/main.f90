program Main
  use Metodos_Computacionais_Fisica, only: muller
  !use muller

  implicit none

  ! Declarando as variáveis necessárias
  integer, parameter :: nnov = 3        ! Número de variáveis do polinômio
  integer, parameter :: nprev = 0       ! Número de raízes conhecidas (opcional)
  integer, parameter :: maxit = 100     ! Número máximo de iteraçsões
  real, parameter :: errabs = 0     ! Erro absoluto
  real, parameter :: errrel = 1.0d-4 !Erro relativo
  complex, dimension(nnov+nprev) :: czeros  ! Vetor para armazenar as raízes
  integer, dimension(nnov) :: iterac    ! Vetor para armazenar o número de iterações (opcional)
  integer :: i !Variável  auxiliar
   czeros = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]

! PRINT*, SIZE(CZEROS)
! CZEROS= (0.0, 0.0)
! PRINT*, CZEROS


  ! Chamando a subrotina "muller"
! PRINT*, 'ENTRA MULLER'
  call muller(fn, nnov, nprev, maxit, errabs, errrel, czeros)
! PRINT*, 'SAIU MULLER'

  ! Imprimindo as raízes na tela
  write(*,*) "Raízes encontradas:"
  do i = 1, nnov
    write(*,*) "Raiz ", i, ": ", czeros(i)!, "errel: ", errrel
  end do


contains

  ! Implementando a função "fn"
  function fn(x)
    complex, intent(in) :: x
    complex :: fn

    ! Inserindo a expressão do polinômio
    fn = x**3 - x  - 1

    return
  end function fn

end program Main


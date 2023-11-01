!********************** MÓDULO METODOS_COMPUTACIONAIS_FISICA ************************
! Contém constantes nomeadas e interfaces das rotinas incluídas na Apostila:
!                         Introdução à Física Computacional                         
! Ref: http://professor.ufrgs.br/rgaelzer/pages/comp-phys
! Autor: Rudi Gaelzer, IF-UFRGS
! Data: Outubro/2019
MODULE Metodos_Computacionais_Fisica
use, intrinsic :: iso_fortran_env
implicit none
real, parameter :: zero=  0.0, half=  0.5, one=  1.0, &
                       r3o2=  1.5, two=   2.0, r5o2= 2.5, &
                       three= 3.0, r7o2=  3.5, four= 4.0, &
                       r9o2=  4.5, five=  5.0, six=  6.0, &
                       seven= 7.0, eight= 8.0, nine= 9.0, &
                       ten=   10.00
real, parameter :: pi=     3.14159265358979323846264338327950288419717
real, parameter :: rtpio2= 0.88622692545275801364908374167057259139877 ! raiz(pi)/2
real, parameter :: pid2=   1.57079632679489661923132169163975144209858 ! pi/2
complex, parameter :: z0= (zero, zero), z1= (one, zero), zi= (zero, one)
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
!***
   subroutine fdfabst(x, fx, dfdx)
   !import :: dp
   real, intent(in) :: x
   real, intent(out) :: fx, dfdx
   end subroutine fdfabst
!***
   subroutine vfdfabst(x, fx, dfdx)
   !import :: dp
   real, intent(in) :: x
   real, dimension(:), intent(in)  :: fx
   real, dimension(:), intent(out) :: dfdx
   end subroutine vfdfabst
end interface
!
interface
   subroutine dfdx_rich(f, x, h_ini, errest, dfdx, err_sai)
   !import :: dp, fabst
    import :: fabst
   real, intent(in)  :: x, h_ini, errest
   real, intent(out) :: dfdx, err_sai
   procedure(fabst) :: f
   end subroutine dfdx_rich
!***
   module function simpson(f, n, a, b)
   real             :: simpson
   integer, intent(in)  :: n
   real, intent(in) ::  a, b
   procedure(fabst)     :: f
   end function simpson
!***
   module function trapez(f, n, a, b)
   real             :: trapez
   integer, intent(in)  :: n
   real, intent(in) ::  a, b
   procedure(fabst)     :: f
   end function trapez
!***
   function trapez_rom(f, a, b, n_ordem)
   !import :: dp, fabst
   import :: fabst
   real             :: trapez_rom
   integer, intent(in)  :: n_ordem
   real, intent(in) :: a, b
   procedure(fabst)     :: f
   end function trapez_rom
!***
   function pntmed_rom(f, a, b, n_ordem)
   !import :: dp, fabst
   import :: fabst
   real             :: pntmed_rom
   integer, intent(in)  :: n_ordem
   real, intent(in) :: a, b
   procedure(fabst)     :: f
   end function pntmed_rom
!***
   function pntmed_inf_rom(f, a, n_ordem)
   import :: fabst
   !import :: dp, fabst
   real             :: pntmed_inf_rom
   integer, intent(in)  :: n_ordem
   real, intent(in) :: a
   procedure(fabst)     :: f
   end function pntmed_inf_rom
!***
   subroutine muller (fn, nnov, nprev, maxit, errabs, &
                     errrel, czeros, fator, fnreal, iterac)
   import :: fzabst
   !import ::  dp, fzabst 
   logical, optional, intent(in)                     :: fnreal
   integer, intent(in)                               :: nnov, nprev, maxit
   integer, dimension(nnov), optional, intent(out)   :: iterac
   real, intent(in)                              :: errabs, errrel
   real, optional, intent(in)                    :: fator
   complex, dimension(nnov+nprev), intent(inout) :: czeros
   procedure(fzabst)                                 :: fn
   end subroutine muller
!***
   subroutine rk4(x, y, h, ysai, derivs)
   !import :: dp, vfdfabst
   import :: vfdfabst
   real, intent(in)                :: x, h
   real, dimension(:), intent(in)  :: y
   real, dimension(:), intent(out) :: ysai
   procedure(vfdfabst)                  :: derivs
   end subroutine rk4
end interface
! Rotinas auxiliares
interface
   module subroutine troca(a, b)
   real, intent(inout) :: a, b
   end subroutine troca
!***
   module subroutine verifica_tamanhos(vint, string)
   character(len=*), intent(in) :: string
   integer, dimension(:), intent(in) :: vint
   end subroutine verifica_tamanhos
!***
   module subroutine xw_Laguerre(n, x_s, w_s)
   integer, intent(in)                       :: n ! 2 <= n <= n_max.
   real(kind= 8), dimension(n), intent(out) :: x_s, w_s
   end subroutine xw_Laguerre
end interface
END MODULE Metodos_Computacionais_Fisica

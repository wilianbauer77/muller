!******************************** SUBMÓDULO MCF_SMOD ********************************
! Contém constantes as definições das interfaces no módulo:
!                             Metodos_Computacionais_Fisica
! Autor: Rudi Gaelzer, IF-UFRGS
! Data: Outubro/2019
SUBMODULE(Metodos_Computacionais_Fisica) MCF_Smod
CONTAINS
   include "simpson.f90"
   include "trapez.f90"
   include "xw_Laguerre.f90"
! Rotinas auxiliares
   include "troca.f90"
   include "verifica_tamanhos.f90"
END SUBMODULE MCF_Smod

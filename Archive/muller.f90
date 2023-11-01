!**************************** SUBROTINA MULLER ******************************
! Encontra as raízes de uma função analítica unívoca pelo Método de Muller.
! Argumentos:
!   fn: Função analítica unívoca f(z) cujas raízes são procuradas (Entrada)
!   nnov: Número total de raízes novas a serem encontradas.       (Entrada)
!   nprev: Número de raízes previamente conhecidas.               (Entrada)
!   maxit: Número maximo de chamadas da função fn(z) por raiz.    (Entrada)
!   errabs: Primeiro critério de parada.                          (Entrada)
!           Iterações são interrompidas se abs(fn(z)) .lt. errabs.
!   errrel: Segundo critério de parada.                           (Entrada)
!           Iterações são interrompidas se abs(h) .lt errrel*abs(z).
!   czeros: Vetor que contém as raízes de fn(z).             (Entrada/Saída)
!       czeros(1), ..., czeros(nprev): raízes previamente conhecidas.
!       czeros(nprev+1), ..., czeros(n): (ent.) aproximações iniciais 
!                                             para raízes.
!                                      (sai.) raízes encontradas.
!   fator: (OPCIONAL) Fator multiplicativo para valores iniciais. (Entrada)
!          Para z= czeros(i), os 3 primeiros pontos para ajustar a parábola
!          são z, z + fator*h e z - fator*h, sendo h um valor fixo.
!   fnreal: (OPCIONAL) Variável lógica para raízes reais.         (Entrada)
!           fnreal = .true. se todas as raízes são reais.
!           fnreal= .false. se há raízes complexas (valor padrão).
!   iterac: (OPCIONAL) Vetor com o número de iterações realizadas
!           para cada raiz.                                       (saída)
! Autor: Rudi Gaelzer, IF-UFRGS.
! Data: Abril/2011 (versão 2).
subroutine muller (fn, nnov, nprev, maxit, errabs, &
                   errrel, czeros, fator, fnreal, iterac)
use Metodos_Computacionais_Fisica, only: zero, half, one, two, &
                                         four, ten, z0, fzabst 
implicit none
! Variaveis mudas.
logical, optional, intent(in)                     :: fnreal
integer, intent(in)                               :: nnov, nprev, maxit
integer, dimension(nnov), optional, intent(out)   :: iterac
real, intent(in)                              :: errabs, errrel
real, optional, intent(in)                    :: fator
complex, dimension(nnov+nprev), intent(inout) :: czeros
procedure(fzabst)                                 :: fn
! Variaveis locais.
logical           :: tes_maxit, tes_g_mxit, tes_real, tes_iter
integer           :: ntot, i, n_it
real    :: eps1, eps2, teste, h_ini= 0.1, tiny1= tiny(one)
real, parameter    :: big1= huge(one)
complex, parameter :: z1o10=(0.1,zero), z2=(two,zero), z4=(four,zero)
complex :: c, den, divdf1, divdf2, dvdf1p, fzr, z_inc, zrprev
complex :: fzrdfl, fzrprv, h, hprev, czero, sqr, him2, him1

! PRINT*, 'ENTROU MULLER'
! PRINT*, NNOV, NPREV, MAXIT
! PRINT*, ERRABS, ERRREL
! PRINT*, CZEROS
! 
! PRINT*, FN(CZEROS(1)), FN(CMPLX(2.0))
! 
! PRINT*, PRESENT(FNREAL), PRESENT(FATOR), PRESENT(ITERAC)
! PRINT*, FATOR
! PRINT*, ITERAC

if(nnov < 1)then
   print*, 'O argumento nnov deve ser >= 1.'
   stop
end if
!                   Inicializações.
tes_real= .false. ; tes_iter= .false. ; tes_g_mxit= .false.

! PRINT*, TES_REAL, TES_ITER, TES_G_MXIT
! PRINT*, 'VAI IMPRIMIR FNREAL'
! PRINT*, FNREAL

if(present(fnreal))tes_real= fnreal
if(present(fator))h_ini= fator*h_ini
if(present(iterac))tes_iter= .true.

! PRINT*, TES_REAL, H_INI, TES_ITER

if(big1*tiny1 < one)tiny1= one/big1
eps1 = max(errrel, ten*epsilon(one)) ! Erro de primeira espécie.
eps2 = max(errabs, ten*tiny1)        ! Erro de segunda espécie.
ntot= nnov + nprev


! PRINT*, EPS1, EPS2, NTOT

!
l_roots: do  i = nprev + 1, ntot
   n_it = 0
   tes_maxit= .false.
!                   Calcule os três primeiros valores da i-ésima raiz como
!                       czeros(i) + h, czeros(i)- h, czeros(i)
   czero = czeros(i)
   h = z1o10*cmplx(h_ini, kind = 8)
   if(abs(czero) > h_ini) h= z1o10*czero
   him2= h
   call dflac(czero, him2, i, fzr, dvdf1p)  ! f(czero + h)
   him1= -h
   call dflac(czero, him1, i, fzr, fzrprv)  ! f(czero - h)
   hprev = him1 - him2
   zrprev= czero + him1
   dvdf1p = (fzrprv - dvdf1p)/hprev
   l_iter: do
      l_div: do
         z_inc= z0
         call dflac(czero, z_inc, i, fzr, fzrdfl)
         czero= czero + z_inc
         if (tes_maxit)then
            tes_g_mxit= .true.
            exit l_iter
         end if
         h= czero - zrprev
!                           Testa convergência da segunda espécie.
         if (max(abs(fzr),abs(fzrdfl)) < eps2)exit l_iter
!                           Testa convergência da primeira espécie.
         teste= abs(h) - eps1*abs(czero)
         if (teste < zero)exit l_iter
!                           Verifique se valor iterado diverge da raiz.
         if (abs(fzrdfl) < ten*abs(fzrprv))exit l_div
         h = cmplx(half, kind = 8)*h
         czero = czero - h
      end do l_div
!                           Inicia algoritmo principal.
      divdf1 = (fzrdfl - fzrprv)/h
      divdf2 = (divdf1 - dvdf1p)/(h + hprev)
      hprev = h ; zrprev = czero
      dvdf1p = divdf1
      c = divdf1 + h*divdf2
      sqr = c*c - z4*fzrdfl*divdf2
      if (tes_real .and. (real(sqr) < zero)) sqr = z0
      sqr = sqrt(sqr)
      teste= sign(one, real(c)*real(sqr)+aimag(c)*aimag(sqr))
      den = c + cmplx(teste, kind = 8)*sqr
      h = -z2*fzrdfl/den
      fzrprv = fzrdfl
      czero = czero + h
   end do l_iter
   czeros(i) = czero
   if(tes_iter)iterac(i - nprev)= n_it
end do l_roots
if(tes_g_mxit)print*, 'Método nao convergiu em', maxit, &
                      'passos para 1 ou mais raízes.'
return
CONTAINS
   subroutine dflac(zero, z_inc, i, fzero, fzrdfl)
   integer, intent(in)              :: i
   complex, intent(in)    :: zero
   complex, intent(out)   :: fzero, fzrdfl
   complex, intent(inout) :: z_inc
   logical           :: t_den= .true.
   integer           :: j
   real, dimension(2) :: v_sft
   complex :: root, den
 !
   l_den: do
      if (n_it == maxit)then
         tes_maxit= .true.
         return
      end if
      n_it = n_it + 1
      root= zero + z_inc
      fzero = fn(root)
      fzrdfl = fzero
      l_deflac: do j = 2, i
         den = root - czeros(j-1)
! Teste para evitar singularidade ou overflow.
         if ((abs(den) < tiny1) .or. &
             (big1*min(abs(den),one) <= abs(fzrdfl))) then
! Desloca o ponto aleatoriamente.
            call random_number(v_sft)
            v_sft= ten*(v_sft - half)*eps1
            z_inc= z_inc + cmplx(v_sft(1), v_sft(2), kind= 8)
            t_den= .false.
            exit l_deflac
         else
            fzrdfl = fzrdfl/den
            t_den= .true.
         end if
      end do l_deflac
      if(t_den) exit l_den
   end do l_den
   return
   end subroutine dflac
end subroutine muller

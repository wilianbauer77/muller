SUBROUTINE MULLER(fn,roots,n,maxiter,epsilon1,epsilon2)
!fn: nome da subrotina (função) de forma fn(x,fx). Devolve um complexo.
!fx: função para um x complexo (entrada)
!n: número de raízes procuradas
!roots(1), ..., roots(n): primeira estimativa das n raízes. Se você não conhece nenhuma, 0 éuma boa escolha.
!maxiter: número máximo de iterações permitidas por raiz (entrada)
!epsilon1: tamanho máximo de h permitido
!epsilon2: tamanho máximo de |f(x)| se o critério de epsilon1 não for atendido.
!roots(1), ..., roots(n): n raízes calculadas

COMPLEX :: roots(n), x, fx, cprim, cprev, c, a, b, hprev, h,  dd, ddp, sqr, den, den1
EXTERNAL fn !disp_det
DO i = 1,n
   kount = 0;
   x=roots(i)
   CALL DEFLATE(fn,x+0.5,i,kount,fx,cprim,roots)
   CALL DEFLATE(fn,x-0.5,i,kount,fx,cprev,roots)
   CALL DEFLATE(fn,x,i,kount,fx,c,roots)
   hprev=-1.0
   h=0.5
   ddp=(cprev-cprim)/hprev
   10 dd=(c-cprev)/h
   a=(dd-ddp)/(hprev+h)
   b=a*h+dd
   sqr=CSQRT(b**2-4.0*a*c)
   den=b+sqr
   den1=b-sqr
   IF(CABS(den1)>CABS(den)) den=den1
   IF(CABS(den)<=0.) den=1.0
   h = -2*c/den
   x=x+h
   cprev=c
   hprev=h
   IF(kount>maxiter) THEN
      PRINT*, 'Máximo de iterações alcançado'
      RETURN
   ENDIF
   20 CALL DEFLATE(fn,x,i,kount,fx,c,roots)
   IF(CABS(h)<epsilon1*CABS(x)) GOTO 30
   IF(MAX(CABS(fx),CABS(c))<epsilon2) GOTO 30
   IF (CABS(c)>10*CABS(cprev)) THEN
      h=0.5*h; x=x-h; GOTO 20
      ELSE
      GOTO 10
   END IF
   30 roots(i)=x
ENDDO
RETURN
END SUBROUTINE MULLER

!*****************************************************

SUBROUTINE DEFLATE(fn,x,i,kount,fx,fxdfl,roots)
COMPLEX :: x, fx, fxdfl, den, roots(i)
EXTERNAL fn
kount = kount+1
CALL fn(x,fx)
fxdfl=fx
IF(i==1) RETURN
DO j=2,i
   den=x-roots(j-1)
   IF(CABS(den)==0.0) THEN
      roots(i)=x*1.001
      RETURN
   ELSE
      fxdfl=fxdfl/den
   END IF
END DO
RETURN
END SUBROUTINE DEFLATE

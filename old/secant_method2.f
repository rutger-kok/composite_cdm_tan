      PROGRAM SECANT_METHOD
 
        IMPLICIT NONE
        REAL*8 rtsec, f, root
        
        root = rtsec(f,-3,5,0.0001)
        PRINT*, 'Root=', root
        PAUSE "Press RETURN to end program"
      END PROGRAM SECANT_METHOD

      FUNCTION f(x)
        IMPLICIT NONE
        REAL*8 f,x
        f = x*x*x - 3.0*x*x + 2.0*x
        RETURN
        END


      FUNCTION rtsec(func,x1,x2,xacc)
        INTEGER MAXIT
        REAL rtsec,x1,x2,xacc,func
        EXTERNAL func
        PARAMETER (MAXIT=30) !Maximum allowed number of iterations.
C        Using the secant method, find the root of a function func thought to lie between x1 and
C        x2. The root, returned as rtsec, is refined until its accuracy is ±xacc.
        INTEGER j
        REAL dx,f,fl,swap,xl
        fl=func(x1)
        f=func(x2)
        if(abs(fl).lt.abs(f))then !Pick the bound with the smaller function value as the most recent guess.
          rtsec=x1
          xl=x2
          swap=fl
          fl=f
          f=swap
        else
          xl=x1
          rtsec=x2
        endif
        
        do j=1,MAXIT !Secant loop.
          dx=(xl-rtsec)*f/(f-fl) !Increment with respect to latest value.
          xl=rtsec
          fl=f
          rtsec=rtsec+dx
          f=func(rtsec)
          if(abs(dx).lt.xacc.or.f.eq.0.) return !Convergence.
        enddo
        pause 'rtsec exceed maximum iterations'
        END
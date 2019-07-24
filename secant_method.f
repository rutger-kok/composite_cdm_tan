      PROGRAM SECANT_METHOD
 
        IMPLICIT NONE
        INTEGER i,limit
        REAL*8 d, e, f, x, x1, x2
        
        i = 1
        limit = 100
         
        x1 = -1.0
        x2 = 3.0
        e = 1.d-15
         
        DO 
          IF (i > limit) THEN
            PRINT*,"Function not converging"
            EXIT
          END IF
          d = (x2 - x1) / (f(x2) - f(x1)) * f(x2)
          IF (ABS(d) < e) THEN
            PRINT*,"Root found at x = ", x2
            EXIT    
          END IF
          x1 = x2
          x2 = x2 - d
          i = i + 1
        END DO

        PAUSE "Press RETURN to end program"
      END PROGRAM SECANT_METHOD

      FUNCTION f(x)
        IMPLICIT NONE
        REAL*8 f,x
        f = x*x*x - 3.0*x*x + 2.0*x
        RETURN
        END
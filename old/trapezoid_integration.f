    
      PROGRAM TRAPZD_TEST
        IMPLICIT NONE
        INTEGER k
        REAL*8 f, Y(200), X(200), result

        DO k = 1,200
          X(k) = (k-1)*0.1
          Y(k) = f(X(k))
        ENDDO

        call trapzd(X,Y,result)
        
        PRINT*, 'Result = ', result
        PAUSE "Press RETURN to end program"
        END PROGRAM TRAPZD_TEST
      
      SUBROUTINE trapzd(xarray,yarray,result)
        IMPLICIT NONE        
C       Input variables
        REAL*8 xarray(200),yarray(200)
C       Local variables
        INTEGER i
        REAL*8 s,x1,x2,y1,y2,delta
C       Output variables
        REAL*8 result

        result = 0.d0
        DO i = 1,199
            x1 = xarray(i)
            x2 = xarray(i+1)
            y1 = yarray(i)
            y2 = yarray(i+1)
            delta = x2 - x1
            s=0.5*delta*(y1+y2)
            result = result + s
        ENDDO
        RETURN
      END

      FUNCTION f(x)
        IMPLICIT NONE
        REAL*8 f,x
        f = x*x*x - 3.0*x*x + 2.0*x
        RETURN
        END
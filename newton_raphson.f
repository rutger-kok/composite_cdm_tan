        PROGRAM newton_raphson

        implicit none
        REAL NEWT, FUNC, DFUNC
        REAL  x, result
        PRINT*,'GIVE AN INITIAL GUESS FOR THE ROOT OF THE FUNCTION F(X)'
        READ*, x ! initial estimate
        result = NEWT(x)
        PRINT*,'Final X= ', result
        PAUSE "Press RETURN to end program"
        END
        
        FUNCTION NEWT(init)
        REAL NEWT, xOld, xNew, tol, init, df,dx,f
        xOld = init
        xNew = init
        dx = 100
        tol = 0.0001
10      IF (abs(dx).gt.tol) THEN
          xOld = xNew
          dx=FUNC(xOld)/DFUNC(xOld)
          xNew = xOld-dx
          GOTO 10
        ELSE
          NEWT = xNew
        END IF
        RETURN
        END

        FUNCTION FUNC(xOld)
        REAL FUNC
        FUNC = xOld**3 - 4.0d0
        RETURN
        END

        FUNCTION DFUNC(xOld)
        REAL DFUNC
        DFUNC = 3*xOld**2
        RETURN
        END


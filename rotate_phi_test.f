      PROGRAM rotate_phi_test
      IMPLICIT NONE

      REAL*8 trialStress(6), phiC, XC, NEWT, FUNC, DFUNC
      REAL*8 trialStressP(6)

      XC = 1.2001
      phiC = 0.116264127346
      trialStress(1) = -1.2001
      trialStress(2) = 0.0
      trialStress(3) = 0.0
      trialStress(4) = 0.0
      trialStress(5) = 0.0
      trialStress(6) = 0.0
    
      call ROTATE_PHI(trialStress, phiC, XC, trialStressP)

      PRINT*,'Trial Stress = ', trialStressP
      PAUSE "Press RETURN to end program"
      END

* <<<<<<<<<<<<<<<<<<<<<< SUBROUTINE ROTATE_PHI>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* ROTATION OF STRESSES TO THE MISALIGNMENT COORDINATE FRAME *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE ROTATE_PHI(trialStress, phiC, XC, trialStressP)
        IMPLICIT NONE
        ! input variables
        REAL*8 trialStress(6), phiC, XC, NEWT
        ! local variables
        REAL*8 theta, phi, X, m, n, u, v
        REAL*8 trialStressT(6)
        ! output variables
        REAL*8 trialStressP(6)

        ! first determine fracture plane angle theta (fiber kinking)   
        IF ((trialStress(4).eq.0.d0).AND.(trialStress(6).eq.0.d0)) THEN
            IF ((trialStress(2)-trialStress(3)).eq.0.d0) THEN
                theta = atan(1.0d0) ! pi/4
            ELSE
                theta = 0.5d0*atan((2.0d0*trialStress(5))
     1                  /(trialStress(2)-trialStress(3))) !Eq.55 CLN
            ENDIF 
        ELSE
            IF (trialStress(4).eq.0.d0) THEN
                theta = 2.0d0*atan(1.0d0) ! pi/2
            ELSE 
                theta = atan(trialStress(6)/trialStress(4)) !Eq. 56 CLN
            ENDIF
        END IF
        PRINT*,'Theta = ', theta
        ! determine the misalignment angle phi
        X = (sin(2.0d0*phiC)*XC)/(2.0d0*phiC) ! Eq. 86 CLN

        IF (trialStress(4).gt.0) THEN
            phi = NEWT(0.1, trialStress, X) ! initial value of 0.1
        ELSE 
            phi = -1.0d0*NEWT(0.1, trialStress, X)
        END IF

        PRINT*,'Phi = ', phi

        m = cos(theta)
        n = sin(theta)
        ! Rotate stresses by angle theta
        trialStressT(1) = trialStress(1)
        trialStressT(2) = trialStress(2)*m**2 
     1                    + 2.0d0*trialStress(5)*m*n 
     2                    + trialStress(3)*n**2
        trialStressT(3) = trialStress(3)*m**2 
     1                    - 2.0d0*trialStress(5)*m*n 
     2                    + trialStress(2)*n**2
        trialStressT(4) = trialStress(4)*m + trialStress(6)*n
        trialStressT(5) = trialStress(5)*(m**2 - n**2)
     1                    - trialStress(2)*n*m 
     2                    + trialStress(3)*n*m
        trialStressT(6) = trialStress(6)*m - trialStress(4)*n

        ! Rotate stresses by angle phi
        u = cos(phi)
        v = sin(phi)
        trialStressP(1) = trialStressT(1)*u**2 
     1                    + 2.0d0*trialStressT(4)*u*v 
     2                    + trialStressT(2)*v**2
        trialStressP(2) = trialStressT(2)*u**2 
     1                    - 2.0d0*trialStressT(4)*v*u 
     2                    + trialStressT(1)*v**2     
        trialStressP(3) = trialStressT(3)      
        trialStressP(4) = trialStressT(4)*(u**2 -v**2) 
     1                    + trialStressT(2)*v*u 
     2                    - trialStressT(1)*v*u
        trialStressP(5) = trialStressT(5)*u - trialStressT(6)*v     
        trialStressP(6) = trialStressT(6)*u + trialStressT(5)*v

        PRINT*,'Trial Stress T = ', trialStressT(1)
        PRINT*,'Trial Stress P = ', trialStressP(1)

        RETURN
        END

* <<<<<<<<<<<<<<<<<<<<<<<<< FUNCTION NEWT >>>>>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Newton Raphson method to determine miaslignment angle phi
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      FUNCTION NEWT(init, trialStress, X)
        REAL*8 NEWT, xOld, xNew, tol, init, df,dx,f, trialStress(6), X
        xOld = init
        xNew = init
        dx = 100
        tol = 0.0001
10      IF (abs(dx).gt.tol) THEN
            xOld = xNew
            dx = FUNC(xOld, trialStress, X)/DFUNC(xOld, trialStress, X)
            xNew = xOld-dx
            GOTO 10
        ELSE
            NEWT = xNew
        END IF
        RETURN
        END

      ! Function used in the Newton Raphson iteration method
      FUNCTION FUNC(gammaOld, trialStress, X)
        REAL*8 FUNC, X, gammaOld, trialStress(6)
        ! Eq. 88 CLN
        FUNC = X*gammaOld + 0.5d0*(trialStress(1) 
     1         - trialStress(2))*sin(2.0d0*gammaOld)
     2         - abs(trialStress(4))*cos(2.0d0*gammaOld) 
      RETURN
      END

      ! Derivative function of FUNC
      FUNCTION DFUNC(gammaOld, trialStress, X)
        REAL DFUNC, X, gammaOld, trialStress(6)
        ! Eq. 89 CLN
        DFUNC = X + (trialStress(1)-trialStress(2))*cos(2.0d0*gammaOld)
     1          + 2.0d0*abs(trialStress(4))*sin(2.0d0*gammaOld) 
      RETURN
      END
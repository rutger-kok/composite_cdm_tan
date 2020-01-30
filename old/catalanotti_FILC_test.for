C----------------------------------------------------------------------C
C     ABAQUS VUMAT USER SUBROUTINE: catalanotti_FC.for                 C
C	    Author(s): Rutger Kok, Francisca Martinez-Hergueta               C
C     Date: 28/06/2019                                           	     C
C     Version 1.0                                                      C
C----------------------------------------------------------------------C

C User subroutine VUMAT
      PROGRAM catalanotti_FC_test
      IMPLICIT NONE
      REAL*8 trialStress(6),trialStrain(6)
      REAL*8 d1Plus,d1Minus,d2Plus,d2Minus,d6, stressPower
      REAL*8 rNew1,rNew2,rNew3,rNew4,rNew5,rNew6
      REAL*8 rOld1,rOld2,rOld3,rOld4,rOld5,rOld6
      REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
      REAL*8 nu21,nu31,nu32,delta
      REAL*8 d11,d22,d33,d12,d13,d23,d44,d55,d66
      REAL*8 XT,XC,YT,YC,SL,ST
      REAL*8 G1plus,G1minus,G2plus,G2minus,G6
      REAL*8 A1plus,A1minus,A2plus,A2minus,A6
      REAL*8 alpha0,etaT,etaL,phiC,omegaValue,kappa,lambda
      REAL*8 FI_LT,FI_LC,FI_MT,FI_MC,theta,phi
      REAL*8 initialcharLength,thickness,traceStrain,meanDamage,expo,A
      REAL*8 trialStressP(6), X, NEWT, FUNC, DFUNC, eps, maxFail
      REAL*8 xOmega1,xOmega2,xOmega3,xOmega4,xOmega5,xOmega6,aminDamage
      INTEGER k1, k2
      CHARACTER*80 cmname
      PARAMETER (eps=1.d-8)
      OPEN(1, file = 
     1 'C:\Workspace\catalanotti_failure_criteria\fi_lc_test.txt',
     2 status = 'unknown')
      
C Elastic constants orthotropic ply

      e11 = 161.0
      e22 = 11.4
      e33 = 11.4
      nu12 = 0.32
      nu13 = 0.32
      nu23 = 0.436
      g12 = 5.29
      g13 = 5.29
      g23 = 3.98

C Ply strength

      XT = 2.3235
      XC = 1.2001
      YT = 0.1602
      YC = 0.198
      SL = 0.1302

C Fracture Angle

      alpha0 = 0.9250245036 !53 degrees

C Fracture toughness

      G1plus = 0.1
      G1minus = 0.1
      G2plus = 0.00075
      G2minus = 0.0025
      G6 = 0.0035

C Initial values

      etaL = 0.5
      phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)*((SL/XC)+etaL))) !Eq.12
     1       /(2.0d0*((SL/XC)+etaL)))
      ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*SL)/  !Eq.12 (CLN)
     1     (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL))
      etaT = (etaL*ST)/SL                      !Eq.10 CLN
      kappa = (ST**2.0d0-YT**2.0)/(ST*YT)  !Eq.43 CLN
      lambda = ((2.0*etaL*ST)/SL)-kappa  !Eq.45 CLN 
      omegaValue = -0.120382            !Eq.20

C Stiffness matrix orthotropic material

      nu21 = nu12*(e22/e11)
      nu31 = nu13*(e33/e11)
      nu32 = nu23*(e33/e22)

      delta = (1.0-nu12*nu21-nu23*nu32-nu13*nu31-2.0*nu21*nu32*nu13)
     1        /(e11*e22*e33)

      d11 = (1-nu23*nu32)/(e22*e33*delta)
      d22 = (1-nu13*nu31)/(e11*e33*delta)
      d33 = (1-nu12*nu21)/(e11*e22*delta)
      d12 = (nu21+nu23*nu31)/(e22*e33*delta)
      d13 = (nu31+nu21*nu32)/(e22*e33*delta)
      d23 = (nu32+nu12*nu31)/(e11*e33*delta)
      d44 = 2.0*g12
      d55 = 2.0*g23
      d66 = 2.0*g13       

100   FORMAT (E,A,E,A,E,A,E,A,E)
      WRITE(1,*) 'trialStress1',',','FI_LT',',','FI_LC',',','theta',
     1           ',','phi'

      DO k1 = 1, 1000
          trialStress(1) = -1.250 + k1*0.0035236
          trialStress(2) = 0.0
          trialStress(3) = 0.0
          trialStress(4) = 0.0
          trialStress(5) = 0.0
          trialStress(6) = 0.0

          ! Evaluation of the damage activation functions

          ! longitudinal failure criteria
          IF (trialStress(1).gt.0.d0) THEN
            FI_LT = trialStress(1)/XT ! Eq. 54 CLN
          ELSEIF (trialStress(1).lt.0.d0) THEN
            call ROTATE_PHI(trialStress,phiC,XC,trialStressP,theta,phi)
            call FAIL_CLN(trialStressP,ST,YT,SL,etaL,etaT,lambda,kappa,
     1                    FI_LT,FI_LC)
          ENDIF
        
          WRITE(1,100) trialStress(1),',',FI_LT,',',FI_LC,',',theta,
     1                 ',',phi
          
      ENDDO
      END

* <<<<<<<<<<<<<<<<<<<<<< SUBROUTINE FAIL_CLN >>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Catalanotti failure criteria
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE FAIL_CLN(trialStress,ST,YT,SL,etaL,etaT,lambda,kappa,
     1                    FI_MT,FI_MC)
        IMPLICIT NONE
        ! input variables
        REAL*8 trialStress(6), ST, YT, SL, etaL, etaT, lambda, kappa
        ! local variables
        REAL*8 a(31), b(31), pi, tN, tT, tL
        REAL*8 trialFI_MT, trialFI_MC, aFail_MC, aFail_MT
        INTEGER i, p
        ! output variables
        REAL*8 FI_MT, FI_MC
       
        pi = 4*atan(1.0d0) ! determine value of pi
        a = (/ (i, i=0,30) /) ! create array of integers from 0 to 30
 
        DO p = 1,31
            b(p) = a(p)*(pi/60) ! create array of angles from 0 to pi/2
        END DO
       
        FI_MT = 0.0d0 ! initialize failure criteria
        FI_MC = 0.0d0

        DO p = 1, 31 ! iterate over angles
            ! Eq 3 CLN (expanded in 59-61)
            tN = trialStress(2)*cos(b(p))**2 + 2.0d0*trialStress(5)
     1          *sin(b(p))*cos(b(p)) + trialStress(3)*sin(b(p))**2
            tT = -1.0*cos(b(p))*sin(b(p))
     1           *(trialStress(2)-trialStress(3))
     2           +(trialStress(5)*(cos(b(p))**2.0 - sin(b(p))**2.0))
            tL = trialStress(4)*cos(b(p)) + trialStress(6)*sin(b(p))

            IF (tN.ge.0.0d0) THEN
                trialFI_MT = (tN/ST)**2 + (tL/SL)**2 + (tT/ST)**2 
     1                       + lambda*(tN/ST)*(tL/SL)**2 
     2                       + kappa*(tN/ST) ! Eq. 42 CLN
                trialFI_MC = 0.0
            ELSE
                trialFI_MC = (tL/(SL-etaL*tN))**2 
     1                      + (tT/(ST-etaT*tN))**2 ! Eq. 5 CLN
                trialFI_MT = 0.0
            ENDIF

            IF (trialFI_MT.gt.FI_MT) THEN
                FI_MT = trialFI_MT
                aFail_MT = b(p) ! record failure plane 
            END IF
            IF (trialFI_MC.gt.FI_MC) THEN
                FI_MC = trialFI_MC
                aFail_MC = b(p) ! record failure plane 
            END IF

        END DO

      RETURN
      END  

* <<<<<<<<<<<<<<<<<<<<<< SUBROUTINE ROTATE_PHI>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* ROTATION OF STRESSES TO THE MISALIGNMENT COORDINATE FRAME *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE ROTATE_PHI(trialStress, phiC, XC, trialStressP, theta,
     1                      phi)
        IMPLICIT NONE
        ! input variables
        REAL*8 trialStress(6), phiC, XC, NEWT, eps
        ! local variables
        REAL*8 theta, phi, X, m, n, u, v, gamma0
        REAL*8 trialStressT(6)
        PARAMETER (eps=1.d-8)
        LOGICAL tS4EQ0, tS6EQ0
        ! output variables
        REAL*8 trialStressP(6)

        ! first determine fracture plane angle theta (fiber kinking)
        tS4EQ0 = (abs(trialStress(4)-0.d0).lt.eps)
        tS6EQ0 = (abs(trialStress(6)-0.d0).lt.eps)
        IF (tS4EQ0.AND.tS6EQ0) THEN
            IF (abs(trialStress(2)-trialStress(3)).lt.eps) THEN
                theta = atan(1.0d0) ! pi/4
            ELSE
                theta = 0.5d0*atan((2.0d0*trialStress(5))
     1                  /(trialStress(2)-trialStress(3))) !Eq.55 CLN
            ENDIF 
        ELSE
            IF (abs(trialStress(4)-0.d0).lt.eps) THEN
                theta = 2.0d0*atan(1.0d0) ! pi/2
            ELSE 
                theta = atan(trialStress(6)/trialStress(4)) !Eq. 56 CLN
            ENDIF
        END IF
        
        ! determine the misalignment angle phi
        X = (sin(2.0d0*phiC)*XC)/(2.0d0*phiC) ! Eq. 86 CLN
        gamma0 = 0.1
        IF (trialStress(4).ge.0) THEN
            phi = NEWT(gamma0, trialStress, X) ! initial value of 0.1
        ELSE 
            phi = -1.0d0*NEWT(gamma0, trialStress, X)
        END IF

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

        RETURN
        END

* <<<<<<<<<<<<<<<<<<<<<<<<< FUNCTION NEWT >>>>>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Newton Raphson method to determine miaslignment angle phi
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      FUNCTION NEWT(gamma0, trialStress, X)
        IMPLICIT NONE
        REAL*8 NEWT, gamma, gamma0, d, trialStress(6), X
        REAL*8 FUNC, DFUNC, fx, fxprime, tol
        INTEGER k, maxiter
        PARAMETER (maxiter = 100, tol = 1.d-6)
        
        ! intial guess
        gamma = gamma0
!        WRITE(1,11) gamma
11      FORMAT('Initial guess: gamma = ', e22.15)
        ! Newton iteration to find a zero of f(x) 
        DO k = 1,maxiter
          ! evaluate function and its derivative:
          fxprime = DFUNC(gamma,trialStress,X)  
          fx = FUNC(gamma,trialStress,X)
!          WRITE(1,*) 'fxprime', fxprime
!          WRITE(1,*) 'fx', fx
          
          IF (abs(fx) < tol) THEN
              EXIT  ! jump out of do loop
          ENDIF
          ! compute Newton increment x:
          d = fx/fxprime
          ! update x:
          gamma = gamma - d
!          WRITE(1,12) k,gamma
12        FORMAT('After', i3, ' iterations, x = ', e22.15)
        ENDDO
        
        NEWT = gamma
        
        RETURN
        END

      ! Function used in the Newton Raphson iteration method
      FUNCTION FUNC(gamma1, trialStress, X)
        IMPLICIT NONE
        REAL*8 FUNC, X, gamma1, trialStress(6)
        ! Eq. 88 CLN
!        WRITE(1,*) 'gamma1', gamma1
!        WRITE(1,*) 'trialStress1', trialStress
        FUNC = X*gamma1 + 0.5d0*(trialStress(1) 
     1         - trialStress(2))*sin(2.0d0*gamma1)
     2         - abs(trialStress(4))*cos(2.0d0*gamma1) 
      RETURN
      END

      ! Derivative function of FUNC
      FUNCTION DFUNC(gamma2, trialStress, X)
        IMPLICIT NONE
        REAL*8 DFUNC, X, gamma2, trialStress(6)
        ! Eq. 89 CLN
!        WRITE(1,*) 'gamma2', gamma2
!        WRITE(1,*) 'trialStress2', trialStress
        DFUNC = X + (trialStress(1)-trialStress(2))*cos(2.0d0*gamma2)
     1          + 2.0d0*abs(trialStress(4))*sin(2.0d0*gamma2) 
      RETURN
      END

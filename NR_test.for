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
      REAL*8 FI_LT,FI_LC,FI_MT,FI_MC, gamma0
      REAL*8 initialcharLength,thickness,traceStrain,meanDamage,expo,A
      REAL*8 trialStressP(6), X, NEWT, FUNC, DFUNC, phi
      REAL*8 xOmega1,xOmega2,xOmega3,xOmega4,xOmega5,xOmega6,aminDamage
      CHARACTER*80 cmname

      OPEN(1, file = 
     1 'C:\Workspace\catalanotti_failure_criteria\NR_data.txt',
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
     
      trialStress(1) = -1.2001
      trialStress(2) = 0.0
      trialStress(3) = 0.0
      trialStress(4) = 0.0
      trialStress(5) = 0.0
      trialStress(6) = 0.0

      gamma0 = 0.1
      X = (sin(2.0d0*phiC)*XC)/(2.0*phiC)
      phi = NEWT(gamma0,trialStress,X)
      WRITE(1,*) 'phi = ', phi

      CLOSE(1)
      END

* <<<<<<<<<<<<<<<<<<<<<<<<< FUNCTION NEWT >>>>>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Newton Raphson method to determine miaslignment angle phi
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      FUNCTION NEWT(gamma0, trialStress, X)
        IMPLICIT NONE
        REAL*8 NEWT, gamma, gamma0, delta, trialStress(6), X
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
          delta = fx/fxprime
          ! update x:
          gamma = gamma - delta
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

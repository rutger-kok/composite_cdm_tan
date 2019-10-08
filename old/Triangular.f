C----------------------------------------------------------------------C
C     ABAQUS VUMAT USER SUBROUTINE: catalanotti_FC.for                 C
C     Author(s): Rutger Kok, Francisca Martinez-Hergueta               C
C     Date: 28/06/2019                                                 C
C     Version 1.0                                                      C
C----------------------------------------------------------------------C

C User subroutine VUMAT
      SUBROUTINE vumat (
C Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )

      INCLUDE 'vaba_param.inc'

      dimension coordMp(nblock,*), charLength(nblock), props(nprops),
     1     density(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr), 
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
     
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
      REAL*8 FI_LT,FI_LC,FI_MT,FI_MC
      REAL*8 initialcharLength,thickness,traceStrain,meanDamage,expo,A
      REAL*8 trialStressP(6), X, NEWT, FUNC, DFUNC
      REAL*8 xOmega1,xOmega2,xOmega3,xOmega4,xOmega5,xOmega6,aminDamage
      REAL*8 delta_0, delta_c
      CHARACTER*80 cmname
c    !   OPEN(1, file = 
c    !  1 'C:\Workspace\catalanotti_failure_criteria\FI_data.txt',
c    !  2 status = 'unknown')

C Elastic constants orthotropic ply

      e11 = props(1)
      e22 = props(2)
      e33 = props(3)
      nu12 = props(4)
      nu13 = props(5)
      nu23 = props(6)
      g12 = props(7)
      g13 = props(8)
      g23 = props(9)

C Ply strength

      XT = props(10)
      XC = props(11)
      YT = props(12)
      YC = props(13)
      SL = props(14)

C Fracture Angle

      alpha0 = props(15)*0.017453292519943295 !converts to radians

C Fracture toughness

      G1plus = props(16)
      G1minus = props(17)
      G2plus = props(18)
      G2minus = props(19)
      G6 = props(20)

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

C Loop through the gauss points

      IF (stepTime.eq.0) THEN     

      ! Initial elastic step, for Abaqus tests
      DO k = 1, nblock
      
       ! Initialisation of state variables
       DO k1 = 1,nstatev
            stateNew(k,k1) = 0.d0
       ENDDO

       DO i = 1,6
           trialStrain(i)=strainInc(k,i)
       ENDDO

       trialStress(1) = d11*trialStrain(1)+d12*trialStrain(2)
     1                  +d13*trialStrain(3)
       trialStress(2) = d12*trialStrain(1)+d22*trialStrain(2)
     1                  +d23*trialStrain(3)
       trialStress(3) = d13*trialStrain(1)+d23*trialStrain(2)
     1                  +d33*trialStrain(3)
       trialStress(4) = d44*trialStrain(4)
       trialStress(5) = d55*trialStrain(5)
       trialStress(6) = d66*trialStrain(6)

       DO i = 1,6
           stressNew(k,i)=trialStress(i)
       ENDDO

      ENDDO

      ELSE

       ! Constitutive model definition

       ! Update of the failure thresholds (r values)
       DO k = 1,nblock
        rOld1 = stateOld(k,7)
        rOld2 = stateOld(k,8)
        rOld3 = stateOld(k,9)
        rOld4 = stateOld(k,10)
        rOld5 = stateOld(k,11)
        rOld6 = stateOld(k,12)

        ! Computation of the total strain
        DO i = 1,6
            stateNew(k,i)=stateOld(k,i)+strainInc(k,i)
        ENDDO

        DO i = 1,6
            trialStrain(i)=stateNew(k,i)
        ENDDO

        !Trial stress
        trialStress(1) = d11*trialStrain(1)+d12*trialStrain(2)
     1                   +d13*trialStrain(3)
        trialStress(2) = d12*trialStrain(1)+d22*trialStrain(2)
     1                   +d23*trialStrain(3)
        trialStress(3) = d13*trialStrain(1)+d23*trialStrain(2)
     1                   +d33*trialStrain(3)
        trialStress(4) = d44*trialStrain(4)
        trialStress(5) = d55*trialStrain(5)
        trialStress(6) = d66*trialStrain(6)

        DO i = 1,6
            stressNew(k,i)=trialStress(i)
        ENDDO

        ! Evaluation of the damage activation functions

        ! longitudinal failure criteria
        IF (trialStress(1).gt.0.d0) THEN
            FI_LT = trialStrain(1)/(XT/e11) ! Eq. 54 CLN
        ELSEIF (trialStress(1).lt.0.d0) THEN
            call ROTATE_PHI(trialStress,phiC,XC,trialStressP)
            call FAIL_CLN(trialStressP,ST,YT,SL,etaL,etaT,lambda,kappa,
     1                    FI_LT,FI_LC)
        ENDIF
    
        ! transverse failure criteria
        call FAIL_CLN(trialStress,ST,YT,SL,etaL,etaT,lambda,kappa,
     1                    FI_MT,FI_MC)

! 100     FORMAT (E,A,E,A,E,A,E,A,E,A,E,A,E,A,E,A,E,A,E)
!         WRITE(1,100) trialStress(1), ',',
!      1  trialStress(2), ',', trialStress(3), ',',
!      2  trialStress(4), ',', trialStress(5), ',',
!      3  trialStress(6), ',', FI_LT, ',', FI_LC, ',',
!      4  FI_MT, ',', FI_MC

        ! Update of the damage thresholds
        rNew1 = max(1.0d0,max(FI_LT,rOld1),max(FI_LC,rOld2))    !Eq.26
        rNew2 = max(1.0d0,max(FI_LC,rOld2))
        rNew3 = max(1.0d0,max(FI_MT,rOld3),max(FI_MC,rOld4))
        rNew4 = max(1.0d0,max(FI_MC,rOld4))
        rNew5 = 1.0d0
        rNew6 = 1.0d0

        ! Softening parameters
c        initialcharLength = 1.0
        initialcharLength = charLength(k)
        
        A1plus = 2.0d0*initialcharLength*XT*XT  !Second part. eq. B.2       
     1           /(2.0d0*e11*G1plus-initialcharLength*XT*XT)
        
        A1minus = 2.0d0*initialcharLength  ! Second part. Appendix B
     1            *(1.2671444971783234*XC)**2.0d0 
     2            /(2.0d0*(1.173191594009027*e11)
     3            *G1minus-initialcharLength
     4            *(1.145996547107106*XC)**2.0d0)
      
        A2plus = 2.0d0*initialcharLength*(0.8308920721648098*YT)**2.0d0
     1           /(2.0d0*0.9159840495203727*e22*G2plus
     2           -initialcharLength*(0.9266465446084761*YT)**2.0d0)

        A2minus = 2.0d0*initialcharLength*(0.709543*YC)**2.0d0
     1            /(2.0d0*0.7881*e22*G2minus-initialcharLength
     2            *(0.802572*YC)**2.0d0)

        A6 = 2.0d0*initialcharLength*SL*SL
     1       /(2.0d0*g12*G6-initialcharLength*SL*SL)
     
        A6 = A2plus
 
        !Eq.5 (from part II of the Maimi paper)
        d1Plus = 1.0d0-1.0d0/rNew1*exp(A1plus*(1.0d0-rNew1)) 
c        d1Minus = 1.0d0-1.0d0/rNew2*exp(A1minus*(1.0d0-rNew2))
        d2Plus = 1.0d0-1.0d0/rNew3*exp(A2plus*(1.0d0-rNew3))
        d2Minus = 1.0d0-1.0d0/rNew4*exp(A2minus*(1.0d0-rNew4))
        d6 = 1.0d0-1.0d0/rNew3*exp(A6*(1.0d0-rNew3))*(1.0d0-d1Plus)

        if (rNew2.gt.1.0d0) then
                delta_0 = XC/e11
                delta_c = (2.0d0*G1minus/XC*1.0d0)
                if (delta_c.lt.delta_0)then               !!Needed for d.nq.0-0
                   delta_c = 1.5*delta_0             !! 1.1 adjustable for stability
                endif
                d1Minus = (delta_c*(dabs(trialStrain(1))-delta_0))
     *            /(dabs(trialStrain(1))*(delta_c-delta_0))
                d1Minus = min(0.999,d1Minus)
        endif

        ! Damage variables

        IF (trialStress(1).gt.0) THEN      !Eq.6
            xOmega1 = d1Plus
        ELSE
            xOmega1 = d1Minus
        ENDIF

        IF (trialStress(2).gt.0) THEN
            xOmega2 = d2Plus
        ELSE
            xOmega2 = d2Minus
        ENDIF

c        xOmega4 = d6
        xOmega4 = 0.0d0
        xOmega3 = 0.0d0
        xOmega5 = 0.0d0
        xOmega6 = 0.0d0 

        ! Stiffness tensor + damage
        nu21 = nu12*(e22/e11)
        nu31 = nu13*(e33/e11)
        nu32 = nu23*(e33/e22)
 
        delta = 1.0d0/(1.0d0-nu12*nu21*(1.0d0-xOmega1)*(1.0d0-xOmega2)
     1        -nu32*nu23*(1.0d0-xOmega2)*(1.0d0-xOmega3)
     2        -nu13*nu31*(1.0d0-xOmega1)*(1.0d0-xOmega3)
     3        -2.0d0*nu12*nu23*nu31*(1.0d0-xOmega1)
     4        *(1.0d0-xOmega2)*(1.0d0-xOmega3))

        d11 = e11*(1.0d0-xOmega1)*(1.0d0-nu23*nu32*(1.0d0-xOmega2)*
     1             (1.0d0-xOmega3))*delta
        d22 = e22*(1.0d0-xOmega2)*(1.0d0-nu13*nu31*(1.0d0-xOmega1)*
     1             (1.0d0-xOmega3))*delta 
        d33 = e33*(1.0d0-xOmega3)*(1.0d0-nu12*nu21*(1.0d0-xOmega1)*
     1             (1.0d0-xOmega2))*delta
        d12 = e11*(1.0d0-xOmega1)*(1.0d0-xOmega2)*(nu21+nu31*nu23*
     1             (1.0d0-xOmega3))*delta
        d13 = e11*(1.0d0-xOmega1)*(1.0d0-xOmega3)*(nu31+nu21*nu32*
     1             (1.0d0-xOmega2))*delta
        d23 = e22*(1.0d0-xOmega2)*(1.0d0-xOmega3)*(nu32+nu12*nu31*
     1             (1.0d0-xOmega1))*delta
        d44 = 2.0d0*g12*(1.0d0-xOmega4)
        d55 = 2.0d0*g23*(1.0d0-xOmega5)
        d66 = 2.0d0*g13*(1.0d0-xOmega6)

        trialStress(1) = d11*trialStrain(1)+d12*trialStrain(2)
     1                   +d13*trialStrain(3)
        trialStress(2) = d12*trialStrain(1)+d22*trialStrain(2)
     1                   +d23*trialStrain(3)
        trialStress(3) = d13*trialStrain(1)+d23*trialStrain(2)
     1                   +d33*trialStrain(3)
        trialStress(4) = d44*trialStrain(4)
        trialStress(5) = d55*trialStrain(5)
        trialStress(6) = d66*trialStrain(6)

        DO i = 1,6
            stressNew(k,i) = trialStress(i)
        ENDDO

        ! Energy
        stressPower = 0.5d0*(
     1              (stressNew(k,1)+stressOld(k,1))*strainInc(k,1)+
     2              (stressNew(k,2)+stressOld(k,2))*strainInc(k,2)+
     3              (stressNew(k,3)+stressOld(k,3))*strainInc(k,3)+
     4              2.0*(stressNew(k,4)+stressOld(k,4))*strainInc(k,4)+
     5              2.0*(stressNew(k,5)+stressOld(k,5))*strainInc(k,5)+
     6              2.0*(stressNew(k,6)+stressOld(k,6))*strainInc(k,6))

        enerInternNew(k) = enerInternOld(k)+stressPower/density(k)

        stateNew(k,7) = rNew1
        stateNew(k,8) = rNew2
        stateNew(k,9) = rNew3
        stateNew(k,10) = rNew4
        stateNew(k,11) = rNew5
        stateNew(k,12) = rNew6
        stateNew(k,13) = xOmega1
        stateNew(k,14) = xOmega2
        stateNew(k,15) = xOmega3
        stateNew(k,16) = xOmega4
        stateNew(k,17) = delta_c
        stateNew(k,18) = delta_0

        aminDamage = min(xOmega1,xOmega2,xOmega4)
        stateNew(k,19) = aminDamage

        ! single line if (so no THEN & no ENDIF)
        IF (aminDamage.gt.0.999) stateNew(k,20) = 0.d0

       ENDDO ! this ends the do loop statement at line ~157

      ENDIF ! this ends the if statement at line ~122
      
      RETURN

      CLOSE(1)

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

c    !   CLOSE(1)

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

  
* <<<<<<<<<<<<<<<<<<<<<< SUBROUTINE ROTATE_PHI>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* ROTATION OF STRESSES TO THE MISALIGNMENT COORDINATE FRAME *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE ROTATE_PHI(trialStress, phiC, XC, trialStressP)
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
        PRINT*,'Theta = ', theta
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
  

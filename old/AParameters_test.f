C----------------------------------------------------------------------C
C     ABAQUS VUMAT USER SUBROUTINE: catalanotti_FC.for                 C
C     Author(s): Rutger Kok, Francisca Martinez-Hergueta               C
C     Date: 28/06/2019                                                 C
C     Version 1.0                                                      C
C----------------------------------------------------------------------C

C User subroutine VUMAT
      PROGRAM AParameters_test

      REAL*8 trialStress(6),trialStrain(6),trialStressP(6)
      REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
      REAL*8 nu21,nu31,nu32,delta
      REAL*8 C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9)
      REAL*8 XT,XC,YT,YC,SL,ST
      REAL*8 G1plus,G1minus,G2plus,G2minus,G6
      REAL*8 A1plus,A1minus,A2plus,A2minus,A6
      REAL*8 alpha0,etaT,etaL,phiC,kappa,lambda,initialcharLength
      REAL*8 failStrain,tStrainArray(50),tStressArray(50),energy
      INTEGER n,i
      OPEN(1, file = 
     1 'C:\Workspace\catalanotti_failure_criteria\A_parameter.txt',
     2 status = 'unknown')
      WRITE(1,*) 'trialStress(1)', ',',
     1  'trialStress(2)', ',', 'trialStress(3)', ',',
     2  'trialStress(4)', ',', 'trialStress(5)', ',',
     3  'trialStress(6)', ',', 'trialStrain(1)', ',',
     4  'trialStrain(2)', ',', 'trialStrain(3)', ',',
     5  'trialStrain(4)', ',', 'trialStrain(5)', ',',
     6  'trialStrain(6)'

C Elastic constants orthotropic ply

      e11 = 154.0
      e22 = 8.5
      e33 = 8.5
      nu12 = 0.309053
      nu13 = 0.309053
      nu23 = 0.33
      g12 = 4.2
      g13 = 4.2
      g23 = 3.2

C Ply strength

      XT = 2.61
      XC = 1.759
      YT = 0.055
      YC = 0.285
      SL = 0.105

C Fracture Angle

      alpha0 = 53.0*0.017453292519943295 !converts to radians

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

C Stiffness matrix orthotropic material

      nu21 = nu12*(e22/e11)
      nu31 = nu13*(e33/e11)
      nu32 = nu23*(e33/e22)

      delta = (1.0-nu12*nu21-nu23*nu32-nu13*nu31-2.0*nu21*nu32*nu13)
     1        /(e11*e22*e33)

      C(1) = (1-nu23*nu32)/(e22*e33*delta)
      C(2) = (1-nu13*nu31)/(e11*e33*delta)
      C(3) = (1-nu12*nu21)/(e11*e22*delta)
      C(4) = (nu21+nu23*nu31)/(e22*e33*delta)
      C(5) = (nu31+nu21*nu32)/(e22*e33*delta)
      C(6) = (nu32+nu12*nu31)/(e11*e33*delta)
      C(7) = 2.0*g12
      C(8) = 2.0*g23
      C(9) = 2.0*g13 

C Calculate A parameters
  
      initialcharLength = 1.0
        
      ! initial approximation
      A1plus = 2.0d0*initialcharLength*XT*XT  !Second part. eq. B.2       
     1         /(2.0d0*e11*G1plus-initialcharLength*XT*XT)

      A1minus = 2.0d0*initialcharLength*XC*XC        
     1          /(2.0d0*e11*G1minus-initialcharLength*XC*XC)

      A2plus = 2.0d0*initialcharLength*YT*YT         
     1         /(2.0d0*e22*G2plus-initialcharLength*YT*YT)

      A2minus = 2.0d0*initialcharLength*YC*YC      
     1         /(2.0d0*e22*G2minus-initialcharLength*YC*YC)

      A6 = 2.0d0*initialcharLength*SL*SL
     1     /(2.0d0*g12*G6-initialcharLength*SL*SL)

      n = 50
      failStrain = 0.16
      DO i = 0,n-1
        trialStrain(1) = (failStrain/n)*i
        trialStrain(2) = 0.d0
        trialStrain(3) = 0.d0
        trialStrain(4) = 0.d0
        trialStrain(5) = 0.d0
        trialStrain(6) = 0.d0
        tStrainArray(i) = trialStrain(1)

        call CDM(trialStrain,C,ST,YT,SL,etaL,etaT,lambda, kappa,
     1               XC,XT,YC,A1Plus,A1Minus,A2Plus,A2Minus,A6,
     2               phiC,e11,nu12,nu13,nu23,g12,g13,g23,trialStress)
100     FORMAT (E,A,E,A,E,A,E,A,E,A,E,A,E,A,E,A,E,A,E,A,E,A,E)
        WRITE(1,100) trialStress(1), ',',
     1  trialStress(2), ',', trialStress(3), ',',
     2  trialStress(4), ',', trialStress(5), ',',
     3  trialStress(6), ',', trialStrain(1), ',',
     4  trialStrain(2), ',', trialStrain(3), ',',
     5  trialStrain(4), ',', trialStrain(5), ',',
     6  trialStrain(6)

        tStressArray(i) = trialStress(1)
      ENDDO 

      call TRAPZD(tStrainArray,tStressArray,energy)
      PRINT*, 'Energy:' energy

      CLOSE(1)
      PAUSE "Press RETURN to end program"
      END PROGRAM AParameters_test

* <<<<<<<<<<<<<<<<<<<<<< SUBROUTINE CDM >>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* CDM from Maimi et al. used to determine A parameters
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE CDM(trialStrain,C,ST,YT,SL,etaL,etaT,lambda, kappa,
     1               XC,XT,YC,A1Plus,A1Minus,A2Plus,A2Minus,A6,
     2               phiC,e11,nu12,nu13,nu23,g12,g13,g23,trialStress)
        IMPLICIT NONE
        ! input variables
        REAL*8 trialStrain(6), ST, YT, SL, etaL, etaT, lambda, kappa
        REAL*8 C(9), XT, XC, YC, phiC, e11, e22, e33, nu12, nu13, nu23
        REAL*8 g12, g13, g23
        ! local variables
        REAL*8 FI_MT,FI_MC,FI_LT,FI_LC, delta, trialStressP(6)
        REAL*8 d1Plus,d1Minus,d2Plus,d2Minus,d6
        REAL*8 rNew1,rNew2,rNew3,rNew4,rNew5,rNew6
        REAL*8 xOmega1,xOmega2,xOmega3,xOmega4,xOmega5,xOmega6
        REAL*8 nu21, nu32, nu31
        ! output variables
        REAL*8 trialStress(6)

        !Trial stress
        trialStress(1) = C(1)*trialStrain(1)+C(4)*trialStrain(2)
     1                   +C(5)*trialStrain(3)
        trialStress(2) = C(4)*trialStrain(1)+C(2)*trialStrain(2)
     1                   +C(6)*trialStrain(3)
        trialStress(3) = C(5)*trialStrain(1)+C(6)*trialStrain(2)
     1                   +C(3)*trialStrain(3)
        trialStress(4) = C(7)*trialStrain(4)
        trialStress(5) = C(8)*trialStrain(5)
        trialStress(6) = C(9)*trialStrain(6)

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

        ! Update of the damage thresholds
        rNew1 = max(1.0,max(FI_LT,FI_LC))    !Eq.26
        rNew2 = max(1.0,FI_LC)
        rNew3 = max(1.0,max(FI_MT,FI_MC))
        rNew4 = max(1.0,FI_MC)
        rNew5 = 1.0d0
        rNew6 = 1.0d0
 
        !Eq.5 (from part II of the Maimi paper)
        d1Plus = 1.0d0-1.0d0/rNew1*exp(A1plus*(1.0d0-rNew1)) 
        d1Minus = 1.0d0-1.0d0/rNew2*exp(A1minus*(1.0d0-rNew2))
        d2Plus = 1.0d0-1.0d0/rNew3*exp(A2plus*(1.0d0-rNew3))
        d2Minus = 1.0d0-1.0d0/rNew4*exp(A2minus*(1.0d0-rNew4))
        d6 = 1.0d0-1.0d0/rNew3*exp(A6*(1.0d0-rNew3))*(1.0d0-d1Plus)

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

        xOmega4 = d6
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

        C(1) = e11*(1.0d0-xOmega1)*(1.0d0-nu23*nu32*(1.0d0-xOmega2)*
     1             (1.0d0-xOmega3))*delta
        C(2) = e22*(1.0d0-xOmega2)*(1.0d0-nu13*nu31*(1.0d0-xOmega1)*
     1             (1.0d0-xOmega3))*delta 
        C(3) = e33*(1.0d0-xOmega3)*(1.0d0-nu12*nu21*(1.0d0-xOmega1)*
     1             (1.0d0-xOmega2))*delta
        C(4) = e11*(1.0d0-xOmega1)*(1.0d0-xOmega2)*(nu21+nu31*nu23*
     1             (1.0d0-xOmega3))*delta
        C(5) = e11*(1.0d0-xOmega1)*(1.0d0-xOmega3)*(nu31+nu21*nu32*
     1             (1.0d0-xOmega2))*delta
        C(6) = e22*(1.0d0-xOmega2)*(1.0d0-xOmega3)*(nu32+nu12*nu31*
     1             (1.0d0-xOmega1))*delta
        C(7) = 2.0d0*g12*(1.0d0-xOmega4)
        C(8) = 2.0d0*g23*(1.0d0-xOmega5)
        C(9) = 2.0d0*g13*(1.0d0-xOmega6)

        trialStress(1) = C(1)*trialStrain(1)+C(4)*trialStrain(2)
     1                   +C(5)*trialStrain(3)
        trialStress(2) = C(4)*trialStrain(1)+C(2)*trialStrain(2)
     1                   +C(6)*trialStrain(3)
        trialStress(3) = C(5)*trialStrain(1)+C(6)*trialStrain(2)
     1                   +C(3)*trialStrain(3)
        trialStress(4) = C(7)*trialStrain(4)
        trialStress(5) = C(8)*trialStrain(5)
        trialStress(6) = C(9)*trialStrain(6)
      
      RETURN
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

* <<<<<<<<<<<<<<<<<<<<<<< SUBROUTINE TRAPZD >>>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Trapezoidal integration subroutine
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE TRAPZD(xarray,yarray,result)
        IMPLICIT NONE        
C       Input variables
        REAL*8 xarray(50),yarray(50)
C       Local variables
        INTEGER i
        REAL*8 s,x1,x2,y1,y2,delta
C       Output variables
        REAL*8 result

        result = 0.d0
        DO i = 1,49
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
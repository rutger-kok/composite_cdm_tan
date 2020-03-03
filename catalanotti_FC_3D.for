!----------------------------------------------------------------------!
!     ABAQUS VUMAT USER SUBROUTINE: catalanotti_FC_v3.for              !
!     Author(s): Rutger Kok, Francisca Martinez-Hergueta               !
!     Date: 31/10/2019                                                 !
!     Version 1.0                                                      !
!----------------------------------------------------------------------!

! User subroutine VUMAT

      SUBROUTINE vumat (
      ! Read only -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
      ! Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew )

      INCLUDE 'vaba_param.inc'
      
      DIMENSION props(nprops), density(nblock), coordMp(nblock),
     1 charLength(nblock), strainInc(nblock, ndir+nshr),
     2 relSpinInc(nblock, nshr), tempOld(nblock),
     3 stretchOld(nblock, ndir+nshr),defgradOld(nblock,ndir+nshr+nshr),
     4 fieldOld(nblock, nfieldv), stressOld(nblock, ndir+nshr),
     5 stateOld(nblock, nstatev), enerInternOld(nblock),
     6 enerInelasOld(nblock), tempNew(nblock),
     7 stretchNew(nblock, ndir+nshr),defgradNew(nblock,ndir+nshr+nshr),
     8 fieldNew(nblock, nfieldv), stressNew(nblock,ndir+nshr),
     9 stateNew(nblock, nstatev), enerInternNew(nblock),
     1 enerInelasNew(nblock)
      CHARACTER*80 cmname
      
      REAL*8, DIMENSION(6) :: trialStress,trialStrain
      REAL*8, DIMENSION(6,6) :: C
      REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23 ! elastic constants
      REAL*8 nu21,nu31,nu32,delta
      REAL*8 XT,XC,YT,YC,SL,alpha0,stressPower
      REAL*8 G1plus,G1minus,G2plus,G2minus,G6 ! fracture energies
      REAL*8, PARAMETER :: tol=1.0d-5 ! for Brent method convergence
      INTEGER k,k1,i,j,f,myTmpVar
    !   OPEN(105,file='C:\Workspace\catalanotti_final\test.txt',status='unknown',
    !  1       action='readwrite',position='append')

      ! Elastic constants orthotropic ply

      e11 = props(1)
      e22 = props(2)
      e33 = props(3)
      nu12 = props(4)
      nu13 = props(5)
      nu23 = props(6)
      g12 = props(7)
      g13 = props(8)
      g23 = props(9)

      ! Ply strengths

      XT = props(10)
      XC = props(11)
      YT = props(12)
      YC = props(13)
      SL = props(14)

      ! Fracture Angle

      alpha0 = props(15)*0.017453292519943295 !converts to radians

      ! Fracture toughness

      G1plus = props(16)
      G1minus = props(17)
      G2plus = props(18)
      G2minus = props(19)
      G6 = props(20)

      ! Stiffness matrix orthotropic material

      nu21 = nu12*(e22/e11)
      nu31 = nu13*(e33/e11)
      nu32 = nu23*(e33/e22)
      
      !DO WHILE(myTmpVar.ne.999)
        !myTmpVar = 1
      !END DO 

      delta = 1.0d0/(e22**2.0d0*nu12**2.0d0 +
     1        2.0d0*e33*e22*nu12*nu13*nu23 +
     1        e33*e22*nu13**2.0d0 - e11*e22 + e11*e33*nu23**2.0d0)

      C = 0.0d0 ! set all elements in array equal to zero

      C(1,1) = -(e11**2.0d0*(-e33*nu23**2.0d0 + e22))*delta
      C(2,1) = -(e11*e22*(e22*nu12 + e33*nu13*nu23))*delta
      C(3,1) = -(e11*e22*e33*(nu13 + nu12*nu23))*delta
      C(1,2) = -(e11*e22*(e22*nu12 + e33*nu13*nu23))*delta
      C(2,2) = -(e22**2.0d0*(-e33*nu13**2.0d0 + e11))*delta
      C(3,2) = -(e22*e33*(e11*nu23 + e22*nu12*nu13))*delta
      C(1,3) = -(e11*e22*e33*(nu13 + nu12*nu23))*delta
      C(2,3) = -(e22*e33*(e11*nu23 + e22*nu12*nu13))*delta
      C(3,3) = -(e22*e33*(-e22*nu12**2.0d0 + e11))*delta
      C(4,4) = g12
      C(5,5) = g23
      C(6,6) = g13
 
      ! Loop through the gauss points
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

          trialStress = matmul(C(1:6,1:6), trialStrain(1:6))

          DO i = 1,6
              stressNew(k,i)=trialStress(i)
          ENDDO

        ENDDO

      ELSE

        DO k = 1,nblock
          
            DO i = 1,6
              stateNew(k,i)=stateOld(k,i)+strainInc(k,i)
          ENDDO

          DO i = 1,6
              trialStrain(i)=stateNew(k,i)
          ENDDO

          call cdm(k,nblock,nstatev,trialStrain,stateOld,charLength(k),
     1             C,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2             XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,G2minus,G6,
     3             alpha0,trialStress,stateNew)

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
    
        ENDDO
      ENDIF
    !   CLOSE(105)

      END SUBROUTINE vumat

!----------------------------------------------------------------------!
! Subroutine cdm: Failure criteria and damage evolution.               !
! Returns a trialStress matrix                                         !
!----------------------------------------------------------------------!

      SUBROUTINE cdm(k,nblock,nstatev,trialStrain,stateOld,l_ch,C,
     1               e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2               XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,G2minus,G6,
     3               alpha0,trialStress,stateNew)
        IMPLICIT NONE
        ! input variables
        INTEGER, INTENT(IN) :: nblock,nstatev,k
        REAL*8, DIMENSION(nblock,nstatev), INTENT(IN) :: stateOld
        REAL*8, DIMENSION(6,6), INTENT(INOUT) :: C
        REAL*8, DIMENSION(6), INTENT(IN) :: trialStrain
        REAL*8, INTENT(IN) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
        REAL*8, INTENT(IN) :: XT,XC,YT,YC,SL,l_ch
        REAL*8, INTENT(IN) :: G1plus,G1minus,G2plus,G2minus,G6
        ! local variables
        REAL*8, DIMENSION(6) :: trialStressP
        REAL*8 d1Plus,d1Minus,d2Plus,d2Minus,d3Plus,d3Minus
        REAL*8 d1,d2,d3,d6
        REAL*8 nu21,nu31,nu32,delta,aminDamage
        REAL*8 alpha0,etaT,etaL,phiC,kappa,lambda,ST
        REAL*8 FI_LT,FI_LC,FI_MT,FI_MC
        REAL*8 xOmega1,xOmega2,xOmega3,xOmega4
        REAL*8 xOmega1Old,xOmega2Old,xOmega3Old,xOmega4Old
        ! output variables
        REAL*8, DIMENSION(6), INTENT(OUT) :: trialStress
        REAL*8, DIMENSION(nblock,nstatev), INTENT(INOUT) :: stateNew
        ! external functions
        REAL*8, EXTERNAL :: triangular, mccauley

        ! Initial values
        etaL = 0.280622004043348
        phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)*((SL/XC)+etaL))) 
     1         /(2.0d0*((SL/XC)+etaL)))  !Eq.12
        ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*SL)/  
     1       (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL)) !Eq.12 CLN
        etaT = (etaL*ST)/SL  !Eq.10 CLN
        kappa = (ST**2.0d0-YT**2.0)/(ST*YT)  !Eq.43 CLN
        lambda = ((2.0*etaL*ST)/SL)-kappa  !Eq.45 CLN 

        !Old xOmega values
        xOmega1Old = stateOld(k,13)
        xOmega2Old = stateOld(k,14)
        xOmega3Old = stateOld(k,15)
        xOmega4Old = stateOld(k,16)

        !Trial stress
        trialStress = matmul(C(1:6,1:6), trialStrain(1:6))

        ! Evaluation of the damage activation functions
        FI_LT = 0.0d0
        FI_LC = 0.0d0
        FI_MT = 0.0d0
        FI_MC = 0.0d0

        ! longitudinal failure criteria
        IF (trialStress(1).gt.0.d0) THEN
            FI_LT = trialStrain(1)/(XT/e11) ! Eq. 54 CLN
        ELSEIF (trialStress(1).lt.0.d0) THEN
            call rotate_stress(trialStress,phiC,XC,g12,trialStressP)
            call fail_cln(trialStressP,ST,SL,etaL,etaT,lambda,kappa,
     1                    FI_LT,FI_LC)
        ENDIF

        ! transverse failure criteria
        call fail_cln(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                FI_MT,FI_MC)

        IF (FI_LT.gt.1.0d0) THEN !Triangular damage dissipation
            d1Plus= triangular(XT,e11,G1Plus,trialStrain(1),l_ch)
        ELSE
            d1Plus = 0.0d0
        ENDIF

        IF (FI_LC.gt.1.0d0) THEN !Triangular damage dissipation
            d1Minus = triangular(XC,e11,G1minus,trialStrain(1),l_ch)
        ELSE
            d1Minus = 0.0d0
        ENDIF

        IF (FI_MT.gt.1.0d0) THEN !Triangular damage dissipation
            d2Plus = triangular(YT,e22,G2Plus,trialStrain(2),l_ch)
        ELSE
            d2Plus = 0.0d0
        ENDIF

        IF (FI_MC.gt.1.0d0) THEN
            d2Minus = triangular(YC,e22,G2minus,trialStrain(2),l_ch)
        ELSE
            d2Minus = 0.0d0
        ENDIF

        IF (FI_MT.gt.1.0d0) THEN !Triangular damage dissipation
            d3Plus = triangular(YT,e33,G2Plus,trialStrain(3),l_ch)
        ELSE
            d3Plus = 0.0d0
        ENDIF

        IF (FI_MC.gt.1.0d0) THEN
            d3Minus = triangular(YC,e33,G2minus,trialStrain(3),l_ch)
        ELSE
            d3Minus = 0.0d0
        ENDIF

        IF ((FI_MT.gt.1.0d0).or.(FI_MC.gt.1.0d0)) THEN !Triangular damage dissipation
            d6 = triangular(SL,e22,G6,trialStrain(4),l_ch)
        ELSE
            d6 = 0.0d0
        ENDIF

        ! Damage variables
        IF (trialStress(1).ne.0.0d0) THEN
          d1 = d1Plus*(mccauley(trialStress(1))/abs(trialStress(1))) + 
     1         d1Minus*(mccauley(-trialStress(1))/abs(trialStress(1)))
        ELSE
          d1 = 0.0d0
        END IF
        
        IF (trialStress(2).ne.0.0d0) THEN
          d2 = d2Plus*(mccauley(trialStress(2))/abs(trialStress(2))) + 
     1         d2Minus*(mccauley(-trialStress(2))/abs(trialStress(2)))
        ELSE
          d2 = 0.0d0
        END IF 
        
        IF (trialStress(3).ne.0.0d0) THEN
          d3 = d3Plus*(mccauley(trialStress(3))/abs(trialStress(3))) + 
     1         d3Minus*(mccauley(-trialStress(3))/abs(trialStress(3)))
        ELSE 
          d3 = 0.0d0 
        END IF

        xOmega1 = min(0.999,max(d1,xOmega1Old))
        xOmega2 = min(0.999,max(d2,xOmega2Old))
        xOmega3 = min(0.999,max(d3,xOmega3Old))
        xOmega4 = min(0.999,max(d6,xOmega4Old))

        ! Stiffness tensor + damage
        nu21 = nu12*(e22/e11)
        nu31 = nu13*(e33/e11)
        nu32 = nu23*(e33/e22)
 
        delta = 1.0d0/(e11*e22 - e22**2.0d0*nu12**2.0d0 + 
     1          e22**2.0d0*xOmega1*nu12**2.0d0 +
     2          e22**2.0d0*xOmega2*nu12**2.0d0 - 
     3          e22*e33*nu13**2.0d0 - e11*e33*nu23**2.0d0 +
     4          e22*e33*xOmega1*nu13**2.0d0 + 
     5          e22*e33*xOmega3*nu13**2.0d0 + 
     6          e11*e33*xOmega2*nu23**2.0d0 +
     7          e11*e33*xOmega3*nu23**2.0d0 -
     8          e22**2.0d0*xOmega1*xOmega2*nu12**2.0d0 -
     9          2.0d0*e22*e33*nu12*nu13*nu23 -
     1          e22*e33*xOmega1*xOmega3*nu13**2.0d0 -
     2          e11*e33*xOmega2*xOmega3*nu23**2.0d0 +
     3          2.0d0*e22*e33*xOmega1*nu12*nu13*nu23 +
     4          2.0d0*e22*e33*xOmega2*nu12*nu13*nu23 +
     5          2.0d0*e22*e33*xOmega3*nu12*nu13*nu23 -
     6          2.0d0*e22*e33*xOmega1*xOmega2*nu12*nu13*nu23 -
     7          2.0d0*e22*e33*xOmega1*xOmega3*nu12*nu13*nu23 -
     8          2.0d0*e22*e33*xOmega2*xOmega3*nu12*nu13*nu23 +
     9          2.0d0*e22*e33*xOmega1*xOmega2*xOmega3*nu12*nu13*nu23)

        C = 0.0d0 ! set all elements in array equal to zero

        C(1,1) = -(e11**2.0d0*(xOmega1 - 1.0d0)*(e22 - e33*nu23**2.0d0 +
     1           e33*xOmega2*nu23**2.0d0 + e33*xOmega3*nu23**2.0d0 -
     2           e33*xOmega2*xOmega3*nu23**2.0d0))*delta
        C(2,1) = (e11*e22*(xOmega1 - 1.0d0)*(xOmega2 - 1.0d0)*
     1           (e22*nu12 + e33*nu13*nu23 - e33*xOmega3*nu13*nu23))
     2           *delta
        C(3,1) = (e11*e22*e33*(xOmega1 - 1.0d0)*(xOmega3 - 1.0d0)*
     1           (nu13 + nu12*nu23 - xOmega2*nu12*nu23))*delta
        C(1,2) = (e11*e22*(xOmega1 - 1.0d0)*(xOmega2 - 1.0d0)*
     1           (e22*nu12 + e33*nu13*nu23 - e33*xOmega3*nu13*nu23))
     2           *delta
        C(2,2) = -(e22**2.0d0*(xOmega2 - 1.0d0)*(e11 - e33*nu13**2.0d0 +
     1           e33*xOmega1*nu13**2.0d0 + e33*xOmega3*nu13**2.0d0 -
     2           e33*xOmega1*xOmega3*nu13**2.0d0))*delta
        C(3,2) = (e22*e33*(xOmega2 - 1.0d0)*(xOmega3 - 1.0d0)*
     1           (e11*nu23 + e22*nu12*nu13 - e22*xOmega1*nu12*nu13))
     2           *delta
        C(1,3) = (e11*e22*e33*(xOmega1 - 1.0d0)*(xOmega3 - 1.0d0)*
     1           (nu13 + nu12*nu23 - xOmega2*nu12*nu23))*delta
        C(2,3) = (e22*e33*(xOmega2 - 1.0d0)*(xOmega3 - 1.0d0)*
     1           (e11*nu23 + e22*nu12*nu13 - e22*xOmega1*nu12*nu13))
     2           *delta
        C(3,3) = -(e22*e33*(xOmega3 - 1.0d0)*(e11 - e22*nu12**2.0d0 +
     1           e22*xOmega1*nu12**2.0d0 + e22*xOmega2*nu12**2.0d0 - 
     2           e22*xOmega1*xOmega2*nu12**2.0d0))*delta
        C(4,4) = -(g12*(xOmega4 - 1.0d0))
        C(5,5) = g23
        C(6,6) = -(g13*(xOmega4 - 1.0d0))
        
        trialStress = matmul(C(1:6,1:6), trialStrain(1:6))

        stateNew(k,7) = C(1,1)
        stateNew(k,8) = C(2,1)
        stateNew(k,9) = C(3,1)
        stateNew(k,10) = C(1,2)
        stateNew(k,11) = C(2,2)
        stateNew(k,12) = C(3,2)
        stateNew(k,13) = xOmega1
        stateNew(k,14) = xOmega2
        stateNew(k,15) = xOmega3
        stateNew(k,16) = xOmega4
        stateNew(k,17) = 0.0d0
        stateNew(k,18) = 0.0d0
        stateNew(k,21) = FI_LT
        stateNew(k,22) = FI_LC
        stateNew(k,23) = FI_MT
        stateNew(k,24) = FI_MC

        aminDamage = min(xOmega1,xOmega2,xOmega3,xOmega4)
        ! aminDamage = min(xOmega1,xOmega2)
        stateNew(k,19) = aminDamage

        IF (aminDamage.gt.0.99) stateNew(k,20) = 0.d0

      RETURN
      END SUBROUTINE cdm

  
!----------------------------------------------------------------------!
! Function rotate_stress: rotates stresses into misaligned             !
! coordinates system                                                   !
!----------------------------------------------------------------------!

      SUBROUTINE rotate_stress(trialStress, phiC, XC, g12,trialStressP)
        IMPLICIT NONE
        ! input variables
        REAL*8, INTENT(IN) :: phiC,XC,g12
        REAL*8, DIMENSION(6), INTENT(IN) :: trialStress
        ! local variables
        REAL*8 X,m,n,u,v,gamma0,phi,theta
        REAL*8 trialStressT(6), gammaMC, gammaM, phi0
        REAL*8, PARAMETER :: eps=1.d-8
        LOGICAL tS4EQ0, tS6EQ0
        ! output variables
        REAL*8, DIMENSION(6), INTENT(OUT) :: trialStressP

        ! first determine G plane angle theta (fiber kinking)
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

        gammaMC = (sin(2*phiC)*XC)/(2*g12)  ! eq 74
        phi0 = phiC - gammaMC  ! eq 75
        gammaM = ((phi0*g12 + abs(trialStressT(4)))/
     1           (g12+trialStressT(1)-trialStressT(2)) - phi0)  ! eq 81
        IF (trialStress(4).ge.0.d0) THEN ! eq 77
            phi = phi0 + gammaM
        ELSE
            phi = -1.0*(phi0+gammaM)
        ENDIF

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
      END SUBROUTINE rotate_stress


!----------------------------------------------------------------------!
! Function triangular: implements triangular damage dissipation        !
!                                                                      !
!----------------------------------------------------------------------!

      REAL*8 FUNCTION triangular(X,E,G,tStrain,l_ch)
        IMPLICIT NONE
          ! input variables
          REAL*8, INTENT(IN) :: X,E,G,tStrain,l_ch
          ! local variables
          REAL*8 delta_0, delta_c,dam

          IF (tStrain.eq.0.0d0) THEN ! no strain = no damage (avoid div. by 0)
            dam = 0.0d0 
          ELSE
            delta_0 = X/E
            delta_c = (2.0d0*G)/(X*l_ch)
            IF (delta_c.lt.delta_0) THEN  !Needed for d.nq.0-0
                delta_c = 1.1*delta_0   !Adjustable for stability
            ENDIF
            dam = (delta_c*(abs(tStrain)-delta_0))
     1            /(abs(tStrain)*(delta_c-delta_0))
            IF (dam.lt.0.0d0) THEN
                dam = 0.0d0
            ENDIF
          END IF
          triangular = min(0.999,dam)
      END FUNCTION triangular


!----------------------------------------------------------------------!
! Function mccauley: Implements the McCauley operator                  !
!                                                                      !
!----------------------------------------------------------------------!

      REAL*8 FUNCTION mccauley(value)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: value 
        mccauley = (value+abs(value))/2.0d0
        RETURN
      END FUNCTION mccauley

!----------------------------------------------------------------------!
! Subroutine FAIl_CLN: Catalanotti failure criteria                    !
!                                                                      !
!----------------------------------------------------------------------!

      SUBROUTINE fail_cln(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                    FI_MT,FI_MC)
        IMPLICIT NONE
        ! input variables
        REAL*8, DIMENSION(6), INTENT(IN) :: trialStress
        REAL*8, INTENT(IN) :: ST,SL,etaL,etaT,lambda,kappa
        ! local variables
        INTEGER, PARAMETER :: j = 51
        REAL*8, DIMENSION(j) :: b 
        REAL*8 a(0:j-1)
        REAL*8 trialFI_MT,trialFI_MC,aFail_MC,aFail_MT,pi,tN,tT,tL
        INTEGER i, p
        ! output variables
        REAL*8, INTENT(OUT) :: FI_MT, FI_MC
       
        pi = 4*atan(1.0d0) ! determine value of pi

        DO i =0,j-1
            a(i) = i 
        END DO

        ! a = (/ (i, i=0,j-1) /) ! create array of integers from 0 to 30
 
        DO p = 1,j
            b(p) = a(p)*(pi/(j-1.0d0)) ! create angles from 0 to pi
        END DO
       
        FI_MT = 0.0d0 ! initialize failure criteria
        FI_MC = 0.0d0
        DO p = 1,j ! iterate over angles
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
      END SUBROUTINE fail_cln



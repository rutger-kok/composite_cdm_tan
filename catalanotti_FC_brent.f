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
      REAL*8, DIMENSION(9) :: C
      REAL*8, DIMENSION(2) :: A,A_init
      REAL*8, DIMENSION(nblock,3) :: ATable
      REAL*8, DIMENSION(nblock,nstatev) :: sN,sO
      REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23 ! elastic constants
      REAL*8 nu21,nu31,nu32,delta
      REAL*8 XT,XC,YT,YC,SL,alpha0,stressPower
      REAL*8 G1plus,G1minus,G2plus,G2minus,G6 ! fracture energies
      REAL*8, PARAMETER :: tol=1.0d-5 ! for Brent method convergence
      INTEGER k,k1,i,j,f
      OPEN(105,file='C:\Workspace\catalanotti_final\test.txt',status='unknown',
     1       action='readwrite',position='append')
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

      delta = (1.0-nu12*nu21-nu23*nu32-nu13*nu31-2.0*nu21*nu32*nu13)
     1        /(e11*e22*e33)

      C(1) = (1-nu23*nu32)/(e22*e33*delta) !d11
      C(2) = (1-nu13*nu31)/(e11*e33*delta) !d22
      C(3) = (1-nu12*nu21)/(e11*e22*delta) !d33
      C(4) = (nu21+nu23*nu31)/(e22*e33*delta) !d12
      C(5) = (nu31+nu21*nu32)/(e22*e33*delta) !d13
      C(6) = (nu32+nu12*nu31)/(e11*e33*delta) !d23
      C(7) = 2.0*g12 !d44
      C(8) = 2.0*g23 !d55
      C(9) = 2.0*g13 !d66

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

          trialStress(1) = C(1)*trialStrain(1)+C(4)*trialStrain(2)
     1                     +C(5)*trialStrain(3)
          trialStress(2) = C(4)*trialStrain(1)+C(2)*trialStrain(2)
     1                     +C(6)*trialStrain(3)
          trialStress(3) = C(5)*trialStrain(1)+C(6)*trialStrain(2)
     1                     +C(3)*trialStrain(3)
          trialStress(4) = C(7)*trialStrain(4)
          trialStress(5) = C(8)*trialStrain(5)
          trialStress(6) = C(9)*trialStrain(6)

          DO i = 1,6
              stressNew(k,i)=trialStress(i)
          ENDDO

          DO k1 = 1,nstatev
            sO(k,k1) = stateNew(k,k1) ! copy state vars
          ENDDO  

          ! Calculate A Parameters using Brent method
          ! Initial approximations
          A_init(1) = 2.0d0*charLength(k)*XT*XT/
     1              (2.0d0*e11*G1plus-charLength(k)*XT*XT)
          A_init(2) = 2.0d0*charLength(k)*YT*YT/
     1              (2.0d0*e22*G2plus-charLength(k)*YT*YT)

          j = 1 ! sets case (case=1 to determine A1Plus)
          A(1) = brent(diff,A_init(1),tol,j)
        !   PRINT*, 'A1Plus = ', A(1)

          DO k1 = 1,nstatev
            sN(k,k1) = stateNew(k,k1) ! reset state vars
          ENDDO  

          j = 2 ! sets case (case=2 to determine A2Plus)
          A(2) = brent(diff,A_init(2),tol,j)
        !   PRINT*, 'A2Plus = ', A(2)

        WRITE(105,*) k, A(1), A(2)

        ENDDO

      ELSE

        REWIND(105)
        READ(105,*) ATable

        DO k = 1,nblock

          DO i = 1,6
              stateNew(k,i)=stateOld(k,i)+strainInc(k,i)
          ENDDO

          DO i = 1,6
              trialStrain(i)=stateNew(k,i)
          ENDDO

          stateOld(k,17) = ATable(k,2)
          stateOld(k,18) = ATable(k,3)

          call cdm(k,nblock,nstatev,trialStrain,stateOld,charLength(k),C,
     1               e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2               XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,G2minus,G6,
     3               alpha0,trialStress,stateNew)

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
      CLOSE(105)

      CONTAINS

!----------------------------------------------------------------------!
! Function brent: Brent method (root finding algorithm) to determine   !
! the A parameter for correct energy dissipation.                      !
!----------------------------------------------------------------------!

        REAL*8 FUNCTION brent(func,A_init,tol,case)
          IMPLICIT NONE
          !input variables
          REAL*8, INTENT(IN) :: tol, A_init!tolerance
          INTEGER, INTENT(IN) :: case !to determine A1Plus or A2Plus
          ! local variables
          INTEGER itmax,iter
          REAL*8  a,b,cc,dd,e,fa,fb,fc,p,q,r,s,tol1,xm,eps
          PARAMETER (itmax=100,eps=3.e-8)
          REAL*8, EXTERNAL :: func

          ! Brent's method to find the root of a function between x1 and x2
          ! Parameters: Maximum allowed number of iterations, machine 
          ! floating-point precision.
          ! Adapted from Numerical Recipes in FORTRAN 77 - Press et al.

          a = A_init/10000
          b = A_init*10000
          fa = func(a,case)
          fb = func(b,case)
          IF ((fa>0..AND.fb>0.).OR.(fa<0..AND.fb<0.)) THEN
              PRINT*, 'f(x1) = ',fa,'f(x2) = ',fb
              PRINT*, charlength(k), case
              PAUSE 'Root must be bracketed for BRENT'
              call XPLB_EXIT
          END IF
          cc = b
          fc = fb
          DO iter = 1,itmax
            IF((fb>0..AND.fc>0.).OR.(fb<0..AND.fc<0.)) THEN
              cc = a !Rename a, b, cc and adjust bounding interval dd.
              fc = fa
              dd = b-a
              e = dd
            END IF
            IF (abs(fc)<abs(fb)) THEN
              a = b
              b = cc
              cc = a
              fa = fb
              fb = fc
              fc = fa
            END IF
            tol1 = 2.*eps*abs(b)+0.5*tol !Convergence check.
            xm = .5*(cc-b)
            IF (abs(xm)<=tol1.OR.fb==0.) THEN
              brent = b
              RETURN
            END IF
            IF (abs(e)>=tol1.AND.abs(fa)>abs(fb)) THEN
              s = fb/fa !Attempt inverse quadratic interpolation.
              IF (a==cc) THEN
                p = 2.*xm*s
                q = 1.-s
              ELSE
                q = fa/fc
                r = fb/fc
                p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
                q = (q-1.)*(r-1.)*(s-1.)
              END IF
              IF (p>0.) q = -q !Check whether in bounds.
              p = abs(p)
              IF (2.*p < min(3.*xm*q-abs(tol1*q),abs(e*q))) THEN
                e = dd !Accept interpolation.
                dd = p/q
              ELSE
                dd = xm !Interpolation failed, use bisection.
                e = dd
              END IF
            ELSE !Bounds decreasing too slowly, use bisection.
              dd = xm
              e = dd
            END IF
            a = b !Move last best guess to a.
            fa = fb
            IF (abs(dd)>tol1) THEN !Evaluate new trial root.
              b = b+dd
            ELSE
              b = b+sign(tol1,xm)
            END IF
            fb = func(b,case)
          END DO
          brent = b
          PRINT*, 'BRENT exceeding maximum iterations'
          call XPLB_EXIT ! Exit analysis if Brent does not converge
        END FUNCTION brent


  !----------------------------------------------------------------------!
  ! Function diff: Difference between energy dissipated in element and   !
  ! material fracture energy                                             !
  !----------------------------------------------------------------------!

        REAL*8 FUNCTION diff(A_current,case)
          IMPLICIT NONE
          !input variables
          INTEGER, INTENT(IN) :: case
          REAL*8, INTENT(IN) :: A_current
          !local variables 
          REAL*8 gM,G_M
          IF (case.eq.1) THEN
            G_M = G1plus
          ELSE IF (case.eq.2) THEN
            G_M = G2plus
          END IF
          gM = integrate_gm(A_current,case)
          diff = gM - (G_M/charLength(k))
        END FUNCTION diff

  !----------------------------------------------------------------------!
  ! Function trapz: integrate array numerically using the trapezoid rule !
  !                                                                      !
  !----------------------------------------------------------------------!

        REAL*8 FUNCTION trapz(xarray,yarray)
          IMPLICIT NONE        
          ! input variables
          REAL*8, DIMENSION(1000), INTENT(IN) :: xarray, yarray
          ! local variables
          INTEGER i
          REAL*8 s,x1,x2,y1,y2,ddd,intg
          
          intg = 0.d0
          DO i = 1,999
              x1 = xarray(i)
              x2 = xarray(i+1)
              y1 = yarray(i)
              y2 = yarray(i+1)
              ddd = x2 - x1
              s=0.5*ddd*(y1+y2)
              intg = intg + s
          END DO
          trapz = intg
        END FUNCTION trapz

  !----------------------------------------------------------------------!
  ! Function integrate_gm: determine dissipated energy by integrating    !
  !                                                                      !
  !----------------------------------------------------------------------!

        REAL*8 FUNCTION integrate_gm(A_current,case)
          IMPLICIT NONE
          ! input variables
          INTEGER, INTENT(IN) :: case
          REAL*8, INTENT(IN) :: A_current
          ! local variables
          INTEGER, PARAMETER :: n=1000
          REAL*8, DIMENSION(6) :: tStrain,tStress
          REAL*8, DIMENSION(1000) :: tStrainArray,tStressArray
          REAL*8, DIMENSION(2) :: A_test
          REAL*8 m,failStress,failStrain,prevStress
          INTEGER i,counter
          m = 50.0d0
          counter = 0
          ! initialize strains
          tStrain(1) = 0.0d0
          tStrain(2) = 0.0d0
          tStrain(3) = 0.0d0
          tStrain(4) = 0.0d0
          tStrain(5) = 0.0d0
          tStrain(6) = 0.0d0
          IF (case.eq.1) THEN          
            ! define strain at initial failure
            tStrain(1) = XT/e11
            prevStress = XT
            sO(k,17) = A_current
            sO(k,18) = A_init(2)
            ! determine strain at final failure (when stress < XT/50)
            DO 
              IF (counter<1000) THEN
                call cdm(k,nblock,nstatev,trialStrain,stateOld,
     1               charlength(k),C,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2               XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,G2minus,G6,
     3               alpha0,trialStress,stateNew)
                IF (trialStress(1)<prevStress) THEN
                  tStrain(1) = tStrain(1) + 0.00025
                  counter = counter + 1
                  prevStress = trialStress(1)
                ELSE
                  EXIT
                END IF
              ELSE 
                PRINT*, 'Max failStress iterations exceeded',trialStress
     1                   ,tStrain
                call XPLB_EXIT
              END IF
            END DO
            failStrain = tStrain(case)
          ELSE IF (case.eq.2) THEN
            tStrain(2) = YT/e22
            prevStress = YT
            sO(k,17) = A_init(1)
            sO(k,18) = A_current
            ! determine strain at final failure (when stress < YT/50)
            DO 
              IF (counter<10000) THEN
                call cdm(k,nblock,nstatev,tStrain,stateOld,
     1               charLength(k),C,e11,e22,e33,nu12,nu13,nu23,g12,g13,
     2               g23,A_test,XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,
     3               G2minus,G6,alpha0,trialStress,stateNew)
                IF (trialStress(2)>failStress) THEN
                  tStrain(2) = tStrain(2) + 0.0001
                  counter = counter + 1
                  prevStress = trialStress(2)
                ELSE
                  EXIT
                END IF
              ELSE 
                PRINT*, 'Max failStress iterations exceeded',trialStress
     1                   ,tStrain
                call XPLB_EXIT
              END IF
            END DO
            failStrain = tStrain(case)
          END IF

          ! iterate over strain values from 0 to failure strain
          tStrain(1) = 0.0d0
          tStrain(2) = 0.0d0
          tStrain(3) = 0.0d0
          tStrain(4) = 0.0d0
          tStrain(5) = 0.0d0
          tStrain(6) = 0.0d0
          DO i = 0,n-1
            tStrain(case) = (failStrain/n)*i
            tStrainArray(i+1) = tStrain(case)    
            call cdm(k,nblock,nstatev,tStrain,sO,
     1               charLength(k),C,e11,e22,e33,nu12,nu13,nu23,g12,g13,
     2               g23,XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,
     3               G2minus,G6,alpha0,tStress,sN)
            tStressArray(i+1) = tStress(case)

          END DO
          
          ! integrate stress strain curve to obtain dissipated energy
          integrate_gm = trapz(tStrainArray,tStressArray)

        END FUNCTION integrate_gm

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
        REAL*8, DIMENSION(9), INTENT(INOUT) :: C
        REAL*8, DIMENSION(6), INTENT(IN) :: trialStrain
        REAL*8, INTENT(IN) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
        REAL*8, INTENT(IN) :: XT,XC,YT,YC,SL,l_ch
        REAL*8, INTENT(IN) :: G1plus,G1minus,G2plus,G2minus,G6
        ! local variables
        REAL*8, DIMENSION(6) :: trialStressP
        REAL*8 d1Plus,d1Minus,d2Plus,d2Minus,d6
        REAL*8 rNew1,rNew2,rNew3,rNew4,rNew5,rNew6
        REAL*8 rOld1,rOld2,rOld3,rOld4,rOld5,rOld6
        REAL*8 nu21,nu31,nu32,delta,aminDamage
        REAL*8 alpha0,etaT,etaL,phiC,kappa,lambda,ST
        REAL*8 FI_LT,FI_LC,FI_MT,FI_MC,A1Plus,A2Plus
        REAL*8 xOmega1,xOmega2,xOmega3,xOmega4,xOmega5,xOmega6
        ! output variables
        REAL*8, DIMENSION(6), INTENT(OUT) :: trialStress
        REAL*8, DIMENSION(nblock,nstatev), INTENT(INOUT) :: stateNew
        ! external functions
        REAL*8, EXTERNAL :: triangular

        ! Initial values
        etaL = 0.280622004043348
        phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)*((SL/XC)+etaL))) 
     1         /(2.0d0*((SL/XC)+etaL)))  !Eq.12
        ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*SL)/  
     1       (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL)) !Eq.12 CLN
        etaT = (etaL*ST)/SL  !Eq.10 CLN
        kappa = (ST**2.0d0-YT**2.0)/(ST*YT)  !Eq.43 CLN
        lambda = ((2.0*etaL*ST)/SL)-kappa  !Eq.45 CLN 

        ! Update of the failure thresholds (r values)
        rOld1 = stateOld(k,7)
        rOld2 = stateOld(k,8)
        rOld3 = stateOld(k,9)
        rOld4 = stateOld(k,10)
        rOld5 = stateOld(k,11)
        rOld6 = stateOld(k,12)
        A1Plus = stateOld(k,17)
        A2Plus = stateOld(k,18)

        ! IF ((A1Plus.lt.0.5).or.(A2Plus.lt.0.3)) THEN
        !     PRINT*, A1Plus,A2Plus
        !     call XPLB_EXIT ! Exit analysis if Brent does not converge
        ! END IF

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
            call rotate_stress(trialStress,phiC,XC,g12,trialStressP)
            call fail_cln(trialStressP,ST,SL,etaL,etaT,lambda,kappa,
     1                    FI_LT,FI_LC)
        ENDIF

        ! transverse failure criteria
        call fail_cln(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                FI_MT,FI_MC)

        ! Update of the damage thresholds
        rNew1 = max(1.0,max(FI_LT,rOld1),max(FI_LC,rOld2))    !Eq.26
        rNew2 = max(1.0,max(FI_LC,rOld2))
        rNew3 = max(1.0,max(FI_MT,rOld3),max(FI_MC,rOld4))
        rNew4 = max(1.0,max(FI_MC,rOld4))
        rNew5 = 1.0d0
        rNew6 = 1.0d0
 
        !Eq.5 (from part II of the Maimi paper)
        d1Plus = 1.0d0-1.0d0/rNew1*exp(A1Plus*(1.0d0-rNew1)) 

        IF (rNew2.gt.1.0d0) THEN !Triangular damage dissipation
          d1Minus = triangular(XC,e11,G1minus,trialStrain(1),l_ch)
        ELSE
          d1Minus = 0.0d0
        ENDIF

        d2Plus = 1.0d0-1.0d0/rNew3*exp(A2Plus*(1.0d0-rNew3))

        IF (rNew4.gt.1.0d0) THEN
          d2Minus = triangular(YC,e22,G2minus,trialStrain(2),l_ch)
        ELSE
          d2Minus = 0.0d0
        ENDIF

        IF (rNew3.gt.1.0d0) THEN !Triangular damage dissipation
          d6 = triangular(SL,e22,G6,trialStrain(4),l_ch)
        ELSE
          d6 = 0.0d0
        ENDIF

        ! Damage variables

        IF (trialStress(1).gt.0) THEN      !Eq.6
            xOmega1 = min(0.999,d1Plus)
        ELSE
            xOmega1 = d1Minus
        ENDIF

        IF (trialStress(2).gt.0) THEN
            xOmega2 = min(0.999, d2Plus)
        ELSE
            xOmega2 = d2Minus
        ENDIF

        xOmega4 = min(0.999,d6)
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


        ! IF (d1Plus.gt.0.935d0) THEN
        !     PRINT*, rNew1
        !     PRINT*, A1Plus,A2Plus
        !     PRINT*, d1Plus
        !     print*, trialStrain(1)
        !     print*, trialStress(1)
        !     print*, '-----------------------'
        ! ENDIF 

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
        stateNew(k,17) = A1Plus
        stateNew(k,18) = A2Plus
        stateNew(k,21) = FI_LT
        stateNew(k,22) = FI_LC
        stateNew(k,23) = FI_MT
        stateNew(k,24) = FI_MC

        aminDamage = min(xOmega1,xOmega2,xOmega4)
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

          delta_0 = X/E
          delta_c = (2.0d0*G/E*l_ch)
          IF (delta_c.lt.delta_0) THEN  !Needed for d.nq.0-0
              delta_c = 1.1*delta_0   !Adjustable for stability
          ENDIF
          dam = (delta_c*(abs(tStrain)-delta_0))
     1              /(abs(tStrain)*(delta_c-delta_0))
          IF (dam.lt.0.0d0) THEN
              dam = 0.0d0
          ENDIF
          triangular = min(0.999,dam)
      END FUNCTION triangular


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
        INTEGER, PARAMETER :: j = 31
        REAL*8, DIMENSION(j) :: a,b 
        REAL*8 trialFI_MT,trialFI_MC,aFail_MC,aFail_MT,pi,tN,tT,tL
        INTEGER i, p
        ! output variables
        REAL*8, INTENT(OUT) :: FI_MT, FI_MC
       
        pi = 4*atan(1.0d0) ! determine value of pi
        a = (/ (i, i=0,j-1) /) ! create array of integers from 0 to 30
 
        DO p = 1,j
            b(p) = a(p)*(pi/((j-1)*2)) ! create angles from 0 to pi/2
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


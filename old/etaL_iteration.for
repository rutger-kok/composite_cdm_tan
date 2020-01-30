C----------------------------------------------------------------------C
C     A Parameter Calculation Program                                  C
C     Author(s): Rutger Kok                                            C
C     Date: 06/08/2019                                                 C
C     Version 1.0                                                      C
C----------------------------------------------------------------------C

      PROGRAM etaL_iteration
      IMPLICIT NONE
      
C Main program variables and external routines

      REAL*8 etaLFinal,tol,BRENT
      EXTERNAL f,BRENT,FAIL_CLN,ROTATE_PHI

C Declare common block for material properties

      REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
      REAL*8 XT,XC,YT,YC,SL,ST
      REAL*8 G1plus,G1minus,G2plus,G2minus,G6
      REAL*8 alpha0,etaT,etaL0,phiC,kappa,lambda
      COMMON /mat_props/ e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     1                   XT,XC,YT,YC,SL,ST,
     2                   G1plus,G1minus,G2plus,G2minus,G6,
     3                   alpha0,etaT,phiC,kappa,lambda
      REAL*8 lch
      COMMON /char_length/ lch

C Open file to write stress strain data to (for plotting curves)   
      
C      OPEN(1, file = 
C     1 'C:\Workspace\A_parameter.txt',
C     2 status = 'unknown')
C200   FORMAT (A,A,A,A,A,A,A,A,A)
C      WRITE(1,200) 'trialStress(2)', ',','trialStrain(2)', ',','FI_MT',
C     1             ',', 'r2Plus', ',', 'd2Plus'  

C --------------------------------------------------------------------
C Define material properties

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

C Energies

        G1Plus = 0.1 !G1c+
        G1Minus = 0.1 !G1c-
        G2Plus = 0.00075 !G2c+
        G2Minus = 0.0025 !G2c-
        G6 = 0.0035 !G6

C --------------------------------------------------------------------
C Start main program

C Initial values

C Initial approximation etaL
      etaL0 = 0.280  
      tol = 0.000001 ! tolerance for Brent method convergence

      lch = 1.0
      phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)*((SL/XC)+etaL0)))
     1       /(2.0d0*((SL/XC)+etaL0))) !Eq.12
      ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*SL)/  !Eq.12 (CLN)
     1     (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL0))
      etaT = (etaL0*ST)/SL                      !Eq.10 CLN
      kappa = (ST**2.0d0-YT**2.0)/(ST*YT)  !Eq.43 CLN
      lambda = ((2.0*etaL0*ST)/SL)-kappa  !Eq.45 CLN 
          
      etaLFinal = BRENT(f,etaL0,tol) ! iteratively find etaL parameter
      
C      CLOSE(1)

      PRINT*, 'etaL Final:', etaLFinal
      PAUSE "Press RETURN to end program"
      END PROGRAM etaL_iteration

* <<<<<<<<<<<<<<<<<<<<<<<< FUNCTION BRENT >>>>>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Brent method to find roots of a function f(x)
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      FUNCTION f(etaLTest)
          IMPLICIT NONE
          REAL*8 f,gM,G_M,lch,etaLTest,trialStress(6),FI_MC,FI_MT
          REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
          REAL*8 XT,XC,YT,YC,SL,ST
          REAL*8 G1plus,G1minus,G2plus,G2minus,G6
          REAL*8 alpha0,etaT,etaL,phiC,kappa,lambda
          COMMON /char_length/ lch
          COMMON /mat_props/ e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     1                     XT,XC,YT,YC,SL,ST,
     2                     G1plus,G1minus,G2plus,G2minus,G6,
     3                     alpha0,etaT,phiC,kappa,lambda
          EXTERNAL FAIL_CLN

          lch = 1.0
          phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)
     1           *((SL/XC)+etaLTest)))/(2.0d0*((SL/XC)+etaLTest)))
          ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*SL)/  !Eq.12 (CLN)
     1         (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaLTest))
          etaT = (etaLTest*ST)/SL                      !Eq.10 CLN
          kappa = (ST**2.0d0-YT**2.0)/(ST*YT)  !Eq.43 CLN
          lambda = ((2.0*etaLTest*ST)/SL)-kappa  !Eq.45 CLN 

          trialStress(1) = 0.0
          trialStress(2) = -0.285
          trialStress(3) = 0.0
          trialStress(4) = 0.0
          trialStress(5) = 0.0
          trialStress(6) = 0.0

          call FAIL_CLN(trialStress,ST,YT,SL,etaLTest,etaT,lambda,kappa,
     1                  FI_MT,FI_MC)
          f = FI_MC - 1.0
          RETURN
      END
        
      FUNCTION BRENT(func,etaL0,tol)
        IMPLICIT NONE
        INTEGER ITMAX
        REAL*8 BRENT,tol,x1,x2,func,EPS,etaL0
        EXTERNAL func
        PARAMETER (ITMAX=100,EPS=3.e-8)
        ! Brent's method to find the root of a function between x1 and x2
        ! Parameters: Maximum allowed number of iterations, and machine 
        ! floating-point precision.
        ! Adapted from Numerical Recipes in FORTRAN 77 - Press et al.
        INTEGER iter
        REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        
        x1 = etaL0/10.0d0 !initial etaL variable
        x2 = etaL0*10.0d0
        
        a=x1
        b=x2
        fa=func(a) !passes x1 to function f (above)
        fb=func(b) !passes x2 to function f (above)
        IF ((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) THEN
            PRINT*, 'f(x1)=',fa,'f(x2)=',fb
            PAUSE 'Root must be bracketed for BRENT'
        END IF
        c=b
        fc=fb
        DO iter=1,ITMAX
          IF((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))THEN
            c=a !Rename a, b, c and adjust bounding interval d.
            fc=fa
            d=b-a
            e=d
          END IF
          IF (abs(fc).lt.abs(fb)) THEN
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
          END IF
          tol1=2.*EPS*abs(b)+0.5*tol !Convergence check.
          xm=.5*(c-b)
          IF (abs(xm).le.tol1 .or. fb.eq.0.) THEN
            BRENT=b
            RETURN
          END IF
          IF (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) THEN
            s=fb/fa !Attempt inverse quadratic interpolation.
            IF (a.eq.c) THEN
              p=2.*xm*s
              q=1.-s
            ELSE
              q=fa/fc
              r=fb/fc
              p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
              q=(q-1.)*(r-1.)*(s-1.)
            END IF
            IF (p.gt.0.) q=-q !Check whether in bounds.
            p=abs(p)
            IF (2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) THEN
              e=d !Accept interpolation.
              d=p/q
            ELSE
              d=xm !Interpolation failed, use bisection.
              e=d
            END IF
          ELSE !Bounds decreasing too slowly, use bisection.
            d=xm
            e=d
          END IF
          a=b !Move last best guess to a.
          fa=fb
          IF (abs(d).gt.tol1) THEN !Evaluate new trial root.
            b=b+d
          ELSE
            b=b+sign(tol1,xm)
          END IF
          fb=func(b)
        END DO
        PAUSE 'BRENT exceeding maximum iterations'
        BRENT=b
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
      SUBROUTINE ROTATE_PHI(trialStress, phiC, XC, g12, trialStressP)
        IMPLICIT NONE
        ! input variables
        REAL*8 trialStress(6), phiC, XC, NEWT, eps, g12
        ! local variables
        REAL*8 X, m, n, u, v, gamma0
        REAL*8 trialStressT(6), gammaMC, gammaM, phi0
        PARAMETER (eps=1.d-8)
        LOGICAL tS4EQ0, tS6EQ0
        ! output variables
        REAL*8 trialStressP(6), theta, phi

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
        END


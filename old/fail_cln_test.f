      PROGRAM fail_cln_test
      IMPLICIT NONE

      REAL*8 trialStress(6),XT,XC,YT,YC,SL,ST,etaL,alpha0,etaT
      REAL*8 kappa, lambda, FI_MT, FI_MC

      !IM78552 taken from Catalanotti (2013)
      XT = 2323.5
      XC = 1200.1
      YT = 160.2
      YC = 198.0
      SL = 130.2
      etaL = 0.5
      alpha0 = 0.9250245036

      ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*SL)
     1     /(((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL))  ! eq 12
      etaT = (etaL*ST)/SL  ! eq 10
      kappa = (ST**2.0-YT**2.0)/(ST*YT)  ! eq 43
      lambda = ((2.0*etaL*ST)/SL)-kappa  ! eq 45 

      !define trial stresses
      trialStress(1) = 0.0
      trialStress(2) = -130.0
      trialStress(3) = 0.0
      trialStress(4) = 120.0
      trialStress(5) = 0.0
      trialStress(6) = 0.0
  
      call FAIL_CLN(trialStress,ST,YT,SL,etaL,etaT,lambda,kappa,
     1              FI_MT,FI_MC)

      PRINT*,'FI_MT = ', FI_MT
      PRINT*,'FI_MC = ', FI_MC
      PAUSE "Press RETURN to end program"
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
      
C-----------------------------------------------------------------------C
C     ABAQUS VUMAT USER SUBROUTINE: catalanotti_criteria.for            C
C	  Author(s): Rutger Kok, Francisca Martinez-Hergueta                 C
C     Date: 24/06/2019                                           	    C
C     Version 1.0                                                       C
C-----------------------------------------------------------------------C

C
C User subroutine VUMAT
      subroutine vumat (
C Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
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
C     
      real*8 trialStress(6),trialStrain(6)
      real*8 xomega(6)
      real*8 d1Plus,d1Minus,d2Plus,d2Minus,d6
      real*8 rNew1,rNew2,rNew3,rNew4,rNew5,rNew6
      real*8 dummy, dummysalida
      real*8 e11,e22,e33,xnu12,xnu13,xnu23,xmu12,xmu13,xmu23
      real*8 XT,XC,YT,YC,SL
      real*8 G1plus,G1minus,G2plus,G2minus,G6
      real*8 A1plus,A1minus,A2plus,A2minus,A6
      real*8 alpha0,etaT,etaL,phic,ST,theta,omegaValue
      real*8 phiLT,phiLC_MC, phiLC_MT, phiLC,phiMT,phiMC
      real*8 phiMC2,phiMC3, phiMT2, phiMT3, UNUSED 
      real*8 initialcharLength,thickness,traceStrain,meanDamage,expo, A
      real*8 trialStressP(6), trialStressT(6), phi, X, NEWT
      real*8 FUNC, DFUNC
      character*80 cmname
      parameter ( zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.0d0,
     *     third = 1.d0 / 3.d0, half = 0.5d0, six = 6.0d0)


C Elastic constants orthotropic ply

      e11   = props(1)
      e22   = props(2)
      e33   = props(3)
      xnu12 = props(4)
      xnu13 = props(5)
      xnu23 = props(6)
      xmu12 = props(7)
      xmu13 = props(8)
      xmu23 = props(9)

C Ply strength

      XT=props(10)
      XC=props(11)
      YT=props(12)
      YC=props(13)
      SL=props(14)

C Fracture Angle

      alpha0=props(15)*0.017453292519943295

C Fracture toughness

      G1plus=props(16)
      G1minus=props(17)
      G2plus=props(18)
      G2minus=props(19)
      G6=props(20)

C Initial values

      etaL=-1.0d0*SL*cos(2.0d0*alpha0)/(YC*cos(alpha0)**2.0d0)     !Eq.10

      dummy=SL/XC
      phic=ATan((1.0d0-Sqrt(1.0d0-4.0d0*dummy*(dummy+etaL)))         !Eq.12
     +/(2.0d0*(dummy+etaL)))

      ST = (0.5*(((2*sin(alpha0)**2.0)-1.0)*Sl_is)/    !Eq 12 Catalanotti (CLN)
     +   (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL))
      etaT=-1.0d0/tan(2.0d0*alpha0)                      !Eq.17

      omegaValue=-0.120382            !Eq.20

C Stiffness matrix orthotropic material

      xnu21=xnu12*e22/e11
      xnu31=xnu13*e33/e11
      xnu32=xnu23*e33/e22

      xnu=1.0/(1.0-xnu12*xnu21-xnu32*xnu23-xnu13*xnu31
     *-2.0*xnu21*xnu13*xnu32)

      d11=e11*(1.0-xnu23*xnu32)*xnu
      d22=e22*(1.0-xnu13*xnu31)*xnu 
      d33=e33*(1.0-xnu12*xnu21)*xnu
      d12=e11*(xnu21+xnu31*xnu23)*xnu
      d13=e11*(xnu31+xnu21*xnu32)*xnu
      d23=e22*(xnu32+xnu12*xnu31)*xnu
      d44=2.0*xmu12
      d55=2.0*xmu23
      d66=2.0*xmu13   

C Loop through the gauss points

      if(stepTime.eq.0) then     

C Initial elastic step, for Abaqus tests

      do k = 1, nblock
      
C Initialisation of state variables
       do k1=1,nstatev
        stateNew(k,k1)=0.d0
       end do

       do i=1,6
       trialStrain(i)=strainInc(k,i)
       enddo


      trialStress(1)=d11*trialStrain(1)+d12*trialStrain(2)
     *+d13*trialStrain(3)
      trialStress(2)=d12*trialStrain(1)+d22*trialStrain(2)
     *+d23*trialStrain(3)
      trialStress(3)=d13*trialStrain(1)+d23*trialStrain(2)
     *+d33*trialStrain(3)
      trialStress(4)=d44*trialStrain(4)
      trialStress(5)=d55*trialStrain(5)
      trialStress(6)=d66*trialStrain(6)

      do i=1,6
      stressNew(k,i)=trialStress(i)
      enddo

      end do

C Constitutive model      
      else

       do k=1,nblock

C Update of the failure thresholds (r values)

      rOld1=stateOld(k,7)
      rOld2=stateOld(k,8)
      rOld3=stateOld(k,9)
      rOld4=stateOld(k,10)
      rOld5=stateOld(k,11)
      rOld6=stateOld(k,12)

C Computation of the total strain

      do i=1,6
      stateNew(k,i)=stateOld(k,i)+strainInc(k,i)
      enddo

       do i=1,6
        trialStrain(i)=stateNew(k,i)
       enddo

C Trial stress

      trialStress(1)=d11*trialStrain(1)+d12*trialStrain(2)
     *+d13*trialStrain(3)
      trialStress(2)=d12*trialStrain(1)+d22*trialStrain(2)
     *+d23*trialStrain(3)
      trialStress(3)=d13*trialStrain(1)+d23*trialStrain(2)
     *+d33*trialStrain(3)
      trialStress(4)=d44*trialStrain(4)
      trialStress(5)=d55*trialStrain(5)
      trialStress(6)=d66*trialStrain(6)

      do i=1,6
       stressNew(k,i)=trialStress(i)
      enddo

C Evaluation of the damage activation functions

c determine fracture plane angle theta (for fiber kinking)   
      
      IF ((trialStress(4).eq.0).AND.(trialStress(6).eq.0)) THEN
        theta = 0.5d0*ATAN((2.0d0*trialStress(5))
     *          /(trialStress(2)-trialStress(3))) ! Eq. 55 CLN
      ELSE
        theta = ATAN(trialStress(6)/trialStress(4)) ! Eq. 56 CLN
      END IF

c determine the misalignment angle phi
      X = (sin(2.0d0*phic)*XC)/(2.0d0*phic) ! Eq. 86 CLN

      IF (trialStress(4).gt.0) THEN
        phi = NEWT(0.1, trialStress, X) ! specified initial value of 0.1
      ELSE 
        phi = -1.0d0*NEWT(0.1, trialStress, X)
      END IF

c Rotate stresses by angle theta
      trialStressT(1) = trialStress(1)
      trialStressT(4) = trialStress(4)*cos(theta) 
     *                  + trialStress(6)*sin(theta)
      trialStressT(6) = trialStress(6)*cos(theta) 
     *                  - trialStress(4)*sin(theta)
      trialStressT(2) = trialStress(2)*cos(theta)**2 
     *                  + 2.0d0*trialStress(5)*cos(theta)*sin(theta) 
     *                  + trialStress(3)*sin(theta)**2
      trialStressT(5) = trialStress(5)*(cos(theta)**2 - sin(theta)**2)
     *                  - trialStress(2)*sin(theta)*cos(theta) 
     *                  + trialStress(3)*sin(theta)*cos(theta)
      trialStressT(3) = trialStress(3)*cos(theta)**2 
     *                  - 2.0d0*trialStress(5)*cos(theta)*sin(theta) 
     *                  + trialStress(2)*sin(theta)**2

c NOTE: phi is also defined as a local variable in the Softening
c parameter subroutine but this is not called in the main program
c Rotate stresses by angle phi
      trialStressP(1) = trialStressT(1)*cos(phi)**2 
     *                  + 2.0d0*trialStressT(4)*cos(phi)*sin(phi) 
     *                  + trialStressT(2)*sin(phi)**2
      trialStressP(4) = trialStressT(4)*(cos(phi)**2 -sin(phi)**2) 
     *                  + trialStressT(2)*sin(phi)*cos(phi) 
     *                  - trialStressT(1)*sin(phi)*cos(phi)
      trialStressP(6) = trialStressT(6)*cos(phi) 
     *                  + trialStressT(5)*sin(phi)
      trialStressP(2) = trialStressT(2)*cos(phi)**2 
     *                  - 2.0d0*trialStressT(4)*sin(phi)*cos(phi) 
     *                  + trialStressT(1)*sin(phi)**2 
      trialStressP(5) = trialStressT(5)*cos(phi) 
     *                  - trialStressT(6)*sin(phi)
      trialStressP(3) = trialStressT(3)

c call subroutine to determine longitudinal compressive failure criteria
      IF (trialStress(1).gt.0.d0) THEN
        phiLT = trialStrain(1)/(XT/e11) ! Eq. 54 CLN
      ELSEIF (trialStress(1).lt.0.d0) THEN
        call FAIL(trialStressP,ST,YT,SL,etaL,etaT, phiMC, phiMT)
        phiLC_MC = phiMC
        phiLC_MT = phiMC
        phiLC = max(phiLC_MC, phiLC_MT)
      ENDIF
    
C call subroutine to determine transverse compressive and tensile 
C failure criteria
C in the 2-direction
      IF (trialStress(2).ge.0.d0) THEN
        call FAIL(trialStress,ST,YT,SL,etaL,etaT, phiMC, phiMT)
        UNUSED = phiMC
        phiMT2 = phiMC
      ELSEIF (trialStress(2).lt.0.d0) THEN
        call FAIL(trialStress,ST,YT,SL,etaL,etaT, phiMC, phiMT)
        phiMC2 = phiMC
        UNUSED = phiMC
      ENDIF 

C in the 3-direction
      IF (trialStress(3).gt.0.d0) THEN
        call FAIL(trialStress,ST,YT,SL,etaL,etaT, phiMC, phiMT)
        UNUSED = phiMC
        phiMT3 = phiMC
      ELSEIF (trialStress(3).lt.0.d0) THEN
        call FAIL(trialStress,ST,YT,SL,etaL,etaT, phiMC, phiMT)
        phiMC3 = phiMC
        UNUSED = phiMC
      ENDIF

C transverse failure criteria is the max of both 2 and 3 direction
      phiMT = max(phiMT2,phiMT3)
      phiMC = max(phiMC2, phiMC3) 

C Update of the damage thresholds

       rNew1=max(1.0,max(phiLT,rOld1),max(phiK,rOld2))    !Eq.26
       rNew2=max(1.0,max(phiLC,rOld2))
       rNew3=max(1.0,max(phiMT,rOld3),max(phiMC,rOld4))
       rNew4=max(1.0,max(phiMC,rOld4))
       rNew5=1.0d0
       rNew6=1.0d0

C Softening parameter A

c	initialcharLength=charLength(k)
      initialcharLength=1.0
      
       A1plus=2.0d0*initialcharLength*XT*XT                    !Second part. eq. B.2
     +/(2.0d0*e11*G1plus-initialcharLength*XT*XT)
      
      A1minus=2.0d0*initialcharLength*(1.2671444971783234*XC)**2.0d0   ! Second part. Appendix B
     +/(2.0d0*(1.173191594009027*e11)*G1minus-initialcharLength
     +*(1.145996547107106*XC)**2.0d0)
      
       A2plus=2.0d0*initialcharLength*(0.8308920721648098*YT)**2.0d0
     +/(2.0d0*0.9159840495203727*e22*G2plus-initialcharLength
     +*(0.9266465446084761*YT)**2.0d0)

       A2minus=2.0d0*initialcharLength*(0.709543*YC)**2.0d0
     +/(2.0d0*0.7881*e22*G2minus-initialcharLength*(0.802572*YC)**2.0d0)

c      A2minus=5.0

      A6=2.0d0*initialcharLength*SL*SL
     +/(2.0d0*xmu12*G6-initialcharLength*SL*SL)
     
      A6=A2plus
      
c      write(6,*) "A1 coefficients", A1plus,A1minus
c      write(6,*) "A2 coefficients", A2plus,A2minus
c      write(6,*) "A6 coefficients", A6
 
      d1Plus=1.0d0-1.0d0/rNew1*dexp(A1plus*(1.0d0-rNew1))         !Eq.28
      d1Minus=1.0d0-1.0d0/rNew2*dexp(A1minus*(1.0d0-rNew2))
      d2Plus=1.0d0-1.0d0/rNew3*dexp(A2plus*(1.0d0-rNew3))
      d2Minus=1.0d0-1.0d0/rNew4*dexp(A2minus*(1.0d0-rNew4))
      d6=1.0d0-1.0d0/rNew3*dexp(A6*(1.0d0-rNew3))*(1.0d0-d1Plus)


C Damage variables

       if(trialStress(1).gt.0) then      !Eq.6
          xOmega1=d1Plus
       else
          xOmega1=d1Minus
       endif

       if(trialStress(2).gt.0) then
          xOmega2=d2Plus
       else
          xOmega2=d2Minus
       endif

c	xOmega2=0.d0

        xOmega4=d6
        xOmega3=0.0d0
        xOmega5=0.0d0
        xOmega6=0.0d0 

C Stiffness tensor + damage

        xnu21=xnu12*e22/e11
        xnu31=xnu13*e33/e11
        xnu32=xnu23*e33/e22
 
        xnu=1.0d0/(1.0d0-xnu12*xnu21*(1.0d0-xOmega1)*(1.0d0-xOmega2)
     *-xnu32*xnu23*(1.0d0-xOmega2)*(1.0d0-xOmega3)
     *-xnu13*xnu31*(1.0d0-xOmega1)*(1.0d0-xOmega3)
     *-2.0d0*xnu12*xnu23*xnu31*(1.0d0-xOmega1)
     **(1.0d0-xOmega2)*(1.0d0-xOmega3))

      d11=e11*(1.0d0-xOmega1)*(1.0d0-xnu23*xnu32*(1.0d0-xOmega2)*
     *           (1.0d0-xOmega3))*xnu
      d22=e22*(1.0d0-xOmega2)*(1.0d0-xnu13*xnu31*(1.0d0-xOmega1)*
     *           (1.0d0-xOmega3))*xnu 
      d33=e33*(1.0d0-xOmega3)*(1.0d0-xnu12*xnu21*(1.0d0-xOmega1)*
     *           (1.0d0-xOmega2))*xnu
      d12=e11*(1.0d0-xOmega1)*(1.0d0-xOmega2)*(xnu21+xnu31*xnu23*
     +           (1.0d0-xOmega3))*xnu
      d13=e11*(1.0d0-xOmega1)*(1.0d0-xOmega3)*(xnu31+xnu21*xnu32*
     +           (1.0d0-xOmega2))*xnu
      d23=e22*(1.0d0-xOmega2)*(1.0d0-xOmega3)*(xnu32+xnu12*xnu31*
     +           (1.0d0-xOmega1))*xnu
      d44=2.0d0*xmu12*(1.0d0-xOmega4)
      d55=2.0d0*xmu23*(1.0d0-xOmega5)
      d66=2.0d0*xmu13*(1.0d0-xOmega6)

      trialStress(1)=d11*trialStrain(1)+d12*trialStrain(2)
     *+d13*trialStrain(3)
      trialStress(2)=d12*trialStrain(1)+d22*trialStrain(2)
     *+d23*trialStrain(3)
      trialStress(3)=d13*trialStrain(1)+d23*trialStrain(2)
     *+d33*trialStrain(3)
      trialStress(4)=d44*trialStrain(4)
      trialStress(5)=d55*trialStrain(5)
      trialStress(6)=d66*trialStrain(6)

      do i=1,6
        stressNew(k,i)=trialStress(i)
      enddo

C Energy

      stressPower=half*(
     +(stressNew(k,1)+stressOld(k,1))*strainInc(k,1)+
     +(stressNew(k,2)+stressOld(k,2))*strainInc(k,2)+
     + (stressNew(k,3)+stressOld(k,3))*strainInc(k,3)+
     + 2.0*(stressNew(k,4)+stressOld(k,4))*strainInc(k,4)+
     + 2.0*(stressNew(k,5)+stressOld(k,5))*strainInc(k,5)+
     + 2.0*(stressNew(k,6)+stressOld(k,6))*strainInc(k,6))

      enerInternNew(k)=enerInternOld(k)+stressPower/density(k)

        stateNew(k,7)=rNew1
        stateNew(k,8)=rNew2
        stateNew(k,9)=rNew3
        stateNew(k,10)=rNew4
        stateNew(k,11)=rNew5
        stateNew(k,12)=rNew6
        stateNew(k,13)=xOmega1
        stateNew(k,14)=xOmega2
        stateNew(k,15)=xOmega3
        stateNew(k,16)=xOmega4
        stateNew(k,17)=xOmega5
        stateNew(k,18)=xOmega6

      aminDamage=min(xOmega1,xOmega2,xOmega4)
        stateNew(k,19)=aminDamage

       if(aminDamage.gt.0.999) stateNew(k,20)=0.d0

       enddo

       endif

c	end if !at the beginning
      
      return
        end


      subroutine McAulay(valor,salida)
        real*8 valor,salida
        if(valor.gt.0.0d0) then 
        salida=valor
        else
        salida=0.0d0
        endif
        return
        end


C Damage varable, triangular damage dissipation       
       subroutine Damage(X,E,G,trialStrain,lch,dNew)
        real*8 X,E,G,trialStrain,lch,dNew
        real*8 delta_0, delta_c
           delta_0 = X/E
           delta_c = (2.0d0*G/X*lch)
           if (delta_c.lt.delta_0)then       !!Needed for d.nq.0-0
                   delta_c = 1.1*delta_0     !! 1.1 adjustable for stability
           endif
           dNew = (delta_c*(abs(trialStrain)-delta_0))
     *            /(abs(trialStrain)*(delta_c-delta_0))
        return
       end


c Subroutine to determine failure in transverse tension and compression
c and in longitudinal compression (with rotated coord system)
       SUBROUTINE FAIL(trialStress,ST,YT,SL,etaL,etaT, phiMC, phiMT)
        IMPLICIT NONE
c input variables
        REAL*8 trialStress(6), ST, YT, SL, etaL, etaT
c local variables
        REAL*8 lambda, kappa, a(31), b(31), pi, tN, tT, tL
        REAL*8 trialphiMC, trialphiMT, aFailC, aFailT
        INTEGER i, p
c output variables
        REAL*8 phiMC, phiMT
        
        pi = 4*atan(1.0d0) ! determine value of pi
        a = (/ (i, i=0,30) /) ! create array of integers from 0 to 30
  
        DO p = 1,31
          b(p) = a(p)*(pi/60) ! create array of angles from 0 to pi/2
        END DO
        
        phiMC = 0.0d0 ! initialize max failure criteria
        phiMT = 0.0d0

        DO p = 1, 31 ! iterate over angles
            ! Eq 3 CLN (expanded in 59-61)
            tN = trialStress(2)*cos(b(p))**2 + 2.0d0*trialStress(5)
     *          *sin(b(p))*cos(b(p)) + trialStress(3)*sin(b(p))**2
            tT = -1.0d0*trialStress(2)*cos(b(p))*sin(b(p)) 
     *          + trialStress(3)*sin(b(p))*cos(b(p)) - trialStress(5)
     *          *(cos(b(p))**2 - sin(b(p))**2)
            tL = trialStress(4)*cos(b(p)) + trialStress(6)*sin(b(p))
      
            kappa = (ST**2 - YT**2)/(ST*YT)             ! Eq 43 CLN
            lambda = (2.0d0*etaL*ST)/SL - kappa         ! Eq 45 CLN
      
            trialphiMC = (tL/(SL-etaL*tN))**2 + (tT/(ST-etaT*tN))**2 ! Eq. 5 CLN
            trialphiMT = (tN/ST)**2 + (tL/SL)**2 + (tT/ST)**2 
     *          + lambda*(tN/ST)*(tL/SL)**2 + kappa*(tN/ST) ! Eq. 42 CLN
            IF (trialphiMC.gt.phiMC) THEN ! update if criteria at b(p) is max 
                phiMC = trialphiMC
                aFailC = b(p) ! record failure plane at max fail criteria
            END IF
            IF (trialphiMT.gt.phiMT) THEN
                phiMT = trialphiMT
                aFailT = b(p)
            END IF
        END DO
       RETURN
      END  
      
c Function to determine the misalignment angle phi using the Newton-
c Raphson method (uses functions FUNC and DFUNC as inputs)
      FUNCTION NEWT(init, trialStress, X)
      REAL*8 NEWT, xOld, xNew, tol, init, df,dx,f, trialStress(6), X
      xOld = init
      xNew = init
      dx = 100
      tol = 0.0001
10    IF (abs(dx).gt.tol) THEN
        xOld = xNew
        dx=FUNC(xOld, trialStress, X)/DFUNC(xOld, trialStress, X)
        xNew = xOld-dx
        GOTO 10
      ELSE
        NEWT = xNew
      END IF
      RETURN
      END

c Definition of function used in the Newton Raphson iteration method
      FUNCTION FUNC(gammaOld, trialStress, X)
      REAL*8 FUNC, X, gammaOld, trialStress(6)
      ! Eq. 88 CLN
      FUNC = X*gammaOld + 0.5d0*(trialStress(1) - trialStress(2))
     *      *sin(2.0d0*gammaOld) - abs(trialStress(4))*cos(2.0d0*gammaOld) 
      RETURN
      END
  
c Definiton of the derivative function of FUNC
      FUNCTION DFUNC(gammaOld, trialStress, X)
      REAL DFUNC, X, gammaOld, trialStress(6)
      ! Eq. 89 CLN
      DFUNC = X + (trialStress(1) - trialStress(2))*cos(2.0d0*gammaOld)
     *      + 2.0d0*abs(trialStress(4))*sin(2.0d0*gammaOld) 
      RETURN
      END
  
  
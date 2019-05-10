C-----------------------------------------------------------------------C
C     ABAQUS VUMAT USER SUBROUTINE: vumat_LARC04.for                    C
C																	    C
C     Version 1.0                                                       C
C     Release 1                                                         C
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
      real*8 s11,s22,s12,s11m,s22m,s12m
      real*8 phi1Plus,phi1Minus,phi2Plus,phi2Minus
      real*8 initialcharLength,thickness,traceStrain,meanDamage,expo, A
      character*80 cmname
      parameter ( zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.0d0,
     *     third = 1.d0 / 3.d0, half = 0.5d0, six = 6.0d0 )


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

      ST=YC*cos(alpha0)*(sin(alpha0)+cos(alpha0)/tan(2.0d0*alpha0))
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

      enddo

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

      s11=trialStress(1)
      s22=trialStress(2)
      s12=trialStress(4)

C Fibre direction

       phi1Plus=0.0d0
       phi1Minus=0.0d0

       if (s11.gt.0.0) then

       phi1Plus=trialStrain(1)*e11/XT      !Eq.8

       else 

       s22m=s22*Cos(phic)**2.0d0-2.0*Abs(s12)*Cos(phic)*Sin(phic)   !Eq.11
     ++s11*Sin(phic)**2

      s12m=(-s11+s22)*Cos(phic)*Sin(phic)            !Eq. 11
     ++Abs(s12)*(Cos(phic)**2.0d0-Sin(phic)**2.0d0)

      dummy=0.0
      dummysalida=0.0

       dummy=abs(s12m)+etaL*s22m
       call McAulay(dummy,dummysalida)

      phi1Minus=min(dummysalida/SL,(etaL-1.0d0)*s11/(2.0d0*SL))       !Eq.9,21

      endif

C Matrix direction

      phi2Plus=0.0d0
      phi2Minus=0.0d0


      if(s22.gt.0) then
        phi2Plus=sqrt((s22/YT)**2.0+(s12/SL)**2.0)      !Eq.13
      else
        dummy=abs(s12)+etaL*s22                         !Eq.13
        call McAulay(dummy,dummysalida)
        phi2Plus=dummysalida/SL
      endif

        if(s22.ne.0) then
        theta=atan(-abs(s12)/(s22*sin(alpha0)))       !Eq.16
        else
        theta=1.570796327
        endif
        dummy=-s22*cos(alpha0)*(sin(alpha0)-etaT*cos(alpha0)*cos(theta))   !Eq.15
        call McAulay(dummy,tauT)
        dummy=cos(alpha0)*(abs(s12)+etaL*s22*cos(alpha0)*sin(theta))    !Eq.15
        call McAulay(dummy,tauL)

       if(s22.lt.0) then
       phi2Minus=min(sqrt((tauT/ST)**2.0+(tauL/SL)**2.0),s22/omegaValue)  !Eq14,20
       endif

C Update of the damage thresholds

       rNew1=max(1.0,max(phi1Plus,rOld1),max(phi1Minus,rOld2))    !Eq.26
       rNew2=max(1.0,max(phi1Minus,rOld2))
       rNew3=max(1.0,max(phi2Plus,rOld3),max(phi2Minus,rOld4))
       rNew4=max(1.0,max(phi2Minus,rOld4))
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

       if(s11.gt.0) then      !Eq.6
          xOmega1=d1Plus
       else
          xOmega1=d1Minus
       endif

       if(s22.gt.0) then
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

! Determination of the softening parameter A  (exponential)      

      subroutine SofteningParam(lch,X,E,G,cte,trialStrain,k,A)     !It will only be used when phi=1.0d0
        real*8 E,X,cte,lch,G,trialStrain
        real*8 Strain(1000), dam(1000), stress(1000),phi(1000)
        real*8 Lim, interval, Error
        real*8 G2, G1, GArea
        real*8 A
        integer i, j, length_integral, n, k

           Lim = 100.0*trialStrain
           length_integral = 1000
           interval = (Lim-trialStrain)/length_integral
           G1 = 0.5*abs(trialStrain)*X*lch
           if (G1.gt.G)then
                   G=1.5*G1
             WRITE(*,*) 'Real energy, G in element', G, k
           endif
           A = 0.01+(2*lch*X*X/((2*E*G)-(lch*X*X)))
           
           do i = 1,length_integral
              Strain(i) = trialStrain +i*interval
              phi(i) = dsqrt((E*Strain(i)/X)**2.0d0+cte)
           enddo

           Error = 1.0d0 
           n = 0

666        if (Error.gt.0.01)then
              A = A - 0.001
              G2 = 0.0d0
              do j = 1,length_integral
                 dam(j) = 1.0d0-((1.0d0/phi(j))*dexp(A*(1.0d0-phi(j))))
                 stress(j) = Strain(j)*E*(1.0d0-dam(j))
                 G2 = G2 + stress(j)* interval
              enddo
              GArea = (G2 *lch) + G1
              Error = G-GArea
              Error = abs(Error)
              n = n +1
              if (n.gt.9999)then
                      A = (2*lch*X*X/((2*E*G)-(lch*X*X)))
                      goto 777
              endif
              if (A.lt.0.0d0)then
                      A = (2*lch*X*X/((2*E*G)-(lch*X*X)))
          WRITE(*,*) 'Negative A' 
                      goto 777
              endif
           else
                   goto 777
           endif   

           goto 666

777        continue           

        return

       end

C Damage varable, triangular damage dissipation       
       subroutine Damage(X,E,G,trialStrain,lch,dNew)
        real*8 X,E,G,trialStrain,lch,dNew
        real*8 delta_0, delta_c
           delta_0 = X/E
           delta_c = (2.0d0*G/X*lch)
           if (delta_c.lt.delta_0)then               !!Needed for d.nq.0-0
                   delta_c = 1.1*delta_0             !! 1.1 adjustable for stability
           endif
           dNew = (delta_c*(abs(trialStrain)-delta_0))
     *            /(abs(trialStrain)*(delta_c-delta_0))
        return
       end


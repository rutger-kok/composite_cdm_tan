      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE ’ABA_PARAM.INC’
C
      common/crdflg/lrdflg
C
      CHARACTER*80 CMNAME,ORNAME,CMNAME1
      CHARACTER*3 FLGRAY(15)
      CHARACTER xoutdir*255, xfname*80
      CHARACTER dmkname*255, FNAMEX*80
      DIMENSION UVAR(*),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C
      DIMENSION stress(6)
C The dimensions of the variables FLGRAY, ARRAY and JARRAY
C must be set equal to or greater than 15.
      double precision alphao,alphamem(1),psimem(1),thetamem(1),
     1      omega(1),lambda(1),
     2      fmat(1),fkink(1),fft(1),epsmato(1),
     3      sigmato(1),epskinko(1),sigkinko(1),epsfto(1),
     4      sig1(1),sig2(1),sig3(1),tau12(1),tau23(1),tau31(1),
     5      eps1(1),eps2(1),eps3(1),eps12(1),eps23(1),eps31(1),
     6      s12,s23,fio,beta
C
      integer lft,llt
C
      pi=dacos(-1.d0)
      degtorad=pi/180.d0
C
      do i=1,nuvarm
       uvar(i) = 0.d0
      enddo
C----------------------------------------------------------------------
C Open and read input file with material properties:
C directory/jobname.mt
C----------------------------------------------------------------------
      lxfname = 0
      lxoutdir = 0
      xfname =’ ’
      xoutdir =’ ’
C
      call getjobname(xfname,lxfname) ! input file name
      call getoutdir(xoutdir,lxoutdir) ! output directory
C
      if(lrdflg.ne.1) then
       fnamex=dmkname(xfname(1:lxfname),xoutdir(1:lxoutdir),’.mt’)
       open(unit=17,file=fnamex,status=’old’)
       lrdflg = 1
      endif
C
      read (17,*)
      read (17,*) klarc
C
      CMNAME1=’**dummy_name**’
      do while(CMNAME1.NE.CMNAME) ! search for material type
       read (17,*)
       read (17,*) CMNAME1
       if(CMNAME1.EQ.CMNAME) then
        read (17,*)
        read (17,*) ym1, ym2, ym3, nu21, nu31, nu32
        read (17,*)
        read (17,*) g12, g23, g31, xt, xc, yt, yc, s12
        read (17,*)
        read (17,*) alphao, beta, g, slis
       else
        do i=1,6
         read(17,*)
        enddo
       endif
      enddo
C
      rewind 17
C----------------------------------------------------------------------
C Compute derived material properties
C----------------------------------------------------------------------
      alphao=alphao*degtorad
      st=yc*dcos(alphao)*(dsin(alphao)+dcos(alphao)/dtan(2.d0*alphao))
      s23=st
      sl=s12
      etat=-1.d0/dtan(2.d0*alphao)
      etal=-s12*dcos(2*alphao)/(yc*dcos(alphao)*dcos(alphao))
C----------------------------------------------------------------------
C Read stress tensor from current increment
C----------------------------------------------------------------------
      CALL GETVRM(’S’,ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,
     1      LACCFLA)
C
      if(klarc.eq.3) then ! LaRC03 failure criteria
       stress(1) = array(1)
       stress(2) = array(2)
       stress(3) = array(3)
       stress(4) = array(4)
       stress(5) = array(5)
       stress(6) = array(6)
C
       call larc03(stress(1),stress(2),stress(3),stress(4),
     1      stress(5),stress(6),XT,XC,YT,YC,
     2      SL,SLIS,ST,G,G12,ETAL,ETAT,NDIM,UVAR,ANGLES,
     3      NOUT,NUVARM)
      endif
C----------------------------------------------------------------------
*
* End of main program
*
C----------------------------------------------------------------------
      RETURN
      END

* <<<<<<<<<<<<<<<<<<<<<<<< SUBROUTINE LARC03 >>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* LaRC03 failure criteria *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE LaRC03(S11,S22,S33,S12,S13,S23,XT,XC,YT,YC,
     1      SL,SL_IS,ST,G,G12,ETA_L,ETA_T,NDIM,FI,ANGLES,
     2      NOUT,NUVARM)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION S11,S22,S33,S12,S13,S23,XT,XC,YT,YC,
     1      SL,SL_IS,ST,G,G12,ETA_L,ETA_T,FI(*),ANGLES,
     2      PI,S11_M,S22_M,S33_M,S12_M,S13_M,S23_M,FLARC03,
     3      FMCCAULEY,ALPHA,FIP(7)
C
      INTEGER NDIM,NOUT,I,NUVARM
C
      PI = DACOS(-1.D0)
C
      do i=1,7
       fip(i) = 0.d0
      enddo
C-----------------------------------------------------------------------
*
* Transverse (matrix)
*
C----------------------------------------------------------------------
      IF(S22.GT.0.D0) THEN ! matrix tension
       IF(S11.LT.0.D0.AND.DABS(S11).LT.XC/2) THEN
        CALL ROTATE_PHI(SL,XC,ETA_L,S11,S22,S12,0.D0,0.D0,
     1      S11_M,S22_M,S12_M,S13_M,S23_M,G12,NDIM)
        FIP(1) = (1-G)*S22_M/YT+G*S22_M/YT*S22_M/YT+S12_M/SL_IS*
     1      S12_M/SL_IS
       ELSE
        FIP(2) = (1-G)*S22/YT+G*S22/YT*S22/YT+S12/SL_IS*S12/SL_IS
       ENDIF
C
      ELSE ! matrix compression
       IF(S11.GE.-YC) THEN
        FIP(3) = FLaRC03(ALPHA,S22,S12,ETA_L,ETA_T,SL_IS,ST,PI)
       ELSE
        CALL ROTATE_PHI(SL,XC,ETA_L,S11,S22,S12,0.D0,0.D0,
     1      S11_M,S22_M,S12_M,S13_M,S23_M,G12,NDIM)
       FIP(4) = FLaRC03(ALPHA,S22_M,S12_M,ETA_L,ETA_T,SL_IS,ST,PI)
       ENDIF
      ENDIF
C-----------------------------------------------------------------------
*
* Longitudinal (fibre)
* 
C----------------------------------------------------------------------
      IF(S11.GE.0.D0) THEN ! fibre tension
       FIP(5) = S11/XT
C
      ELSE ! fibre compression
       CALL ROTATE_PHI(SL,XC,ETA_L,S11,S22,S12,0.D0,0.D0,
     1      S11_M,S22_M,S12_M,S13_M,S23_M,G12,NDIM)
       IF(S22_M.LT.0.D0) THEN
        FIP(6) = FMcCAULEY((DABS(S12_M)+ETA_L*S22_M)/SL_IS) ! LaRC#4
       ELSEIF (DABS(S11).GE.XC/2.D0) THEN
        FIP(7) = (1-G)*S22_M/YT+G*S22_M/YT*S22_M/YT+S12_M/SL_IS
     1      *S12_M/SL_IS
       ENDIF
      ENDIF
C
      FI(1) = MAX(FIP(1),FIP(2)) ! Transverse with S22>0
C
      FI(2) = MAX(FIP(3),FIP(4)) ! Transverse with S22<0
C
      FI(3) = FIP(5) ! Longitudinal with S11>0
C
      FI(4) = MAX(FIP(6),FIP(7)) ! Longitudinal with S11<0
C----------------------------------------------------------------------
      RETURN
      END
* <<<<<<<<<<<<<<<<<<<<<<<<< FUNCTION FLaRC03 >>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* MATRIX COMPRESSION FAILURE CRITERION (LaRC03) *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      REAL*8 FUNCTION FLaRC03(ALPHA,S22,S12,ETAL,ETAT,SL_IS,ST,PI)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION ALPHA,S22,S12,ETAL,ETAT,SL_IS,ST,PI,THETA,
     1      TAUT_EFF, FMCCAULEY,TAUL_EFF,FIC,ALPHA1,taul,taut
      INTEGER I
C Cycle over possible fracture angles
      FLaRC03=0.d0
      DO i=0,56 ! Determination of the fracture angle
       ALPHA1 = i*PI/180.D0
       IF(ALPHA1.EQ.0.D0.OR.S22.EQ.0.D0) THEN ! Avoids divisions by zero
        THETA = PI/2.D0
       ELSE
        THETA = DATAN(-dabs(S12)/(S22*DSIN(ALPHA1)))
       ENDIF
c
       TAUT_EFF = FMcCAULEY(-S22*dcos(alpha1)*(dsin(alpha1)-etat*
     1       dcos(alpha1)*dcos(theta)))
c
       TAUL_EFF = FMcCAULEY(dcos(alpha1)*(dabs(s12)+etal*s22*
     1      dcos(alpha1)*dsin(theta)))
c
       FIC = (TAUT_EFF/ST)*(TAUT_EFF/ST)+
     1       (TAUL_EFF/SL_IS)*(TAUL_EFF/SL_IS)
C
       FLaRC03 = max(FLaRC03,FIC)
      ENDDO
C
      RETURN
      END
* <<<<<<<<<<<<<<<<<<<<<< SUBROUTINE ROTATE_PHI>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* ROTATION OF STRESSES TO THE MISALIGNMENT COORDINATE FRAME *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      SUBROUTINE ROTATE_PHI(SL,XC,ETAL,S11,S22,S12,S13,S23,
     1      S11T,S22T,S12T,S13T,S23T,G12,NDI)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION SL,XC,ETAL,S11,S22,S12,S13,S23,
     1      S11T,S22T,S12T,S13T,S23T,G12,aa,cc,phiC,
     2      PHI,sqr,phi0,cp,ss,c2,s2
C
      INTEGER NDI
C Calculate fiber misalignment angle (Linear shear law and
C small angle approximations)
      cc = dABS(SL/XC)
      aa = cc+ETAL
      sqr = dsqrt(1.d0-4.0d0*aa*cc)
      phiC = datan((1.d0-sqr)/(2.0d0*aa)) ! select smallest root
C
      phi0 = (dabs(S12)+(G12-XC)*phiC)/(G12+S11-S22)
c
      cp = dcos(phi0)
      ss = dsin(phi0)
      c2 = cp*cp
      s2 = ss*ss
C Calculate stresses in misalignment coordinate frame C
      S11T = S11*c2+S22*s2+2.0d0*cp*ss*DABS(S12)
      S22T = S11*s2+S22*c2-2.0d0*cp*ss*DABS(S12)
      S12T = -ss*cp*S11+ss*cp*S22+(c2-s2)*DABS(S12)
C
      RETURN
      END
* <<<<<<<<<<<<<<<<<<<<<<< FUNCTION FMcCAULEY >>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* McCAULEY OPERATOR *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      REAL*8 FUNCTION FMcCAULEY(X)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION X
C
      IF(X.LE.0.D0) THEN
       FMcCAULEY = 0.D0
      ELSE
       FMcCAULEY = X
      ENDIF
C
      RETURN
      END
* <<<<<<<<<<<<<<<<<<<<<<<< FUNCTION DMKNAME >>>>>>>>>>>>>>>>>>>>>>>>> *
* *
* Compose a filename directory/jobname.exten *
* *
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *
      character*(*) function dmkname(fname,dname,exten)
C
      character*(*) fname,dname,exten
C C     fname I jobname C dname I directory C exten I
Cextension C dmkname O directory/jobname.exten C
      ltot = len(fname)
      lf = 0
      do k1 = ltot,2,-1
       if (lf.eq.0.and.fname(k1:k1).ne.’ ’) lf = k1
      end do
C
      ltot = len(dname)
      ld = 0
      do k1 = ltot,2,-1
       if (ld.eq.0.and.dname(k1:k1).ne.’ ’) ld = k1
      end do
C
      ltot = len(exten)
      le = 0
      do k1 = ltot,2,-1
       if (le.eq.0.and.exten(k1:k1).ne.’ ’) le = k1
      end do
C
      if ((lf + ld + le) .le. len(dmkname)) then
       dmkname = dname(1:ld)//’/’//fname(1:lf)
       ltot = ld + lf + 1
       if ( le.gt.0) then
        dmkname = dmkname(1:ltot)//exten(1:le)
       end if
      end if
C
      return
      end
C=======================================================================C
C ==== end of program ====
C=======================================================================C

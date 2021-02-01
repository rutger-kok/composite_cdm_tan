!----------------------------------------------------------------------!
!     ABAQUS VUMAT USER SUBROUTINE: catalanotti_FC_3D.for              !
!     Author(s): Rutger Kok, Francisca Martinez-Hergueta               !
!     Last Updated: 01/02/2021                                         !
!     Version 1.0                                                      !
!----------------------------------------------------------------------!

! This is a VUMAT subroutine implementing a continuum damage mechanics
! framework for composite materials in Abaqus. The subroutine is based
! on the work of Catalanotti et al. (failure criteria) [1], Maimi
! et al. [2][3], and Camanho et al. [4]. 

! [1] G. Catalanotti, P.P. Camanho, A.T. Marques
! Three-dimensional failure criteria for fiber-reinforced laminates
! Composite Structures 95 (2013) 63–79
! http://dx.doi.org/10.1016/j.compstruct.2012.07.016

! [2] P. Maimi, P.P. Camanho, J.A. Mayugo, C.G. Davila
! A continuum damage model for composite laminates: Part I –
! Constitutive model
! Mechanics of Materials 39 (2007) 897–908
! http://dx.doi.org/10.1016/j.mechmat.2007.03.005

! [3] P. Maimi, P.P. Camanho, J.A. Mayugo, C.G. Davila
! A continuum damage model for composite laminates: Part II –
! Computational implementation and validation
! Mechanics of Materials 39 (2007) 909–919
! http://dx.doi.org/10.1016/j.mechmat.2007.03.006

! [4] P.P. Camanho, M.A. Bessa, G. Catalanotti, M. Vogler, R. Rolfes
! Modeling the inelastic deformation and fracture of polymer
! composites – Part II: Smeared crack model
! Mechanics of Materials 59 (2013) 36–49
! http://dx.doi.org/10.1016/j.mechmat.2012.12.001

! ----------------------------------------------------------------------
! VARIABLE DICTIONARY (main program)
! For a full rundown of each Abaqus variable see the Abaqus VUMAT
! documentation. 
! nblock = Number of material points to be processed in this call to
!          VUMAT. In other words, the number of integration points,
!          which depends on the element choice. C3D8R elements have 1
!          integration point, so nblock=number of elements.
! props = material property array
! charLength = characteristic element length
! stepTime = time since step began
! nstatev = number of user-defined state variables that are associated
!           with this material type, 24 in this case (not all state
!           variable slots are used)
! strainInc = strain increment tensor at each material point
! stressOld = stress tensor at each material point at the beginning of
!             the increment.
! stateOld = state variables at each material point at the beginning
!            of the increment.
! stressNew = stress tensor at each material point at the end of the
!             increment.
! stateNew = state variables at each material point at the end of the
!            increment
! enerInternNew = internal energy per unit mass at each material point
!                 at the end of the increment.
! C = stiffness matrix (6x6)
! k = material point index
! i = integer used for indexing arrays in a number of do loops
! strain = strain tensor in Voigt notation, see Abaqus documentation
!           for more info (strain is diff. in implicit and explicit)
! stress = stress tensor in Voigt notation
! stressPower = internal energy per unit mass
!----------------------------------------------------------------------!

! User subroutine VUMAT

      SUBROUTINE vumat (
      ! Read only - input variables from Abaqus
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
      ! Write only - outputs Abaqus needs from subroutine
     7     stressNew, stateNew, enerInternNew, enerInelasNew )

      INCLUDE 'vaba_param.inc'
      
      ! Dimension the Abaqus input variables
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
      
      ! Declare variables used in main part of script
      REAL*8, DIMENSION(6) :: stress,strain
      REAL*8, DIMENSION(6,6) :: C
      REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23 ! elastic constants
      REAL*8 nu21,nu31,nu32,delta
      REAL*8 G1plus,G1minus,G2plus,G2minus,G6 ! fracture energies
      REAL*8 XT,XC,YT,YC,SL,alpha0,stressPower
      INTEGER k,i

      ! Elastic constants orthotropic ply
      e11 = props(1) ! stiffness fiber direction
      e22 = props(2) ! stiffness transverse direction (in-plane)
      e33 = props(3) ! stiffness transverse direction (out-of-plane)
      nu12 = props(4) ! Poisson's ratio 12 direction
      nu13 = props(5) ! Poisson's ratio 13 direction
      nu23 = props(6) ! Poisson's ratio 23 direction
      g12 = props(7) ! shear modulus 12 direction
      g13 = props(8) ! shear modulus 13 direction
      g23 = props(9) ! shear modulus 23 direction

      ! Ply strengths
      XT = props(10) ! tensile strength fiber direction
      XC = props(11) ! compressive strength fiber direction
      YT = props(12) ! tensile strength transverse direction
      YC = props(13) ! compressive strength transverse direction
      SL = props(14) ! shear strength

      ! Fracture toughnesses
      G1plus = props(16) ! tensile frac. toughness fiber direction
      G1minus = props(17) ! compressive frac. toughness fiber direction
      G2plus = props(18) ! tensile frac. toughness transverse direction
      G2minus = props(19) ! compressive frac. toughness transverse direction
      G6 = props(20) ! shear frac. toughness

      ! Calculate remaining Poisson's ratios
      nu21 = nu12*(e22/e11) ! Poisson's ratio 21 direction
      nu31 = nu13*(e33/e11) ! Poisson's ratio 31 direction
      nu32 = nu23*(e33/e22) ! Poisson's ratio 32 direction

      ! Calculate undamaged stiffness matrix
      delta = 1.0d0/(e22**2.0d0*nu12**2.0d0 +
     1        2.0d0*e33*e22*nu12*nu13*nu23 +
     1        e33*e22*nu13**2.0d0 - e11*e22 + e11*e33*nu23**2.0d0)

      C = 0.0d0 ! initialize stiffness matrix, set all elements = 0
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
 
      ! Initial elastic step, for Abaqus tests 
      IF (stepTime.eq.0) THEN        
        DO k = 1, nblock
          ! Initialize state variables
          DO i = 1,nstatev
            stateNew(k,i) = 0.d0
          ENDDO
          ! Initial strain
          DO i = 1,6
            strain(i) = strainInc(k,i)
          ENDDO
          ! Calculate initial stress (multiply stiffness and strain)
          stress = matmul(C(1:6,1:6), strain(1:6))
          ! Assign calculated stresses to stressNew array
          DO i = 1,6
            stressNew(k,i) = stress(i)
          ENDDO
        ENDDO
      ELSE
        ! If not initial step, calc. stresses according to CDM model
        DO k = 1,nblock
          ! Calc. new strain and assign to state variable array (1-6)
          DO i = 1,6
            stateNew(k,i) = stateOld(k,i) + strainInc(k,i)
          ENDDO
          ! Assign strain values to strain array
          DO i = 1,6
            strain(i) = stateNew(k,i)
          ENDDO
          ! Call CDM subroutine
          call cdm(k,nblock,nstatev,strain,stateOld,charLength(k),
     1             C,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2             XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,G2minus,G6,
     3             alpha0,stress,stateNew)
          ! Assign stresses output by cdm subroutine to stressNew array
          DO i = 1,6
            stressNew(k,i) = stress(i)
          ENDDO
          ! Calculate internal energy
          stressPower = 0.5d0*(
     1              (stressNew(k,1)+stressOld(k,1))*strainInc(k,1)+
     2              (stressNew(k,2)+stressOld(k,2))*strainInc(k,2)+
     3              (stressNew(k,3)+stressOld(k,3))*strainInc(k,3)+
     4              2.0*(stressNew(k,4)+stressOld(k,4))*strainInc(k,4)+
     5              2.0*(stressNew(k,5)+stressOld(k,5))*strainInc(k,5)+
     6              2.0*(stressNew(k,6)+stressOld(k,6))*strainInc(k,6))
          ! Assign energies to enerInternNew array
          enerInternNew(k) = enerInternOld(k)+stressPower/density(k)
        ENDDO
      ENDIF
      END SUBROUTINE vumat

!----------------------------------------------------------------------!
! Subroutine cdm: Failure criteria and damage evolution.               !
! Returns a stress matrix                                              !
!----------------------------------------------------------------------!
! VARIABLE DICTIONARY (cdm subroutine)
! l_ch = characteristic element length, note: for a single element!
! trialStress = stress calculated using the undamaged stiffness tensor.
! d1Plus,d1Minus,d2Plus,d2Minus,d3Plus,d3Minus = scalar damage
!     variables quantifying damage in the fiber (1) and transverse
!     (2,3) directions. Minus and Plus refer to tensile and
!     compressive loading respectively.
! d1,d2,d3,d6 = scalar damage variables quantifying damage in the
!     fiber (1), transverse (2,3) and in-plane shear (6) directions.
! dState1,dState2,dState3,dState6 = maximum scalar damage variables
!     for the 1,2,3,6 directions. These ensure damage state retention
!     as loads fluctuate.
! dStateOld1, dStateOld2, dStateOld3, dStateOld6 = maximum scalar
!     damage variables from the previous iteration.
! FI_LT = longitudinal (fiber direction) tensile failure index
! FI_LC = longitudinal (fiber direction) compressive failure index
! FI_MT = longitudinal (transvere direction) tensile failure index
! FI_MC = longitudinal (transverse direction) compressive failure index
! minDamage = element deletion flag. If element is damaged in all
!     directions the element is flagged for deletion.
!----------------------------------------------------------------------!

      SUBROUTINE cdm(k,nblock,nstatev,strain,stateOld,l_ch,C,
     1               e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2               XT,XC,YT,YC,SL,G1plus,G1minus,G2plus,G2minus,G6,
     3               alpha0,stress,stateNew)
        IMPLICIT NONE
        ! input variables
        INTEGER, INTENT(IN) :: nblock,nstatev,k
        REAL*8, DIMENSION(nblock,nstatev), INTENT(IN) :: stateOld
        REAL*8, DIMENSION(6,6), INTENT(INOUT) :: C
        REAL*8, DIMENSION(6), INTENT(IN) :: strain
        REAL*8, INTENT(IN) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
        REAL*8, INTENT(IN) :: XT,XC,YT,YC,SL,l_ch
        REAL*8, INTENT(IN) :: G1plus,G1minus,G2plus,G2minus,G6
        ! local variables
        REAL*8, DIMENSION(6) :: trialStress, trialStressP
        REAL*8 d1Plus,d1Minus,d2Plus,d2Minus,d3Plus,d3Minus
        REAL*8 d1,d2,d3,d6
        REAL*8 nu21,nu31,nu32,delta,minDamage
        REAL*8 FI_LT,FI_LC,FI_MT,FI_MC
        REAL*8 dState1,dState2,dState3,dState6
        REAL*8 dState1Old,dState2Old,dState3Old,dState6Old
        ! output variables
        REAL*8, DIMENSION(6), INTENT(OUT) :: stress
        REAL*8, DIMENSION(nblock,nstatev), INTENT(INOUT) :: stateNew
        ! external functions
        REAL*8, EXTERNAL :: triangular, mccauley

        ! Load old dState values
        dState1Old = stateOld(k,13)
        dState2Old = stateOld(k,14)
        dState3Old = stateOld(k,15)
        dState6Old = stateOld(k,16)

        ! Calculate stress using an undamaged stiffness matrix
        trialStress = matmul(C(1:6,1:6), strain(1:6))

        ! Initialize damage activation functions (failure indices)
        FI_LT = 0.0d0
        FI_LC = 0.0d0
        FI_MT = 0.0d0
        FI_MC = 0.0d0

        ! Calculate failure indices
        IF (trialStress(1).gt.0.d0) THEN
            ! If stress in fiber direction is positive calculate the
            ! longitudinal failure index in tension
            FI_LT = trialStress(1)/XT ! max stress criterion
        ELSEIF (trialStress(1).lt.0.d0) THEN
            ! Else if stress in fiber direction is negative calculate the
            ! longitudinal failure index in compression
            FI_LC = trialStress(1)/XT
        ENDIF

        IF (trialStress(2).gt.0.d0) THEN
            ! If stress in transverse direction is positive calculate the
            ! transverse failure index in tension
            FI_MT = trialStress(1)/YT ! max stress criterion
        ELSEIF (trialStress(2).lt.0.d0) THEN
            ! Else if stress in transverse direction is negative calculate the
            ! transverse failure index in compression
            FI_MC = trialStress(1)/YC
        ENDIF

        ! Calculate scalar damage variables assuming triangular
        ! damage dissipation.
        IF (FI_LT.gt.1.0d0) THEN ! Longitudinal tensile damage variable
            d1Plus= triangular(XT,e11,G1Plus,strain(1),l_ch)
        ELSE
            d1Plus = 0.0d0
        ENDIF

        IF (FI_LC.gt.1.0d0) THEN ! Longitudinal compressive damage variable
            d1Minus = triangular(XC,e11,G1minus,strain(1),l_ch)
        ELSE
            d1Minus = 0.0d0
        ENDIF

        IF (FI_MT.gt.1.0d0) THEN ! Transverse (2) tensile damage variable
            d2Plus = triangular(YT,e22,G2Plus,strain(2),l_ch)
        ELSE
            d2Plus = 0.0d0
        ENDIF

        IF (FI_MC.gt.1.0d0) THEN ! Transverse (2) compressive damage variable
            d2Minus = triangular(YC,e22,G2minus,strain(2),l_ch)
        ELSE
            d2Minus = 0.0d0
        ENDIF

        IF (FI_MT.gt.1.0d0) THEN ! Transverse (3) compressive damage variable
            d3Plus = triangular(YT,e33,G2Plus,strain(3),l_ch)
        ELSE
            d3Plus = 0.0d0
        ENDIF

        IF (FI_MC.gt.1.0d0) THEN ! Transverse (3) compressive damage variable
            d3Minus = triangular(YC,e33,G2minus,strain(3),l_ch)
        ELSE
            d3Minus = 0.0d0
        ENDIF

        IF ((FI_MT.gt.1.0d0).or.(FI_MC.gt.1.0d0)) THEN ! Shear damage variable
            d6 = triangular(SL,e22,G6,strain(4),l_ch)
        ELSE
            d6 = 0.0d0
        ENDIF

        ! Calculate damage variables
        IF (trialStress(1).ne.0.0d0) THEN ! Longitudinal damage variable
          d1 = d1Plus*(mccauley(trialStress(1))/abs(trialStress(1))) + 
     1         d1Minus*(mccauley(-trialStress(1))/abs(trialStress(1)))
        ELSE
          d1 = 0.0d0
        END IF
        
        IF (trialStress(2).ne.0.0d0) THEN ! Longitudinal damage variable
          d2 = d2Plus*(mccauley(trialStress(2))/abs(trialStress(2))) + 
     1         d2Minus*(mccauley(-trialStress(2))/abs(trialStress(2)))
        ELSE
          d2 = 0.0d0
        END IF 
        
        IF (trialStress(3).ne.0.0d0) THEN ! Longitudinal damage variable
          d3 = d3Plus*(mccauley(trialStress(3))/abs(trialStress(3))) + 
     1         d3Minus*(mccauley(-trialStress(3))/abs(trialStress(3)))
        ELSE 
          d3 = 0.0d0 
        END IF

        ! Calculate damage state variables. Damage state must be
        ! retained between iterations.
        dState1 = min(0.999,max(d1,dState1Old))
        dState2 = min(0.999,max(d2,dState2Old))
        dState3 = min(0.999,max(d3,dState3Old))
        dState6 = min(0.999,max(d6,dState6Old))

        ! Calculate damaged stiffness tensor
        ! Recalculate additional Poisson's ratio terms (not passed as args.)
        nu21 = nu12*(e22/e11)
        nu31 = nu13*(e33/e11)
        nu32 = nu23*(e33/e22)
 
        delta = 1.0d0/(e11*e22 - e22**2.0d0*nu12**2.0d0 + 
     1          e22**2.0d0*dState1*nu12**2.0d0 +
     2          e22**2.0d0*dState2*nu12**2.0d0 - 
     3          e22*e33*nu13**2.0d0 - e11*e33*nu23**2.0d0 +
     4          e22*e33*dState1*nu13**2.0d0 + 
     5          e22*e33*dState3*nu13**2.0d0 + 
     6          e11*e33*dState2*nu23**2.0d0 +
     7          e11*e33*dState3*nu23**2.0d0 -
     8          e22**2.0d0*dState1*dState2*nu12**2.0d0 -
     9          2.0d0*e22*e33*nu12*nu13*nu23 -
     1          e22*e33*dState1*dState3*nu13**2.0d0 -
     2          e11*e33*dState2*dState3*nu23**2.0d0 +
     3          2.0d0*e22*e33*dState1*nu12*nu13*nu23 +
     4          2.0d0*e22*e33*dState2*nu12*nu13*nu23 +
     5          2.0d0*e22*e33*dState3*nu12*nu13*nu23 -
     6          2.0d0*e22*e33*dState1*dState2*nu12*nu13*nu23 -
     7          2.0d0*e22*e33*dState1*dState3*nu12*nu13*nu23 -
     8          2.0d0*e22*e33*dState2*dState3*nu12*nu13*nu23 +
     9          2.0d0*e22*e33*dState1*dState2*dState3*nu12*nu13*nu23)

        C = 0.0d0 ! set all elements in array equal to zero
        C(1,1) = -(e11**2.0d0*(dState1 - 1.0d0)*(e22 - e33*nu23**2.0d0 +
     1           e33*dState2*nu23**2.0d0 + e33*dState3*nu23**2.0d0 -
     2           e33*dState2*dState3*nu23**2.0d0))*delta
        C(2,1) = (e11*e22*(dState1 - 1.0d0)*(dState2 - 1.0d0)*
     1           (e22*nu12 + e33*nu13*nu23 - e33*dState3*nu13*nu23))
     2           *delta
        C(3,1) = (e11*e22*e33*(dState1 - 1.0d0)*(dState3 - 1.0d0)*
     1           (nu13 + nu12*nu23 - dState2*nu12*nu23))*delta
        C(1,2) = (e11*e22*(dState1 - 1.0d0)*(dState2 - 1.0d0)*
     1           (e22*nu12 + e33*nu13*nu23 - e33*dState3*nu13*nu23))
     2           *delta
        C(2,2) = -(e22**2.0d0*(dState2 - 1.0d0)*(e11 - e33*nu13**2.0d0 +
     1           e33*dState1*nu13**2.0d0 + e33*dState3*nu13**2.0d0 -
     2           e33*dState1*dState3*nu13**2.0d0))*delta
        C(3,2) = (e22*e33*(dState2 - 1.0d0)*(dState3 - 1.0d0)*
     1           (e11*nu23 + e22*nu12*nu13 - e22*dState1*nu12*nu13))
     2           *delta
        C(1,3) = (e11*e22*e33*(dState1 - 1.0d0)*(dState3 - 1.0d0)*
     1           (nu13 + nu12*nu23 - dState2*nu12*nu23))*delta
        C(2,3) = (e22*e33*(dState2 - 1.0d0)*(dState3 - 1.0d0)*
     1           (e11*nu23 + e22*nu12*nu13 - e22*dState1*nu12*nu13))
     2           *delta
        C(3,3) = -(e22*e33*(dState3 - 1.0d0)*(e11 - e22*nu12**2.0d0 +
     1           e22*dState1*nu12**2.0d0 + e22*dState2*nu12**2.0d0 - 
     2           e22*dState1*dState2*nu12**2.0d0))*delta
        C(4,4) = -(g12*(dState6 - 1.0d0))
        C(5,5) = g23
        C(6,6) = -(g13*(dState6 - 1.0d0))
        
        ! Calculate new stress value using damaged stiffness tensor
        stress = matmul(C(1:6,1:6), strain(1:6))

        ! Record scalar damage variables state array
        stateNew(k,13) = dState1
        stateNew(k,14) = dState2
        stateNew(k,15) = dState3
        stateNew(k,16) = dState6
        ! Record failure indices in state array
        stateNew(k,21) = FI_LT
        stateNew(k,22) = FI_LC
        stateNew(k,23) = FI_MT
        stateNew(k,24) = FI_MC

        ! Check damage state of each element. If element is damaged
        ! completely in each direction it is flagged for deletion
        minDamage = min(dState1,dState2,dState3,dState6)
        stateNew(k,19) = minDamage
        IF (minDamage.gt.0.99) stateNew(k,20) = 0.d0

      RETURN
      END SUBROUTINE cdm

!----------------------------------------------------------------------!
! Function triangular: implements triangular damage dissipation        !
!                                                                      !
!----------------------------------------------------------------------!
! VARIABLE DICTIONARY (rotate_stress function)
! X = failure stress (loading and material direction unspecified)
! E = modulus
! G = energy release rate (NOT shear modulus)
! tStrain = scalar trial strain value
! l_ch = characteristic element length
! delta_0 = strain at yield
! delta_c = strain at final failure
! dam = scalar damage variable
!----------------------------------------------------------------------!
      REAL*8 FUNCTION triangular(X,E,G,tStrain,l_ch)
        IMPLICIT NONE
          ! input variables
          REAL*8, INTENT(IN) :: X,E,G,tStrain,l_ch
          ! local variables
          REAL*8 delta_0, delta_c,dam
          ! Calculate scalar damage variable assuming triangular dissipation
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
C VARIABLE DICTIONARY (mccauley function)
C value = variable to execute mccauley brackets operation on
C----------------------------------------------------------------------C
      REAL*8 FUNCTION mccauley(value)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: value 
        mccauley = (value+abs(value))/2.0d0
        RETURN
      END FUNCTION mccauley



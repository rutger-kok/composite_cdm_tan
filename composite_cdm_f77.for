! This is a VUMAT subroutine implementing a continuum damage mechanics
! framework for composite materials in Abaqus.
! Copyright (C) 2021 Rutger Wouter Kok

! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
! 02110-1301 USA


      subroutine vumat(
        ! Read only - input variables from Abaqus
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
        ! Write only - outputs Abaqus needs from subroutine
     7  stressNew, stateNew, enerInternNew, enerInelasNew)

        include 'vaba_param.inc'
        
        ! Dimension the Abaqus input variables
        dimension props(nprops), density(nblock), coordMp(nblock),
     1    charLength(nblock), strainInc(nblock, ndir+nshr),
     2    relSpinInc(nblock, nshr), tempOld(nblock),
     3    stretchOld(nblock, ndir+nshr), fieldOld(nblock, nfieldv),
     4    defgradOld(nblock,ndir+nshr+nshr), stateOld(nblock, nstatev),
     5    stressOld(nblock, ndir+nshr), enerInternOld(nblock),
     6    enerInelasOld(nblock), tempNew(nblock),
     7    stretchNew(nblock, ndir+nshr), fieldNew(nblock, nfieldv),
     8    defgradNew(nblock,ndir+nshr+nshr), enerInelasNew(nblock),
     9    stressNew(nblock,ndir+nshr), stateNew(nblock, nstatev),
     1    enerInternNew(nblock)     
        character*80 cmname

        ! Purpose: This is a VUMAT subroutine implementing a CDM 
        ! framework for composite materials in Abaqus. The subroutine
        ! is based on the work of Catalanotti et al. [1] (for the
        ! failure criteria), Maimi et al. [2][3], and Camanho et al.
        ! [4] (damage evolution). 

        ! [1] G. Catalanotti, P.P. Camanho, A.T. Marques
        ! Three-dimensional failure criteria for fiber-reinforced
        ! laminates
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

        ! [4] P.P. Camanho, M.A. Bessa, G. Catalanotti, M. Vogler,
        ! R. Rolfes
        ! Modeling the inelastic deformation and fracture of polymer
        ! composites – Part II: Smeared crack model
        ! Mechanics of Materials 59 (2013) 36–49
        ! http://dx.doi.org/10.1016/j.mechmat.2012.12.001

        ! Variable dictionary:
        ! For a full rundown of each Abaqus variable see the Abaqus
        ! VUMAT documentation. 
        ! nblock = Number of material points to be processed in this
        !   call to VUMAT. In other words, the number of integration
        !   points, which depends on the element choice. C3D8R elements
        !   have 1 integration point, so nblock=number of elements.
        ! props = material property array
        ! charLength = characteristic element length array
        ! stepTime = time since step began
        ! nstatev = number of user-defined state variables that are
        !   associated with this material type, 24 in this case
        !   (not all state variable slots are used)
        ! strainInc = strain increment tensor at each material point
        ! stressOld = stress tensor at each material point at the
        !   beginning of the increment
        ! stateOld = state variables at each material point at the
        !   beginning of the increment
        ! stressNew = stress tensor at each material point at the end
        !   of the increment.
        ! stateNew = state variables at each material point at the end
        !   of the increment
        ! enerInternNew = internal energy per unit mass at each
        !   material point at the end of the increment.
        ! C = stiffness matrix (6x6)
        ! k = material point index
        ! i = integer used for indexing arrays in a number of do loops
        ! strain = strain tensor in Voigt notation, see Abaqus
        !   documentation for more info (strain is different in
        !   implicit and explicit)
        ! stress = stress tensor in Voigt notation
        ! stressPower = internal energy per unit mass

        ! Declare variables used in main part of script
        real*8, dimension(6) :: stress,strain
        real*8, dimension(6,6) :: C
        real*8 :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
        real*8 :: nu21,nu31,nu32,delta
        real*8 :: G1plus,G1minus,G2plus,G2minus,G6
        real*8 :: XT,XC,YT,YC,SL,alpha0,stressPower
        integer :: k,i,debugLoop

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

        ! Fracture angle
        alpha0 = props(15)*0.017453292519943295d0 ! converts to radians

        ! Fracture toughnesses
        G1plus = props(16) ! tensile frac. toughness fiber direction
        G1minus = props(17) ! comp. frac. toughness fiber direction
        G2plus = props(18) ! tensile frac. toughness trans. direction
        G2minus = props(19) ! comp. frac. toughness trans. direction
        G6 = props(20) ! shear frac. toughness

        ! Calculate remaining Poisson's ratios
        nu21 = nu12*(e22/e11) ! Poisson's ratio 21 direction
        nu31 = nu13*(e33/e11) ! Poisson's ratio 31 direction
        nu32 = nu23*(e33/e22) ! Poisson's ratio 32 direction

        ! Calculate undamaged stiffness matrix
        delta = 1.0d0/(e22**2.0d0*nu12**2.0d0 + 2.0d0*e33*e22*
     1          nu12*nu13*nu23 + e33*e22*nu13**2.0d0 - e11*e22 +
     2          e11*e33*nu23**2.0d0)

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
        if (stepTime == 0) then        
          do k = 1, nblock
            ! Initialize state variables
            do i = 1,nstatev
              stateNew(k,i) = 0.d0
            end do
            ! Initial strain
            do i = 1,6
              strain(i) = strainInc(k,i)
            end do
            ! Calculate initial stress (multiply stiffness and strain)
            stress = matmul(C(1:6,1:6), strain(1:6))
            ! Assign calculated stresses to stressNew array
            do i = 1,6
              stressNew(k,i) = stress(i)
            end do
          end do
        else
          ! infinite loop for debugging with Visual Studio
          !do while (debugLoop /= 999)
            !debugLoop = 1
          !end do
          ! If not initial step, calc. stresses according to CDM model
          do k = 1,nblock
            ! Calc. new strain and assign to state variable array (1-6)
            do i = 1,6
              stateNew(k,i) = stateOld(k,i) + strainInc(k,i)
            end do
            ! Assign strain values to strain array
            do i = 1,6
              strain(i) = stateNew(k,i)
            end do
            ! Call CDM subroutine
            call cdm(k,nblock,nstatev,strain,stateOld,charLength(k),C,
     1                e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,XT,XC,YT,
     2                YC,SL, G1plus,G1minus,G2plus,G2minus,G6,alpha0,
     3                stress,stateNew)
            ! Assign stresses from cdm subroutine to stressNew array
            do i = 1,6
              stressNew(k,i) = stress(i)
            end do
            ! Calculate internal energy
            stressPower = 0.5d0*((stressNew(k,1)+stressOld(k,1))*
     1                    strainInc(k,1)+(stressNew(k,2)+
     2                    stressOld(k,2))*strainInc(k,2)+
     3                    (stressNew(k,3)+stressOld(k,3))*
     4                    strainInc(k,3)+2.0d0*(stressNew(k,4)+
     5                    stressOld(k,4))*strainInc(k,4)+
     6                    2.0d0*(stressNew(k,5)+stressOld(k,5))*
     7                    strainInc(k,5)+2.0d0*(stressNew(k,6)+
     8                    stressOld(k,6))*strainInc(k,6))
            ! Assign energies to enerInternNew array
            enerInternNew(k) = enerInternOld(k)+stressPower/density(k)
          end do
        end if
      end subroutine vumat

      subroutine cdm(k,nblock,nstatev,strain,stateOld,lch,C,e11,e22,
     1               e33,nu12,nu13,nu23,g12,g13,g23,XT,XC,YT,YC,SL,
     2               G1plus,G1minus,G2plus,G2minus,G6,alpha0,stress,
     3               stateNew)
        ! Purpose: Implements continuum damage mechanics framework
        !   from Maimi et al. (2007)
        ! Variable Dictionary:
        ! lch = characteristic element length
        ! trialStress = stress calculated using the undamaged
        !   stiffness tensor.
        ! trialStressP = trialStress rotated into the misalignment
        !   coordinate frame to calculate the failure index for
        !   longitudinal compression.
        ! d1Plus,d1Minus,d2Plus,d2Minus,d3Plus,d3Minus = scalar
        !   damage variables quantifying damage in the fiber (1) and
        !   transverse (2,3) directions. Minus and Plus refer to
        !   tensile and compressive loading respectively.
        ! d1,d2,d3,d6 = scalar damage variables quantifying damage in
        !   the fiber (1), transverse (2,3) and in-plane shear (6)
        !   directions.
        ! dState1,dState2,dState3,dState6 = maximum scalar damage
        !   variables for the 1,2,3,6 directions. These ensure damage
        !   state retention as loads fluctuate.
        ! dStateOld1, dStateOld2, dStateOld3, dStateOld6 = maximum
        !   scalar damage variables from the previous iteration.
        ! FI_LT = fiber direction tensile failure index
        ! FI_LC = fiber direction compressive failure index
        ! FI_MT = transvere direction tensile failure index
        ! FI_MC = transverse direction compressive failure index
        ! etaT = friction coefficient in the transverse direction
        ! etaL = friction coefficient in the longitudinal direction
        ! phiC = initial misalignment angle for compressive failure
        !   calculation
        ! kappa = parameter used to calculate failure indices
        ! lambda = parameter used to calculate failure indices
        ! ST = in-situ (ideally) transverse shear strength
        ! minDamage = element deletion flag. If element is damaged
        !   in all directions the element is flagged for deletion.

        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev,k
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld
        real*8, dimension(6,6), intent(inout) :: C
        real*8, dimension(6), intent(in) :: strain
        real*8, intent(in) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
        real*8, intent(in) :: XT,XC,YT,YC,SL,lch,alpha0
        real*8, intent(in) :: G1plus,G1minus,G2plus,G2minus,G6
        ! local variables
        real*8, dimension(6) :: trialStress, trialStressP, e0
        real*8, dimension(3) :: tbar_cr
        real*8 d1Plus,d1Minus,d2Plus,d2Minus,d3Plus,d3Minus
        real*8 d1,d2,d3,d6
        real*8 nu21,nu31,nu32,delta,minDamage,crack_initiation
        real*8 etaT,etaL,phiC,kappa,lambda,ST,a,tN,tT,tL
        real*8 FI_LT,FI_LC,FI_MT,FI_MC
        real*8 dState1,dState2,dState3,dState6
        real*8 dState1Old,dState2Old,dState3Old,dState6Old
        integer :: i,j
        ! output variables
        real*8, dimension(6), intent(inout) :: stress
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew
        ! external functions
        real*8, external :: triangular,mccauley

        ! Parameters needed for failure criteria calculation
        etaL = 0.280622004043348d0
        phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)*((SL/XC)+etaL)))/
     1         (2.0d0*((SL/XC)+etaL)))  !Eq.12 Maimi
        ST = (0.5d0*(((2.0d0*sin(alpha0)**2.0d0)-1.0d0)*SL)/
     1       (((1.0d0-sin(alpha0)**2.0d0)**0.5d0)*sin(alpha0)*etaL)) !Eq.12 CLN
        etaT = (etaL*ST)/SL  !Eq.10 CLN
        kappa = (ST**2.0d0-YT**2.0d0)/(ST*YT)  !Eq.43 CLN
        lambda = ((2.0d0*etaL*ST)/SL)-kappa  !Eq.45 CLN 

        ! Load old dState values
        dState1Old = stateOld(k,13)
        dState2Old = stateOld(k,14)

        ! Calculate stress using an undamaged stiffness matrix
        trialStress = matmul(C(1:6,1:6), strain(1:6))

        ! Initialize damage activation functions (failure indices)
        FI_LT = 0.0d0
        FI_LC = 0.0d0
        FI_MT = 0.0d0
        FI_MC = 0.0d0

        ! Calculate longitudinal failure indices
        if (trialStress(1) > 0.d0) then
          ! If stress in fiber direction is positive calculate the
          ! longitudinal failure index
          FI_LT = strain(1)/(XT/e11) ! Eq. 54 CLN
        else if (trialStress(1) < 0.d0) then
          ! If stress is negative first rotate stresses into
          ! misalignment coordinate frame using rotate_stress function
          call rotate_stress(trialStress,phiC,XC,g12,trialStressP)
          ! Call fail_cln function with rotated stresses to determine
          ! longitudinal compressive failure index
          call catalanotti(trialStressP,ST,SL,etaL,etaT,lambda,kappa,
     1                     FI_LT,FI_LC,a,tN,tT,tL)
        end if

        ! Calculate scalar damage variables assuming triangular
        ! damage dissipation.
        if (FI_LT > 1.0d0) then ! Longitudinal tensile damage
          d1Plus= triangular(XT,e11,G1Plus,strain(1),lch)
        else
          d1Plus = 0.0d0
        end if

        if (FI_LC > 1.0d0) then ! Longitudinal compressive damage
            d1Minus = triangular(XC,e11,G1minus,strain(1),lch)
        else
            d1Minus = 0.0d0
        end if

        ! Calculate damage variables
        if (trialStress(1) /= 0.0d0) then ! Longitudinal damage
          d1 = d1Plus*(mccauley(trialStress(1))/abs(trialStress(1))) +
     1         d1Minus*(mccauley(-trialStress(1))/abs(trialStress(1)))
        else
          d1 = 0.0d0
        end if

        ! Calculate damage state variables. Damage state must be
        ! retained between iterations.
        dState1 = min(0.999d0,max(d1,dState1Old))

        ! Calculate damaged stiffness tensor
        ! Recalculate additional Poisson's ratio terms
        nu21 = nu12*(e22/e11)
        nu31 = nu13*(e33/e11)
        nu32 = nu23*(e33/e22)
  
        delta = 1.0d0/(e11*e22 - e22**2.0d0*nu12**2.0d0 +
     1          dState1*e22**2.0d0*nu12**2.0d0 - e11*e33*nu23**2.0d0 -
     2          e22*e33*nu13**2.0d0 + dState1*e22*e33*nu13**2.0d0 -
     3          2.0d0*e22*e33*nu12*nu13*nu23 + 2.0d0*dState1*e22*e33*
     4          nu12*nu13*nu23)

        C = 0.0d0 ! set all elements in array equal to zero
        C(1,1) = -((- e33*e11**2.0d0*nu23**2.0d0 + e22*e11**2.0d0)*
     1           (dState1 - 1.0d0))*delta
        C(2,1) = -((e11*nu12*e22**2.0d0 + e11*e33*nu13*nu23*e22)*
     1           (dState1 - 1.0d0))*delta
        C(3,1) = -((dState1 - 1.0d0)*(e11*e22*e33*nu13 +
     1           e11*e22*e33*nu12*nu23))*delta
        C(1,2) = (e22*(e11*e22*nu12 + e11*e33*nu13*nu23 -
     1           dState1*e11*e22*nu12 - dState1*e11*e33*nu13*nu23))*
     2           delta
        C(2,2) = (e22*(e11*e22 - e22*e33*nu13**2.0d0 +
     1           dState1*e22*e33*nu13**2.0d0))*delta
        C(3,2) = (e22*(e11*e33*nu23 + e22*e33*nu12*nu13 -
     1           dState1*e22*e33*nu12*nu13))*delta
        C(1,3) = (e22*e33*(e11*nu13 - dState1*e11*nu13 +
     1           e11*nu12*nu23 - dState1*e11*nu12*nu23))*delta
        C(2,3) = (e22*e33*(e11*nu23 + e22*nu12*nu13 -
     1           dState1*e22*nu12*nu13))*delta
        C(3,3) = (e22*e33*(e11 - e22*nu12**2.0d0 +
     1           dState1*e22*nu12**2.0d0))*delta
        C(4,4) = g12
        C(5,5) = g23
        C(6,6) = g13

        ! Calculate transverse failure indices using fail_cln function
        call catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                   FI_MT,FI_MC,a,tN,tT,tL)

        crack_initiation = stateOld(k,18)
        if ((FI_MT > 1.0d0).or.(FI_MC > 1.0d0)) then
          if (crack_initiation == 0.0d0) then
            stateNew(k,18) = 1.0d0
            do i = 7,12
              j = i - 6
              stateNew(k,i) = strain(j)
            end do
            
            stateNew(k,15) = tL
            stateNew(k,16) = tN
            stateNew(k,17) = tT
            
            stress = matmul(C(1:6,1:6), strain(1:6))
            
          else ! dont call smeared_crack at crack initiation
            
            do i = 7,12
              j = i - 6
              e0(j) = stateOld(k,i)
            end do

            tbar_cr(1) = stateOld(k,15)
            tbar_cr(2) = stateOld(k,16)
            tbar_cr(3) = stateOld(k,17)

            call smeared_crack(strain,e0,a,tbar_cr,C,G1Plus,G2Plus,lch,
     1                         e22,stress,dState2)
          end if

        else

          stress = matmul(C(1:6,1:6), strain(1:6))

        end if

        ! Record scalar damage variables state array
        stateNew(k,13) = dState1
        stateNew(k,14) = dState2
        ! Record failure indices in state array
        stateNew(k,21) = FI_LT
        stateNew(k,22) = FI_LC
        stateNew(k,23) = FI_MT
        stateNew(k,24) = FI_MC

        ! Check damage state of each element. If element is damaged
        ! completely in each direction it is flagged for deletion
        minDamage = min(dState1,dState2)
        stateNew(k,19) = minDamage
        if (minDamage > 0.99d0) stateNew(k,20) = 0.0d0
      end subroutine cdm

      subroutine rotate_stress(trialStress,phiC,XC,g12,trialStressP)
        ! Purpose: Rotate stresses to misaligned coordinate for
        !   calculation of longitudinal compressive failure index.
        ! Variable dictionary:
        ! trialStress = stress calculated using the undamaged stiffness
        !   tensor.
        ! trialStressP = trialStress rotated into the misalignment
        !   coordinate frame to calculate the failure index for
        !   longitudinal compression.
        ! trialStressT = trialStress rotated by in-plane angle theta.
        ! phiC = initial misalignment angle for compressive failure
        !   calculation
        ! g12 = longitudinal compressive failure stress
        ! phi = angle of kink band
        ! theta = kinking plane angle
        ! phi0 = initial misalignment angle
        ! m = cosine of theta
        ! n = sine of theta
        ! u = cosine of phi
        ! v = sine of phi
        ! tS4EQ0 = boolean - checks if trialStress(4) = 0
        ! tS6EQ0 = boolean - checks if trialStress(6) = 0
        ! gammaM = shear stress induced misalignment angle
        ! gammaMC = gammaM under axial compression loading
        ! eps = tolerance for conditional statements

        implicit none
        ! input variables
        real*8, intent(in) :: phiC,XC,g12
        real*8, dimension(6), intent(in) :: trialStress
        ! local variables
        real*8, dimension(6) :: trialStressT
        real*8  m,n,u,v,gamma0,phi,theta
        real*8  gammaMC, gammaM, phi0
        real*8, PARAMETER :: eps=1.d-8
        LOGICAL tS4EQ0, tS6EQ0
        ! output variables
        real*8, dimension(6), intent(out) :: trialStressP

        ! first determine kink plane angle theta
        tS4EQ0 = (abs(trialStress(4)-0.d0) < eps)
        tS6EQ0 = (abs(trialStress(6)-0.d0) < eps)

        if (tS4EQ0.and.tS6EQ0) then
          if (abs(trialStress(2)-trialStress(3)) < eps) then
            theta = atan(1.0d0) ! pi/4
          else
            theta = 0.5d0*atan((2.0d0*trialStress(5))/(trialStress(2)-
     1              trialStress(3))) !Eq.55 CLN
          end if 
        else
          if (abs(trialStress(4)-0.d0) < eps) then
            theta = 2.0d0*atan(1.0d0) ! pi/2
          else 
            theta = atan(trialStress(6)/trialStress(4)) !Eq. 56 CLN
          end if
        end if

        ! Rotate stresses by angle theta
        m = cos(theta)
        n = sin(theta)
        trialStressT(1) = trialStress(1)
        trialStressT(2) = trialStress(2)*m**2 +
     1                    2.0d0*trialStress(5)*m*n +
     2                    trialStress(3)*n**2
        trialStressT(3) = trialStress(3)*m**2 -
     1                    2.0d0*trialStress(5)*m*n +
     2                    trialStress(2)*n**2
        trialStressT(4) = trialStress(4)*m + trialStress(6)*n
        trialStressT(5) = trialStress(5)*(m**2 - n**2) -
     1                    trialStress(2)*n*m + trialStress(3)*n*m
        trialStressT(6) = trialStress(6)*m - trialStress(4)*n

        ! Determine kink band angle phi
        gammaMC = (sin(2*phiC)*XC)/(2*g12)  ! eq 74
        phi0 = phiC - gammaMC  ! eq 75
        gammaM = ((phi0*g12 + abs(trialStressT(4)))/
     1           (g12+trialStressT(1)-trialStressT(2)) - phi0)  ! eq 81

        if (trialStress(4) >= 0.d0) then ! eq 77
            phi = phi0 + gammaM
        else
            phi = -1.0*(phi0+gammaM)
        end if

        ! Rotate stresses by angle phi
        u = cos(phi)
        v = sin(phi)
        trialStressP(1) = trialStressT(1)*u**2 +
     1                    2.0d0*trialStressT(4)*u*v +
     2                    trialStressT(2)*v**2
        trialStressP(2) = trialStressT(2)*u**2 -
     1                    2.0d0*trialStressT(4)*v*u +
     2                    trialStressT(1)*v**2     
        trialStressP(3) = trialStressT(3)      
        trialStressP(4) = trialStressT(4)*(u**2 -v**2) +
     1                    trialStressT(2)*v*u - trialStressT(1)*v*u
        trialStressP(5) = trialStressT(5)*u - trialStressT(6)*v     
        trialStressP(6) = trialStressT(6)*u + trialStressT(5)*v

      end subroutine rotate_stress

      pure function triangular(X,E,G,tStrain,lch) result(dam)
        ! Purpose: Implement triangular damage dissipation
        ! Variable dictionary:
        ! X = failure stress (material direction unspecified)
        ! E = modulus
        ! G = energy release rate (NOT shear modulus)
        ! tStrain = scalar trial strain value
        ! lch = characteristic element length
        ! delta_0 = strain at yield
        ! delta_c = strain at final failure
        ! dam = scalar damage variable

        implicit none
        ! input variables
        real*8, intent(in) :: X,E,G,tStrain,lch
        ! local variables
        real*8 delta_0, delta_c,dam

        ! Calculate scalar damage variable
        if (tStrain == 0.0d0) then ! no strain = no damage
          dam = 0.0d0 
        else
          delta_0 = X/E
          delta_c = (2.0d0*G)/(X*lch)
          if (delta_c < delta_0) then  !Needed for d.nq.0-0
            delta_c = 1.1d0*delta_0   !Adjustable for stability
          end if
          dam = (delta_c*(abs(tStrain)-delta_0))/
     1          (abs(tStrain)*(delta_c-delta_0))
          if (dam < 0.0d0) then
              dam = 0.0d0
          end if
        end if
        dam = min(0.999,dam)

      end function triangular

      pure function mccauley(value)
        ! Purpose: implements the McCauley operator
        ! Variable dictionary:
        ! value = variable to execute mccauley brackets operation on
        implicit none
        real*8, intent(in) :: value
        real*8 :: mccauley
        mccauley = (value+abs(value))/2.0d0
      end function mccauley

      subroutine catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                       FI_MT,FI_MC,aFail,tN,tT,tL)
        ! Purpose: Implements transverse (and longitudinal compression)
        !   failure indices from Catalanotti et al. (2013)
        ! Variable Dictionary:
        ! trialStress = stress calculated using the undamaged
        !   stiffness tensor.
        ! FI_MT = tensile failure index
        ! FI_MC = compressive failure index
        ! etaT = friction coefficient in the transverse direction
        ! etaL = friction coefficient in the longitudinal direction
        ! kappa = parameter used to calculate failure indices
        ! lambda = parameter used to calculate failure indices
        ! ST = in-situ (ideally) transverse shear strength
        ! SL = in-situ (ideally) longitudinal shear strength
        ! trialFI_MT = stores trial values of tensile failure index
        ! trialFI_MC = stores trial values of compression failure index
        ! aFail_MT = angle of failure plane in tension
        ! aFail_MC = angle of failure plane in compression
        ! pi = variable to store the value of pi
        ! tN = traction normal to failure plane
        ! tT = traction transverse to failure plane
        ! tL = traction aligned with failure plane
        ! aR = angle in radians
        ! a = angle in degrees

        implicit none
        ! input variables
        real*8, dimension(6), intent(in) :: trialStress
        real*8, intent(in) :: ST,SL,etaL,etaT,lambda,kappa
        ! local variables
        real*8 :: trialFI_MT,trialFI_MC,pi,aR
        integer :: a
        ! output variables
        real*8, intent(out) :: FI_MT,FI_MC,aFail,tN,tT,tL
    
        pi = 4.0d0*atan(1.0d0) ! determine value of pi
    
        FI_MT = 0.0d0 ! initialize failure indices
        FI_MC = 0.0d0
        a = 0 ! initial failure plane angle

        ! iterate over angles to determine angle which maximizes FIs
        do while (a <= 180) ! iterate over angles from 0 to 180 degrees
          aR= a*(pi/180.0d0)  ! angle in radians
          ! Calculate tractions on failure plane: Eq 3 CLN (and 59-61)
          tN = trialStress(2)*cos(aR)**2.0d0 + 2.0d0*trialStress(5)*
     1         sin(aR)*cos(aR) + trialStress(3)*sin(aR)**2.0d0
          tT = -1.0d0*cos(aR)*sin(aR)*(trialStress(2)-trialStress(3))+
     1         (trialStress(5)*(cos(aR)**2.0d0 - sin(aR)**2.0d0))
          tL = trialStress(4)*cos(aR) + trialStress(6)*sin(aR)
      
          ! Calculate value of failure indices at current angle
          if (tN >= 0.0d0) then
            trialFI_MT = (tN/ST)**2.0d0 +
     1                   (tL/SL)**2.0d0 + (tT/ST)**2.0d0 + lambda*
     2                   (tN/ST)*(tL/SL)**2.0d0 + kappa*(tN/ST) ! Eq. 42 CLN
            trialFI_MC = 0.0d0
          else
            trialFI_MC = (tL/(SL-etaL*tN))**2.0d0 +
     1                   (tT/(ST-etaT*tN))**2.0d0 ! Eq. 5 CLN
            trialFI_MT = 0.0d0
          end if

          ! Update max values if current value > max value
          if (trialFI_MT > FI_MT) then
            FI_MT = trialFI_MT
            aFail = aR ! record failure plane 
          end if
          if (trialFI_MC > FI_MC) then
            FI_MC = trialFI_MC
            aFail = aR ! record failure plane 
          end if
      
          ! Update angle
          a = a + 1
        end do
      end subroutine catalanotti

      subroutine smeared_crack(eTotal,e0,alpha,tbar_cr,C,GIc,GIIc,lch,
     1                         E22,stress,d)
        implicit none
        ! input variables
        real*8, dimension(6,6), intent(in) :: C
        real*8, dimension(6), intent(in) :: eTotal, e0
        real*16, dimension(3), intent(in) :: tbar_cr
        real*8, intent(in) :: alpha, GIc, GIIc, lch, E22
        ! local variables
        real*16, dimension(6) :: ec_cr, eTotal_cr
        real*16, dimension(3,3) :: M, Minv
        real*8, dimension(6,6) :: R_gc, R_cg
        real*8, dimension(3,6) :: H
        real*8, dimension(6) :: ec_cr_new, ec, e0_cr
        real*8, dimension(6) :: sc, sTotal, s
        real*8, dimension(3,3) ::  R, RT
        real*8, dimension(3) :: n1, n2, n3, w_cr, f, res, ec_cr_delta
        real*8, dimension(3) :: t_cr, t
        real*8 :: A, B, diff, tol, beta, lambda, wmf, tbar_cr_norm
        real*8 :: w_cr_s, tbar_cr_s, p, q, eta
        integer i
        ! output variables
        real*8, dimension(6), intent(inout) :: stress
        real*8, intent(out) :: d
        ! external functions
        real*8, external :: mccauley, kdelta
        interface
          function matinv3(A) result (B)
            real*16, intent(in) :: A(3,3)   ! Matrix
            real*16             :: B(3,3)   ! Inverse matrix
          end function
        end interface
        interface
          function jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,
     1                       E22,eTotal_cr) result(A)
            real*8, dimension(6,6), intent(in) :: C
            real*16, dimension(6), intent(in) :: ec_cr, eTotal_cr
            real*16, dimension(3), intent(in) :: tbar_cr
            real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
            real*16, dimension(3,3) :: A
          end function
        end interface
        
        interface
          function rot(angle) result(A)
            real*8, intent(in) :: angle
            real*8, dimension(6,6) :: A
          end function
        end interface

        ! mixed mode interaction parameter
        eta = 1.45d0
        
        ! crack plane vectors (n2 is normal to crack plane)
        p = cos(alpha)
        q = sin(alpha)
        n1 = (/1.0d0, 0.0d0, 0.0d0 /)
        n2 = (/0.0d0, p, q /)
        n3 = (/0.0d0, -q, p /)

        write(*, *) tbar_cr
        ! tbar_cr = tractions on the failure plane at damage onset
        tbar_cr_s = (tbar_cr(1)**2 + tbar_cr(2)**2)**0.5
        ! compute the 2 norm of the tbar vector
        tbar_cr_norm = (tbar_cr(1)**2 + tbar_cr(2)**2 +
     1                 tbar_cr(3)**2.0d0)**0.5d0
        
        ! define 6x6 rotation matrices
        R_cg = rot(alpha) ! crack to global
        R_gc = rot(-alpha) ! global to crack
        ! define 3x3 rotation matrices
        R = 0.0d0
        R(1,1) = 1.0d0
        R(2,2) = p
        R(2,3) = q
        R(3,2) = -q
        R(3,3) = p
        RT = transpose(R)

        ! define matrix to calculate tractions on failure plane
        H = 0.0d0
        H(1,4) = p
        H(1,6) = q
        H(2,2) = p
        H(2,5) = q
        H(3,3) = q
        H(3,5) = p

        ! rotate total strain tensor to crack coordinate system
        eTotal_cr = matmul(R_gc(1:6,1:6), eTotal(1:6))
        ! rotate strain at crack initiation to crack coord system
        e0_cr = matmul(R_gc(1:6,1:6), e0(1:6))

        ! calculate initial trial strain
        ec_cr = eTotal_cr - e0_cr
        i = 0 ! iteration counter
        diff = 1.0d8 ! initial value for while loop
        tol = 1.0d-10 ! convergence criterion
        do while (diff.gt.tol)
          ! calculate first term of residual
          ! compute total stress
          sTotal = matmul(C(1:6,1:6), eTotal(1:6))
          ! rotate cracking strains to global coord system
          ec = matmul(R_cg(1:6,1:6), ec_cr(1:6))
          ! compute cracking stress in global coords
          sc = matmul(C(1:6,1:6), ec(1:6))
          ! subtract cracking stress from total stress
          s = sTotal - sc
          ! compute tractions
          t = matmul(H(1:3,1:6), s(1:6))
          ! compute f term of the residual
          f = matmul(RT(1:3,1:3), t(1:3))

          ! estimate displacement jump at crack
          w_cr(1) = 2.0d0*ec_cr(4)*lch ! e_12
          w_cr(2) = ec_cr(2)*lch ! e_22
          w_cr(3) = 2.0d0*ec_cr(5)*lch ! e_23

          ! second term of the residual according to the cohesive law
          ! compute lambda 
          lambda = (w_cr(1)**2.0d0 + mccauley(w_cr(2))**2.0d0 +
     1             w_cr(3)**2.0d0)**0.5d0
          ! compute the shear components of the displacement jump
          w_cr_s = (w_cr(1)**2.0d0 + w_cr(3)**2.0d0)**0.5d0  
          beta = w_cr_s/(w_cr_s+w_cr(2))
          ! compute Benzeggagh-Kenane parameters A and B
          A = GIIc - GIc
          B = (tbar_cr_s*beta)/(beta*(tbar_cr_s-tbar_cr(2))+tbar_cr(2))
          wmf = (2.0d0*(GIc + A*B**eta))/tbar_cr_norm

          ! compute the damage parameter d
          d = max(0.0d0,min((lambda/wmf),1.0d0))

          ! compute the tractions on the fracture plane
          do i = 1,3
            t_cr(i) = ((1.0d0-d)/d)*(w_cr(i)/wmf)*tbar_cr(i)-
     1                kdelta(i,2)*(mccauley(-w_cr(2))/w_cr(2))*
     2                (((1.0d0-d)/d)*(w_cr(i)/wmf)*tbar_cr(i)-
     3                E22*(eTotal_cr(4)- ec_cr(4)))
          end do

          ! calculate the residual res 
          res = f - t_cr
          ! calculate the Jacobian matrix M
          M = jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,E22,
     1                 eTotal_cr)
          ! invert M and solve for ec_cr_delta
          Minv = matinv3(M)
          ec_cr_delta = matmul(-Minv(1:3,1:3), res(1:3))
          ! calculate new crack strain values and determine error
          ec_cr_new(1) = ec_cr(1)
          ec_cr_new(2) = ec_cr(2) + ec_cr_delta(2)
          ec_cr_new(3) = ec_cr(3)
          ec_cr_new(4) = ec_cr(4) + ec_cr_delta(1)
          ec_cr_new(5) = ec_cr(5) + ec_cr_delta(3)
          ec_cr_new(6) = ec_cr(6)
          diff = maxval(abs(ec_cr_new - ec_cr))
          ! increment counter and reassign crack strain
          ec_cr = ec_cr_new
          i = i + 1
        end do

        sTotal = matmul(C(1:6,1:6), eTotal(1:6))
        ! rotate cracking strains to global coord system
        ec = matmul(R_cg(1:6,1:6), ec_cr(1:6))
        ! compute cracking stress in global coords
        sc = matmul(C(1:6,1:6), ec(1:6))
        ! subtract cracking stress from total stress
        stress = sTotal - sc
      end subroutine smeared_crack

      pure function matinv3(A) result(B)
        ! Adapted from the matrix inversion page on the Fortran Wiki
        ! see: http://fortranwiki.org/fortran/show/Matrix+inversion
        ! Performs a direct calculation of the inverse of a 3×3 matrix.
        implicit none
        real*16, intent(in) :: A(3,3)   ! Matrix
        real*16             :: B(3,3)   ! Inverse matrix
        real*16             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1.0d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     1           - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     2           + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
      end function matinv3

      pure function rot(angle) result(A)
        implicit none
        ! input variables
        real*8, intent(in) :: angle
        ! local variables
        real*8 :: m,n
        ! output variables
        real*8, dimension(6,6) :: A
        ! calculate sine and cosine of angle
        m = cos(angle)
        n = sin(angle)
        ! define rotation matrix
        A = 0.0d0
        A(1,1) = 1.0d0
        A(2,1) = m**2
        A(2,2) = n**2
        A(2,5) = 2*m*n
        A(3,2) = n**2
        A(3,3) = m**2
        A(3,5) = -2*m*n
        A(4,4) = m 
        A(4,6) = n
        A(5,2) = -m*n
        A(5,3) = m*n
        A(5,5) = m**2 - n**2
        A(6,4) = -n
        A(6,6) = m
      end function rot

      pure function kDelta(i,j)
        implicit none
        integer, intent(in) :: i, j
        real*8 :: kDelta
        if (i /= j) then
          kDelta = 0.0d0
        else
          kDelta = 1.0d0
        end if
      end function kDelta

      pure function jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,
     1                       E22,eTotal_cr) result(A)
        implicit none
        ! input variables
        real*8, dimension(6,6), intent(in) :: C
        real*16, dimension(6), intent(in) :: ec_cr, eTotal_cr
        real*16, dimension(3), intent(in) :: tbar_cr
        real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
        ! local variables
        real*16 :: m,n, term1, term2, term3, term4, term5, term6, term7
        ! output variables
        real*16, dimension(3,3) :: A
        ! define cosine and sine of fracture plane angle for brevity
        m = cos(alpha)
        n = sin(alpha)
        ! terms1-7 are repeated parts of the Jacobian matrix. They have
        ! been defined here to reduce the length of the equations.
        term1 = (GIc + ((tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     3    tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     5    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)
        term2 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*
     1    ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*
     2    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
     3    sign(1.0d0,ec_cr(2)*lch)**2.0d0)
        term3 = ((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     1    lch) + 1.0d0)**2.0d0*(tbar_cr(1)**2.0d0*
     2    sign(1.0d0,tbar_cr(2))**2.0d0*
     3    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     4    sign(1.0d0,tbar_cr(1))**2.0d0*
     5    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     6    sign(1.0d0,tbar_cr(1))**2.0d0*
     7    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
     8    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     9    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1*
     1    term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*(ec_cr(2)*
     2    tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     4    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     5    2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     8    (GIIc - GIc)*(tbar_cr(1)**2.0d0*
     9    sign(1.0d0,tbar_cr(2))**2.0d0*
     1    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     2    sign(1.0d0,tbar_cr(1))**2.0d0*
     3    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     4    sign(1.0d0,tbar_cr(1))**2.0d0*
     5    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*
     6    term2**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0*
     7    tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*ec_cr(4)**2.0d0*
     1    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
     3    ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     5    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*
     6    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
     7    ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
     8    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
     9    (sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     1    sign(1.0d0,tbar_cr(3))*term1**2.0d0*(4.0d0*
     2    ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
     3    tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     4    tbar_cr(1)**2.0d0)**2.0d0))
        term4 = ((0.25d0*(tbar_cr(1)**2.0d0*
     1    sign(1.0d0,tbar_cr(2))**2.0d0*
     2    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     3    sign(1.0d0,tbar_cr(1))**2.0d0*
     4    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     5    sign(1.0d0,tbar_cr(1))**2.0d0*
     6    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*
     7    term2**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     8    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1) - 1.0d0)
        term5 = ((4.0d0*ec_cr(4)*lch**2.0d0*
     1    (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 +
     2    abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0))/ (term1*(16.0d0*
     3    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     4    2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*
     5    ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     6    lch)**2.0d0)**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     7    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     9    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     1    2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     3    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     4    (GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
     5    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     6    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     7    ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     8    lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*
     9    lch*abs(ec_cr(2)*lch))**(1.0d0/2.0d0)*(2.0d0*
     1    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
     4    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     6    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*
     7    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
     8    ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
     9    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
     1    (term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 +
     2    4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
     3    ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0))
        term6 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     1    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     2    abs(ec_cr(2)*lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 +
     3    abs(ec_cr(2)*lch)**2.0d0)
        term7 = ((0.25d0*(abs(tbar_cr(1))**2.0d0 +
     1    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     2    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     3    ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     4    lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*
     5    lch*abs(ec_cr(2)*lch))**(1.0d0/2.0d0))/term1 - 1.0d0)
        A(1,1) = C(6,6)*(m**2.0d0 - 1.0d0) - C(4,4)*m**2.0d0 +
     1    (4.0d0*lch*tbar_cr(1)*term7)/term6**0.5d0 - (64.0d0*
     2    ec_cr(4)**2.0d0*lch**3.0d0*tbar_cr(1)*term7)/
     3    term6**1.5d0 + (4.0d0*ec_cr(4)*lch*tbar_cr(1)*term5)/
     4    term6**0.5d0
        A(1,2) = (ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
     1    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
     2    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     3    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     4    sign(1.0d0,tbar_cr(1))**2.0d0*
     5    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     6    sign(1.0d0,tbar_cr(1))**2.0d0*
     7    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
     8    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     9    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1*
     1    term2) - (4.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
     2    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/
     3    (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0) - (4.0d0*
     4    ec_cr(4)**3.0d0*eta*lch*tbar_cr(1)*tbar_cr(2)*(ec_cr(2)*
     5    tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     2    (GIIc - GIc)*(tbar_cr(1)**2.0d0*
     3    sign(1.0d0,tbar_cr(2))**2.0d0*
     4    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     5    sign(1.0d0,tbar_cr(1))**2.0d0*
     6    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     7    sign(1.0d0,tbar_cr(1))**2.0d0*
     8    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(2.0d0*
     9    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
     1    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
     3    4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - ec_cr(2)**2.0d0*
     6    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) +
     9    8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 +
     1    2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) +
     2    8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
     3    tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*
     4    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*
     5    term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 +
     6    4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
     7    ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)
        A(1,3) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*
     1    tbar_cr(1))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 +
     2    16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     3    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*
     4    ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0
        A(2,1) = (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term5)/
     1    term6**0.5d0 - ((0.5d0*abs(ec_cr(2)*lch) - 0.5d0*
     2    ec_cr(2)*lch)*((2.0d0*ec_cr(2)*lch*tbar_cr(2)*term5)/
     3    term6**0.5d0 - (32*ec_cr(4)*ec_cr(2)*lch**3.0d0*
     4    tbar_cr(2)*term7)/term6**1.5d0))/(ec_cr(2)*lch) - (32*
     5    ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*term7)/
     6    term6**1.5d0
        A(2,2) = n**2.0d0*(C(2,3)*m**2.0d0 + 2.0d0*C(5,5)*
     1    m**2.0d0 + C(3,3)*n**2.0d0) - m**2.0d0*(C(2,2) - C(2,2)*
     2    n**2.0d0 + C(2,3)*n**2.0d0 + 2.0d0*C(5,5)*n**2.0d0) +
     3    (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*
     4    (ec_cr(2) - eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*
     5    m**2.0d0 + 2.0d0*eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*
     6    tbar_cr(2)*term4)/term2**0.5d0))/ec_cr(2) - (0.5d0*
     7    (sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22 - (2.0d0*lch*
     8    tbar_cr(2)*term4)/term2**0.5d0 - (2.0d0*ec_cr(2)*lch*
     9    tbar_cr(2)*term3)/term2**0.5d0 + (2.0d0*ec_cr(2)**2.0d0*
     1    lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) +
     2    1.0d0)**2.0d0*term4)/(sign(1.0d0,ec_cr(2)*lch)*
     3    term2**1.5d0)))/sign(1.0d0,ec_cr(2)*lch) + (2.0d0*lch*
     4    tbar_cr(2)*term4)/term2**0.5d0 + (0.5d0*
     5    (sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) -
     6    eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*m**2.0d0 + 2.0d0*
     7    eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*
     8    term4)/term2**0.5d0))/(ec_cr(2)*sign(1.0d0,ec_cr(2)*
     9    lch)) + (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term3)/
     1    term2**0.5d0 - (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*
     2    tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
     3    term4)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
        A(2,3) = 0.5d0*C(3,3)*2.0d0*m*n - 0.5d0*C(2,3)*2.0d0*m*
     1    n + 0.5d0*C(2,2)*2.0d0*m*n*m**2.0d0 - 0.5d0*C(3,3)*2.0d0*
     2    m*n*m**2.0d0 + (48*ec_cr(2)*ec_cr(5)*lch**3.0d0*
     3    tbar_cr(2))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 +
     4    16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     5    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*
     6    ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0 - (16.0d0*
     7    ec_cr(5)*lch**2.0d0*tbar_cr(2)*abs(ec_cr(2))*abs(lch))/
     8    (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     9    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     1    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     2    abs(ec_cr(2))*abs(lch))**1.5d0
        A(3,1) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*
     1    tbar_cr(3))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     2    ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     3    lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*
     4    lch*abs(ec_cr(2)*lch))**1.5d0 + (4.0d0*ec_cr(4)*ec_cr(2)*
     5    ec_cr(5)*eta*lch*tbar_cr(2)*tbar_cr(3)*(ec_cr(2)*
     6    tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     9    2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     3    (GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
     4    abs(tbar_cr(2))**2.0d0 +
     5    abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0)*(2.0d0*
     6    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
     9    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) +
     3    8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 +
     4    2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) +
     5    8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
     6    tbar_cr(2)))/(term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*
     7    tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
     8    tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     9    tbar_cr(1)**2.0d0)**2.0d0)
        A(3,2) = (4.0d0*ec_cr(5)*lch*tbar_cr(3)*term3)/
     1    term2**0.5d0 - m*n*(C(2,2) - C(2,2)*n**2.0d0 + C(2,3)*
     2    n**2.0d0 + 2.0d0*C(5,5)*n**2.0d0) - m*n*(C(2,3)*
     3    m**2.0d0 + 2.0d0*C(5,5)*m**2.0d0 + C(3,3)*n**2.0d0) -
     4    (4.0d0*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(3)*
     5    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/
     6    (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
        A(3,3) = (64.0d0*ec_cr(5)**2.0d0*lch**3.0d0*tbar_cr(3))/
     1    (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     3    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     4    abs(ec_cr(2))*abs(lch))**1.5d0 - C(5,5) - (4.0d0*lch*
     5    tbar_cr(3))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 +
     6    16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     7    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*
     8    ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**0.5d0 - C(2,2)*
     9    n**2.0d0*(n**2.0d0 - 1.0d0) + C(3,3)*n**2.0d0*(n**2.0d0 -
     1    1.0d0) + 2.0d0*C(5,5)*n**2.0d0 + (lch*tbar_cr(3)*
     2    (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 +
     3    abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0))/(GIc + (GIIc*
     4    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*
     5    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta))/
     6    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta -
     9    (GIc*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*
     1    eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*
     2    eta))/(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
     4    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)
      end function jacobian
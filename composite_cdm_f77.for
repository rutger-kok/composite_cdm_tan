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
        real*8, dimension(6) :: stress, strain
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
          do while (debugLoop /= 999)
            debugLoop = 1
          end do
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
        real*8, dimension(6) :: trialStress, trialStressP,e0
        real*16, dimension(3) :: tbar_cr
        real*8 d1Plus,d1Minus, d1
        real*8 nu21,nu31,nu32,delta,minDamage,crack_initiation
        real*8 etaT,etaL,phiC,kappa,lambda,ST,a,tN,tT,tL
        real*8 FI_LT,FI_LC,FI_MT,FI_MC
        real*8 dState1,dState2,dState1Old,dState2Old
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

            call smeared_crack(strain,e0,a,tbar_cr,C,G2Plus,G6,lch,
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
        real*16, dimension(3), intent(inout) :: tbar_cr
        real*8, intent(in) :: alpha, GIc, GIIc, lch, E22
        ! local variables
        real*16, dimension(6) :: ec_cr, eTotal_cr, e0_cr, ec
        real*16, dimension(3,3) :: M, Minv
        real*8, dimension(6,6) :: R_gc, R_cg
        real*8, dimension(3,6) :: H
        real*8, dimension(6) :: ec_cr_new
        real*8, dimension(6) :: sc, sTotal, s
        real*8, dimension(3,3) ::  R, RT
        real*8, dimension(3) :: n1, n2, n3, w_cr, f, res, ec_cr_delta
        real*8, dimension(3) :: t_cr, t
        real*16 :: tbar_cr_norm, tbar_cr_s, wmf, B
        real*8 :: A, diff, tol, beta, lambda
        real*8 :: w_cr_s, p, q, eta
        integer i, j
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

          function jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,
     1                      E22,eTotal_cr) result(M)
            real*8, dimension(6,6), intent(in) :: C
            real*16, dimension(6), intent(in) :: ec_cr, eTotal_cr
            real*16, dimension(3), intent(in) :: tbar_cr
            real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
            real*16, dimension(3,3) :: M
          end function

          function rot(angle) result(A)
            real*8, intent(in) :: angle
            real*8, dimension(6,6) :: A
          end function
        end interface

        ! mixed mode interaction parameter
        eta = 1.45d0
        ! compute Benzeggagh-Kenane parameter A
        A = GIIc - GIc
        
        ! crack plane vectors (n2 is normal to crack plane)
        p = cos(alpha)
        q = sin(alpha)
        n1 = (/1.0d0, 0.0d0, 0.0d0 /)
        n2 = (/0.0d0, p, q /)
        n3 = (/0.0d0, -q, p /)
        ! adjust tbar_cr for stability
        call stability_adjust(tbar_cr)
        ! tbar_cr = tractions on the failure plane at damage onset
        tbar_cr_s = (tbar_cr(1)**2.0d0 + tbar_cr(3)**2.0d0)**0.5d0
        ! compute the 2 norm of the tbar vector
        tbar_cr_norm = (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0 +
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
          ! adjust for stability
          call stability_adjust(ec)
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
          ! compute Benzeggagh-Kenane parameter B
          B = (tbar_cr_s*beta)/(beta*(tbar_cr_s-tbar_cr(2))+tbar_cr(2))
          wmf = (2.0d0*(GIc + A*B**eta))/tbar_cr_norm

          ! compute the damage parameter d
          d = max(0.0d0,min((lambda/wmf),1.0d0))

          ! compute the tractions on the fracture plane
          do j = 1,3
            t_cr(j) = ((1.0d0-d)/d)*(w_cr(j)/wmf)*tbar_cr(j)-
     1                kdelta(j,2)*(mccauley(-w_cr(2))/w_cr(2))*
     2                (((1.0d0-d)/d)*(w_cr(j)/wmf)*tbar_cr(j)-
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
        
      contains
      
        subroutine stability_adjust(A)
          ! Purpose: adjust values for stability. This function will
          ! replace all values in array A that are smaller than tol
          ! with 0.0d0. Note: adjusts arrays in-place.
          ! Variable dictionary
          ! tol = tolerance for adjustment
          ! n = size of input array (must be 1 dimensional)
          ! A = input array to adjust for numerical stability
          ! i = loop counter
          implicit none
          ! input variables
          real*16, intent(inout) :: A(:)
          ! local variables
          real*8, parameter :: tol = 1.0d-8
          integer :: n, i
          n = size(A)
          !output variables
          do i = 1,n
            if (abs(A(i)) < tol) then
              A(i) = 1.0d-8
            end if
          end do
        end subroutine stability_adjust
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

      function jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,
     1                  E22,eTotal_cr) result(M)
        implicit none
        ! input variables
        real*8, dimension(6,6), intent(in) :: C
        real*16, dimension(6), intent(in) :: ec_cr, eTotal_cr
        real*16, dimension(3), intent(in) :: tbar_cr
        real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
        ! local variables
        real*16 :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
        real*16 :: t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
        real*16 :: t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42
        real*16 :: t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55
        real*16 :: t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68
        real*16 :: t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81
        real*16 :: t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94
        real*16 :: t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105
        real*16 :: t106,t107,t108,t109,t110,t111,t112,t113,t114,t115
        real*16 :: t116,t117,t118,t119,t120,t121
        ! output variables
        real*16, dimension(3,3) :: M
        
        t2 = cos(alpha)
        t3 = sin(alpha)
        t4 = lch**2.0D0
        t5 = ec_cr(4)**2.0D0
        t6 = t4*t5*4.0D0
        t7 = ec_cr(2)**2.0D0
        t8 = t4*t7
        t9 = t6+t8
        t10 = sqrt(t9)
        t11 = ec_cr(2)*lch
        t12 = t10+t11
        t13 = 1.0D0/t12
        t14 = tbar_cr(1)**2.0D0
        t15 = tbar_cr(2)**2.0D0
        t16 = t14+t15
        t17 = sqrt(t16)
        t22 = abs(t11)
        t23 = t22/2.0D0
        t24 = (ec_cr(2)*lch)/2.0D0
        t18 = t23+t24
        t19 = abs(tbar_cr(1))
        t20 = abs(tbar_cr(2))
        t21 = abs(tbar_cr(3))
        t25 = t18**2.0D0
        t26 = ec_cr(5)**2.0D0
        t27 = t4*t26*4.0D0
        t28 = t6+t25+t27
        t29 = GIc*2.0D0
        t30 = GIIc-GIc
        t46 = t17-tbar_cr(2)
        t31 = t10*t13*t46
        t32 = t31+tbar_cr(2)
        t33 = 1.0D0/t32
        t34 = t10*t13*t17*t33
        t35 = t34**eta
        t36 = t30*t35*2.0D0
        t37 = t29+t36
        t38 = 1.0D0/t37
        t39 = 1.0D0/sqrt(t28)
        t40 = t19**2.0D0
        t41 = t20**2.0D0
        t42 = t21**2.0D0
        t43 = t40+t41+t42
        t44 = sqrt(t43)
        t45 = sqrt(t28)
        t47 = 1.0D0/t12**2.0D0
        t48 = 1.0D0/sqrt(t9)
        t49 = t10*t13*(t17-tbar_cr(2))
        t50 = t49+tbar_cr(2)
        t51 = 1.0D0/t50
        t52 = t10*t13*t17*t51
        t53 = t52**eta
        t54 = t30*t53*2.0D0
        t55 = t29+t54
        t56 = 1.0D0/t55
        t57 = ec_cr(2)*t4*t48
        t58 = lch+t57
        t59 = eta-1.0D0
        t60 = t52**t59
        t61 = lch/2.0D0
        t62 = (t11/abs(t11))
        t63 = (lch*t62)/2.0D0
        t64 = t61+t63
        t65 = t44*t45*t56
        t66 = t65-1.0D0
        t67 = 1.0D0/t28**(3.0D0/2.0D0)
        t68 = 1.0D0/t55**2.0D0
        t69 = ec_cr(4)*t4*t46*t47*4.0D0
        t70 = 1.0D0/t50**2.0D0
        t74 = ec_cr(4)*t4*t13*t46*t48*4.0D0
        t71 = t69-t74
        t72 = ec_cr(4)*t4*t13*t17*t48*t51*4.0D0
        t73 = ec_cr(4)*t4*t39*t44*t56*4.0D0
        t75 = t10*t13*t17*t70*t71
        t103 = ec_cr(4)*t4*t17*t47*t51*4.0D0
        t76 = t72+t75-t103
        t104 = eta*t30*t44*t45*t60*t68*t76*2.0D0
        t77 = t73-t104
        t78 = ec_cr(2)*lch*t39*t77*tbar_cr(2)
        t79 = t2**2.0D0
        t80 = t3**2.0D0
        t81 = 1.0D0/ec_cr(2)
        t82 = 1.0D0/lch
        t83 = t23-t24
        t84 = lch*t39*t66*tbar_cr(2)
        t85 = t18*t39*t44*t56*t64
        t86 = t10*t46*t47*t58
        t96 = ec_cr(2)*t4*t13*t46*t48
        t87 = t86-t96
        t88 = t10*t13*t17*t70*t87
        t89 = ec_cr(2)*t4*t13*t17*t48*t51
        t97 = t10*t17*t47*t51*t58
        t90 = t88+t89-t97
        t98 = eta*t30*t44*t45*t60*t68*t90*2.0D0
        t91 = t85-t98
        t92 = eTotal_cr(5)*t2*t3*2.0D0
        t93 = ec_cr(2)+t92-eTotal_cr(2)*t79-eTotal_cr(3)*t80
        t94 = E22*t93
        t95 = t94-ec_cr(2)*lch*t39*t66*tbar_cr(2)
        t99 = ec_cr(2)*lch*t18*t64*t66*t67*tbar_cr(2)
        t100 = t79-t80
        t101 = 1.0D0/t28
        t102 = ec_cr(2)*ec_cr(5)*lch*t4*t66*t67*tbar_cr(2)*4.0D0
        t105 = C(2,2)*t79
        t106 = C(2,3)*t80
        t107 = t105+t106
        t108 = t2*t107
        t109 = C(5,5)*t2*t80*2.0D0
        t110 = t108+t109
        t111 = C(2,3)*t79
        t112 = C(3,3)*t80
        t113 = t111+t112
        t114 = t3*t113
        t115 = C(5,5)*t3*t79*2.0D0
        t116 = t114+t115
        t117 = C(2,3)*t2*t3
        t118 = t117-C(3,3)*t2*t3
        t119 = t3*t118
        t120 = t119-C(5,5)*t2*t100
        t121 = C(2,2)*t2*t3

        M(1,1) = -C(4,4)*t79-C(6,6)*t80+lch*t39*tbar_cr(1)*
     1    (t38*t44*t45-1.0D0)*2.0D0+ec_cr(4)*lch*t39*tbar_cr(1)*
     2    (ec_cr(4)*t4*t38*t39*t44*4.0D0-eta*t30*1.0D0/t37**2.0d0*t44*
     3    t45*t60*(t72+t10*t13*t17*1.0D0/t32**2.0D0*t71-ec_cr(4)*t4*t17*
     4    t47*t51*4.0D0)*2.0D0)*2.0D0-lch*t4*t5*t66*t67*tbar_cr(1)*8.0D0
        M(1,2) = ec_cr(4)*lch*t39*t91*tbar_cr(1)*2.0D0-ec_cr(4)*lch*
     1    t18*t64*t66*t67*tbar_cr(1)*2.0D0
        M(1,3) = ec_cr(4)*ec_cr(5)*lch*t4*t66*t67*tbar_cr(1)*(-8.0D0)+
     1    ec_cr(4)*ec_cr(5)*lch*t4*t44*t56*t101*tbar_cr(1)*8.0D0
        M(2,1) = t78-t81*t82*t83*(t78-ec_cr(4)*ec_cr(2)*lch*t4*t66*t67*
     1   tbar_cr(2)*4.0D0)-ec_cr(4)*ec_cr(2)*lch*t4*t66*t67*tbar_cr(2)*
     2   4.0D0
        M(2,2) = t84-t99-t2*t110+t3*t116+t81*t82*t83*(E22-t84+t99-
     1   ec_cr(2)*lch*t39*t91*tbar_cr(2))-t81*t82*t95*(t61-t63)-1.0D0/
     2   ec_cr(2)**2.0D0*t82*t83*t95+ec_cr(2)*lch*t39*t91*tbar_cr(2)
        M(2,3) = -t102-t3*t120+t2*(t2*(t121-C(2,3)*t2*t3)-C(5,5)*t3*
     1    t100)+t81*t82*t83*(t102-ec_cr(2)*ec_cr(5)*lch*t4*t44*t56*t101*
     2    tbar_cr(2)*4.0D0)+ec_cr(2)*ec_cr(5)*lch*t4*t44*t56*t101*
     3    tbar_cr(2)*4.0D0
        M(3,1) = ec_cr(5)*lch*t39*t77*tbar_cr(3)*2.0D0-ec_cr(4)*
     1    ec_cr(5)*lch*t4*t66*t67*tbar_cr(3)*8.0D0
        M(3,2) = -t3*t110-t2*t116+ec_cr(5)*lch*t39*t91*tbar_cr(3)*
     1    2.0D0-ec_cr(5)*lch*t18*t64*t66*t67*tbar_cr(3)*2.0D0
        M(3,3) = -t3*(t2*(t117-t121)+C(5,5)*t3*t100)+t2*t120+lch*t39*
     1    t66*tbar_cr(3)*2.0D0-lch*t4*t26*t66*t67*tbar_cr(3)*8.0D0+lch*
     2    t4*t26*t44*t56*t101*tbar_cr(3)*8.0D0
      end function jacobian
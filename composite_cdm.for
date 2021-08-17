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
        ! framework for composite materials in Abaqus. Damage initiation
        ! is based on the work of Catalanotti et al. [1] while damage
        ! evolution is defined according to an amalgamation of the
        ! models developed by Maimi et al. [2][3] (for longitudinal
        ! tensile damage) and Tan et al. [4][5] (for matrix damage). 

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

        ! [4] W. Tan, B. G. Falzon, L. N. S. Chiu, and M. Price
        ! Predicting low velocity impact damage and
        ! Compression-After-Impact (CAI) behaviour of composite
        ! laminates
        ! Composites Part A 71 (2015) 212–226.
        ! http://doi.org/10.1016/j.compositesa.2015.01.025

        ! [5] B. G. Falzon, H. Liu, and W. Tan
        ! Comment on ‘A tensorial based progressive damage model for
        ! fibre reinforced polymers’
        ! Composite Structures 176 (2017) 877–882.
        ! http://doi.org/10.1016/j.compstruct.2017.06.011
  
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
        real*8 :: nu21,nu31,nu32,psi
        real*8 :: GL1Plus,GE1Plus,G1Minus,G2Plus,G6
        real*8 :: XT,XPO,XC,YT,YC,SL,ST,alpha0,stressPower,etaL,lch,t
        integer :: k,i,debugLoop

        ! Elastic constants orthotropic ply
        e11 = props(1) ! stiffness fiber 
        e22 = props(2) ! stiffness transverse  (in-plane)
        e33 = props(3) ! stiffness transverse  (out-of-plane)
        nu12 = props(4) ! Poisson's ratio 12 
        nu13 = props(5) ! Poisson's ratio 13 
        nu23 = props(6) ! Poisson's ratio 23 
        ! shear moduli multiplied by two (tensorial-engineering strain)
        g12 = 2.0d0*props(7) ! shear modulus 12  
        g13 = 2.0d0*props(8) ! shear modulus 13 
        g23 = 2.0d0*props(9) ! shear modulus 23 

        ! Ply strengths
        XT = props(10) ! tensile strength fiber 
        XPO = props(11) ! tensile pull-out strength
        XC = props(12) ! compressive strength fiber 
        YT = props(13) ! in-situ tensile strength transverse 
        YC = props(14) ! in-situ compressive strength transverse 
        SL = props(15) ! in-situ longitudinal shear strength
        ST = props(16) ! in-situ transverse shear strength

        alpha0 = props(17)*0.017453292519943295d0 ! Fracture angle
        etaL = props(18) ! longitudinal friction coefficient

        ! Fracture toughnesses (Gc)
        GL1Plus = props(19) ! tensile Gc fiber (linear)
        GE1Plus = props(20) ! tensile Gc fiber (exponential)
        G1Minus = props(21) ! comp. Gc fiber
        G2Plus = props(22) ! tensile Gc trans. 
        G6 = props(23) ! shear Gc

        ! Calculate remaining Poisson's ratios
        nu21 = nu12*(e22/e11) ! Poisson's ratio 21
        nu31 = nu13*(e33/e11) ! Poisson's ratio 31
        nu32 = nu23*(e33/e22) ! Poisson's ratio 32

        ! Calculate undamaged stiffness matrix (Eq. 2 [4])
        psi = (1.0d0 - nu12*nu21 - nu23*nu32 - nu31*nu13 -
     1         2.0d0*nu12*nu23*nu31)/(e11*e22*e33)

        C = 0.0d0 ! initialize stiffness matrix, set all elements = 0
        C(1,1) = (1.0d0 - nu23*nu32)/(e22*e33*psi)
        C(1,2) = (nu21 + nu31*nu23)/(e22*e33*psi)
        C(1,3) = (nu31 + nu21*nu32)/(e22*e33*psi)
        C(2,1) = (nu12 + nu13*nu32)/(e11*e33*psi)
        C(2,2) = (1.0d0 - nu31*nu13)/(e11*e33*psi)
        C(2,3) = (nu32 + nu31*nu12)/(e11*e33*psi)
        C(3,1) = (nu13 + nu12*nu23)/(e11*e22*psi)
        C(3,2) = (nu23 + nu13*nu21)/(e11*e22*psi)
        C(3,3) = (1.0d0 - nu12*nu21)/(e11*e22*psi)
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
            stress = matmul(C, strain)
            ! Assign calculated stresses to stressNew array
            do i = 1,6
              stressNew(k,i) = stress(i)
            end do

          end do
        else
        ! enabled only for debugging in Visual Studio
        !   t = stepTime - dt
        !   do while ((debugLoop.ne.999).and.(t == 0.0d0))
        !     debugLoop = 1
        !   end do

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

            lch = charLength(k)
            
            ! Call CDM subroutine
            call cdm(k,nblock,nstatev,strain,stateOld,C,e11,e22,e33,lch,
     1               nu12,nu13,nu23,g12,g13,g23,XT,XPO,XC,YT,YC,SL,ST,
     2               GL1Plus,GE1Plus,G1Minus,G2Plus,G6,etaL,defGradNew,
     3               stress,stateNew)

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
        close(95)
      end subroutine vumat

      subroutine cdm(k,nblock,nstatev,strain,stateOld,C,e11,e22,e33,lch,
     1               nu12,nu13,nu23,g12,g13,g23,XT,XPO,XC,YT,YC,SL,ST,
     2               GL1Plus,GE1Plus,G1Minus,G2Plus,G6,etaL,defGradNew,
     3               stress,stateNew)
        ! Purpose: Implements continuum damage mechanics framework
        !   from Maimi et al. (2007) for longitudinal tensile damage and
        !   from Tan et al. (2015) for transverse and longitudinal
        !   compressive damage.
        ! Variable Dictionary:
        ! lch = current characteristic length
        ! lch_fib = longitudinal damage characteristic length
        ! lch_mat = matrix damage characteristic element length
        ! trialStress = stress calculated using the undamaged
        !   stiffness tensor.
        ! stress_np = nonphysical stress calculated using Poisson
        !   degraded stiffness matrix
        ! stress_np_crack = stress_np roatted to crack coordinate frame
        ! strain_0_crack = strains at matrix damage onset in crack CSYS
        ! stress_0_crack = stress at matrix damage onset in crack CSYS
        ! stress_crack = final degraded stresses in crack CSYS
        ! dE1Plus,dL1Plus,dF1Plus,d1Plus = damage variables for
        !   longitudinal tension (see Maimi et al. [2][3])
        ! d1Minus = damage variable longitudinal compression
        ! dN (N=1,2,3,4,5,6) = damage variables quantifying
        !   damage in the 11,22,33,12,23,13 directions
        ! dNState (N=1,2,3,4,5,6) = state damage variables.
        !   Ensure damage state retention as loads fluctuate.
        ! dMat/dMatState = matrix damage variable
        ! D = damage variable matrix
        ! r1Plus,rE1Plus,rF1Plus,rL1Plus = elastic limits for
        !   longitudinal tensile damage (see Maimi [2][3])
        ! K1 = slope of the linear softening portion of the longitudinal
        !   tensile material response.
        ! A1Plus = parameter to ensure correct dissipation of energy
        !   with longitudinal tensile damage.
        ! FI_LT = fiber direction tensile failure index
        ! FI_LC = fiber direction compressive failure index
        ! FI_M = matrix failure index
        ! etaT = friction coefficient in the transverse direction
        ! kappa = parameter used to calculate failure indices
        ! lambda = parameter used to calculate failure indices
        ! R_GC,R_CG = transformation matrices (global to crack CSYS and
        !   crack to global CSYS)
        ! nuXY_d (X,Y = 1,2,3) = degraded Poisson's ratios
        ! a_fail = angle of matrix damage failure plane
        ! epsR_0 (equivalent strain at matrix failure onset)
        ! tN_0,tT_0,tL_0 = tractions on failure plane at damage onset
        ! epsN_0,epsT_0,epsL_0 = strains on fail plane at damage onset
        ! tR_0 (equivalent traction at matrix failure onset)
        ! Gamma (mixed mode critical strain energy release rate)
        ! epsR_f (equivalent strain at ultimate failure)
        ! Variables containing 'Old' are state variables from the
        !   previous iteration
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev,k
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld
        real*8, dimension(nblock,9), intent(in) :: defGradNew ! hardcoded dimensions assume 3D
        real*8, dimension(6,6), intent(inout) :: C
        real*8, dimension(6), intent(in) :: strain
        real*8, intent(in) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
        real*8, intent(in) :: XT,XPO,XC,YT,YC,SL,ST,etaL,lch
        real*8, intent(in) :: GL1Plus,GE1Plus,G1Minus,G2Plus,G6
        ! local variables
        real*8, dimension(6,6) :: R_GC,R_CG,D
        real*8, dimension(6) :: trialStress, stress_np,stress_np_crack
        real*8, dimension(6) :: strain_0_crack,stress_0_crack
        real*8, dimension(6) :: stress_crack
        real*8 dE1Plus,dL1Plus,dF1Plus,d1Plus,d1Minus,d1,d2,d3,d4,d5,d6
        real*8 d1State,d1StateOld,dMat,dMatState,dMatStateOld,psi
        real*8 r1Plus,rE1Plus,rF1Plus,rL1Plus,K1,A1Plus,r1Minus,A1Minus
        real*8 nu21,nu31,nu32,nu12_d,nu21_d,nu13_d,nu31_d,nu23_d,nu32_d
        real*8 etaT,kappa,lambda,a_fail,lch_fib,lch_mat,epsR_f,F
        real*8 Gamma,g0,tR_0,epsR_0,tN_0,tL_0,tT_0,epsT_0,epsN_0,epsL_0
        real*8 FI_LT,FI_LC,FI_M,FI_LT_old,FI_LC_old,FI_LT_max,FI_LC_max
        integer i
        ! output variables
        real*8, dimension(6), intent(inout) :: stress
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew
        ! external functions
        real*8, external :: linear_damage,mccauley,det_3x3

        ! Recalculate additional Poisson's ratio terms
        nu21 = nu12*(e22/e11) ! Poisson's ratio 21
        nu31 = nu13*(e33/e11) ! Poisson's ratio 31
        nu32 = nu23*(e33/e22) ! Poisson's ratio 32

        ! Parameters needed for failure criteria calculation
        etaT = (etaL*ST)/SL  !Eq. 10 [1]

        kappa = (ST**2.0d0-YT**2.0d0)/(ST*YT)  !Eq. 43 [1]
        lambda = ((2.0d0*etaL*ST)/SL)-kappa  !Eq. 45 [1] 

        ! Load old state variables
        FI_LT_old = stateOld(k,7)
        FI_LC_old = stateOld(k,8)

        ! Calculate stress using an undamaged stiffness matrix
        trialStress = matmul(C, strain)

        ! Initialize damage activation functions (failure indices)
        FI_LT = 0.0d0
        FI_LC = 0.0d0

        ! Calculate longitudinal tensile failure indices
        if (trialStress(1) >= 0.d0) then
          FI_LT = strain(1)/(XT/e11) ! Eq. 54 [1]
        else if (trialStress(1) < 0.d0) then
          FI_LC = abs(strain(1))/(XC/e11) ! Eq. 6 [4]
        end if

        ! Longitudinal tensile damage
        ! Compute max failure indices
        FI_LT_max = max(FI_LT,FI_LT_old)
        FI_LC_max = max(FI_LC,FI_LC_old)
        ! Record failure indices in state array
        stateNew(k,7) = FI_LT_max
        stateNew(k,8) = FI_LC_max   

        if ((FI_LT_max >= 1.0d0).or.(FI_LC_max >= 1.0d0)) then
          ! Check if damage has already initiated
          if (stateOld(k,9) == 0.0d0) then
            stateNew(k,9) = 1.0d0 ! longitudinal failure flag
            stateNew(k,10) = lch ! record char. length at onset
            lch_fib = lch
          else
            ! load fiber direction characteristic length
            lch_fib = stateOld(k,10)
            stateNew(k,9) = 1.0d0 ! ensure state retention
            stateNew(k,10) = lch_fib
          end if
          ! calculate linear-exponential softening response (Eq. 12 [3])
          K1 = (lch_fib*XT*e11*(XT-XPO))/(2.0d0*GL1Plus*e11-lch_fib*XT*
     1         (XT-XPO)) 
          rF1Plus = 1.0d0 + (e11/K1)*(1.0d0-XPO/XT)
          dF1Plus = 1.0d0 + (K1/e11)-(K1/e11+1.0d0)*(1.0d0/rF1Plus)
          A1Plus = (2.0d0*lch_fib*XPO**2.0d0)/(2.0d0*(1.0d0-dF1Plus)*
     1              e11*GE1Plus-lch_fib*XPO**2.0d0)
          r1Plus = max(1.0d0, FI_LT_max, FI_LC_max)
          rL1Plus = max(1.0d0, min(r1Plus, rF1Plus)) ! Eq. 8 [3]
          rE1Plus = max(1.0d0, (1.0d0-dF1Plus)*(XT/XPO)*r1Plus) 
          dL1Plus = 1.0d0 + (K1/e11) - (K1/e11 + 1.0d0)*(1.0d0/rL1Plus) ! Eq. 7 [3]
          dE1Plus = 1.0d0 - (1.0d0/rE1Plus)*exp(A1Plus*(1.0d0-rE1Plus))
          d1Plus = 1.0d0 - (1.0d0-dL1Plus)*(1.0d0-dE1Plus) ! Eq. 5 [3]
        else
          d1Plus = 0.0d0
        end if

        ! Longitudinal compressive damage
        if (FI_LC_max >= 1.0d0) then
          ! Check if damage has already initiated
          if (stateOld(k,9) == 0.0d0) then
            stateNew(k,9) = 1.0d0 ! longitudinal failure flag
            stateNew(k,10) = lch ! record char. length at onset
          else
            lch_fib = stateOld(k,10) ! load char. length
            do i = 9,10 ! ensure state variables are retained
              stateNew(k,i) = stateOld(k,i)
            end do
          end if
          r1Minus = FI_LC_max
          A1Minus = (2.0d0*lch_fib*XC**2.0d0)/
     1              (2.0d0*e11*G1Minus - lch_fib*XC**2.0d0)
          d1Minus = 1.0d0 - (1.0d0/r1Minus)*exp(A1Minus*(1.0d0-r1Minus))
        else
          d1Minus = 0.0d0
        end if
        
        ! Calculate damage variables
        if (trialStress(1) /= 0.0d0) then ! longitudinal damage
          d1 = d1Plus*(mccauley(trialStress(1))/abs(trialStress(1))) +
     1         d1Minus*(mccauley(-trialStress(1))/abs(trialStress(1)))
        else
          d1 = 0.0d0
        end if
  
        ! Matrix damage
        FI_M = 0.0d0 ! initialize variables

        ! Check if damage has already initiated
        if (stateOld(k,11) == 0.0d0) then
          ! if not, check if failure occurs in current step
          call fail_initiation(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                         FI_M,a_fail)
          if (FI_M >= 1.0d0) then
            ! if damage is triggered
            stateNew(k,11) = 1.0d0 ! matrix damage onset flag
            stateNew(k,12) = lch ! record char. length at onset
            stateNew(k,13) = a_fail ! record failure plane
            ! Rotate onset strains to failure plane
            call rot_matrix(a_fail,R_GC) ! transform global to crack
            strain_0_crack = matmul(R_GC, strain)
            ! Assign onset strains
            epsN_0 = strain_0_crack(2)
            epsT_0 = strain_0_crack(5)*2.0d0
            epsL_0 = strain_0_crack(4)*2.0d0
            ! Rotate onset stresses to failure plane
            stress_0_crack = matmul(R_GC, trialStress)
            ! Assign onset stresses
            tN_0 = stress_0_crack(2)
            tT_0 = stress_0_crack(5)
            tL_0 = stress_0_crack(4)
            ! calculuate equivalent strain at onset (Eq. 17 [4])
            epsR_0 = (epsT_0**2.0d0 + mccauley(epsN_0)**2.0d0 +
     1               epsL_0**2.0d0)**0.5d0
            ! Calculate equivalent stress at onset (Eq. 16 [4])
            tR_0 = (tT_0**2.0d0 + mccauley(tN_0)**2.0d0 +
     1             tL_0**2.0d0)**0.5d0
            ! Calculate mixed-mode critical strain energy release rate (Eq. 25 [4])
            Gamma = G6*(tT_0/tR_0)**2.0d0 +
     1              G6*(tL_0/tR_0)**2.0d0 +
     2              G2Plus*(mccauley(tN_0)/tR_0)**2.0d0
            ! calculate mixed mode ultimate failure strain (Eq. 22 [4])
            epsR_f = (2.0d0 * Gamma) / (tR_0 * lch)
            if (epsR_f <= epsR_0) then
              print *, 'Characteristic element length too large'
              print *, a_fail, tR_0, Gamma, g0, epsR_0, lch
              print *, trialStress
              call xplb_exit
            end if
            ! Record variables at onset
            stateNew(k,14) = epsR_0
            stateNew(k,15) = tR_0
            stateNew(k,16) = Gamma
            stateNew(k,17) = g0
            stateNew(k,18) = epsR_f
            call matrix_damage(strain,a_fail,epsR_0,epsR_f,dMat)
          else
            dMat = 0.0d0
            a_fail = 0.0d0
          end if
        else
          ! load onset variables
          lch_mat = stateOld(k,12)
          a_fail = stateOld(k,13) ! load failure plane
          epsR_0 = stateOld(k,14)
          tR_0 = stateOld(k,15)
          Gamma = stateOld(k,16)
          g0 = stateOld(k,17)
          epsR_f = stateOld(k,18)
          call matrix_damage(strain,a_fail,epsR_0,epsR_f,dMat)
          do i = 11,18
            stateNew(k,i) = stateOld(k,i)
          end do
        end if
  
        ! Record fiber direction damage state
        d1StateOld = stateOld(k,20)
        d1State = min(0.9d0,max(d1,d1StateOld))
        stateNew(k,20) = d1State
        ! Record matrix damage state
        dMatStateOld = stateOld(k,19)
        dMatState = min(0.9d0,max(dMat,dMatStateOld))
        stateNew(k,19) = dMatState

        ! Degrade Poisson ratios (Eq. 3 [4])
        nu12_d = nu12*(1.0d0 - d1State)
        nu13_d = nu13*(1.0d0 - d1State)
        nu21_d = nu21*(1.0d0 - dMatState)
        nu23_d = nu23*(1.0d0 - dMatState)
        nu31_d = nu31*(1.0d0 - dMatState)
        nu32_d = nu32*(1.0d0 - dMatState)

        ! Calculate degraded Poisson ratio stiffness matrix
        psi = (1.0d0 - nu12_d*nu21_d - nu23_d*nu32_d - nu31_d*nu13_d -
     1         2.0d0*nu12_d*nu23_d*nu31_d)/(e11*e22*e33)

        C = 0.0d0 ! initialize stiffness matrix, set all elements = 0
        C(1,1) = (1.0d0 - nu23_d*nu32_d)/(e22*e33*psi)
        C(1,2) = (nu21_d + nu31_d*nu23_d)/(e22*e33*psi)
        C(1,3) = (nu31_d + nu21_d*nu32_d)/(e22*e33*psi)
        C(2,1) = (nu12_d + nu13_d*nu32_d)/(e11*e33*psi)
        C(2,2) = (1.0d0 - nu31_d*nu13_d)/(e11*e33*psi)
        C(2,3) = (nu32_d + nu31_d*nu12_d)/(e11*e33*psi)
        C(3,1) = (nu13_d + nu12_d*nu23_d)/(e11*e22*psi)
        C(3,2) = (nu23_d + nu13_d*nu21_d)/(e11*e22*psi)
        C(3,3) = (1.0d0 - nu12_d*nu21_d)/(e11*e22*psi)
        C(4,4) = g12
        C(5,5) = g23
        C(6,6) = g13

        ! Calculate non-physical stresses
        stress_np = matmul(C, strain)

        ! Define rotation matrices to transform to and from crack CSYS
        call rot_matrix(a_fail,R_GC) ! transform global to crack
        call rot_matrix(-a_fail,R_CG) ! transform crack to global

        ! Rotate non-physical stresses to crack coordinate system
        stress_np_crack = matmul(R_GC, stress_np)

        ! Define damage matrix D
        D = 0.0d0
        D(1,1) = 1.0d0 - d1State
        if (stress_np_crack(2) == 0.0d0) then
          D(2,2) = 1.0d0 ! avoid division by zero
        else
          D(2,2) = 1.0d0 - dMatState*(mccauley(stress_np_crack(2))/
     1             stress_np_crack(2))
        end if
        D(3,3) = 1.0d0
        D(4,4) = 1.0d0 - dMatState
        D(5,5) = 1.0d0 - dMatState
        D(6,6) = 1.0d0

        ! Calculate degraded stresses using damage matrix (Eq. 1 [4])
        stress_crack = matmul(D, stress_np_crack)

        ! transform back to material coordinate system
        stress = matmul(R_CG, stress_crack)
        
        ! define global scalar damage variables (for tracking ONLY)
        d2 = 1.0d0 - abs(stress(2)/trialStress(2)) ! damage global 2
        d3 = 1.0d0 - abs(stress(3)/trialStress(3)) ! damage global 3
        d4 = 1.0d0 - abs(stress(4)/trialStress(4)) ! damage global 12
        d5 = 1.0d0 - abs(stress(5)/trialStress(5)) ! damage global 23
        d6 = 1.0d0 - abs(stress(6)/trialStress(6)) ! damage global 13
        ! Set as state variables
        stateNew(k,21) = d2
        stateNew(k,22) = d3
        stateNew(k,23) = d4
        stateNew(k,24) = d5
        stateNew(k,25) = d6
        ! Element deletion criteria
        ! calculate the determinant of the defomartion gradient
        F = det_3x3(defGradNew(k,:))
        stateNew(k,26) = F
        if (((F > 0.0d0).and.(F < 0.5d0)).or.
     1      (F >= 2.5d0).or.
     2      (d1State >= 0.9d0).or.
     3      (min(strain(1), strain(2)) <= -1.0d0)) then
          stateNew(k,27) = 0.0d0  ! flag element for deletion
        end if
      end subroutine cdm

      subroutine matrix_damage(strain,aR,epsR_0,epsR_f,dM)
        ! Purpose: Calculate matrix damage evolution according to
        ! the methodlogy proposed in Tan et al. [4]
        ! Variable dictionary
        ! strain = global strain vector
        ! epsR_0 = equivalent strain at damage onset
        ! tR_0 = equivalent tractions at damage onset
        ! aR = angle of failure plane
        ! R_GC = transformation matrix - global to failure CSYS
        ! strain_crack = strain in failure plane coordinate system
        ! epsN, epsT, epsL = strains in crack coordinate system
        !   corresponding to tractions tN, tT, tL
        ! epsR_f = equivalent strain at final failure
        ! d/dM = matrix scalar damage variable
        implicit none
        ! input variables
        real*8, dimension(6), intent(in) :: strain
        real*8, intent(in) :: epsR_0,aR,epsR_f
        ! local variables
        real*8, dimension(6,6) :: R_GC
        real*8, dimension(6) :: strain_crack
        real*8 epsN,epsT,epsL,epsR,d
        integer ios
        ! output variables
        real*8, intent(out) :: dM
        ! external functions
        real*8, external :: mccauley
        ! Rotate current strains to failure plane
        call rot_matrix(aR,R_GC) ! transform global to crack
        strain_crack = matmul(R_GC, strain)
        epsN = strain_crack(2)
        epsT = strain_crack(5)*2.0d0
        epsL = strain_crack(4)*2.0d0
        ! Calculate equivalent strain (Eq. 17 [4])
        epsR = (epsT**2.0d0 + mccauley(epsN)**2.0d0 + 
     1          epsL**2.0d0)**0.5d0
        ! calculate damage variable
        d = (epsR_f*(epsR-epsR_0))/(epsR*(epsR_f-epsR_0)) ! Eq. 19 [4]
        dM = max(0.0d0,d)
      end subroutine matrix_damage

      pure function mccauley(val)
        ! Purpose: implements the McCauley operator
        ! Variable dictionary:
        ! val = variable to execute mccauley brackets operation on
        implicit none
        real*8, intent(in) :: val
        real*8 :: mccauley
        if (val >= 0.0d0) then
            mccauley = val
        else
            mccauley = 0.0d0
        end if
      end function mccauley
        
      subroutine fail_initiation(trialStress,ST,SL,etaL,etaT,lambda,
     1                           kappa,FI,a_fail)
        ! Purpose: Iterates over failure plane angles from 0 to 180
        !   degrees to determine the failure angle and max failure
        !   index (using Catalanotti failure criteria)
        ! Variable Dictionary:
        ! FI = failure index
        ! a = angle of failure plane in degrees
        ! aR = angle of failure plane in radians

        implicit none
        ! input variables
        real*8, dimension(6), intent(in) :: trialStress
        real*8, intent(in) :: ST,SL,etaL,etaT,lambda,kappa
        ! local variables
        real*8 :: pi,aR,trialFI
        integer :: a
        ! output variables
        real*8, intent(out) :: FI,a_fail
    
        pi = 4.0d0*atan(1.0d0) ! determine value of pi
    
        FI = 0.0d0 ! initialize failure index
        a = 0 ! initial failure plane angle
        ! iterate over angles to determine angle which maximizes FIs
        do while (a <= 180) ! iterate over angles from 0 to 180 degrees
          aR = a*(pi/180.0d0)  ! angle in radians
          call catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                     aR,trialFI)
          ! Update max value if current value > max value
          if (trialFI > FI) then
            FI = trialFI
            a_fail = aR ! record failure plane 
          end if
          ! Update angle
          a = a + 1
        end do
      end subroutine fail_initiation

      subroutine catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                       aR,trialFI)
        ! Purpose: Failure criteria from Catalanotti et al. (2013)
        ! Variable Dictionary:
        ! trialStress = stress calculated using the undamaged
        !   stiffness tensor.
        ! trialStress_crack = trialStress in crack coordinate system
        ! etaT = friction coefficient in the transverse direction
        ! etaL = friction coefficient in the longitudinal direction
        ! kappa = parameter used to calculate failure indices
        ! lambda = parameter used to calculate failure indices
        ! ST = in-situ transverse shear strength
        ! SL = in-situ longitudinal shear strength
        ! trialFI_MT = stores trial values of tensile failure index
        ! trialFI_MC = stores trial values of compression failure index
        ! trialFI = max of tensile and compressive failure indices
        ! tN, tT, tL = tractions on failure plane
        ! aR = angle in radians
        ! R_GC = transformation matrix from global to crack (failure
        !   plane) coordinate system
        implicit none
        ! input variables
        real*8, dimension(6), intent(in) :: trialStress
        real*8, intent(in) :: ST,SL,etaL,etaT,lambda,kappa,aR
        ! local variables
        real*8 :: tN,tT,tL,trialFI_MT,trialFI_MC
        real*8, dimension(6,6) :: R_GC
        real*8, dimension(6) :: trialStress_crack
        ! output variables
        real*8, intent(out) :: trialFI
        call rot_matrix(aR, R_GC)
        ! rotate stresses to crack coordinate frame
        trialStress_crack = matmul(R_GC(1:6,1:6), trialStress(1:6))
        ! define tractions
        tN = trialStress_crack(2)
        tT = trialStress_crack(5)
        tL = trialStress_crack(4)
        ! Calculate value of failure indices at current angle
        if (tN >= 0.0d0) then
          trialFI_MT = (tN/ST)**2.0d0 + (tL/SL)**2.0d0 +
     1                 (tT/ST)**2.0d0 + lambda*(tN/ST)*(tL/SL)**2.0d0 +
     2                 kappa*(tN/ST) ! Eq. 42 [1]
          trialFI_MC = 0.0d0
        else
          trialFI_MC = (tL/(SL-etaL*tN))**2.0d0 +
     1                 (tT/(ST-etaT*tN))**2.0d0 ! Eq. 5 [1]
          trialFI_MT = 0.0d0
        end if
        ! Return the maximum trial failure index
        trialFI = max(trialFI_MT,trialFI_MC)
      end subroutine catalanotti

      subroutine rot_matrix(angle, R)
        ! Purpose: defines transformation matrix to rotate from failure
        !   plane coordinates system and back.
        ! Variable Dictionary:
        ! angle = angle through which to rotate
        ! R = transformation matrix
        ! m = cosine of angle
        ! n = sine of angle
        implicit none
        ! input variables
        real*8, intent(in) :: angle
        ! local variables
        real*8 m,n
        ! output variables
        real*8, dimension(6,6), intent(out) :: R
        m = cos(angle)
        n = sin(angle)
        R = 0.0d0
        R(1,1) = 1.0d0
        R(2,2) = m**2.0d0
        R(2,3) = n**2.0d0
        R(2,5) = 2.0d0*m*n
        R(3,2) = n**2.0d0
        R(3,3) = m**2.0d0
        R(3,5) = -2.0d0*m*n
        R(4,4) = m
        R(4,6) = n
        R(5,2) = -m*n
        R(5,3) = m*n
        R(5,5) = (m**2.0d0) - (n**2.0d0)
        R(6,4) = n
        R(6,6) = m
      end subroutine rot_matrix

      function det_3x3(A) result(det)
        ! Purpose: Compute the determinant of a 3x3 matrix (input as a
        !   9 element vector)
        ! Variable Dictionary:
        ! A = input matrix (9x1)
        ! det = determinant of A
        implicit none
        real*8, dimension(9), intent(in)  :: A
        real*8 det
        det = A(1)*(A(2)*A(3) - A(5)*A(8)) -
     1        A(4)*(A(7)*A(3) - A(5)*A(6)) +
     2        A(9)*(A(7)*A(8) - A(2)*A(6))
      end function det_3x3
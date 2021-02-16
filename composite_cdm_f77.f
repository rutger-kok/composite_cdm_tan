!----------------------------------------------------------------------!
!     ABAQUS VUMAT USER subroutine: catalanotti_FC_3D.for              !
!     Author(s): Rutger Kok, Francisca Martinez-Hergueta               !
!     Last Updated: 18/01/2021                                         !
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
! Variable Dictionary:
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
! ! = stiffness matrix (6x6)
! k = material point index
! i = integer used for indexing arrays in a number of do loops
! strain = strain tensor in Voigt notation, see Abaqus documentation
!           for more info (strain is diff. in implicit and explicit)
! stress = stress tensor in Voigt notation
! stressPower = internal energy per unit mass
!----------------------------------------------------------------------!

! User subroutine VUMAT
      include 'transverse_damage.f'
      include 'damage_evolution.f'
      include 'failure_criteria.f'

! User subroutine VUMAT
      subroutine vumat (
! Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'

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

        ! Declare variables used in main part of script
        real*8, dimension(6) :: stress,strain
        real*8, dimension(6,6) :: C
        real*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23 ! elastic const.
        real*8 nu21,nu31,nu32,delta
        real*8 XT,XC,YT,YC,SL,alpha0,stressPower
        real*8 G1plus,G1minus,G2plus,G2minus,G6 ! fracture energies
        integer k,i

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




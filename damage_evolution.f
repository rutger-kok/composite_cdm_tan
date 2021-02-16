!DIR$ FREEFORM

module damage_evolution

  ! no variables

contains 

  subroutine cdm(k,nblock,nstatev,strain,stateOld,lch,C,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,XT,XC,YT,YC,SL, &
                G1plus,G1minus,G2plus,G2minus,G6,alpha0,stress,stateNew)
    ! Purpose: Implements continuum damage mechanics framework from Maimi et al. (2007)
    ! Variable Dictionary:
    ! lch = characteristic element length, note: for a single element!
    ! trialStress = stress calculated using the undamaged stiffness tensor.
    ! trialStressP = trialStress rotated into the misalignment coordinate
    !     frame to calculate the failure index for longitudinal compression.
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
    ! etaT = friction coefficient in the transverse direction
    ! etaL = friction coefficient in the longitudinal direction
    ! phiC = initial misalignment angle for compressive failure calculation
    ! kappa = parameter used to calculate failure indices
    ! lambda = parameter used to calculate failure indices
    ! ST = in-situ (ideally) transverse shear strength
    ! minDamage = element deletion flag. If element is damaged in all
    !     directions the element is flagged for deletion.

    implicit none
    ! input variables
    integer, intent(in) :: nblock,nstatev,k
    real*8, dimension(nblock,nstatev), intent(in) :: stateOld
    real*8, dimension(6,6), intent(inout) :: C
    real*8, dimension(6), intent(in) :: strain
    real*8, intent(in) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23
    real*8, intent(in) :: XT,XC,YT,YC,SL,lch
    real*8, intent(in) :: G1plus,G1minus,G2plus,G2minus,G6
    ! local variables
    real*8, dimension(6) :: trialStress, trialStressP, e0
    real*8, dimension(3) :: tbar_cr
    real*8 d1Plus,d1Minus,d2Plus,d2Minus,d3Plus,d3Minus
    real*8 d1,d2,d3,d6
    real*8 nu21,nu31,nu32,delta,minDamage
    real*8 alpha0,etaT,etaL,phiC,kappa,lambda,ST,a,tN,tT,tL,crack_initiation
    real*8 FI_LT,FI_LC,FI_MT,FI_MC
    real*8 dState1,dState2,dState3,dState6
    real*8 dState1Old,dState2Old,dState3Old,dState6Old
    integer :: i,j
    ! output variables
    real*8, dimension(6), intent(inout) :: stress
    real*8, dimension(nblock,nstatev), intent(inout) :: stateNew

    ! Calculate parameters needed for failure criteria calculation
    etaL = 0.280622004043348d0
    phiC = atan((1.0d0-sqrt(1.0d0-4.0d0*(SL/XC)*((SL/XC)+etaL)))/(2.0d0*((SL/XC)+etaL)))  !Eq.12 Maimi
    ST = (0.5d0*(((2.0d0*sin(alpha0)**2.0d0)-1.0d0)*SL)/(((1.0d0-sin(alpha0)**2.0d0)**0.5d0)*sin(alpha0)*etaL)) !Eq.12 CLN
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
      call catalanotti(trialStressP,ST,SL,etaL,etaT,lambda,kappa,FI_LT,FI_LC,a,tN,tT,tL)
    end if

    ! Calculate scalar damage variables assuming triangular
    ! damage dissipation.
    if (FI_LT > 1.0d0) then ! Longitudinal tensile damage variable
      d1Plus= triangular(XT,e11,G1Plus,strain(1),lch)
    else
      d1Plus = 0.0d0
    end if

    if (FI_LC > 1.0d0) then ! Longitudinal compressive damage variable
        d1Minus = triangular(XC,e11,G1minus,strain(1),lch)
    else
        d1Minus = 0.0d0
    end if

    ! Calculate damage variables
    if (trialStress(1) /= 0.0d0) then ! Longitudinal damage variable
      d1 = d1Plus*(mccauley(trialStress(1))/abs(trialStress(1))) + d1Minus*(mccauley(-trialStress(1))/abs(trialStress(1)))
    else
      d1 = 0.0d0
    end if

    ! Calculate damage state variables. Damage state must be
    ! retained between iterations.
    dState1 = min(0.999d0,max(d1,dState1Old))

    ! Calculate damaged stiffness tensor
    ! Recalculate additional Poisson's ratio terms (not passed as args.)
    nu21 = nu12*(e22/e11)
    nu31 = nu13*(e33/e11)
    nu32 = nu23*(e33/e22)
  
    delta = 1.0d0/(e11*e22 - e22**2.0d0*nu12**2.0d0 + dState1*e22**2.0d0*nu12**2.0d0 - e11*e33*nu23**2.0d0 - e22*e33*nu13**2.0d0 + &
            dState1*e22*e33*nu13**2.0d0 - 2.0d0*e22*e33*nu12*nu13*nu23 + 2.0d0*dState1*e22*e33*nu12*nu13*nu23)

    C = 0.0d0 ! set all elements in array equal to zero
    C(1,1) = -((- e33*e11**2.0d0*nu23**2.0d0 + e22*e11**2.0d0)*(dState1 - 1.0d0))*delta
    C(2,1) = -((e11*nu12*e22**2.0d0 + e11*e33*nu13*nu23*e22)*(dState1 - 1.0d0))*delta
    C(3,1) = -((dState1 - 1.0d0)*(e11*e22*e33*nu13 + e11*e22*e33*nu12*nu23))*delta
    C(1,2) = (e22*(e11*e22*nu12 + e11*e33*nu13*nu23 - dState1*e11*e22*nu12 - dState1*e11*e33*nu13*nu23))*delta
    C(2,2) = (e22*(e11*e22 - e22*e33*nu13**2.0d0 + dState1*e22*e33*nu13**2.0d0))*delta
    C(3,2) = (e22*(e11*e33*nu23 + e22*e33*nu12*nu13 - dState1*e22*e33*nu12*nu13))*delta
    C(1,3) = (e22*e33*(e11*nu13 - dState1*e11*nu13 + e11*nu12*nu23 - dState1*e11*nu12*nu23))*delta
    C(2,3) = (e22*e33*(e11*nu23 + e22*nu12*nu13 - dState1*e22*nu12*nu13))*delta
    C(3,3) = (e22*e33*(e11 - e22*nu12**2.0d0 + dState1*e22*nu12**2.0d0))*delta
    C(4,4) = g12
    C(5,5) = g23
    C(6,6) = g13

    ! Calculate transverse failure indices using fail_cln function
    call catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,FI_MT,FI_MC,a,tN,tT,tL)

    crack_initiation = stateOld(k,18)
    if ((FI_MT > 1.0d0).or.(FI_MC > 1.0d0)) then ! transverse damage
      if (crack_initiation == 0.0d0) then
        stateNew(k,18) = 1.0d0
        do i = 7,12
          stateNew(k,i) = strain(i)
        end do
        stateNew(k,15) = tL
        stateNew(k,16) = tN
        stateNew(k,17) = tT
      end if

      do i = 7,12
        j = i - 6
        e0(j) = stateNew(k,i)
      end do

      tbar_cr(1) = stateNew(k,15)
      tbar_cr(2) = stateNew(k,16)
      tbar_cr(3) = stateNew(k,17)

      call smeared_crack(strain, e0, a, tbar_cr, C, G1Plus, G2Plus, lch, e22, stress, dState2)

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
    if (minDamage > 0.99) stateNew(k,20) = 0.d0
  end subroutine cdm

  subroutine rotate_stress(trialStress,phiC,XC,g12,trialStressP)
    ! Purpose: Rotate stresses to misaligned coordinate for calculation of longitudinal compressive failure index.
    ! Variable dictionary
    ! trialStress = stress calculated using the undamaged stiffness tensor.
    ! trialStressP = trialStress rotated into the misalignment coordinate frame to calculate the failure index for longitudinal
    !                compression.
    ! trialStressT = trialStress rotated by in-plane angle theta.
    ! phiC = initial misalignment angle for compressive failure calculation
    ! g12 = longitudinal compressive failure stress
    ! phi = angle of kink band
    ! theta = kinking plane angle
    ! phi0 = initial misalignment angle (manufacturing defects etc.)
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
        theta = 0.5d0*atan((2.0d0*trialStress(5))/(trialStress(2)-trialStress(3))) !Eq.55 CLN
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
    trialStressT(2) = trialStress(2)*m**2 + 2.0d0*trialStress(5)*m*n + trialStress(3)*n**2
    trialStressT(3) = trialStress(3)*m**2 - 2.0d0*trialStress(5)*m*n + trialStress(2)*n**2
    trialStressT(4) = trialStress(4)*m + trialStress(6)*n
    trialStressT(5) = trialStress(5)*(m**2 - n**2) - trialStress(2)*n*m + trialStress(3)*n*m
    trialStressT(6) = trialStress(6)*m - trialStress(4)*n

    ! Determine kink band angle phi
    gammaMC = (sin(2*phiC)*XC)/(2*g12)  ! eq 74
    phi0 = phiC - gammaMC  ! eq 75
    gammaM = ((phi0*g12 + abs(trialStressT(4)))/(g12+trialStressT(1)-trialStressT(2)) - phi0)  ! eq 81

    if (trialStress(4) >= 0.d0) then ! eq 77
        phi = phi0 + gammaM
    else
        phi = -1.0*(phi0+gammaM)
    end if

    ! Rotate stresses by angle phi
    u = cos(phi)
    v = sin(phi)
    trialStressP(1) = trialStressT(1)*u**2 + 2.0d0*trialStressT(4)*u*v + trialStressT(2)*v**2
    trialStressP(2) = trialStressT(2)*u**2 - 2.0d0*trialStressT(4)*v*u + trialStressT(1)*v**2     
    trialStressP(3) = trialStressT(3)      
    trialStressP(4) = trialStressT(4)*(u**2 -v**2) + trialStressT(2)*v*u - trialStressT(1)*v*u
    trialStressP(5) = trialStressT(5)*u - trialStressT(6)*v     
    trialStressP(6) = trialStressT(6)*u + trialStressT(5)*v

  end subroutine rotate_stress

  pure function triangular(X,E,G,tStrain,lch) result(dam)
    ! Purpose: Implement triangular damage dissipation
    ! Variable dictionary:
    ! X = failure stress (loading and material direction unspecified)
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

    ! Calculate scalar damage variable assuming triangular dissipation
    if (tStrain == 0.0d0) then ! no strain = no damage (avoid div. by 0)
      dam = 0.0d0 
    else
      delta_0 = X/E
      delta_c = (2.0d0*G)/(X*lch)
      if (delta_c < delta_0) then  !Needed for d.nq.0-0
        delta_c = 1.1d0*delta_0   !Adjustable for stability
      end if
      dam = (delta_c*(abs(tStrain)-delta_0))/(abs(tStrain)*(delta_c-delta_0))
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
end module damage_evolution
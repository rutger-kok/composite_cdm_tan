!DIR$ FREEFORM

module failure_criteria

  ! no variables

contains

  subroutine catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,FI_MT,FI_MC,a,tN,tT,tL)
    ! Purpose: Implements transverse (and longitudinal compression) failure
    !          indices from Catalanotti et al. (2013)
    ! Variable Dictionary:
    ! trialStress = stress calculated using the undamaged stiffness tensor.
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
    real*8 trialFI_MT,trialFI_MC,aFail_MC,aFail_MT,pi,tN,tT,tL,aR
    integer a
    ! output variables
    real*8, intent(out) :: FI_MT, FI_MC
    
    pi = 4.0d0*atan(1.0d0) ! determine value of pi
    
    FI_MT = 0.0d0 ! initialize failure indices
    FI_MC = 0.0d0
    a = 0 ! initial failure plane angle

    ! iterate over angles to determine angle which maximizes FIs
    do while (a <= 180) ! iterate over angles from 0 to 180 degrees
      aR= a*(pi/180.0d0)  ! angle in radians
      ! Calculate tractions on failure plane: Eq 3 CLN (and 59-61)
      tN = trialStress(2)*cos(aR)**2.0d0 + 2.0d0*trialStress(5)*sin(aR)*cos(aR) + trialStress(3)*sin(aR)**2.0d0
      tT = -1.0d0*cos(aR)*sin(aR)*(trialStress(2)-trialStress(3))+(trialStress(5)*(cos(aR)**2.0d0 - sin(aR)**2.0d0))
      tL = trialStress(4)*cos(aR) + trialStress(6)*sin(aR)
      
      ! Calculate value of failure indices at current angle
      if (tN >= 0.0d0) then
        trialFI_MT = (tN/ST)**2.0d0 + (tL/SL)**2.0d0 + (tT/ST)**2.0d0 + lambda*(tN/ST)*(tL/SL)**2.0d0 + kappa*(tN/ST) ! Eq. 42 CLN
        trialFI_MC = 0.0d0
      else
        trialFI_MC = (tL/(SL-etaL*tN))**2.0d0 + (tT/(ST-etaT*tN))**2.0d0 ! Eq. 5 CLN
        trialFI_MT = 0.0d0
      end if

      ! Update max values if current value > max value
      if (trialFI_MT > FI_MT) then
        FI_MT = trialFI_MT
        aFail_MT = aR ! record failure plane 
      end if
      if (trialFI_MC > FI_MC) then
        FI_MC = trialFI_MC
        aFail_MC = aR ! record failure plane 
      end if
      
      ! Update angle
      a = a + 1
    end do
  end subroutine catalanotti
end module failure_criteria

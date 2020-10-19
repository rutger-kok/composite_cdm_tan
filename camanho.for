!----------------------------------------------------------------------!
!     ABAQUS VUMAT USER SUBROUTINE: catalanotti_FC_v3.for              !
!     Author(s): Rutger Kok, Francisca Martinez-Hergueta               !
!     Date: 31/10/2019                                                 !
!     Version 1.0                                                      !
!----------------------------------------------------------------------!

! Back calculate stresses in undulation regions using the Chou
! homogenization, then implement failure criteria on transformed stresses
! transverse fracture predicted by Camanho model while longitudinal CDM
! an adaptation of Maimi et al.

      SUBROUTINE camanho(strain_voigt, alpha, C_voigt, init_strain_voigt)
      ! input variables
      REAL*8, DIMENSION(6), INTENT(IN) :: strain_voigt, init_strain_voigt
      REAL*8, DIMENSION(6,6) :: C_voigt
      REAL*8, INTENT(IN) :: alpha
      ! local variables
      REAL*8, DIMENSION(3,3) :: trialStrain,totalStrain,initStrain
      REAL*8, DIMENSION(3,3,3,3) :: C
      
      totalStrain(1,1) = strain_voigt(1)
      totalStrain(2,2) = strain_voigt(2)
      totalStrain(3,3) = strain_voigt(3)
      totalStrain(2,3) = strain_voigt(4)
      totalStrain(1,3) = strain_voigt(5)
      totalStrain(1,2) = strain_voigt(6)

      initStrain(1,1) = init_strain_voigt(1)
      initStrain(2,2) = init_strain_voigt(2)
      initStrain(3,3) = init_strain_voigt(3)
      initStrain(2,3) = init_strain_voigt(4)
      initStrain(1,3) = init_strain_voigt(5)
      initStrain(1,2) = init_strain_voigt(6)

      n2 = (/1.0, cos(alpha), sin(alpha)/)

      ! define rotation matrices
      R = ROTX(alpha0)
      RT = TRANSPOSE(R)
      ! rotate total strain tensor to crack coordinate system
      totalStrain_cr = dot22(dot22(R, totalStrain), RT)
      ! rotate strain at crack initiation to crack coord system
      initStrain_cr = dot22(dot22(R, initStrain), RT)

      ! calculate initial trial strain
      trialStrain_cr = rotTotalStrain - initStrain

      DO WHILE (error.gt.tol)

        ! calculate first term of residual
  
        ! convert voigt stiffness matrix to 4th order tensor
        C = voigtToTensor(C_voigt)
        ! compute elastic stress
        elasticStress = ddot42(C, totalStrain)
        ! rotate trial strain to global coord system
        trialStrain = dot22(dot22(R, trialStrain_cr), RT)
        ! compute cracking stress 
        crackingStress = ddot42(C, trialStrain_cr)
        ! subtract cracking stress from elastic stress
        sub = elasticStress - crackingStress
        ! compute f term
        f = dot21(dot22(RT, sub), n2)

        ! define crack coordinate system normals
        n1_cr = (/1.0, 0.0, 0.0/)
        n1_cr_tp = TRANSPOSE(n1_cr)
        n2_cr = (/0.0, 1.0, 0.0/)
        n2_cr_tp = TRANSPOSE(n2_cr)
        n3_cr = (/0.0, 0.0, 1.0/)
        n3_cr_tp = TRANPOSE(n3_cr)

        ! estimate displacement jump at crack
        ! note the use of matmul(a,b) in this case is the dyadic product of
        ! the first order tensors n1,n2,n3
        w_cr1 = (2*ddot22(trialStrain_cr, MATMUL(n1_cr_tp, n2_cr)))*n1_cr
        w_cr2 = (ddot22(trialStrain_cr, MATMUL(n2_cr_tp, n2_cr)))*n2_cr
        w_cr3 = (2*ddot22(trialStrain_cr, MATMUL(n2_cr_tp, n3_cr)))*n3_cr
        w_cr = (w_cr1+w_cr2+w_cr3)*lch

        ! calculate second term of the residual according to the
        ! cohesive law

        ! compute lambda 
        lambda = (w_cr(1)**2 + mccauley(w_cr(2))**2 + w_cr(3)**2)**0.5

        ! compute tbar_cr, the tractions on the failure plane at damage onset
        tbar_cr_s = (tbar_cr(1)**2 + tbar_cr(2)**2)**0.5
        ! compute the shear components of the displacement jump
        w_cr_s = (w_cr(1)**2 + w_cr(2)**2)**0.5   
        beta = w_cr_s/(w_cr_s+w_cr(2))
        ! compute Benzeggagh-Kenane parameters A and B
        A = GIIc - GIc
        B = (tbar_cr_s*beta)/(beta*(tbar_cr_s-tbar_cr(2))+tbar_cr(2))
        ! compute the 2 norm of the tbar vector
        tbar_cr_norm = (tbar_cr(1)**2 + tbar_cr(2)**2
      1               + tbar_cr(3)**2)**0.5
        wmf = (2*(GIc + A*B**eta))/tbar_cr

        ! compute the damage parameter d
        d = max(0.0d0, min((lambda/wmf),1.0d0))

        ! compute the tractions on the fracture plane
        t_cr(i) = ((1-d)/d)*(w_cr(i)/wmf)*tbar_cr(i)-kdelta(i,2)*
      1           (mccauley(-w_cr(2))/w_cr(2))*(((1-d)/d)*(w_cr(i)/wmf)*
      2           tbar_cr(i) - e22*(totalStrain_cr(2,2)-
      3           trialStrain_cr(2,2)))

        g = t_cr
        ! calculate the residual
        r = f - g

        ! calculate the Jacobian matrix A

        ! solve for strainDelta_cr



      END SUBROUTINE camanho

!----------------------------------------------------------------------!
! Function ROTX: rotates 2nd order tensor by angle about the xx axis   !
!----------------------------------------------------------------------!

      FUNCTION ROTX(angle)
        IMPLICIT NONE
        ! input variables
        REAL*8, INTENT(IN) :: angle
        ! local variables
        REAL*8 m,n
        REAL*8, PARAMETER :: eps=1.d-8
        ! output variables
        REAL*8, DIMENSION(3,3) :: ROTX

        m = cos(theta)
        n = sin(theta)
        ROTX = 0.0d0
        ! Rotate stresses by angle theta
        ROTX(1,1) = 1.0d0
        ROTX(2,2) = m
        ROTX(2,3) = n
        ROTX(3,2) = -n
        ROTX(3,3) = m

        RETURN
      END FUNCTION ROTX


!----------------------------------------------------------------------!
! Function mccauley: Implements the McCauley operator                  !
!                                                                      !
!----------------------------------------------------------------------!

      REAL*8 FUNCTION mccauley(value)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: value 
        mccauley = (value+abs(value))/2.0d0
        RETURN
      END FUNCTION mccauley

!----------------------------------------------------------------------!
! Function ddot42: Double contraction of 4th and 2nd order tensors     !
!                                                                      !
!----------------------------------------------------------------------!

      FUNCTION ddot42(A,B)
        IMPLICIT NONE
        ! input variables
        REAL*8, DIMENSION(3,3,3,3), INTENT(IN) :: A
        REAL*8, DIMENSION(3,3), INTENT(IN) :: B
        ! local variables
        INTEGER :: i,j,k,l
        ! output variables
        REAL*8, DIMENSION(3,3) :: ddot42

        ddot42 = 0.0d0
        DO i = 1,3
          DO j = 1,3
            DO k = 1,3
              DO l = 1,3
                ddot42(i,j) = ddot42(i,j) + A(i,j,k,l)*B(k,l)
              END DO
            END DO
          END DO
        END DO
      
        RETURN
      END FUNCTION ddot42

!----------------------------------------------------------------------!
! Function ddot22: Double contraction of two 2nd order tensors         !
!                                                                      !
!----------------------------------------------------------------------!

      FUNCTION ddot22(A,B)
        IMPLICIT NONE
        ! input variables
        REAL*8, DIMENSION(3,3), INTENT(IN) :: A
        REAL*8, DIMENSION(3,3), INTENT(IN) :: B
        ! local variables
        INTEGER :: i,j
        ! output variables
        REAL*8 :: ddot22  ! scalar

        ddot22 = 0.0d0
        DO i = 1,3
          DO j = 1,3
            ddot22 = ddot22 + A(i,j)*B(i,j)
          END DO
        END DO

        RETURN
      END FUNCTION ddot22

!----------------------------------------------------------------------!
! Function dot22: Single contraction of two 2nd order tensors         !
!                                                                      !
!----------------------------------------------------------------------!

      FUNCTION dot22(A,B)
        IMPLICIT NONE
        ! input variables
        REAL*8, DIMENSION(3,3), INTENT(IN) :: A,B
        ! local variables
        INTEGER :: i,j,k
        ! output variables
        REAL*8, DIMENSION(3,3) :: dot22 ! result is a 2nd order tensor

        dot22 = 0.0d0
        DO i = 1,3
          DO j = 1,3
            DO k = 1,3
                dot22(i,j) = dot22(i,j) + A(i,k)*B(k,j)
            END DO
          END DO
        END DO
      
        RETURN
      END FUNCTION dot22

!----------------------------------------------------------------------!
! Function dot21: Single contraction of a 2nd and 1st order tensor     !
!                                                                      !
!----------------------------------------------------------------------!

      FUNCTION dot21(A,B)
        IMPLICIT NONE
        ! input variables
        REAL*8, DIMENSION(3,3), INTENT(IN) :: A
        REAL*8, DIMENSION(3,1), INTENT(IN) :: B
        ! local variables
        INTEGER :: i,k
        ! output variables
        REAL*8, DIMENSION(3) :: dot21

        dot21 = 0.0d0
        DO i = 1,3
          DO k = 1,3
              dot21(i) = dot21(i) + A(i,k)*B(k)
          END DO
        END DO
      
        RETURN
      END FUNCTION dot21

!----------------------------------------------------------------------!
! Function voigtToTensor: Converts a Voigt notation 6x6 (stiffness)    !
! matrix to a fourth order stiffness tensor                            !
!                                                                      !
!----------------------------------------------------------------------!

      FUNCTION voigtToTensor(A)
        IMPLICIT NONE
        ! input variables 
        REAL*8, DIMENSION(6,6) :: A
        ! local variables
        INTEGER :: i,j,k,l,s1,s2,idx1,idx2
        REAL*8 :: term
        ! output variables
        REAL*8, DIMENSION(3,3,3,3) :: voigtToTensor

        DO i = 1,3
          DO j = 1,3
            DO k = 1,3
              DO l = 1,3
                IF (i.eq.j) THEN
                  IF (i.eq.1) THEN
                    idx1 = 1
                  ELSE IF (i.eq.2) THEN
                    idx1 = 2
                  ELSE IF (i.eq.3) THEN
                    idx1 = 3
                  END IF
                ELSE
                    s1 = i + j
                  IF (s1.eq.5) THEN
                    idx1 = 4
                  ELSE IF (s1.eq.4) THEN
                    idx1 = 5
                  ELSE IF (s1.eq.3) THEN
                    idx1 = 6
                  END IF
                END IF
                IF (k.eq.l) THEN
                  IF (k.eq.1) THEN
                    idx2 = 1
                  ELSE IF (k.eq.2) THEN
                    idx2 = 2
                  ELSE IF (k.eq.3) THEN
                    idx2 = 3
                  END IF
                ELSE
                    s2 = k + l
                  IF (s2.eq.5) THEN
                    idx2 = 4
                  ELSE IF (s2.eq.4) THEN
                    idx2 = 5
                  ELSE IF (s2.eq.3) THEN
                    idx2 = 6
                  END IF
                END IF

                term = A(idx1,idx2)
                voigtToTensor(i,j,k,l) = term
              END DO
            END DO
          END DO
        END DO
        RETURN
      END FUNCTION voigtToTensor



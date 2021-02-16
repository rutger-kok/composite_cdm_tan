!DIR$ FREEFORM

module transverse_damage

  ! no variables
  
contains

  subroutine smeared_crack(eTotal, e0, alpha, tbar_cr, C, GIc, GIIc, lch, E22, stress, d)
    implicit none
    ! input variables
    real*8, dimension(6,6), intent(in) :: C
    real*8, dimension(6), intent(in) :: eTotal, e0
    real*8, dimension(3), intent(in) :: tbar_cr
    real*8, intent(in) :: alpha, GIc, GIIc, lch, E22
    ! local variables
    real*8, dimension(6,6) :: R_gc, R_cg
    real*8, dimension(3,6) :: H
    real*8, dimension(6) :: ec_cr, eTotal_cr, ec_cr_new, ec, e0_cr, sc, sTotal, s
    real*8, dimension(3,3) :: M, Minv, R, RT
    real*8, dimension(3) :: n1, n2, n3, w_cr, f, res, ec_cr_delta, t_cr, t
    real*8 :: A, B, diff, tol, beta, lambda, wmf, d, tbar_cr_norm, w_cr_s, tbar_cr_s, p, q, eta
    integer i
    ! output variables
    real*8, dimension(6) :: stress

    ! crack plane vectors (n2 is normal to crack plane)
    p = cos(alpha)
    q = sin(alpha)
    n1 = (/1.0d0, 0.0d0, 0.0d0 /)
    n2 = (/0.0d0, p, q /)
    n3 = (/0.0d0, -q, p /)

    ! tbar_cr = tractions on the failure plane at damage onset
    tbar_cr_s = (tbar_cr(1)**2 + tbar_cr(2)**2)**0.5

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
    diff = 1.0d8 ! initial value for difference between ec_cr_new and ec_cr
    tol = 1.0d-6 ! convergence criterion
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
      lambda = (w_cr(1)**2.0d0 + mccauley(w_cr(2))**2.0d0 + w_cr(3)**2.0d0)**0.5d0
      ! compute the shear components of the displacement jump
      w_cr_s = (w_cr(1)**2.0d0 + w_cr(2)**2.0d0)**0.5d0  
      beta = w_cr_s/(w_cr_s+w_cr(2))
      ! compute Benzeggagh-Kenane parameters A and B
      A = GIIc - GIc
      B = (tbar_cr_s*beta)/(beta*(tbar_cr_s-tbar_cr(2))+tbar_cr(2))
      ! compute the 2 norm of the tbar vector
      tbar_cr_norm = (tbar_cr(1)**2 + tbar_cr(2)**2 + tbar_cr(3)**2.0d0)**0.5d0
      wmf = (2.0d0*(GIc + A*B**eta))/tbar_cr_norm

      ! compute the damage parameter d
      d = max(0.0d0,min((lambda/wmf),1.0d0))

      ! compute the tractions on the fracture plane
      t_cr(i) = ((1.0d0-d)/d)*(w_cr(i)/wmf)*tbar_cr(i)-kdelta(i,2)*(mccauley(-w_cr(2))/w_cr(2))*(((1.0d0-d)/d)*(w_cr(i)/wmf)* &
                tbar_cr(i)- E22*(eTotal_cr(4)- ec_cr(4)))

      ! calculate the residual res 
      res = f - t_cr
      ! calculate the Jacobian matrix M
      M = jacobian(C, ec_cr, tbar_cr, lch, alpha, eta, GIc, GIIc, E22, eTotal_cr)
      ! invert M and solve for ec_cr_delta
      Minv = matinv3(M)
      ec_cr_delta = matmul(-Minv(1:3,1:3), res(1:3))
      ! calculate new crack strain values and determine error
      ec_cr_new(4) = ec_cr(4) + ec_cr_delta(1)
      ec_cr_new(2) = ec_cr(2) + ec_cr_delta(2)
      ec_cr_new(5) = ec_cr(5) + ec_cr_delta(3)
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
    real*8, intent(in) :: A(3,3)   ! Matrix
    real*8             :: B(3,3)   ! Inverse matrix
    real*8             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

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

  pure function mccauley(value)
    implicit none
    real*8, intent(in) :: value
    real*8 :: mccauley
    mccauley = (value+abs(value))/2.0d0
  end function mccauley

  pure function jacobian(C, ec_cr, tbar_cr, lch, alpha, eta, GIc, GIIc, E22, eTotal_cr) result(A)
    implicit none
    ! input variables
    real*8, dimension(6,6), intent(in) :: C
    real*8, dimension(6), intent(in) :: ec_cr, eTotal_cr
    real*8, dimension(3), intent(in) :: tbar_cr
    real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
    ! local variables
    real*8 :: m,n, term1, term2, term3, term4, term5, term6, term7
    ! output variables
    real*8, dimension(3,3) :: A
    ! define cosine and sine of fracture plane angle for brevity
    m = cos(alpha)
    n = sin(alpha)
    ! terms1-7 are repeated parts of the Jacobian matrix. They have been defined here to reduce the length of the equations.
    term1 = (GIc + ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)* &
            (GIIc - GIc))/(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)
    term2 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0* &
            (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*lch)**2.0d0)
    term3 = ((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*(tbar_cr(1)**2.0d0* &
            sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0* &
            sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0* &
            sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,ec_cr(2)*lch)* &
            sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1*term2**0.5d0) - &
            (ec_cr(4)**2.0d0*eta*tbar_cr(2)*(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)* &
            (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)* &
            (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*(tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0* &
            sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + &
            tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*term2**(1.0d0/2.0d0)* &
            (2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0* &
            (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 + &
            tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - ec_cr(2)**2.0d0* &
            tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + &
            2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/ &
            (sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0* &
            tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0))
    term4 = ((0.25d0*(tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0* &
            sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0* &
            sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*term2**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))* &
            sign(1.0d0,tbar_cr(3))*term1) - 1.0d0)
    term5 = ((4.0d0*ec_cr(4)*lch**2.0d0*(abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 + &
            abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0))/ (term1*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + &
            2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0)**0.5d0) + &
            (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)* &
            (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + &
            tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)* &
            (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0)* &
            (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + &
            abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0* &
            tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - &
            ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0* &
            tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/(term1**2.0d0* &
            (4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0* &
            tbar_cr(1)**2.0d0)**2.0d0))
    term6 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + &
            16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0)
    term7 = ((0.25d0*(abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0)* &
            (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + &
            abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/2.0d0))/term1 - 1.0d0)
    A(1,1) = C(6,6)*(m**2.0d0 - 1.0d0) - C(4,4)*m**2.0d0 + (4.0d0*lch*tbar_cr(1)*term7)/term6**0.5d0 - &
            (64.0d0*ec_cr(4)**2.0d0*lch**3.0d0*tbar_cr(1)*term7)/term6**1.5d0 + (4.0d0*ec_cr(4)*lch*tbar_cr(1)*term5)/ &
            term6**0.5d0
    A(1,2) = (ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*(sign(1.0d0,ec_cr(2)*lch) + &
            1.0d0)**2.0d0*(tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + &
            tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + &
            tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/ &
            (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1*term2) - &
            (4.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/ &
            (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0) - (4.0d0*ec_cr(4)**3.0d0*eta*lch*tbar_cr(1)*tbar_cr(2)*(ec_cr(2)*tbar_cr(2) + &
            (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/ &
            2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - &
            GIc)*(tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + &
            tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + &
            tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/ &
            2.0d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + &
            tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - &
            4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/ &
            2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + &
            2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/ &
            (sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1**2.0d0* &
            (4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + &
            ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)
    A(1,3) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*tbar_cr(1))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + &
            16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + &
            2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0
    A(2,1) = (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term5)/term6**0.5d0 - ((0.5d0*abs(ec_cr(2)*lch) - &
            0.5d0*ec_cr(2)*lch)*((2.0d0*ec_cr(2)*lch*tbar_cr(2)*term5)/term6**0.5d0 - &
            (32*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*term7)/term6**1.5d0))/(ec_cr(2)*lch) - &
            (32*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*term7)/term6**1.5d0
    A(2,2) = n**2.0d0*(C(2,3)*m**2.0d0 + 2.0d0*C(5,5)*m**2.0d0 + C(3,3)*n**2.0d0) - m**2.0d0*(C(2,2) - C(2,2)*n**2.0d0 + &
            C(2,3)*n**2.0d0 + 2.0d0*C(5,5)*n**2.0d0) + (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) - &
            eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*m**2.0d0 + 2.0d0*eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term4)/ &
            term2**0.5d0))/ec_cr(2) - (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22 - (2.0d0*lch*tbar_cr(2)*term4)/term2**0.5d0 - &
            (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term3)/term2**0.5d0 + &
            (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/ &
            (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)))/sign(1.0d0,ec_cr(2)*lch) + (2.0d0*lch*tbar_cr(2)*term4)/term2**0.5d0 + &
            (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) - eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*m**2.0d0 + &
            2.0d0*eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term4)/term2**0.5d0))/(ec_cr(2)*sign(1.0d0,ec_cr(2)*lch)) + &
            (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term3)/term2**0.5d0 - &
            (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/ &
            (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
    A(2,3) = 0.5d0*C(3,3)*2.0d0*m*n - 0.5d0*C(2,3)*2.0d0*m*n + 0.5d0*C(2,2)*2.0d0*m*n*m**2.0d0 - &
            0.5d0*C(3,3)*2.0d0*m*n*m**2.0d0 + (48*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(2))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + &
            16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + &
            2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0 - (16.0d0*ec_cr(5)*lch**2.0d0*tbar_cr(2)*abs(ec_cr(2))*abs(lch))/ &
            (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + &
            16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0
    A(3,1) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*tbar_cr(3))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + &
            ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + &
            2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**1.5d0 + &
            (4.0d0*ec_cr(4)*ec_cr(2)*ec_cr(5)*eta*lch*tbar_cr(2)*tbar_cr(3)*(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 + &
            tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - &
            eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 + &
            abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 + &
            tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - &
            ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + &
            2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/ &
            (term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + &
            ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)
    A(3,2) = (4.0d0*ec_cr(5)*lch*tbar_cr(3)*term3)/term2**0.5d0 - m*n*(C(2,2) - C(2,2)*n**2.0d0 + C(2,3)*n**2.0d0 + &
            2.0d0*C(5,5)*n**2.0d0) - m*n*(C(2,3)*m**2.0d0 + 2.0d0*C(5,5)*m**2.0d0 + C(3,3)*n**2.0d0) - &
            (4.0d0*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(3)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/ &
            (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
    A(3,3) = (64.0d0*ec_cr(5)**2.0d0*lch**3.0d0*tbar_cr(3))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + &
            16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + &
            2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0 - C(5,5) - (4.0d0*lch*tbar_cr(3))/ &
            (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + &
            16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**0.5d0 - C(2,2)*n**2.0d0*(n**2.0d0 - &
            1.0d0) + C(3,3)*n**2.0d0*(n**2.0d0 - 1.0d0) + 2.0d0*C(5,5)*n**2.0d0 + (lch*tbar_cr(3)*(abs(tbar_cr(1))**2.0d0 + &
            abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0))/(GIc + (GIIc*(tbar_cr(1)**2.0d0 + &
            tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta))/(ec_cr(2)*tbar_cr(2) + &
            (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/ &
            2.0d0))**eta - (GIc*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 + &
            ec_cr(2)**2.0d0)**(0.5d0*eta))/(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/ &
            2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)
  end function jacobian
end module transverse_damage
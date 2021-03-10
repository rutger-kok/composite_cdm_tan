      module jacobian_f77_noterms

        ! no variables

      contains
      
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
        A(1.0d0,1.0d0) = C(6,6)*(m**2.0d0 - 1.0d0) - C(4,4)*m**2.0d0 +
     1    (4.0d0*lch*tbar_cr(1)*((0.25d0*(abs(tbar_cr(1))**2.0d0 +
     2    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     3    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     4    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     5    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     6    2.0d0))/(GIc + ((tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     9    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     1    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     2    2.0d0))**eta) - 1.0d0))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     3    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     4    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     5    lch)**2.0d0)**0.5d0 - (64.0d0*ec_cr(4)**2.0d0*lch**3.0d0*
     6    tbar_cr(1)*((0.25d0*(abs(tbar_cr(1))**2.0d0 +
     7    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     8    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     9    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     1    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     2    2.0d0))/(GIc + ((tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     5    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     6    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     7    2.0d0))**eta) - 1.0d0))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     8    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     9    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     1    lch)**2.0d0)**1.5d0 + (4.0d0*ec_cr(4)*lch*tbar_cr(1)*((4.0d0*
     2    ec_cr(4)*lch**2.0d0*(abs(tbar_cr(1))**2.0d0 +
     3    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     4    2.0d0))/((GIc + ((tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     7    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     9    2.0d0))**eta)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     1    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     2    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     3    lch)**2.0d0)**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     4    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*
     7    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     8    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     9    1.5d0)*(GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
     1    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     2    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     3    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     4    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     5    2.0d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 -
     6    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
     8    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     2    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     3    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     4    tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + ((tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     7    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     9    2.0d0))**eta)**2.0d0*(4.0d0*ec_cr(4)**2.0d0*
     1    tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
     2    ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)))/(16.0d0*
     3    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     4    2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*
     5    ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0)**0.5d0
        A(1.0d0,2.0d0) = (ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
     1    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*(tbar_cr(1)**2.0d0*
     2    sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     3    tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     4    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     5    sign(1.0d0,tbar_cr(1))**2.0d0*
     6    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
     7    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     8    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     9    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     2    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)*(16.0d0*
     5    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     6    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     7    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*lch)**2.0d0)) -
     8    (4.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
     9    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*((0.25d0*
     1    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     2    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     3    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     4    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     5    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     6    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     7    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     8    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     9    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     1    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     2    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     3    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     4    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/
     7    (sign(1.0d0,ec_cr(2)*lch)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     8    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*
     9    lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
     1    sign(1.0d0,ec_cr(2)*lch)**2.0d0)**1.5d0) - (4.0d0*
     2    ec_cr(4)**3.0d0*eta*lch*tbar_cr(1)*tbar_cr(2)*(ec_cr(2)*
     3    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     4    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     5    2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     8    (GIIc - GIc)*(tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     9    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     1    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     2    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     3    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(2.0d0*
     4    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
     5    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     6    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*
     7    ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - ec_cr(2)**2.0d0*
     1    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     4    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     5    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     6    tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*
     7    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     8    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     9    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     1    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)**2.0d0*(4.0d0*
     4    ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
     5    tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     6    tbar_cr(1)**2.0d0)**2.0d0)
        A(1.0d0,3.0d0) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*
     1    tbar_cr(1))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     3    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     4    abs(ec_cr(2))*abs(lch))**1.5d0
        A(2.0d0,1.0d0) = (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((4.0d0*
     1    ec_cr(4)*lch**2.0d0*(abs(tbar_cr(1))**2.0d0 +
     2    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     3    2.0d0))/((GIc + ((tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     6    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     7    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0))**eta)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     9    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     1    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     2    lch)**2.0d0)**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     3    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*
     6    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     7    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     8    1.5d0)*(GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
     9    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     1    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     2    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     3    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     4    2.0d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 -
     5    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     6    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
     7    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     1    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     2    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     3    tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + ((tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     6    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     7    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0))**eta)**2.0d0*(4.0d0*ec_cr(4)**2.0d0*
     9    tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
     1    ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)))/(16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     3    2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*
     4    ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     5    lch)**2.0d0)**0.5d0 - ((0.5d0*abs(ec_cr(2)*lch) - 0.5d0*
     6    ec_cr(2)*lch)*((2.0d0*ec_cr(2)*lch*tbar_cr(2)*((4.0d0*
     7    ec_cr(4)*lch**2.0d0*(abs(tbar_cr(1))**2.0d0 +
     8    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     9    2.0d0))/((GIc + ((tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     3    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     4    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     5    2.0d0))**eta)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     6    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     7    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     8    lch)**2.0d0)**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     9    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*
     3    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     4    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     5    1.5d0)*(GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
     6    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     7    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     8    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     9    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     1    2.0d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 -
     2    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     3    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
     4    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     7    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     8    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     9    tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + ((tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     3    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     4    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     5    2.0d0))**eta)**2.0d0*(4.0d0*ec_cr(4)**2.0d0*
     6    tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
     7    ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)))/(16.0d0*
     8    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     9    2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*
     1    ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     2    lch)**2.0d0)**0.5d0 - (32.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*
     3    tbar_cr(2)*((0.25d0*(abs(tbar_cr(1))**2.0d0 +
     4    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     5    2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     6    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     7    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     8    2.0d0))/(GIc + ((tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     2    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     3    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     4    2.0d0))**eta) - 1.0d0))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     5    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     6    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     7    lch)**2.0d0)**1.5d0))/(ec_cr(2)*lch) - (32.0d0*ec_cr(4)*
     8    ec_cr(2)*lch**3.0d0*tbar_cr(2)*((0.25d0*
     9    (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 +
     1    abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     3    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     4    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))**(1.0d0/
     5    2.0d0))/(GIc + ((tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     7    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     8    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     9    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     1    2.0d0))**eta) - 1.0d0))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     2    ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*
     3    lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     4    lch)**2.0d0)**1.5d0
        A(2.0d0,2.0d0) = n**2.0d0*(C(2,3)*m**2.0d0 + 2.0d0*C(5,5)*
     1    m**2.0d0 + C(3,3)*n**2.0d0) - m**2.0d0*(C(2,2) - C(2,2)*
     2    n**2.0d0 + C(2,3)*n**2.0d0 + 2.0d0*C(5,5)*n**2.0d0) + (0.5d0*
     3    (sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) -
     4    eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*m**2.0d0 + 2.0d0*
     5    eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((0.25d0*
     6    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     7    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     8    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     9    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     1    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     3    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     4    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     5    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     6    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     7    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     9    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/(16.0d0*
     3    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     4    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     5    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     6    lch)**2.0d0)**0.5d0))/ec_cr(2) - (0.5d0*(sign(1.0d0,ec_cr(2)*
     7    lch) - 1.0d0)*(E22 - (2.0d0*lch*tbar_cr(2)*((0.25d0*
     8    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     9    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     1    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     2    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     3    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     4    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     5    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     6    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     7    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     8    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     9    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     2    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/(16.0d0*
     5    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     6    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     7    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     8    lch)**2.0d0)**0.5d0 - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((0.25d0*
     9    ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
     1    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     2    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     3    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     4    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     5    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
     6    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     7    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     8    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     9    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     1    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)*(16.0d0*
     4    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     5    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     6    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     7    lch)**2.0d0)**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*
     8    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*
     2    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     3    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     4    1.5d0)*(GIIc - GIc)*(tbar_cr(1)**2.0d0*
     5    sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     6    tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     7    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     8    sign(1.0d0,tbar_cr(1))**2.0d0*
     9    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     1    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     2    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     3    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     4    lch)**2.0d0)**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0*
     5    tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     7    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*ec_cr(4)**2.0d0*
     8    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     9    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
     1    ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     4    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     5    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     6    tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*
     7    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     8    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     9    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     1    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)**2.0d0*(4.0d0*
     4    ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
     5    tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     6    tbar_cr(1)**2.0d0)**2.0d0)))/(16.0d0*ec_cr(4)**2.0d0*
     7    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 +
     8    (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) +
     9    1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*lch)**2.0d0)**0.5d0 +
     1    (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*tbar_cr(2)*
     2    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*((0.25d0*
     3    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     4    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     5    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     6    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     7    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     8    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     9    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     1    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     2    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     3    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     4    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     5    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     6    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/
     9    (sign(1.0d0,ec_cr(2)*lch)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     1    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*
     2    lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
     3    sign(1.0d0,ec_cr(2)*lch)**2.0d0)**1.5d0)))/
     4    sign(1.0d0,ec_cr(2)*lch) + (2.0d0*lch*tbar_cr(2)*((0.25d0*
     5    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     6    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     7    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     8    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     9    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     1    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     2    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     3    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     4    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     5    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     6    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     8    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/(16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     3    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     4    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     5    lch)**2.0d0)**0.5d0 + (0.5d0*(sign(1.0d0,ec_cr(2)*lch) -
     6    1.0d0)*(E22*(ec_cr(2) - eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*
     7    m**2.0d0 + 2.0d0*eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*
     8    tbar_cr(2)*((0.25d0*(tbar_cr(1)**2.0d0*
     9    sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     1    tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     2    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     3    sign(1.0d0,tbar_cr(1))**2.0d0*
     4    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     5    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     6    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     7    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     8    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     9    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     1    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     3    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/(16.0d0*
     6    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     7    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     8    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     9    lch)**2.0d0)**0.5d0))/(ec_cr(2)*sign(1.0d0,ec_cr(2)*lch)) +
     1    (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((0.25d0*ec_cr(2)*lch**2.0d0*
     2    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*(tbar_cr(1)**2.0d0*
     3    sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     4    tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     5    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     6    sign(1.0d0,tbar_cr(1))**2.0d0*
     7    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
     8    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     9    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     1    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     3    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)*(16.0d0*
     6    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     7    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     8    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     9    lch)**2.0d0)**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*
     1    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     2    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     3    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*
     4    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     5    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     6    1.5d0)*(GIIc - GIc)*(tbar_cr(1)**2.0d0*
     7    sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     8    tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     9    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     1    sign(1.0d0,tbar_cr(1))**2.0d0*
     2    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     3    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     4    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     5    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     6    lch)**2.0d0)**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0*
     7    tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*ec_cr(4)**2.0d0*
     1    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
     3    ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     6    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     7    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     8    tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*
     9    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     1    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     3    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)**2.0d0*(4.0d0*
     6    ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
     7    tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     8    tbar_cr(1)**2.0d0)**2.0d0)))/(16.0d0*ec_cr(4)**2.0d0*
     9    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 +
     1    (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) +
     2    1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*lch)**2.0d0)**0.5d0 -
     3    (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*tbar_cr(2)*
     4    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*((0.25d0*
     5    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     6    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     7    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     8    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     9    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     1    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     2    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     3    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     4    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     5    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     6    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     8    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/
     2    (sign(1.0d0,ec_cr(2)*lch)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     3    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*
     4    lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
     5    sign(1.0d0,ec_cr(2)*lch)**2.0d0)**1.5d0)
        A(2.0d0,3.0d0) = 0.5d0*C(3,3)*2*m*n - 0.5d0*C(2,3)*2*m*n +
     1    0.5d0*C(2,2)*2*m*n*m**2.0d0 - 0.5d0*C(3,3)*2*m*n*m**2.0d0 +
     2    (48*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(2))/
     3    (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     4    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     5    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     6    abs(ec_cr(2))*abs(lch))**1.5d0 - (16.0d0*ec_cr(5)*lch**2.0d0*
     7    tbar_cr(2)*abs(ec_cr(2))*abs(lch))/(abs(ec_cr(2))**2.0d0*
     8    abs(lch)**2.0d0 + 16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     9    ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     1    lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*
     2    abs(lch))**1.5d0
        A(3.0d0,1.0d0) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*
     1    tbar_cr(3))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     2    ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     3    lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*
     4    abs(ec_cr(2)*lch))**1.5d0 + (4.0d0*ec_cr(4)*ec_cr(2)*ec_cr(5)*
     5    eta*lch*tbar_cr(2)*tbar_cr(3)*(ec_cr(2)*tbar_cr(2) +
     6    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*
     7    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     2    (GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
     3    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     4    2.0d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 -
     5    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     6    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
     7    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     1    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     2    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     3    tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + ((tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
     6    tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/
     7    2.0d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0))**eta)**2.0d0*(4.0d0*ec_cr(4)**2.0d0*
     9    tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
     1    ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)
        A(3.0d0,2.0d0) = (4.0d0*ec_cr(5)*lch*tbar_cr(3)*((0.25d0*
     1    ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
     2    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     3    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     4    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     5    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     6    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
     7    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     8    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     9    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     2    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)*(16.0d0*
     5    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     6    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     7    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     8    lch)**2.0d0)**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*
     9    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**(1.0d0 - eta)*
     3    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     4    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     5    1.5d0)*(GIIc - GIc)*(tbar_cr(1)**2.0d0*
     6    sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     7    tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     8    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
     9    sign(1.0d0,tbar_cr(1))**2.0d0*
     1    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     3    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     4    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     5    lch)**2.0d0)**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0*
     6    tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*ec_cr(4)**2.0d0*
     9    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
     2    ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*ec_cr(4)**2.0d0*
     5    ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*
     6    tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     7    tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*
     8    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     9    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     2    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)**2.0d0*(4.0d0*
     5    ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
     6    tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     7    tbar_cr(1)**2.0d0)**2.0d0)))/(16.0d0*ec_cr(4)**2.0d0*
     8    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 +
     9    (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) +
     1    1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*lch)**2.0d0)**0.5d0 - m*n*
     2    (C(2,2) - C(2,2)*n**2.0d0 + C(2,3)*n**2.0d0 + 2.0d0*C(5,5)*
     3    n**2.0d0) - m*n*(C(2,3)*m**2.0d0 + 2.0d0*C(5,5)*m**2.0d0 +
     4    C(3,3)*n**2.0d0) - (4.0d0*ec_cr(2)*ec_cr(5)*lch**3.0d0*
     5    tbar_cr(3)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*((0.25d0*
     6    (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     7    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     8    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     9    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     1    sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
     3    lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
     4    lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*
     5    lch)**2.0d0)**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
     6    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc +
     7    ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/
     9    (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)) - 1.0d0))/
     3    (sign(1.0d0,ec_cr(2)*lch)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
     4    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*
     5    lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
     6    sign(1.0d0,ec_cr(2)*lch)**2.0d0)**1.5d0)
        A(3.0d0,3.0d0) = (64.0d0*ec_cr(5)**2.0d0*lch**3.0d0*
     1    tbar_cr(3))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     2    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     3    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     4    abs(ec_cr(2))*abs(lch))**1.5d0 - C(5,5) - (4.0d0*lch*
     5    tbar_cr(3))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     6    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     7    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     8    abs(ec_cr(2))*abs(lch))**0.5d0 - C(2,2)*n**2.0d0*(n**2.0d0 -
     9    1.0d0) + C(3,3)*n**2.0d0*(n**2.0d0 - 1.0d0) + 2.0d0*C(5,5)*
     1    n**2.0d0 + (lch*tbar_cr(3)*(abs(tbar_cr(1))**2.0d0 +
     2    abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
     3    2.0d0))/(GIc + (GIIc*(tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**(0.5d0*eta))/(ec_cr(2)*tbar_cr(2) +
     6    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*
     7    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     8    2.0d0))**eta - (GIc*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(0.5d0*eta))/(ec_cr(2)*tbar_cr(2) +
     2    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*
     3    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
     4    2.0d0))**eta)
      end function jacobian

      end module jacobian_f77_noterms
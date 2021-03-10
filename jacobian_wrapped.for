        term1 = (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(0.5d0))
        term2 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*
     1    ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*
     2    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
     3    sign(1.0d0,ec_cr(2)*lch)**2.0d0)
        term3 = ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*
     1    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc -
     2    GIc))
        term4 = (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
     1    sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
     2    sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 +
     3    tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*
     4    sign(1.0d0,tbar_cr(2))**2.0d0)
        term5 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     1    lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*
     2    ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0)
        term6 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
     1    lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
     2    lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))
        term7 = ((0.25d0*term4**(0.5d0)*term2**(0.5d0))/
     1    (sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     2    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)) - 1.0d0)
        term8 = (4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*
     1    ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
     2    tbar_cr(1)**2.0d0)
        term9 = (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 +
     1    abs(tbar_cr(3))**2.0d0)
        term10 = (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
     1    ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
     2    16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
     3    abs(ec_cr(2))*abs(lch))
        M(1,1) = C(6,6)*(p**2.0d0 - 1.0d0) - C(4,4)*p**2.0d0 + (4.0d0*
     1    lch*tbar_cr(1)*((0.25d0*term9**(0.5d0)*term6**(0.5d0))/(GIc +
     2    term3/term1**eta) - 1.0d0))/term5**0.5d0 - (64.0d0*
     3    ec_cr(4)**2.0d0*lch**3.0d0*tbar_cr(1)*((0.25d0*term9**(0.5d0)*
     4    term6**(0.5d0))/(GIc + term3/term1**eta) - 1.0d0))/
     5    term5**1.5d0 + (4.0d0*ec_cr(4)*lch*tbar_cr(1)*((4.0d0*
     6    ec_cr(4)*lch**2.0d0*term9**(0.5d0))/((GIc + term3/term1**eta)*
     7    term5**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     8    term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     2    (GIIc - GIc)*term9**(0.5d0)*term6**(0.5d0)*(2.0d0*
     3    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     5    ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*
     6    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*
     7    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*
     8    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
     9    ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
     1    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
     2    ((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term5**0.5d0
        M(1,2) = (ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
     1    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4**(0.5d0))/
     2    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     3    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/
     4    term1**eta)*term2) - (4.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*
     5    tbar_cr(1)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term7)/
     6    (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0) - (4.0d0*
     7    ec_cr(4)**3.0d0*eta*lch*tbar_cr(1)*tbar_cr(2)*term1**(1.0d0 -
     8    eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta -
     9    0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*
     1    eta - 1.5d0)*(GIIc - GIc)*term4**(0.5d0)*(2.0d0*
     2    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
     3    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     4    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*
     5    ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     7    ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*
     8    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 +
     1    ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     2    tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*
     3    tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
     4    tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     5    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*
     6    term8**2.0d0)
        M(1,3) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*tbar_cr(1))/
     1    term10**1.5d0
        M(2,1) = (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((4.0d0*ec_cr(4)*
     1    lch**2.0d0*term9**(0.5d0))/((GIc + term3/term1**eta)*
     2    term5**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     3    term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     4    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     5    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     6    (GIIc - GIc)*term9**(0.5d0)*term6**(0.5d0)*(2.0d0*
     7    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*
     1    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*
     2    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*
     3    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
     4    ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
     5    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
     6    ((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/
     7    term5**0.5d0 - ((0.5d0*abs(ec_cr(2)*lch) - 0.5d0*ec_cr(2)*
     8    lch)*((2.0d0*ec_cr(2)*lch*tbar_cr(2)*((4.0d0*ec_cr(4)*
     9    lch**2.0d0*term9**(0.5d0))/((GIc + term3/term1**eta)*
     1    term5**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
     2    term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     4    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     5    (GIIc - GIc)*term9**(0.5d0)*term6**(0.5d0)*(2.0d0*
     6    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*
     9    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*
     2    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
     3    ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
     4    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
     5    ((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/
     6    term5**0.5d0 - (32.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*
     7    tbar_cr(2)*((0.25d0*term9**(0.5d0)*term6**(0.5d0))/(GIc +
     8    term3/term1**eta) - 1.0d0))/term5**1.5d0))/(ec_cr(2)*lch) -
     9    (32.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*((0.25d0*
     1    term9**(0.5d0)*term6**(0.5d0))/(GIc + term3/term1**eta) -
     2    1.0d0))/term5**1.5d0
        M(2,2) = q**2.0d0*(C(2,3)*p**2.0d0 + 2.0d0*C(5,5)*p**2.0d0 +
     1    C(3,3)*q**2.0d0) - p**2.0d0*(C(2,2) - C(2,2)*q**2.0d0 +
     2    C(2,3)*q**2.0d0 + 2.0d0*C(5,5)*q**2.0d0) + (0.5d0*
     3    (sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) -
     4    eTotal_cr(3)*q**2.0d0 - eTotal_cr(2)*p**2.0d0 + 2.0d0*
     5    eTotal_cr(5)*p*q) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term7)/
     6    term2**0.5d0))/ec_cr(2) - (0.5d0*(sign(1.0d0,ec_cr(2)*lch) -
     7    1.0d0)*(E22 - (2.0d0*lch*tbar_cr(2)*term7)/term2**0.5d0 -
     8    (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((0.25d0*ec_cr(2)*lch**2.0d0*
     9    (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4**(0.5d0))/
     1    (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
     2    sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/
     3    term1**eta)*term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*
     4    term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     6    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     7    (GIIc - GIc)*term4**(0.5d0)*term2**(0.5d0)*(2.0d0*
     8    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
     9    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*
     2    ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
     3    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     4    ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*
     5    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     6    tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 +
     7    ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     8    tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*
     9    tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
     1    tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     2    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*
     3    term8**2.0d0)))/term2**0.5d0 + (2.0d0*ec_cr(2)**2.0d0*
     4    lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) +
     5    1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*
     6    term2**1.5d0)))/sign(1.0d0,ec_cr(2)*lch) + (2.0d0*lch*
     7    tbar_cr(2)*term7)/term2**0.5d0 + (0.5d0*(sign(1.0d0,ec_cr(2)*
     8    lch) - 1.0d0)*(E22*(ec_cr(2) - eTotal_cr(3)*q**2.0d0 -
     9    eTotal_cr(2)*p**2.0d0 + 2.0d0*eTotal_cr(5)*p*q) - (2.0d0*
     1    ec_cr(2)*lch*tbar_cr(2)*term7)/term2**0.5d0))/(ec_cr(2)*
     2    sign(1.0d0,ec_cr(2)*lch)) + (2.0d0*ec_cr(2)*lch*tbar_cr(2)*
     3    ((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) +
     4    1.0d0)**2.0d0*term4**(0.5d0))/(sign(1.0d0,ec_cr(2)*lch)*
     5    sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     6    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)*
     7    term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*
     8    term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     9    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     2    (GIIc - GIc)*term4**(0.5d0)*term2**(0.5d0)*(2.0d0*
     3    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
     4    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     5    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*
     6    ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*
     9    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     1    tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 +
     2    ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     3    tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*
     4    tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
     5    tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     6    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*
     7    term8**2.0d0)))/term2**0.5d0 - (2.0d0*ec_cr(2)**2.0d0*
     8    lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) +
     9    1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
        M(2,3) = 0.5d0*C(3,3)*2*p*q - 0.5d0*C(2,3)*2*p*q + 0.5d0*
     1    C(2,2)*2*p*q*p**2.0d0 - 0.5d0*C(3,3)*2*p*q*p**2.0d0 + (48*
     2    ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(2))/term10**1.5d0 -
     3    (16.0d0*ec_cr(5)*lch**2.0d0*tbar_cr(2)*abs(ec_cr(2))*
     4    abs(lch))/term10**1.5d0
        M(3,1) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*tbar_cr(3))/
     1    term6**1.5d0 + (4.0d0*ec_cr(4)*ec_cr(2)*ec_cr(5)*eta*lch*
     2    tbar_cr(2)*tbar_cr(3)*term1**(1.0d0 - eta)*
     3    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*
     4    (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta -
     5    1.5d0)*(GIIc - GIc)*term9**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*
     6    tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     8    ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*
     9    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*
     1    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*
     2    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
     3    ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
     4    ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
     5    ((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)
        M(3,2) = (4.0d0*ec_cr(5)*lch*tbar_cr(3)*((0.25d0*ec_cr(2)*
     1    lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
     2    term4**(0.5d0))/(sign(1.0d0,ec_cr(2)*lch)*
     3    sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     4    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)*
     5    term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*
     6    term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
     7    tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
     9    (GIIc - GIc)*term4**(0.5d0)*term2**(0.5d0)*(2.0d0*
     1    ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
     2    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
     3    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*
     4    ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*
     7    tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
     8    tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 +
     9    ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*
     1    tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*
     2    tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
     3    tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
     4    sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*
     5    term8**2.0d0)))/term2**0.5d0 - p*q*(C(2,2) - C(2,2)*q**2.0d0 +
     6    C(2,3)*q**2.0d0 + 2.0d0*C(5,5)*q**2.0d0) - p*q*(C(2,3)*
     7    p**2.0d0 + 2.0d0*C(5,5)*p**2.0d0 + C(3,3)*q**2.0d0) - (4.0d0*
     8    ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(3)*(sign(1.0d0,ec_cr(2)*
     9    lch) + 1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*
     1    term2**1.5d0)
        M(3,3) = (64.0d0*ec_cr(5)**2.0d0*lch**3.0d0*tbar_cr(3))/
     1    term10**1.5d0 - C(5,5) - (4.0d0*lch*tbar_cr(3))/
     2    term10**0.5d0 - C(2,2)*q**2.0d0*(q**2.0d0 - 1.0d0) + C(3,3)*
     3    q**2.0d0*(q**2.0d0 - 1.0d0) + 2.0d0*C(5,5)*q**2.0d0 + (lch*
     4    tbar_cr(3)*term9**(0.5d0))/(GIc + (GIIc*(tbar_cr(1)**2.0d0 +
     5    tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
     6    ec_cr(2)**2.0d0)**(0.5d0*eta))/term1**eta - (GIc*
     7    (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*
     8    ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta))/term1**eta)

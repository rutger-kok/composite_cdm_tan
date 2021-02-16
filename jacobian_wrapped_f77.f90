       term1 = (GIc + ((tbar_cr(1)**2.0d0 +
      1   tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 +
      2   ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))/(ec_cr(2)*
      3  tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      4   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      5  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)
       term2 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*
      1  ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*
      2  (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/
      3  sign(1.0d0,ec_cr(2)*lch)**2.0d0)
       term3 = ((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*
      1  lch) + 1.0d0)**2.0d0*(tbar_cr(1)**2.0d0*
      2  sign(1.0d0,tbar_cr(2))**2.0d0*
      3  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
      4  sign(1.0d0,tbar_cr(1))**2.0d0*
      5  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
      6  sign(1.0d0,tbar_cr(1))**2.0d0*
      7  sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
      8  (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
      9  sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1*
      1  term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*(ec_cr(2)*
      2  tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      3   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      4  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
      5  2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
      6   tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
      7  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
      8  (GIIc - GIc)*(tbar_cr(1)**2.0d0*
      9  sign(1.0d0,tbar_cr(2))**2.0d0*
      1  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
      2  sign(1.0d0,tbar_cr(1))**2.0d0*
      3  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
      4  sign(1.0d0,tbar_cr(1))**2.0d0*
      5  sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*
      6  term2**(1.0d0/2.0d0)*(2.0d0*ec_cr(2)**3.0d0*
      7  tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
      8   tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
      9   ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - 4.0d0*ec_cr(4)**2.0d0*
      1  (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
      2  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
      3   ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
      4   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      5  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*
      6  ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
      7  ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
      8  ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
      9  (sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*
      1  sign(1.0d0,tbar_cr(3))*term1**2.0d0*(4.0d0*
      2  ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
      3  tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
      4  tbar_cr(1)**2.0d0)**2.0d0))
       term4 = ((0.25d0*(tbar_cr(1)**2.0d0*
      1  sign(1.0d0,tbar_cr(2))**2.0d0*
      2  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
      3  sign(1.0d0,tbar_cr(1))**2.0d0*
      4  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
      5  sign(1.0d0,tbar_cr(1))**2.0d0*
      6  sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*
      7  term2**(1.0d0/2.0d0))/(sign(1.0d0,tbar_cr(1))*
      8  sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1) -
      9   1.0d0)
       term5 = ((4.0d0*ec_cr(4)*lch**2.0d0*
      1  (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 +
      2   abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0))/ (term1*(16.0d0*
      3  ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
      4   2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*
      5  ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*
      6  lch)**2.0d0)**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*
      7  (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      8   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      9  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
      1  2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
      2   tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
      3  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
      4  (GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
      5   abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
      6  2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
      7   ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
      8  lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*
      9  lch*abs(ec_cr(2)*lch))**(1.0d0/2.0d0)*(2.0d0*
      1  ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
      2   tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
      3   ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
      4  tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
      5   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      6  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) + 8.0d0*
      7  ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*
      8  ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*
      9  ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/
      1  (term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 +
      2   4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
      3   ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0))
       term6 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
      1   ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
      2  abs(ec_cr(2)*lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 +
      3   abs(ec_cr(2)*lch)**2.0d0)
       term7 = ((0.25d0*(abs(tbar_cr(1))**2.0d0 +
      1   abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)**(1.0d0/
      2  2.0d0)*(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
      3   ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
      4  lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*
      5  lch*abs(ec_cr(2)*lch))**(1.0d0/2.0d0))/term1 - 1.0d0)
       A(1,1) = C(6,6)*(m**2.0d0 - 1.0d0) - C(4,4)*m**2.0d0 +
      1   (4.0d0*lch*tbar_cr(1)*term7)/term6**0.5d0 - (64.0d0*
      2  ec_cr(4)**2.0d0*lch**3.0d0*tbar_cr(1)*term7)/
      3  term6**1.5d0 + (4.0d0*ec_cr(4)*lch*tbar_cr(1)*term5)/
      4  term6**0.5d0
       A(1,2) = (ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
      1  (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
      2  (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*
      3  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
      4  sign(1.0d0,tbar_cr(1))**2.0d0*
      5  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
      6  sign(1.0d0,tbar_cr(1))**2.0d0*
      7  sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0))/
      8  (sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*
      9  sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1*
      1  term2) - (4.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*
      2  (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/
      3  (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0) - (4.0d0*
      4  ec_cr(4)**3.0d0*eta*lch*tbar_cr(1)*tbar_cr(2)*(ec_cr(2)*
      5  tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      6   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      7  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
      8  2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
      9   tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
      1  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
      2  (GIIc - GIc)*(tbar_cr(1)**2.0d0*
      3  sign(1.0d0,tbar_cr(2))**2.0d0*
      4  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*
      5  sign(1.0d0,tbar_cr(1))**2.0d0*
      6  sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*
      7  sign(1.0d0,tbar_cr(1))**2.0d0*
      8  sign(1.0d0,tbar_cr(2))**2.0d0)**(1.0d0/2.0d0)*(2.0d0*
      9  ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*
      1  (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*
      2  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) -
      3   4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 +
      4   tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
      5   ec_cr(2)**2.0d0)**(1.0d0/2.0d0) - ec_cr(2)**2.0d0*
      6  tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
      7   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      8  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) +
      9   8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 +
      1   2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) +
      2   8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
      3  tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*
      4  sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*
      5  term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 +
      6   4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 +
      7   ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)
       A(1,3) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*
      1  tbar_cr(1))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 +
      2   16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
      3  lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*
      4  ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0
       A(2,1) = (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term5)/
      1  term6**0.5d0 - ((0.5d0*abs(ec_cr(2)*lch) - 0.5d0*
      2  ec_cr(2)*lch)*((2.0d0*ec_cr(2)*lch*tbar_cr(2)*term5)/
      3  term6**0.5d0 - (32*ec_cr(4)*ec_cr(2)*lch**3.0d0*
      4  tbar_cr(2)*term7)/term6**1.5d0))/(ec_cr(2)*lch) - (32*
      5  ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*term7)/
      6  term6**1.5d0
       A(2,2) = n**2.0d0*(C(2,3)*m**2.0d0 + 2.0d0*C(5,5)*
      1  m**2.0d0 + C(3,3)*n**2.0d0) - m**2.0d0*(C(2,2) - C(2,2)*
      2  n**2.0d0 + C(2,3)*n**2.0d0 + 2.0d0*C(5,5)*n**2.0d0) +
      3   (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*
      4  (ec_cr(2) - eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*
      5  m**2.0d0 + 2.0d0*eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*
      6  tbar_cr(2)*term4)/term2**0.5d0))/ec_cr(2) - (0.5d0*
      7  (sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22 - (2.0d0*lch*
      8  tbar_cr(2)*term4)/term2**0.5d0 - (2.0d0*ec_cr(2)*lch*
      9  tbar_cr(2)*term3)/term2**0.5d0 + (2.0d0*ec_cr(2)**2.0d0*
      1  lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) +
      2   1.0d0)**2.0d0*term4)/(sign(1.0d0,ec_cr(2)*lch)*
      3  term2**1.5d0)))/sign(1.0d0,ec_cr(2)*lch) + (2.0d0*lch*
      4  tbar_cr(2)*term4)/term2**0.5d0 + (0.5d0*
      5  (sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) -
      6   eTotal_cr(3)*n**2.0d0 - eTotal_cr(2)*m**2.0d0 + 2.0d0*
      7  eTotal_cr(5)*m*n) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*
      8  term4)/term2**0.5d0))/(ec_cr(2)*sign(1.0d0,ec_cr(2)*
      9  lch)) + (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term3)/
      1  term2**0.5d0 - (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*
      2  tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*
      3  term4)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
       A(2,3) = 0.5d0*C(3,3)*2.0d0*m*n - 0.5d0*C(2,3)*2.0d0*m*
      1  n + 0.5d0*C(2,2)*2.0d0*m*n*m**2.0d0 - 0.5d0*C(3,3)*2.0d0*
      2  m*n*m**2.0d0 + (48*ec_cr(2)*ec_cr(5)*lch**3.0d0*
      3  tbar_cr(2))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 +
      4   16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
      5  lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*
      6  ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**1.5d0 - (16.0d0*
      7  ec_cr(5)*lch**2.0d0*tbar_cr(2)*abs(ec_cr(2))*abs(lch))/
      8  (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
      9  ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
      1   16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
      2  abs(ec_cr(2))*abs(lch))**1.5d0
       A(3,1) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*
      1  tbar_cr(3))/(16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 +
      2   ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*
      3  lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*
      4  lch*abs(ec_cr(2)*lch))**1.5d0 + (4.0d0*ec_cr(4)*ec_cr(2)*
      5  ec_cr(5)*eta*lch*tbar_cr(2)*tbar_cr(3)*(ec_cr(2)*
      6  tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      7   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      8  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/
      9  2.0d0))**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 +
      1   tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*
      2  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*
      3  (GIIc - GIc)*(abs(tbar_cr(1))**2.0d0 +
      4   abs(tbar_cr(2))**2.0d0 +
      5   abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0)*(2.0d0*
      6  ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 +
      7   tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 +
      8   ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*
      9  tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 +
      1   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      2  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0) +
      3   8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 +
      4   2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) +
      5   8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*
      6  tbar_cr(2)))/(term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*
      7  tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*
      8  tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*
      9  tbar_cr(1)**2.0d0)**2.0d0)
       A(3,2) = (4.0d0*ec_cr(5)*lch*tbar_cr(3)*term3)/
      1  term2**0.5d0 - m*n*(C(2,2) - C(2,2)*n**2.0d0 + C(2,3)*
      2  n**2.0d0 + 2.0d0*C(5,5)*n**2.0d0) - m*n*(C(2,3)*
      3  m**2.0d0 + 2.0d0*C(5,5)*m**2.0d0 + C(3,3)*n**2.0d0) -
      4   (4.0d0*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(3)*
      5  (sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4)/
      6  (sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
       A(3,3) = (64.0d0*ec_cr(5)**2.0d0*lch**3.0d0*tbar_cr(3))/
      1  (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*
      2  ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 +
      3   16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*
      4  abs(ec_cr(2))*abs(lch))**1.5d0 - C(5,5) - (4.0d0*lch*
      5  tbar_cr(3))/(abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 +
      6   16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*
      7  lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*
      8  ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))**0.5d0 - C(2,2)*
      9  n**2.0d0*(n**2.0d0 - 1.0d0) + C(3,3)*n**2.0d0*(n**2.0d0 -
      1   1.0d0) + 2.0d0*C(5,5)*n**2.0d0 + (lch*tbar_cr(3)*
      2  (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 +
      3   abs(tbar_cr(3))**2.0d0)**(1.0d0/2.0d0))/(GIc + (GIIc*
      4  (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*
      5  (4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta))/
      6  (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      7   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      8  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta -
      9   (GIc*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*
      1  eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*
      2  eta))/(ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 +
      3   tbar_cr(2)**2.0d0)**(1.0d0/2.0d0)*(4.0d0*
      4  ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(1.0d0/2.0d0))**eta)

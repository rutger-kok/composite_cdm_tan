term1 = (ec_cr(2)*tbar_cr(2) + (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0))
term2 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + (ec_cr(2)**2.0d0*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0)/sign(1.0d0,ec_cr(2)*lch)**2.0d0)
term3 = ((tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta)*(GIIc - GIc))
term4 = (tbar_cr(1)**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(2)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(3))**2.0d0 + tbar_cr(3)**2.0d0*sign(1.0d0,tbar_cr(1))**2.0d0*sign(1.0d0,tbar_cr(2))**2.0d0)
term5 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch) + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0)
term6 = (16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + abs(ec_cr(2)*lch)**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2)*lch))
term7 = ((0.25d0*term4**(0.5d0)*term2**(0.5d0))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)) - 1.0d0)
term8 = (4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + 4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)
term9 = (abs(tbar_cr(1))**2.0d0 + abs(tbar_cr(2))**2.0d0 + abs(tbar_cr(3))**2.0d0)
term10 = (abs(ec_cr(2))**2.0d0*abs(lch)**2.0d0 + 16.0d0*ec_cr(4)**2.0d0*lch**2.0d0 + ec_cr(2)**2.0d0*lch**2.0d0 + 16.0d0*ec_cr(5)**2.0d0*lch**2.0d0 + 2.0d0*ec_cr(2)*lch*abs(ec_cr(2))*abs(lch))
M(1,1) = C(6,6)*(p**2.0d0 - 1.0d0) - C(4,4)*p**2.0d0 + (4.0d0*lch*tbar_cr(1)*((0.25d0*term9**(0.5d0)*term6**(0.5d0))/(GIc + term3/term1**eta) - 1.0d0))/term5**0.5d0 - (64.0d0*ec_cr(4)**2.0d0*lch**3.0d0*tbar_cr(1)*((0.25d0*term9**(0.5d0)*term6**(0.5d0))/(GIc + term3/term1**eta) - 1.0d0))/term5**1.5d0 + (4.0d0*ec_cr(4)*lch*tbar_cr(1)*((4.0d0*ec_cr(4)*lch**2.0d0*term9**(0.5d0))/((GIc + term3/term1**eta)*term5**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term9**(0.5d0)*term6**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term5**0.5d0
M(1,2) = (ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4**(0.5d0))/(sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)*term2) - (4.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(1)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0) - (4.0d0*ec_cr(4)**3.0d0*eta*lch*tbar_cr(1)*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term4**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*term8**2.0d0)
M(1,3) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*tbar_cr(1))/term10**1.5d0
M(2,1) = (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((4.0d0*ec_cr(4)*lch**2.0d0*term9**(0.5d0))/((GIc + term3/term1**eta)*term5**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term9**(0.5d0)*term6**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term5**0.5d0 - ((0.5d0*abs(ec_cr(2)*lch) - 0.5d0*ec_cr(2)*lch)*((2.0d0*ec_cr(2)*lch*tbar_cr(2)*((4.0d0*ec_cr(4)*lch**2.0d0*term9**(0.5d0))/((GIc + term3/term1**eta)*term5**0.5d0) + (ec_cr(4)*ec_cr(2)*eta*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term9**(0.5d0)*term6**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term5**0.5d0 - (32.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*((0.25d0*term9**(0.5d0)*term6**(0.5d0))/(GIc + term3/term1**eta) - 1.0d0))/term5**1.5d0))/(ec_cr(2)*lch) - (32.0d0*ec_cr(4)*ec_cr(2)*lch**3.0d0*tbar_cr(2)*((0.25d0*term9**(0.5d0)*term6**(0.5d0))/(GIc + term3/term1**eta) - 1.0d0))/term5**1.5d0
M(2,2) = q**2.0d0*(C(2,3)*p**2.0d0 + 2.0d0*C(5,5)*p**2.0d0 + C(3,3)*q**2.0d0) - p**2.0d0*(C(2,2) - C(2,2)*q**2.0d0 + C(2,3)*q**2.0d0 + 2.0d0*C(5,5)*q**2.0d0) + (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) - eTotal_cr(3)*q**2.0d0 - eTotal_cr(2)*p**2.0d0 + 2.0d0*eTotal_cr(5)*p*q) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term7)/term2**0.5d0))/ec_cr(2) - (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22 - (2.0d0*lch*tbar_cr(2)*term7)/term2**0.5d0 - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4**(0.5d0))/(sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)*term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term4**(0.5d0)*term2**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term2**0.5d0 + (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)))/sign(1.0d0,ec_cr(2)*lch) + (2.0d0*lch*tbar_cr(2)*term7)/term2**0.5d0 + (0.5d0*(sign(1.0d0,ec_cr(2)*lch) - 1.0d0)*(E22*(ec_cr(2) - eTotal_cr(3)*q**2.0d0 - eTotal_cr(2)*p**2.0d0 + 2.0d0*eTotal_cr(5)*p*q) - (2.0d0*ec_cr(2)*lch*tbar_cr(2)*term7)/term2**0.5d0))/(ec_cr(2)*sign(1.0d0,ec_cr(2)*lch)) + (2.0d0*ec_cr(2)*lch*tbar_cr(2)*((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4**(0.5d0))/(sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)*term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term4**(0.5d0)*term2**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term2**0.5d0 - (2.0d0*ec_cr(2)**2.0d0*lch**3.0d0*tbar_cr(2)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
M(2,3) = 0.5d0*C(3,3)*2*p*q - 0.5d0*C(2,3)*2*p*q + 0.5d0*C(2,2)*2*p*q*p**2.0d0 - 0.5d0*C(3,3)*2*p*q*p**2.0d0 + (48*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(2))/term10**1.5d0 - (16.0d0*ec_cr(5)*lch**2.0d0*tbar_cr(2)*abs(ec_cr(2))*abs(lch))/term10**1.5d0
M(3,1) = (64.0d0*ec_cr(4)*ec_cr(5)*lch**3.0d0*tbar_cr(3))/term6**1.5d0 + (4.0d0*ec_cr(4)*ec_cr(2)*ec_cr(5)*eta*lch*tbar_cr(2)*tbar_cr(3)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term9**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - (tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**1.5d0 - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/((GIc + term3/term1**eta)**2.0d0*term8**2.0d0)
M(3,2) = (4.0d0*ec_cr(5)*lch*tbar_cr(3)*((0.25d0*ec_cr(2)*lch**2.0d0*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term4**(0.5d0))/(sign(1.0d0,ec_cr(2)*lch)*sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)*term2**0.5d0) - (ec_cr(4)**2.0d0*eta*tbar_cr(2)*term1**(1.0d0 - eta)*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta - 0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta - 1.5d0)*(GIIc - GIc)*term4**(0.5d0)*term2**(0.5d0)*(2.0d0*ec_cr(2)**3.0d0*tbar_cr(2)**3.0d0 - ec_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - 4.0d0*ec_cr(4)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**1.5d0*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) - ec_cr(2)**2.0d0*tbar_cr(2)**2.0d0*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(2)**3.0d0 + 2.0d0*ec_cr(2)**3.0d0*tbar_cr(1)**2.0d0*tbar_cr(2) + 8.0d0*ec_cr(4)**2.0d0*ec_cr(2)*tbar_cr(1)**2.0d0*tbar_cr(2)))/(sign(1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*(GIc + term3/term1**eta)**2.0d0*term8**2.0d0)))/term2**0.5d0 - p*q*(C(2,2) - C(2,2)*q**2.0d0 + C(2,3)*q**2.0d0 + 2.0d0*C(5,5)*q**2.0d0) - p*q*(C(2,3)*p**2.0d0 + 2.0d0*C(5,5)*p**2.0d0 + C(3,3)*q**2.0d0) - (4.0d0*ec_cr(2)*ec_cr(5)*lch**3.0d0*tbar_cr(3)*(sign(1.0d0,ec_cr(2)*lch) + 1.0d0)**2.0d0*term7)/(sign(1.0d0,ec_cr(2)*lch)*term2**1.5d0)
M(3,3) = (64.0d0*ec_cr(5)**2.0d0*lch**3.0d0*tbar_cr(3))/term10**1.5d0 - C(5,5) - (4.0d0*lch*tbar_cr(3))/term10**0.5d0 - C(2,2)*q**2.0d0*(q**2.0d0 - 1.0d0) + C(3,3)*q**2.0d0*(q**2.0d0 - 1.0d0) + 2.0d0*C(5,5)*q**2.0d0 + (lch*tbar_cr(3)*term9**(0.5d0))/(GIc + (GIIc*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta))/term1**eta - (GIc*(tbar_cr(1)**2.0d0 + tbar_cr(2)**2.0d0)**(0.5d0*eta)*(4.0d0*ec_cr(4)**2.0d0 + ec_cr(2)**2.0d0)**(0.5d0*eta))/term1**eta)

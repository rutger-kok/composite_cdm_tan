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
         1.0d0,tbar_cr(1))*sign(1.0d0,tbar_cr(2))*sign(1.0d0,tbar_cr(3))*term1**2.0d0*(4.0d0*ec_cr(4)**2.0d0*tbar_cr(1)**2.0d0 + &
         4.0d0*ec_cr(4)**2.0d0*tbar_cr(2)**2.0d0 + ec_cr(2)**2.0d0*tbar_cr(1)**2.0d0)**2.0d0)
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

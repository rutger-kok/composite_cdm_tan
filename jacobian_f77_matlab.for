      module jacobian_f77_matlab

        ! no variables

      contains
  
      function jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,
     1                  E22,eTotal_cr) result(M)
        implicit none
        ! input variables
        real*8, dimension(6,6), intent(in) :: C
        real*8, dimension(6), intent(in) :: ec_cr, eTotal_cr
        real*8, dimension(3), intent(in) :: tbar_cr
        real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
        ! local variables
        real*8 :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
        real*8 :: t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
        real*8 :: t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42
        real*8 :: t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55
        real*8 :: t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68
        real*8 :: t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81
        real*8 :: t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94
        real*8 :: t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105
        real*8 :: t106,t107,t108,t109,t110,t111,t112,t113,t114,t115
        real*8 :: t116,t117,t118,t119,t120,t121,t122,t123,t124,t125
        real*8 :: t126,t127,t128,t129,t130,t131
        ! output variables
        real*8, dimension(3,3) :: M
        
        t2 = cos(alpha)
        t3 = sin(alpha)
        t4 = lch**2.0D0
        t5 = ec_cr(4)**2.0D0
        t6 = t4*t5*4.0D0
        t7 = ec_cr(5)**2.0D0
        t8 = t4*t7*4.0D0
        t9 = t6+t8
        t10 = sqrt(t9)
        t11 = ec_cr(2)*lch
        t12 = t10+t11
        t13 = 1.0D0/t12
        t14 = tbar_cr(1)**2.0D0
        t15 = tbar_cr(3)**2.0D0
        t16 = t14+t15
        t17 = sqrt(t16)
        t22 = abs(t11)
        t23 = t22/2.0D0
        t24 = (ec_cr(2)*lch)/2.0D0
        t18 = t23+t24
        t19 = abs(tbar_cr(1))
        t20 = abs(tbar_cr(2))
        t21 = abs(tbar_cr(3))
        t25 = t18**2.0D0
        t26 = t6+t8+t25
        t27 = GIc*2.0D0
        t28 = GIIc-GIc
        t44 = t17-tbar_cr(2)
        t29 = t10*t13*t44
        t30 = t29+tbar_cr(2)
        t31 = 1.0D0/t30
        t32 = t10*t13*t17*t31
        t33 = t32**eta
        t34 = t28*t33*2.0D0
        t35 = t27+t34
        t36 = 1.0D0/t35
        t37 = 1.0D0/sqrt(t26)
        t38 = t19**2.0D0
        t39 = t20**2.0D0
        t40 = t21**2.0D0
        t41 = t38+t39+t40
        t42 = sqrt(t41)
        t43 = sqrt(t26)
        t45 = 1.0D0/t12**2.0D0
        t46 = 1.0D0/sqrt(t9)
        t47 = t10*t13*(t17-tbar_cr(2))
        t48 = t47+tbar_cr(2)
        t49 = 1.0D0/t48
        t50 = t10*t13*t17*t49
        t51 = t50**eta
        t52 = t28*t51*2.0D0
        t53 = t27+t52
        t54 = 1.0D0/t53
        t55 = -t17+tbar_cr(2)
        t57 = t10*t13*t55
        t56 = -t57+tbar_cr(2)
        t58 = 1.0D0/t56
        t59 = t10*t13*t17*t58
        t60 = eta-1.0D0
        t61 = lch/2.0D0
        t62 = (t11/abs(t11))
        t63 = (lch*t62)/2.0D0
        t64 = t61+t63
        t65 = t59**eta
        t66 = t28*t65*2.0D0
        t67 = t27+t66
        t68 = 1.0D0/t26**(3.0D0/2.0D0)
        t69 = 1.0D0/t67
        t70 = 1.0D0/t67**2.0D0
        t71 = 1.0D0/t56**2.0D0
        t72 = t59**t60
        t73 = t42*t43*t69
        t74 = t73-1.0D0
        t75 = ec_cr(4)*t4*t37*t42*t69*4.0D0
        t76 = ec_cr(4)*t4*t45*t55*4.0D0
        t113 = ec_cr(4)*t4*t13*t46*t55*4.0D0
        t77 = t76-t113
        t78 = t10*t13*t17*t71*t77
        t79 = ec_cr(4)*t4*t17*t45*t58*4.0D0
        t114 = ec_cr(4)*t4*t13*t17*t46*t58*4.0D0
        t80 = t78+t79-t114
        t81 = eta*t28*t42*t43*t70*t72*t80*2.0D0
        t82 = t75+t81
        t83 = ec_cr(2)*lch*t37*t82*tbar_cr(2)
        t84 = t2**2.0D0
        t85 = t3**2.0D0
        t86 = 1.0D0/ec_cr(2)
        t87 = 1.0D0/lch
        t88 = t23-t24
        t89 = 1.0D0/t12**3
        t90 = lch*t9*t17*t55*t71*t89
        t91 = eTotal_cr(5)*t2*t3*2.0D0
        t92 = ec_cr(2)+t91-eTotal_cr(2)*t84-eTotal_cr(3)*t85
        t93 = E22*t92
        t94 = t93-ec_cr(2)*lch*t37*t74*tbar_cr(2)
        t95 = t18*t37*t42*t64*t69
        t96 = lch*t10*t17*t45*t58
        t97 = t90+t96
        t98 = eta*t28*t42*t43*t70*t72*t97*2.0D0
        t99 = t95+t98
        t100 = ec_cr(2)*lch*t18*t64*t68*t74*tbar_cr(2)
        t101 = t84-t85
        t102 = ec_cr(5)*t4*t37*t42*t69*4.0D0
        t103 = ec_cr(5)*t4*t45*t55*4.0D0
        t110 = ec_cr(5)*t4*t13*t46*t55*4.0D0
        t104 = t103-t110
        t105 = t10*t13*t17*t71*t104
        t106 = ec_cr(5)*t4*t17*t45*t58*4.0D0
        t111 = ec_cr(5)*t4*t13*t17*t46*t58*4.0D0
        t107 = t105+t106-t111
        t108 = eta*t28*t42*t43*t70*t72*t107*2.0D0
        t109 = t102+t108
        t112 = ec_cr(2)*lch*t37*t109*tbar_cr(2)
        t115 = C(2,2)*t84
        t116 = C(2,3)*t85
        t117 = t115+t116
        t118 = t2*t117
        t119 = C(5,5)*t2*t85*2.0D0
        t120 = t118+t119
        t121 = C(2,3)*t84
        t122 = C(3,3)*t85
        t123 = t121+t122
        t124 = t3*t123
        t125 = C(5,5)*t3*t84*2.0D0
        t126 = t124+t125
        t127 = C(2,3)*t2*t3
        t128 = t127-C(3,3)*t2*t3
        t129 = t3*t128
        t130 = t129-C(5,5)*t2*t101
        t131 = C(2,2)*t2*t3
        M(1,1) = -C(4,4)*t84-C(6,6)*t85+lch*t37*tbar_cr(1)*(t36*t42*
     1    t43-1.0D0)*2.0D0+ec_cr(4)*lch*t37*tbar_cr(1)*(ec_cr(4)*t4*t36*
     2    t37*t42*4.0D0-eta*t28*1.0D0/t35**2.0D0*t42*t43*t50**t60*
     3    (ec_cr(4)*t4*t17*t45*t49*(-4.0D0)+t10*t13*t17*1.0D0/
     4    t30**2.0D0*(ec_cr(4)*t4*t44*t45*4.0D0-ec_cr(4)*t4*t13*t44*t46*
     5    4.0D0)+ec_cr(4)*t4*t13*t17*t46*t49*4.0D0)*2.0D0)*2.0D0-lch*t4*
     6    t5*t68*tbar_cr(1)*(t42*t43*t54-1.0D0)*8.0D0
        M(1,2) = ec_cr(4)*lch*t37*tbar_cr(1)*(t18*t37*t42*t54*t64+eta*
     1    t28*t42*t43*t70*t72*(t90+lch*t10*t17*t45*t49)*2.0D0)*2.0D0-
     2    ec_cr(4)*lch*t18*t64*t68*t74*tbar_cr(1)*2.0D0
        M(1,3) = ec_cr(4)*lch*t37*t109*tbar_cr(1)*2.0D0-ec_cr(4)*
     1    ec_cr(5)*lch*t4*t68*t74*tbar_cr(1)*8.0D0
        M(2,1) = t83-t86*t87*t88*(t83-ec_cr(4)*ec_cr(2)*lch*t4*t68*t74*
     1    tbar_cr(2)*4.0D0)-ec_cr(4)*ec_cr(2)*lch*t4*t68*t74*tbar_cr(2)*
     2    4.0D0
        M(2,2) = -t100-t2*t120+t3*t126+lch*t37*t74*tbar_cr(2)-t86*t87*
     1    t94*(t61-t63)-1.0D0/ec_cr(2)**2.0D0*t87*t88*t94+t86*t87*t88*
     2    (E22+t100-lch*t37*t74*tbar_cr(2)-ec_cr(2)*lch*t37*t99*
     3    tbar_cr(2))+ec_cr(2)*lch*t37*t99*tbar_cr(2)
        M(2,3) = t112-t3*t130+t2*(t2*(t131-C(2,3)*t2*t3)-C(5,5)*t3*
     1    t101)-t86*t87*t88*(t112-ec_cr(2)*ec_cr(5)*lch*t4*t68*t74*
     2    tbar_cr(2)*4.0D0)-ec_cr(2)*ec_cr(5)*lch*t4*t68*t74*tbar_cr(2)*
     3    4.0D0
        M(3,1) = ec_cr(5)*lch*t37*t82*tbar_cr(3)*2.0D0-ec_cr(4)*
     1    ec_cr(5)*lch*t4*t68*t74*tbar_cr(3)*8.0D0
        M(3,2) = -t3*t120-t2*t126+ec_cr(5)*lch*t37*t99*tbar_cr(3)*
     1    2.0D0-ec_cr(5)*lch*t18*t64*t68*t74*tbar_cr(3)*2.0D0
        M(3,3) = -t3*(t2*(t127-t131)+C(5,5)*t3*t101)+t2*t130+lch*t37*
     1    t74*tbar_cr(3)*2.0D0+ec_cr(5)*lch*t37*t109*tbar_cr(3)*2.0D0-
     2    lch*t4*t7*t68*t74*tbar_cr(3)*8.0D0
  
      end function jacobian

      end module jacobian_f77_matlab
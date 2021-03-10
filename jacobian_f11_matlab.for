      module jacobian_f77_matlab

        ! no variables

      contains
  
      function jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,
     1                  E22,eTotal_cr) result(M)
        implicit none
        ! input variables
        real*8, dimension(6,6), intent(in) :: C
        real*16, dimension(6), intent(in) :: ec_cr, eTotal_cr
        real*16, dimension(3), intent(in) :: tbar_cr
        real*8, intent(in) :: alpha, lch, eta, GIc, GIIc, E22
        ! local variables
        real*16 :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
        real*16 :: t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
        real*16 :: t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42
        real*16 :: t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55
        real*16 :: t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68
        real*16 :: t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81
        real*16 :: t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94
        real*16 :: t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105
        real*16 :: t106,t107,t108,t109,t110,t111,t112,t113,t114,t115
        real*16 :: t116,t117,t118,t119,t120,t121
        ! output variables
        real*16, dimension(3,3) :: M
        
        t2 = cos(alpha)
        t3 = sin(alpha)
        t4 = lch**2.0D0
        t5 = ec_cr(4)**2.0D0
        t6 = t4*t5*4.0D0
        t7 = ec_cr(2)**2.0D0
        t8 = t4*t7
        t9 = t6+t8
        t10 = sqrt(t9)
        t11 = ec_cr(2)*lch
        t12 = t10+t11
        t13 = 1.0D0/t12
        t14 = tbar_cr(1)**2.0D0
        t15 = tbar_cr(2)**2.0D0
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
        t26 = ec_cr(5)**2.0D0
        t27 = t4*t26*4.0D0
        t28 = t6+t25+t27
        t29 = GIc*2.0D0
        t30 = GIIc-GIc
        t46 = t17-tbar_cr(2)
        t31 = t10*t13*t46
        t32 = t31+tbar_cr(2)
        t33 = 1.0D0/t32
        t34 = t10*t13*t17*t33
        t35 = t34**eta
        t36 = t30*t35*2.0D0
        t37 = t29+t36
        t38 = 1.0D0/t37
        t39 = 1.0D0/sqrt(t28)
        t40 = t19**2.0D0
        t41 = t20**2.0D0
        t42 = t21**2.0D0
        t43 = t40+t41+t42
        t44 = sqrt(t43)
        t45 = sqrt(t28)
        t47 = 1.0D0/t12**2.0D0
        t48 = 1.0D0/sqrt(t9)
        t49 = t10*t13*(t17-tbar_cr(2))
        t50 = t49+tbar_cr(2)
        t51 = 1.0D0/t50
        t52 = t10*t13*t17*t51
        t53 = t52**eta
        t54 = t30*t53*2.0D0
        t55 = t29+t54
        t56 = 1.0D0/t55
        t57 = ec_cr(2)*t4*t48
        t58 = lch+t57
        t59 = eta-1.0D0
        t60 = t52**t59
        t61 = lch/2.0D0
        t62 = (t11/abs(t11))
        t63 = (lch*t62)/2.0D0
        t64 = t61+t63
        t65 = t44*t45*t56
        t66 = t65-1.0D0
        t67 = 1.0D0/t28**(3.0D0/2.0D0)
        t68 = 1.0D0/t55**2.0D0
        t69 = ec_cr(4)*t4*t46*t47*4.0D0
        t70 = 1.0D0/t50**2.0D0
        t74 = ec_cr(4)*t4*t13*t46*t48*4.0D0
        t71 = t69-t74
        t72 = ec_cr(4)*t4*t13*t17*t48*t51*4.0D0
        t73 = ec_cr(4)*t4*t39*t44*t56*4.0D0
        t75 = t10*t13*t17*t70*t71
        t103 = ec_cr(4)*t4*t17*t47*t51*4.0D0
        t76 = t72+t75-t103
        t104 = eta*t30*t44*t45*t60*t68*t76*2.0D0
        t77 = t73-t104
        t78 = ec_cr(2)*lch*t39*t77*tbar_cr(2)
        t79 = t2**2.0D0
        t80 = t3**2.0D0
        t81 = 1.0D0/ec_cr(2)
        t82 = 1.0D0/lch
        t83 = t23-t24
        t84 = lch*t39*t66*tbar_cr(2)
        t85 = t18*t39*t44*t56*t64
        t86 = t10*t46*t47*t58
        t96 = ec_cr(2)*t4*t13*t46*t48
        t87 = t86-t96
        t88 = t10*t13*t17*t70*t87
        t89 = ec_cr(2)*t4*t13*t17*t48*t51
        t97 = t10*t17*t47*t51*t58
        t90 = t88+t89-t97
        t98 = eta*t30*t44*t45*t60*t68*t90*2.0D0
        t91 = t85-t98
        t92 = eTotal_cr(5)*t2*t3*2.0D0
        t93 = ec_cr(2)+t92-eTotal_cr(2)*t79-eTotal_cr(3)*t80
        t94 = E22*t93
        t95 = t94-ec_cr(2)*lch*t39*t66*tbar_cr(2)
        t99 = ec_cr(2)*lch*t18*t64*t66*t67*tbar_cr(2)
        t100 = t79-t80
        t101 = 1.0D0/t28
        t102 = ec_cr(2)*ec_cr(5)*lch*t4*t66*t67*tbar_cr(2)*4.0D0
        t105 = C(2,2)*t79
        t106 = C(2,3)*t80
        t107 = t105+t106
        t108 = t2*t107
        t109 = C(5,5)*t2*t80*2.0D0
        t110 = t108+t109
        t111 = C(2,3)*t79
        t112 = C(3,3)*t80
        t113 = t111+t112
        t114 = t3*t113
        t115 = C(5,5)*t3*t79*2.0D0
        t116 = t114+t115
        t117 = C(2,3)*t2*t3
        t118 = t117-C(3,3)*t2*t3
        t119 = t3*t118
        t120 = t119-C(5,5)*t2*t100
        t121 = C(2,2)*t2*t3

        M(1,1) = -C(4,4)*t79-C(6,6)*t80+lch*t39*tbar_cr(1)*
     1    (t38*t44*t45-1.0D0)*2.0D0+ec_cr(4)*lch*t39*tbar_cr(1)*
     2    (ec_cr(4)*t4*t38*t39*t44*4.0D0-eta*t30*1.0D0/t37**2.0d0*t44*
     3    t45*t60*(t72+t10*t13*t17*1.0D0/t32**2.0D0*t71-ec_cr(4)*t4*t17*
     4    t47*t51*4.0D0)*2.0D0)*2.0D0-lch*t4*t5*t66*t67*tbar_cr(1)*8.0D0
        M(1,2) = ec_cr(4)*lch*t39*t91*tbar_cr(1)*2.0D0-ec_cr(4)*lch*
     1    t18*t64*t66*t67*tbar_cr(1)*2.0D0
        M(1,3) = ec_cr(4)*ec_cr(5)*lch*t4*t66*t67*tbar_cr(1)*(-8.0D0)+
     1    ec_cr(4)*ec_cr(5)*lch*t4*t44*t56*t101*tbar_cr(1)*8.0D0
        M(2,1) = t78-t81*t82*t83*(t78-ec_cr(4)*ec_cr(2)*lch*t4*t66*t67*
     1   tbar_cr(2)*4.0D0)-ec_cr(4)*ec_cr(2)*lch*t4*t66*t67*tbar_cr(2)*
     2   4.0D0
        M(2,2) = t84-t99-t2*t110+t3*t116+t81*t82*t83*(E22-t84+t99-
     1   ec_cr(2)*lch*t39*t91*tbar_cr(2))-t81*t82*t95*(t61-t63)-1.0D0/
     2   ec_cr(2)**2.0D0*t82*t83*t95+ec_cr(2)*lch*t39*t91*tbar_cr(2)
        M(2,3) = -t102-t3*t120+t2*(t2*(t121-C(2,3)*t2*t3)-C(5,5)*t3*
     1    t100)+t81*t82*t83*(t102-ec_cr(2)*ec_cr(5)*lch*t4*t44*t56*t101*
     2    tbar_cr(2)*4.0D0)+ec_cr(2)*ec_cr(5)*lch*t4*t44*t56*t101*
     3    tbar_cr(2)*4.0D0
        M(3,1) = ec_cr(5)*lch*t39*t77*tbar_cr(3)*2.0D0-ec_cr(4)*
     1    ec_cr(5)*lch*t4*t66*t67*tbar_cr(3)*8.0D0
        M(3,2) = -t3*t110-t2*t116+ec_cr(5)*lch*t39*t91*tbar_cr(3)*
     1    2.0D0-ec_cr(5)*lch*t18*t64*t66*t67*tbar_cr(3)*2.0D0
        M(3,3) = -t3*(t2*(t117-t121)+C(5,5)*t3*t100)+t2*t120+lch*t39*
     1    t66*tbar_cr(3)*2.0D0-lch*t4*t26*t66*t67*tbar_cr(3)*8.0D0+lch*
     2    t4*t26*t44*t56*t101*tbar_cr(3)*8.0D0
      end function jacobian

      end module jacobian_f77_matlab
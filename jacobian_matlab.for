      t2 = cos(x)
      t3 = sin(x)
      t4 = lch**2
      t5 = ec12**2
      t6 = t4*t5*4.0D0
      t7 = ec22**2
      t8 = t4*t7
      t9 = t6+t8
      t10 = sqrt(t9)
      t11 = ec22*lch
      t12 = t10+t11
      t13 = 1.0D0/t12
      t14 = tbar_1**2
      t15 = tbar_2**2
      t16 = t14+t15
      t17 = sqrt(t16)
      t22 = abs(t11)
      t23 = t22/2.0D0
      t24 = (ec22*lch)/2.0D0
      t18 = t23+t24
      t19 = abs(tbar_1)
      t20 = abs(tbar_2)
      t21 = abs(tbar_3)
      t25 = t18**2
      t26 = ec23**2
      t27 = t4*t26*4.0D0
      t28 = t6+t25+t27
      t29 = GIc*2.0D0
      t30 = GIIc-GIc
      t46 = t17-tbar_2
      t31 = t10*t13*t46
      t32 = t31+tbar_2
      t33 = 1.0D0/t32
      t34 = t10*t13*t17*t33
      t35 = t34**eta
      t36 = t30*t35*2.0D0
      t37 = t29+t36
      t38 = 1.0D0/t37
      t39 = 1.0D0/sqrt(t28)
      t40 = t19**2
      t41 = t20**2
      t42 = t21**2
      t43 = t40+t41+t42
      t44 = sqrt(t43)
      t45 = sqrt(t28)
      t47 = 1.0D0/t12**2
      t48 = 1.0D0/sqrt(t9)
      t49 = t10*t13*(t17-tbar_2)
      t50 = t49+tbar_2
      t51 = 1.0D0/t50
      t52 = t10*t13*t17*t51
      t53 = t52**eta
      t54 = t30*t53*2.0D0
      t55 = t29+t54
      t56 = 1.0D0/t55
      t57 = ec22*t4*t48
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
      t68 = 1.0D0/t55**2
      t69 = ec12*t4*t46*t47*4.0D0
      t70 = 1.0D0/t50**2
      t74 = ec12*t4*t13*t46*t48*4.0D0
      t71 = t69-t74
      t72 = ec12*t4*t13*t17*t48*t51*4.0D0
      t73 = ec12*t4*t39*t44*t56*4.0D0
      t75 = t10*t13*t17*t70*t71
      t103 = ec12*t4*t17*t47*t51*4.0D0
      t76 = t72+t75-t103
      t104 = eta*t30*t44*t45*t60*t68*t76*2.0D0
      t77 = t73-t104
      t78 = ec22*lch*t39*t77*tbar_2
      t79 = t2**2
      t80 = t3**2
      t81 = 1.0D0/ec22
      t82 = 1.0D0/lch
      t83 = t23-t24
      t84 = lch*t39*t66*tbar_2
      t85 = t18*t39*t44*t56*t64
      t86 = t10*t46*t47*t58
      t96 = ec22*t4*t13*t46*t48
      t87 = t86-t96
      t88 = t10*t13*t17*t70*t87
      t89 = ec22*t4*t13*t17*t48*t51
      t97 = t10*t17*t47*t51*t58
      t90 = t88+t89-t97
      t98 = eta*t30*t44*t45*t60*t68*t90*2.0D0
      t91 = t85-t98
      t92 = e23*t2*t3*2.0D0
      t93 = ec22+t92-e22*t79-e33*t80
      t94 = estiff22*t93
      t95 = t94-ec22*lch*t39*t66*tbar_2
      t99 = ec22*lch*t18*t64*t66*t67*tbar_2
      t100 = t79-t80
      t101 = 1.0D0/t28
      t102 = ec22*ec23*lch*t4*t66*t67*tbar_2*4.0D0
      t105 = C22*t79
      t106 = C23*t80
      t107 = t105+t106
      t108 = t2*t107
      t109 = C55*t2*t80*2.0D0
      t110 = t108+t109
      t111 = C23*t79
      t112 = C33*t80
      t113 = t111+t112
      t114 = t3*t113
      t115 = C55*t3*t79*2.0D0
      t116 = t114+t115
      t117 = C23*t2*t3
      t118 = t117-C33*t2*t3
      t119 = t3*t118
      t120 = t119-C55*t2*t100
      t121 = C22*t2*t3
      A0(1,1) = -C44*t79-C66*t80+lch*t39*tbar_1*(t38*t44*t45-1.0D0)*2.0D
     &0+ec12*lch*t39*tbar_1*(ec12*t4*t38*t39*t44*4.0D0-eta*t30*1.0D0/t37
     &**2*t44*t45*t60*(t72+t10*t13*t17*1.0D0/t32**2*t71-ec12*t4*t17*t47*
     &t51*4.0D0)*2.0D0)*2.0D0-lch*t4*t5*t66*t67*tbar_1*8.0D0
      A0(1,2) = ec12*lch*t39*t91*tbar_1*2.0D0-ec12*lch*t18*t64*t66*t67*t
     &bar_1*2.0D0
      A0(1,3) = ec12*ec23*lch*t4*t66*t67*tbar_1*(-8.0D0)+ec12*ec23*lch*t
     &4*t44*t56*t101*tbar_1*8.0D0
      A0(2,1) = t78-t81*t82*t83*(t78-ec12*ec22*lch*t4*t66*t67*tbar_2*4.0
     &D0)-ec12*ec22*lch*t4*t66*t67*tbar_2*4.0D0
      A0(2,2) = t84-t99-t2*t110+t3*t116+t81*t82*t83*(estiff22-t84+t99-ec
     &22*lch*t39*t91*tbar_2)-t81*t82*t95*(t61-t63)-1.0D0/ec22**2*t82*t83
     &*t95+ec22*lch*t39*t91*tbar_2
      A0(2,3) = -t102-t3*t120+t2*(t2*(t121-C23*t2*t3)-C55*t3*t100)+t81*t
     &82*t83*(t102-ec22*ec23*lch*t4*t44*t56*t101*tbar_2*4.0D0)+ec22*ec23
     &*lch*t4*t44*t56*t101*tbar_2*4.0D0
      A0(3,1) = ec23*lch*t39*t77*tbar_3*2.0D0-ec12*ec23*lch*t4*t66*t67*t
     &bar_3*8.0D0
      A0(3,2) = -t3*t110-t2*t116+ec23*lch*t39*t91*tbar_3*2.0D0-ec23*lch*
     &t18*t64*t66*t67*tbar_3*2.0D0
      A0(3,3) = -t3*(t2*(t117-t121)+C55*t3*t100)+t2*t120+lch*t39*t66*tba
     &r_3*2.0D0-lch*t4*t26*t66*t67*tbar_3*8.0D0+lch*t4*t26*t44*t56*t101*
     &tbar_3*8.0D0

program jacobian_test

  use jacobian_f77

  implicit none
  integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: dp = selected_real_kind(15, 307)
  real*8, dimension(6,6) :: C
  real*16, dimension(6) :: eTotal_cr, ec_cr
  real*16, dimension(3) :: tbar_cr
  real*16, dimension(3,3) :: result
  real*8 :: e22, x, alpha, lch, GIIc, GIc, eta

  C = 0.0d0
  C(1,1) = 83.4654904587297d0
  C(1,2) = 3.11992755363710d0
  C(1,3) = 3.11992755363710d0
  C(2,2) = 9.99982098790910d0
  C(2,3) = 4.76172574981386d0
  C(3,3) = 9.99982098790910d0
  C(4,4) = 6.47000000000000d0
  C(5,5) = 3.5d0
  C(6,6) = 6.47d0
  eta = 1.45d0
  alpha = 1.37881010907552d0
  lch = 1.00502458838063d0
  GIc = 0.1d0
  GIIc = 3.0d-4
  e22 = 7.7d0
  eTotal_cr(1) = 2.604713340655021d-002
  eTotal_cr(2) = -5.505555810078023d-003
  eTotal_cr(3) = -5.505555810078023d-003
  eTotal_cr(4) = -4.235821198543642d-022
  eTotal_cr(5) = 1.045701209432247d-021
  eTotal_cr(6) = 1.820553897296562d-021
  ec_cr = 0.0d0
  ec_cr(2) = 1.075982801850006d-007
  ec_cr(4) = -5.302694461785174d-022
  ec_cr(5) = 1.083949357049625d-018
  tbar_cr(1) = 4.184479366209085d-021
  tbar_cr(2) = 3.430434285843303d-002
  tbar_cr(3) = 2.705608424223070d-021

  result = jacobian(C,ec_cr,tbar_cr,lch,alpha,eta,GIc,GIIc,e22,eTotal_cr)

  call printarray(result)

    contains

      subroutine printarray(array)
        implicit none
        real*16, dimension(3,3), intent(in) :: array
        integer :: i
        do i = 1, ubound(array, 1)
          print *, array(i, :)
        end do
      end subroutine printarray

end program jacobian_test


program camanho test
    implicit none
    REAL*8 e11,e22,e33,nu12,nu13,nu23,g12,g13,g23 ! elastic constants
    REAL*8 nu21,nu31,nu32,delta,alpha0
    
    ! Elastic constants orthotropic ply

    e11 = 146.8
    e22 = 11.6
    e33 = 11.6
    nu12 = 0.32
    nu13 = 0.32
    nu23 = 0.45
    g12 = 6.46
    g13 = 6.46
    g23 = 4.38

    ! Fracture Angle

    alpha0 = 53.0*0.017453292519943295 !converts to radians

    ! Fracture toughness

    G1plus = props(16)
    G1minus = props(17)
    G2plus = props(18)
    G2minus = props(19)
    G6 = props(20)

    ! Stiffness matrix orthotropic material

    nu21 = nu12*(e22/e11)
    nu31 = nu13*(e33/e11)
    nu32 = nu23*(e33/e22)

    delta = 1.0d0/(e22**2.0d0*nu12**2.0d0 +
    1        2.0d0*e33*e22*nu12*nu13*nu23 +
    1        e33*e22*nu13**2.0d0 - e11*e22 + e11*e33*nu23**2.0d0)

    C_voigt = 0.0d0 ! set all elements in array equal to zero

    C_voigt(1,1) = -(e11**2.0d0*(-e33*nu23**2.0d0 + e22))*delta
    C_voigt(2,1) = -(e11*e22*(e22*nu12 + e33*nu13*nu23))*delta
    C_voigt(3,1) = -(e11*e22*e33*(nu13 + nu12*nu23))*delta
    C_voigt(1,2) = -(e11*e22*(e22*nu12 + e33*nu13*nu23))*delta
    C_voigt(2,2) = -(e22**2.0d0*(-e33*nu13**2.0d0 + e11))*delta
    C_voigt(3,2) = -(e22*e33*(e11*nu23 + e22*nu12*nu13))*delta
    C_voigt(1,3) = -(e11*e22*e33*(nu13 + nu12*nu23))*delta
    C_voigt(2,3) = -(e22*e33*(e11*nu23 + e22*nu12*nu13))*delta
    C_voigt(3,3) = -(e22*e33*(-e22*nu12**2.0d0 + e11))*delta
    C_voigt(4,4) = g12
    C_voigt(5,5) = g23
    C_voigt(6,6) = g13 

    result = ddot(A,B)

    print *, 'Result =', result
    pause 'Press RETURN to exit'

    contains

        real*8 function ddot(A,B)
            implicit none
            real*8, intent(in) :: A(:, :), B(:, :)
            integer i,j,dim1,dim2
            dim1 = UBOUND(A,DIM = 1)  !Ascertain the bounds of the incoming arrays.
            dim2 = UBOUND(A,DIM = 2)

            ddot = 0.0d0
            do i = 1,dim1
                do j = 1,dim2
                    ddot = ddot + A(i,j) * B(i,j)
                    print *, ddot
                end do
            end do
        end function ddot

end program matmul_test



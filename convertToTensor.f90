program tensorProduct
    implicit none
    real*8, dimension(3) :: n1, n2
    real*8, dimension(3,3) :: result

    n1(1) = 1.0d0
    n1(2) = 2.0d0
    n1(3) = 3.0d0

    n2(1) = 4.0d0
    n2(2) = 5.0d0
    n2(3) = 6.0d0

    result = tensor_product(n1,n2)

    print *, 'Result =', result
    pause 'Press RETURN to exit'

    contains

      FUNCTION convertToTensor(A)
        IMPLICIT NONE
        ! input variables 
        REAL*8, DIMENSION(6,6) :: A
        ! local variables
        INTEGER :: i,j,k,l,s1,s2,idx1,idx2
        REAL*8 :: term
        ! output variables
        REAL*8, DIMENSION(3,3,3,3) :: convertToTensor

        DO i = 1,3
          DO j = 1,3
            DO k = 1,3
              DO l = 1,3
                IF (i.eq.j) THEN
                  IF (i.eq.1) THEN
                    idx1 = 1
                  ELSE IF (i.eq.2) THEN
                    idx1 = 2
                  ELSE IF (i.eq.3) THEN
                    idx1 = 3
                  END IF
                ELSE
                    s1 = i + j
                  IF (s1.eq.5) THEN
                    idx1 = 4
                  ELSE IF (s1.eq.4) THEN
                    idx1 = 5
                  ELSE IF (s1.eq.3) THEN
                    idx1 = 6
                  END IF
                END IF
                IF (k.eq.l) THEN
                  IF (k.eq.1) THEN
                    idx2 = 1
                  ELSE IF (k.eq.2) THEN
                    idx2 = 2
                  ELSE IF (k.eq.3) THEN
                    idx2 = 3
                  END IF
                ELSE
                    s2 = k + l
                  IF (s2.eq.5) THEN
                    idx2 = 4
                  ELSE IF (s2.eq.4) THEN
                    idx2 = 5
                  ELSE IF (s2.eq.3) THEN
                    idx2 = 6
                  END IF
                END IF

                term = A(idx1,idx2)
                convertToTensor(i,j,k,l) = term
              END DO
            END DO
          END DO
        END DO
        RETURN
      END FUNCTION convertToTensor

end program tensorProduct



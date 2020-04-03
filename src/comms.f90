!------------------------------------------------------------------------------
MODULE comms
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP), PARAMETER :: zero = 0.d0
  COMPLEX(DP), PARAMETER :: cone = (1.d0, 0.d0)
  COMPLEX(DP), PARAMETER :: czero = (0.d0, 0.d0)
  COMPLEX(DP), PARAMETER :: ci = (0.d0, 1.d0)
  REAL(DP), PARAMETER :: PI = 3.1415926535897932384626433832795028841971694_dp
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE io_error(error_msg)
  !----------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: error_msg
    WRITE(*,*) error_msg
    STOP
  END SUBROUTINE io_error
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE utility_zgemm(a, b, c, beta_opt, transa_opt, transb_opt)
  !----------------------------------------------------------------------------
  !!
  !! Wrapper for LAPACK ZGEMM subroutine
  !!
  !! Return matrix product of complex matrices a and b: C = beta * Op(A) Op(B)
  !!
  !! Default: beta = (0.d0, 0.d0), Op(A) = A, Op(B) = B
  !!
  !! transa = 'N'  ==> Op(X) = X
  !! transa = 'T'  ==> Op(X) = transpose(X)
  !! transa = 'C'  ==> Op(X) = congj(transpose(X))
  !!
  !! Adapted from utility_zgemm_new in wannier90/src/utility.F90
  !!
  !----------------------------------------------------------------------------
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: a(:, :)
    COMPLEX(DP), INTENT(IN) :: b(:, :)
    COMPLEX(DP), INTENT(OUT) :: c(:, :)
    COMPLEX(DP), INTENT(IN), OPTIONAL :: beta_opt
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL  :: transa_opt
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL  :: transb_opt
    !
    INTEGER :: m
    !! number of rows in Op(A) and C
    INTEGER :: n
    !! number of columns in Op(B) and C
    INTEGER :: k
    !! number of columns in Op(A) resp. rows in Op(B)
    COMPLEX(DP) :: beta
    CHARACTER(LEN=1) :: transa, transb
    !
    transa = 'N'
    transb = 'N'
    beta = czero
    IF (PRESENT(transa_opt)) transa = transa_opt
    IF (PRESENT(transb_opt)) transb = transb_opt
    IF (PRESENT(beta_opt)) beta = beta_opt
    !
    m = SIZE(c, 1)
    n = SIZE(c, 2)
    !
    IF (transa /= 'N') THEN
      k = SIZE(a, 1)
    ELSE
      k = SIZE(a, 2)
    ENDIF
    !
    CALL ZGEMM(transa, transb, m, n, k, &
      cone, a, SIZE(a, 1), b, size(b, 1), beta, c, m)
    !
  END SUBROUTINE utility_zgemm
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE inv_omega_minus_mat(ndim, mat_in, omega, mat_out, flag)
  !----------------------------------------------------------------------------
  !! Compute mat_out = (omega - mat_in)^-1
  !----------------------------------------------------------------------------
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ndim
    COMPLEX(DP), INTENT(IN) :: omega
    COMPLEX(DP), INTENT(IN) :: mat_in(ndim, ndim)
    COMPLEX(DP), INTENT(OUT) :: mat_out(ndim, ndim)
    CHARACTER(LEN=*), INTENT(IN) :: flag
    !! Description to be included in the error message when inversion fails
!f2py depend(ndim) :: mat_in, mat_out
    !
    COMPLEX(DP), ALLOCATABLE :: mat_temp(:, :)
    INTEGER, ALLOCATABLE :: ipiv(:)
    INTEGER :: i
    INTEGER :: info
    !
    ! Set mat_temp = omega - mat_in
    ALLOCATE(mat_temp(ndim, ndim))
    mat_temp = - mat_in
    mat_out = czero
    DO i = 1, ndim
        mat_temp(i, i) = omega + mat_temp(i, i)
        mat_out(i, i) = cone
    END DO
    !
    ! Invert matrix
    ALLOCATE(ipiv(ndim))
    CALL ZGESV(ndim,ndim, mat_temp, ndim, ipiv, mat_out, ndim, info)
    IF (info /= 0) CALL io_error('Error in ZGESV: &
        &matinv inversion in ' // flag)
    !
    DEALLOCATE(mat_temp)
    DEALLOCATE(ipiv)
    !
  END SUBROUTINE inv_omega_minus_mat
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE util_sum_abs(arr, n1, n2, out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, n2
    COMPLEX(DP), INTENT(IN) :: arr(n1, n2)
    REAL(DP), INTENT(OUT) :: out
!f2py intent(hide) :: n1, n2
    out = SUM(ABS(arr))
  END SUBROUTINE util_sum_abs
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE util_trace(arr, n, out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    COMPLEX(DP), INTENT(IN) :: arr(n, n)
    COMPLEX(DP), INTENT(OUT) :: out
!f2py intent(hide) :: n
    INTEGER :: i
    out = czero
    DO i = 1, n
      out = out + arr(i, i)
    ENDDO
  END SUBROUTINE util_trace
  !----------------------------------------------------------------------------
END MODULE comms
!------------------------------------------------------------------------------

SUBROUTINE get_hk(kvec, nw, nrpts, hr, rvec, ndegen, hk)
  !
  USE kinds, ONLY : DP
  USE comms, ONLY : ci, PI
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: kvec(3)
  INTEGER, INTENT(IN) :: nw
  INTEGER, INTENT(IN) :: nrpts
  COMPLEX(DP), INTENT(IN) :: hr(nrpts, nw, nw)
  INTEGER, INTENT(IN) :: rvec(3, nrpts)
  INTEGER, INTENT(IN) :: ndegen(nrpts)
  COMPLEX(DP), INTENT(OUT) :: hk(nw, nw)
!f2py depend(nw) :: hk, hr
!f2py depend(nrpts) :: hk, hr, rvec, ndegen
  !
  INTEGER :: ir
  REAL(DP) :: herm_check, rdotk
  COMPLEX(DP) :: coeff
  !
  DO ir = 1, nrpts
    rdotk = SUM(kvec * REAL(rvec(:, ir), DP))
    coeff = EXP((0.d0, 2.d0) * PI * rdotk) / REAL(ndegen(ir), DP)
    hk = hk + hr(ir, :, :) * coeff
  ENDDO
  !
  herm_check = SQRT(SUM(ABS(hk - TRANSPOSE(CONJG(hk))) ** 2))
  IF (herm_check > 1.d-6) THEN
    print*, "WARNING: Hermiticity error: ", herm_check
  ENDIF
  !
END SUBROUTINE get_hk

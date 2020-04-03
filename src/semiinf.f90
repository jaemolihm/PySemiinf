SUBROUTINE update_ideal(nbulk, omega, e_b, e_s, a_b, b_b)
!----------------------------------------------------------------------------
!! use recurrence relations to update matrices
!----------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE comms, ONLY : cone, inv_omega_minus_mat, utility_zgemm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nbulk
  COMPLEX(DP), INTENT(IN) :: omega
  COMPLEX(DP), INTENT(INOUT) :: e_b(nbulk, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: e_s(nbulk, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: a_b(nbulk, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: b_b(nbulk, nbulk)
!f2py depend(nbulk) :: e_b, e_s, a_b, b_b
!f2py intent(in, out) :: e_b, e_s, a_b, b_b
  !
  COMPLEX(DP) :: inv_e_b(nbulk, nbulk)
  COMPLEX(DP) :: a_b_prev(nbulk, nbulk)
  COMPLEX(DP) :: b_b_prev(nbulk, nbulk)
  COMPLEX(DP) :: temp_mat(nbulk, nbulk)
  !
  ! inv_e_b = (omgea - e_b)^-1
  CALL inv_omega_minus_mat(nbulk, e_b, omega, inv_e_b, &
      'update_ideal for inv_e_b')
  !
  CALL ZCOPY(nbulk * nbulk, a_b, 1, a_b_prev, 1)
  CALL ZCOPY(nbulk * nbulk, b_b, 1, b_b_prev, 1)
  !
  CALL utility_zgemm(a_b_prev, inv_e_b, temp_mat)
  CALL utility_zgemm(temp_mat, b_b_prev, e_s, cone)
  CALL utility_zgemm(temp_mat, b_b_prev, e_b, cone)
  CALL utility_zgemm(temp_mat, a_b_prev, a_b)
  !
  CALL utility_zgemm(b_b_prev, inv_e_b, temp_mat)
  CALL utility_zgemm(temp_mat, b_b_prev, e_b, cone)
  CALL utility_zgemm(temp_mat, b_b_prev, b_b)
!------------------------------------------------------------------------------
END SUBROUTINE update_ideal
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE update_nonideal(nsurf, nbulk, omega, e_s0, e_s1, e_b, a_s, a_b, &
    b_s, b_b)
!------------------------------------------------------------------------------
!! use recurrence relation to update matrices
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE comms, ONLY : cone, inv_omega_minus_mat, utility_zgemm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nsurf
  INTEGER, INTENT(IN) :: nbulk
  COMPLEX(DP), INTENT(IN) :: omega
  COMPLEX(DP), INTENT(INOUT) :: e_s0(nsurf, nsurf)
  COMPLEX(DP), INTENT(INOUT) :: e_s1(nbulk, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: e_b(nbulk, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: a_s(nsurf, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: a_b(nbulk, nbulk)
  COMPLEX(DP), INTENT(INOUT) :: b_s(nbulk, nsurf)
  COMPLEX(DP), INTENT(INOUT) :: b_b(nbulk, nbulk)
!f2py depend(nsurf) :: e_s0, a_s, b_s
!f2py depend(nbulk) :: e_s1, e_b, a_s, a_b, b_s, b_b
!f2py intent(in, out) :: e_s0, e_s1, e_b, a_s, a_b, b_s, b_b
  !
  COMPLEX(DP) :: inv_e_b(nbulk, nbulk)
  COMPLEX(DP) :: a_s_prev(nsurf, nbulk)
  COMPLEX(DP) :: b_s_prev(nbulk, nsurf)
  COMPLEX(DP) :: a_b_prev(nbulk, nbulk)
  COMPLEX(DP) :: b_b_prev(nbulk, nbulk)
  COMPLEX(DP) :: temp_mat(nbulk, nbulk)
  COMPLEX(DP) :: temp_mat_sb(nsurf, nbulk)
  !
  ! inv_e_b = (omgea - e_b)^-1
  CALL inv_omega_minus_mat(nbulk, e_b, omega, inv_e_b, &
      'slab_update for inv_e_b')
  !
  CALL ZCOPY(nsurf*nbulk, a_s, 1, a_s_prev, 1)
  CALL ZCOPY(nbulk*nsurf, b_s, 1, b_s_prev, 1)
  CALL ZCOPY(nbulk*nbulk, a_b, 1, a_b_prev, 1)
  CALL ZCOPY(nbulk*nbulk, b_b, 1, b_b_prev, 1)
  !
  CALL utility_zgemm(inv_e_b, a_b, temp_mat)
  CALL utility_zgemm(a_s_prev, temp_mat, a_s)
  CALL utility_zgemm(a_b_prev, temp_mat, a_b)
  CALL utility_zgemm(b_b_prev, temp_mat, e_b, cone)
  !
  CALL utility_zgemm(b_b_prev, inv_e_b, temp_mat)
  CALL utility_zgemm(temp_mat, b_s_prev, b_s)
  CALL utility_zgemm(temp_mat, b_b_prev, b_b)
  !
  CALL utility_zgemm(inv_e_b, b_b_prev, temp_mat)
  CALL utility_zgemm(a_b_prev, temp_mat, e_b, cone)
  CALL utility_zgemm(a_b_prev, temp_mat, e_s1, cone)
  !
  CALL utility_zgemm(a_s_prev, inv_e_b, temp_mat_sb)
  CALL utility_zgemm(temp_mat_sb, b_s_prev, e_s0, cone)
  !
!------------------------------------------------------------------------
END SUBROUTINE update_nonideal
!------------------------------------------------------------------------

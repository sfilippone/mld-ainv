subroutine mld_d_ainvk_solver_check(sv,info)
  

  use psb_base_mod
  use mld_d_ainvk_solver, mld_protect_name => mld_d_ainvk_solver_check

  Implicit None

  ! Arguments
  class(mld_d_ainvk_solver_type), intent(inout) :: sv
  integer, intent(out)                   :: info
  Integer           :: err_act
  character(len=20) :: name='d_ainvk_solver_check'

  call psb_erractionsave(err_act)
  info = psb_success_

  call mld_check_def(sv%fill_in,&
       & 'Level',0,is_legal_ml_lev)
  call mld_check_def(sv%inv_fill,&
       & 'Level',0,is_legal_ml_lev)
  call mld_check_def(sv%thresh,&
       & 'Eps',dzero,is_legal_d_fact_thrs)

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_ainvk_solver_check

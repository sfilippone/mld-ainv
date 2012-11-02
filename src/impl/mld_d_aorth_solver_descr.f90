subroutine mld_d_ainv_solver_descr(sv,info,iout,coarse)
  

  use psb_base_mod
  use mld_d_ainv_solver, mld_protect_name => mld_d_ainv_solver_descr

  Implicit None

  ! Arguments
  class(mld_d_ainv_solver_type), intent(in) :: sv
  integer, intent(out)                     :: info
  integer, intent(in), optional            :: iout
  logical, intent(in), optional       :: coarse

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_d_ainv_solver_descr'
  integer :: iout_

  call psb_erractionsave(err_act)
  info = psb_success_
  if (present(iout)) then 
    iout_ = iout 
  else
    iout_ = 6
  endif

  write(iout_,*) '  AINV: Approximate Inverse with sparse boconjugation '
  write(iout_,*) '  Algoritm variant       :',sv%alg    
  write(iout_,*) '  Fill level             :',sv%fill_in
  write(iout_,*) '  Fill threshold         :',sv%thresh

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_ainv_solver_descr

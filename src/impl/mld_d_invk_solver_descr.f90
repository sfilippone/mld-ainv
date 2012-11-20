subroutine mld_d_invk_solver_descr(sv,info,iout,coarse)
  

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_d_invk_solver_descr

  Implicit None

  ! Arguments
  class(mld_d_invk_solver_type), intent(in) :: sv
  integer, intent(out)                     :: info
  integer, intent(in), optional            :: iout
  logical, intent(in), optional       :: coarse

  ! Local variables
  integer      :: err_act
  integer      :: ictxt, me, np
  character(len=20), parameter :: name='mld_d_invk_solver_descr'
  integer :: iout_

  call psb_erractionsave(err_act)
  info = psb_success_
  if (present(iout)) then 
    iout_ = iout 
  else
    iout_ = 6
  endif

  write(iout_,*) '  INVK Approximate Inverse with ILU(N) '
  write(iout_,*) '  Fill level             :',sv%fill_in
  write(iout_,*) '  Inverse fill level     :',sv%inv_fill

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_invk_solver_descr
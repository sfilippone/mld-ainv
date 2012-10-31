subroutine mld_d_invk_solver_seti(sv,what,val,info)
  

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_d_invk_solver_seti

  Implicit None

  ! Arguments
  class(mld_d_invk_solver_type), intent(inout) :: sv 
  integer, intent(in)                    :: what 
  integer, intent(in)                    :: val
  integer, intent(out)                   :: info
  Integer :: err_act
  character(len=20)  :: name='d_invk_solver_seti'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(what) 
  case(mld_sub_fillin_)
    sv%fill_in   = val
  case(mld_inv_fillin_)
    sv%inv_fill  = val
  case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_invk_solver_seti

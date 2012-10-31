subroutine mld_d_ainvk_solver_setr(sv,what,val,info)
  

  use psb_base_mod
  use mld_d_ainvk_solver, mld_protect_name => mld_d_ainvk_solver_setr

  Implicit None

  ! Arguments
  class(mld_d_ainvk_solver_type), intent(inout) :: sv 
  integer, intent(in)                    :: what 
  real(psb_dpk_), intent(in)             :: val
  integer, intent(out)                   :: info
  Integer :: err_act
  character(len=20)  :: name='d_ainvk_solver_setr'

  call psb_erractionsave(err_act)
  info = psb_success_

  select case(what)
  case(mld_sub_iluthrs_) 
    sv%thresh = val
  case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
!!$      goto 9999
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
end subroutine mld_d_ainvk_solver_setr

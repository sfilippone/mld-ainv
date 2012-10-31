subroutine mld_d_ainvk_solver_free(sv,info)
  

  use psb_base_mod
  use mld_d_ainvk_solver, mld_protect_name => mld_d_ainvk_solver_free

  Implicit None

  ! Arguments
  class(mld_d_ainvk_solver_type), intent(inout) :: sv
  integer, intent(out)                       :: info
  Integer :: err_act
  character(len=20)  :: name='d_ainvk_solver_free'

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(sv%d)) then 
    deallocate(sv%d,stat=info)
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999 
    end if
  end if
  call sv%l%free()
  call sv%u%free()
  call sv%dv%free(info)
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_ainvk_solver_free

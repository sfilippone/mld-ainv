subroutine  psb_d_dsc_trim(a)
  
  use psb_realloc_mod
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_trim
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  Integer :: err_act, info, nz, n
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$  n   = a%get_ncols()
!!$  nz  = a%get_nzeros()
!!$  if (info == psb_success_) call psb_realloc(n+1,a%icp,info)
!!$  if (info == psb_success_) call psb_realloc(nz,a%ia,info)
!!$  if (info == psb_success_) call psb_realloc(nz,a%val,info)
!!$
!!$  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_trim

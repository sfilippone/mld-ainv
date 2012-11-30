subroutine  psb_d_dsc_reallocate_nz(nz,a) 
  
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_reallocate_nz
  implicit none 
  integer, intent(in) :: nz
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='d_dsc_reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  call psb_realloc(nz,a%ia,info)
!!$  if (info == psb_success_) call psb_realloc(nz,a%val,info)
!!$  if (info == psb_success_) call psb_realloc(max(nz,a%get_nrows()+1,a%get_ncols()+1),a%icp,info)
!!$  if (info /= psb_success_) then 
!!$    call psb_errpush(psb_err_alloc_dealloc_,name)
!!$    goto 9999
!!$  end if
!!$
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_reallocate_nz

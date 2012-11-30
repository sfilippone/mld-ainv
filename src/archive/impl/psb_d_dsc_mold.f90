subroutine psb_d_dsc_mold(a,b,info) 
  
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_mold
  use psb_error_mod
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer, intent(out)                    :: info
  Integer :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  allocate(psb_d_dsc_sparse_mat :: b, stat=info)

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_dsc_mold

subroutine psb_d_dsc_cp_from(a,b)
  
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_cp_from
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  type(psb_d_dsc_sparse_mat), intent(in)   :: b


  Integer :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_success_

  call a%allocate(b%get_nrows(),b%get_ncols(),b%get_nzeros())
  call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
  call psb_safe_ab_cpy( b%cols, a%cols , info)

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_dsc_cp_from

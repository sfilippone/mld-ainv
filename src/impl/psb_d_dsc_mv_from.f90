subroutine psb_d_dsc_mv_from(a,b)
  
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_mv_from
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout)  :: a
  type(psb_d_dsc_sparse_mat), intent(inout) :: b


  Integer :: err_act, info
  character(len=20)  :: name='mv_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
  call move_alloc(b%cols, a%cols)
  call b%free()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_dsc_mv_from

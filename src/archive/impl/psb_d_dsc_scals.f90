subroutine psb_d_dsc_scals(d,a,info) 
  
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_scals
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer, intent(out)            :: info

  Integer :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  do i=1,a%get_nzeros()
!!$    a%val(i) = a%val(i) * d
!!$  enddo
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

end subroutine psb_d_dsc_scals

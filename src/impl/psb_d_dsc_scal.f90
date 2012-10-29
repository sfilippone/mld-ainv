subroutine psb_d_dsc_scal(d,a,info,side) 
  
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_scal
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer, intent(out)            :: info
  character, intent(in), optional :: side

  Integer :: err_act,mnm, i, j, n
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  n = a%get_ncols()
!!$  if (size(d) < n) then 
!!$    info=psb_err_input_asize_invalid_i_
!!$    call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$  do i=1, n
!!$    do j = a%icp(i), a%icp(i+1) -1 
!!$      a%val(j) = a%val(j) * d(a%ia(j))
!!$    end do
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

end subroutine psb_d_dsc_scal

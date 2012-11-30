subroutine psb_d_dsc_get_diag(a,d,info) 
  
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_get_diag
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act, mnm, i, j, k
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999

  mnm = min(a%get_nrows(),a%get_ncols())
!!$  if (size(d) < mnm) then 
!!$    info=psb_err_input_asize_invalid_i_
!!$    call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  if (a%is_triangle().and.a%is_unit()) then 
!!$    d(1:mnm) = done 
!!$  else
!!$    do i=1, mnm
!!$      d(i) = dzero
!!$      do k=a%icp(i),a%icp(i+1)-1
!!$        j=a%ia(k)
!!$        if ((j == i) .and.(j <= mnm )) then 
!!$          d(i) = a%val(k)
!!$        endif
!!$      enddo
!!$    end do
!!$  endif
!!$  do i=mnm+1,size(d) 
!!$    d(i) = dzero
!!$  end do
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_get_diag

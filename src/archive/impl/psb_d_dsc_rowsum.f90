subroutine psb_d_dsc_rowsum(d,a) 
  
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_rowsum
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  m = a%get_ncols()
!!$  n = a%get_nrows()
!!$  if (size(d) < n) then 
!!$    info=psb_err_input_asize_small_i_
!!$    int_err(1) = 1
!!$    int_err(2) = size(d)
!!$    int_err(3) = n
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  d   = dzero
!!$
!!$  do i=1, m
!!$    do j=a%icp(i),a%icp(i+1)-1
!!$      k = a%ia(j)
!!$      d(k) = d(k) + (a%val(k))
!!$    end do
!!$  end do
!!$
  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_rowsum

function psb_d_dsc_csnmi(a) result(res)
  
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csnmi
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nr, ir, jc, nc, info
  real(psb_dpk_), allocatable  :: acc(:) 
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_csnmi'
  logical, parameter :: debug=.false.


  res = dzero 
!!$  nr = a%get_nrows()
!!$  nc = a%get_ncols()
!!$  allocate(acc(nr),stat=info)
!!$  if (info /= psb_success_) then 
!!$    return
!!$  end if
!!$  acc(:) = dzero
!!$  do i=1, nc
!!$    do j=a%icp(i),a%icp(i+1)-1  
!!$      acc(a%ia(j)) = acc(a%ia(j)) + abs(a%val(j))
!!$    end do
!!$  end do
!!$  do i=1, nr
!!$    res = max(res,acc(i))
!!$  end do
!!$  deallocate(acc)

end function psb_d_dsc_csnmi

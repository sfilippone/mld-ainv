function psb_d_dsc_csnm1(a) result(res)
  
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csnm1

  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_dsc_csnm1'
  logical, parameter :: debug=.false.


  res = dzero 
!!$  m = a%get_nrows()
!!$  n = a%get_ncols()
!!$  do j=1, n
!!$    acc = dzero 
!!$    do k=a%icp(j),a%icp(j+1)-1
!!$      acc = acc + abs(a%val(k))
!!$    end do
!!$    res = max(res,acc)
!!$  end do
!!$  
!!$  return
!!$
end function psb_d_dsc_csnm1

subroutine mld_d_base_ainv_update_a(sv,x,desc_data,info)
  use mld_d_base_ainv_mod, mld_protect_name => mld_d_base_ainv_update_a
  use psb_base_mod

  implicit none 

  type(psb_desc_type), intent(in)      :: desc_data
  class(mld_d_base_ainv_solver_type), intent(inout) :: sv
  real(psb_dpk_),intent(in)            :: x(:)
  integer, intent(out)                 :: info

  ! Local variables
  real(psb_dpk_), allocatable   :: dd(:), ee(:)
  integer(psb_ipk_) :: nrows, ncols, i, j, k, nzr, nzc, ir, ic
  
  type(psb_d_csc_sparse_mat) :: ac
  type(psb_d_csr_sparse_mat) :: ar
  
  dd = sv%dv%get_vect()
  if (size(x)<size(dd)) then 
    write(0,*) 'Wrong size in update',size(x),size(dd)
    info = -1
    return
  end if
  
  ee = x(1:size(dd))
  ! Now compute the update
  call sv%z%cp_to(ac)
  call sv%w%cp_to(ar)
  ! Compute diag(W ee Z)
  ! ee Z
  call ac%scal(ee,info,side='L')
  ! Compute diag(W ee Z)
  do i=1,ar%get_nrows()
    ir  = ar%irp(i)
    nzr = ar%irp(i+1)-ar%irp(i)
    ic  = ac%icp(i)
    nzc = ac%icp(i+1)-ac%icp(i)
    ee(i) = psb_spdot_srtd(nzr,ar%ja(ir:ir+nzr-1),ar%val(ir:ir+nzr-1),&
         & nzc,ac%ia(ic:ic+nzc-1),ac%val(ic:ic+nzc-1))    
  end do
  ! Need to invert
  dd = done/dd
  dd = dd + ee
  dd = done/dd
  sv%d = dd
  call sv%dv%bld(dd,mold=sv%dv%v)

end subroutine mld_d_base_ainv_update_a

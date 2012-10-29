subroutine psb_d_cp_dsc_from_coo(a,b,info) 
  
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_from_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  type(psb_d_coo_sparse_mat)   :: tmp
  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  ! This is to have fix_coo called behind the scenes
  call tmp%cp_from_coo(b,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

end subroutine psb_d_cp_dsc_from_coo

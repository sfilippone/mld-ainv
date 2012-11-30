subroutine psb_d_mv_dsc_to_coo(a,b,info) 
  
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_mv_dsc_to_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout)   :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  
  ! Too tired to figure out a smart mv right now. 
  call a%cp_to_coo(b,info)
  call a%free()

end subroutine psb_d_mv_dsc_to_coo

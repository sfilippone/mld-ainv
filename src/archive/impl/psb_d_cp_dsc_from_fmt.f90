subroutine psb_d_cp_dsc_from_fmt(a,b,info) 
  
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_from_fmt
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_d_dsc_sparse_mat) 
    call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
    call psb_safe_ab_cpy( b%cols, a%cols , info)

  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_cp_dsc_from_fmt

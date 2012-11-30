subroutine psb_d_cp_dsc_to_coo(a,b,info) 
  
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_to_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(in)  :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                      :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,k,irw, err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)

  if (allocated(a%cols)) then 
    k = 0
    do i=1, nc
      do j=1,a%cols(i)%nz
        k = k + 1
        b%ia(k)  = a%cols(i)%idx(j)
        b%ja(k)  = i
        b%val(k) = a%cols(i)%val(j)
      end do
    end do
  end if
  call b%set_nzeros(nza)
  call b%fix(info)

end subroutine psb_d_cp_dsc_to_coo

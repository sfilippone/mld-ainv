subroutine mld_d_invk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,ktrw,trw,info,sign,inlevs)
  

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_d_invk_copyin

  implicit none

  ! Arguments 
  type(psb_d_csr_sparse_mat), intent(in)    :: a
  type(psb_d_coo_sparse_mat), intent(inout) :: trw
  integer, intent(in)                  :: i,m,jmin,jmax
  integer, intent(inout)               :: ktrw,info
  integer, intent(inout)               :: rowlevs(:)
  real(psb_dpk_), intent(inout)        :: row(:)
  type(psb_int_heap), intent(inout)    :: heap
  real(psb_dpk_), optional, intent(in) :: sign
  integer, intent(in), optional        :: inlevs(:)

  ! Local variables
  integer             :: k,j,irb,err_act, nz
  integer, parameter  :: nrb=16
  real(psb_dpk_)      :: sign_
  character(len=20), parameter  :: name='invk_copyin'
  character(len=20)             :: ch_err

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  call psb_init_heap(heap,info) 
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_init_heap')
    goto 9999
  end if

  if (present(sign)) then 
    sign_ = sign
  else
    sign_ = done
  end if


  !
  ! Take a fast shortcut if the matrix is stored in CSR format
  !
  if (present(inlevs)) then 
    do j = a%irp(i), a%irp(i+1) - 1
      k          = a%ja(j)
      if ((jmin<=k).and.(k<=jmax)) then 
        row(k)     = sign_ * a%val(j)
        rowlevs(k) = inlevs(j)
        call psb_insert_heap(k,heap,info)
      end if
    end do
  else
    do j = a%irp(i), a%irp(i+1) - 1
      k          = a%ja(j)
      if ((jmin<=k).and.(k<=jmax)) then 
        row(k)     = sign_ * a%val(j)
        rowlevs(k) = 0
        call psb_insert_heap(k,heap,info)
      end if
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_d_invk_copyin

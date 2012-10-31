subroutine mld_dinvt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
  
     & irwt,ktrw,trw,info,sign)
  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_dinvt_copyin
  implicit none 
  type(psb_d_csr_sparse_mat), intent(in)    :: a
  type(psb_d_coo_sparse_mat), intent(inout) :: trw
  integer, intent(in)                  :: i, m,jmin,jmax,jd
  integer, intent(inout)               :: ktrw,nlw,nup,jmaxup,info
  integer, intent(inout)               :: irwt(:)
  real(psb_dpk_), intent(inout)        :: nrmi,row(:)
  type(psb_int_heap), intent(inout)    :: heap
  real(psb_dpk_), intent(in), optional :: sign

  integer                     :: k,j,irb,kin,nz, err_act
  integer, parameter          :: nrb=16
  real(psb_dpk_)              :: dmaxup, sign_
  real(psb_dpk_), external    :: dnrm2
  character(len=20), parameter  :: name='invt_copyin'

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)

  call psb_init_heap(heap,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_init_heap')
    goto 9999
  end if
  sign_ = done
  if (present(sign)) sign_ = sign
  !
  ! nrmi is the norm of the current sparse row (for the time being,
  ! we use the 2-norm).
  ! NOTE: the 2-norm below includes also elements that are outside
  ! [jmin:jmax] strictly. Is this really important? TO BE CHECKED.
  !

  nlw    = 0
  nup    = 0
  jmaxup = 0
  dmaxup = dzero
  nrmi   = dzero

  do j = a%irp(i), a%irp(i+1) - 1
    k = a%ja(j)
    if ((jmin<=k).and.(k<=jmax)) then 
      row(k)     = sign_ * a%val(j)
      call psb_insert_heap(k,heap,info)
      irwt(k) = 1
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if
    end if
    if (k<jd) nlw = nlw + 1 
    if (k>jd) then 
      nup = nup + 1
      if (abs(row(k))>dmaxup) then 
        jmaxup = k
        dmaxup = abs(row(k))
      end if
    end if
  end do
  nz   = a%irp(i+1) - a%irp(i)
  nrmi = dnrm2(nz,a%val(a%irp(i):),ione)

  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_dinvt_copyin

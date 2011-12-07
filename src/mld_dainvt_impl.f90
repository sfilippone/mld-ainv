
subroutine mld_dinvt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
     & irwt,ktrw,trw,info,sign)
  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dinvt_copyin
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

subroutine mld_dinvt(thres,i,nrmi,row,heap,irwt,ja,irp,val,nidx,idxs,info)

  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dinvt

  implicit none 

  ! Arguments
  type(psb_int_heap), intent(inout)   :: heap 
  integer, intent(in)                 :: i
  integer, intent(inout)              :: nidx,info
  integer, intent(inout)              :: irwt(:) 
  real(psb_dpk_), intent(in)          :: thres,nrmi
  integer, allocatable, intent(inout) :: idxs(:)
  integer, intent(in)                 :: ja(:),irp(:)
  real(psb_dpk_), intent(in)          :: val(:)
  real(psb_dpk_), intent(inout)       :: row(:)

  ! Local Variables
  integer               :: k,j,jj,lastk,iret
  real(psb_dpk_)      :: rwk

  info  = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i
  irwt(i) = 1 

  !
  ! Do while there are indices to be processed
  !
  do

    call psb_heap_get_first(k,heap,iret) 
    if (iret < 0) exit

    ! 
    ! An index may have been put on the heap more than once.
    !
    if (k == lastk) cycle

    lastk = k 

    !
    ! Dropping rule based on the threshold: compare the absolute
    ! value of each updated entry of row with thres * 2-norm of row.
    !
    rwk    = row(k)

    if (abs(rwk) < thres*nrmi) then
      ! 
      ! Drop the entry.
      !
      row(k) = dzero
      cycle
    else
      !
      ! Note: since U is scaled while copying it out (see ilut_copyout),
      ! we can use rwk in the update below.
      !           
      do jj=irp(k),irp(k+1)-1
        j = ja(jj)
        if (j<=k) then 
          info = -i 
          return
        endif
        !
        ! Update row(j) and, if it is not to be discarded, insert
        ! its index into the heap for further processing.
        !
        row(j)     = row(j) - rwk * val(jj)
        if (abs(row(j)) < thres*nrmi) then
          ! 
          ! Drop the entry.
          !
          row(j) = dzero
        else
          !
          ! Do the insertion.
          !
          if (irwt(j) == 0) then 
            call psb_insert_heap(j,heap,info)
            if (info /= psb_success_) return
            irwt(j) = 1
          end if
        endif
      end do
    end if

    !
    ! If we get here it is an index we need to keep on copyout.
    !
    nidx       = nidx + 1
    call psb_ensure_size(nidx,idxs,info,addsz=psb_heap_resize)      
    if (info /= psb_success_) return
    idxs(nidx) = k
    irwt(k)    = 0
  end do
  irwt(i) = 0

end subroutine mld_dinvt

subroutine mld_dinvt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
     & nidx,idxs,l2,ja,irp,val,info)

  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dinvt_copyout

  implicit none 

  ! Arguments
  integer, intent(in)                       :: fill_in,i,m,nidx,nlw,nup,jmaxup
  integer, intent(in)                       :: idxs(:)
  integer, intent(inout)                    :: l2, info
  integer, allocatable, intent(inout)       :: ja(:),irp(:)
  real(psb_dpk_), intent(in)                :: thres,nrmi
  real(psb_dpk_),allocatable, intent(inout) :: val(:)
  real(psb_dpk_), intent(inout)             :: row(:)

  ! Local variables
  real(psb_dpk_),allocatable   :: xw(:)
  integer, allocatable         :: xwid(:), indx(:)
  real(psb_dpk_)               :: witem, wmin
  integer                      :: widx
  integer                      :: k,isz,err_act,int_err(5),idxp, nz
  type(psb_double_idx_heap)    :: heap
  character(len=20), parameter :: name='invt_copyout'
  character(len=20)            :: ch_err
  logical                      :: fndmaxup

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)

  !
  ! Here we need to apply also the dropping rule base on the fill-in. 
  ! We do it by putting into a heap the elements that are not dropped
  ! by using the 2-norm rule, and then copying them out. 
  !
  ! The heap goes down on the entry absolute value, so the first item
  ! is the largest absolute value. 
  !

  call psb_init_heap(heap,info,dir=psb_asort_down_)

  if (info == psb_success_) allocate(xwid(nidx),xw(nidx),indx(nidx),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/3*nidx,0,0,0,0/),&
         & a_err='real(psb_dpk_)')
    goto 9999      
  end if

  !
  ! First the lower part
  !

  nz   = 0
  idxp = 0

  do

    idxp = idxp + 1
    if (idxp > nidx) exit
    if (idxs(idxp) >= i) exit
    widx      = idxs(idxp)
    witem     = row(widx)
    !
    ! Dropping rule based on the 2-norm
    !
    if (abs(witem) < thres*nrmi) cycle

    nz       = nz + 1 
    xw(nz)   = witem 
    xwid(nz) = widx
    call psb_insert_heap(witem,widx,heap,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_insert_heap')
      goto 9999
    end if
  end do

  if (nz > 1) then 
    write(psb_err_unit,*) 'Warning: lower triangle from ainvt???? '
  end if


  if (idxp <= size(idxs)) then 
    if (idxs(idxp) < i) then 
      do 
        idxp = idxp + 1
        if (idxp > nidx) exit
        if (idxs(idxp) >= i) exit
      end do
    end if
  end if
  idxp = idxp - 1 
  nz   = 0
  wmin=HUGE(wmin)
  do

    idxp = idxp + 1
    if (idxp > nidx) exit
    widx      = idxs(idxp)
    if (widx < i) then 
      write(psb_err_unit,*) 'Warning: lower triangle in upper copy',widx,i,idxp,idxs(idxp)
      cycle
    end if
    if (widx > m) then 
      cycle
    end if
    witem     = row(widx)
    !
    ! Dropping rule based on the 2-norm. But keep the jmaxup-th entry anyway.
    !
    if ((widx /= jmaxup) .and. (widx /= i) .and. (abs(witem) < thres*nrmi)) then 
      cycle 
    end if
    if ((widx/=jmaxup).and.(nz > nup+fill_in)) then
      if (abs(witem) < wmin) cycle
    endif
    wmin = min(abs(witem),wmin)
    nz       = nz + 1
    xw(nz)   = witem 
    xwid(nz) = widx
    call psb_insert_heap(witem,widx,heap,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_insert_heap')
      goto 9999
    end if

  end do

  !
  ! Now we have to take out the first nup-fill_in entries. But make sure
  ! we include entry jmaxup.
  !
  if (nz <= nup+fill_in) then
    ! 
    ! Just copy everything from xw
    !
    fndmaxup=.true.
  else
    fndmaxup = .false.
    nz = nup+fill_in
    do k=1,nz
      call psb_heap_get_first(witem,widx,heap,info)
      xw(k)   = witem
      xwid(k) = widx
      if (widx == jmaxup) fndmaxup=.true.
    end do
  end if
  if ((i<jmaxup).and.(jmaxup<=m)) then 
    if (.not.fndmaxup) then 
      ! 
      ! Include entry jmaxup, if it is not already there.
      ! Put it in the place of the smallest coefficient. 
      !
      xw(nz)   = row(jmaxup) 
      xwid(nz) = jmaxup
    endif
  end if

  !
  ! Now we put things back into ascending column order
  !
  call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)

  !
  ! Copy out the upper part of the row
  !
  do k=1,nz
    l2     = l2 + 1 
    if (size(val) < l2) then
      ! 
      ! Figure out a good reallocation size!
      ! 
      isz  = max(int(1.2*l2),l2+100)
      call psb_realloc(isz,val,info) 
      if (info == psb_success_) call psb_realloc(isz,ja,info) 
      if (info /= psb_success_) then 
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999
      end if
    end if
    ja(l2)   = xwid(k)
    val(l2)  = xw(indx(k))
  end do

  !
  ! Set row to zero
  !
  do idxp=1,nidx
    row(idxs(idxp)) = dzero
  end do

  irp(i+1) = l2 + 1

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine mld_dinvt_copyout



subroutine mld_dsparse_ainvt(n,a,z,nzrmax,sp_thresh,info)
  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dsparse_ainvt

  implicit none 
  integer, intent(in)                  :: n
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_dspmat_type), intent(inout) :: z
  integer, intent(in)                  :: nzrmax
  real(psb_dpk_), intent(in)           :: sp_thresh
  integer, intent(out)                 :: info

  integer :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2
  integer, allocatable        :: ia(:), ja(:), iz(:),jz(:) 
  real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
  integer, allocatable        :: uplevs(:), rowlevs(:),idxs(:)
  real(psb_dpk_), allocatable :: row(:)
  type(psb_d_coo_sparse_mat)  :: trw
  type(psb_d_csr_sparse_mat)  :: acsr, zcsr
  integer                  :: ktrw, nidx, nlw,nup,jmaxup
  type(psb_int_heap)       :: heap
  real(psb_dpk_)     :: alpha, nrmi
  character(len=20)  :: name='mld_sp_ainvt'


  if(psb_get_errstatus() /= psb_success_) return 
  info = psb_success_
  call psb_erractionsave(err_act)

  if (.not.(a%is_triangle().and.a%is_unit().and.a%is_upper())) then 
    write(psb_err_unit,*) 'Wrong A ' 
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='wrong A')
    goto 9999      
  end if
  call a%cp_to(acsr)
  call trw%allocate(0,0,1)
  if (info == psb_success_) allocate(zw(n),iz(n),valz(n),&
       & row(n),rowlevs(n),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  call zcsr%allocate(n,n,n*nzrmax)
  call zcsr%set_triangle()
  call zcsr%set_unit(.false.)
  call zcsr%set_upper()
  ! 
  !
  nzz        = 0
  row(:)     = dzero 
  rowlevs(:) = 0
  l2         = 0
  zcsr%irp(1) = 1 
  
  outer: do i = 1, n-1
    ! ZW = e_i
    call mld_invt_copyin(i,n,acsr,i,1,n,nlw,nup,jmaxup,nrmi,row,&
         & heap,rowlevs,ktrw,trw,info,sign=-done)
    if (info /= 0) exit
    row(i) = done
    ! Adjust norm
    if (nrmi < done) then 
      nrmi = sqrt(done + nrmi**2)
    else 
      nrmi = nrmi*sqrt(done+done/(nrmi**2))
    end if

    call mld_invt(sp_thresh,i,nrmi,row,heap,rowlevs,&
         & acsr%ja,acsr%irp,acsr%val,nidx,idxs,info)
    if (info /= 0) exit
!!$    write(0,*) 'Calling copyout ',nzrmax,nlw,nup,nidx,l2
    call mld_invt_copyout(nzrmax,sp_thresh,i,n,nlw,nup,jmaxup,nrmi,row,&
         & nidx,idxs,l2,zcsr%ja,zcsr%irp,zcsr%val,info)
    if (info /= 0) exit
    nzz = l2
  end do outer
  if (info /= psb_success_) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='mainloop')
    goto 9999
  end if

  ipz1 = nzz+1
  call psb_ensure_size(ipz1,zcsr%val,info)
  call psb_ensure_size(ipz1,zcsr%ja,info)
  zcsr%val(ipz1) = done
  zcsr%ja(ipz1)  = n
  zcsr%irp(n+1)  = ipz1+1 
  
  call z%mv_from(zcsr)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_dsparse_ainvt


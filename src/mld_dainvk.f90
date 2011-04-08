
subroutine mld_dinvk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,ktrw,trw,info,sign,inlevs)

  use psb_base_mod
  use mld_dainv_mod, mld_protect_name => mld_dinvk_copyin

  implicit none

  ! Arguments 
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_dspmat_type), intent(inout) :: trw
  integer, intent(in)                  :: i,m,jmin,jmax
  integer, intent(inout)               :: ktrw,info
  integer, intent(inout)               :: rowlevs(:)
  real(psb_dpk_), intent(inout)        :: row(:)
  type(psb_int_heap), intent(inout)    :: heap
  real(psb_dpk_), optional, intent(in) :: sign
  integer, intent(in), optional        :: inlevs(:)

  ! Local variables
  integer             :: k,j,irb,err_act
  integer, parameter  :: nrb=16
  real(psb_dpk_)      :: sign_
  character(len=20), parameter  :: name='invk_copyin'
  character(len=20)             :: ch_err

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  call psb_init_heap(heap,info) 

  if (present(sign)) then 
    sign_ = sign
  else
    sign_ = done
  end if

  if (psb_toupper(a%fida) == 'CSR') then

    !
    ! Take a fast shortcut if the matrix is stored in CSR format
    !
    if (present(inlevs)) then 
      do j = a%ia2(i), a%ia2(i+1) - 1
        k          = a%ia1(j)
        if ((jmin<=k).and.(k<=jmax)) then 
          row(k)     = sign_ * a%aspk(j)
          rowlevs(k) = inlevs(j)
          call psb_insert_heap(k,heap,info)
        end if
      end do
    else
      do j = a%ia2(i), a%ia2(i+1) - 1
        k          = a%ia1(j)
        if ((jmin<=k).and.(k<=jmax)) then 
          row(k)     = sign_ * a%aspk(j)
          rowlevs(k) = 0
          call psb_insert_heap(k,heap,info)
        end if
      end do
    end if
  else

    !
    ! Otherwise use psb_sp_getblk, slower but able (in principle) of 
    ! handling any format. In this case, a block of rows is extracted
    ! instead of a single row, for performance reasons, and these
    ! rows are copied one by one into the array row, through successive
    ! calls to invk_copyin.
    !

    if ((mod(i,nrb) == 1).or.(nrb == 1)) then 
      irb = min(m-i+1,nrb)
      call psb_sp_getblk(i,a,trw,info,lrw=i+irb-1)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='psb_sp_getblk'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      ktrw=1
    end if

    do 
      if (ktrw > trw%infoa(psb_nnz_)) exit
      if (trw%ia1(ktrw) > i) exit
      k          = trw%ia2(ktrw)
      if ((jmin<=k).and.(k<=jmax)) then 
        row(k)     = sign_*trw%aspk(ktrw)
        rowlevs(k) = 0
        call psb_insert_heap(k,heap,info)
      end if
      ktrw       = ktrw + 1
    enddo
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

end subroutine mld_dinvk_copyin

subroutine mld_dinvk(fill_in,i,row,rowlevs,heap,uia1,uia2,uaspk,uplevs,nidx,idxs,info)

  use psb_base_mod
  use mld_dainv_mod, mld_protect_name => mld_dinvk

  implicit none 

  ! Arguments
  type(psb_int_heap), intent(inout)    :: heap 
  integer, intent(in)                  :: i, fill_in
  integer, intent(inout)               :: nidx,info
  integer, intent(inout)               :: rowlevs(:)
  integer, allocatable, intent(inout)  :: idxs(:)
  integer, intent(in)                  :: uia1(:),uia2(:),uplevs(:)
  real(psb_dpk_), intent(in)           :: uaspk(:)
  real(psb_dpk_), intent(inout)        :: row(:)

  ! Local variables
  integer             :: k,j,lrwk,jj,lastk, iret
  real(psb_dpk_)      :: rwk


  info = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i

  !
  ! Do while there are indices to be processed
  !
  do
    ! Beware: (iret < 0) means that the heap is empty, not an error.
    call psb_heap_get_first(k,heap,iret) 
    if (iret < 0) then 
!!$        write(psb_err_unit,*) 'IINVK: ',i,' returning at ',lastk
      return
    end if

    ! 
    ! Just in case an index has been put on the heap more than once.
    !
    if (k == lastk) cycle

    lastk = k 
    nidx = nidx + 1
    if (nidx>size(idxs)) then 
      call psb_realloc(nidx+psb_heap_resize,idxs,info)
      if (info /= psb_success_) return
    end if
    idxs(nidx) = k

    if ((row(k) /= dzero).and.(rowlevs(k) <= fill_in)) then 
      !
      ! Note: since U is scaled while copying it out (see iluk_copyout),
      ! we can use rwk in the update below
      ! 
      rwk    = row(k)
      lrwk   = rowlevs(k)

      do jj=uia2(k),uia2(k+1)-1
        j = uia1(jj)
        if (j<=k) then 
          info = -i
          return
        endif
        !
        ! Insert the index into the heap for further processing.
        ! The fill levels are initialized to a negative value. If we find
        ! one, it means that it is an as yet untouched index, so we need
        ! to insert it; otherwise it is already on the heap, there is no
        ! need to insert it more than once. 
        !
        if (rowlevs(j)<0) then 
          call psb_insert_heap(j,heap,info)
          if (info /= psb_success_) return
          rowlevs(j) = abs(rowlevs(j))
        end if
        !
        ! Update row(j) and the corresponding fill level
        !
        row(j)     = row(j) - rwk * uaspk(jj)
        rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
      end do

    end if
  end do

end subroutine mld_dinvk

subroutine mld_dinvk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
     &  l2,uia1,uia2,uaspk,info)

  use psb_base_mod
  use mld_dainv_mod, mld_protect_name => mld_dinvk_copyout

  implicit none 

  ! Arguments
  integer, intent(in)                        :: fill_in, i, m, nidx
  integer, intent(inout)                     :: l2, info
  integer, intent(inout)                     :: rowlevs(:), idxs(:)
  integer, allocatable, intent(inout)        :: uia1(:), uia2(:)
  real(psb_dpk_), allocatable, intent(inout) :: uaspk(:)
  real(psb_dpk_), intent(inout)              :: row(:)

  ! Local variables
  integer               :: j,isz,err_act,int_err(5),idxp
  character(len=20), parameter  :: name='mld_diluk_factint'
  character(len=20)             :: ch_err

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)


  do idxp=1,nidx

    j = idxs(idxp)


    if (j>=i) then 
      !
      ! Copy the upper part of the row
      ! 
      if (rowlevs(j) <= fill_in) then 
        l2     = l2 + 1 
        if (size(uaspk) < l2) then 
          ! 
          ! Figure out a good reallocation size!
          !
          isz  = max((l2/i)*m,int(1.2*l2),l2+100)
          call psb_realloc(isz,uaspk,info) 
          if (info == psb_success_) call psb_realloc(isz,uia1,info) 
          if (info /= psb_success_) then 
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='Allocate')
            goto 9999
          end if
        end if
        uia1(l2)   = j
        uaspk(l2)  = row(j)
      end if
      !
      ! Re-initialize row(j) and rowlevs(j)
      !
      row(j)     = dzero
      rowlevs(j) = -(m+1)
    end if
  end do

  uia2(i+1) = l2 + 1

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine mld_dinvk_copyout

subroutine mld_dsparse_ainvk(n,a,z,fill_in,sp_thresh,info,inlevs)
  use psb_base_mod
  use mld_dainv_mod, mld_protect_name => mld_dsparse_ainvk

  integer, intent(in)                  :: n
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_dspmat_type), intent(inout) :: z
  integer, intent(in)                  :: fill_in
  real(psb_dpk_), intent(in)           :: sp_thresh
  integer, intent(out)                 :: info
  integer, intent(in), optional        :: inlevs(:)

  integer :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2 
  integer, allocatable        :: ia(:), ja(:), iz(:),jz(:) 
  real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
  integer, allocatable        :: uplevs(:), rowlevs(:),idxs(:)
  real(psb_dpk_), allocatable :: row(:)
  type(psb_dspmat_type)    :: trw
  integer                  :: ktrw, nidx
  type(psb_int_heap)       :: heap

  real(psb_dpk_)     :: alpha
  character(len=20)  :: name='mld_sp_ainvk'


  if(psb_get_errstatus() /= psb_success_) return 
  info = psb_success_
  call psb_erractionsave(err_act)


  if (psb_tolower(a%fida) /= 'csr') then 
    write(psb_err_unit,*) 'AFMT wrong ',a%fida
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='Afmt wrong')
    goto 9999      
  end if
  if (psb_tolower(a%descra(1:3)) /= 'tuu') then 
    write(psb_err_unit,*) 'DESCR wrong ',a%descra
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='descra wrong')
    goto 9999      
  end if

  allocate(uplevs(size(a%aspk)),rowlevs(n),row(n),stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='Allocate')
    goto 9999
  end if
  uplevs(:)  = 0
  row(:)     = dzero
  rowlevs(:) = -(n+1)

  allocate(zw(n),iz(n),valz(n),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  call psb_ensure_size(n+1, z%ia2,  info)
  call psb_ensure_size(n+1, z%ia1,  info)
  call psb_ensure_size(n+1, z%aspk, info)
  call psb_ensure_size(n+1, idxs,  info)


  ! 
  !
  z%descra  = 'GUN'
  z%fida    = 'CSR'
  z%m       = n
  z%k       = n
  z%ia2(1)  = 1
  nzz       = 0

  l2 = 0
  outer: do i = 1, n-1
    ! ZW = e_i
    call mld_invk_copyin(i,n,a,1,n,row,rowlevs,heap,ktrw,trw,info,&
         & sign=-done,inlevs=inlevs)
    row(i)     = done
    rowlevs(i) = 0
!!$      call psb_insert_heap(i,heap,info) ! No we don't want to put I in. 

    ! Update loop
    call mld_invk(fill_in,i,row,rowlevs,heap,a%ia1,a%ia2,a%aspk,uplevs,nidx,idxs,info)

    call mld_invk_copyout(fill_in,i,n,row,rowlevs,nidx,idxs,&
         & l2,z%ia1,z%ia2,z%aspk,info)

    nzz = l2
  end do outer
  if (info /= psb_success_) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='mainloop')
    goto 9999
  end if
  ipz1 = nzz+1
  z%aspk(ipz1) = done
  z%ia1(ipz1)  = n
  z%ia2(n+1)   = ipz1+1 
  call psb_spcnv(z,info,afmt='CSR')
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine mld_dsparse_ainvk



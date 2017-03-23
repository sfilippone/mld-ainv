!  
!   
!                       MLD-AINV: Approximate Inverse plugin for
!                             MLD2P4  version 2.0
!    
!    (C) Copyright 2012
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine mld_d_invt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
     & nidx,idxs,l2,ja,irp,val,info)

  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_d_invt_copyout

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
  type(psb_d_idx_heap)         :: heap
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
!!$  write(0,*) 'invt_copyout ',nidx,nup+fill_in
  call heap%init(info,dir=psb_asort_down_)

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
    call heap%insert(witem,widx,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_insert_heap')
      goto 9999
    end if
  end do

  if (nz > 1) then 
    write(psb_err_unit,*) 'Warning: lower triangle from invt???? '
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
  if (.false.) then 
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
      call heap%insert(witem,widx,info)
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
        call heap%get_first(witem,widx,info)
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

  else if (.true.) then 

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
      if ((widx /= i) .and. (abs(witem) < thres*nrmi)) then 
        cycle 
      end if
      if (nz > nup+fill_in) then
        if (abs(witem) < wmin) cycle
      endif
      wmin = min(abs(witem),wmin)
      nz       = nz + 1
      xw(nz)   = witem 
      xwid(nz) = widx
      call heap%insert(witem,widx,info)
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
    if (nz >  nup+fill_in) then
      nz = nup+fill_in
      do k=1,nz
        call heap%get_first(witem,widx,info)
        xw(k)   = witem
        xwid(k) = widx
      end do
    end if
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

end subroutine mld_d_invt_copyout


subroutine mld_d_sparsify(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,info,istart,iheap,ikr)
  use psb_base_mod
  implicit none 

  real(psb_dpk_), intent(in)  :: sp_thresh
  integer, intent(in)         :: idiag, n, nzrmax
  real(psb_dpk_), intent(inout)  :: zw(:)
  integer, intent(out)        :: nz
  integer, intent(out)        :: iz(:)
  real(psb_dpk_), intent(out) :: valz(:)
  integer, intent(out)        :: info
  integer, intent(in), optional :: istart
  type(psb_int_heap), optional :: iheap
  integer, optional            :: ikr(:)

  integer :: i, istart_, last_i, iret,k
  real(psb_dpk_)     :: witem
  integer            :: widx
  real(psb_dpk_), allocatable :: xw(:)
  integer, allocatable        :: xwid(:), indx(:)
  type(psb_dreal_idx_heap)    :: heap


  info = psb_success_
  istart_ = 1
  if (present(istart)) istart_ = max(1,istart)
  if (.false.) then 
    nz = 0
    do i=istart_, n
      if ((i == idiag).or.(abs(zw(i)) >= sp_thresh)) then 
        nz       = nz + 1 
        iz(nz)   = i
        valz(nz) = zw(i) 
      end if
    end do

  else

    allocate(xw(nzrmax),xwid(nzrmax),indx(nzrmax),stat=info)
    if (info /= psb_success_) then 
      return
    end if

    call psb_init_heap(heap,info,dir=psb_asort_down_)

    ! Keep at least the diagonal
    nz = 0 

    if (present(iheap)) then 
      if (.not.(present(ikr))) then 
        write(psb_err_unit,*) 'Error: if IHEAP then also IKR'
        info = -1
        return
      end if
      last_i = -1
      do 
        call psb_heap_get_first(i,iheap,iret) 
        if (iret < 0) exit
        ! An index may have been put on the heap more than once.
        if (i == last_i) cycle
        last_i = i 
        if (i == idiag) then 
          xw(1)   = zw(i)
          xwid(1) = i
        else if (abs(zw(i)) >= sp_thresh) then 
          call psb_insert_heap(zw(i),i,heap,info)
        end if
        zw(i)  = dzero
        ikr(i) = 0
      end do

    else

      do i=istart_, n
        if (i == idiag) then 
          xw(1)   = zw(i)
          xwid(1) = i
        else if (abs(zw(i)) >= sp_thresh) then 
          call psb_insert_heap(zw(i),i,heap,info)
        end if
        zw(i) = dzero
      end do
    end if

    k = 1
    do 
      if (k == nzrmax) exit 
      call psb_heap_get_first(witem,widx,heap,info)
      if (info == -1) then 
        info = psb_success_
        exit 
      endif
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        return
      end if
      k = k + 1 
      xw(k)   = witem
      xwid(k) = widx
    end do
    call psb_free_heap(heap,info)
    nz = k 
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)
    do i=1, nz
      valz(i) = xw(indx(i))
      iz(i)   = xwid(i)
    end do
!!$    write(0,*) 'sparsify output for idiag ',idiag,' :',nz

  end if

  return

end subroutine mld_d_sparsify


subroutine mld_d_sparsify_list(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,lhead,listv,ikr,info)
  use psb_base_mod
  implicit none 

  real(psb_dpk_), intent(in)  :: sp_thresh
  integer, intent(in)         :: idiag, n, nzrmax
  real(psb_dpk_), intent(inout)  :: zw(:)
  integer, intent(out)        :: nz
  integer, intent(out)        :: iz(:)
  real(psb_dpk_), intent(out) :: valz(:)
  integer, intent(out)        :: info
  integer, intent(inout)      :: lhead, listv(:)
  integer                     :: ikr(:)

  integer :: i, istart_, last_i, iret,k,current, next
  real(psb_dpk_)     :: witem
  integer            :: widx
  real(psb_dpk_), allocatable :: xw(:)
  integer, allocatable        :: xwid(:), indx(:)


  info = psb_success_
  istart_ = 1
  allocate(xw(n),xwid(n),indx(n),stat=info)

  current = lhead 
  lhead = -1 
  i = 0
  do while (current >0) 
    i = i + 1 
    xw(i)   = zw(current) 
    xwid(i) = current 
    
    if (current == idiag) then 
      ! Bring the diagona into first position
      witem = xw(1) 
      widx  = xwid(1)
      xw(1)   = xw(i) 
      xwid(1) = xwid(i)
      xw(i)   = witem
      xwid(i) = widx
    end if            

    zw(current)  = dzero
    ikr(current) = 0 
    
    next = listv(current)
    listv(current) = -1
    current = next
  end do
  nz = i
  if (nz > 2) call psb_msort(xw(2:nz),ix=xwid(2:nz),dir=psb_asort_down_,flag=psb_sort_keep_idx_)
!!$  write(0,*) 'Done first msort '
  i = 2
  do while (i<=nz)
    if (abs(xw(i)) < sp_thresh) exit
    i = i + 1 
  end do
!!$  write(0,*) 'NZ ',nz, i, nzrmax
  nz = max(1,min(i-1,nzrmax))
  call psb_msort(xwid(1:nz),ix=indx(1:nz),dir=psb_sort_up_)
!!$  write(0,*) 'Done second msort '
  
  do i=1, nz
    valz(i) = xw(indx(i))
    iz(i)   = xwid(i)
  end do
!!$  write(0,*) 'sparsify output for idiag ',idiag,' :',nz

  return

end subroutine mld_d_sparsify_list


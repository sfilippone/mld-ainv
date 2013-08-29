subroutine mld_d_sp_drop(idiag,nzrmax,sp_thresh,nz,iz,valz,info)
  use psb_base_mod
  implicit none 
  real(psb_dpk_), intent(in)    :: sp_thresh
  integer, intent(in)           :: idiag, nzrmax
  integer, intent(inout)        :: nz
  integer, intent(inout)        :: iz(:)
  real(psb_dpk_), intent(inout) :: valz(:)
  integer, intent(out)          :: info

  integer :: i, j, idf, nw
  real(psb_dpk_)     :: witem
  integer            :: widx
  real(psb_dpk_), allocatable :: xw(:)
  integer, allocatable        :: xwid(:), indx(:)


  info = psb_success_

  if (nz > min(size(iz),size(valz))) then 
    write(0,*) 'Serious size problem ',nz,size(iz),size(valz)
    info = -2
    return
  end if
  allocate(xw(nz),xwid(nz),indx(nz),stat=info) 
  if (info /= psb_success_) then 
    write(psb_err_unit,*) ' Memory allocation failure in sp_drop',nz,info
    return
  endif

  ! Always keep the diagonal element
  idf = -1 
  do i=1, nz
    if (iz(i) == idiag) then 
      idf     = i
      witem   = valz(i)
      widx    = iz(i)
      valz(i) = valz(1) 
      iz(i)   = iz(1) 
      valz(1) = witem
      iz(1)   = widx
      exit
    end if
  end do

  if (idf == -1) then

    xw(1:nz) = valz(1:nz)
    call psb_qsort(xw(1:nz),indx(1:nz),dir=psb_asort_down_)
    do i=1, nz
      xwid(i) = iz(indx(i))
    end do
    nw = min(nz,nzrmax)
    do 
      if (nw <= 1) exit
      if (abs(xw(nw)) < sp_thresh) then 
        nw = nw - 1
      else 
        exit
      end if
    end do
    nw = max(nw, 1)

  else

    nw = nz-1

    xw(1:nw) = valz(2:nz)

    call psb_qsort(xw(1:nw),indx(1:nw),dir=psb_asort_down_)
    nw = min(nw,nzrmax-1)
    do 
      if (nw <= 1) exit
      if (abs(xw(nw)) < sp_thresh) then 
        nw = nw - 1
      else 
        exit
      end if
    end do

    do i=1, nw
      xwid(i) = iz(1+indx(i))
    end do
    nw       = nw + 1 
    xw(nw)   = valz(1)
    xwid(nw) = iz(1)
  end if

  call psb_msort(xwid(1:nw),indx(1:nw),dir=psb_sort_up_)

  do i=1, nw
    valz(i) = xw(indx(i))
    iz(i)   = xwid(i)
  end do
  nz = nw
  if (nz>nzrmax) write(0,*) 'in sp_drop: ',nw,nzrmax,nz
  deallocate(xw,xwid,indx,stat=info) 
  if (info /= psb_success_) then 
    write(psb_err_unit,*) ' Memory deallocation failure in sp_drop',info
    return
  endif
  return
end subroutine mld_d_sp_drop


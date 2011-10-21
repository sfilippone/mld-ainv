!
module mld_base_ainv_mod
  
  use mld_base_prec_type

  interface sp_drop
    module procedure d_sp_drop
  end interface

  interface rwclip
    module procedure drwclip
  end interface
  
  interface sparsify
    module procedure d_sparsify
  end interface
  
  integer, parameter   :: mld_inv_fillin_ = mld_ifpsz_ + 1
  integer, parameter   :: mld_ainv_alg_   = mld_inv_fillin_+1
  integer, parameter   :: mld_ainv_orth1_ = mld_max_sub_solve_+1
  integer, parameter   :: mld_ainv_orth2_ = mld_ainv_orth1_+1
  integer, parameter   :: mld_ainv_orth3_ = mld_ainv_orth2_+1
  integer, parameter   :: mld_ainv_orth4_ = mld_ainv_orth3_+1
  integer, parameter   :: mld_inv_thresh_ = mld_ainv_orth4_ + 1

contains

  subroutine drwclip(nz,ia,ja,val,imin,imax,jmin,jmax)
    use psb_base_mod
    implicit none 
    integer, intent(inout) :: nz
    integer, intent(inout) :: ia(*), ja(*)
    real(psb_dpk_), intent(inout) :: val(*)
    integer, intent(in)    :: imin,imax,jmin,jmax

    integer :: i,j 

    j = 0
    do i=1, nz
      if ((imin <= ia(i)).and.&
           & (ia(i) <= imax).and.&
           & (jmin <= ja(i)).and.&
           & (ja(i) <= jmax) ) then 
        j = j + 1 
        ia(j) = ia(i) 
        ja(j) = ja(i)
        val(j) = val(i)
      end if
    end do
    nz = j 
  end subroutine drwclip


  subroutine d_sparsify(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,info,istart,iheap,ikr)
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
    type(psb_double_idx_heap)   :: heap


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

    end if

    return

  end subroutine d_sparsify

  subroutine d_sp_drop(idiag,nzrmax,sp_thresh,nz,iz,valz,info)
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
    if (nz <= nzrmax) return 
!!$    write(0,*) 'sp_drop start',nz
!!$    write(0,*) 'sp_drop start',nz,size(iz),size(valz)

    if (nz > min(size(iz),size(valz))) then 
      write(0,*) 'Serious size problem ',nz,size(iz),size(valz)
      info = -2
      return
    end if
!!$    write(0,*) 'sp_drop allocation',nz
    allocate(xw(nz),xwid(nz),indx(nz),stat=info) 
!!$    write(0,*) 'sp_drop allocation',nz
    if (info /= psb_success_) then 
      write(psb_err_unit,*) ' Memory allocation in sp_drop'
      return
    endif

    ! Always keep the diagonal element
!!$    write(0,*) 'sp_drop looking for diag ',idiag
    call flush(0)
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
!!$    write(0,*) 'sp_drop diag found :',idf

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
      nw = min(nw,nzrmax)
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
!!$    write(0,*) 'sp_drop into msort ',nw

    call psb_msort(xwid(1:nw),indx(1:nw),dir=psb_sort_up_)
    do i=1, nw
      valz(i) = xw(indx(i))
      iz(i)   = xwid(i)
    end do
    nz = nw
    deallocate(xw,xwid,indx,stat=info) 
    if (info /= psb_success_) then 
      write(psb_err_unit,*) ' Memory allocation in sp_drop'
      return
    endif

    return
  end subroutine d_sp_drop


end module mld_base_ainv_mod


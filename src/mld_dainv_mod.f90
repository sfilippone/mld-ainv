!!$ 
!!$ 
!!$                           MLD2P4  version 1.2
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.3.1)
!!$  
!!$  (C) Copyright 2008,2009
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MLD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$
! File: mld_dainv_bld.f90
!
! Subroutine: mld_dainv_bld
! Version:    real
!
!  For details on the above factorizations see
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!
!
! Arguments:
!
module mld_dainv_mod

  interface mld_sparse_orthbase
    subroutine mld_dsparse_orthbase(alg,n,a,p,z,nzrmax,sp_thresh,info)
      use psb_base_mod
      integer, intent(in)                  :: alg,n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: nzrmax
      real(psb_dpk_), intent(in)           :: sp_thresh
      real(psb_dpk_), intent(out)          :: p(:)
      integer, intent(out)                 :: info

    end subroutine mld_dsparse_orthbase
  end interface

  interface sp_drop
    module procedure d_sp_drop
!!$    subroutine d_sp_drop(idiag,nzrmax,sp_thresh,nz,iz,valz,info)
!!$      use psb_base_mod
!!$      real(psb_dpk_), intent(in)  :: sp_thresh
!!$      integer, intent(in)         :: idiag, nzrmax
!!$      integer, intent(inout)        :: nz
!!$      integer, intent(inout)        :: iz(:)
!!$      real(psb_dpk_), intent(inout) :: valz(:)
!!$      integer, intent(out)        :: info
!!$
!!$    end subroutine d_sp_drop
  end interface

  interface mld_invt_copyin
    subroutine mld_dinvt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
         & irwt,ktrw,trw,info,sign)
      use psb_base_mod
      implicit none 
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_d_coo_sparse_mat), intent(inout) :: trw
      integer, intent(in)                  :: i, m,jmin,jmax,jd
      integer, intent(inout)               :: ktrw,nlw,nup,jmaxup,info
      integer, intent(inout)               :: irwt(:)
      real(psb_dpk_), intent(inout)        :: nrmi,row(:)
      type(psb_int_heap), intent(inout)    :: heap
      real(psb_dpk_), intent(in), optional :: sign

    end subroutine mld_dinvt_copyin
  end interface

  interface mld_invt
    subroutine mld_dinvt(thres,i,nrmi,row,heap,irwt,ja,irp,val,nidx,idxs,info)

      use psb_base_mod
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

    end subroutine mld_dinvt
  end interface

  interface mld_invt_copyout
    subroutine mld_dinvt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
         & nidx,idxs,l2,ja,irp,val,info)

      use psb_base_mod

      implicit none 

      ! Arguments
      integer, intent(in)                       :: fill_in,i,m,nidx,nlw,nup,jmaxup
      integer, intent(in)                       :: idxs(:)
      integer, intent(inout)                    :: l2, info
      integer, allocatable, intent(inout)       :: ja(:),irp(:)
      real(psb_dpk_), intent(in)                :: thres,nrmi
      real(psb_dpk_),allocatable, intent(inout) :: val(:)
      real(psb_dpk_), intent(inout)             :: row(:)

    end subroutine mld_dinvt_copyout
  end interface

  interface  mld_sparse_ainvt
    subroutine mld_dsparse_ainvt(n,a,z,nzrmax,sp_thresh,info)
      use psb_base_mod
      implicit none 
      integer, intent(in)                  :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: nzrmax
      real(psb_dpk_), intent(in)           :: sp_thresh
      integer, intent(out)                 :: info

    end subroutine mld_dsparse_ainvt
  end interface

  interface  mld_invk_copyin
    subroutine mld_dinvk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,&
         & ktrw,trw,info,sign,inlevs)

      use psb_base_mod
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

    end subroutine mld_dinvk_copyin
  end interface

  interface mld_invk
    subroutine mld_dinvk(fill_in,i,row,rowlevs,heap,ja,irp,val,uplevs,nidx,idxs,info)

      use psb_base_mod
      implicit none 

      ! Arguments
      type(psb_int_heap), intent(inout)    :: heap 
      integer, intent(in)                  :: i, fill_in
      integer, intent(inout)               :: nidx,info
      integer, intent(inout)               :: rowlevs(:)
      integer, allocatable, intent(inout)  :: idxs(:)
      integer, intent(in)                  :: ja(:),irp(:), uplevs(:)
      real(psb_dpk_), intent(in)           :: val(:)
      real(psb_dpk_), intent(inout)        :: row(:)


    end subroutine mld_dinvk
  end interface

  interface mld_invk_copyout
    subroutine mld_dinvk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
         &  l2,ja,irp,val,info)

      use psb_base_mod

      implicit none 

      ! Arguments
      integer, intent(in)                        :: fill_in, i, m, nidx
      integer, intent(inout)                     :: l2, info
      integer, intent(inout)                     :: rowlevs(:), idxs(:)
      integer, allocatable, intent(inout)        :: ja(:), irp(:)
      real(psb_dpk_), allocatable, intent(inout) :: val(:)
      real(psb_dpk_), intent(inout)              :: row(:)

    end subroutine mld_dinvk_copyout
  end interface
  
  interface mld_sparse_ainvk
    subroutine mld_dsparse_ainvk(n,a,z,fill_in,sp_thresh,info,inlevs)
      use psb_base_mod
      integer, intent(in)                  :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: fill_in
      real(psb_dpk_), intent(in)           :: sp_thresh
      integer, intent(out)                 :: info
      integer, intent(in), optional        :: inlevs(:)

    end subroutine mld_dsparse_ainvk
  end interface

  interface rwclip
    module procedure drwclip
  end interface
  
  interface sparsify
    module procedure d_sparsify
  end interface

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
    real(psb_dpk_), intent(in)  :: sp_thresh
    integer, intent(in)         :: idiag, nzrmax
    integer, intent(inout)        :: nz
    integer, intent(inout)        :: iz(:)
    real(psb_dpk_), intent(inout) :: valz(:)
    integer, intent(out)        :: info

    integer :: i, j, idf, nw
    real(psb_dpk_)     :: witem
    integer            :: widx
    real(psb_dpk_), allocatable :: xw(:)
    integer, allocatable        :: xwid(:), indx(:)


    info = psb_success_

    allocate(xw(nz),xwid(nz),indx(nz),stat=info) 
    if (info /= psb_success_) then 
      write(psb_err_unit,*) ' Memory allocation in sp_drop'
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
      nw = min(nw,nzrmax)
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

    call psb_msort(xwid(1:nw),indx(1:nw),dir=psb_sort_up_)
    do i=1, nw
      valz(i) = xw(indx(i))
      iz(i)   = xwid(i)
    end do
    nz = nw

    return
  end subroutine d_sp_drop


end module mld_dainv_mod


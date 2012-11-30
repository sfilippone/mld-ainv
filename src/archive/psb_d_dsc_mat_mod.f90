!!$
!!$ 
!!$                     MLD-AINV: Approximate Inverse plugin for
!!$                           MLD2P4  version 2.0
!!$  
!!$  (C) Copyright 2012
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
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
!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!
! package: psb_d_dsc_mat_mod
!
! This module contains the definition of the psb_d_dsc_sparse_mat type
! which implements an actual storage format (the DSC in this case) for
! a sparse matrix as well as the related methods (those who are
! specific to the type and could not be defined higher in the
! hierarchy). We are at the bottom level of the inheritance chain.
! 
module psb_d_spvect_mod

  use psb_base_mod, only : psb_safe_cpy, psb_dpk_

  type psb_d_spvect
    integer :: nz=0
    integer, allocatable        :: idx(:)
    real(psb_dpk_), allocatable :: val(:)
  end type psb_d_spvect

  interface psb_realloc
    module procedure psb_reallocate_d_spvect
  end interface psb_realloc

  interface psb_safe_ab_cpy
    module procedure psb_d_spvect_ab_copy
  end interface psb_safe_ab_cpy

contains

  Subroutine psb_reallocate_d_spvect(len,rrax,info,pad,lb)
    use psb_base_mod

    ! ...Subroutine Arguments  
    Integer,intent(in) :: len
    type(psb_d_spvect),allocatable, intent(inout) :: rrax(:)
    integer :: info
    type(psb_d_spvect), optional, intent(in) :: pad
    integer, optional, intent(in) :: lb

    ! ...Local Variables
    type(psb_d_spvect),allocatable  :: tmp(:)
    Integer :: dim,err_act,err, lb_, lbi,ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1d'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (debug) write(psb_err_unit,*) 'reallocate D',len

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name,i_err=(/len,0,0,0,0/),a_err='type(psb_d_spvect)')
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name,i_err=(/len,0,0,0,0/),a_err='type(psb_d_spvect)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call move_alloc(tmp,rrax)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name,i_err=(/len,0,0,0,0/),a_err='type(psb_d_spvect)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_reallocate_d_spvect

  subroutine psb_d_spvect_ab_copy(vin,vout,info)
    use psb_base_mod
    type(psb_d_spvect), intent(in), allocatable :: vin(:)
    type(psb_d_spvect), intent(out), allocatable :: vout(:)
    integer :: info
    !
    Integer :: isz,err_act,lb,i
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    
    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error()
    end if
    return

  end subroutine psb_d_spvect_ab_copy

end module psb_d_spvect_mod

module psb_d_dsc_mat_mod

  use psb_d_base_mat_mod
  use psb_d_spvect_mod

  type, extends(psb_d_base_sparse_mat) :: psb_d_dsc_sparse_mat
    
    type(psb_d_spvect), allocatable :: cols(:) 
    
  contains
    procedure, pass(a) :: get_size     => d_dsc_get_size
    procedure, pass(a) :: get_nzeros   => d_dsc_get_nzeros
    procedure, nopass  :: get_fmt      => d_dsc_get_fmt
    procedure, pass(a) :: sizeof       => d_dsc_sizeof
    procedure, pass(a) :: d_csmm       => psb_d_dsc_csmm
    procedure, pass(a) :: d_csmv       => psb_d_dsc_csmv
    procedure, pass(a) :: d_inner_cssm => psb_d_dsc_cssm
    procedure, pass(a) :: d_inner_cssv => psb_d_dsc_cssv
    procedure, pass(a) :: d_scals      => psb_d_dsc_scals
    procedure, pass(a) :: d_scal       => psb_d_dsc_scal
    procedure, pass(a) :: csnmi        => psb_d_dsc_csnmi
    procedure, pass(a) :: csnm1        => psb_d_dsc_csnm1
    procedure, pass(a) :: rowsum       => psb_d_dsc_rowsum
    procedure, pass(a) :: arwsum       => psb_d_dsc_arwsum
    procedure, pass(a) :: colsum       => psb_d_dsc_colsum
    procedure, pass(a) :: aclsum       => psb_d_dsc_aclsum
    procedure, pass(a) :: reallocate_nz => psb_d_dsc_reallocate_nz
    procedure, pass(a) :: allocate_mnnz => psb_d_dsc_allocate_mnnz
    procedure, pass(a) :: cp_to_coo    => psb_d_cp_dsc_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_d_cp_dsc_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_d_cp_dsc_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_d_cp_dsc_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_d_mv_dsc_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_d_mv_dsc_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_d_mv_dsc_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_d_mv_dsc_from_fmt
    procedure, pass(a) :: csput        => psb_d_dsc_csput
    procedure, pass(a) :: get_diag     => psb_d_dsc_get_diag
    procedure, pass(a) :: csgetptn     => psb_d_dsc_csgetptn
    procedure, pass(a) :: d_csgetrow   => psb_d_dsc_csgetrow
    procedure, pass(a) :: get_nz_col   => d_dsc_get_nz_col
    procedure, pass(a) :: reinit       => psb_d_dsc_reinit
    procedure, pass(a) :: trim         => psb_d_dsc_trim
    procedure, pass(a) :: print        => psb_d_dsc_print
    procedure, pass(a) :: free         => d_dsc_free
    procedure, pass(a) :: mold         => psb_d_dsc_mold
    procedure, pass(a) :: psb_d_dsc_cp_from
    generic, public    :: cp_from => psb_d_dsc_cp_from
    procedure, pass(a) :: psb_d_dsc_mv_from
    generic, public    :: mv_from => psb_d_dsc_mv_from

  end type psb_d_dsc_sparse_mat

 private :: d_dsc_get_nzeros, d_dsc_free,  d_dsc_get_fmt, &
       & d_dsc_get_size, d_dsc_sizeof, d_dsc_get_nz_col

  interface
    subroutine  psb_d_dsc_reallocate_nz(nz,a) 
      import :: psb_d_dsc_sparse_mat
      integer, intent(in) :: nz
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
    end subroutine psb_d_dsc_reallocate_nz
  end interface
  
  interface 
    subroutine psb_d_dsc_reinit(a,clear)
      import :: psb_d_dsc_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_d_dsc_reinit
  end interface
  
  interface
    subroutine  psb_d_dsc_trim(a)
      import :: psb_d_dsc_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
    end subroutine psb_d_dsc_trim
  end interface
  
  interface
    subroutine  psb_d_dsc_allocate_mnnz(m,n,a,nz) 
      import :: psb_d_dsc_sparse_mat
      integer, intent(in) :: m,n
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      integer, intent(in), optional :: nz
    end subroutine psb_d_dsc_allocate_mnnz
  end interface

  interface 
    subroutine psb_d_dsc_mold(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_base_sparse_mat, psb_long_int_k_
      class(psb_d_dsc_sparse_mat), intent(in)               :: a
      class(psb_d_base_sparse_mat), intent(out), allocatable :: b
      integer, intent(out)                                 :: info
    end subroutine psb_d_dsc_mold
  end interface

  interface
    subroutine psb_d_dsc_print(iout,a,iv,head,ivr,ivc)
      import :: psb_d_dsc_sparse_mat
      integer, intent(in)               :: iout
      class(psb_d_dsc_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_d_dsc_print
  end interface
  
  interface 
    subroutine psb_d_cp_dsc_to_coo(a,b,info) 
      import :: psb_d_coo_sparse_mat, psb_d_dsc_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)            :: info
    end subroutine psb_d_cp_dsc_to_coo
  end interface
  
  interface 
    subroutine psb_d_cp_dsc_from_coo(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_cp_dsc_from_coo
  end interface
  
  interface 
    subroutine psb_d_cp_dsc_to_fmt(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(in)   :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                       :: info
    end subroutine psb_d_cp_dsc_to_fmt
  end interface
  
  interface 
    subroutine psb_d_cp_dsc_from_fmt(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(in)   :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_cp_dsc_from_fmt
  end interface
  
  interface 
    subroutine psb_d_mv_dsc_to_coo(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout)   :: b
      integer, intent(out)            :: info
    end subroutine psb_d_mv_dsc_to_coo
  end interface
  
  interface 
    subroutine psb_d_mv_dsc_from_coo(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_coo_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_mv_dsc_from_coo
  end interface
  
  interface 
    subroutine psb_d_mv_dsc_to_fmt(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(inout)  :: b
      integer, intent(out)                        :: info
    end subroutine psb_d_mv_dsc_to_fmt
  end interface
  
  interface 
    subroutine psb_d_mv_dsc_from_fmt(a,b,info) 
      import :: psb_d_dsc_sparse_mat, psb_d_base_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(inout)  :: a
      class(psb_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine psb_d_mv_dsc_from_fmt
  end interface
  
  interface 
    subroutine psb_d_dsc_cp_from(a,b)
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      type(psb_d_dsc_sparse_mat), intent(in)   :: b
    end subroutine psb_d_dsc_cp_from
  end interface
  
  interface 
    subroutine psb_d_dsc_mv_from(a,b)
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(inout)  :: a
      type(psb_d_dsc_sparse_mat), intent(inout) :: b
    end subroutine psb_d_dsc_mv_from
  end interface
  
  
  interface 
    subroutine psb_d_dsc_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz,ia(:), ja(:),&
           &  imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psb_d_dsc_csput
  end interface
  
  interface 
    subroutine psb_d_dsc_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_dsc_csgetptn
  end interface
  
  interface 
    subroutine psb_d_dsc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_dsc_csgetrow
  end interface

  interface 
    subroutine psb_d_dsc_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
      import :: psb_d_dsc_sparse_mat, psb_dpk_, psb_d_coo_sparse_mat
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      class(psb_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(in)                  :: imin,imax
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_d_dsc_csgetblk
  end interface
    
  interface 
    subroutine psb_d_dsc_cssv(alpha,a,x,beta,y,info,trans) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_dsc_cssv
    subroutine psb_d_dsc_cssm(alpha,a,x,beta,y,info,trans) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_dsc_cssm
  end interface
  
  interface 
    subroutine psb_d_dsc_csmv(alpha,a,x,beta,y,info,trans) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_dsc_csmv
    subroutine psb_d_dsc_csmm(alpha,a,x,beta,y,info,trans) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine psb_d_dsc_csmm
  end interface
  
  
  interface 
    function psb_d_dsc_csnmi(a) result(res)
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_dsc_csnmi
  end interface
  
  interface 
    function psb_d_dsc_csnm1(a) result(res)
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function psb_d_dsc_csnm1
  end interface

  interface 
    subroutine psb_d_dsc_rowsum(d,a) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_dsc_rowsum
  end interface

  interface 
    subroutine psb_d_dsc_arwsum(d,a) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_dsc_arwsum
  end interface
  
  interface 
    subroutine psb_d_dsc_colsum(d,a) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_dsc_colsum
  end interface

  interface 
    subroutine psb_d_dsc_aclsum(d,a) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)              :: d(:)
    end subroutine psb_d_dsc_aclsum
  end interface
    
  interface 
    subroutine psb_d_dsc_get_diag(a,d,info) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
    end subroutine psb_d_dsc_get_diag
  end interface
  
  interface 
    subroutine psb_d_dsc_scal(d,a,info,side) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d(:)
      integer, intent(out)            :: info
      character, intent(in), optional :: side
    end subroutine psb_d_dsc_scal
  end interface
  
  interface
    subroutine psb_d_dsc_scals(d,a,info) 
      import :: psb_d_dsc_sparse_mat, psb_dpk_
      class(psb_d_dsc_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: d
      integer, intent(out)            :: info
    end subroutine psb_d_dsc_scals
  end interface
  

contains 

  ! == ===================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  ! == ===================================

  
  function d_dsc_sizeof(a) result(res)
    implicit none 
    class(psb_d_dsc_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    integer :: i

    res = 12
    
    if (allocated(a%cols)) then     
      do i=1, size(a%cols)
        res = res + psb_sizeof_int + psb_sizeof_dp  * size(a%cols(i)%val) &
             &  + psb_sizeof_int * size(a%cols(i)%idx)
      end do
    end if
    
  end function d_dsc_sizeof

  function d_dsc_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'DSC'
  end function d_dsc_get_fmt
  
  function d_dsc_get_nzeros(a) result(res)
    implicit none 
    class(psb_d_dsc_sparse_mat), intent(in) :: a
    integer :: res
    if (allocated(a%cols)) then 
      res = sum(a%cols(:)%nz)
    end if
  end function d_dsc_get_nzeros

  function d_dsc_get_size(a) result(res)
    implicit none 
    class(psb_d_dsc_sparse_mat), intent(in) :: a
    integer :: res
    integer :: i, szt

    res = 0
    if (allocated(a%cols)) then 
      do i=1, size(a%cols)
        if (allocated(a%cols(i)%idx).and.allocated(a%cols(i)%val)) then 
          res = res + min(size(a%cols(i)%idx),size(a%cols(i)%val))
        end if
      end do
    end if

  end function d_dsc_get_size



  function  d_dsc_get_nz_col(idx,a) result(res)
    use psb_const_mod
    implicit none
    
    class(psb_d_dsc_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_ncols()).and.allocated(a%cols)) &
         & res = a%cols(idx)%nz
    
  end function d_dsc_get_nz_col



  ! == ===================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  ! == ===================================  


  subroutine  d_dsc_free(a) 
    implicit none 

    class(psb_d_dsc_sparse_mat), intent(inout) :: a
    integer :: info
    
    if (allocated(a%cols)) &
         & deallocate(a%cols,stat=info) 
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return

  end subroutine d_dsc_free

end module psb_d_dsc_mat_mod

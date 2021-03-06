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
!
!
!
!

module mld_d_invt_solver

  use mld_d_base_solver_mod
  use mld_d_base_ainv_mod 
  use psb_base_mod, only : psb_d_vect_type

  type, extends(mld_d_base_ainv_solver_type) :: mld_d_invt_solver_type
    integer(psb_ipk_)   :: fill_in, inv_fill
    real(psb_dpk_)      :: thresh, inv_thresh
  contains
    procedure, pass(sv) :: check   => mld_d_invt_solver_check
    procedure, pass(sv) :: clone   => mld_d_invt_solver_clone
    procedure, pass(sv) :: clone_settings => mld_d_invt_solver_clone_settings
    procedure, pass(sv) :: build   => mld_d_invt_solver_bld
    procedure, pass(sv) :: cseti   => mld_d_invt_solver_cseti
    procedure, pass(sv) :: csetr   => mld_d_invt_solver_csetr
    procedure, pass(sv) :: descr   => mld_d_invt_solver_descr
    procedure, pass(sv) :: default => d_invt_solver_default
  end type mld_d_invt_solver_type


  private :: d_invt_solver_default


  interface 
    subroutine mld_d_invt_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
       & mld_d_base_solver_type, psb_dpk_, mld_d_invt_solver_type, psb_ipk_
      Implicit None
      class(mld_d_invt_solver_type), intent(inout)              :: sv
      class(mld_d_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)                            :: info
    end subroutine mld_d_invt_solver_clone
  end interface

  interface 
    subroutine mld_d_invt_solver_clone_settings(sv,svout,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
       & mld_d_base_solver_type, psb_dpk_, mld_d_invt_solver_type, psb_ipk_
      Implicit None
      class(mld_d_invt_solver_type), intent(inout) :: sv
      class(mld_d_base_solver_type), intent(inout) :: svout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_d_invt_solver_clone_settings
  end interface
  
  interface 
    subroutine mld_d_invt_solver_bld(a,desc_a,sv,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
       & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_,&
       & mld_d_invt_solver_type, psb_i_base_vect_type, psb_ipk_
      
      Implicit None
      
      ! Arguments
      type(psb_dspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(inout)                  :: desc_a 
      class(mld_d_invt_solver_type), intent(inout)        :: sv
      integer(psb_ipk_), intent(out)                      :: info
      type(psb_dspmat_type), intent(in), target, optional :: b
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_d_invt_solver_bld
  end interface
  
  interface 
    subroutine mld_d_invt_solver_check(sv,info)
      import :: mld_d_invt_solver_type, psb_ipk_

      Implicit None
      
      ! Arguments
      class(mld_d_invt_solver_type), intent(inout) :: sv
      integer(psb_ipk_), intent(out)               :: info
    end subroutine mld_d_invt_solver_check
  end interface
  
  interface 
    subroutine mld_d_invt_solver_cseti(sv,what,val,info,idx)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat,  psb_ipk_,&
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_invt_solver_type
      Implicit None
      ! Arguments
      class(mld_d_invt_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                 :: what 
      integer(psb_ipk_), intent(in)                :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine mld_d_invt_solver_cseti
  end interface
  
  interface 
    subroutine mld_d_invt_solver_csetr(sv,what,val,info,idx)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat,  psb_ipk_,&
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_invt_solver_type
      Implicit None
      ! Arguments
      class(mld_d_invt_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                 :: what 
      real(psb_dpk_), intent(in)                   :: val
      integer(psb_ipk_), intent(out)               :: info
      integer(psb_ipk_), intent(in), optional      :: idx
    end subroutine mld_d_invt_solver_csetr
  end interface 

  interface
    subroutine mld_d_invt_solver_descr(sv,info,iout,coarse)
      import :: psb_dpk_, mld_d_invt_solver_type, psb_ipk_
      
      Implicit None
      
      ! Arguments
      class(mld_d_invt_solver_type), intent(in) :: sv
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), intent(in), optional   :: iout
      logical, intent(in), optional             :: coarse

    end subroutine mld_d_invt_solver_descr
  end interface 

  interface mld_invt_bld
    subroutine mld_d_invt_bld(a,fillin,invfill,thresh,invthresh,&
         & lmat,d,umat,desc,info,blck)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_dpk_, psb_ipk_      
      
      implicit none
      
      ! Arguments                                                     
      type(psb_dspmat_type), intent(in), target   :: a
      integer(psb_ipk_), intent(in)               :: fillin,invfill
      real(psb_dpk_), intent(in)                  :: thresh
      real(psb_dpk_), intent(in)                  :: invthresh
      type(psb_dspmat_type), intent(inout)        :: lmat, umat
      real(psb_dpk_), allocatable                 :: d(:)
      Type(psb_desc_type), Intent(inout)          :: desc
      integer(psb_ipk_), intent(out)              :: info
      type(psb_dspmat_type), intent(in), optional :: blck
    end subroutine mld_d_invt_bld
  end interface
  

  interface mld_invt_copyin
    subroutine mld_d_invt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
         & irwt,ktrw,trw,info,sign)
      use psb_base_mod, only : psb_d_csr_sparse_mat, psb_d_coo_sparse_mat,&
           & psb_dpk_, psb_i_heap, psb_ipk_
      implicit none 
      type(psb_d_csr_sparse_mat), intent(in)    :: a
      type(psb_d_coo_sparse_mat), intent(inout) :: trw
      integer(psb_ipk_), intent(in)        :: i, m,jmin,jmax,jd
      integer(psb_ipk_), intent(inout)     :: ktrw,nlw,nup,jmaxup,info
      integer(psb_ipk_), intent(inout)     :: irwt(:)
      real(psb_dpk_), intent(inout)        :: nrmi,row(:)
      type(psb_i_heap), intent(inout)      :: heap
      real(psb_dpk_), intent(in), optional :: sign

    end subroutine mld_d_invt_copyin
  end interface

  interface mld_invt_inv
    subroutine mld_d_invt_inv(thres,i,nrmi,row,heap,irwt,uia1,uia2,uaspk,&
         & nidx,idxs,info)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_i_heap, psb_ipk_
      implicit none 
      ! Arguments
      type(psb_i_heap), intent(inout)     :: heap 
      integer(psb_ipk_), intent(in)       :: i
      integer(psb_ipk_), intent(inout)    :: nidx,info
      integer(psb_ipk_), intent(inout)    :: irwt(:) 
      real(psb_dpk_), intent(in)          :: thres,nrmi
      integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
      integer(psb_ipk_), intent(in)       :: uia1(:),uia2(:)
      real(psb_dpk_), intent(in)          :: uaspk(:)
      real(psb_dpk_), intent(inout)       :: row(:)

    end subroutine mld_d_invt_inv
  end interface

  interface mld_invt_copyout
    subroutine mld_d_invt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
         & nidx,idxs,l2,uia1,uia2,uaspk,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_i_heap, psb_ipk_

      implicit none 

      ! Arguments
      integer(psb_ipk_), intent(in)                 :: fill_in,i,m,nidx,nlw,nup,jmaxup
      integer(psb_ipk_), intent(in)                 :: idxs(:)
      integer(psb_ipk_), intent(inout)              :: l2, info
      integer(psb_ipk_), allocatable, intent(inout) :: uia1(:),uia2(:)
      real(psb_dpk_), intent(in)                    :: thres,nrmi
      real(psb_dpk_),allocatable, intent(inout)     :: uaspk(:)
      real(psb_dpk_), intent(inout)                 :: row(:)

    end subroutine mld_d_invt_copyout
  end interface

  interface  mld_sparse_invt
    subroutine mld_dsparse_invt(n,a,z,nzrmax,sp_thresh,info)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_ipk_
      implicit none 
      integer(psb_ipk_), intent(in)        :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer(psb_ipk_), intent(in)        :: nzrmax
      real(psb_dpk_), intent(in)           :: sp_thresh
      integer(psb_ipk_), intent(out)       :: info
    end subroutine mld_dsparse_invt
  end interface

  
contains

  subroutine d_invt_solver_default(sv)

    !use psb_base_mod
    
    Implicit None

    ! Arguments
    class(mld_d_invt_solver_type), intent(inout) :: sv

    sv%fill_in    = 0
    sv%inv_fill   = 0
    sv%thresh     = dzero
    sv%inv_thresh = dzero

    return
  end subroutine d_invt_solver_default

end module mld_d_invt_solver

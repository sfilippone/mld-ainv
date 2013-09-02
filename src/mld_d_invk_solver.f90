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
!
!
!
!
!

module mld_d_invk_solver

  use mld_d_base_solver_mod
  use mld_d_base_ainv_mod 
  use psb_base_mod, only : psb_d_vect_type

  type, extends(mld_d_base_ainv_solver_type) :: mld_d_invk_solver_type
    integer                     :: fill_in, inv_fill
  contains
    procedure, pass(sv) :: clone   => mld_d_invk_solver_clone
    procedure, pass(sv) :: build   => mld_d_invk_solver_bld
    procedure, pass(sv) :: seti    => mld_d_invk_solver_seti
    procedure, pass(sv) :: cseti   => mld_d_invk_solver_cseti
    procedure, pass(sv) :: descr   => mld_d_invk_solver_descr
    procedure, pass(sv) :: default => d_invk_solver_default
  end type mld_d_invk_solver_type


  private ::  d_invk_solver_default


  interface 
    subroutine mld_d_invk_solver_clone(sv,svout,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
       & mld_d_base_solver_type, psb_dpk_, mld_d_invk_solver_type, psb_ipk_
      Implicit None
      class(mld_d_invk_solver_type), intent(inout)              :: sv
      class(mld_d_base_solver_type), allocatable, intent(inout) :: svout
      integer(psb_ipk_), intent(out)                            :: info
    end subroutine mld_d_invk_solver_clone
  end interface
  
  interface 
    subroutine mld_d_invk_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold,imold)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
       & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, &
       & mld_d_invk_solver_type, psb_i_base_vect_type
      
      Implicit None
      
      ! Arguments
      type(psb_dspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(in)                     :: desc_a 
      class(mld_d_invk_solver_type), intent(inout)        :: sv
      character, intent(in)                               :: upd
      integer, intent(out)                                :: info
      type(psb_dspmat_type), intent(in), target, optional :: b
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_d_invk_solver_bld
  end interface
  
  interface 
    subroutine mld_d_invk_solver_check(sv,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_invk_solver_type

      Implicit None
      
      ! Arguments
      class(mld_d_invk_solver_type), intent(inout) :: sv
      integer, intent(out)                   :: info
    end subroutine mld_d_invk_solver_check
  end interface
  
  interface 
    subroutine mld_d_invk_solver_seti(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_invk_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_invk_solver_type), intent(inout) :: sv 
      integer, intent(in)                          :: what 
      integer, intent(in)                          :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_invk_solver_seti
  end interface
  
  interface 
    subroutine mld_d_invk_solver_cseti(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_invk_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_invk_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                 :: what 
      integer, intent(in)                          :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_invk_solver_cseti
  end interface
  
  interface
    subroutine mld_d_invk_solver_descr(sv,info,iout,coarse)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_invk_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_invk_solver_type), intent(in) :: sv
      integer, intent(out)                      :: info
      integer, intent(in), optional             :: iout
      logical, intent(in), optional             :: coarse

    end subroutine mld_d_invk_solver_descr
  end interface 
  

  interface mld_invk_bld
    subroutine mld_d_invk_bld(a,fill1, fill2,lmat,d,umat,desc,info,blck)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_dpk_ 
      
      implicit none
      
      ! Arguments                                                     
      type(psb_dspmat_type), intent(in), target   :: a
      integer, intent(in)                         :: fill1, fill2 
      type(psb_dspmat_type), intent(inout)        :: lmat, umat
      real(psb_dpk_), allocatable                 :: d(:)
      Type(psb_desc_type), Intent(in)             :: desc
      integer, intent(out)                        :: info
      type(psb_dspmat_type), intent(in), optional :: blck
    end subroutine mld_d_invk_bld
  end interface
  
  interface  mld_invk_copyin
    subroutine mld_d_invk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,&
         & ktrw,trw,info,sign,inlevs)

      use psb_base_mod, only : psb_d_csr_sparse_mat, psb_d_coo_sparse_mat,&
           & psb_dpk_, psb_int_heap
      implicit none

      ! Arguments 
      type(psb_d_csr_sparse_mat), intent(in)    :: a
      type(psb_d_coo_sparse_mat), intent(inout) :: trw
      integer, intent(in)                  :: i,m,jmin,jmax
      integer, intent(inout)               :: ktrw,info
      integer, intent(inout)               :: rowlevs(:)
      real(psb_dpk_), intent(inout)        :: row(:)
      type(psb_int_heap), intent(inout)    :: heap
      real(psb_dpk_), optional, intent(in) :: sign
      integer, intent(in), optional        :: inlevs(:)

    end subroutine mld_d_invk_copyin
  end interface

  interface mld_invk_inv
    subroutine mld_d_invk_inv(fill_in,i,row,rowlevs,heap,uia1,uia2,uaspk,uplevs,&
         & nidx,idxs,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
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


    end subroutine mld_d_invk_inv
  end interface

  interface mld_invk_copyout
    subroutine mld_d_invk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
         &  l2,uia1,uia2,uaspk,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap

      implicit none 

      ! Arguments
      integer, intent(in)                        :: fill_in, i, m, nidx
      integer, intent(inout)                     :: l2, info
      integer, intent(inout)                     :: rowlevs(:), idxs(:)
      integer, allocatable, intent(inout)        :: uia1(:), uia2(:)
      real(psb_dpk_), allocatable, intent(inout) :: uaspk(:)
      real(psb_dpk_), intent(inout)              :: row(:)

    end subroutine mld_d_invk_copyout
  end interface
  
  interface mld_sparse_invk
    subroutine mld_dsparse_invk(n,a,z,fill_in,info,inlevs)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      integer, intent(in)                  :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: fill_in
      integer, intent(out)                 :: info
      integer, intent(in), optional        :: inlevs(:)

    end subroutine mld_dsparse_invk
  end interface

  
contains

  subroutine d_invk_solver_default(sv)

    !use psb_base_mod
    
    Implicit None

    ! Arguments
    class(mld_d_invk_solver_type), intent(inout) :: sv

    sv%fill_in    = 0
    sv%inv_fill   = 0

    return
  end subroutine d_invk_solver_default

end module mld_d_invk_solver

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
module mld_d_base_ainv_mod
  
  use mld_base_ainv_mod
  use mld_d_prec_type
  use mld_d_ainv_tools_mod
  use psb_base_mod, only : psb_d_vect_type, psb_dpk_,psb_epk_

  

  type, extends(mld_d_base_solver_type) :: mld_d_base_ainv_solver_type
    ! 
    !  Compute an approximate factorization
    !      A^-1 = Z D^-1 W^T
    !  Note that here W is going to be transposed explicitly,
    !  so that the component w will in the end contain W^T.     
    !
    type(psb_dspmat_type)       :: w, z
    type(psb_d_vect_type)       :: dv
    real(psb_dpk_), allocatable :: d(:)

  contains
    procedure, pass(sv) :: cnv     => mld_d_base_ainv_solver_cnv
    procedure, pass(sv) :: dump    => mld_d_base_ainv_solver_dmp
    procedure, pass(sv) :: apply_v => mld_d_base_ainv_solver_apply_vect
    procedure, pass(sv) :: apply_a => mld_d_base_ainv_solver_apply
    procedure, pass(sv) :: free    => mld_d_base_ainv_solver_free
    procedure, pass(sv) :: sizeof  => d_base_ainv_solver_sizeof
    procedure, pass(sv) :: get_nzeros => d_base_ainv_get_nzeros
    procedure, nopass   :: get_wrksz => d_base_ainv_get_wrksize
    procedure, pass(sv) :: update_a => mld_d_base_ainv_update_a
    generic, public     :: update => update_a
  end type mld_d_base_ainv_solver_type

  private ::  d_base_ainv_solver_sizeof, &
       & d_base_ainv_get_nzeros, d_base_ainv_solver_get_wrksize

  
  interface 
    subroutine mld_d_base_ainv_solver_cnv(sv,info,amold,vmold,imold)
      import :: psb_d_base_sparse_mat, psb_d_base_vect_type, psb_dpk_, &
       & mld_d_base_ainv_solver_type, psb_ipk_, psb_i_base_vect_type      
      Implicit None      
      ! Arguments
      class(mld_d_base_ainv_solver_type), intent(inout)  :: sv
      integer(psb_ipk_), intent(out)                     :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine mld_d_base_ainv_solver_cnv
  end interface

  interface 
    subroutine mld_d_base_ainv_update_a(sv,x,desc_data,info)
      import :: psb_desc_type, psb_dpk_,mld_d_base_ainv_solver_type, psb_d_vect_type
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_d_base_ainv_solver_type), intent(inout) :: sv
      real(psb_dpk_),intent(in)            :: x(:)
      integer, intent(out)                 :: info
    end subroutine mld_d_base_ainv_update_a
  end interface

  interface 
    subroutine mld_d_base_ainv_solver_apply(alpha,sv,x,beta,y,desc_data,&
         & trans,work,info,init,initu)
      import :: psb_desc_type, psb_dpk_,mld_d_base_ainv_solver_type
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_d_base_ainv_solver_type), intent(inout) :: sv
      real(psb_dpk_),intent(inout)         :: x(:)
      real(psb_dpk_),intent(inout)         :: y(:)
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer, intent(out)                 :: info
      character, intent(in), optional       :: init
      real(psb_dpk_),intent(inout), optional :: initu(:)
    end subroutine mld_d_base_ainv_solver_apply
  end interface

  interface 
    subroutine mld_d_base_ainv_solver_apply_vect(alpha,sv,x,beta,y,desc_data,&
         &  trans,work,wv,info,init,initu)
      import :: psb_desc_type, psb_dpk_,mld_d_base_ainv_solver_type, psb_d_vect_type
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_d_base_ainv_solver_type), intent(inout) :: sv
      type(psb_d_vect_type),intent(inout)  :: x
      type(psb_d_vect_type),intent(inout)  :: y
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      type(psb_d_vect_type),intent(inout)  :: wv(:)
      integer, intent(out)                 :: info
      character, intent(in), optional                :: init
      type(psb_d_vect_type),intent(inout), optional   :: initu
    end subroutine mld_d_base_ainv_solver_apply_vect
  end interface
 
  
  interface
    subroutine mld_d_base_ainv_solver_free(sv,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_base_ainv_solver_type
      Implicit None
      
      ! Arguments
      class(mld_d_base_ainv_solver_type), intent(inout) :: sv
      integer, intent(out)                         :: info
    end subroutine mld_d_base_ainv_solver_free
  end interface
  
  
  interface 
    subroutine mld_d_base_ainv_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_base_ainv_solver_type
      
      implicit none 
      class(mld_d_base_ainv_solver_type), intent(in) :: sv
      integer, intent(in)              :: ictxt,level
      integer, intent(out)             :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: solver
    end subroutine mld_d_base_ainv_solver_dmp
  end interface



contains


  function d_base_ainv_get_nzeros(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_d_base_ainv_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer             :: i

    val = 0
    val = val + sv%dv%get_nrows()
    val = val + sv%w%get_nzeros()
    val = val + sv%z%get_nzeros()

    return
  end function d_base_ainv_get_nzeros

  function d_base_ainv_solver_sizeof(sv) result(val)
    implicit none 
    ! Arguments
    class(mld_d_base_ainv_solver_type), intent(in) :: sv
    integer(psb_epk_) :: val
    integer             :: i

    val = 2*psb_sizeof_ip + psb_sizeof_dp
    val = val + sv%dv%sizeof()
    val = val + sv%w%sizeof()
    val = val + sv%z%sizeof()
    
    return
  end function d_base_ainv_solver_sizeof

  function d_base_ainv_get_wrksize() result(val)
    implicit none 
    integer(psb_ipk_)  :: val

    val = 2
  end function d_base_ainv_get_wrksize

end module mld_d_base_ainv_mod


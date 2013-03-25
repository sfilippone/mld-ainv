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
module mld_d_ainv_solver

  use mld_d_prec_type
  use mld_base_ainv_mod 
  use psb_base_mod, only : psb_d_vect_type

  type, extends(mld_d_base_ainv_solver_type) :: mld_d_ainv_solver_type
    ! 
    !  Compute an approximate factorization
    !      A^-1 = Z D^-1 W^T
    !  Note that here W is going to be transposed explicitly,
    !  so that the component w will in the end contain W^T.     
    !
    integer                     :: alg, fill_in
    real(psb_dpk_)              :: thresh
  contains
    procedure, pass(sv) :: dump    => mld_d_ainv_solver_dmp
    procedure, pass(sv) :: build   => mld_d_ainv_solver_bld
    procedure, pass(sv) :: seti    => mld_d_ainv_solver_seti
    procedure, pass(sv) :: setc    => mld_d_ainv_solver_setc
    procedure, pass(sv) :: setr    => mld_d_ainv_solver_setr
    procedure, pass(sv) :: cseti   => mld_d_ainv_solver_cseti
    procedure, pass(sv) :: csetc   => mld_d_ainv_solver_csetc
    procedure, pass(sv) :: csetr   => mld_d_ainv_solver_csetr
    procedure, pass(sv) :: descr   => mld_d_ainv_solver_descr
    procedure, pass(sv) :: default => d_ainv_solver_default
    procedure, nopass   :: stringval  => d_ainv_stringval
    procedure, nopass   :: algname => d_ainv_algname
  end type mld_d_ainv_solver_type


  private :: d_ainv_stringval, d_ainv_solver_default, &
       &  d_ainv_algname

  interface 
    subroutine mld_d_ainv_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, psb_dpk_,mld_d_ainv_solver_type
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_d_ainv_solver_type), intent(in) :: sv
      real(psb_dpk_),intent(inout)         :: x(:)
      real(psb_dpk_),intent(inout)         :: y(:)
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer, intent(out)                 :: info
    end subroutine mld_d_ainv_solver_apply
  end interface

  interface 
    subroutine mld_d_ainv_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
      import :: psb_desc_type, psb_dpk_,mld_d_ainv_solver_type, psb_d_vect_type
      type(psb_desc_type), intent(in)      :: desc_data
      class(mld_d_ainv_solver_type), intent(inout) :: sv
      type(psb_d_vect_type),intent(inout)  :: x
      type(psb_d_vect_type),intent(inout)  :: y
      real(psb_dpk_),intent(in)            :: alpha,beta
      character(len=1),intent(in)          :: trans
      real(psb_dpk_),target, intent(inout) :: work(:)
      integer, intent(out)                 :: info
    end subroutine mld_d_ainv_solver_apply_vect
  end interface
 
  
  interface 
    subroutine mld_d_ainv_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
       & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      
      Implicit None
      
      ! Arguments
      type(psb_dspmat_type), intent(in), target           :: a
      Type(psb_desc_type), Intent(in)                     :: desc_a 
      class(mld_d_ainv_solver_type), intent(inout)        :: sv
      character, intent(in)                               :: upd
      integer, intent(out)                                :: info
      type(psb_dspmat_type), intent(in), target, optional :: b
      class(psb_d_base_sparse_mat), intent(in), optional  :: amold
      class(psb_d_base_vect_type), intent(in), optional   :: vmold
    end subroutine mld_d_ainv_solver_bld
  end interface
  
  interface 
    subroutine mld_d_ainv_solver_check(sv,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type

      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv
      integer, intent(out)                   :: info
    end subroutine mld_d_ainv_solver_check
  end interface
  
  interface 
    subroutine mld_d_ainv_solver_seti(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv 
      integer, intent(in)                          :: what 
      integer, intent(in)                          :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_seti
  end interface
  
  interface 
    subroutine mld_d_ainv_solver_setc(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv
      integer, intent(in)                          :: what 
      character(len=*), intent(in)                 :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_setc
  end interface 
  
  interface 
    subroutine mld_d_ainv_solver_setr(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
            
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv 
      integer, intent(in)                          :: what 
      real(psb_dpk_), intent(in)                   :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_setr
  end interface 
  
  interface 
    subroutine mld_d_ainv_solver_cseti(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                 :: what 
      integer, intent(in)                          :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_cseti
  end interface 
  
  
  interface 
    subroutine mld_d_ainv_solver_csetc(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                 :: what 
      character(len=*), intent(in)                 :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_csetc
  end interface 
  
  interface 
    subroutine mld_d_ainv_solver_csetr(sv,what,val,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
            
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv 
      character(len=*), intent(in)                 :: what 
      real(psb_dpk_), intent(in)                   :: val
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_csetr
  end interface 
 
  
  interface
    subroutine mld_d_ainv_solver_free(sv,info)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(inout) :: sv
      integer, intent(out)                         :: info
    end subroutine mld_d_ainv_solver_free
  end interface
  
  interface
    subroutine mld_d_ainv_solver_descr(sv,info,iout,coarse)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      
      Implicit None
      
      ! Arguments
      class(mld_d_ainv_solver_type), intent(in) :: sv
      integer, intent(out)                      :: info
      integer, intent(in), optional             :: iout
      logical, intent(in), optional             :: coarse

    end subroutine mld_d_ainv_solver_descr
  end interface 
  
  interface 
    subroutine mld_d_ainv_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_, mld_d_ainv_solver_type
      
      implicit none 
      class(mld_d_ainv_solver_type), intent(in) :: sv
      integer, intent(in)              :: ictxt,level
      integer, intent(out)             :: info
      character(len=*), intent(in), optional :: prefix, head
      logical, optional, intent(in)    :: solver
    end subroutine mld_d_ainv_solver_dmp
  end interface
  
  interface  mld_ainv_bld
    subroutine mld_d_ainv_bld(a,alg,fillin,thresh,wmat,d,zmat,desc,info,blck,iscale)
      import :: psb_desc_type, psb_dspmat_type,  psb_d_base_sparse_mat, &
           & psb_d_vect_type, psb_d_base_vect_type, psb_dpk_
      implicit none
      type(psb_dspmat_type), intent(in), target   :: a
      integer, intent(in)                         :: fillin,alg
      real(psb_dpk_), intent(in)                  :: thresh
      type(psb_dspmat_type), intent(inout)        :: wmat, zmat
      real(psb_dpk_), allocatable                 :: d(:)
      Type(psb_desc_type), Intent(in)             :: desc
      integer, intent(out)                        :: info
      type(psb_dspmat_type), intent(in), optional :: blck
      integer, intent(in), optional               :: iscale
    end subroutine mld_d_ainv_bld
  end interface


contains

  subroutine d_ainv_solver_default(sv)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainv_solver_type), intent(inout) :: sv

    sv%alg     = mld_ainv_llk_
    sv%fill_in = 0
    sv%thresh  = dzero

    return
  end subroutine d_ainv_solver_default

  function d_ainv_stringval(string) result(val)
    use psb_base_mod, only : psb_ipk_,psb_toupper
    implicit none 
  ! Arguments
    character(len=*), intent(in) :: string
    integer(psb_ipk_) :: val 
    character(len=*), parameter :: name='d_ainv_stringval'
    
    select case(psb_toupper(trim(string)))
    case('LLK')
      val = mld_ainv_llk_
    case('S-LLK')
      val = mld_ainv_s_llk_
#if defined(HAVE_TUMA_SAINV)
    case('SAINV-TUMA')
      val = mld_ainv_s_tuma_ 
#endif 
    case default
      val  = mld_stringval(string)
    end select
  end function d_ainv_stringval


  function d_ainv_algname(ialg) result(val)
    integer(psb_ipk_), intent(in) :: ialg 
    character(len=40) :: val
    
    character(len=*), parameter :: llkname   = 'Left-looking '
    character(len=*), parameter :: sllkname  = 'Symmetric Left-looking '
    character(len=*), parameter :: sainvname = 'SAINV (Benzi & Tuma) '
    character(len=*), parameter :: defname   = 'Unknown alg variant '
    
    select case (ialg)
    case(mld_ainv_llk_)
      val = llkname
    case(mld_ainv_s_llk_)
      val = sllkname
#if defined(HAVE_TUMA_SAINV)
    case(mld_ainv_s_tuma_ )
      val = sainvname
#endif 
    case default
      val = defname
    end select
    
  end function d_ainv_algname

end module mld_d_ainv_solver

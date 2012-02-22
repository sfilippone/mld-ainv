!!$
!!$ 
!!$                           MLD2P4  version 2.0
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 3.0)
!!$  
!!$  (C) Copyright 2008,2009,2010, 2010
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      CNRS-IRIT, Toulouse
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
!
!
!
!
!
!

module mld_d_ainvk_solver

  use mld_d_prec_type
  use mld_base_ainv_mod
  use psb_base_mod, only : psb_d_vect_type

  type, extends(mld_d_base_solver_type) :: mld_d_ainvk_solver_type
    type(psb_dspmat_type)       :: l, u
    type(psb_d_vect_type)       :: dv
    real(psb_dpk_), allocatable :: d(:)
    integer                     :: fill_in, inv_fill
    real(psb_dpk_)              :: thresh
  contains
    procedure, pass(sv) :: dump    => d_ainvk_solver_dmp
    procedure, pass(sv) :: build   => d_ainvk_solver_bld
    procedure, pass(sv) :: apply_v => d_ainvk_solver_apply_vect
    procedure, pass(sv) :: apply_a => d_ainvk_solver_apply
    procedure, pass(sv) :: free    => d_ainvk_solver_free
    procedure, pass(sv) :: seti    => d_ainvk_solver_seti
    procedure, pass(sv) :: setc    => d_ainvk_solver_setc
    procedure, pass(sv) :: setr    => d_ainvk_solver_setr
    procedure, pass(sv) :: descr   => d_ainvk_solver_descr
    procedure, pass(sv) :: sizeof  => d_ainvk_solver_sizeof
    procedure, pass(sv) :: default => d_ainvk_solver_default
    procedure, pass(sv) :: get_nzeros => d_ainvk_get_nzeros
  end type mld_d_ainvk_solver_type


  private :: d_ainvk_solver_bld, d_ainvk_solver_apply, &
       &  d_ainvk_solver_free,   d_ainvk_solver_seti, &
       &  d_ainvk_solver_setc,   d_ainvk_solver_setr,&
       &  d_ainvk_solver_descr,  d_ainvk_solver_sizeof, &
       &  d_ainvk_solver_default, d_ainvk_solver_dmp,&
       &  d_ainvk_solver_apply_vect,  d_ainvk_get_nzeros


contains

  subroutine d_ainvk_solver_default(sv)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(inout) :: sv

    sv%fill_in   = 0
    sv%inv_fill  = 0
    sv%thresh    = dzero

    return
  end subroutine d_ainvk_solver_default

  subroutine d_ainvk_solver_check(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(inout) :: sv
    integer, intent(out)                   :: info
    Integer           :: err_act
    character(len=20) :: name='d_ainvk_solver_check'

    call psb_erractionsave(err_act)
    info = psb_success_

    call mld_check_def(sv%fill_in,&
         & 'Level',0,is_legal_ml_lev)
    call mld_check_def(sv%inv_fill,&
         & 'Level',0,is_legal_ml_lev)
    call mld_check_def(sv%thresh,&
         & 'Eps',dzero,is_legal_d_fact_thrs)
    
    if (info /= psb_success_) goto 9999
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_check


  subroutine d_ainvk_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_d_ainvk_solver_type), intent(in) :: sv
    real(psb_dpk_),intent(inout)         :: x(:)
    real(psb_dpk_),intent(inout)         :: y(:)
    real(psb_dpk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    real(psb_dpk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row,n_col
    real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='d_ainvk_solver_apply'

    call psb_erractionsave(err_act)

    info = psb_success_

    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N')
    case('T','C')
    case default
      call psb_errpush(psb_err_iarg_invalid_i_,name)
      goto 9999
    end select

    n_row = psb_cd_get_local_rows(desc_data)
    n_col = psb_cd_get_local_cols(desc_data)

    if (n_col <= size(work)) then 
      ww => work(1:n_col)
      if ((4*n_col+n_col) <= size(work)) then 
        aux => work(n_col+1:)
      else
        allocate(aux(4*n_col),stat=info)
        if (info /= psb_success_) then 
          info=psb_err_alloc_request_
          call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
               & a_err='real(psb_dpk_)')
          goto 9999      
        end if
      endif
    else
      allocate(ww(n_col),aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if
    endif

    select case(trans_)
    case('N')
      call psb_spmm(done,sv%l,x,dzero,ww,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)
      ww(1:n_row) = ww(1:n_row) * sv%d(1:n_row)
      if (info == psb_success_) &
           & call psb_spmm(alpha,sv%u,ww,beta,y,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)
      
    case('T','C')
      call psb_spmm(done,sv%u,x,dzero,ww,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)
      ww(1:n_row) = ww(1:n_row) * sv%d(1:n_row)
      if (info == psb_success_) &
           & call psb_spmm(alpha,sv%l,ww,beta,y,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)

    case default
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Invalid TRANS in AINVK subsolve')
      goto 9999
    end select


    if (info /= psb_success_) then

      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error in subsolve')
      goto 9999
    endif

    if (n_col <= size(work)) then 
      if ((4*n_col+n_col) <= size(work)) then 
      else
        deallocate(aux,stat=info)
      endif
    else
      deallocate(ww,aux,stat=info)
    endif

    if (info /= psb_success_) then

      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Deallocate')
      goto 9999
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_ainvk_solver_apply


  subroutine d_ainvk_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
    use psb_base_mod
    type(psb_desc_type), intent(in)      :: desc_data
    class(mld_d_ainvk_solver_type), intent(inout) :: sv
    type(psb_d_vect_type), intent(inout) :: x
    type(psb_d_vect_type), intent(inout) :: y
    real(psb_dpk_),intent(in)            :: alpha,beta
    character(len=1),intent(in)          :: trans
    real(psb_dpk_),target, intent(inout) :: work(:)
    integer, intent(out)                 :: info

    integer    :: n_row,n_col
    real(psb_dpk_), pointer :: ww(:), aux(:)
    type(psb_d_vect_type)   :: tx,ty
    integer    :: ictxt,np,me,i, err_act
    character          :: trans_
    character(len=20)  :: name='d_ainvk_solver_apply'

    call psb_erractionsave(err_act)

    info = psb_success_

    trans_ = psb_toupper(trans)
    select case(trans_)
    case('N')
    case('T','C')
    case default
      call psb_errpush(psb_err_iarg_invalid_i_,name)
      goto 9999
    end select

    n_row = psb_cd_get_local_rows(desc_data)
    n_col = psb_cd_get_local_cols(desc_data)

    if (n_col <= size(work)) then 
      ww => work(1:n_col)
      if ((4*n_col+n_col) <= size(work)) then 
        aux => work(n_col+1:)
      else
        allocate(aux(4*n_col),stat=info)
        if (info /= psb_success_) then 
          info=psb_err_alloc_request_
          call psb_errpush(info,name,i_err=(/4*n_col,0,0,0,0/),&
               & a_err='real(psb_dpk_)')
          goto 9999      
        end if
      endif
    else
      allocate(ww(n_col),aux(4*n_col),stat=info)
      if (info /= psb_success_) then 
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/5*n_col,0,0,0,0/),&
             & a_err='real(psb_dpk_)')
        goto 9999      
      end if
    endif
    
    call tx%bld(x%get_nrows(),mold=x%v)
    call ty%bld(x%get_nrows(),mold=x%v)

    select case(trans_)
    case('N')
      call psb_spmm(done,sv%l,x,dzero,tx,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)
      if (info == psb_success_) call ty%mlt(done,sv%dv,tx,dzero,info)
      if (info == psb_success_) &
           & call psb_spmm(alpha,sv%u,ty,beta,y,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)
      
    case('T','C')
      call psb_spmm(done,sv%u,x,dzero,tx,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)
      if (info == psb_success_) call ty%mlt(done,sv%dv,tx,dzero,info)
      if (info == psb_success_) &
           & call psb_spmm(alpha,sv%l,ty,beta,y,desc_data,info,&
           & trans=trans_,work=aux,doswap=.false.)

    case default
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Invalid TRANS in AINVK subsolve')
      goto 9999
    end select


    if (info /= psb_success_) then

      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error in subsolve')
      goto 9999
    endif


    call tx%free(info) 
    if (info == psb_success_) call ty%free(info)
    if (n_col <= size(work)) then 
      if ((4*n_col+n_col) <= size(work)) then 
      else
        deallocate(aux,stat=info)
      endif
    else
      deallocate(ww,aux,stat=info)
    endif

    if (info /= psb_success_) then

      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Deallocate')
      goto 9999
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_ainvk_solver_apply_vect


  subroutine d_ainvk_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)

    use psb_base_mod
    use mld_d_ainv_bld_mod
    Implicit None

    ! Arguments
    type(psb_dspmat_type), intent(in), target  :: a
    Type(psb_desc_type), Intent(in)            :: desc_a
    class(mld_d_ainvk_solver_type), intent(inout) :: sv
    character, intent(in)                      :: upd
    integer, intent(out)                       :: info
    type(psb_dspmat_type), intent(in), target, optional :: b
    class(psb_d_base_sparse_mat), intent(in), optional  :: amold
    class(psb_d_base_vect_type), intent(in), optional   :: vmold

    ! Local variables
    integer :: n_row,n_col, nrow_a, nztota
    real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
    integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
    character(len=20)  :: name='d_ainvk_solver_bld', ch_err
    
    info=psb_success_
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'


    call mld_ainv_invk_bld(a,sv%fill_in,sv%inv_fill,sv%thresh,&
         & sv%l,sv%d,sv%u,desc_a,info,b)    
    
    
    if ((info == psb_success_) .and.present(amold)) then 
      call sv%l%cscnv(info,mold=amold)
      if (info == psb_success_) &
           & call sv%u%cscnv(info,mold=amold)
    end if

    if (info == psb_success_) then 
      call sv%dv%bld(sv%d,mold=vmold)
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' end'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_bld


  subroutine d_ainvk_solver_seti(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    integer, intent(in)                    :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='d_ainvk_solver_seti'

    info = psb_success_
    call psb_erractionsave(err_act)

    select case(what) 
    case(mld_sub_fillin_)
      sv%fill_in   = val
    case(mld_inv_fillin_)
      sv%inv_fill  = val
    case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_seti

  subroutine d_ainvk_solver_setc(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(inout) :: sv
    integer, intent(in)                    :: what 
    character(len=*), intent(in)           :: val
    integer, intent(out)                   :: info
    Integer :: err_act, ival
    character(len=20)  :: name='d_ainvk_solver_setc'

    info = psb_success_
    call psb_erractionsave(err_act)
    
    call mld_stringval(val,ival,info)

    if (info == psb_success_) call sv%set(what,ival,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_setc
  
  subroutine d_ainvk_solver_setr(sv,what,val,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(inout) :: sv 
    integer, intent(in)                    :: what 
    real(psb_dpk_), intent(in)             :: val
    integer, intent(out)                   :: info
    Integer :: err_act
    character(len=20)  :: name='d_ainvk_solver_setr'

    call psb_erractionsave(err_act)
    info = psb_success_

    select case(what)
    case(mld_sub_iluthrs_) 
      sv%thresh = val
    case default
!!$      write(0,*) name,': Error: invalid WHAT'
!!$      info = -2
!!$      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_setr

  subroutine d_ainvk_solver_free(sv,info)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(inout) :: sv
    integer, intent(out)                       :: info
    Integer :: err_act
    character(len=20)  :: name='d_ainvk_solver_free'

    call psb_erractionsave(err_act)
    info = psb_success_
    
    if (allocated(sv%d)) then 
      deallocate(sv%d,stat=info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        goto 9999 
      end if
    end if
    call sv%l%free()
    call sv%u%free()
    call sv%dv%free(info)
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_free

  subroutine d_ainvk_solver_descr(sv,info,iout,coarse)

    use psb_base_mod

    Implicit None

    ! Arguments
    class(mld_d_ainvk_solver_type), intent(in) :: sv
    integer, intent(out)                     :: info
    integer, intent(in), optional            :: iout
    logical, intent(in), optional       :: coarse

    ! Local variables
    integer      :: err_act
    integer      :: ictxt, me, np
    character(len=20), parameter :: name='mld_d_ainvk_solver_descr'
    integer :: iout_

    call psb_erractionsave(err_act)
    info = psb_success_
    if (present(iout)) then 
      iout_ = iout 
    else
      iout_ = 6
    endif
    
    write(iout_,*) '  AINVK Approximate Inverse with ILU(N) '
    write(iout_,*) '  Fill level             :',sv%fill_in
    write(iout_,*) '  Inverse fill level     :',sv%inv_fill

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_ainvk_solver_descr

  function d_ainvk_get_nzeros(sv) result(val)
    use psb_base_mod, only : psb_long_int_k_
    implicit none 
    ! Arguments
    class(mld_d_ainvk_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = 0
    val = val + sv%dv%get_nrows()
    val = val + sv%l%get_nzeros()
    val = val + sv%u%get_nzeros()

    return
  end function d_ainvk_get_nzeros


  function d_ainvk_solver_sizeof(sv) result(val)
    use psb_base_mod, only : psb_long_int_k_
    implicit none 
    ! Arguments
    class(mld_d_ainvk_solver_type), intent(in) :: sv
    integer(psb_long_int_k_) :: val
    integer             :: i

    val = 2*psb_sizeof_int + psb_sizeof_dp
    val = val + sv%dv%sizeof()
    val = val + sv%l%sizeof()
    val = val + sv%u%sizeof()

    return
  end function d_ainvk_solver_sizeof

  subroutine d_ainvk_solver_dmp(sv,ictxt,level,info,prefix,head,solver)
    use psb_base_mod
    implicit none 
    class(mld_d_ainvk_solver_type), intent(in) :: sv
    integer, intent(in)              :: ictxt,level
    integer, intent(out)             :: info
    character(len=*), intent(in), optional :: prefix, head
    logical, optional, intent(in)    :: solver
    integer :: i, j, il1, iln, lname, lev
    integer :: icontxt,iam, np
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    logical :: solver_
    !  len of prefix_ 

    info = 0

    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_ainvk_d"
    end if

    call psb_info(ictxt,iam,np)

    if (present(solver)) then 
      solver_ = solver
    else
      solver_ = .false. 
    end if
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5

    if (solver_) then 
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_lower.mtx'
      if (sv%l%is_asb()) &
           & call sv%l%print(fname,head=head)
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_diag.mtx'
      if (allocated(sv%d)) &
           & call psb_geprt(fname,sv%d,head=head)
      write(fname(lname+1:),'(a,i3.3,a)')'_l',level,'_upper.mtx'
      if (sv%u%is_asb()) &
           & call sv%u%print(fname,head=head)

    end if

  end subroutine d_ainvk_solver_dmp

end module mld_d_ainvk_solver

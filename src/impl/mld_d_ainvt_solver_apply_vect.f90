subroutine mld_d_ainvt_solver_apply_vect(alpha,sv,x,beta,y,desc_data,trans,work,info)
  
  use psb_base_mod
  use mld_d_ainvt_solver, mld_protect_name => mld_d_ainvt_solver_apply_vect
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(mld_d_ainvt_solver_type), intent(inout) :: sv
  type(psb_d_vect_type),intent(inout)  :: x
  type(psb_d_vect_type),intent(inout)  :: y
  real(psb_dpk_),intent(in)            :: alpha,beta
  character(len=1),intent(in)          :: trans
  real(psb_dpk_),target, intent(inout) :: work(:)
  integer, intent(out)                 :: info

  integer    :: n_row,n_col
  real(psb_dpk_), pointer :: ww(:), aux(:)
  type(psb_d_vect_type) :: tx,ty
  integer    :: ictxt,np,me,i, err_act
  character          :: trans_
  character(len=20)  :: name='d_ainvt_solver_apply'

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

  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

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
         & a_err='Invalid TRANS in AINVT subsolve')
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

end subroutine mld_d_ainvt_solver_apply_vect

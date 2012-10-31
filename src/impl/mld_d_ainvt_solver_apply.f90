subroutine mld_d_invt_solver_apply(alpha,sv,x,beta,y,desc_data,trans,work,info)
  
  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_d_invt_solver_apply
  implicit none 
  type(psb_desc_type), intent(in)      :: desc_data
  class(mld_d_invt_solver_type), intent(in) :: sv
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
  character(len=20)  :: name='d_invt_solver_apply'

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
         & a_err='Invalid TRANS in invt subsolve')
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

end subroutine mld_d_invt_solver_apply

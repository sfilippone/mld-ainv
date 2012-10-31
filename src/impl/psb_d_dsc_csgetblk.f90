subroutine psb_d_dsc_csgetblk(imin,imax,a,b,info,&
  
     & jmin,jmax,iren,append,rscale,dscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csgetblk
  implicit none

  class(psb_d_dsc_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(in)                  :: imin,imax
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,dscale
  Integer :: err_act, nzin, nzout
  character(len=20)  :: name='csget'
  logical :: append_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  if (present(append)) then 
!!$    append_ = append
!!$  else
!!$    append_ = .false.
!!$  endif
!!$  if (append_) then 
!!$    nzin = a%get_nzeros()
!!$  else
!!$    nzin = 0
!!$  endif
!!$
!!$  call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
!!$       & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
!!$       & nzin=nzin, rscale=rscale, dscale=dscale)
!!$
!!$  if (info /= psb_success_) goto 9999
!!$
!!$  call b%set_nzeros(nzin+nzout)
!!$  call b%fix(info)
!!$  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_csgetblk

subroutine psb_d_dsc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
  
     & jmin,jmax,iren,append,nzin,rscale,dscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csgetrow
  implicit none

  class(psb_d_dsc_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,dscale

  logical :: append_, rscale_, dscale_ 
  integer :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999

!!$  if (present(jmin)) then
!!$    jmin_ = jmin
!!$  else
!!$    jmin_ = 1
!!$  endif
!!$  if (present(jmax)) then
!!$    jmax_ = jmax
!!$  else
!!$    jmax_ = a%get_ncols()
!!$  endif
!!$
!!$  if ((imax<imin).or.(jmax_<jmin_)) return
!!$
!!$  if (present(append)) then
!!$    append_=append
!!$  else
!!$    append_=.false.
!!$  endif
!!$  if ((append_).and.(present(nzin))) then 
!!$    nzin_ = nzin
!!$  else
!!$    nzin_ = 0
!!$  endif
!!$  if (present(rscale)) then 
!!$    rscale_ = rscale
!!$  else
!!$    rscale_ = .false.
!!$  endif
!!$  if (present(dscale)) then 
!!$    dscale_ = dscale
!!$  else
!!$    dscale_ = .false.
!!$  endif
!!$  if ((rscale_.or.dscale_).and.(present(iren))) then 
!!$    info = psb_err_many_optional_arg_
!!$    call psb_errpush(info,name,a_err='iren (rscale.or.dscale)')
!!$    goto 9999
!!$  end if
!!$
!!$  call dsc_getrow(imin,imax,jmin_,jmax_,a,nz,ia,ja,val,nzin_,append_,info,&
!!$       & iren)
!!$  
!!$  if (rscale_) then 
!!$    do i=nzin_+1, nzin_+nz
!!$      ia(i) = ia(i) - imin + 1
!!$    end do
!!$  end if
!!$  if (dscale_) then 
!!$    do i=nzin_+1, nzin_+nz
!!$      ja(i) = ja(i) - jmin_ + 1
!!$    end do
!!$  end if
!!$
!!$  if (info /= psb_success_) goto 9999
!!$
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
!!$
!!$contains
!!$
!!$  subroutine dsc_getrow(imin,imax,jmin,jmax,a,nz,ia,ja,val,nzin,append,info,&
!!$       & iren)
!!$
!!$    use psb_const_mod
!!$    use psb_error_mod
!!$    use psb_realloc_mod
!!$    use psb_sort_mod
!!$    implicit none
!!$
!!$    class(psb_d_dsc_sparse_mat), intent(in)    :: a
!!$    integer                              :: imin,imax,jmin,jmax
!!$    integer, intent(out)                 :: nz
!!$    integer, allocatable, intent(inout)  :: ia(:), ja(:)
!!$    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
!!$    integer, intent(in)                  :: nzin
!!$    logical, intent(in)                  :: append
!!$    integer                              :: info
!!$    integer, optional                    :: iren(:)
!!$    integer  :: nzin_, nza, idx,i,j,k, nzt, irw, lrw, icl, lcl,isz,m
!!$    integer  :: debug_level, debug_unit
!!$    character(len=20) :: name='coo_getrow'
!!$
!!$    debug_unit  = psb_get_debug_unit()
!!$    debug_level = psb_get_debug_level()
!!$
!!$    m = a%get_nrows()
!!$    nza = a%get_nzeros()
!!$    irw = imin
!!$    lrw = min(imax,a%get_nrows())
!!$    icl = jmin
!!$    lcl = min(jmax,a%get_ncols())
!!$    if (irw<0) then 
!!$      info = psb_err_pivot_too_small_
!!$      return
!!$    end if
!!$
!!$    if (append) then 
!!$      nzin_ = nzin
!!$    else
!!$      nzin_ = 0
!!$    endif
!!$
!!$    nzt = min((a%icp(lcl+1)-a%icp(icl)),&
!!$         & ((nza*(lrw+1-irw))/m) )
!!$    nz = 0 
!!$
!!$
!!$    call psb_ensure_size(nzin_+nzt,ia,info)
!!$    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)
!!$    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,val,info)
!!$
!!$    if (info /= psb_success_) return
!!$    isz = min(size(ia),size(ja),size(val))
!!$    if (present(iren)) then 
!!$      do i=icl, lcl
!!$        do j=a%icp(i), a%icp(i+1) - 1
!!$          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
!!$            nzin_ = nzin_ + 1
!!$            if (nzin_>isz) then 
!!$              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
!!$              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
!!$              call psb_ensure_size(int(1.25*nzin_)+1,val,info)
!!$              isz = min(size(ia),size(ja),size(val))
!!$            end if
!!$            nz    = nz + 1
!!$            val(nzin_) = a%val(j)
!!$            ia(nzin_)  = iren(a%ia(j))
!!$            ja(nzin_)  = iren(i)
!!$          end if
!!$        enddo
!!$      end do
!!$    else
!!$      do i=icl, lcl
!!$        do j=a%icp(i), a%icp(i+1) - 1
!!$          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
!!$            nzin_ = nzin_ + 1
!!$            if (nzin_>isz) then 
!!$              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
!!$              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
!!$              call psb_ensure_size(int(1.25*nzin_)+1,val,info)
!!$              isz = min(size(ia),size(ja),size(val))
!!$            end if
!!$            nz    = nz + 1
!!$            val(nzin_) = a%val(j)
!!$            ia(nzin_)  = (a%ia(j))
!!$            ja(nzin_)  = (i)
!!$          end if
!!$        enddo
!!$      end do
!!$    end if
!!$  end subroutine dsc_getrow

end subroutine psb_d_dsc_csgetrow

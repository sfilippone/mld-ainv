subroutine psb_d_dsc_cssm(alpha,a,x,beta,y,info,trans) 
  
  use psb_error_mod
  use psb_string_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_cssm
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:,:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_base_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999

!!$  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
!!$  m   = a%get_nrows()
!!$
!!$  if (size(x,1)<m) then 
!!$    info = 36
!!$    call psb_errpush(info,name,i_err=(/3,n,0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$  if (size(y,1)<m) then 
!!$    info = 36
!!$    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$  nc  = min(size(x,2) , size(y,2)) 
!!$
!!$  if (.not. (a%is_triangle())) then 
!!$    info = psb_err_invalid_mat_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  if (alpha == dzero) then
!!$    if (beta == dzero) then
!!$      do i = 1, m
!!$        y(i,:) = dzero
!!$      enddo
!!$    else
!!$      do  i = 1, m
!!$        y(i,:) = beta*y(i,:)
!!$      end do
!!$    endif
!!$    return
!!$  end if
!!$
!!$  if (beta == dzero) then 
!!$    call inner_dscsm(tra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
!!$         & a%icp,a%ia,a%val,x,size(x,1),y,size(y,1),info) 
!!$    do  i = 1, m
!!$      y(i,1:nc) = alpha*y(i,1:nc)
!!$    end do
!!$  else 
!!$    allocate(tmp(m,nc), stat=info) 
!!$    if(info /= psb_success_) then
!!$      info=psb_err_from_subroutine_
!!$      call psb_errpush(info,name,a_err='allocate')
!!$      goto 9999
!!$    end if
!!$
!!$    tmp(1:m,:) = x(1:m,1:nc)
!!$    call inner_dscsm(tra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
!!$         & a%icp,a%ia,a%val,tmp,size(tmp,1),y,size(y,1),info) 
!!$    do  i = 1, m
!!$      y(i,1:nc) = alpha*tmp(i,1:nc) + beta*y(i,1:nc)
!!$    end do
!!$  end if
!!$
!!$  if(info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    call psb_errpush(info,name,a_err='inner_dscsm')
!!$    goto 9999
!!$  end if
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
!!$  subroutine inner_dscsm(tra,lower,unit,nr,nc,&
!!$       & icp,ia,val,x,ldx,y,ldy,info) 
!!$    implicit none 
!!$    logical, intent(in)                 :: tra,lower,unit
!!$    integer, intent(in)                 :: nr,nc,ldx,ldy,icp(*),ia(*)
!!$    real(psb_dpk_), intent(in)          :: val(*), x(ldx,*)
!!$    real(psb_dpk_), intent(out)         :: y(ldy,*)
!!$    integer, intent(out)                :: info
!!$    integer :: i,j,k,m, ir, jc
!!$    real(psb_dpk_), allocatable  :: acc(:)
!!$
!!$    info = psb_success_
!!$    allocate(acc(nc), stat=info)
!!$    if(info /= psb_success_) then
!!$      info=psb_err_from_subroutine_
!!$      return
!!$    end if
!!$
!!$
!!$    if (tra) then 
!!$
!!$      if (lower) then 
!!$        if (unit) then 
!!$          do i=nr, 1, -1 
!!$            acc = dzero 
!!$            do j=a%icp(i), a%icp(i+1)-1
!!$              acc = acc + a%val(j)*y(a%ia(j),1:nc)
!!$            end do
!!$            y(i,1:nc) = x(i,1:nc) - acc
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=nr, 1, -1 
!!$            acc = dzero 
!!$            do j=a%icp(i)+1, a%icp(i+1)-1
!!$              acc = acc + a%val(j)*y(a%ia(j),1:nc)
!!$            end do
!!$            y(i,1:nc) = (x(i,1:nc) - acc)/a%val(a%icp(i))
!!$          end do
!!$        end if
!!$
!!$      else if (.not.lower) then 
!!$
!!$        if (unit) then 
!!$          do i=1, nr
!!$            acc = dzero 
!!$            do j=a%icp(i), a%icp(i+1)-1
!!$              acc = acc + a%val(j)*y(a%ia(j),1:nc)
!!$            end do
!!$            y(i,1:nc) = x(i,1:nc) - acc
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=1, nr
!!$            acc = dzero 
!!$            do j=a%icp(i), a%icp(i+1)-2
!!$              acc = acc + a%val(j)*y(a%ia(j),1:nc)
!!$            end do
!!$            y(i,1:nc) = (x(i,1:nc) - acc)/a%val(a%icp(i+1)-1)
!!$          end do
!!$        end if
!!$
!!$      end if
!!$
!!$    else if (.not.tra) then 
!!$
!!$      do i=1, nr
!!$        y(i,1:nc) = x(i,1:nc)
!!$      end do
!!$
!!$      if (lower) then
!!$
!!$        if (unit) then  
!!$          do i=1, nr
!!$            acc  = y(i,1:nc) 
!!$            do j=a%icp(i), a%icp(i+1)-1
!!$              jc    = a%ia(j)
!!$              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
!!$            end do
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=1, nr
!!$            y(i,1:nc) = y(i,1:nc)/a%val(a%icp(i))
!!$            acc    = y(i,1:nc) 
!!$            do j=a%icp(i)+1, a%icp(i+1)-1
!!$              jc      = a%ia(j)
!!$              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
!!$            end do
!!$          end do
!!$        end if
!!$
!!$      else if (.not.lower) then 
!!$
!!$        if (unit) then 
!!$          do i=nr, 1, -1
!!$            acc = y(i,1:nc) 
!!$            do j=a%icp(i), a%icp(i+1)-1
!!$              jc    = a%ia(j)
!!$              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
!!$            end do
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=nr, 1, -1
!!$            y(i,1:nc) = y(i,1:nc)/a%val(a%icp(i+1)-1)
!!$            acc  = y(i,1:nc) 
!!$            do j=a%icp(i), a%icp(i+1)-2
!!$              jc    = a%ia(j)
!!$              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
!!$            end do
!!$          end do
!!$        end if
!!$
!!$      end if
!!$    end if
!!$  end subroutine inner_dscsm

end subroutine psb_d_dsc_cssm

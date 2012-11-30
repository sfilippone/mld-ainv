subroutine psb_d_dsc_cssv(alpha,a,x,beta,y,info,trans) 
  
  use psb_error_mod
  use psb_string_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_cssv
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m, nnz, ir, jc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_dsc_cssv'
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

!!$
!!$  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
!!$  m = a%get_nrows()
!!$
!!$  if (.not. (a%is_triangle())) then 
!!$    info = psb_err_invalid_mat_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  end if
!!$
!!$  if (size(x,1)<m) then 
!!$    info = 36
!!$    call psb_errpush(info,name,i_err=(/3,m,0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$  if (size(y,1)<m) then 
!!$    info = 36
!!$    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  
!!$  if (alpha == dzero) then
!!$    if (beta == dzero) then
!!$      do i = 1, m
!!$        y(i) = dzero
!!$      enddo
!!$    else
!!$      do  i = 1, m
!!$        y(i) = beta*y(i)
!!$      end do
!!$    endif
!!$    return
!!$  end if
!!$
!!$  if (beta == dzero) then 
!!$    call inner_dscsv(tra,a%is_lower(),a%is_unit(),a%get_nrows(),&
!!$         & a%icp,a%ia,a%val,x,y) 
!!$    if (alpha == done) then 
!!$      ! do nothing
!!$    else if (alpha == -done) then 
!!$      do  i = 1, m
!!$        y(i) = -y(i)
!!$      end do
!!$    else
!!$      do  i = 1, m
!!$        y(i) = alpha*y(i)
!!$      end do
!!$    end if
!!$  else 
!!$    allocate(tmp(m), stat=info) 
!!$    if (info /= psb_success_) then 
!!$      return
!!$    end if
!!$    tmp(1:m) = x(1:m)
!!$    call inner_dscsv(tra,a%is_lower(),a%is_unit(),a%get_nrows(),&
!!$         & a%icp,a%ia,a%val,tmp,y) 
!!$    do  i = 1, m
!!$      y(i) = alpha*tmp(i) + beta*y(i)
!!$    end do
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

!!$contains 
!!$
!!$  subroutine inner_dscsv(tra,lower,unit,n,icp,ia,val,x,y) 
!!$    implicit none 
!!$    logical, intent(in)                 :: tra,lower,unit  
!!$    integer, intent(in)                 :: icp(*), ia(*),n
!!$    real(psb_dpk_), intent(in)          :: val(*)
!!$    real(psb_dpk_), intent(in)          :: x(*)
!!$    real(psb_dpk_), intent(out)         :: y(*)
!!$
!!$    integer :: i,j,k,m, ir, jc
!!$    real(psb_dpk_) :: acc
!!$
!!$    if (tra) then 
!!$
!!$      if (lower) then 
!!$        if (unit) then 
!!$          do i=n, 1, -1 
!!$            acc = dzero 
!!$            do j=icp(i), icp(i+1)-1
!!$              acc = acc + val(j)*y(ia(j))
!!$            end do
!!$            y(i) = x(i) - acc
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=n, 1, -1 
!!$            acc = dzero 
!!$            do j=icp(i)+1, icp(i+1)-1
!!$              acc = acc + val(j)*y(ia(j))
!!$            end do
!!$            y(i) = (x(i) - acc)/val(icp(i))
!!$          end do
!!$        end if
!!$
!!$      else if (.not.lower) then 
!!$
!!$        if (unit) then 
!!$          do i=1, n
!!$            acc = dzero 
!!$            do j=icp(i), icp(i+1)-1
!!$              acc = acc + val(j)*y(ia(j))
!!$            end do
!!$            y(i) = x(i) - acc
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=1, n
!!$            acc = dzero 
!!$            do j=icp(i), icp(i+1)-2
!!$              acc = acc + val(j)*y(ia(j))
!!$            end do
!!$            y(i) = (x(i) - acc)/val(icp(i+1)-1)
!!$          end do
!!$        end if
!!$
!!$      end if
!!$
!!$    else if (.not.tra) then 
!!$
!!$      do i=1, n
!!$        y(i) = x(i)
!!$      end do
!!$
!!$      if (lower) then 
!!$
!!$        if (unit) then 
!!$          do i=1, n
!!$            acc  = y(i) 
!!$            do j=icp(i), icp(i+1)-1
!!$              jc    = ia(j)
!!$              y(jc) = y(jc) - val(j)*acc 
!!$            end do
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=1, n
!!$            y(i) = y(i)/val(icp(i))
!!$            acc  = y(i) 
!!$            do j=icp(i)+1, icp(i+1)-1
!!$              jc    = ia(j)
!!$              y(jc) = y(jc) - val(j)*acc 
!!$            end do
!!$          end do
!!$        end if
!!$
!!$      else if (.not.lower) then 
!!$
!!$        if (unit) then 
!!$          do i=n, 1, -1
!!$            acc = y(i) 
!!$            do j=icp(i), icp(i+1)-1
!!$              jc    = ia(j)
!!$              y(jc) = y(jc) - val(j)*acc 
!!$            end do
!!$          end do
!!$        else if (.not.unit) then 
!!$          do i=n, 1, -1
!!$            y(i) = y(i)/val(icp(i+1)-1)
!!$            acc  = y(i) 
!!$            do j=icp(i), icp(i+1)-2
!!$              jc    = ia(j)
!!$              y(jc) = y(jc) - val(j)*acc 
!!$            end do
!!$          end do
!!$        end if
!!$
!!$      end if
!!$    end if
!!$  end subroutine inner_dscsv
!!$
end subroutine psb_d_dsc_cssv

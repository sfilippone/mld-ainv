
! == ===================================
!
!
!
! Computational routines
!
!
!
!
!
!
! == ===================================

subroutine psb_d_dsc_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csmv
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_dsc_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

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

  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
!!$
!!$  if (tra) then 
!!$    m = a%get_ncols()
!!$    n = a%get_nrows()
!!$  else
!!$    n = a%get_ncols()
!!$    m = a%get_nrows()
!!$  end if
!!$
!!$
!!$  if (size(x,1)<n) then 
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
!!$  if (tra) then 
!!$
!!$    if (beta == dzero) then 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = -acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$
!!$    else if (beta == done) then 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = y(i) + acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = y(i) -acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = y(i) + alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$    else if (beta == -done) then 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = -y(i) + acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = -y(i) -acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = -y(i) + alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$    else 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = beta*y(i) + acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = beta*y(i) - acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j))          
!!$          enddo
!!$          y(i) = beta*y(i) + alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$    end if
!!$
!!$  else if (.not.tra) then 
!!$
!!$    if (beta == dzero) then 
!!$      do i=1, m
!!$        y(i) = dzero
!!$      end do
!!$    else if (beta == done) then 
!!$      ! Do nothing
!!$    else if (beta == -done) then 
!!$      do i=1, m
!!$        y(i) = -y(i) 
!!$      end do
!!$    else
!!$      do i=1, m
!!$        y(i) = beta*y(i) 
!!$      end do
!!$    end if
!!$
!!$    if (alpha == done) then
!!$
!!$      do i=1,n
!!$        do j=a%icp(i), a%icp(i+1)-1
!!$          ir = a%ia(j)
!!$          y(ir) = y(ir) +  a%val(j)*x(i)
!!$        end do
!!$      enddo
!!$
!!$    else if (alpha == -done) then
!!$
!!$      do i=1,n
!!$        do j=a%icp(i), a%icp(i+1)-1
!!$          ir = a%ia(j)
!!$          y(ir) = y(ir) -  a%val(j)*x(i)
!!$        end do
!!$      enddo
!!$
!!$    else                    
!!$
!!$      do i=1,n
!!$        do j=a%icp(i), a%icp(i+1)-1
!!$          ir = a%ia(j)
!!$          y(ir) = y(ir) + alpha*a%val(j)*x(i)
!!$        end do
!!$      enddo
!!$
!!$    end if
!!$
!!$  endif
!!$
!!$  if (a%is_triangle().and.a%is_unit()) then 
!!$    do i=1, min(m,n)
!!$      y(i) = y(i) + alpha*x(i)
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

end subroutine psb_d_dsc_csmv

subroutine psb_d_dsc_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csmm
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_dsc_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
!!$  if (.not.a%is_asb()) then 
!!$    info = psb_err_invalid_mat_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$
!!$  if (tra) then 
!!$    m = a%get_ncols()
!!$    n = a%get_nrows()
!!$  else
!!$    n = a%get_ncols()
!!$    m = a%get_nrows()
!!$  end if
!!$
!!$  if (size(x,1)<n) then 
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
!!$
!!$  nc = min(size(x,2) , size(y,2) )
!!$
!!$  allocate(acc(nc), stat=info)
!!$  if(info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    call psb_errpush(info,name,a_err='allocate')
!!$    goto 9999
!!$  end if
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
!!$  if (tra) then 
!!$
!!$    if (beta == dzero) then 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = -acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$
!!$    else if (beta == done) then 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = y(i,:) + acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = y(i,:) -acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = y(i,:) + alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$    else if (beta == -done) then 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = -y(i,:) + acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = -y(i,:) -acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = -y(i,:) + alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$    else 
!!$
!!$      if (alpha == done) then 
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = beta*y(i,:) + acc
!!$        end do
!!$
!!$      else if (alpha == -done) then 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = beta*y(i,:) - acc
!!$        end do
!!$
!!$      else 
!!$
!!$        do i=1,m 
!!$          acc  = dzero
!!$          do j=a%icp(i), a%icp(i+1)-1
!!$            acc  = acc + a%val(j) * x(a%ia(j),:)          
!!$          enddo
!!$          y(i,:) = beta*y(i,:) + alpha*acc
!!$        end do
!!$
!!$      end if
!!$
!!$    end if
!!$
!!$  else if (.not.tra) then 
!!$
!!$    if (beta == dzero) then 
!!$      do i=1, m
!!$        y(i,:) = dzero
!!$      end do
!!$    else if (beta == done) then 
!!$      ! Do nothing
!!$    else if (beta == -done) then 
!!$      do i=1, m
!!$        y(i,:) = -y(i,:) 
!!$      end do
!!$    else
!!$      do i=1, m
!!$        y(i,:) = beta*y(i,:) 
!!$      end do
!!$    end if
!!$
!!$    if (alpha == done) then
!!$
!!$      do i=1,n
!!$        do j=a%icp(i), a%icp(i+1)-1
!!$          ir = a%ia(j)
!!$          y(ir,:) = y(ir,:) +  a%val(j)*x(i,:)
!!$        end do
!!$      enddo
!!$
!!$    else if (alpha == -done) then
!!$
!!$      do i=1,n
!!$        do j=a%icp(i), a%icp(i+1)-1
!!$          ir = a%ia(j)
!!$          y(ir,:) = y(ir,:) -  a%val(j)*x(i,:)
!!$        end do
!!$      enddo
!!$
!!$    else                    
!!$
!!$      do i=1,n
!!$        do j=a%icp(i), a%icp(i+1)-1
!!$          ir = a%ia(j)
!!$          y(ir,:) = y(ir,:) + alpha*a%val(j)*x(i,:)
!!$        end do
!!$      enddo
!!$
!!$    end if
!!$
!!$  endif
!!$
!!$  if (a%is_triangle().and.a%is_unit()) then 
!!$    do i=1, min(m,n)
!!$      y(i,:) = y(i,:) + alpha*x(i,:)
!!$    end do
!!$  end if

  call psb_erractionrestore(err_act)
  return
9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_csmm


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

function psb_d_dsc_csnmi(a) result(res)
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csnmi
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nr, ir, jc, nc, info
  real(psb_dpk_), allocatable  :: acc(:) 
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_csnmi'
  logical, parameter :: debug=.false.


  res = dzero 
!!$  nr = a%get_nrows()
!!$  nc = a%get_ncols()
!!$  allocate(acc(nr),stat=info)
!!$  if (info /= psb_success_) then 
!!$    return
!!$  end if
!!$  acc(:) = dzero
!!$  do i=1, nc
!!$    do j=a%icp(i),a%icp(i+1)-1  
!!$      acc(a%ia(j)) = acc(a%ia(j)) + abs(a%val(j))
!!$    end do
!!$  end do
!!$  do i=1, nr
!!$    res = max(res,acc(i))
!!$  end do
!!$  deallocate(acc)

end function psb_d_dsc_csnmi


function psb_d_dsc_csnm1(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csnm1

  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_dsc_csnm1'
  logical, parameter :: debug=.false.


  res = dzero 
!!$  m = a%get_nrows()
!!$  n = a%get_ncols()
!!$  do j=1, n
!!$    acc = dzero 
!!$    do k=a%icp(j),a%icp(j+1)-1
!!$      acc = acc + abs(a%val(k))
!!$    end do
!!$    res = max(res,acc)
!!$  end do
!!$  
!!$  return
!!$
end function psb_d_dsc_csnm1

subroutine psb_d_dsc_colsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_colsum
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999

!!$  m = a%get_ncols()
!!$  if (size(d) < m) then 
!!$    info=psb_err_input_asize_small_i_
!!$    int_err(1) = 1
!!$    int_err(2) = size(d)
!!$    int_err(3) = m
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  do i = 1, a%get_ncols()
!!$    d(i) = dzero
!!$    do j=a%icp(i),a%icp(i+1)-1  
!!$      d(i) = d(i) + (a%val(j))
!!$    end do
!!$  end do
!!$
  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_colsum

subroutine psb_d_dsc_aclsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_aclsum
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  m = a%get_ncols()
!!$  if (size(d) < m) then 
!!$    info=psb_err_input_asize_small_i_
!!$    int_err(1) = 1
!!$    int_err(2) = size(d)
!!$    int_err(3) = m
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  do i = 1, a%get_ncols()
!!$    d(i) = dzero
!!$    do j=a%icp(i),a%icp(i+1)-1  
!!$      d(i) = d(i) + abs(a%val(j))
!!$    end do
!!$  end do
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

end subroutine psb_d_dsc_aclsum

subroutine psb_d_dsc_rowsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_rowsum
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  m = a%get_ncols()
!!$  n = a%get_nrows()
!!$  if (size(d) < n) then 
!!$    info=psb_err_input_asize_small_i_
!!$    int_err(1) = 1
!!$    int_err(2) = size(d)
!!$    int_err(3) = n
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  d   = dzero
!!$
!!$  do i=1, m
!!$    do j=a%icp(i),a%icp(i+1)-1
!!$      k = a%ia(j)
!!$      d(k) = d(k) + (a%val(k))
!!$    end do
!!$  end do
!!$
  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_rowsum

subroutine psb_d_dsc_arwsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_arwsum
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='arwsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$  m = a%get_ncols()
!!$  n = a%get_nrows()
!!$  if (size(d) < n) then 
!!$    info=psb_err_input_asize_small_i_
!!$    int_err(1) = 1
!!$    int_err(2) = size(d)
!!$    int_err(3) = n
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  d   = dzero
!!$
!!$  do i=1, m
!!$    do j=a%icp(i),a%icp(i+1)-1
!!$      k = a%ia(j)
!!$      d(k) = d(k) + abs(a%val(k))
!!$    end do
!!$  end do
!!$
  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_arwsum


subroutine psb_d_dsc_get_diag(a,d,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_get_diag
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act, mnm, i, j, k
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999

  mnm = min(a%get_nrows(),a%get_ncols())
!!$  if (size(d) < mnm) then 
!!$    info=psb_err_input_asize_invalid_i_
!!$    call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$
!!$  if (a%is_triangle().and.a%is_unit()) then 
!!$    d(1:mnm) = done 
!!$  else
!!$    do i=1, mnm
!!$      d(i) = dzero
!!$      do k=a%icp(i),a%icp(i+1)-1
!!$        j=a%ia(k)
!!$        if ((j == i) .and.(j <= mnm )) then 
!!$          d(i) = a%val(k)
!!$        endif
!!$      enddo
!!$    end do
!!$  endif
!!$  do i=mnm+1,size(d) 
!!$    d(i) = dzero
!!$  end do
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_get_diag


subroutine psb_d_dsc_scal(d,a,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_scal
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act,mnm, i, j, n
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  n = a%get_ncols()
!!$  if (size(d) < n) then 
!!$    info=psb_err_input_asize_invalid_i_
!!$    call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
!!$    goto 9999
!!$  end if
!!$
!!$  do i=1, n
!!$    do j = a%icp(i), a%icp(i+1) -1 
!!$      a%val(j) = a%val(j) * d(a%ia(j))
!!$    end do
!!$  enddo
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

end subroutine psb_d_dsc_scal


subroutine psb_d_dsc_scals(d,a,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_scals
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer, intent(out)            :: info

  Integer :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  do i=1,a%get_nzeros()
!!$    a%val(i) = a%val(i) * d
!!$  enddo
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

end subroutine psb_d_dsc_scals


! == =================================== 
!
!
!
! Data management
!
!
!
!
!
! == ===================================   

subroutine psb_d_dsc_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,dscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csgetptn
  implicit none

  class(psb_d_dsc_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
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
!!$
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
!!$  call dsc_getptn(imin,imax,jmin_,jmax_,a,nz,ia,ja,nzin_,append_,info,iren)
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

!!$contains
!!$
!!$  subroutine dsc_getptn(imin,imax,jmin,jmax,a,nz,ia,ja,nzin,append,info,&
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
!!$    integer, intent(in)                  :: nzin
!!$    logical, intent(in)                  :: append
!!$    integer                              :: info
!!$    integer, optional                    :: iren(:)
!!$    integer  :: nzin_, nza, idx,i,j,k, nzt, irw, lrw, icl, lcl,m,isz
!!$    integer  :: debug_level, debug_unit
!!$    character(len=20) :: name='coo_getrow'
!!$
!!$    debug_unit  = psb_get_debug_unit()
!!$    debug_level = psb_get_debug_level()
!!$
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
!!$
!!$    nzt = min((a%icp(lcl+1)-a%icp(icl)),&
!!$         & ((nza*(lcl+1-icl))/a%get_ncols()) )
!!$    nz = 0 
!!$
!!$
!!$    call psb_ensure_size(nzin_+nzt,ia,info)
!!$    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)
!!$
!!$    if (info /= psb_success_) return
!!$    isz = min(size(ia),size(ja))
!!$    if (present(iren)) then 
!!$      do i=icl, lcl
!!$        do j=a%icp(i), a%icp(i+1) - 1
!!$          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
!!$            nzin_ = nzin_ + 1
!!$            if (nzin_>isz) then 
!!$              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
!!$              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
!!$              isz = min(size(ia),size(ja))
!!$            end if
!!$            nz    = nz + 1
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
!!$              isz = min(size(ia),size(ja))
!!$            end if
!!$            nz    = nz + 1
!!$            ia(nzin_)  = (a%ia(j))
!!$            ja(nzin_)  = (i)
!!$          end if
!!$        enddo
!!$      end do
!!$    end if
!!$    
!!$  end subroutine dsc_getptn
  
end subroutine psb_d_dsc_csgetptn


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



subroutine psb_d_dsc_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csput
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='d_dsc_csput'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)

  call psb_erractionsave(err_act)
  info = psb_success_


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$  if (nz <= 0) then 
!!$    info = psb_err_iarg_neg_
!!$    int_err(1)=1
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$  if (size(ia) < nz) then 
!!$    info = psb_err_input_asize_invalid_i_
!!$    int_err(1)=2
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  if (size(ja) < nz) then 
!!$    info = psb_err_input_asize_invalid_i_
!!$    int_err(1)=3
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$  if (size(val) < nz) then 
!!$    info = psb_err_input_asize_invalid_i_
!!$    int_err(1)=4
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  if (nz == 0) return
!!$
!!$  nza  = a%get_nzeros()
!!$
!!$  if (a%is_bld()) then 
!!$    ! Build phase should only ever be in COO
!!$    info = psb_err_invalid_mat_state_
!!$
!!$  else  if (a%is_upd()) then 
!!$    call  psb_d_dsc_srch_upd(nz,ia,ja,val,a,&
!!$         & imin,imax,jmin,jmax,info,gtl)
!!$
!!$    if (info /= psb_success_) then  
!!$
!!$      info = psb_err_invalid_mat_state_
!!$    end if
!!$
!!$  else 
!!$    ! State is wrong.
!!$    info = psb_err_invalid_mat_state_
!!$  end if
!!$  if (info /= psb_success_) then
!!$    call psb_errpush(info,name)
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
!!$
!!$contains
!!$
!!$  subroutine psb_d_dsc_srch_upd(nz,ia,ja,val,a,&
!!$       & imin,imax,jmin,jmax,info,gtl)
!!$
!!$    use psb_const_mod
!!$    use psb_realloc_mod
!!$    use psb_string_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$
!!$    class(psb_d_dsc_sparse_mat), intent(inout) :: a
!!$    integer, intent(in) :: nz, imin,imax,jmin,jmax
!!$    integer, intent(in) :: ia(:),ja(:)
!!$    real(psb_dpk_), intent(in) :: val(:)
!!$    integer, intent(out) :: info
!!$    integer, intent(in), optional  :: gtl(:)
!!$    integer  :: i,ir,ic, ilr, ilc, ip, &
!!$         & i1,i2,nr,nc,nnz,dupl,ng, nar, nac
!!$    integer              :: debug_level, debug_unit
!!$    character(len=20)    :: name='d_dsc_srch_upd'
!!$
!!$    info = psb_success_
!!$    debug_unit  = psb_get_debug_unit()
!!$    debug_level = psb_get_debug_level()
!!$
!!$    dupl = a%get_dupl()
!!$
!!$    if (.not.a%is_sorted()) then 
!!$      info = -4
!!$      return
!!$    end if
!!$
!!$    ilr = -1 
!!$    ilc = -1 
!!$    nnz = a%get_nzeros()
!!$    nar = a%get_nrows()
!!$    nac = a%get_ncols()
!!$
!!$    if (present(gtl)) then 
!!$      ng = size(gtl)
!!$
!!$      select case(dupl)
!!$      case(psb_dupl_ovwrt_,psb_dupl_err_)
!!$        ! Overwrite.
!!$        ! Cannot test for error, should have been caught earlier.
!!$
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
!!$            ir = gtl(ir)
!!$            ic = gtl(ic)
!!$            if ((ic > 0).and.(ic <= nac)) then 
!!$              i1 = a%icp(ic)
!!$              i2 = a%icp(ic+1)
!!$              nr=i2-i1
!!$
!!$              ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$              if (ip>0) then 
!!$                a%val(i1+ip-1) = val(i)
!!$              else
!!$                if (debug_level >= psb_debug_serial_) &
!!$                     & write(debug_unit,*) trim(name),&
!!$                     & ': Was searching ',ir,' in: ',i1,i2,&
!!$                     & ' : ',a%ia(i1:i2-1)
!!$                info = i
!!$                return
!!$              end if
!!$
!!$            else
!!$
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Discarding row that does not belong to us.'
!!$            end if
!!$          end if
!!$        end do
!!$
!!$      case(psb_dupl_add_)
!!$        ! Add
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
!!$            ir = gtl(ir)
!!$            ic = gtl(ic)
!!$            if ((ic > 0).and.(ic <= nac)) then 
!!$              i1 = a%icp(ic)
!!$              i2 = a%icp(ic+1)
!!$              nr=i2-i1
!!$
!!$              ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$              if (ip>0) then 
!!$                a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
!!$              else
!!$                if (debug_level >= psb_debug_serial_) &
!!$                     & write(debug_unit,*) trim(name),&
!!$                     & ': Was searching ',ir,' in: ',i1,i2,&
!!$                     & ' : ',a%ia(i1:i2-1)
!!$                info = i
!!$                return
!!$              end if
!!$            else
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Discarding row that does not belong to us.'
!!$            end if
!!$
!!$          end if
!!$        end do
!!$
!!$      case default
!!$        info = -3
!!$        if (debug_level >= psb_debug_serial_) &
!!$             & write(debug_unit,*) trim(name),&
!!$             & ': Duplicate handling: ',dupl
!!$      end select
!!$
!!$    else
!!$
!!$      select case(dupl)
!!$      case(psb_dupl_ovwrt_,psb_dupl_err_)
!!$        ! Overwrite.
!!$        ! Cannot test for error, should have been caught earlier.
!!$
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$
!!$          if ((ic > 0).and.(ic <= nac)) then 
!!$            i1 = a%icp(ic)
!!$            i2 = a%icp(ic+1)
!!$            nr=i2-i1
!!$
!!$            ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$            if (ip>0) then 
!!$              a%val(i1+ip-1) = val(i)
!!$            else
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Was searching ',ir,' in: ',i1,i2,&
!!$                   & ' : ',a%ia(i1:i2-1)
!!$              info = i
!!$              return
!!$            end if
!!$
!!$          else
!!$            if (debug_level >= psb_debug_serial_) &
!!$                 & write(debug_unit,*) trim(name),&
!!$                 & ': Discarding col that does not belong to us.'
!!$          end if
!!$
!!$        end do
!!$
!!$      case(psb_dupl_add_)
!!$        ! Add
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$          if ((ic > 0).and.(ic <= nac)) then 
!!$            i1 = a%icp(ic)
!!$            i2 = a%icp(ic+1)
!!$            nr=i2-i1
!!$
!!$            ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$            if (ip>0) then 
!!$              a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
!!$            else
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Was searching ',ir,' in: ',i1,i2,&
!!$                   & ' : ',a%ia(i1:i2-1)
!!$              info = i
!!$              return
!!$            end if
!!$          else
!!$            if (debug_level >= psb_debug_serial_) &
!!$                 & write(debug_unit,*) trim(name),&
!!$                 & ': Discarding col that does not belong to us.'
!!$          end if
!!$        end do
!!$
!!$      case default
!!$        info = -3
!!$        if (debug_level >= psb_debug_serial_) &
!!$             & write(debug_unit,*) trim(name),&
!!$             & ': Duplicate handling: ',dupl
!!$      end select
!!$
!!$    end if
!!$
!!$  end subroutine psb_d_dsc_srch_upd

end subroutine psb_d_dsc_csput



subroutine psb_d_cp_dsc_from_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_from_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  type(psb_d_coo_sparse_mat)   :: tmp
  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  ! This is to have fix_coo called behind the scenes
  call tmp%cp_from_coo(b,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

end subroutine psb_d_cp_dsc_from_coo



subroutine psb_d_cp_dsc_to_coo(a,b,info) 
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_to_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(in)  :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                      :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,k,irw, err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)

  if (allocated(a%cols)) then 
    k = 0
    do i=1, nc
      do j=1,a%cols(i)%nz
        k = k + 1
        b%ia(k)  = a%cols(i)%idx(j)
        b%ja(k)  = i
        b%val(k) = a%cols(i)%val(j)
      end do
    end do
  end if
  call b%set_nzeros(nza)
  call b%fix(info)

end subroutine psb_d_cp_dsc_to_coo


subroutine psb_d_mv_dsc_to_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_mv_dsc_to_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout)   :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  
  ! Too tired to figure out a smart mv right now. 
  call a%cp_to_coo(b,info)
  call a%free()

end subroutine psb_d_mv_dsc_to_coo



subroutine psb_d_mv_dsc_from_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_mv_dsc_from_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc, icl
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ! Too tired to figure out a smart mv right now. 
  call a%cp_from_coo(b,info)
  call b%free()
!!$
!!$  call b%fix(info, idir=1)
!!$  if (info /= psb_success_) return
!!$
!!$  nr  = b%get_nrows()
!!$  nc  = b%get_ncols()
!!$  nza = b%get_nzeros()
!!$  
!!$  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
!!$
!!$  ! Dirty trick: call move_alloc to have the new data allocated just once.
!!$  call move_alloc(b%ja,itemp)
!!$  call move_alloc(b%ia,a%ia)
!!$  call move_alloc(b%val,a%val)
!!$  call psb_realloc(max(nr+1,nc+1),a%icp,info)
!!$  call b%free()
!!$
!!$  if (nza <= 0) then 
!!$    a%icp(:) = 1
!!$  else
!!$    a%icp(1) = 1
!!$    if (nc < itemp(nza)) then 
!!$      write(debug_unit,*) trim(name),': CLSHR=.false. : ',&
!!$           &nc,itemp(nza),' Expect trouble!'
!!$      info = 12
!!$    end if
!!$
!!$    j = 1 
!!$    i = 1
!!$    icl = itemp(j) 
!!$
!!$    outer: do 
!!$      inner: do 
!!$        if (i >= icl) exit inner
!!$        if (i > nc) then 
!!$          write(debug_unit,*) trim(name),&
!!$               & 'Strange situation: i>nr ',i,nc,j,nza,icl
!!$          exit outer
!!$        end if
!!$        a%icp(i+1) = a%icp(i) 
!!$        i = i + 1
!!$      end do inner
!!$      j = j + 1
!!$      if (j > nza) exit
!!$      if (itemp(j) /= icl) then 
!!$        a%icp(i+1) = j
!!$        icl = itemp(j) 
!!$        i = i + 1
!!$      endif
!!$      if (i > nc) exit
!!$    enddo outer
!!$    !
!!$    ! Cleanup empty rows at the end
!!$    !
!!$    if (j /= (nza+1)) then 
!!$      write(debug_unit,*) trim(name),': Problem from loop :',j,nza
!!$      info = 13
!!$    endif
!!$    do 
!!$      if (i > nc) exit
!!$      a%icp(i+1) = j
!!$      i = i + 1
!!$    end do
!!$
!!$  endif
!!$

end subroutine psb_d_mv_dsc_from_coo


subroutine psb_d_mv_dsc_to_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_mv_dsc_to_fmt
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout)  :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%mv_to_coo(b,info)
    ! Need to fix trivial copies! 
  type is (psb_d_dsc_sparse_mat) 
    call b%psb_d_base_sparse_mat%mv_from(a%psb_d_base_sparse_mat)
    call move_alloc(a%cols, b%cols)
    call a%free()
    
  class default
    call a%mv_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_mv_dsc_to_fmt
!!$

subroutine psb_d_cp_dsc_to_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_to_fmt
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(in)   :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                       :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_


  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_to_coo(b,info)

  type is (psb_d_dsc_sparse_mat) 
    call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)
    call psb_safe_ab_cpy( a%cols, b%cols , info)

  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_cp_dsc_to_fmt


subroutine psb_d_mv_dsc_from_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_mv_dsc_from_fmt
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout)  :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                         :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

  type is (psb_d_dsc_sparse_mat) 
    call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
    call move_alloc(b%cols, a%cols)
    call b%free()

  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_mv_dsc_from_fmt



subroutine psb_d_cp_dsc_from_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_cp_dsc_from_fmt
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_d_dsc_sparse_mat) 
    call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
    call psb_safe_ab_cpy( b%cols, a%cols , info)

  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_cp_dsc_from_fmt

subroutine psb_d_dsc_mold(a,b,info) 
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_mold
  use psb_error_mod
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer, intent(out)                    :: info
  Integer :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  allocate(psb_d_dsc_sparse_mat :: b, stat=info)

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_dsc_mold

!!$
subroutine  psb_d_dsc_reallocate_nz(nz,a) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_reallocate_nz
  implicit none 
  integer, intent(in) :: nz
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='d_dsc_reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$
!!$  call psb_realloc(nz,a%ia,info)
!!$  if (info == psb_success_) call psb_realloc(nz,a%val,info)
!!$  if (info == psb_success_) call psb_realloc(max(nz,a%get_nrows()+1,a%get_ncols()+1),a%icp,info)
!!$  if (info /= psb_success_) then 
!!$    call psb_errpush(psb_err_alloc_dealloc_,name)
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

end subroutine psb_d_dsc_reallocate_nz



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

subroutine psb_d_dsc_reinit(a,clear)
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_reinit
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  Integer :: err_act, info,i
  character(len=20)  :: name='reinit'
  logical  :: clear_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_


  if (present(clear)) then 
    clear_ = clear
  else
    clear_ = .true.
  end if

  if (a%is_bld() .or. a%is_upd()) then 
    ! do nothing
    return
  else if (a%is_asb()) then 
    if (clear_) then 
      if (allocated(a%cols)) then 
        do i=1,size(a%cols)
          a%cols(i)%val(:) = dzero
        end do
      end if
    end if
    call a%set_upd()
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
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

end subroutine psb_d_dsc_reinit

subroutine  psb_d_dsc_trim(a)
  use psb_realloc_mod
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_trim
  implicit none 
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  Integer :: err_act, info, nz, n
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$  n   = a%get_ncols()
!!$  nz  = a%get_nzeros()
!!$  if (info == psb_success_) call psb_realloc(n+1,a%icp,info)
!!$  if (info == psb_success_) call psb_realloc(nz,a%ia,info)
!!$  if (info == psb_success_) call psb_realloc(nz,a%val,info)
!!$
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

end subroutine psb_d_dsc_trim

subroutine  psb_d_dsc_allocate_mnnz(m,n,a,nz) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_allocate_mnnz
  implicit none 
  integer, intent(in) :: m,n
  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  integer, intent(in), optional :: nz
  Integer :: err_act, info, nz_,i
  character(len=20)  :: name='allocate_mnz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (m < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/1,0,0,0,0/))
    goto 9999
  endif
  if (n < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/2,0,0,0,0/))
    goto 9999
  endif
  if (present(nz)) then 
    nz_ = nz
  else
    nz_ = max(7*m,7*n,1)
  end if
  if (nz_ < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/3,0,0,0,0/))
    goto 9999
  endif

  if (info == psb_success_) call psb_realloc(n,a%cols,info)
  if (info == psb_success_) then 
    nz_ = (nz_ + n - 1 )/n
    do i=1,n
      call psb_realloc(nz_,a%cols(i)%idx,info)
      if (info == psb_success_) call psb_realloc(nz_,a%cols(i)%val,info)
      if (info /= psb_success_) exit
    end do
  end if
      
  if (info == psb_success_) then 
    call a%set_nrows(m)
    call a%set_ncols(n)
    call a%set_bld()
    call a%set_triangle(.false.)
    call a%set_unit(.false.)
    call a%set_dupl(psb_dupl_def_)
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

end subroutine psb_d_dsc_allocate_mnnz

subroutine psb_d_dsc_print(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_string_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_print
  implicit none 

  integer, intent(in)               :: iout
  class(psb_d_dsc_sparse_mat), intent(in) :: a   
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  Integer :: err_act
  character(len=20)  :: name='d_dsc_print'
  logical, parameter :: debug=.false.

  character(len=80)                 :: frmtv 
  integer  :: irs,ics,i,j, nmx, ni, nr, nc, nz

  if (present(eirs)) then 
    irs = eirs
  else
    irs = 0
  endif
  if (present(eics)) then 
    ics = eics
  else
    ics = 0
  endif

  if (present(head)) then 
    write(iout,'(a)') '%%MatrixMarket matrix coordinate real general'
    write(iout,'(a,a)') '% ',head 
    write(iout,'(a)') '%'    
    write(iout,'(a,a)') '% COO'
  endif

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nz  = a%get_nzeros()
  nmx = max(nr,nc,1)
  ni  = floor(log10(1.0*nmx)) + 1

  write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),es26.18,1x,2(i',ni,',1x))'
  write(iout,*) nr, nc, nz 
!!$  if(present(iv)) then 
!!$    do i=1, nc
!!$      do j=a%icp(i),a%icp(i+1)-1 
!!$        write(iout,frmtv) iv(a%ia(j)),iv(i),a%val(j)
!!$      end do
!!$    enddo
!!$  else      
!!$    if (present(ivr).and..not.present(ivc)) then 
!!$      do i=1, nc
!!$        do j=a%icp(i),a%icp(i+1)-1 
!!$          write(iout,frmtv) ivr(a%ia(j)),i,a%val(j)
!!$        end do
!!$      enddo
!!$    else if (present(ivr).and.present(ivc)) then 
!!$      do i=1, nc
!!$        do j=a%icp(i),a%icp(i+1)-1 
!!$          write(iout,frmtv) ivr(a%ia(j)),ivc(i),a%val(j)
!!$        end do
!!$      enddo
!!$    else if (.not.present(ivr).and.present(ivc)) then 
!!$      do i=1, nc
!!$        do j=a%icp(i),a%icp(i+1)-1 
!!$          write(iout,frmtv) (a%ia(j)),ivc(i),a%val(j)
!!$        end do
!!$      enddo
!!$    else if (.not.present(ivr).and..not.present(ivc)) then 
!!$      do i=1, nc
!!$        do j=a%icp(i),a%icp(i+1)-1 
!!$          write(iout,frmtv) (a%ia(j)),(i),a%val(j)
!!$        end do
!!$      enddo
!!$    endif
!!$  endif

end subroutine psb_d_dsc_print

subroutine psb_d_dsc_cp_from(a,b)
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_cp_from
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  type(psb_d_dsc_sparse_mat), intent(in)   :: b


  Integer :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_success_

  call a%allocate(b%get_nrows(),b%get_ncols(),b%get_nzeros())
  call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
  call psb_safe_ab_cpy( b%cols, a%cols , info)

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_dsc_cp_from

subroutine psb_d_dsc_mv_from(a,b)
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_mv_from
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout)  :: a
  type(psb_d_dsc_sparse_mat), intent(inout) :: b


  Integer :: err_act, info
  character(len=20)  :: name='mv_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
  call move_alloc(b%cols, a%cols)
  call b%free()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_dsc_mv_from




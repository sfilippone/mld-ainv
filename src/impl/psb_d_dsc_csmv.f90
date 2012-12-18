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

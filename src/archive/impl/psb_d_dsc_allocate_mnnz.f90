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

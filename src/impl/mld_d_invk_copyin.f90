!  
!   
!                       MLD-AINV: Approximate Inverse plugin for
!                             MLD2P4  version 2.0
!    
!    (C) Copyright 2012
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine mld_d_invk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,ktrw,trw,info,sign,inlevs)
  

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_d_invk_copyin

  implicit none

  ! Arguments 
  type(psb_d_csr_sparse_mat), intent(in)    :: a
  type(psb_d_coo_sparse_mat), intent(inout) :: trw
  integer, intent(in)                  :: i,m,jmin,jmax
  integer, intent(inout)               :: ktrw,info
  integer, intent(inout)               :: rowlevs(:)
  real(psb_dpk_), intent(inout)        :: row(:)
  type(psb_i_heap), intent(inout)      :: heap
  real(psb_dpk_), optional, intent(in) :: sign
  integer, intent(in), optional        :: inlevs(:)

  ! Local variables
  integer             :: k,j,irb,err_act, nz
  integer, parameter  :: nrb=16
  real(psb_dpk_)      :: sign_
  character(len=20), parameter  :: name='invk_copyin'
  character(len=20)             :: ch_err

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  call heap%init(info) 
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_init_heap')
    goto 9999
  end if

  if (present(sign)) then 
    sign_ = sign
  else
    sign_ = done
  end if


  !
  ! Take a fast shortcut if the matrix is stored in CSR format
  !
  if (present(inlevs)) then 
    do j = a%irp(i), a%irp(i+1) - 1
      k          = a%ja(j)
      if ((jmin<=k).and.(k<=jmax)) then 
        row(k)     = sign_ * a%val(j)
        rowlevs(k) = inlevs(j)
        call heap%insert(k,info)
      end if
    end do
  else
    do j = a%irp(i), a%irp(i+1) - 1
      k          = a%ja(j)
      if ((jmin<=k).and.(k<=jmax)) then 
        row(k)     = sign_ * a%val(j)
        rowlevs(k) = 0
        call heap%insert(k,info)
      end if
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_d_invk_copyin

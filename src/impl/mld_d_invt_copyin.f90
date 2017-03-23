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
subroutine mld_d_invt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
     & irwt,ktrw,trw,info,sign)
  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_d_invt_copyin
  implicit none 
  type(psb_d_csr_sparse_mat), intent(in)    :: a
  type(psb_d_coo_sparse_mat), intent(inout) :: trw
  integer, intent(in)                  :: i, m,jmin,jmax,jd
  integer, intent(inout)               :: ktrw,nlw,nup,jmaxup,info
  integer, intent(inout)               :: irwt(:)
  real(psb_dpk_), intent(inout)        :: nrmi,row(:)
  type(psb_i_heap), intent(inout)      :: heap
  real(psb_dpk_), intent(in), optional :: sign

  integer                     :: k,j,irb,kin,nz, err_act
  integer, parameter          :: nrb=16
  real(psb_dpk_)              :: dmaxup, sign_
  real(psb_dpk_), external    :: dnrm2
  character(len=20), parameter  :: name='invt_copyin'

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)

  call heap%init(info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_init_heap')
    goto 9999
  end if
  sign_ = done
  if (present(sign)) sign_ = sign
  !
  ! nrmi is the norm of the current sparse row (for the time being,
  ! we use the 2-norm).
  ! NOTE: the 2-norm below includes also elements that are outside
  ! [jmin:jmax] strictly. Is this really important? TO BE CHECKED.
  !

  nlw    = 0
  nup    = 0
  jmaxup = 0
  dmaxup = dzero
  nrmi   = dzero

  do j = a%irp(i), a%irp(i+1) - 1
    k = a%ja(j)
    if ((jmin<=k).and.(k<=jmax)) then 
      row(k)     = sign_ * a%val(j)
      call heap%insert(k,info)
      irwt(k) = 1
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if
    end if
    if (k<jd) nlw = nlw + 1 
    if (k>jd) then 
      nup = nup + 1
      if (abs(row(k))>dmaxup) then 
        jmaxup = k
        dmaxup = abs(row(k))
      end if
    end if
  end do
  nz   = a%irp(i+1) - a%irp(i)
  nrmi = dnrm2(nz,a%val(a%irp(i):),ione)

  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine mld_d_invt_copyin

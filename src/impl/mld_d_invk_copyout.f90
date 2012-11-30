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
subroutine mld_d_invk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
     &  l2,uia1,uia2,uaspk,info)

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_d_invk_copyout

  implicit none 

  ! Arguments
  integer, intent(in)                        :: fill_in, i, m, nidx
  integer, intent(inout)                     :: l2, info
  integer, intent(inout)                     :: rowlevs(:), idxs(:)
  integer, allocatable, intent(inout)        :: uia1(:), uia2(:)
  real(psb_dpk_), allocatable, intent(inout) :: uaspk(:)
  real(psb_dpk_), intent(inout)              :: row(:)

  ! Local variables
  integer               :: j,isz,err_act,int_err(5),idxp
  character(len=20), parameter  :: name='mld_diluk_factint'
  character(len=20)             :: ch_err

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)


  do idxp=1,nidx

    j = idxs(idxp)


    if (j>=i) then 
      !
      ! Copy the upper part of the row
      ! 
      if (rowlevs(j) <= fill_in) then 
        l2     = l2 + 1 
        if (size(uaspk) < l2) then 
          ! 
          ! Figure out a good reallocation size!
          !
          isz  = max(int(1.2*l2),l2+100)
          call psb_realloc(isz,uaspk,info) 
          if (info == psb_success_) call psb_realloc(isz,uia1,info) 
          if (info /= psb_success_) then 
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='Allocate')
            goto 9999
          end if
        end if
        uia1(l2)   = j
        uaspk(l2)  = row(j)
      end if
      !
      ! Re-initialize row(j) and rowlevs(j)
      !
      row(j)     = dzero
      rowlevs(j) = -(m+1)
    end if
  end do

  uia2(i+1) = l2 + 1

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine mld_d_invk_copyout

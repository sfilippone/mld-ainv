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
subroutine mld_d_ainv_solver_csetc(sv,what,val,info,idx)
  

  use psb_base_mod
  use mld_d_ainv_solver, mld_protect_name => mld_d_ainv_solver_csetc

  Implicit None

  ! Arguments
  class(mld_d_ainv_solver_type), intent(inout) :: sv
  character(len=*), intent(in)                 :: what 
  character(len=*), intent(in)                 :: val
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), intent(in), optional      :: idx
  !
  integer(psb_ipk_ ) :: err_act, ival
  character(len=20)  :: name='mld_d_ainv_solver_setc'

  info = psb_success_
  call psb_erractionsave(err_act)

  select case(psb_toupper(trim(what)))
  case('AINV_ALG')
    sv%alg  = sv%stringval(val)
  case default
    call sv%mld_d_base_solver_type%set(what,val,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info, name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine mld_d_ainv_solver_csetc

!  
!   
!                             MLD2P4  version 2.0
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.0)
!    
!    (C) Copyright 2008,2009,2010,2010,2012
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!                        Alfredo Buttari      CNRS-IRIT, Toulouse
!                        Pasqua D'Ambra       ICAR-CNR, Naples
!                        Daniela di Serafino  Second University of Naples
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
subroutine mld_d_invt_solver_clone(sv,svout,info)
  
  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_d_invt_solver_clone

  Implicit None

  ! Arguments
  class(mld_d_invt_solver_type), intent(inout)              :: sv
  class(mld_d_base_solver_type), allocatable, intent(inout) :: svout
  integer(psb_ipk_), intent(out)                            :: info
  ! Local variables
  integer(psb_ipk_) :: err_act


  info=psb_success_
  call psb_erractionsave(err_act)

  if (allocated(svout)) then
    call svout%free(info)
    if (info == psb_success_) deallocate(svout, stat=info)
  end if
  if (info == psb_success_) &
       & allocate(mld_d_invt_solver_type :: svout, stat=info)
  if (info /= 0) then 
    info = psb_err_alloc_dealloc_
    goto 9999 
  end if
  select type(svo => svout)
  type is (mld_d_invt_solver_type)
    svo%fill_in    = sv%fill_in
    svo%inv_fill   = sv%inv_fill
    svo%thresh     = sv%thresh
    svo%inv_thresh = sv%inv_thresh
    call psb_safe_ab_cpy(sv%d,svo%d,info)
    if (info == psb_success_) &
         & call sv%dv%clone(svo%dv,info)
    if (info == psb_success_) &
         & call sv%w%clone(svo%w,info)
    if (info == psb_success_) &
         & call sv%z%clone(svo%z,info)
    
  class default
    info = psb_err_internal_error_
  end select

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_invt_solver_clone

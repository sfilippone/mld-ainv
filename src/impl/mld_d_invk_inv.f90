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
subroutine mld_d_invk_inv(fill_in,i,row,rowlevs,heap,ja,irp,val,uplevs,&
     & nidx,idxs,info)
  
  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_d_invk_inv

  implicit none 

  ! Arguments
  type(psb_i_heap), intent(inout)     :: heap 
  integer, intent(in)                 :: i, fill_in
  integer, intent(inout)              :: nidx,info
  integer, intent(inout)              :: rowlevs(:)
  integer, allocatable, intent(inout) :: idxs(:)
  integer, intent(in)                 :: ja(:),irp(:),uplevs(:)
  real(psb_dpk_), intent(in)          :: val(:)
  real(psb_dpk_), intent(inout)       :: row(:)

  ! Local variables
  integer             :: k,j,lrwk,jj,lastk, iret
  real(psb_dpk_)      :: rwk


  info = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i

  !
  ! Do while there are indices to be processed
  !
  do
    ! Beware: (iret < 0) means that the heap is empty, not an error.
    call heap%get_first(k,iret) 
    if (iret < 0) then 
!!$        write(psb_err_unit,*) 'IINVK: ',i,' returning at ',lastk
      return
    end if

    ! 
    ! Just in case an index has been put on the heap more than once.
    !
    if (k == lastk) cycle

    lastk = k 
    nidx = nidx + 1
    if (nidx>size(idxs)) then 
      call psb_realloc(nidx+psb_heap_resize,idxs,info)
      if (info /= psb_success_) return
    end if
    idxs(nidx) = k

    if ((row(k) /= dzero).and.(rowlevs(k) <= fill_in)) then 
      !
      ! Note: since U is scaled while copying it out (see iluk_copyout),
      ! we can use rwk in the update below
      ! 
      rwk    = row(k)
      lrwk   = rowlevs(k)

      do jj=irp(k),irp(k+1)-1
        j = ja(jj)
        if (j<=k) then 
          info = -i
          return
        endif
        !
        ! Insert the index into the heap for further processing.
        ! The fill levels are initialized to a negative value. If we find
        ! one, it means that it is an as yet untouched index, so we need
        ! to insert it; otherwise it is already on the heap, there is no
        ! need to insert it more than once. 
        !
        if (rowlevs(j)<0) then 
          call heap%insert(j,info)
          if (info /= psb_success_) return
          rowlevs(j) = abs(rowlevs(j))
        end if
        !
        ! Update row(j) and the corresponding fill level
        !
        row(j)     = row(j) - rwk * val(jj)
        rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
      end do

    end if
  end do

end subroutine mld_d_invk_inv

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
subroutine mld_d_invt_inv(thres,i,nrmi,row,heap,irwt,ja,irp,val,nidx,idxs,info)
 
  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_d_invt_inv

  implicit none 

  ! Arguments
  type(psb_i_heap), intent(inout)     :: heap 
  integer(psb_ipk_), intent(in)       :: i
  integer(psb_ipk_), intent(inout)    :: nidx,info
  integer(psb_ipk_), intent(inout)    :: irwt(:) 
  real(psb_dpk_), intent(in)          :: thres,nrmi
  integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
  integer(psb_ipk_), intent(in)       :: ja(:),irp(:)
  real(psb_dpk_), intent(in)          :: val(:)
  real(psb_dpk_), intent(inout)       :: row(:)

  ! Local Variables
  integer(psb_ipk_) :: k,j,jj,lastk,iret
  real(psb_dpk_)    :: rwk, alpha

  info  = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i
  irwt(i) = 1 
!!$  write(0,*) 'Drop Threshold ',thres*nrmi
  !
  ! Do while there are indices to be processed
  !
  do

    call heap%get_first(k,iret) 
    if (iret < 0)  exit

    ! 
    ! An index may have been put on the heap more than once.
    ! Should not happen, but just in case. 
    !
    if (k == lastk)  cycle
    lastk = k 

    !
    ! Dropping rule based on the threshold: compare the absolute
    ! value of each updated entry of row with thres * 2-norm of row.
    !
    rwk    = row(k)

    if (abs(rwk) < thres*nrmi) then
      ! 
      ! Drop the entry.
      !
      row(k)  = dzero
      irwt(k) = 0
      cycle
    else
      !
      ! Note: since U is scaled while copying it out (see ilut_copyout),
      ! we can use rwk in the update below.
      !           
      do jj=irp(k),irp(k+1)-1
        j = ja(jj)
        if (j<=k) then 
          info = -i 
          return
        endif
        !
        ! Update row(j) and, if it is not to be discarded, insert
        ! its index into the heap for further processing.
        !
        row(j)     = row(j) - rwk * val(jj)
        if (irwt(j) == 0) then 
          if (abs(row(j)) < thres*nrmi) then
            ! 
            ! Drop the entry.
            !
            row(j)  = dzero
          else
            !
            ! Do the insertion.
            !
            call heap%insert(j,info)
            if (info /= psb_success_) return
            irwt(j) = 1
          end if
        end if
      end do
    end if

    !
    ! If we get here it is an index we need to keep on copyout.
    !

    nidx       = nidx + 1
    call psb_ensure_size(nidx,idxs,info,addsz=psb_heap_resize)      
    if (info /= psb_success_) return
    idxs(nidx) = k
    irwt(k)    = 0
  end do

  irwt(i) = 0
end subroutine mld_d_invt_inv

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
subroutine psb_d_dsc_print(iout,a,iv,head,ivr,ivc)
  
  use psb_string_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_print
  implicit none 

  integer, intent(in)               :: iout
  class(psb_d_dsc_sparse_mat), intent(in) :: a   
  integer, intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  Integer :: err_act
  character(len=20)  :: name='d_dsc_print'
  logical, parameter :: debug=.false.

  character(len=80)                 :: frmtv 
  integer  :: irs,ics,i,j, nmx, ni, nr, nc, nz

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

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
subroutine mld_dsparse_invt(n,a,z,nzrmax,sp_thresh,info)
  
  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name => mld_dsparse_invt

  implicit none 
  integer(psb_ipk_), intent(in)        :: n
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_dspmat_type), intent(inout) :: z
  integer(psb_ipk_), intent(in)        :: nzrmax
  real(psb_dpk_), intent(in)           :: sp_thresh
  integer(psb_ipk_), intent(out)       :: info
  !
  integer(psb_ipk_) :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2
  integer(psb_ipk_), allocatable :: ia(:), ja(:), iz(:),jz(:) 
  real(psb_dpk_), allocatable    :: zw(:), val(:), valz(:)
  integer(psb_ipk_), allocatable :: uplevs(:), rowlevs(:),idxs(:)
  real(psb_dpk_), allocatable :: row(:)
  type(psb_d_coo_sparse_mat)  :: trw
  type(psb_d_csr_sparse_mat)  :: acsr, zcsr
  integer(psb_ipk_)           :: ktrw, nidx, nlw,nup,jmaxup
  type(psb_i_heap)            :: heap
  real(psb_dpk_)     :: alpha, nrmi
  character(len=20)  :: name='mld_sp_invt'

  info = psb_success_  
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

  if (.not.(a%is_triangle().and.a%is_unit().and.a%is_upper())) then 
    write(psb_err_unit,*) 'Wrong A ' 
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='wrong A')
    goto 9999      
  end if
  call a%cp_to(acsr)
  call trw%allocate(izero,izero,ione)
  if (info == psb_success_) allocate(zw(n),iz(n),valz(n),&
       & row(n),rowlevs(n),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  call zcsr%allocate(n,n,n*nzrmax)
  call zcsr%set_triangle()
  call zcsr%set_unit(.false.)
  call zcsr%set_upper()
  ! 
  !
  nzz        = 0
  row(:)     = dzero 
  rowlevs(:) = 0
  l2         = 0
  zcsr%irp(1) = 1 

  outer: do i = 1, n-1
    ! ZW = e_i
    call mld_invt_copyin(i,n,acsr,i,ione,n,nlw,nup,jmaxup,nrmi,row,&
         & heap,rowlevs,ktrw,trw,info,sign=-done)
    if (info /= 0) exit
    row(i) = done
    ! Adjust norm
    if (nrmi < done) then 
      nrmi = sqrt(done + nrmi**2)
    else 
      nrmi = nrmi*sqrt(done+done/(nrmi**2))
    end if

    call mld_invt_inv(sp_thresh,i,nrmi,row,heap,rowlevs,&
         & acsr%ja,acsr%irp,acsr%val,nidx,idxs,info)
    if (info /= 0) exit
!!$    write(0,*) 'Calling copyout ',nzrmax,nlw,nup,nidx,l2
    call mld_invt_copyout(nzrmax,sp_thresh,i,n,nlw,nup,jmaxup,nrmi,row,&
         & nidx,idxs,l2,zcsr%ja,zcsr%irp,zcsr%val,info)
    if (info /= 0) exit
    nzz = l2
  end do outer
  if (info /= psb_success_) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='mainloop')
    goto 9999
  end if

  ipz1 = nzz+1
  call psb_ensure_size(ipz1,zcsr%val,info)
  call psb_ensure_size(ipz1,zcsr%ja,info)
  zcsr%val(ipz1) = done
  zcsr%ja(ipz1)  = n
  zcsr%irp(n+1)  = ipz1+1 

  call z%mv_from(zcsr)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine mld_dsparse_invt

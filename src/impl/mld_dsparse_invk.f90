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
subroutine mld_dsparse_invk(n,a,z,fill_in,info,inlevs)
  
  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_dsparse_invk

  integer(psb_ipk_), intent(in)           :: n
  type(psb_dspmat_type), intent(in)       :: a
  type(psb_dspmat_type), intent(inout)    :: z
  integer(psb_ipk_), intent(in)           :: fill_in
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), intent(in), optional :: inlevs(:)
  !
  integer(psb_ipk_) :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2 
  integer(psb_ipk_), allocatable        :: ia(:), ja(:), iz(:), jz(:) 
  real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
  integer(psb_ipk_), allocatable        :: uplevs(:), rowlevs(:), idxs(:)
  real(psb_dpk_), allocatable :: row(:)
  type(psb_d_coo_sparse_mat)  :: trw
  type(psb_d_csr_sparse_mat)  :: acsr, zcsr
  integer(psb_ipk_)           :: ktrw, nidx
  type(psb_i_heap)            :: heap

  real(psb_dpk_)     :: alpha
  character(len=20)  :: name='mld_sp_invk'

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

  allocate(uplevs(acsr%get_nzeros()),stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='Allocate')
    goto 9999
  end if
  uplevs(:)  = 0
  row(:)     = dzero
  rowlevs(:) = -(n+1)

  call zcsr%allocate(n,n,n*fill_in)
  call zcsr%set_triangle()
  call zcsr%set_unit(.false.)
  call zcsr%set_upper()
  call psb_ensure_size(n+1, idxs,  info)


  ! 
  !
  zcsr%irp(1)  = 1
  nzz          = 0

  l2 = 0
  outer: do i = 1, n-1
    ! ZW = e_i
    call mld_invk_copyin(i,n,acsr,ione,n,row,rowlevs,heap,ktrw,trw,info,&
         & sign=-done,inlevs=inlevs)
    row(i)     = done
    rowlevs(i) = 0

    ! Update loop
    call mld_invk_inv(fill_in,i,row,rowlevs,heap,&
         & acsr%ja,acsr%irp,acsr%val,uplevs,nidx,idxs,info)

    call mld_invk_copyout(fill_in,i,n,row,rowlevs,nidx,idxs,&
         & l2,zcsr%ja,zcsr%irp,zcsr%val,info)

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
  call zcsr%set_sorted()  
  call z%mv_from(zcsr)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine mld_dsparse_invk

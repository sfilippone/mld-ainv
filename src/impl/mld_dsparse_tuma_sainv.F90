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
subroutine mld_dsparse_tuma_sainv(n,a,p,z,w,nzrmax,sp_thresh,info)
  use psb_base_mod
  use mld_d_base_ainv_mod
  ! Interface to TUMA's code
  !
  implicit none 
  integer(psb_ipk_), intent(in)             :: n
  type(psb_d_csr_sparse_mat), intent(in)    :: a
  type(psb_d_csc_sparse_mat), intent(inout) :: z,w
  integer(psb_ipk_), intent(in)             :: nzrmax
  real(psb_dpk_), intent(in)                :: sp_thresh
  real(psb_dpk_), intent(out)               :: p(:)
  integer(psb_ipk_), intent(out)            :: info

  ! Locals
  type(psb_d_csr_sparse_mat)  :: ztum
  integer(psb_ipk_), pointer        :: ia(:), ja(:), iz(:),jz(:)
  real(psb_dpk_), pointer :: val(:), valz(:)
  integer(psb_ipk_) :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj, nza,&
       & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz, ipz1, ipz2
  integer(psb_ipk_) ::  msglvl,msgunit,size_r,size_c,size_p
  integer(psb_ipk_) garcol,garrow,droptyp
  integer(psb_ipk_) imodif,diag_one,fill,fillmax,ifillmax
  double precision mi,drfl,diagtol

  real(psb_dpk_)     :: alpha, t0, t1
  character(len=20)  :: name='mld_TUMA_sainv'

  interface 
    subroutine ainvsr2(msglvl,msgunit,n,ia,ja,a,ip,jp,ap,&
         &  size_p,size_c,size_r,diagtol,&
         &  drfl,mi,diag_one,droptyp,imodif,fill,fillmax,&
         &  ifillmax,garrow,garcol,info)
      integer msglvl,msgunit,n,size_r,size_c,size_p
      integer, intent(in) :: ia(*),ja(*)
      double precision, intent(in) :: a(*)
      integer, pointer :: ip(:),jp(:)
      double precision, pointer :: ap(:)
      integer garcol,garrow,droptyp
      integer imodif,diag_one,fill,fillmax,ifillmax,info
      double precision mi,drfl,diagtol
    end subroutine ainvsr2
  end interface

  info = psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

#if defined(HAVE_TUMA_SAINV) && defined(IPK4)
  !
  ! First step. 
  ! 
  nza = a%get_nzeros()
  allocate(ia(n+1),ja(nza),val(nza),stat=info)
  nullify(iz,jz,valz)
  if (info /= 0) then 
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='Allocate')
    goto 9999      
  end if

  msglvl  = 0
  msgunit = psb_err_unit
  drfl    = sp_thresh



  size_p  = nza
  size_r  = nza
  size_c  = nza
  ! These are taken straight from TUMA's code. 
  diagtol  = 1.1d-16
  droptyp  = 0
  mi       = 0.1d0
  diag_one = 1
  t0 = psb_wtime()
  call  ainvsr2(msglvl,msgunit,n,a%irp,a%ja,a%val,iz,jz,valz,&
       &  size_p,size_c,size_r,diagtol,&
       &  drfl,mi,diag_one,droptyp,imodif,fill,fillmax,&
       &  ifillmax,garrow,garcol,info)
  t1 = psb_wtime()
  nz=iz(n+1)-1
  ! !$    write(0,*) 'On output from AINVSR2 ',info,fillmax,a%get_nzeros(),iz(n+1)-1,t1-t0
  if (info /= 0) then 
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='ainvsr2')
    goto 9999      
  end if
  call ztum%allocate(n,n,nz)
  ztum%irp(1:n+1) = iz(1:n+1)
  ztum%ja(1:nz)   = jz(1:nz)
  ztum%val(1:nz)  = valz(1:nz)
  call ztum%transp(w)
  call w%cp_to_fmt(z,info)
  p = done
  deallocate(iz,jz,valz,stat=info)
  if (info /= 0) then 
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='ainvsr2')
    goto 9999      
  end if

#else 
  info = psb_err_from_subroutine_
  call psb_errpush(psb_err_internal_error_,name,a_err='sainv not linked')
  goto 9999      

#endif
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine mld_dsparse_tuma_sainv


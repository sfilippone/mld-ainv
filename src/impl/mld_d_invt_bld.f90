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

subroutine mld_d_invt_bld(a,fillin,invfill,thresh,invthresh,&
     & lmat,d,umat,desc,info,blck)

  use psb_base_mod
  use mld_d_invt_solver, mld_protect_name =>  mld_d_invt_bld
  use mld_d_ilu_fact_mod
  implicit none

  ! Arguments                                                     
  type(psb_dspmat_type), intent(in), target   :: a
  integer, intent(in)                         :: fillin,invfill
  real(psb_dpk_), intent(in)                  :: thresh
  real(psb_dpk_), intent(in)                  :: invthresh
  type(psb_dspmat_type), intent(inout)        :: lmat, umat
  real(psb_dpk_), allocatable                 :: d(:)
  Type(psb_desc_type), Intent(in)             :: desc
  integer, intent(out)                        :: info
  type(psb_dspmat_type), intent(in), optional :: blck
  integer   :: i, nztota, err_act, n_row, nrow_a, n_col
  type(psb_dspmat_type)          :: atmp
  real(psb_dpk_), allocatable :: pq(:), pd(:), w(:)
  integer   :: debug_level, debug_unit
  integer   :: ictxt,np,me
  integer            :: nzrmax
  real(psb_dpk_)     :: sp_thresh

  character(len=20)  :: name, ch_err, fname


  if(psb_get_errstatus() /= psb_success_) return 
  info = psb_success_
  name='mld_dainv_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = psb_cd_get_context(desc)
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  !
  ! Check the memory available to hold the incomplete L and U factors
  ! and allocate it if needed
  !
  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()

  if (present(blck)) then 
    nztota = nztota + blck%get_nzeros()
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ': out get_nnzeros',nrow_a,nztota,&
       & a%get_nrows(),a%get_ncols(),a%get_nzeros()


  n_row  = psb_cd_get_local_rows(desc)
  n_col  = psb_cd_get_local_cols(desc)
  allocate(pd(n_row),w(n_row),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  nzrmax    = fillin
  sp_thresh = thresh

  call lmat%csall(n_row,n_row,info,nz=nztota)
  if (info == psb_success_) call umat%csall(n_row,n_row,info,nz=nztota)

  if (info == 0) call mld_ilut_fact(nzrmax,sp_thresh,&
       & a,lmat,umat,pd,info,blck=blck,iscale=mld_ilu_scale_maxval_)

  if (info == psb_success_) call atmp%csall(n_row,n_row,info,nz=nztota)
  if(info/=0) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (.false.) then 
!!$    if (debug_level >= psb_debug_inner_) then 
    write(fname,'(a,i0,a)') 'invt-lo-',me,'.mtx'
    call lmat%print(fname,head="INVTLOW")
    write(fname,'(a,i0,a)') 'invt-up-',me,'.mtx'
    call umat%print(fname,head="INVTUPP")
  end if

  !
  ! Compute the approx U^-1  and L^-1
  !
  nzrmax    = invfill
  call mld_sparse_invt(n_row,umat,atmp,nzrmax,invthresh,info)
  if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
  if (info == psb_success_) call lmat%transp()
  if (info == psb_success_) call mld_sparse_invt(n_row,lmat,atmp,nzrmax,invthresh,info)
  if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
  if (info == psb_success_) call lmat%transp()
  ! Done. Hopefully.... 

  if (info /= psb_success_) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invt')
    goto 9999
  end if

  call psb_move_alloc(pd,d,info)
  call lmat%set_asb()
  call lmat%trim()
  call umat%set_asb()
  call umat%trim()

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_invt_bld


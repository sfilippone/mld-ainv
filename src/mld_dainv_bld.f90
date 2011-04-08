!!$ 
!!$ 
!!$                           MLD2P4  version 1.2
!!$  MultiLevel Domain Decomposition Parallel Preconditioners Package
!!$             based on PSBLAS (Parallel Sparse BLAS version 2.3.1)
!!$  
!!$  (C) Copyright 2008,2009
!!$
!!$                      Salvatore Filippone  University of Rome Tor Vergata
!!$                      Alfredo Buttari      University of Rome Tor Vergata
!!$                      Pasqua D'Ambra       ICAR-CNR, Naples
!!$                      Daniela di Serafino  Second University of Naples
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
! File: mld_dainv_bld.f90
!
! Subroutine: mld_dainv_bld
! Version:    real
!
!  For details on the above factorizations see
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!
!
! Arguments:
!


subroutine mld_dainv_bld(a,p,upd,info,blck)

  use psb_base_mod
  use mld_inner_mod, mld_protect_name => mld_dainv_bld
  use mld_dainv_mod
  implicit none

  ! Arguments                                                     
  type(psb_dspmat_type), intent(in), target   :: a
  type(mld_dbaseprec_type), intent(inout)     :: p
  character, intent(in)                       :: upd
  integer, intent(out)                        :: info
  type(psb_dspmat_type), intent(in), optional :: blck

  !     Local Variables                       
  integer   :: i, nztota, err_act, n_row, nrow_a, n_col, k
  integer   :: debug_level, debug_unit
  integer   :: ictxt,np,me
  character(len=20)  :: name, ch_err


  if(psb_get_errstatus() /= psb_success_) return 
  info = psb_success_
  name='mld_dainv_bld'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = psb_cd_get_context(p%desc_data)
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'


  select case (p%iprcparm(mld_sub_alg_))
  case (mld_ainv_orth1_,mld_ainv_orth2_,mld_ainv_orth3_,mld_ainv_orth4_)
    call mld_d_ainv_orth_bld(a,p,upd,info,blck)    
  case (mld_ainv_spai_,mld_ainv_invt_)
    call mld_d_ainv_invt_bld(a,p,upd,info,blck)    
  case (mld_ainv_invk_)
    call mld_d_ainv_invk_bld(a,p,upd,info,blck)    
  case default 
    info = psb_err_alloc_dealloc_ 
  end select

  if (info /= psb_success_) then 
    call psb_errpush(info,name)
    goto 9999
  end if

  n_row = psb_cd_get_local_rows(p%desc_data)
  n_col = psb_cd_get_local_cols(p%desc_data)
  nrow_a = a%m 
  ! The following is known to work 
  ! given that the output from CLIP is in COO. 
  call psb_sp_clip(a,p%av(mld_ap_nd_),info,&
       & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
!!$    if (info == psb_success_) call psb_sp_clip(blck_,atmp,info,&
!!$         & jmin=nrow_a+1,rscale=.false.,cscale=.false.)
!!$    if (info == psb_success_) call psb_rwextd(n_row,p%av(mld_ap_nd_),info,b=atmp) 
  if (info == psb_success_) call psb_spcnv(p%av(mld_ap_nd_),info,&
       & afmt='csr',dupl=psb_dupl_add_)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='clip & psb_spcnv csr 4')
    goto 9999
  end if
  
  k = psb_sp_get_nnzeros(p%av(mld_ap_nd_))
  call psb_sum(ictxt,k)
  
  if (k == 0) then 
    !
    ! If the off block-diagonal part is emtpy, there is no point in doing
    ! multiple Jacobi sweeps. This is certain to happen when running
    ! on a single processor.
    !
    p%iprcparm(mld_smoother_sweeps_) = 1
  end if
  if (info/=0) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_sp_free')
      goto 9999
  end if


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

contains

  subroutine mld_d_ainv_orth_bld(a,p,upd,info,blck)

    use psb_base_mod
    use mld_dainv_mod
    implicit none

    ! Arguments                                                     
    type(psb_dspmat_type), intent(in), target   :: a
    type(mld_dbaseprec_type), intent(inout)     :: p
    character, intent(in)                       :: upd
    integer, intent(out)                        :: info
    type(psb_dspmat_type), intent(in), optional :: blck
    integer   :: i, nztota, err_act, n_row, nrow_a
    type(psb_dspmat_type) :: wmat, zmat, atmp
    real(psb_dpk_), allocatable :: pq(:)
    integer   :: debug_level, debug_unit
    integer   :: ictxt,np,me
    integer            :: nzrmax, alg
    real(psb_dpk_)     :: sp_thresh

    character(len=20)  :: name, ch_err


    if(psb_get_errstatus() /= psb_success_) return 
    info = psb_success_
    name='mld_dainv_bld'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(p%desc_data)
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'


    !
    ! Check the memory available to hold the W and Z factors
    ! and allocate it if needed
    !
    nrow_a = psb_sp_get_nrows(a)
    nztota = psb_sp_get_nnzeros(a)
    if (present(blck)) then 
      nztota = nztota + psb_sp_get_nnzeros(blck)
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ': out get_nnzeros',nztota,a%m,a%k,nrow_a


    n_row  = psb_cd_get_local_rows(p%desc_data)
    allocate(pq(n_row),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    call psb_nullify_sp(wmat)
    call psb_nullify_sp(zmat)
    call psb_nullify_sp(atmp)
    call psb_sp_all(n_row,n_row,wmat,nztota,info)
    if (info == psb_success_) call psb_sp_all(n_row,n_row,zmat,nztota,info)
    if (info == psb_success_) call psb_sp_all(n_row,n_row,atmp,nztota,info)
    if(info/=0) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    alg       = p%iprcparm(mld_sub_alg_)
    nzrmax    = p%iprcparm(mld_sub_fillin_)
    sp_thresh = p%rprcparm(mld_sub_iluthrs_)
    !
    ! Ok, let's start first with W. 
    !
    call psb_transp(a,atmp,fmt='CSR') 
    call mld_sparse_orthbase(alg,n_row,atmp,pq,wmat,nzrmax,sp_thresh,info)
    call psb_transp(wmat)
    ! Now for Z
    if (info == psb_success_) call psb_transp(atmp,fmt='CSR') 
    if (info == psb_success_) call mld_sparse_orthbase(alg,n_row,atmp,pq,zmat,nzrmax,sp_thresh,info)

    ! Done. Hopefully.... 

    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='orthbase')
      goto 9999
    end if
    do i=1, n_row
      if (abs(pq(i)) < d_epstol) then
        pq(i) = done
      else
        pq(i) = done/pq(i)
      end if
    end do


    if (allocated(p%av)) then 
      if (size(p%av) < mld_bp_ilu_avsz_) then 
        do i=1,size(p%av) 
          call psb_sp_free(p%av(i),info)
          if (info /= psb_success_) then 
            ! Actually, we don't care here about this. Just let it go.
            ! return
          end if
        enddo
        deallocate(p%av,stat=info)
      endif
    end if


    if (.not.allocated(p%av)) then 
      allocate(p%av(mld_max_avsz_),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(4000,name)
        goto 9999
      end if
    endif

    call psb_move_alloc(pq,p%d,info)
    call psb_move_alloc(wmat,p%av(mld_ainv_w_),info)
    call psb_move_alloc(zmat,p%av(mld_ainv_z_),info)

    if (psb_sp_getifld(psb_upd_,p%av(mld_ainv_z_),info) /= psb_upd_perm_) then
      call psb_sp_trim(p%av(mld_ainv_z_),info)
    endif

    if (psb_sp_getifld(psb_upd_,p%av(mld_ainv_w_),info) /= psb_upd_perm_) then
      call psb_sp_trim(p%av(mld_ainv_w_),info)
    endif

!!$    call psb_csprt('z_orth.mtx',p%av(mld_ainv_z_))
!!$    call psb_csprt('w_orth.mtx',p%av(mld_ainv_w_))


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
  end subroutine mld_d_ainv_orth_bld

  subroutine mld_d_ainv_invt_bld(a,p,upd,info,blck)

    use psb_base_mod
    use mld_dainv_mod
    implicit none

    ! Arguments                                                     
    type(psb_dspmat_type), intent(in), target   :: a
    type(mld_dbaseprec_type), intent(inout)     :: p
    character, intent(in)                       :: upd
    integer, intent(out)                        :: info
    type(psb_dspmat_type), intent(in), optional :: blck
    integer   :: i, nztota, err_act, n_row, nrow_a, n_col
    type(psb_dspmat_type) :: lmat, umat, atmp
    real(psb_dpk_), allocatable :: pq(:), pd(:)
    integer   :: debug_level, debug_unit
    integer   :: ictxt,np,me
    integer            :: nzrmax
    real(psb_dpk_)     :: sp_thresh

    character(len=20)  :: name, ch_err


    if(psb_get_errstatus() /= psb_success_) return 
    info = psb_success_
    name='mld_dainv_bld'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(p%desc_data)
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'

    !
    ! Check the memory available to hold the incomplete L and U factors
    ! and allocate it if needed
    !
    nrow_a = psb_sp_get_nrows(a)
    nztota = psb_sp_get_nnzeros(a)
    if (present(blck)) then 
      nztota = nztota + psb_sp_get_nnzeros(blck)
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ': out get_nnzeros',nztota,a%m,a%k,nrow_a


    n_row  = psb_cd_get_local_rows(p%desc_data)
    n_col  = psb_cd_get_local_cols(p%desc_data)
    allocate(pq(n_row),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    call psb_nullify_sp(lmat)
    call psb_nullify_sp(umat)
    call psb_nullify_sp(atmp)
    nzrmax    = p%iprcparm(mld_sub_fillin_)
    sp_thresh = p%rprcparm(mld_sub_iluthrs_)


    call psb_sp_all(n_row,n_row,lmat,nztota,info)
    if (info == psb_success_) call psb_sp_all(n_row,n_row,umat,nztota,info)

    call mld_ilut_fact(nzrmax,sp_thresh,&
         & a,lmat,umat,pq,info,blck=blck)
!!$
!!$    call mld_iluk_fact(3,mld_ilu_n_,&
!!$         & a,lmat,umat,pq,info,blck=blck)

    if (info == psb_success_) call psb_sp_all(n_row,n_row,atmp,nztota,info)
    if(info/=0) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    !
    ! Compute the aprox U^-1  and L^-1
    !
    if (.true.) then 
      call mld_sparse_ainvt(n_row,umat,atmp,nzrmax,sp_thresh,info)
      if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
      if (info == psb_success_) call psb_transp(lmat)
      if (info == psb_success_) call mld_sparse_ainvt(n_row,lmat,atmp,nzrmax,sp_thresh,info)
      if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
      if (info == psb_success_) call psb_transp(lmat)
    else
      call mld_sparse_ainvk(n_row,umat,atmp,1,sp_thresh,info)
      if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
      if (info == psb_success_) call psb_transp(lmat)
      if (info == psb_success_) call mld_sparse_ainvk(n_row,lmat,atmp,1,sp_thresh,info)
      if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
      if (info == psb_success_) call psb_transp(lmat)
    endif
    ! Done. Hopefully.... 

    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='ainvt')
      goto 9999
    end if

    if (allocated(p%av)) then 
      if (size(p%av) < mld_bp_ilu_avsz_) then 
        do i=1,size(p%av) 
          call psb_sp_free(p%av(i),info)
          if (info /= psb_success_) then 
            ! Actually, we don't care here about this. Just let it go.
            ! return
          end if
        enddo
        deallocate(p%av,stat=info)
      endif
    end if


    if (.not.allocated(p%av)) then 
      allocate(p%av(mld_max_avsz_),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    call psb_move_alloc(pq,p%d,info)
    call psb_move_alloc(lmat,p%av(mld_l_pr_),info)
    call psb_move_alloc(umat,p%av(mld_u_pr_),info)

    if (psb_sp_getifld(psb_upd_,p%av(mld_l_pr_),info) /= psb_upd_perm_) then
      call psb_sp_trim(p%av(mld_l_pr_),info)
    endif

    if (psb_sp_getifld(psb_upd_,p%av(mld_u_pr_),info) /= psb_upd_perm_) then
      call psb_sp_trim(p%av(mld_u_pr_),info)
    endif


    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' end'
!!$
!!$    call psb_csprt('u_spai.mtx',p%av(mld_u_pr_))
!!$    call psb_csprt('l_spai.mtx',p%av(mld_l_pr_))


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mld_d_ainv_invt_bld


  subroutine mld_d_ainv_invk_bld(a,p,upd,info,blck)

    use psb_base_mod
    use mld_dainv_mod
    implicit none

    ! Arguments                                                     
    type(psb_dspmat_type), intent(in), target   :: a
    type(mld_dbaseprec_type), intent(inout)     :: p
    character, intent(in)                       :: upd
    integer, intent(out)                        :: info
    type(psb_dspmat_type), intent(in), optional :: blck
    integer   :: i, nztota, err_act, n_row, nrow_a, n_col, fill_in
    type(psb_dspmat_type) :: lmat, umat, atmp
    real(psb_dpk_), allocatable :: pq(:), pd(:)
    integer, allocatable :: uplevs(:)
    integer   :: debug_level, debug_unit
    integer   :: ictxt,np,me
    integer            :: nzrmax
    real(psb_dpk_)     :: sp_thresh

    character(len=20)  :: name, ch_err


    if(psb_get_errstatus() /= psb_success_) return 
    info = psb_success_
    name='mld_dainv_bld'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(p%desc_data)
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'

    !
    ! Check the memory available to hold the incomplete L and U factors
    ! and allocate it if needed
    !
    nrow_a = psb_sp_get_nrows(a)
    nztota = psb_sp_get_nnzeros(a)
    if (present(blck)) then 
      nztota = nztota + psb_sp_get_nnzeros(blck)
    end if
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ': out get_nnzeros',nztota,a%m,a%k,nrow_a


    n_row  = psb_cd_get_local_rows(p%desc_data)
    n_col  = psb_cd_get_local_cols(p%desc_data)
    allocate(pq(n_row),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    call psb_nullify_sp(lmat)
    call psb_nullify_sp(umat)
    call psb_nullify_sp(atmp)

    call psb_sp_all(n_row,n_row,lmat,nztota,info)
    if (info == psb_success_) call psb_sp_all(n_row,n_row,umat,nztota,info)

    fill_in   = p%iprcparm(mld_sub_fillin_)
    sp_thresh = p%rprcparm(mld_sub_iluthrs_)

    call mld_iluk_fact(fill_in/2,mld_ilu_n_,&
         & a,lmat,umat,pq,info,blck=blck,uplevs=uplevs)

    if (info == psb_success_) call psb_sp_all(n_row,n_row,atmp,nztota,info)
    if(info/=0) then
      info=psb_err_from_subroutine_  
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    

    !
    ! Compute the aprox U^-1  and L^-1
    !
    call mld_sparse_ainvk(n_row,umat,atmp,fill_in,sp_thresh,info,inlevs=uplevs)
    if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
    if (info == psb_success_) call psb_transp(lmat)
    if (info == psb_success_) call mld_sparse_ainvk(n_row,lmat,atmp,fill_in,&
         & sp_thresh,info)
    if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
    if (info == psb_success_) call psb_transp(lmat)
    
    ! Done. Hopefully.... 

    if (info /= psb_success_) then 
      info = psb_err_internal_error_ 
      call psb_errpush(info,name,a_err='ainvt')
      goto 9999
    end if

    if (allocated(p%av)) then 
      if (size(p%av) < mld_bp_ilu_avsz_) then 
        do i=1,size(p%av) 
          call psb_sp_free(p%av(i),info)
          if (info /= psb_success_) then 
            ! Actually, we don't care here about this. Just let it go.
            ! return
          end if
        enddo
        deallocate(p%av,stat=info)
      endif
    end if
    

    if (.not.allocated(p%av)) then 
      allocate(p%av(mld_max_avsz_),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    call psb_move_alloc(pq,p%d,info)
    call psb_move_alloc(lmat,p%av(mld_l_pr_),info)
    call psb_move_alloc(umat,p%av(mld_u_pr_),info)
    
    if (psb_sp_getifld(psb_upd_,p%av(mld_l_pr_),info) /= psb_upd_perm_) then
      call psb_sp_trim(p%av(mld_l_pr_),info)
    endif

    if (psb_sp_getifld(psb_upd_,p%av(mld_u_pr_),info) /= psb_upd_perm_) then
      call psb_sp_trim(p%av(mld_u_pr_),info)
    endif


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
  end subroutine mld_d_ainv_invk_bld



end subroutine mld_dainv_bld



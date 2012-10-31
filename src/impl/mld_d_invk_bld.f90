
subroutine mld_d_invk_bld(a,fill1, fill2,thresh,lmat,d,umat,desc,info,blck)

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name =>  mld_d_invk_bld
  use mld_d_ilu_fact_mod
  implicit none

  ! Arguments                                                     
  type(psb_dspmat_type), intent(in), target   :: a
  integer, intent(in)                         :: fill1, fill2 
  real(psb_dpk_), intent(in)                  :: thresh
  type(psb_dspmat_type), intent(inout)        :: lmat, umat
  real(psb_dpk_), allocatable                 :: d(:)
  Type(psb_desc_type), Intent(in)             :: desc
  integer, intent(out)                        :: info
  type(psb_dspmat_type), intent(in), optional :: blck
  integer   :: i, nztota, err_act, n_row, nrow_a, n_col
  type(psb_dspmat_type)          :: atmp
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
  allocate(pd(n_row),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  sp_thresh = thresh

  call lmat%csall(n_row,n_row,info,nz=nztota)
  if (info == psb_success_) call umat%csall(n_row,n_row,info,nz=nztota)


  call mld_iluk_fact(fill1,mld_ilu_n_,&
       & a,lmat,umat,pd,info,blck=blck)!,uplevs=uplevs)

  if (info == psb_success_) call atmp%csall(n_row,n_row,info,nz=nztota)
  if(info/=0) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  !
  ! Compute the aprox U^-1  and L^-1
  !
  call mld_sparse_invk(n_row,umat,atmp,fill2,sp_thresh,info)
  if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
  if (info == psb_success_) call lmat%transp()
  if (info == psb_success_) call mld_sparse_invk(n_row,lmat,atmp,fill2,sp_thresh,info)
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
end subroutine mld_d_invk_bld

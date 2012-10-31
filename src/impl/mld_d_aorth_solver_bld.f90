subroutine mld_d_aorth_solver_bld(a,desc_a,sv,upd,info,b,amold,vmold)
  

  use psb_base_mod
  use mld_d_aorth_solver, mld_protect_name => mld_d_aorth_solver_bld
  use mld_prec_mod
  use mld_d_aorth_bld_mod
  Implicit None

  ! Arguments
  type(psb_dspmat_type), intent(in), target  :: a
  Type(psb_desc_type), Intent(in)            :: desc_a
  class(mld_d_aorth_solver_type), intent(inout) :: sv
  character, intent(in)                      :: upd
  integer, intent(out)                       :: info
  type(psb_dspmat_type), intent(in), target, optional :: b
  class(psb_d_base_sparse_mat), intent(in), optional  :: amold
  class(psb_d_base_vect_type), intent(in), optional   :: vmold

  ! Local variables
  integer :: n_row,n_col, nrow_a, nztota
  real(psb_dpk_), pointer :: ww(:), aux(:), tx(:),ty(:)
  integer :: ictxt,np,me,i, err_act, debug_unit, debug_level
  character(len=20)  :: name='mld_d_aorth_solver_bld', ch_err

  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt       = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  call mld_ainv_orth_bld(a,sv%alg,sv%fill_in,sv%thresh,&
       & sv%w,sv%d,sv%z,desc_a,info,b,iscale=mld_ilu_scale_maxval_)    
  if ((info == psb_success_) .and.present(amold)) then 
    call sv%w%set_asb()
    call sv%w%trim()
    call sv%z%set_asb()
    call sv%z%trim()
    call sv%w%cscnv(info,mold=amold)
    if (info == psb_success_) &
         & call sv%z%cscnv(info,mold=amold)
  end if

  if (info == psb_success_) &
       & call sv%dv%bld(sv%d,mold=vmold)

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_internal_error_,name) 
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine mld_d_aorth_solver_bld

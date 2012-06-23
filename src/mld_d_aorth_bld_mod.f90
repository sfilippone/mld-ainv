module mld_d_aorth_bld_mod
  
  use mld_base_ainv_mod
  

  use mld_d_orthbase_mod
  use mld_d_biconjg_mod
  interface mld_ainv_orth_bld
    module procedure mld_d_ainv_orth_bld
  end interface mld_ainv_orth_bld


contains

  subroutine mld_d_ainv_orth_bld(a,alg,fillin,thresh,wmat,d,zmat,desc,info,blck,iscale)

    use psb_base_mod
    use mld_prec_mod
    implicit none

    ! Arguments                                                     
    type(psb_dspmat_type), intent(in), target   :: a
    integer, intent(in)                         :: fillin,alg
    real(psb_dpk_), intent(in)                  :: thresh
    type(psb_dspmat_type), intent(inout)        :: wmat, zmat
    real(psb_dpk_), allocatable                 :: d(:)
    Type(psb_desc_type), Intent(in)             :: desc
    integer, intent(out)                        :: info
    type(psb_dspmat_type), intent(in), optional :: blck
    integer, intent(in), optional               :: iscale
    integer   :: i, nztota, err_act, n_row, nrow_a
    type(psb_d_coo_sparse_mat)  :: acoo
    type(psb_d_csr_sparse_mat)  :: acsr
    real(psb_dpk_), allocatable :: pq(:), arws(:), acls(:), ad(:)
    integer            :: debug_level, debug_unit
    integer            :: ictxt,np,me
    integer            :: nzrmax, iscale_
    real(psb_dpk_)     :: sp_thresh, weight
    character(len=20)  :: name, ch_err


    if (psb_get_errstatus() /= psb_success_) return 
    info = psb_success_
    name = 'mld_dainv_bld'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    ictxt       = psb_cd_get_context(desc)
    call psb_info(ictxt, me, np)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'

    iscale_ = mld_ilu_scale_none_
    if (present(iscale)) iscale_ = iscale
    weight = done
    !
    ! Check the memory available to hold the W and Z factors
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


    n_row  = desc%get_local_rows()
    allocate(pq(n_row),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    nzrmax    = fillin
    sp_thresh = thresh

    !
    ! Ok, let's start first with Z (i.e. Upper) 
    !
    call a%csclip(acoo,info,imax=n_row,jmax=n_row)
    call acsr%mv_from_coo(acoo,info)
    select case(iscale_)
    case(mld_ilu_scale_none_) 
      ! Ok, do nothing.

    case(mld_ilu_scale_maxval_) 
      weight = acsr%maxval()
      weight = done/weight
      call acsr%scal(weight,info)

    case(mld_ilu_scale_arcsum_) 
      allocate(arws(n_row),acls(n_row),ad(n_row),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
      call acsr%arwsum(arws)
      call acsr%aclsum(acls)
      ad(1:n_row) = sqrt(sqrt(arws(1:n_row)*acls(1:n_row)))
      ad(1:n_row) = done/ad(1:n_row)
      call acsr%scal(ad,info,side='L')
      call acsr%scal(ad,info,side='R')
    case default
      call psb_errpush(psb_err_from_subroutine_,name,a_err='wrong iscale')
      goto 9999      
    end select

    if (alg /=  mld_ainv_orth3_) then 
      call mld_sparse_orthbase(alg,n_row,acsr,pq,&
           & zmat,nzrmax,sp_thresh,info)
      ! Now for W  (i.e. Lower) 
      if (info == psb_success_) call acsr%transp() 
      if (info == psb_success_) &
           & call mld_sparse_orthbase(alg,n_row,acsr,pq,&
           &   wmat,nzrmax,sp_thresh,info)
      call wmat%transp()
    else
      call mld_sparse_biconjg(alg,n_row,acsr,pq,&
           &   zmat,wmat,nzrmax,sp_thresh,info)
    end if
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

    select case(iscale_)
    case(mld_ilu_scale_none_) 
      ! Ok, do nothing.
    case(mld_ilu_scale_maxval_) 
      pq(:) = pq(:)*weight

    case(mld_ilu_scale_arcsum_) 
      call zmat%scal(ad,info,side='L')
      call wmat%scal(ad,info,side='R') 

    case default
      call psb_errpush(psb_err_from_subroutine_,name,a_err='wrong iscale')
      goto 9999      
    end select

    call psb_move_alloc(pq,d,info)

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

end module mld_d_aorth_bld_mod

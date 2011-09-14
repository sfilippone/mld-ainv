module mld_d_ainv_bld_mod
  
  use mld_base_ainv_mod
  
  use mld_d_ilu_fact_mod
  
  interface mld_ainv_invk_bld
    module procedure mld_d_ainv_invk_bld
  end interface mld_ainv_invk_bld


  interface  mld_invk_copyin
    subroutine mld_dinvk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,&
         & ktrw,trw,info,sign,inlevs)

      use psb_base_mod, only : psb_d_csr_sparse_mat, psb_d_coo_sparse_mat,&
           & psb_dpk_, psb_int_heap
      implicit none

      ! Arguments 
      type(psb_d_csr_sparse_mat), intent(in)    :: a
      type(psb_d_coo_sparse_mat), intent(inout) :: trw
      integer, intent(in)                  :: i,m,jmin,jmax
      integer, intent(inout)               :: ktrw,info
      integer, intent(inout)               :: rowlevs(:)
      real(psb_dpk_), intent(inout)        :: row(:)
      type(psb_int_heap), intent(inout)    :: heap
      real(psb_dpk_), optional, intent(in) :: sign
      integer, intent(in), optional        :: inlevs(:)

    end subroutine mld_dinvk_copyin
  end interface

  interface mld_invk
    subroutine mld_dinvk(fill_in,i,row,rowlevs,heap,uia1,uia2,uaspk,uplevs,nidx,idxs,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      implicit none 

      ! Arguments
      type(psb_int_heap), intent(inout)    :: heap 
      integer, intent(in)                  :: i, fill_in
      integer, intent(inout)               :: nidx,info
      integer, intent(inout)               :: rowlevs(:)
      integer, allocatable, intent(inout)  :: idxs(:)
      integer, intent(in)                  :: uia1(:),uia2(:),uplevs(:)
      real(psb_dpk_), intent(in)           :: uaspk(:)
      real(psb_dpk_), intent(inout)        :: row(:)


    end subroutine mld_dinvk
  end interface

  interface mld_invk_copyout
    subroutine mld_dinvk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
         &  l2,uia1,uia2,uaspk,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap

      implicit none 

      ! Arguments
      integer, intent(in)                        :: fill_in, i, m, nidx
      integer, intent(inout)                     :: l2, info
      integer, intent(inout)                     :: rowlevs(:), idxs(:)
      integer, allocatable, intent(inout)        :: uia1(:), uia2(:)
      real(psb_dpk_), allocatable, intent(inout) :: uaspk(:)
      real(psb_dpk_), intent(inout)              :: row(:)

    end subroutine mld_dinvk_copyout
  end interface
  
  interface mld_sparse_ainvk
    subroutine mld_dsparse_ainvk(n,a,z,fill_in,sp_thresh,info,inlevs)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      integer, intent(in)                  :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: fill_in
      real(psb_dpk_), intent(in)           :: sp_thresh
      integer, intent(out)                 :: info
      integer, intent(in), optional        :: inlevs(:)

    end subroutine mld_dsparse_ainvk
  end interface

  interface mld_ainv_invt_bld
    module procedure mld_d_ainv_invt_bld
  end interface mld_ainv_invt_bld
  

  interface mld_invt_copyin
    subroutine mld_dinvt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
         & irwt,ktrw,trw,info,sign)
      use psb_base_mod, only : psb_d_csr_sparse_mat, psb_d_coo_sparse_mat,&
           & psb_dpk_, psb_int_heap
      implicit none 
      type(psb_d_csr_sparse_mat), intent(in)    :: a
      type(psb_d_coo_sparse_mat), intent(inout) :: trw
      integer, intent(in)                  :: i, m,jmin,jmax,jd
      integer, intent(inout)               :: ktrw,nlw,nup,jmaxup,info
      integer, intent(inout)               :: irwt(:)
      real(psb_dpk_), intent(inout)        :: nrmi,row(:)
      type(psb_int_heap), intent(inout)    :: heap
      real(psb_dpk_), intent(in), optional :: sign

    end subroutine mld_dinvt_copyin
  end interface

  interface mld_invt
    subroutine mld_dinvt(thres,i,nrmi,row,heap,irwt,uia1,uia2,uaspk,nidx,idxs,info)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      implicit none 
      ! Arguments
      type(psb_int_heap), intent(inout)   :: heap 
      integer, intent(in)                 :: i
      integer, intent(inout)              :: nidx,info
      integer, intent(inout)              :: irwt(:) 
      real(psb_dpk_), intent(in)          :: thres,nrmi
      integer, allocatable, intent(inout) :: idxs(:)
      integer, intent(in)                 :: uia1(:),uia2(:)
      real(psb_dpk_), intent(in)          :: uaspk(:)
      real(psb_dpk_), intent(inout)       :: row(:)

    end subroutine mld_dinvt
  end interface

  interface mld_invt_copyout
    subroutine mld_dinvt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
         & nidx,idxs,l2,uia1,uia2,uaspk,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap

      implicit none 

      ! Arguments
      integer, intent(in)                       :: fill_in,i,m,nidx,nlw,nup,jmaxup
      integer, intent(in)                       :: idxs(:)
      integer, intent(inout)                    :: l2, info
      integer, allocatable, intent(inout)       :: uia1(:),uia2(:)
      real(psb_dpk_), intent(in)                :: thres,nrmi
      real(psb_dpk_),allocatable, intent(inout) :: uaspk(:)
      real(psb_dpk_), intent(inout)             :: row(:)

    end subroutine mld_dinvt_copyout
  end interface

  interface  mld_sparse_ainvt
    subroutine mld_dsparse_ainvt(n,a,z,nzrmax,sp_thresh,info)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      implicit none 
      integer, intent(in)                  :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: nzrmax
      real(psb_dpk_), intent(in)           :: sp_thresh
      integer, intent(out)                 :: info
    end subroutine mld_dsparse_ainvt
  end interface



contains

  subroutine mld_d_ainv_invt_bld(a,fillin,invfill,thresh,lmat,d,umat,desc,info,blck)

    use psb_base_mod

    implicit none

    ! Arguments                                                     
    type(psb_dspmat_type), intent(in), target   :: a
    integer, intent(in)                         :: fillin,invfill
    real(psb_dpk_), intent(in)                  :: thresh
    type(psb_dspmat_type), intent(inout)        :: lmat, umat
    real(psb_dpk_), allocatable                 :: d(:)
    Type(psb_desc_type), Intent(in)             :: desc
    integer, intent(out)                        :: info
    type(psb_dspmat_type), intent(in), optional :: blck
    integer   :: i, nztota, err_act, n_row, nrow_a, n_col
    type(psb_dspmat_type)          :: atmp
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

    nzrmax    = fillin
    sp_thresh = thresh

    call lmat%csall(n_row,n_row,info,nz=nztota)
    if (info == psb_success_) call umat%csall(n_row,n_row,info,nz=nztota)

    call mld_ilut_fact(nzrmax,sp_thresh,&
         & a,lmat,umat,pd,info,blck=blck)

    if (info == psb_success_) call atmp%csall(n_row,n_row,info,nz=nztota)
    if(info/=0) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    !
    ! Compute the approx U^-1  and L^-1
    !
    nzrmax    = invfill
    call mld_sparse_ainvt(n_row,umat,atmp,nzrmax,sp_thresh,info)
    if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
    if (info == psb_success_) call lmat%transp()
    if (info == psb_success_) call mld_sparse_ainvt(n_row,lmat,atmp,nzrmax,sp_thresh,info)
    if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
    if (info == psb_success_) call lmat%transp()
    ! Done. Hopefully.... 

    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='ainvt')
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
  end subroutine mld_d_ainv_invt_bld


  subroutine mld_d_ainv_invk_bld(a,fill1, fill2,thresh,lmat,d,umat,desc,info,blck)

    use psb_base_mod

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
    call mld_sparse_ainvk(n_row,umat,atmp,fill2,sp_thresh,info)
    if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
    if (info == psb_success_) call lmat%transp()
    if (info == psb_success_) call mld_sparse_ainvk(n_row,lmat,atmp,fill2,sp_thresh,info)
    if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
    if (info == psb_success_) call lmat%transp()
    ! Done. Hopefully.... 

    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='ainvt')
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
  end subroutine mld_d_ainv_invk_bld


end module mld_d_ainv_bld_mod

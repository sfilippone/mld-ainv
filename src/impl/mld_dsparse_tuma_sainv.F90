  
  subroutine mld_dsparse_tuma_sainv(n,a,p,z,w,nzrmax,sp_thresh,info)
    use psb_base_mod
    use mld_d_base_ainv_mod
    ! Interface to TUMA's code
      !
    implicit none 
    integer, intent(in)                       :: n
    type(psb_d_csr_sparse_mat), intent(in)    :: a
    type(psb_d_csc_sparse_mat), intent(inout) :: z,w
    integer, intent(in)                       :: nzrmax
    real(psb_dpk_), intent(in)                :: sp_thresh
    real(psb_dpk_), intent(out)               :: p(:)
    integer, intent(out)                      :: info

    
    ! Locals
    type(psb_d_csr_sparse_mat)  :: ztum
    integer, pointer        :: ia(:), ja(:), iz(:),jz(:)
    real(psb_dpk_), pointer :: val(:), valz(:)
    integer :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj, nza,&
         & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz, ipz1, ipz2
    integer ::  msglvl,msgunit,size_r,size_c,size_p
    integer garcol,garrow,droptyp
    integer imodif,diag_one,fill,fillmax,ifillmax
    double precision mi,drfl,diagtol

    type(psb_int_heap) :: heap 
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



    if (psb_get_errstatus() /= psb_success_) return 
    info = psb_success_
    call psb_erractionsave(err_act)

#ifdef HAVE_TUMA_SAINV
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

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine mld_dsparse_tuma_sainv


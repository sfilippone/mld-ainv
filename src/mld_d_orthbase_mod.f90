module mld_d_orthbase_mod

  interface mld_sparse_orthbase
    module procedure mld_dsparse_orthbase
  end interface mld_sparse_orthbase

contains

  subroutine mld_dsparse_orthbase(alg,n,a,p,z,nzrmax,sp_thresh,info)
    use psb_base_mod
    use mld_base_ainv_mod
    integer, intent(in)                    :: alg,n
    type(psb_d_csr_sparse_mat), intent(in) :: a
    type(psb_dspmat_type), intent(out)     :: z
    integer, intent(in)                    :: nzrmax
    real(psb_dpk_), intent(in)             :: sp_thresh
    real(psb_dpk_), intent(out)            :: p(:)
    integer, intent(out)                   :: info

    type(psb_d_csc_sparse_mat)             :: zcsc
    integer :: err_act
    character(len=20)  :: name='mld_sp_orthbase'
    integer, parameter :: variant=1


    if (psb_get_errstatus() /= psb_success_) return 
    info = psb_success_
    call psb_erractionsave(err_act)

    if (size(p)<n) then 
      write(psb_err_unit,*) 'Size of P wrong'
      info = psb_err_internal_error_
      call psb_errpush(psb_err_internal_error_,name,a_err='Allocate')
      goto 9999      
    end if

    select case(alg)
    case (mld_ainv_orth1_)
      call mld_dsparse_orth_llk(n,a,p,zcsc,nzrmax,sp_thresh,info)
      if (info /= 0) goto 9999
    case (mld_ainv_orth2_) 
      call mld_dsparse_orth_rlk(n,a,p,zcsc,nzrmax,sp_thresh,info)
!!$    if (info /= 0) goto 9999
!!$  case (mld_ainv_orth3_) 
!!$    call mld_dsparse_orth_3(n,a,p,z,nzrmax,sp_thresh,info)
!!$    if (info /= 0) goto 9999
!!$  case (mld_ainv_orth4_) 
!!$    call mld_dsparse_orth_4(n,a,p,zcsc,nzrmax,sp_thresh,info)
!!$    if (info /= 0) goto 9999
    case default
      info = psb_err_internal_error_
      call psb_errpush(psb_err_internal_error_,name,a_err='Invalid alg')
      goto 9999      
    end select

    call z%mv_from(zcsc)
    call z%cscnv(info,type='CSR')

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine mld_dsparse_orthbase

  subroutine mld_dsparse_orth_llk(n,a,p,z,nzrmax,sp_thresh,info)
    use psb_base_mod
    use mld_base_ainv_mod
    !
    ! Need to figure out how I am doing this.... 
    !
    !
    implicit none 
    integer, intent(in)                       :: n
    type(psb_d_csr_sparse_mat), intent(in)    :: a
    type(psb_d_csc_sparse_mat), intent(inout) :: z
    integer, intent(in)                       :: nzrmax
    real(psb_dpk_), intent(in)                :: sp_thresh
    real(psb_dpk_), intent(out)               :: p(:)
    integer, intent(out)                      :: info

    ! Locals
    integer, allocatable        :: ia(:), ja(:), iz(:),jz(:), icr(:), ikr(:), ljr(:)
    real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
    integer :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
         & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz, ipz1, ipz2
    type(psb_int_heap) :: heap 
    real(psb_dpk_)     :: alpha
    character(len=20)  :: name='mld_orth_llk'

    allocate(zw(n),iz(n),valz(n),icr(n),ikr(n),ljr(n),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      return      
    end if

    !
    ! Init pointers to:
    !  ljr(i): last occupied column index within row  I
    !  icr(i): first occupied row index within column I
    !
    do i=1,n
      icr(i) = n
      ikr(i) = 0
      ljr(n) = 0
    end do
    do i=1, n
      do j=a%irp(i),a%irp(i+1)-1
        k      = a%ja(j)
        if (k<=n) then 
          ljr(i) = max(ljr(i),k)
          icr(k) = min(icr(k),i)
        end if
      end do
    end do

    ! Init z_1=e_1 and p_1=a_11
    p(1) = dzero
    i   = 1
    nz  = a%irp(i+1) - a%irp(i)
    do j=1,nz
      if (a%ja(j) == 1) then 
        p(1) = a%val(j)
        exit
      end if
    end do
    if (abs(p(1)) < d_epstol) &
         & p(1) = 1.d-3 

    ! 
    !
    call z%allocate(n,n,n*nzrmax)

    z%icp(1)  = 1
    z%icp(2)  = 2
    z%ia(1)  = 1
    z%val(1) = done
    nzz       = 1

    do i = 2, n
      ! ZW = e_i
      ! !$        do j=1, i-1
      ! !$          zw(j) = dzero
      ! !$        end do
      zw(i)  = done
      ikr(i) = 1
      ifnz   = i
      call psb_init_heap(heap,info)
      if (info == psb_success_) call psb_insert_heap(i,heap,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_init_heap')
        return
      end if

      ! Update loop
      do j = icr(i), i-1
        if (ljr(j) < ifnz) then 
          ! !$          write(psb_err_unit,*) 'Cycling '
          cycle
        end if
        ip1 = a%irp(j)
        ip2 = a%irp(j+1) - 1
        do 
          if (ip2 < ip1 ) exit
          if (a%ja(ip2) <= n) exit
          ip2 = ip2 -1 
        end do
        nzra = max(0,ip2 - ip1 + 1) 
        p(i) = psb_spge_dot(nzra,a%ja(ip1:ip2),a%val(ip1:ip2),zw)
        ! !$          write(psb_err_unit,*) j,i,p(i)

        ipz1 = z%icp(j) 
        ipz2 = z%icp(j+1) 
        nzrz = ipz2-ipz1
        alpha = (-p(i)/p(j))
        if (abs(alpha) > sp_thresh) then 
          ! !$          call psb_aspxpby(alpha,nzrz,&
          ! !$               & z%ia(ipz1:ipz2-1),z%val(ipz1:ipz2-1), &
          ! !$               & done,zw,info)
          do k=ipz1, ipz2-1
            kr     = z%ia(k)
            zw(kr) = zw(kr) + alpha*z%val(k)
            ifnz   = min(ifnz,kr)
            if (ikr(kr) == 0) then 
              call psb_insert_heap(kr,heap,info) 
              ikr(kr) = 1
            end if
            if (info /= psb_success_) exit
          end do
          if (info /= psb_success_) then
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='psb_insert_heap')
            return
          end if
        end if
      end do
      call a%csget(i,i,nzra,ia,ja,val,info)
      call rwclip(nzra,ia,ja,val,1,n,1,n)      
      p(i) = psb_spge_dot(nzra,ja,val,zw)
      ! !$      write(psb_err_unit,*) i,i,p(i)
      !
      ! Sparsify current ZW and put into ZMAT
      ! 
      call sparsify(i,nzrmax,sp_thresh,n,zw,nzrz,iz,valz,info,iheap=heap,ikr=ikr)
      if (info /= psb_success_) then 
        info = psb_err_internal_error_
        call psb_errpush(info,name,a_err='sparsify')
        return
      end if
      call psb_ensure_size(nzz+nzrz, z%ia,  info)
      call psb_ensure_size(nzz+nzrz, z%val, info)
      ipz1 = z%icp(i)
      do j=1, nzrz
        z%ia(ipz1  + j -1) = iz(j)
        z%val(ipz1 + j -1) = valz(j)
      end do
      z%icp(i+1) = ipz1 + nzrz
      nzz        = nzz + nzrz
    end do

  end subroutine mld_dsparse_orth_llk

  subroutine mld_dsparse_orth_rlk(n,a,p,z,nzrmax,sp_thresh,info)
    use psb_base_mod
    use mld_base_ainv_mod
    use psb_d_dsc_mat_mod
    !
    ! Benzi-Tuma (98): alg biconjugation (section 4).
    !  dds implementation. Is this really what they claim it is?? 
    !

    implicit none 
    integer, intent(in)                  :: n
    type(psb_d_csr_sparse_mat), intent(in)    :: a
    type(psb_d_csc_sparse_mat), intent(inout) :: z
    integer, intent(in)                  :: nzrmax
    real(psb_dpk_), intent(in)           :: sp_thresh
    real(psb_dpk_), intent(out)          :: p(:)
    integer, intent(out)                 :: info
    integer, allocatable        :: ia(:), ja(:), iz(:),jz(:)
    real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
    integer :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
         & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz, ljr
    integer :: debug_unit, debug_level, me
    type(psb_int_heap) :: heap 
    real(psb_dpk_)     :: alpha
    type(psb_d_dsc_sparse_mat) :: zmat
    character(len=20)  :: name='mld_orth_rlk'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    me          = -1 
    ! !$    debug_level = psb_debug_outer_
    allocate(iz(n),jz(n),valz(n),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      return      
    end if

    ! Init z_1=e_1 and p_1=a_11
    
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' start'

    ljr = 0
    call zmat%allocate(n,n,n*nzrmax)
    ! Init Z to identity 
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' Init ZMAT'

    do i=1, n
      zmat%cols(i)%nz = 1
      zmat%cols(i)%idx(1) = i
      zmat%cols(i)%val(1) = done
    end do

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' Start main loop'

    do i = 1, n
      
      ipza = a%irp(i)
      nzra = a%irp(i+1) - ipza

      ! !$      write(psb_err_unit,*) 'Into sparse dot: ',nzra,a%ja(ipza:ipza+nzra-1),a%val(ipza:ipza+nzra-1)
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),' Main loop',i

      ljr = a%ja(a%irp(i+1)-1)

      do j=i, n
        if (ljr < zmat%cols(j)%idx(1)) then 
          p(j) = dzero
          cycle
        end if
        nzrz = zmat%cols(j)%nz
        p(j) = psb_spdot_srtd(nzra,a%ja(ipza:ipza+nzra-1),a%val(ipza:ipza+nzra-1),&
             & nzrz,zmat%cols(j)%idx,zmat%cols(j)%val)
      end do
      if (p(i) == dzero) then 
        write(psb_err_unit,*) 'Breakdown!! '
        p(i) = 1.d-3
      end if

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),' p(i)',p(i)

      nzzi = zmat%cols(i)%nz
      do j=i+1,n 
        ! New kernel psb_aspxpbspy 
        alpha = (-p(j)/p(i))
        if (debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' inner iteration ',j,' :',nzzi,nzzj,alpha,sp_thresh

        if (abs(alpha) > sp_thresh) then 

          nzzj  = zmat%cols(j)%nz
          call psb_nspaxpby(nzrz,iz,valz,&
               & alpha, nzzi, zmat%cols(i)%idx, zmat%cols(i)%val,&
               & done, nzzj, zmat%cols(j)%idx, zmat%cols(j)%val,&
               & info)
          if (debug_level >= psb_debug_outer_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' done nspaxpby',nzrz,size(iz),size(valz)
        
          ! Drop tolerance for new column 
          if (info == psb_success_) call sp_drop(j,nzrmax,sp_thresh,nzrz,iz,valz,info)
          if (debug_level >= psb_debug_outer_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' done drop', nzrz
          ! Copy into znew(j) 
          if (info == psb_success_) call psb_ensure_size(nzrz,zmat%cols(j)%idx,info)
          if (info == psb_success_) call psb_ensure_size(nzrz,zmat%cols(j)%val,info)
          if (info /= psb_success_) then 
            info = psb_err_internal_error_
            call psb_errpush(psb_err_internal_error_,name,a_err='Inner loop ')
            return      
          end if
          zmat%cols(j)%nz          = nzrz
          do k=1,nzrz
            zmat%cols(j)%idx(k) = iz(k)
            zmat%cols(j)%val(k) = valz(k)
          end do
        end if
        if (debug_level >= psb_debug_outer_) &
             & write(debug_unit,*) me,' ',trim(name),' completed inner iteration ',j

        
      end do
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),' completed outer iteration ',i
      if (info /= psb_success_) then 
        info = psb_err_internal_error_
        call psb_errpush(psb_err_internal_error_,name,a_err='Outer loop ')
        return      
      end if
    end do
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' End main loop'

    call zmat%mv_to_fmt(z,info)
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),' end'

  end subroutine mld_dsparse_orth_rlk

!!$  subroutine mld_dsparse_orth_3(n,a,p,z,nzrmax,sp_thresh,info)
!!$    !
!!$    ! Benzi-Tuma (98): alg biconjugation (section 4).
!!$    !  Better implementation, with mixed dinamic data. 
!!$    !
!!$    implicit none 
!!$    integer, intent(in)                  :: n
!!$    type(psb_dspmat_type), intent(in)    :: a
!!$    type(psb_dspmat_type), intent(inout) :: z
!!$    integer, intent(in)                  :: nzrmax
!!$    real(psb_dpk_), intent(in)           :: sp_thresh
!!$    real(psb_dpk_), intent(out)          :: p(:)
!!$    integer, intent(out)                 :: info
!!$    integer, allocatable        :: iv(:), icpz(:), ncz(:),irwz(:)
!!$    real(psb_dpk_), allocatable :: vv(:), valz(:), zw(:)
!!$    integer :: i,j,k, kc, kr, err_act, nz, nv,&
!!$         & ki,kj, nzralc, ipza, nzra
!!$    type(psb_dspmat_type) :: znew, ztmp, acsc
!!$    type(psb_int_heap) :: heap 
!!$    real(psb_dpk_)     :: alpha
!!$
!!$    
!!$    nzralc = nint(2.5*nzrmax)
!!$    write(0,*) 'fill in nzr',nzrmax,nzralc
!!$    allocate(iv(n),vv(n),icpz(n+1),ncz(n),&
!!$         & irwz(n*nzralc),valz(n*nzralc),stat=info)
!!$    if (info /= psb_success_) then 
!!$      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
!!$      return      
!!$    end if
!!$
!!$! !$    call psb_sp_clip(a,acsc,info,imax=n,jmax=n,&
!!$! !$         & rscale=.false.,cscale=.false.)
!!$! !$    call psb_transp(acsc,fmt='CSC')
!!$! !$  
!!$
!!$    !
!!$    ! Init Z to identity.
!!$    !  Each column of Z has room for 2.5 nzrmax, should be plenty
!!$    !  to handle temporaries. 
!!$    !
!!$    j = 1
!!$    do i=1, n
!!$      icpz(i)  = j
!!$      ncz(i)   = 1 
!!$      irwz(j)  = i
!!$      valz(j)  = done
!!$      j        = j + nzralc
!!$    end do
!!$
!!$  
!!$
!!$    do i = 1, n
!!$! !$      ! First, compute v_i=Az_i
!!$      k  = icpz(i)
!!$      nv = 0 
!!$! !$      call psb_d_spmspv(done,acsc,ncz(i),irwz(k:),valz(k:),dzero,nv,iv,vv, info)
!!$! !$      if (info /=0) then 
!!$! !$        write(0,*) 'Internal error from spmspv ',info
!!$! !$        return
!!$! !$      end if
!!$      ! Second, compute the P(J)
!!$! !$      write(psb_err_unit,*) 'Into sparse dot: ',nv,iv(1:nv),vv(1:nv)
!!$      ipza = a%irp(i)
!!$      nzra = a%irp(i+1) - ipza
!!$      k    = icpz(i)
!!$      p(i) = psb_spdot_srtd(nzra,a%ja(ipza:ipza+nzra-1),a%val(ipza:ipza+nzra-1),&
!!$           & ncz(i),irwz(k:),valz(k:))
!!$      
!!$      if (p(i) == dzero) then 
!!$        write(psb_err_unit,*) 'Breakdown!! '
!!$        p(i) = 1.d-3
!!$      end if
!!$
!!$      ! Third, right-looking update. 
!!$      ki = icpz(i) 
!!$    
!!$      do j=i+1,n 
!!$        kj   = icpz(j)
!!$        p(j) = psb_spdot_srtd(nzra,a%ja(ipza:ipza+nzra-1),a%val(ipza:ipza+nzra-1),&
!!$             & ncz(j),irwz(kj:),valz(kj:))
!!$
!!$        alpha = (-p(j)/p(i))
!!$
!!$        if (abs(alpha) > sp_thresh) then 
!!$
!!$          call psb_nspaxpby(nv, iv, vv,&
!!$               & alpha, ncz(i), irwz(ki:), valz(ki:),&
!!$               & done, ncz(j), irwz(kj:), valz(kj:), info)
!!$          
!!$! !$          write(psb_err_unit,*) i,j,alpha,' Inner loop spaxpby',iv(1:nv),vv(1:nv)
!!$          ! Drop tolerance for new column 
!!$          if (info == psb_success_) call sp_drop(j,nzrmax,sp_thresh,nv,iv,vv,info)
!!$          if (nv > nzralc) then 
!!$            write(0,*) 'Error on size estimate: ',nv,nzralc
!!$            info = -996
!!$            return
!!$          end if
!!$          ! Copy into znew(j) 
!!$          ncz(j)           = nv
!!$          irwz(kj:kj+nv-1) = iv(1:nv)
!!$          valz(kj:kj+nv-1) = vv(1:nv)
!!$! !$          write(psb_err_unit,*) 'At step ',i,j,': ',nv
!!$! !$          write(psb_err_unit,*) '       idxs',irwz(kj:kj+nv-1)
!!$! !$          write(psb_err_unit,*) '       vals',valz(kj:kj+nv-1)
!!$        
!!$        end if
!!$      end do
!!$    end do
!!$    call psb_sp_free(acsc,info) 
!!$    ! Now rebuild Z in CSC format from the vectors
!!$    !
!!$    ki = 1  
!!$    kj = 1
!!$    icpz(1) = 1
!!$    do i=1, n 
!!$      nv = ncz(i) 
!!$      do j=0, nv-1
!!$        irwz(ki+j) = irwz(kj+j)
!!$        valz(ki+j) = valz(kj+j)
!!$      end do
!!$      ki        = ki + nv 
!!$      icpz(i+1) = ki
!!$      kj        = kj + nzralc
!!$    end do
!!$    z%descra  = 'GUN'
!!$    z%fida    = 'CSC'
!!$    z%m       = n
!!$    z%k       = n
!!$    z%infoa   = 0
!!$    call move_alloc(icpz,z%icp)
!!$    call move_alloc(irwz,z%ia)
!!$    call move_alloc(valz,z%val)
!!$    call psb_spcnv(z,info,afmt=psb_csr_afmt_)
!!$    call psb_sp_trim(z,info) 
!!$    
!!$  end subroutine mld_dsparse_orth_3


!!$  subroutine mld_dsparse_orth_4(n,a,p,z,nzrmax,sp_thresh,info)
!!$    !
!!$    ! Attempt to get Benzi-Cullum-Tuma (00) alg. 3.1
!!$    !  with mixed dinamic data. 
!!$    !
!!$    implicit none 
!!$    integer, intent(in)                  :: n
!!$    type(psb_d_csr_sparse_mat), intent(in)    :: a
!!$    type(psb_d_csc_sparse_mat), intent(inout) :: z
!!$    integer, intent(in)                  :: nzrmax
!!$    real(psb_dpk_), intent(in)           :: sp_thresh
!!$    real(psb_dpk_), intent(out)          :: p(:)
!!$    integer, intent(out)                 :: info
!!$    integer, allocatable        :: iv(:), icpz(:), ncz(:),irwz(:), ivz(:)
!!$    real(psb_dpk_), allocatable :: vv(:), valz(:), zw(:), vz(:)
!!$    integer :: i,j,k, kc, kr, err_act, nz, nv, nvz,&
!!$         & ki,kj, nzralc, ipza, nzra
!!$    type(psb_dspmat_type) :: znew, ztmp, acsc
!!$    type(psb_int_heap) :: heap 
!!$    real(psb_dpk_)     :: alpha
!!$
!!$    
!!$    nzralc = nint(2.5*nzrmax)
!!$! !$    write(0,*) 'orth 4 fill in nzr',nzrmax,nzralc
!!$    allocate(iv(n),vv(n),ivz(n),vz(n),icpz(n+1),ncz(n),&
!!$         & irwz(n*nzralc),valz(n*nzralc),stat=info)
!!$    if (info /= psb_success_) then 
!!$      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
!!$      return      
!!$    end if
!!$
!!$    call psb_sp_clip(a,acsc,info,imax=n,jmax=n,&
!!$         & rscale=.false.,cscale=.false.)
!!$    call psb_transp(acsc,fmt='CSC')
!!$  
!!$    !
!!$    ! Init Z to identity.
!!$    !  Each column of Z has room for 2.5 nzrmax, should be plenty
!!$    !  to handle temporaries. 
!!$    !
!!$    j = 1
!!$    do i=1, n
!!$      icpz(i)  = j
!!$      ncz(i)   = 1 
!!$      irwz(j)  = i
!!$      valz(j)  = done
!!$      j        = j + nzralc
!!$    end do
!!$
!!$  
!!$
!!$    do i = 1, n
!!$      ! First, compute v_i=Az_i
!!$      k  = icpz(i)
!!$      nv = 0 
!!$      call psb_d_spmspv(done,acsc,ncz(i),irwz(k:),valz(k:),dzero,nv,iv,vv, info)
!!$      if (info /=0) then 
!!$        write(0,*) 'Internal error from spmspv ',info
!!$        return
!!$      end if
!!$      ! Second, compute the P(J)
!!$! !$      write(psb_err_unit,*) 'Into sparse dot: ',nv,iv(1:nv),vv(1:nv)
!!$      ipza = a%irp(i)
!!$      nzra = a%irp(i+1) - ipza
!!$      k    = icpz(i)
!!$      p(i) = psb_spdot_srtd(nv,iv(1:nv),vv(1:nv),&
!!$           & ncz(i),irwz(k:),valz(k:))
!!$      
!!$      if (p(i) == dzero) then 
!!$        write(psb_err_unit,*) 'Breakdown!! '
!!$        p(i) = 1.d-3
!!$      end if
!!$
!!$      ! Third, right-looking update. 
!!$      ki = icpz(i) 
!!$    
!!$      do j=i+1,n 
!!$        kj   = icpz(j)
!!$        p(j) = psb_spdot_srtd(nv,iv(1:nv),vv(1:nv),&
!!$             & ncz(j),irwz(kj:),valz(kj:))
!!$
!!$        alpha = (-p(j)/p(i))
!!$
!!$        if (abs(alpha) > sp_thresh) then 
!!$
!!$          call psb_nspaxpby(nvz, ivz, vz,&
!!$               & alpha, ncz(i), irwz(ki:), valz(ki:),&
!!$               & done, ncz(j), irwz(kj:), valz(kj:), info)
!!$          
!!$! !$          write(psb_err_unit,*) i,j,alpha,' Inner loop spaxpby',iv(1:nv),vv(1:nv)
!!$          ! Drop tolerance for new column 
!!$          if (info == psb_success_) call sp_drop(j,nzrmax,sp_thresh,nvz,ivz,vz,info)
!!$          if (nvz > nzralc) then 
!!$            write(0,*) 'Error on size estimate: ',nvz,nzralc
!!$            info = -996
!!$            return
!!$          end if
!!$          ! Copy into znew(j) 
!!$          ncz(j)            = nvz
!!$          irwz(kj:kj+nv-1) = ivz(1:nvz)
!!$          valz(kj:kj+nv-1) = vz(1:nvz)
!!$! !$          write(psb_err_unit,*) 'At step ',i,j,': ',nv
!!$! !$          write(psb_err_unit,*) '       idxs',irwz(kj:kj+nv-1)
!!$! !$          write(psb_err_unit,*) '       vals',valz(kj:kj+nv-1)
!!$        
!!$        end if
!!$      end do
!!$    end do
!!$    call psb_sp_free(acsc,info) 
!!$    ! Now rebuild Z in CSC format from the vectors
!!$    !
!!$    ki = 1  
!!$    kj = 1
!!$    icpz(1) = 1
!!$    do i=1, n 
!!$      nv = ncz(i) 
!!$      do j=0, nv-1
!!$        irwz(ki+j) = irwz(kj+j)
!!$        valz(ki+j) = valz(kj+j)
!!$      end do
!!$      ki = ki + nv 
!!$      icpz(i+1) = ki
!!$      kj = kj + nzralc
!!$    end do
!!$    z%descra  = 'GUN'
!!$    z%fida    = 'CSC'
!!$    z%m       = n
!!$    z%k       = n
!!$    z%infoa   = 0
!!$    call move_alloc(icpz,z%icp)
!!$    call move_alloc(irwz,z%ia)
!!$    call move_alloc(valz,z%val)
!!$    call psb_spcnv(z,info,afmt=psb_csr_afmt_)
!!$    call psb_sp_trim(z,info) 
!!$    
!!$  end subroutine mld_dsparse_orth_4

  subroutine psb_d_spmspv(alpha,a,nx,ix,vx,beta,ny,iy,vy, info) 
    use psb_base_mod
    implicit none 
    integer, intent(in)           :: nx, ix(:) 
    real(psb_dpk_), intent(in)    :: alpha, beta, vx(:)
    integer, intent(inout)        :: ny, iy(:) 
    real(psb_dpk_), intent(inout) :: vy(:)
    type(psb_d_csc_sparse_mat), intent(in)  :: a
    integer, intent(out)          :: info 

    integer :: i,j,k,m,n, nv, na, iszy
    integer, allocatable        :: iv(:)
    real(psb_dpk_), allocatable :: vv(:)

    info = 0
    if (beta == -done) then 
      do i=1, ny
        vy(i) = -vy(i) 
      end do
    else if (beta == dzero) then 
      do i=1, ny
        vy(i) = dzero
      end do
    else if (beta /= done) then 
      do i=1, ny
        vy(i) = vy(i) * beta
      end do
    end if
    if (alpha == dzero)  return
    iszy = min(size(iy),size(vy))
    m = a%get_nrows()
    n = a%get_ncols()

    if ((ny > m) .or. (nx > n)) then 
      write(0,*) 'Wrong input spmspv rows: ',m,ny,&
           & ' cols: ',n,nx
      info = -4 
      return 
    end if

    allocate(iv(m), vv(m), stat=info) 
    if (info /= 0) then 
      write(0,*) 'Allocation error in spmspv'
      info = -999
      return
    endif

    do i = 1, nx
      j  = ix(i) 
      ! Access column J of A
      k  = a%icp(j)
      na = a%icp(j+1) - a%icp(j)
      call psb_nspaxpby(nv,iv,vv,&
           & alpha, na, a%ia(k:k+na-1), a%val(k:k+na-1),&
           & done, ny, iy, vy, info)

      if (info /= 0) then 
        write(0,*) 'Internal error in spmspv from nspaxpby'
        info = -998 
        return
      endif
      if (nv > iszy) then 
        write(0,*) 'Error in spmspv: out of memory for output' 
        info = -997
        return
      endif
      ny = nv
      iy(1:ny) = iv(1:ny) 
      vy(1:ny) = vv(1:ny) 
    end do
  end subroutine psb_d_spmspv

end module mld_d_orthbase_mod

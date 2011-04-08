
subroutine mld_dsparse_orthbase(alg,n,a,p,z,nzrmax,sp_thresh,info)
  use psb_base_mod
  use mld_base_prec_type
  use mld_dainv_mod, mld_protect_name => mld_dsparse_orthbase
  integer, intent(in)                  :: alg,n
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_dspmat_type), intent(inout) :: z
  integer, intent(in)                  :: nzrmax
  real(psb_dpk_), intent(in)           :: sp_thresh
  real(psb_dpk_), intent(out)          :: p(:)
  integer, intent(out)                 :: info

!!$  integer :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
!!$       & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz
!!$  integer, allocatable        :: ia(:), ja(:), iz(:),jz(:), icr(:), ikr(:), ljr(:)
!!$  real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
!!$  type(psb_dspmat_type) :: znew, ztmp
!!$  type(psb_int_heap) :: heap 
!!$  real(psb_dpk_)     :: alpha
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

  if (psb_tolower(a%fida) /= 'csr') then 
    write(psb_err_unit,*) 'AFMT wrong'
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='Afmt wrong')
    goto 9999      
  end if

  select case(alg)
  case (mld_ainv_orth1_)
    call mld_dsparse_orth_1(n,a,p,z,nzrmax,sp_thresh,info)
    if (info /= 0) goto 9999
  case (mld_ainv_orth2_) 
    call mld_dsparse_orth_2(n,a,p,z,nzrmax,sp_thresh,info)
    if (info /= 0) goto 9999
  case (mld_ainv_orth3_) 
    call mld_dsparse_orth_3(n,a,p,z,nzrmax,sp_thresh,info)
    if (info /= 0) goto 9999
  case (mld_ainv_orth4_) 
    call mld_dsparse_orth_4(n,a,p,z,nzrmax,sp_thresh,info)
    if (info /= 0) goto 9999
  case default
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='Invalid alg')
    goto 9999      
  end select

  call psb_spcnv(z,info,afmt='CSR')
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

  subroutine mld_dsparse_orth_1(n,a,p,z,nzrmax,sp_thresh,info)
    !
    ! Need to figure out how I am doing this.... 
    !
    !
    implicit none 
    integer, intent(in)                  :: n
    type(psb_dspmat_type), intent(in)    :: a
    type(psb_dspmat_type), intent(inout) :: z
    integer, intent(in)                  :: nzrmax
    real(psb_dpk_), intent(in)           :: sp_thresh
    real(psb_dpk_), intent(out)          :: p(:)
    integer, intent(out)                 :: info
    integer, allocatable        :: ia(:), ja(:), iz(:),jz(:), icr(:), ikr(:), ljr(:)
    real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
    integer :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
         & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz, ipz1, ipz2
    type(psb_int_heap) :: heap 
    real(psb_dpk_)     :: alpha

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
      do j=a%ia2(i),a%ia2(i+1)-1
        k      = a%ia1(j)
        if (k<=n) then 
          ljr(i) = max(ljr(i),k)
          icr(k) = min(icr(k),i)
        end if
      end do
    end do

    ! Init z_1=e_1 and p_1=a_11
    p(1) = dzero
    i   = 1
    nz  = a%ia2(i+1) - a%ia2(i)
    do j=1,nz
      if (a%ia1(j) == 1) then 
        p(1) = a%aspk(j)
        exit
      end if
    end do
    if (abs(p(1)) < d_epstol) &
         & p(1) = 1.d-3 

    call psb_ensure_size(n+1, z%ia2,  info)
    call psb_ensure_size(n*nzrmax, z%ia1,  info)
    call psb_ensure_size(n*nzrmax, z%aspk, info)
    ! 
    !
    z%descra  = 'GUN'
    z%fida    = 'CSC'
    z%m       = n
    z%k       = n
    z%ia2(1)  = 1
    z%ia2(2)  = 2
    z%ia1(1)  = 1
    z%aspk(1) = done
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
        if (.false.) then 
          call psb_sp_getrow(j,a,nzra,ia,ja,val,info)
          call rwclip(nzra,ia,ja,val,1,n,1,n)

          p(i) = psb_spge_dot(nzra,ja,val,zw)
        else
          ip1 = a%ia2(j)
          ip2 = a%ia2(j+1) - 1
          do 
            if (ip2 < ip1 ) exit
            if (a%ia1(ip2) <= n) exit
            ip2 = ip2 -1 
          end do
          nzra = max(0,ip2 - ip1 + 1) 
          p(i) = psb_spge_dot(nzra,a%ia1(ip1:ip2),a%aspk(ip1:ip2),zw)
          ! !$          write(psb_err_unit,*) j,i,p(i)
        end if

        ipz1 = z%ia2(j) 
        ipz2 = z%ia2(j+1) 
        nzrz = ipz2-ipz1
        alpha = (-p(i)/p(j))
        if (abs(alpha) > sp_thresh) then 
          ! !$          call psb_aspxpby(alpha,nzrz,&
          ! !$               & z%ia1(ipz1:ipz2-1),z%aspk(ipz1:ipz2-1), &
          ! !$               & done,zw,info)
          do k=ipz1, ipz2-1
            kr     = z%ia1(k)
            zw(kr) = zw(kr) + alpha*z%aspk(k)
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
      call psb_sp_getrow(i,a,nzra,ia,ja,val,info)
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
      call psb_ensure_size(nzz+nzrz, z%ia1,   info)
      call psb_ensure_size(nzz+nzrz, z%aspk, info)
      ipz1 = z%ia2(i)
      do j=1, nzrz
        z%ia1(ipz1  + j -1) = iz(j)
        z%aspk(ipz1 + j -1) = valz(j)
      end do
      z%ia2(i+1) = ipz1 + nzrz
      nzz        = nzz + nzrz
    end do

  end subroutine mld_dsparse_orth_1

  subroutine mld_dsparse_orth_2(n,a,p,z,nzrmax,sp_thresh,info)
    !
    ! Benzi-Tuma (98): alg biconjugation (section 4).
    !  Almost idiotic implementation 
    !

    implicit none 
    integer, intent(in)                  :: n
    type(psb_dspmat_type), intent(in)    :: a
    type(psb_dspmat_type), intent(inout) :: z
    integer, intent(in)                  :: nzrmax
    real(psb_dpk_), intent(in)           :: sp_thresh
    real(psb_dpk_), intent(out)          :: p(:)
    integer, intent(out)                 :: info
    integer, allocatable        :: ia(:), ja(:), iz(:),jz(:)
    real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
    integer :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
         & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn,ifnz
    type(psb_dspmat_type) :: znew, ztmp
    type(psb_int_heap) :: heap 
    real(psb_dpk_)     :: alpha

    allocate(zw(n),iz(n),jz(n),valz(n),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      return      
    end if

    ! Init z_1=e_1 and p_1=a_11
    ! !$    p(1) = dzero
    ! !$    i   = 1
    ! !$    nz  = a%ia2(i+1) - a%ia2(i)
    ! !$    do j=1,nz
    ! !$      if (a%ia1(j) == 1) then 
    ! !$        p(1) = a%aspk(j)
    ! !$        exit
    ! !$      end if
    ! !$    end do
    ! !$    if (abs(p(1)) < d_epstol) &
    ! !$         & p(1) = 1.d-3 

    call psb_ensure_size(n+1, z%ia2,  info)
    call psb_ensure_size(n*nzrmax, z%ia1,  info)
    call psb_ensure_size(n*nzrmax, z%aspk, info)
    ! 
    !
    z%descra  = 'GUN'
    z%fida    = 'CSC'
    z%m       = n
    z%k       = n
    z%ia2(1)  = 1
    ! Init Z to identity 
    do i=1, n
      j           = z%ia2(i) 
      z%ia1(j)    = i
      z%aspk(j)   = done
      z%ia2(i+1)  = j+1
    end do

    call psb_sp_clone(z,znew,info)   
    if (info /= psb_success_) then 
      info = psb_err_internal_error_
      call psb_errpush(psb_err_internal_error_,name,a_err='clone')
      return      
    end if

    nzz = 0 
    do i = 1, n

      ipza = a%ia2(i)
      nzra = a%ia2(i+1) - ipza
!!$      write(psb_err_unit,*) 'Into sparse dot: ',nzra,a%ia1(ipza:ipza+nzra-1),a%aspk(ipza:ipza+nzra-1)
      do j=i, n
        ipzz = z%ia2(j)
        nzrz = z%ia2(j+1) - ipzz
        p(j) = psb_spdot_srtd(nzra,a%ia1(ipza:ipza+nzra-1),a%aspk(ipza:ipza+nzra-1),&
             & nzrz,z%ia1(ipzz:ipzz+nzrz-1),z%aspk(ipzz:ipzz+nzrz-1))
!!$        write(psb_err_unit,*) '     ..........: ',nzrz,z%ia1(ipzz:ipzz+nzrz-1),z%aspk(ipzz:ipzz+nzrz-1)
!!$        write(psb_err_unit,*) i,j,p(i),p(j)
      end do
      if (p(i) == dzero) then 
        write(psb_err_unit,*) 'Breakdown!! '
        p(i) = 1.d-3
      end if

      ipzi = z%ia2(i)
      nzzi = z%ia2(i+1) - ipzi
      do j=i+1,n 
        ! New kernel psb_aspxpbspy 
        ipzj = z%ia2(j)
        nzzj = z%ia2(j+1) - ipzj
        ipzn = znew%ia2(j)
        alpha = (-p(j)/p(i))
        if (abs(alpha) > sp_thresh) then 

          call psb_nspaxpby(nzrz,iz,valz,&
               & alpha, nzzi, z%ia1(ipzi:ipzi+nzzi-1), z%aspk(ipzi:ipzi+nzzi-1),&
               & done, nzzj, z%ia1(ipzj:ipzj+nzzj-1), z%aspk(ipzj:ipzj+nzzj-1),&
               & info)

!!$          write(psb_err_unit,*) i,j,alpha,' Inner loop spaxpby',iz(1:nzrz),valz(1:nzrz)
          ! Drop tolerance for new column 
          if (info == psb_success_) call sp_drop(j,nzrmax,sp_thresh,nzrz,iz,valz,info)
          ! Copy into znew(j) 
          if (info == psb_success_) call psb_ensure_size(ipzn+nzrz-1,znew%ia1,info)
          if (info == psb_success_) call psb_ensure_size(ipzn+nzrz-1,znew%aspk,info)
          if (info /= psb_success_) then 
            info = psb_err_internal_error_
            call psb_errpush(psb_err_internal_error_,name,a_err='Inner loop ')
            return      
          end if
          znew%ia1(ipzn:ipzn+nzrz-1)  = iz(1:nzrz)
          znew%aspk(ipzn:ipzn+nzrz-1) = valz(1:nzrz)
          znew%ia2(j+1) = ipzn + nzrz
        else 
          nzrz = nzzj 
          call psb_ensure_size(ipzn+nzrz-1,znew%ia1,info)
          if (info == psb_success_) call psb_ensure_size(ipzn+nzrz-1,znew%aspk,info)
          znew%ia1(ipzn:ipzn+nzrz-1)  = z%ia1(ipzj:ipzj+nzzj-1)
          znew%aspk(ipzn:ipzn+nzrz-1) = z%aspk(ipzj:ipzj+nzzj-1)
          znew%ia2(j+1) = ipzn + nzrz
        end if
!!$        write(psb_err_unit,*) 'At step ',i,j,': ',nzrz
!!$        write(psb_err_unit,*) '       idxs',znew%ia1(ipzn:ipzn+nzrz-1)
!!$        write(psb_err_unit,*) '       vals',znew%aspk(ipzn:ipzn+nzrz-1)
      end do
      if (i < n) then 
        ! !$        write(psb_err_unit,*)  i,'Copy step ', z%ia2(i+1:i+2),':',znew%ia2(i+1:i+2)
        z%ia2(i+1:i+2) = znew%ia2(i+1:i+2)
        ipzj = z%ia2(i+1)
        nzzj = z%ia2(i+2) - ipzj
        z%ia1(ipzj:ipzj+nzzj-1)  = znew%ia1(ipzj:ipzj+nzzj-1)
        z%aspk(ipzj:ipzj+nzzj-1) = znew%aspk(ipzj:ipzj+nzzj-1) 
        ! !$        write(psb_err_unit,*)  i,'Copy step ', z%ia1(ipzj:ipzj+nzzj-1), z%aspk(ipzj:ipzj+nzzj-1) 
        call psb_move_alloc(znew,ztmp,info)
        if (info == psb_success_) call psb_move_alloc(z,znew,info)
        if (info == psb_success_) call psb_move_alloc(ztmp,z,info)
      end if
      if (info /= psb_success_) then 
        info = psb_err_internal_error_
        call psb_errpush(psb_err_internal_error_,name,a_err='Outer loop ')
        return      
      end if
    end do
  end subroutine mld_dsparse_orth_2

  subroutine mld_dsparse_orth_3(n,a,p,z,nzrmax,sp_thresh,info)
    !
    ! Benzi-Tuma (98): alg biconjugation (section 4).
    !  Better implementation, with mixed dinamic data. 
    !
    implicit none 
    integer, intent(in)                  :: n
    type(psb_dspmat_type), intent(in)    :: a
    type(psb_dspmat_type), intent(inout) :: z
    integer, intent(in)                  :: nzrmax
    real(psb_dpk_), intent(in)           :: sp_thresh
    real(psb_dpk_), intent(out)          :: p(:)
    integer, intent(out)                 :: info
    integer, allocatable        :: iv(:), icpz(:), ncz(:),irwz(:)
    real(psb_dpk_), allocatable :: vv(:), valz(:), zw(:)
    integer :: i,j,k, kc, kr, err_act, nz, nv,&
         & ki,kj, nzralc, ipza, nzra
    type(psb_dspmat_type) :: znew, ztmp, acsc
    type(psb_int_heap) :: heap 
    real(psb_dpk_)     :: alpha

    
    nzralc = nint(2.5*nzrmax)
    write(0,*) 'fill in nzr',nzrmax,nzralc
    allocate(iv(n),vv(n),icpz(n+1),ncz(n),&
         & irwz(n*nzralc),valz(n*nzralc),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      return      
    end if

!!$    call psb_sp_clip(a,acsc,info,imax=n,jmax=n,&
!!$         & rscale=.false.,cscale=.false.)
!!$    call psb_transp(acsc,fmt='CSC')
!!$  

    !
    ! Init Z to identity.
    !  Each column of Z has room for 2.5 nzrmax, should be plenty
    !  to handle temporaries. 
    !
    j = 1
    do i=1, n
      icpz(i)  = j
      ncz(i)   = 1 
      irwz(j)  = i
      valz(j)  = done
      j        = j + nzralc
    end do

  

    do i = 1, n
!!$      ! First, compute v_i=Az_i
      k  = icpz(i)
      nv = 0 
!!$      call psb_d_spmspv(done,acsc,ncz(i),irwz(k:),valz(k:),dzero,nv,iv,vv, info)
!!$      if (info /=0) then 
!!$        write(0,*) 'Internal error from spmspv ',info
!!$        return
!!$      end if
      ! Second, compute the P(J)
!!$      write(psb_err_unit,*) 'Into sparse dot: ',nv,iv(1:nv),vv(1:nv)
      ipza = a%ia2(i)
      nzra = a%ia2(i+1) - ipza
      k    = icpz(i)
      p(i) = psb_spdot_srtd(nzra,a%ia1(ipza:ipza+nzra-1),a%aspk(ipza:ipza+nzra-1),&
           & ncz(i),irwz(k:),valz(k:))
      
      if (p(i) == dzero) then 
        write(psb_err_unit,*) 'Breakdown!! '
        p(i) = 1.d-3
      end if

      ! Third, right-looking update. 
      ki = icpz(i) 
    
      do j=i+1,n 
        kj   = icpz(j)
        p(j) = psb_spdot_srtd(nzra,a%ia1(ipza:ipza+nzra-1),a%aspk(ipza:ipza+nzra-1),&
             & ncz(j),irwz(kj:),valz(kj:))

        alpha = (-p(j)/p(i))

        if (abs(alpha) > sp_thresh) then 

          call psb_nspaxpby(nv, iv, vv,&
               & alpha, ncz(i), irwz(ki:), valz(ki:),&
               & done, ncz(j), irwz(kj:), valz(kj:), info)
          
!!$          write(psb_err_unit,*) i,j,alpha,' Inner loop spaxpby',iv(1:nv),vv(1:nv)
          ! Drop tolerance for new column 
          if (info == psb_success_) call sp_drop(j,nzrmax,sp_thresh,nv,iv,vv,info)
          if (nv > nzralc) then 
            write(0,*) 'Error on size estimate: ',nv,nzralc
            info = -996
            return
          end if
          ! Copy into znew(j) 
          ncz(j)           = nv
          irwz(kj:kj+nv-1) = iv(1:nv)
          valz(kj:kj+nv-1) = vv(1:nv)
!!$          write(psb_err_unit,*) 'At step ',i,j,': ',nv
!!$          write(psb_err_unit,*) '       idxs',irwz(kj:kj+nv-1)
!!$          write(psb_err_unit,*) '       vals',valz(kj:kj+nv-1)
        
        end if
      end do
    end do
    call psb_sp_free(acsc,info) 
    ! Now rebuild Z in CSC format from the vectors
    !
    ki = 1  
    kj = 1
    icpz(1) = 1
    do i=1, n 
      nv = ncz(i) 
      do j=0, nv-1
        irwz(ki+j) = irwz(kj+j)
        valz(ki+j) = valz(kj+j)
      end do
      ki        = ki + nv 
      icpz(i+1) = ki
      kj        = kj + nzralc
    end do
    z%descra  = 'GUN'
    z%fida    = 'CSC'
    z%m       = n
    z%k       = n
    z%infoa   = 0
    call move_alloc(icpz,z%ia2)
    call move_alloc(irwz,z%ia1)
    call move_alloc(valz,z%aspk)
    call psb_spcnv(z,info,afmt=psb_csr_afmt_)
    call psb_sp_trim(z,info) 
    
  end subroutine mld_dsparse_orth_3


  subroutine mld_dsparse_orth_4(n,a,p,z,nzrmax,sp_thresh,info)
    !
    ! Attempt to get Benzi-Cullum-Tuma (00) alg. 3.1
    !  with mixed dinamic data. 
    !
    implicit none 
    integer, intent(in)                  :: n
    type(psb_dspmat_type), intent(in)    :: a
    type(psb_dspmat_type), intent(inout) :: z
    integer, intent(in)                  :: nzrmax
    real(psb_dpk_), intent(in)           :: sp_thresh
    real(psb_dpk_), intent(out)          :: p(:)
    integer, intent(out)                 :: info
    integer, allocatable        :: iv(:), icpz(:), ncz(:),irwz(:), ivz(:)
    real(psb_dpk_), allocatable :: vv(:), valz(:), zw(:), vz(:)
    integer :: i,j,k, kc, kr, err_act, nz, nv, nvz,&
         & ki,kj, nzralc, ipza, nzra
    type(psb_dspmat_type) :: znew, ztmp, acsc
    type(psb_int_heap) :: heap 
    real(psb_dpk_)     :: alpha

    
    nzralc = nint(2.5*nzrmax)
!!$    write(0,*) 'orth 4 fill in nzr',nzrmax,nzralc
    allocate(iv(n),vv(n),ivz(n),vz(n),icpz(n+1),ncz(n),&
         & irwz(n*nzralc),valz(n*nzralc),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      return      
    end if

    call psb_sp_clip(a,acsc,info,imax=n,jmax=n,&
         & rscale=.false.,cscale=.false.)
    call psb_transp(acsc,fmt='CSC')
  
    !
    ! Init Z to identity.
    !  Each column of Z has room for 2.5 nzrmax, should be plenty
    !  to handle temporaries. 
    !
    j = 1
    do i=1, n
      icpz(i)  = j
      ncz(i)   = 1 
      irwz(j)  = i
      valz(j)  = done
      j        = j + nzralc
    end do

  

    do i = 1, n
      ! First, compute v_i=Az_i
      k  = icpz(i)
      nv = 0 
      call psb_d_spmspv(done,acsc,ncz(i),irwz(k:),valz(k:),dzero,nv,iv,vv, info)
      if (info /=0) then 
        write(0,*) 'Internal error from spmspv ',info
        return
      end if
      ! Second, compute the P(J)
!!$      write(psb_err_unit,*) 'Into sparse dot: ',nv,iv(1:nv),vv(1:nv)
      ipza = a%ia2(i)
      nzra = a%ia2(i+1) - ipza
      k    = icpz(i)
      p(i) = psb_spdot_srtd(nv,iv(1:nv),vv(1:nv),&
           & ncz(i),irwz(k:),valz(k:))
      
      if (p(i) == dzero) then 
        write(psb_err_unit,*) 'Breakdown!! '
        p(i) = 1.d-3
      end if

      ! Third, right-looking update. 
      ki = icpz(i) 
    
      do j=i+1,n 
        kj   = icpz(j)
        p(j) = psb_spdot_srtd(nv,iv(1:nv),vv(1:nv),&
             & ncz(j),irwz(kj:),valz(kj:))

        alpha = (-p(j)/p(i))

        if (abs(alpha) > sp_thresh) then 

          call psb_nspaxpby(nvz, ivz, vz,&
               & alpha, ncz(i), irwz(ki:), valz(ki:),&
               & done, ncz(j), irwz(kj:), valz(kj:), info)
          
!!$          write(psb_err_unit,*) i,j,alpha,' Inner loop spaxpby',iv(1:nv),vv(1:nv)
          ! Drop tolerance for new column 
          if (info == psb_success_) call sp_drop(j,nzrmax,sp_thresh,nvz,ivz,vz,info)
          if (nvz > nzralc) then 
            write(0,*) 'Error on size estimate: ',nvz,nzralc
            info = -996
            return
          end if
          ! Copy into znew(j) 
          ncz(j)            = nvz
          irwz(kj:kj+nv-1) = ivz(1:nvz)
          valz(kj:kj+nv-1) = vz(1:nvz)
!!$          write(psb_err_unit,*) 'At step ',i,j,': ',nv
!!$          write(psb_err_unit,*) '       idxs',irwz(kj:kj+nv-1)
!!$          write(psb_err_unit,*) '       vals',valz(kj:kj+nv-1)
        
        end if
      end do
    end do
    call psb_sp_free(acsc,info) 
    ! Now rebuild Z in CSC format from the vectors
    !
    ki = 1  
    kj = 1
    icpz(1) = 1
    do i=1, n 
      nv = ncz(i) 
      do j=0, nv-1
        irwz(ki+j) = irwz(kj+j)
        valz(ki+j) = valz(kj+j)
      end do
      ki = ki + nv 
      icpz(i+1) = ki
      kj = kj + nzralc
    end do
    z%descra  = 'GUN'
    z%fida    = 'CSC'
    z%m       = n
    z%k       = n
    z%infoa   = 0
    call move_alloc(icpz,z%ia2)
    call move_alloc(irwz,z%ia1)
    call move_alloc(valz,z%aspk)
    call psb_spcnv(z,info,afmt=psb_csr_afmt_)
    call psb_sp_trim(z,info) 
    
  end subroutine mld_dsparse_orth_4

  subroutine psb_d_spmspv(alpha,a,nx,ix,vx,beta,ny,iy,vy, info) 
    use psb_base_mod
    implicit none 
    integer, intent(in)           :: nx, ix(:) 
    real(psb_dpk_), intent(in)    :: alpha, beta, vx(:)
    integer, intent(inout)        :: ny, iy(:) 
    real(psb_dpk_), intent(inout) :: vy(:)
    type(psb_dspmat_type), intent(in)  :: a
    integer, intent(out)          :: info 
    
    integer :: i,j,k,m,n, nv, na, iszy
    integer, allocatable        :: iv(:)
    real(psb_dpk_), allocatable :: vv(:)

    info = 0
    if (psb_get_fmt(a)  /= psb_csc_afmt_) then 
      write(0,*) 'Only CSC so far in spmspv'
      info = -2 
      return
    end if
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
    m = psb_sp_get_nrows(a)
    n = psb_sp_get_ncols(a)
    
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
      k  = a%ia2(j)
      na = a%ia2(j+1) - a%ia2(j)
      call psb_nspaxpby(nv,iv,vv,&
               & alpha, na, a%ia1(k:k+na-1), a%aspk(k:k+na-1),&
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

end subroutine mld_dsparse_orthbase

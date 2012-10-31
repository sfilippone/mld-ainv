subroutine mld_dsparse_ainvk(n,a,z,fill_in,sp_thresh,info,inlevs)
  
  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dsparse_ainvk

  integer, intent(in)                  :: n
  type(psb_dspmat_type), intent(in)    :: a
  type(psb_dspmat_type), intent(inout) :: z
  integer, intent(in)                  :: fill_in
  real(psb_dpk_), intent(in)           :: sp_thresh
  integer, intent(out)                 :: info
  integer, intent(in), optional        :: inlevs(:)

  integer :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2 
  integer, allocatable        :: ia(:), ja(:), iz(:), jz(:) 
  real(psb_dpk_), allocatable :: zw(:), val(:), valz(:)
  integer, allocatable        :: uplevs(:), rowlevs(:), idxs(:)
  real(psb_dpk_), allocatable :: row(:)
  type(psb_d_coo_sparse_mat)  :: trw
  type(psb_d_csr_sparse_mat)  :: acsr, zcsr
  integer                  :: ktrw, nidx
  type(psb_int_heap)       :: heap

  real(psb_dpk_)     :: alpha
  character(len=20)  :: name='mld_sp_ainvk'


  if(psb_get_errstatus() /= psb_success_) return 
  info = psb_success_
  call psb_erractionsave(err_act)

  if (.not.(a%is_triangle().and.a%is_unit().and.a%is_upper())) then 
    write(psb_err_unit,*) 'Wrong A ' 
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='wrong A')
    goto 9999      
  end if
  call a%cp_to(acsr)
  call trw%allocate(0,0,1)
  if (info == psb_success_) allocate(zw(n),iz(n),valz(n),&
       & row(n),rowlevs(n),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  allocate(uplevs(acsr%get_nzeros()),stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='Allocate')
    goto 9999
  end if
  uplevs(:)  = 0
  row(:)     = dzero
  rowlevs(:) = -(n+1)

  call zcsr%allocate(n,n,n*fill_in)
  call zcsr%set_triangle()
  call zcsr%set_unit(.false.)
  call zcsr%set_upper()
  call psb_ensure_size(n+1, idxs,  info)


  ! 
  !
  zcsr%irp(1)  = 1
  nzz          = 0

  l2 = 0
  outer: do i = 1, n-1
    ! ZW = e_i
    call mld_invk_copyin(i,n,acsr,1,n,row,rowlevs,heap,ktrw,trw,info,&
         & sign=-done,inlevs=inlevs)
    row(i)     = done
    rowlevs(i) = 0
!!$      call psb_insert_heap(i,heap,info) ! No we don't want to put I in. 

    ! Update loop
    call mld_invk(fill_in,i,row,rowlevs,heap,&
         & acsr%ja,acsr%irp,acsr%val,uplevs,nidx,idxs,info)

    call mld_invk_copyout(fill_in,i,n,row,rowlevs,nidx,idxs,&
         & l2,zcsr%ja,zcsr%irp,zcsr%val,info)

    nzz = l2
  end do outer
  if (info /= psb_success_) then 
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='mainloop')
    goto 9999
  end if
  ipz1 = nzz+1
  call psb_ensure_size(ipz1,zcsr%val,info)
  call psb_ensure_size(ipz1,zcsr%ja,info)
  zcsr%val(ipz1) = done
  zcsr%ja(ipz1)  = n
  zcsr%irp(n+1)  = ipz1+1 
  call zcsr%set_sorted()  
  call z%mv_from(zcsr)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine mld_dsparse_ainvk
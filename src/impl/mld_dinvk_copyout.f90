subroutine mld_dinvk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
  
     &  l2,uia1,uia2,uaspk,info)

  use psb_base_mod
  use mld_d_invk_solver, mld_protect_name => mld_dinvk_copyout

  implicit none 

  ! Arguments
  integer, intent(in)                        :: fill_in, i, m, nidx
  integer, intent(inout)                     :: l2, info
  integer, intent(inout)                     :: rowlevs(:), idxs(:)
  integer, allocatable, intent(inout)        :: uia1(:), uia2(:)
  real(psb_dpk_), allocatable, intent(inout) :: uaspk(:)
  real(psb_dpk_), intent(inout)              :: row(:)

  ! Local variables
  integer               :: j,isz,err_act,int_err(5),idxp
  character(len=20), parameter  :: name='mld_diluk_factint'
  character(len=20)             :: ch_err

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)


  do idxp=1,nidx

    j = idxs(idxp)


    if (j>=i) then 
      !
      ! Copy the upper part of the row
      ! 
      if (rowlevs(j) <= fill_in) then 
        l2     = l2 + 1 
        if (size(uaspk) < l2) then 
          ! 
          ! Figure out a good reallocation size!
          !
          isz  = max(int(1.2*l2),l2+100)
          call psb_realloc(isz,uaspk,info) 
          if (info == psb_success_) call psb_realloc(isz,uia1,info) 
          if (info /= psb_success_) then 
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='Allocate')
            goto 9999
          end if
        end if
        uia1(l2)   = j
        uaspk(l2)  = row(j)
      end if
      !
      ! Re-initialize row(j) and rowlevs(j)
      !
      row(j)     = dzero
      rowlevs(j) = -(m+1)
    end if
  end do

  uia2(i+1) = l2 + 1

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine mld_dinvk_copyout

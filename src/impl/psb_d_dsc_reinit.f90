subroutine psb_d_dsc_reinit(a,clear)
  
  use psb_error_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_reinit
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  Integer :: err_act, info,i
  character(len=20)  :: name='reinit'
  logical  :: clear_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_


  if (present(clear)) then 
    clear_ = clear
  else
    clear_ = .true.
  end if

  if (a%is_bld() .or. a%is_upd()) then 
    ! do nothing
    return
  else if (a%is_asb()) then 
    if (clear_) then 
      if (allocated(a%cols)) then 
        do i=1,size(a%cols)
          a%cols(i)%val(:) = dzero
        end do
      end if
    end if
    call a%set_upd()
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_dsc_reinit

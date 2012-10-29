subroutine psb_d_mv_dsc_from_coo(a,b,info) 
  
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_mv_dsc_from_coo
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc, icl
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ! Too tired to figure out a smart mv right now. 
  call a%cp_from_coo(b,info)
  call b%free()
!!$
!!$  call b%fix(info, idir=1)
!!$  if (info /= psb_success_) return
!!$
!!$  nr  = b%get_nrows()
!!$  nc  = b%get_ncols()
!!$  nza = b%get_nzeros()
!!$  
!!$  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
!!$
!!$  ! Dirty trick: call move_alloc to have the new data allocated just once.
!!$  call move_alloc(b%ja,itemp)
!!$  call move_alloc(b%ia,a%ia)
!!$  call move_alloc(b%val,a%val)
!!$  call psb_realloc(max(nr+1,nc+1),a%icp,info)
!!$  call b%free()
!!$
!!$  if (nza <= 0) then 
!!$    a%icp(:) = 1
!!$  else
!!$    a%icp(1) = 1
!!$    if (nc < itemp(nza)) then 
!!$      write(debug_unit,*) trim(name),': CLSHR=.false. : ',&
!!$           &nc,itemp(nza),' Expect trouble!'
!!$      info = 12
!!$    end if
!!$
!!$    j = 1 
!!$    i = 1
!!$    icl = itemp(j) 
!!$
!!$    outer: do 
!!$      inner: do 
!!$        if (i >= icl) exit inner
!!$        if (i > nc) then 
!!$          write(debug_unit,*) trim(name),&
!!$               & 'Strange situation: i>nr ',i,nc,j,nza,icl
!!$          exit outer
!!$        end if
!!$        a%icp(i+1) = a%icp(i) 
!!$        i = i + 1
!!$      end do inner
!!$      j = j + 1
!!$      if (j > nza) exit
!!$      if (itemp(j) /= icl) then 
!!$        a%icp(i+1) = j
!!$        icl = itemp(j) 
!!$        i = i + 1
!!$      endif
!!$      if (i > nc) exit
!!$    enddo outer
!!$    !
!!$    ! Cleanup empty rows at the end
!!$    !
!!$    if (j /= (nza+1)) then 
!!$      write(debug_unit,*) trim(name),': Problem from loop :',j,nza
!!$      info = 13
!!$    endif
!!$    do 
!!$      if (i > nc) exit
!!$      a%icp(i+1) = j
!!$      i = i + 1
!!$    end do
!!$
!!$  endif
!!$

end subroutine psb_d_mv_dsc_from_coo

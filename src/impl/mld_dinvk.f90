subroutine mld_dinvk_inv(fill_in,i,row,rowlevs,heap,uia1,uia2,uaspk,uplevs,&
     & nidx,idxs,info)
  

  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dinvk_inv

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

  ! Local variables
  integer             :: k,j,lrwk,jj,lastk, iret
  real(psb_dpk_)      :: rwk


  info = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i

  !
  ! Do while there are indices to be processed
  !
  do
    ! Beware: (iret < 0) means that the heap is empty, not an error.
    call psb_heap_get_first(k,heap,iret) 
    if (iret < 0) then 
!!$        write(psb_err_unit,*) 'IINVK: ',i,' returning at ',lastk
      return
    end if

    ! 
    ! Just in case an index has been put on the heap more than once.
    !
    if (k == lastk) cycle

    lastk = k 
    nidx = nidx + 1
    if (nidx>size(idxs)) then 
      call psb_realloc(nidx+psb_heap_resize,idxs,info)
      if (info /= psb_success_) return
    end if
    idxs(nidx) = k

    if ((row(k) /= dzero).and.(rowlevs(k) <= fill_in)) then 
      !
      ! Note: since U is scaled while copying it out (see iluk_copyout),
      ! we can use rwk in the update below
      ! 
      rwk    = row(k)
      lrwk   = rowlevs(k)

      do jj=uia2(k),uia2(k+1)-1
        j = uia1(jj)
        if (j<=k) then 
          info = -i
          return
        endif
        !
        ! Insert the index into the heap for further processing.
        ! The fill levels are initialized to a negative value. If we find
        ! one, it means that it is an as yet untouched index, so we need
        ! to insert it; otherwise it is already on the heap, there is no
        ! need to insert it more than once. 
        !
        if (rowlevs(j)<0) then 
          call psb_insert_heap(j,heap,info)
          if (info /= psb_success_) return
          rowlevs(j) = abs(rowlevs(j))
        end if
        !
        ! Update row(j) and the corresponding fill level
        !
        row(j)     = row(j) - rwk * uaspk(jj)
        rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
      end do

    end if
  end do

end subroutine mld_dinvk_inv

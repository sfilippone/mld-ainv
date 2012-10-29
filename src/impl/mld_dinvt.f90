subroutine mld_dinvt(thres,i,nrmi,row,heap,irwt,ja,irp,val,nidx,idxs,info)
  

  use psb_base_mod
  use mld_d_ainv_bld_mod, mld_protect_name => mld_dinvt

  implicit none 

  ! Arguments
  type(psb_int_heap), intent(inout)   :: heap 
  integer, intent(in)                 :: i
  integer, intent(inout)              :: nidx,info
  integer, intent(inout)              :: irwt(:) 
  real(psb_dpk_), intent(in)          :: thres,nrmi
  integer, allocatable, intent(inout) :: idxs(:)
  integer, intent(in)                 :: ja(:),irp(:)
  real(psb_dpk_), intent(in)          :: val(:)
  real(psb_dpk_), intent(inout)       :: row(:)

  ! Local Variables
  integer             :: k,j,jj,lastk,iret
  real(psb_dpk_)      :: rwk, alpha

  info  = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i
  irwt(i) = 1 
!!$  write(0,*) 'Drop Threshold ',thres*nrmi
  !
  ! Do while there are indices to be processed
  !
  do

    call psb_heap_get_first(k,heap,iret) 
    if (iret < 0)  exit

    ! 
    ! An index may have been put on the heap more than once.
    ! Should not happen, but just in case. 
    !
    if (k == lastk)  cycle
    lastk = k 

    !
    ! Dropping rule based on the threshold: compare the absolute
    ! value of each updated entry of row with thres * 2-norm of row.
    !
    rwk    = row(k)

    if (abs(rwk) < thres*nrmi) then
      ! 
      ! Drop the entry.
      !
      row(k)  = dzero
      irwt(k) = 0
      cycle
    else
      !
      ! Note: since U is scaled while copying it out (see ilut_copyout),
      ! we can use rwk in the update below.
      !           
      do jj=irp(k),irp(k+1)-1
        j = ja(jj)
        if (j<=k) then 
          info = -i 
          return
        endif
        !
        ! Update row(j) and, if it is not to be discarded, insert
        ! its index into the heap for further processing.
        !
        row(j)     = row(j) - rwk * val(jj)
        if (irwt(j) == 0) then 
          if (abs(row(j)) < thres*nrmi) then
            ! 
            ! Drop the entry.
            !
            row(j)  = dzero
          else
            !
            ! Do the insertion.
            !
            call psb_insert_heap(j,heap,info)
            if (info /= psb_success_) return
            irwt(j) = 1
          end if
        end if
      end do
    end if

    !
    ! If we get here it is an index we need to keep on copyout.
    !

    nidx       = nidx + 1
    call psb_ensure_size(nidx,idxs,info,addsz=psb_heap_resize)      
    if (info /= psb_success_) return
    idxs(nidx) = k
    irwt(k)    = 0
  end do

  irwt(i) = 0
end subroutine mld_dinvt

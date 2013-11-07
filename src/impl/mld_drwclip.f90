
subroutine mld_d_rwclip(nz,ia,ja,val,imin,imax,jmin,jmax)
  use psb_base_mod

  implicit none 
  integer, intent(inout) :: nz
  integer, intent(inout) :: ia(*), ja(*)
  real(psb_dpk_), intent(inout) :: val(*)
  integer, intent(in)    :: imin,imax,jmin,jmax

  integer :: i,j 

  j = 0
  do i=1, nz
    if ((imin <= ia(i)).and.&
         & (ia(i) <= imax).and.&
         & (jmin <= ja(i)).and.&
         & (ja(i) <= jmax) ) then 
      j = j + 1 
      ia(j) = ia(i) 
      ja(j) = ja(i)
      val(j) = val(i)
    end if
  end do
  nz = j 
end subroutine mld_d_rwclip


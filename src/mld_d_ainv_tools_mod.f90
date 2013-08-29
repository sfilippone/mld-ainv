module mld_d_ainv_tools_mod

  interface sp_drop
    subroutine mld_d_sp_drop(idiag,nzrmax,sp_thresh,nz,iz,valz,info)
      use psb_base_mod, only : psb_dpk_, psb_ipk_
      implicit none 
      real(psb_dpk_), intent(in)    :: sp_thresh
      integer, intent(in)           :: idiag, nzrmax
      integer, intent(inout)        :: nz
      integer, intent(inout)        :: iz(:)
      real(psb_dpk_), intent(inout) :: valz(:)
      integer, intent(out)          :: info
    end subroutine mld_d_sp_drop
  end interface

  interface rwclip
    subroutine mld_d_rwclip(nz,ia,ja,val,imin,imax,jmin,jmax)
      use psb_base_mod, only : psb_dpk_, psb_ipk_
      
      implicit none 
      integer, intent(inout) :: nz
      integer, intent(inout) :: ia(*), ja(*)
      real(psb_dpk_), intent(inout) :: val(*)
      integer, intent(in)    :: imin,imax,jmin,jmax
    end subroutine mld_d_rwclip
  end interface
  
  interface sparsify
    subroutine mld_d_sparsify(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,info, &
         & istart,iheap,ikr)
      use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_int_heap
      implicit none 
      
      real(psb_dpk_), intent(in)  :: sp_thresh
      integer, intent(in)         :: idiag, n, nzrmax
      real(psb_dpk_), intent(inout)  :: zw(:)
      integer, intent(out)        :: nz
      integer, intent(out)        :: iz(:)
      real(psb_dpk_), intent(out) :: valz(:)
      integer, intent(out)        :: info
      integer, intent(in), optional :: istart
      type(psb_int_heap), optional :: iheap
      integer, optional            :: ikr(:)
    end subroutine mld_d_sparsify
  end interface

end module mld_d_ainv_tools_mod

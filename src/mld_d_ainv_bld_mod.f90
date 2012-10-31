module mld_d_ainv_bld_mod
  
  use mld_base_ainv_mod
  
  use mld_d_ilu_fact_mod
  interface mld_ainv_invt_bld
    module procedure mld_d_ainv_invt_bld
  end interface
  

  interface mld_invt_copyin
    subroutine mld_dinvt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
         & irwt,ktrw,trw,info,sign)
      use psb_base_mod, only : psb_d_csr_sparse_mat, psb_d_coo_sparse_mat,&
           & psb_dpk_, psb_int_heap
      implicit none 
      type(psb_d_csr_sparse_mat), intent(in)    :: a
      type(psb_d_coo_sparse_mat), intent(inout) :: trw
      integer, intent(in)                  :: i, m,jmin,jmax,jd
      integer, intent(inout)               :: ktrw,nlw,nup,jmaxup,info
      integer, intent(inout)               :: irwt(:)
      real(psb_dpk_), intent(inout)        :: nrmi,row(:)
      type(psb_int_heap), intent(inout)    :: heap
      real(psb_dpk_), intent(in), optional :: sign

    end subroutine mld_dinvt_copyin
  end interface

  interface mld_invt_inv
    subroutine mld_dinvt_inv(thres,i,nrmi,row,heap,irwt,uia1,uia2,uaspk,&
         & nidx,idxs,info)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      implicit none 
      ! Arguments
      type(psb_int_heap), intent(inout)   :: heap 
      integer, intent(in)                 :: i
      integer, intent(inout)              :: nidx,info
      integer, intent(inout)              :: irwt(:) 
      real(psb_dpk_), intent(in)          :: thres,nrmi
      integer, allocatable, intent(inout) :: idxs(:)
      integer, intent(in)                 :: uia1(:),uia2(:)
      real(psb_dpk_), intent(in)          :: uaspk(:)
      real(psb_dpk_), intent(inout)       :: row(:)

    end subroutine mld_dinvt_inv
  end interface

  interface mld_invt_copyout
    subroutine mld_dinvt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
         & nidx,idxs,l2,uia1,uia2,uaspk,info)

      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap

      implicit none 

      ! Arguments
      integer, intent(in)                       :: fill_in,i,m,nidx,nlw,nup,jmaxup
      integer, intent(in)                       :: idxs(:)
      integer, intent(inout)                    :: l2, info
      integer, allocatable, intent(inout)       :: uia1(:),uia2(:)
      real(psb_dpk_), intent(in)                :: thres,nrmi
      real(psb_dpk_),allocatable, intent(inout) :: uaspk(:)
      real(psb_dpk_), intent(inout)             :: row(:)

    end subroutine mld_dinvt_copyout
  end interface

  interface  mld_sparse_invt
    subroutine mld_dsparse_invt(n,a,z,nzrmax,sp_thresh,info)
      use psb_base_mod, only : psb_dspmat_type, psb_dpk_, psb_int_heap
      implicit none 
      integer, intent(in)                  :: n
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_dspmat_type), intent(inout) :: z
      integer, intent(in)                  :: nzrmax
      real(psb_dpk_), intent(in)           :: sp_thresh
      integer, intent(out)                 :: info
    end subroutine mld_dsparse_invt
  end interface



contains


end module mld_d_ainv_bld_mod

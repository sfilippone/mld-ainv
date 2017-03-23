!  
!   
!                       MLD-AINV: Approximate Inverse plugin for
!                             MLD2P4  version 2.0
!    
!    (C) Copyright 2012
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
subroutine psb_d_dsc_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_dsc_mat_mod, psb_protect_name => psb_d_dsc_csput
  implicit none 

  class(psb_d_dsc_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='d_dsc_csput'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)

  call psb_erractionsave(err_act)
  info = psb_success_


  info = psb_err_missing_override_method_
  call psb_errpush(info,name)
  goto 9999
!!$  if (nz <= 0) then 
!!$    info = psb_err_iarg_neg_
!!$    int_err(1)=1
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$  if (size(ia) < nz) then 
!!$    info = psb_err_input_asize_invalid_i_
!!$    int_err(1)=2
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  if (size(ja) < nz) then 
!!$    info = psb_err_input_asize_invalid_i_
!!$    int_err(1)=3
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$  if (size(val) < nz) then 
!!$    info = psb_err_input_asize_invalid_i_
!!$    int_err(1)=4
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  end if
!!$
!!$  if (nz == 0) return
!!$
!!$  nza  = a%get_nzeros()
!!$
!!$  if (a%is_bld()) then 
!!$    ! Build phase should only ever be in COO
!!$    info = psb_err_invalid_mat_state_
!!$
!!$  else  if (a%is_upd()) then 
!!$    call  psb_d_dsc_srch_upd(nz,ia,ja,val,a,&
!!$         & imin,imax,jmin,jmax,info,gtl)
!!$
!!$    if (info /= psb_success_) then  
!!$
!!$      info = psb_err_invalid_mat_state_
!!$    end if
!!$
!!$  else 
!!$    ! State is wrong.
!!$    info = psb_err_invalid_mat_state_
!!$  end if
!!$  if (info /= psb_success_) then
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  end if
!!$
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
!!$
!!$
!!$contains
!!$
!!$  subroutine psb_d_dsc_srch_upd(nz,ia,ja,val,a,&
!!$       & imin,imax,jmin,jmax,info,gtl)
!!$
!!$    use psb_const_mod
!!$    use psb_realloc_mod
!!$    use psb_string_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$
!!$    class(psb_d_dsc_sparse_mat), intent(inout) :: a
!!$    integer, intent(in) :: nz, imin,imax,jmin,jmax
!!$    integer, intent(in) :: ia(:),ja(:)
!!$    real(psb_dpk_), intent(in) :: val(:)
!!$    integer, intent(out) :: info
!!$    integer, intent(in), optional  :: gtl(:)
!!$    integer  :: i,ir,ic, ilr, ilc, ip, &
!!$         & i1,i2,nr,nc,nnz,dupl,ng, nar, nac
!!$    integer              :: debug_level, debug_unit
!!$    character(len=20)    :: name='d_dsc_srch_upd'
!!$
!!$    info = psb_success_
!!$    debug_unit  = psb_get_debug_unit()
!!$    debug_level = psb_get_debug_level()
!!$
!!$    dupl = a%get_dupl()
!!$
!!$    if (.not.a%is_sorted()) then 
!!$      info = -4
!!$      return
!!$    end if
!!$
!!$    ilr = -1 
!!$    ilc = -1 
!!$    nnz = a%get_nzeros()
!!$    nar = a%get_nrows()
!!$    nac = a%get_ncols()
!!$
!!$    if (present(gtl)) then 
!!$      ng = size(gtl)
!!$
!!$      select case(dupl)
!!$      case(psb_dupl_ovwrt_,psb_dupl_err_)
!!$        ! Overwrite.
!!$        ! Cannot test for error, should have been caught earlier.
!!$
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
!!$            ir = gtl(ir)
!!$            ic = gtl(ic)
!!$            if ((ic > 0).and.(ic <= nac)) then 
!!$              i1 = a%icp(ic)
!!$              i2 = a%icp(ic+1)
!!$              nr=i2-i1
!!$
!!$              ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$              if (ip>0) then 
!!$                a%val(i1+ip-1) = val(i)
!!$              else
!!$                if (debug_level >= psb_debug_serial_) &
!!$                     & write(debug_unit,*) trim(name),&
!!$                     & ': Was searching ',ir,' in: ',i1,i2,&
!!$                     & ' : ',a%ia(i1:i2-1)
!!$                info = i
!!$                return
!!$              end if
!!$
!!$            else
!!$
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Discarding row that does not belong to us.'
!!$            end if
!!$          end if
!!$        end do
!!$
!!$      case(psb_dupl_add_)
!!$        ! Add
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
!!$            ir = gtl(ir)
!!$            ic = gtl(ic)
!!$            if ((ic > 0).and.(ic <= nac)) then 
!!$              i1 = a%icp(ic)
!!$              i2 = a%icp(ic+1)
!!$              nr=i2-i1
!!$
!!$              ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$              if (ip>0) then 
!!$                a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
!!$              else
!!$                if (debug_level >= psb_debug_serial_) &
!!$                     & write(debug_unit,*) trim(name),&
!!$                     & ': Was searching ',ir,' in: ',i1,i2,&
!!$                     & ' : ',a%ia(i1:i2-1)
!!$                info = i
!!$                return
!!$              end if
!!$            else
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Discarding row that does not belong to us.'
!!$            end if
!!$
!!$          end if
!!$        end do
!!$
!!$      case default
!!$        info = -3
!!$        if (debug_level >= psb_debug_serial_) &
!!$             & write(debug_unit,*) trim(name),&
!!$             & ': Duplicate handling: ',dupl
!!$      end select
!!$
!!$    else
!!$
!!$      select case(dupl)
!!$      case(psb_dupl_ovwrt_,psb_dupl_err_)
!!$        ! Overwrite.
!!$        ! Cannot test for error, should have been caught earlier.
!!$
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$
!!$          if ((ic > 0).and.(ic <= nac)) then 
!!$            i1 = a%icp(ic)
!!$            i2 = a%icp(ic+1)
!!$            nr=i2-i1
!!$
!!$            ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$            if (ip>0) then 
!!$              a%val(i1+ip-1) = val(i)
!!$            else
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Was searching ',ir,' in: ',i1,i2,&
!!$                   & ' : ',a%ia(i1:i2-1)
!!$              info = i
!!$              return
!!$            end if
!!$
!!$          else
!!$            if (debug_level >= psb_debug_serial_) &
!!$                 & write(debug_unit,*) trim(name),&
!!$                 & ': Discarding col that does not belong to us.'
!!$          end if
!!$
!!$        end do
!!$
!!$      case(psb_dupl_add_)
!!$        ! Add
!!$        ilr = -1 
!!$        ilc = -1 
!!$        do i=1, nz
!!$          ir = ia(i)
!!$          ic = ja(i) 
!!$          if ((ic > 0).and.(ic <= nac)) then 
!!$            i1 = a%icp(ic)
!!$            i2 = a%icp(ic+1)
!!$            nr=i2-i1
!!$
!!$            ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
!!$            if (ip>0) then 
!!$              a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
!!$            else
!!$              if (debug_level >= psb_debug_serial_) &
!!$                   & write(debug_unit,*) trim(name),&
!!$                   & ': Was searching ',ir,' in: ',i1,i2,&
!!$                   & ' : ',a%ia(i1:i2-1)
!!$              info = i
!!$              return
!!$            end if
!!$          else
!!$            if (debug_level >= psb_debug_serial_) &
!!$                 & write(debug_unit,*) trim(name),&
!!$                 & ': Discarding col that does not belong to us.'
!!$          end if
!!$        end do
!!$
!!$      case default
!!$        info = -3
!!$        if (debug_level >= psb_debug_serial_) &
!!$             & write(debug_unit,*) trim(name),&
!!$             & ': Duplicate handling: ',dupl
!!$      end select
!!$
!!$    end if
!!$
!!$  end subroutine psb_d_dsc_srch_upd

end subroutine psb_d_dsc_csput

!   
!   
!                             MLD2P4  version 2.0
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.0)
!    
!    (C) Copyright 2008,2009,2010,2012
!  
!                        Salvatore Filippone  University of Rome Tor Vergata
!                        Alfredo Buttari      CNRS-IRIT, Toulouse
!                        Pasqua D'Ambra       ICAR-CNR, Naples
!                        Daniela di Serafino  Second University of Naples
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
! File: ppde2d.f90
!
! Program: ppde2d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs. 
! 
!
! The PDE is a general second order equation in 2d
!
!   a1 dd(u)  a2 dd(u)   b1 d(u)   b2 d(u) 
! -   ------ -  ------   -----  +  ------  + c u = f
!      dxdx     dydy        dx       dy    
!
! with Dirichlet boundary conditions
!   u = g 
!
!  on the unit square  0<=x,y<=1.
!
!
! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
!
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to a BLOCK
! data distribution.
!
program ppde2d
  use psb_base_mod
  use mld_prec_mod
  use psb_krylov_mod
  use psb_util_mod
!  use pde2d_base_mod
  use pde2d_exp_mod
  use mld_d_invt_solver
  use mld_d_invk_solver
  use mld_d_ainv_solver
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype, renum
  character(len=5)  :: afmt
  integer   :: idim

  ! miscellaneous 
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2, tprec, gt1,gt2 


  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
  type(mld_dprec_type)  :: prec
  type(mld_d_invt_solver_type) :: invtsv
  type(mld_d_invk_solver_type) :: invksv
  type(mld_d_ainv_solver_type) :: ainvsv
  type ainvparms 
    character(len=12) :: alg, orth_alg
    integer           :: fill, inv_fill
    real(psb_dpk_)    :: thresh, inv_thresh
  end type ainvparms
  type(ainvparms)     :: parms
  real(psb_dpk_)      :: nu
  ! descriptor
  type(psb_desc_type)   :: desc_a, desc_b
  ! dense matrices
  real(psb_dpk_), allocatable :: b(:), x(:)
  type(psb_d_vect_type)      :: xv,bv, vtst, xg, bg 
  ! blacs parameters
  integer            :: ictxt, iam, np

  ! solver parameters
  integer            :: iter, itmax,itrace, istopc, irst,giter, nr
  integer(psb_long_int_k_) :: amatsize, precsize, descsize, d2size, precnz,amatnz
  real(psb_dpk_)   :: err, eps, gerr, amxval
  integer, allocatable :: perm(:)

  ! other variables
  logical            :: pdump=.false.
  integer            :: info, i
  character(len=20)  :: name,ch_err
  character(len=40)  :: fname

  info=psb_success_

  
  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='ppde2d'
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
  !
  !  get parameters
  !
  call get_parms(ictxt,kmethd,ptype,afmt,renum,&
       & idim,istopc,itmax,itrace,irst,nu,parms)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !
  call pde_set_parm(nu)
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_gen_pde2d(ictxt,idim,a,bv,xv,desc_a,afmt,&
       & a1,a2,b1,b2,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (iam == psb_root_) write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  b=bv%get_vect()
  x=xv%get_vect()
  info = psb_get_errstatus()
  if (info /= 0) call psb_error(ictxt)
  if (np == 1) then 
    if (psb_toupper(trim(renum)) /= 'NONE') then 
      write(0,*) 'Going for renum "',trim(renum),'"'
      nr = a%get_nrows()
      call psb_mat_renum(renum,a,info,perm=perm)
      if (info /= 0) then 
        write(0,*) 'Error from RENUM',info
        goto 9999
      end if
      if (allocated(perm)) then
      else
        write(0,*) allocated(perm)
        goto 9999
      end if
      call psb_gelp('N',perm(1:nr),b(1:nr),info)
    end if
  end if
  info = psb_get_errstatus()
  if (info /= 0) call psb_error(ictxt)
  t1 = psb_wtime()
  call bv%bld(b)


  if (iam == psb_root_) write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  !
  !  prepare the preconditioner.
  !  
  if(iam == psb_root_) write(psb_out_unit,'("Setting preconditioner to : ",a)')ptype
  if (psb_toupper(trim(ptype)) == "AINV") then 
    call mld_precinit(prec,"BJAC",info)
  else
    call mld_precinit(prec,ptype,info)
  end if

  info = psb_get_errstatus()
  if (info /= 0) call psb_error(ictxt)

  if (psb_toupper(ptype) == 'DIAG') then 
    call prec%set(mld_smoother_sweeps_,1,info)
  end if
  if (psb_toupper(ptype) == 'AINV') then 
    select case (psb_toupper(parms%alg)) 
    case ('INVK') 
      call prec%set(invksv,info) 
    case ('INVT') 
      call prec%set(invtsv,info) 
    case ('AINV') 
      call prec%set(ainvsv,info) 
    end select
    call prec%set('ainv_alg',  parms%orth_alg,  info)
    call prec%set('sub_fillin',  parms%fill,    info)
    call prec%set('sub_iluthrs', parms%thresh,  info)
    call prec%set('inv_fillin', parms%inv_fill, info)
    call prec%set('inv_thresh', parms%inv_thresh, info)
  end if
  info = psb_get_errstatus()
  if (info /= 0) call psb_error(ictxt)

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call mld_precbld(a,desc_a,prec,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_precbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  tprec = psb_wtime()-t1

  call psb_amx(ictxt,tprec)


  if (pdump) then 
    call prec%dump(info,prefix='mlainv-',head='AINV prec',&
         & solver=.true., smoother=.false.)
    call a%print('amat.mtx',head='ELL mat')
  end if

  if (iam == psb_root_) write(psb_out_unit,'("Preconditioner time : ",es12.5)')tprec
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  amatsize = psb_sizeof(a)
  descsize = psb_sizeof(desc_a)
  precsize = mld_sizeof(prec)
  amatnz   = a%get_nzeros()
  precnz   = prec%get_nzeros()
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  call psb_sum(ictxt,amatnz)
  call psb_sum(ictxt,precnz)
  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Renumbering algorithm  :      ",a)') trim(renum)
    write(psb_out_unit,'("Total nonzeros          for A:      ",i12)')amatnz
    write(psb_out_unit,'("Total nonzeros          for PREC:   ",i12)')precnz      

    write(psb_out_unit,'("Total memory occupation for A:      ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for PREC:   ",i12)')precsize    
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(psb_out_unit,'("Storage type for DESC_A: ",a)') desc_a%indxmap%get_fmt()
    write(psb_out_unit,'("Storage type for DESC_B: ",a)') desc_b%indxmap%get_fmt()
    write(psb_out_unit,'(" ")')
    call mld_precdescr(prec,info) 
  end if
  eps = 1.d-7
  call psb_barrier(ictxt)
  !call psb_set_debug_level(psb_debug_ext_)
  t1 = psb_wtime()  
  call psb_krylov(kmethd,a,prec,bv,xv,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)     
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  call psb_set_debug_level(0)
  call psb_amx(ictxt,t2)
  call psb_barrier(ictxt)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to solve matrix          : ",es12.5)')t2
    write(psb_out_unit,'("Time per iteration            : ",es12.5)')t2/iter
    write(psb_out_unit,'("Number of iterations          : ",i0)')iter
    write(psb_out_unit,'("Convergence indicator on exit : ",es12.5)')err
    write(psb_out_unit,'("Info  on exit                 : ",i0)')info
  end if

  !  
  !  cleanup storage and exit
  !
  call psb_gefree(bv,desc_a,info)
  call psb_gefree(xv,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call mld_precfree(prec,info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

9999 continue
  if(info /= psb_success_) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,kmethd,ptype,afmt,renum,&
       & idim,istopc,itmax,itrace,irst,nu,parms)
    integer      :: ictxt
    character(len=*) :: kmethd, ptype, afmt,renum
    integer      :: idim, istopc,itmax,itrace,irst
    real(psb_dpk_)  :: nu
    type(ainvparms) :: parms
    integer      :: np, iam
    integer      :: intbuf(10), ip

    call psb_info(ictxt, iam, np)

    if (iam == 0) then
      read(psb_inp_unit,*) kmethd
      read(psb_inp_unit,*) ptype
      read(psb_inp_unit,*) afmt
      read(psb_inp_unit,*) renum
      read(psb_inp_unit,*) idim
      read(psb_inp_unit,*) nu
      read(psb_inp_unit,*) istopc
      read(psb_inp_unit,*) itmax
      read(psb_inp_unit,*) itrace
      read(psb_inp_unit,*) irst
      read(psb_inp_unit,*) parms%alg
      read(psb_inp_unit,*) parms%fill
      read(psb_inp_unit,*) parms%inv_fill
      read(psb_inp_unit,*) parms%thresh
      read(psb_inp_unit,*) parms%inv_thresh
      read(psb_inp_unit,*) parms%orth_alg
      
      write(psb_out_unit,'("Solving matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4)')idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)')np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'("Preconditioner       : ",a)') ptype
      write(psb_out_unit,'("Iterative method     : ",a)') kmethd
      write(psb_out_unit,'(" ")')
    end if
    call psb_bcast(ictxt,kmethd)
    call psb_bcast(ictxt,afmt)
    call psb_bcast(ictxt,renum)
    call psb_bcast(ictxt,ptype)
    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,nu)
    call psb_bcast(ictxt,istopc)
    call psb_bcast(ictxt,itmax)
    call psb_bcast(ictxt,itrace)
    call psb_bcast(ictxt,irst)
    call psb_bcast(ictxt,parms%alg)
    call psb_bcast(ictxt,parms%fill)
    call psb_bcast(ictxt,parms%inv_fill)
    call psb_bcast(ictxt,parms%thresh)
    call psb_bcast(ictxt,parms%inv_thresh)
    call psb_bcast(ictxt,parms%orth_alg)

    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  pde90 methd prec dim &
         &[istop itmax itrace]'  
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl' 
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '  
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage
end program ppde2d

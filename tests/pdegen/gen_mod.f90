module uniform
contains
function r8_uniform_01(seed)

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
end module uniform

module normal_rnd
use uniform
contains
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
end


end module normal_rnd


module psb_gen_min_mod
use normal_rnd
contains

subroutine psb_gen_min(ictxt,idim,a,bv,xv,desc_a,afmt,&
     & a1,a2,a3,b1,b2,b3,c,g,info,f,amold,vmold,imold,nrl)
  use psb_base_mod
  use genpde_mod, psb_protect_name => psb_d_gen_pde3d
  !
  !   Discretizes the partial differential equation
  ! 
  !   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
  ! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
  !      dxdx     dydy       dzdz        dx       dy         dz   
  !
  ! with Dirichlet boundary conditions
  !   u = g 
  !
  !  on the unit cube  0<=x,y,z<=1.
  !
  !
  ! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
  !
  implicit none
  procedure(d_func_3d)  :: b1,b2,b3,c,a1,a2,a3,g
  integer(psb_ipk_)     :: idim
  type(psb_dspmat_type) :: a
  type(psb_d_vect_type) :: xv,bv
  type(psb_desc_type)   :: desc_a
  integer(psb_ipk_)     :: ictxt, info
  character(len=*)      :: afmt
  procedure(d_func_3d), optional :: f
  class(psb_d_base_sparse_mat), optional :: amold
  class(psb_d_base_vect_type), optional :: vmold
  class(psb_i_base_vect_type), optional :: imold
  integer(psb_ipk_), optional :: nrl

  ! Local variables.

  integer(psb_ipk_), parameter :: nb=20
  type(psb_d_csc_sparse_mat)  :: acsc
  type(psb_d_coo_sparse_mat)  :: acoo
  type(psb_d_csr_sparse_mat)  :: acsr
  real(psb_dpk_)           :: zt(nb),x,y,z
  integer(psb_ipk_) :: m,n,nnz,glob_row,nlr,i,ii,ib,k
  integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
  integer(psb_ipk_) :: np, iam, nr, nt
  integer(psb_ipk_) :: icoeff
  integer(psb_ipk_), allocatable     :: irow(:),icol(:),myidx(:)
  real(psb_dpk_), allocatable :: val(:)
  ! deltah dimension of each grid cell
  ! deltat discretization time
  real(psb_dpk_)            :: deltah, sqdeltah, deltah2
  real(psb_dpk_), parameter :: rhs=0.e0,one=1.e0,zero=0.e0
  real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
  integer(psb_ipk_) :: err_act
  procedure(d_func_3d), pointer :: f_
  character(len=20)  :: name, ch_err,tmpfmt
  !mini app variables
  real(psb_dpk_), allocatable :: TX(:,:,:),TY(:,:,:), TZ(:,:,:), K_rnd(:,:,:,:) 
  real(psb_dpk_),allocatable :: x1(:), y1(:), z1(:), x2(:), y2(:), z2(:)
  integer(psb_ipk_) :: idx, mm, jj, kk, m_problem, k_eff, seed=4053

  info = psb_success_
  name = 'create_matrix'
  call psb_erractionsave(err_act)

  call psb_info(ictxt, iam, np)

  
  if (present(f)) then 
    f_ => f
  else
    f_ => d_null_func_3d
  end if

  deltah   = 1.e0/(idim)
  sqdeltah = deltah*deltah
  deltah2  = 2.e0* deltah

  ! initialize array descriptor and sparse matrix storage. provide an
  ! estimate of the number of non zeroes 

  m   = idim*idim*idim
  n   = m
  nnz = ((n*9)/(np))
  if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n

  if (present(nrl)) then 
    nr = nrl
  else
    !
    ! Using a simple BLOCK distribution.
    !
    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))
  end if
  nt = nr
  call psb_sum(ictxt,nt) 
  if (nt /= m) then 
    write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
    info = -1
    call psb_barrier(ictxt)
    call psb_abort(ictxt)
    return    
  end if
    
  call psb_barrier(ictxt)
  t0 = psb_wtime()
  call psb_cdall(ictxt,desc_a,info,nl=nr)
  if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
  ! define  rhs from boundary conditions; also build initial guess 
  if (info == psb_success_) call psb_geall(xv,desc_a,info)
  if (info == psb_success_) call psb_geall(bv,desc_a,info)

  call psb_barrier(ictxt)
  talc = psb_wtime()-t0

  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='allocation rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! we build an auxiliary matrix consisting of one row at a
  ! time; just a small matrix. might be extended to generate 
  ! a bunch of rows per call. 
  ! 
  allocate(val(20*nb),irow(20*nb),&
       &icol(20*nb),stat=info)
  if (info /= psb_success_ ) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  myidx = desc_a%get_global_indices()
  nlr = size(myidx)
  !Create TX, TY, TX, x1, y1, z1

  if (allocated(TX)) deallocate(TX)
  if (allocated(TY)) deallocate(TY)
  if (allocated(TZ)) deallocate(TZ)
  
  allocate(TX(idim+1,idim,idim), TY(idim,idim+1,idim), TZ(idim,idim,idim+1))
  TX=dzero
  TY=dzero
  TZ=dzero
  if (allocated(x1)) deallocate(x1)
  if (allocated(y1)) deallocate(y1)
  if (allocated(z1)) deallocate(z1)
  if (allocated(x2)) deallocate(x2)
  if (allocated(y2)) deallocate(y2)
  if (allocated(z2)) deallocate(z2)

  !Case 1: K is uniform = 1

  !TX(2:idim,:,:)=deltah
  !TY(:,2:idim,:)=deltah
  !TZ(:,:,2:idim)=deltah

  !Case 2: K is generated according to lognorm dist
  allocate(K_rnd(3,idim,idim,idim))
  do ii=1,idim
    do jj=1,idim
      do kk=1,idim
        do mm=1,3
          K_rnd(mm,kk,jj,ii)=exp(r8_normal_01(seed))
        enddo
      enddo
    enddo
  enddo

  TX(2:idim,:,:)=deltah2/(1./K_rnd(1,1:idim-1,:,:)+1./K_rnd(1,2:idim,:,:))
  TY(:,2:idim,:)=deltah2/(1./K_rnd(2,:,1:idim-1,:)+1./K_rnd(2,:,2:idim,:))
  TZ(:,:,2:idim)=deltah2/(1./K_rnd(3,:,:,1:idim-1)+1./K_rnd(3,:,:,2:idim))
 

  allocate(x1(m),y1(m),z1(m))
  allocate(x2(m),y2(m),z2(m))
  x1=reshape(TX(1:idim,:,:),(/m/))
  y1=reshape(TY(:,1:idim,:),(/m/))
  z1=reshape(TZ(:,:,1:idim),(/m/))
  x2=reshape(TX(2:idim+1,:,:),(/m/))
  y2=reshape(TY(:,2:idim+1,:),(/m/))
  z2=reshape(TZ(:,:,2:idim+1),(/m/))

  ! loop over rows belonging to current process in a block
  ! distribution.

  call psb_barrier(ictxt)
  t1 = psb_wtime()

  !rhs and xv
  zt(:)=0.e0
  do ii=1,nlr,nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    do k=1,ib
      idx = ii -1 + k
      !local matrix pointer
      glob_row=myidx(idx)
      if (glob_row == 1) zt(k)=1
      if (glob_row == m) zt(k)=-1
    enddo
    if(info /= psb_success_) exit
    call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
    if(info /= psb_success_) exit
    zt(:)=0.e0
    call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
    if(info /= psb_success_) exit
  enddo
  !Main diagonal
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      glob_row=myidx(idx)
      val(k) = x1(glob_row)+x2(glob_row)+y1(glob_row)+y2(glob_row)+z1(glob_row)+z2(glob_row)
      if (glob_row == 1 ) val(k) = val(k) + sum(K_rnd(:,1,1,1))
      irow(k)= glob_row
      icol(k)= glob_row
      icoeff = icoeff+1
    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
  enddo
  
  !Diagonal -1 
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    k_eff=0
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      glob_row=myidx(idx)
      if (glob_row == 64) print*,'hey sista'
      if (glob_row > m) cycle
      if (glob_row -1 < 1) cycle
      k_eff=k_eff + 1
      val(k_eff) = -x2(glob_row-1)
      irow(k_eff)= glob_row
      icol(k_eff)= glob_row-1
      icoeff = icoeff+1
      !if (iam==0) print*, glob_row, glob_row-1,val(k_eff), idx, k_eff

      !if (iam==0) print*,'pppp',icoeff-1,glob_row, glob_row-1, val(k) 
    end do
    call psb_spins(icoeff-1,irow(1:icoeff-1),icol(1:icoeff-1),val(1:icoeff-1),a,desc_a,info)
    !todo
  enddo

  !Diagonal 1 
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    k_eff = 0
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      glob_row=myidx(idx)
      if ((glob_row + 1)>m) cycle
      k_eff=k_eff + 1
      val(k_eff) = -x1(glob_row+1)
      irow(k_eff)= glob_row
      icol(k_eff)= glob_row+1
      icoeff = icoeff+1
      !if (iam==0) print*,'jjjj',nlr, idx, glob_row, glob_row + 1
    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
  enddo


  !Diagonal -idim
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    k_eff=0
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      glob_row=myidx(idx)
      !if ((glob_row + idim)>m) cycle
      if (glob_row > m) cycle
      if (glob_row -idim < 1) cycle
      k_eff=k_eff + 1
      val(k_eff) = -y2(glob_row-idim)
      irow(k_eff)= glob_row
      icol(k_eff)= glob_row-idim
      !if (iam==0) print*,'éééééééé', glob_row, glob_row-idim,val(k_eff), idx, k_eff
      icoeff = icoeff+1
    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
  enddo

  !Diagonal idim
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    k_eff=0
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      glob_row=myidx(idx)
      if ((glob_row + idim)>m) cycle
      k_eff=k_eff + 1
      val(k_eff) = -y1(glob_row+idim)
      irow(k_eff)= glob_row
      icol(k_eff)= glob_row+idim 
      icoeff = icoeff+1
    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
  enddo
  !Diagonal -idim^2
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    k_eff=0
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      !print*,idx,nb,nlr-idim*idim, ii
      glob_row=myidx(idx)
      if (glob_row > m) cycle
      if (glob_row -idim*idim < 1) cycle
      !if ((glob_row + idim*idim)>m) cycle
      k_eff=k_eff + 1
      val(k_eff) = -z2(glob_row-idim*idim)
      irow(k_eff)= glob_row
      icol(k_eff)= glob_row-idim*idim
      icoeff = icoeff+1
      !if (iam==0) print*,'éééééééé', glob_row, glob_row-idim*idim,val(k_eff), idx, k_eff
      !print*, iam, idx, glob_row, glob_row - idim*idim
    end do
    call psb_spins(icoeff-1,irow(1:icoeff-1),icol(1:icoeff-1),val(1:icoeff-1),a,desc_a,info)
  enddo

  !Diagonal idim^2
  do ii =1,nlr, nb
    ib = min(nb, nlr-ii+1)
    icoeff=1
    k_eff=0
    do k=1, ib
      idx=ii-1+k
      !local matrix pointer
      glob_row=myidx(idx)
      if ((glob_row + idim*idim)>m) cycle
      k_eff=k_eff + 1
      val(k_eff) = -z1(glob_row+idim*idim)
      irow(k_eff)= glob_row
      icol(k_eff)= glob_row+idim*idim
      icoeff = icoeff+1
      !if (iam == 0) print*,'ppppp',idx,ii,nlr-idim*idim, glob_row, glob_row+idim*idim
    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
  enddo

  tgen = psb_wtime()-t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='insert rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  deallocate(val,irow,icol)

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_cdasb(desc_a,info,mold=imold)
  tcdasb = psb_wtime()-t1
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  if (info == psb_success_) then 
    if (present(amold)) then 
      call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,mold=amold)
    else
      call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    end if
    call a%clean_zeros(info)
  end if
  call psb_barrier(ictxt)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (info == psb_success_) call psb_geasb(xv,desc_a,info,mold=vmold)
  if (info == psb_success_) call psb_geasb(bv,desc_a,info,mold=vmold)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  tasb = psb_wtime()-t1
  call psb_barrier(ictxt)
  ttot = psb_wtime() - t0 

  call psb_amx(ictxt,talc)
  call psb_amx(ictxt,tgen)
  call psb_amx(ictxt,tasb)
  call psb_amx(ictxt,ttot)
  if(iam == psb_root_) then
    tmpfmt = a%get_fmt()
    write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")')&
         &   tmpfmt
    write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
    write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
    write(psb_out_unit,'("-desc asbly  time : ",es12.5)') tcdasb
    write(psb_out_unit,'("- mat asbly  time : ",es12.5)') tasb
    write(psb_out_unit,'("-total       time : ",es12.5)') ttot

  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_gen_min

end module psb_gen_min_mod


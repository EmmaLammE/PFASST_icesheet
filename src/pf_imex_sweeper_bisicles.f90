!!  IMEX Sweeper Module
!
! This file is part of LIBPFASST.
!
!>  Module for the IMEX Sweeper  of the  the derived sweeper class for doing IMEX sweeps for an equation of the form
!!         $$   y' = f_1(y) + f_2(y)  $$
!!  The \(f_1\) piece is treated explicitly and \(f_2\) implicitly
!!  Afer this sweeper is initialized (usually in main), the logical flags can be changed if desired
module pf_mod_imex_sweeper_bisicles
  use pf_mod_dtype
  use pf_mod_utils
  use pf_mod_timer

  implicit none

  !>  IMEX SDC sweeper type, extends abstract sweeper
  type, extends(pf_sweeper_t), abstract :: pf_imex_sweeper_bisicles_t
     real(pfdp), allocatable :: QtilE(:,:)   !!  Approximate explicit quadrature rule
     real(pfdp), allocatable :: QtilI(:,:)   !!  Approximate implicit quadrature rule
     real(pfdp), allocatable :: dtsdc(:)     !!  SDC step sizes
     real(pfdp), allocatable :: QdiffE(:,:)  !!  qmat-QtilE
     real(pfdp), allocatable :: QdiffI(:,:)  !!  qmat-QtilI

     logical    :: explicit  !!  True if there is an explicit piece (must set in derived sweeper)
     logical    :: implicit  !!  True if there an implicit piece (must set in derived sweeper)
     integer    :: m_sub     !!  Substep loop variable (useful in the function evaluation routines in derived sweepers)
     class(pf_encap_t), allocatable :: rhs   !! holds rhs for implicit solve

   contains
     procedure(pf_f_eval_p), deferred :: f_eval   !!  RHS function evaluations
     procedure(pf_f_comp_p), deferred :: f_comp   !!  Implicit solver
     !>  Set the generic functions
     procedure :: sweep      => imex_bisicles_sweep
     procedure :: initialize => imex_bisicles_initialize
     procedure :: evaluate   => imex_bisicles_evaluate
     procedure :: integrate  => imex_bisicles_integrate
     procedure :: residual   => imex_bisicles_residual
     procedure :: spreadq0   => imex_bisicles_spreadq0
     procedure :: compute_dt => imex_bisicles_compute_dt
     procedure :: evaluate_all => imex_bisicles_evaluate_all
     procedure :: destroy   => imex_bisicles_destroy
     procedure :: imex_bisicles_destroy
     procedure :: imex_bisicles_initialize
  end type pf_imex_sweeper_bisicles_t

  interface
     !>  The interface to the routine to compute the RHS function values
     !>  Evaluate f_piece(y), where piece is one or two
     subroutine pf_f_eval_p(this,y, t, level_index, f, c_AmrIceHolderPtr)
       use :: iso_c_binding
       import pf_imex_sweeper_bisicles_t, pf_encap_t, pfdp
       class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
       class(pf_encap_t), intent(in   )  :: y           !!  Argument for evaluation
       real(pfdp),        intent(in   )  :: t           !!  Time at evaluation
       integer,    intent(in   )         :: level_index !!  Level index
       class(pf_encap_t), intent(inout)  :: f           !!  RHS function value
       !integer,    intent(in   )         :: piece       !!  Which piece to evaluate
       type(c_ptr), intent(in)           :: c_AmrIceHolderPtr
     end subroutine pf_f_eval_p

     !>  The interface to the routine to do implicit solve 
     !>  i.e, solve the equation y - dtq*f_2(y) =rhs
     subroutine pf_f_comp_p(this,y, t, dtq, rhs, level_index, f, piece)
       import pf_imex_sweeper_bisicles_t, pf_encap_t, pfdp
       class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout)  :: y           !!  Solution of implicit solve
       real(pfdp),        intent(in   )  :: t           !!  Time of solve
       real(pfdp),        intent(in   )  :: dtq         !!  dt*quadrature weight
       class(pf_encap_t), intent(in   )  :: rhs         !!  RHS for solve
       integer,    intent(in   )         :: level_index !!  Level index
       class(pf_encap_t), intent(inout) :: f            !!  f_2 of solution y
       integer,    intent(in   ) :: piece               !!  Which piece to evaluate
     end subroutine pf_f_comp_p
  end interface

contains

  !> Perform nsweeps SDC sweeps on level level_index and set qend appropriately.
  subroutine imex_bisicles_sweep(this, pf, level_index, t0, dt,nsweeps, flags)
    use pf_mod_hooks

    !>  Inputs
    class(pf_imex_sweeper_bisicles_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to sweep
    real(pfdp),        intent(in   ) :: t0           !!  time at beginning of time step
    real(pfdp),        intent(in   ) :: dt           !!  time step size
    integer,           intent(in)    :: nsweeps      !!  number of sweeps to do
    integer, optional, intent(in   ) :: flags        !!  sweep specific flags

    !>  Local variables
    type(pf_level_t), pointer :: lev    !!  points to current level
    type(c_ptr)      :: c_AmrIceHolderPtr
    !type(bisicles_holder_ptr), pointer           :: pf_bisicles

    integer     :: m,n,k   !!  Loop variables
    real(pfdp)  :: t        !!  Time at nodes

    lev => pf%levels(level_index)   !  Assign level pointer
    
    c_AmrIceHolderPtr=pf%cptr_AmrIceHolder
    !pf_bisicles => cast_as_pf_bisicles_t(pf)

    !print *,'pf_imex_bisicles_sweeper.f90 0000 c_AmrIceHolderPtr ', c_AmrIceHolderPtr
    !call ABORT
    !print *, '------ adding rhs to q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1'
    sweeps: do k = 1,nsweeps   !!  Loop over sweeps
       pf%state%sweep=k
       call call_hooks(pf, level_index, PF_PRE_SWEEP)
       call pf_start_timer(pf, T_SWEEP,level_index)

       !  Add terms from previous iteration  (not passing CI tests)
       !do m = 1, lev%nnodes-1
       !   call lev%I(m)%setval(0.0_pfdp)
       !end do

        !if (this%explicit) call pf_apply_mat(lev%I, dt, this%QdiffE, lev%F(:,1), .false.)             
       !if (this%implicit) call pf_apply_mat(lev%I, dt, this%QdiffI, lev%F(:,2), .false.)
       ! compute integrals and add fas correction
       do m = 1, lev%nnodes-1
          call lev%I(m)%setval(0.0_pfdp)
          if (this%explicit) then
             do n = 1, lev%nnodes
                call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1))
             end do
          end if
          if (this%implicit) then
             do n = 1, lev%nnodes
                call lev%I(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2))
             end do
          end if
       end do

       !  Add the tau FAS correction
       if (level_index < pf%state%finest_level) then
          do m = 1, lev%nnodes-1
             call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
             if (m>1 .and. pf%use_Sform) then
                call lev%I(m)%axpy(-1.0_pfdp, lev%tauQ(m-1))
             end if
          end do
       end if
       
       !  Recompute the first function value if this is first sweep
       if (k .eq. 1) then
          !print *,'imex_sweeper_bisicles 1111111 '
          call lev%Q(1)%copy(lev%q0)
          if (this%explicit) then
             call pf_start_timer(pf,T_FEVAL,level_index)
             !print *,'imex_sweeper_bisicles 222222222 '
             call this%f_eval(lev%Q(1), t0, level_index, lev%F(1,1),c_AmrIceHolderPtr)
             call pf_stop_timer(pf,T_FEVAL,level_index)
          end if
          
          if (this%implicit) then
             call pf_start_timer(pf,T_FEVAL,level_index)
              call this%f_eval(lev%Q(1), t0, level_index, lev%F(1,2),c_AmrIceHolderPtr)
             call pf_stop_timer(pf,T_FEVAL,level_index)
          end if
       end if

       t = t0
       ! do the sub-stepping in sweep

       substeps: do m = 1, lev%nnodes-1  !!  Loop over substeps
       !print *, '------ adding rhs to q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2 m',m
          t = t + dt*this%dtsdc(m)

          this%m_sub=m
          !>  Accumulate rhs
          call this%rhs%setval(0.0_pfdp)
          do n = 1, m
             if (this%explicit) &
                  !print *, ' before accumulatin rhs, a ',dt*this%QtilE(m,n)
                  !call lev%F(n,1)%eprint()
                  call this%rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1))
             if (this%implicit) &
                  call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2))
          end do
          !>  Add the integral term
          !print *, '    before adding lev%I '
          !call lev%I(m)%eprint()
          !call this%rhs%eprint()
          call this%rhs%axpy(1.0_pfdp, lev%I(m))
          !call this%rhs%eprint()

          !>  Add the starting value
          if (pf%use_Sform) then
             call this%rhs%axpy(1.0_pfdp, lev%Q(m))
          else
             call this%rhs%axpy(1.0_pfdp, lev%Q(1))
          end if

          !>  Solve for the implicit piece
          if (this%implicit) then
          !print *, '------ adding rhs to q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 3 implicit'
             call pf_start_timer(pf,T_FCOMP,level_index)
             call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, level_index,lev%F(m+1,2),2)
             call pf_stop_timer(pf,T_FCOMP,level_index)
          else
             !print *,'imex_sweeper_bisicles 22222 '
             !print *, '------ adding rhs to q !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 3 explicit'
             !print *, 'before adding, rhs     &   Q(m+1)'
             !call this%rhs%eprint()
             !call lev%Q(m+1)%eprint()
             call lev%Q(m+1)%copy(this%rhs)
             !print *, 'after adding, Q(m+1)'
             !call lev%Q(m+1)%eprint()
             !call ABORT
          end if
          !>  Compute explicit function on new value
          if (this%explicit) then
             !print *,'imex_sweeper_bisicles 22222 '
             call pf_start_timer(pf,T_FEVAL,level_index)
             call this%f_eval(lev%Q(m+1), t, level_index, lev%F(m+1,1),c_AmrIceHolderPtr)
             call pf_stop_timer(pf,T_FEVAL,level_index)
          end if
            

       end do substeps !!  End substep loop
       
       call pf_residual(pf, level_index, dt)
       !print *,'imex_sweeper_bisicles 33333 '
       call lev%qend%copy(lev%Q(lev%nnodes))
       call pf_stop_timer(pf, T_SWEEP,level_index)

       call call_hooks(pf, level_index, PF_POST_SWEEP)
    end do sweeps  !  End loop on sweeps

  end subroutine imex_bisicles_sweep

  !> Subroutine to initialize matrices and space for sweeper
  subroutine imex_bisicles_initialize(this, pf, level_index)
    use pf_mod_quadrature
    class(pf_imex_sweeper_bisicles_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize

    type(pf_level_t), pointer  :: lev    !  Current level

    integer    ::  nnodes,ierr
    lev => pf%levels(level_index)   !  Assign level pointer
    this%npieces = 2

    !  The default is to use both pieces, but can be overriddent in local sweeper
    this%explicit=.TRUE.
    this%implicit=.TRUE.
    
    nnodes = lev%nnodes
    allocate(this%QdiffE(nnodes-1,nnodes),stat=ierr)
    if (ierr /=0) stop "allocate fail in imex_initialize for QdiffE"
    allocate(this%QdiffI(nnodes-1,nnodes),stat=ierr)
    if (ierr /=0) stop "allocate fail in imex_initialize for QdiffI"
    allocate(this%QtilE(nnodes-1,nnodes),stat=ierr)
    if (ierr /=0) stop "allocate fail in imex_initialize for QtilE"
    allocate(this%QtilI(nnodes-1,nnodes),stat=ierr)
    if (ierr /=0) stop "allocate fail in imex_initialize for QtilI"
    allocate(this%dtsdc(nnodes-1),stat=ierr)
    if (ierr /=0) stop "allocate fail in imex_initialize for dtsdc"

    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    !>  Array of substep sizes
    this%dtsdc = lev%sdcmats%qnodes(2:nnodes) - lev%sdcmats%qnodes(1:nnodes-1)

    ! Implicit matrix
    if (this%use_LUq) then
       this%QtilI = lev%sdcmats%qmatLU
    else
       this%QtilI =  lev%sdcmats%qmatBE
    end if

    ! Explicit matrix
    this%QtilE =  lev%sdcmats%qmatFE

    this%QdiffE = lev%sdcmats%qmat-this%QtilE
    this%QdiffI = lev%sdcmats%qmat-this%QtilI

    if (pf%use_Sform) then
          this%QdiffE(2:nnodes-1,:) = this%QdiffE(2:nnodes-1,:)- this%QdiffE(1:nnodes-2,:)
          this%QdiffI(2:nnodes-1,:) = this%QdiffI(2:nnodes-1,:)- this%QdiffI(1:nnodes-2,:)
          this%QtilE(2:nnodes-1,:) = this%QtilE(2:nnodes-1,:)- this%QtilE(1:nnodes-2,:)
          this%QtilI(2:nnodes-1,:) = this%QtilI(2:nnodes-1,:)- this%QtilI(1:nnodes-2,:)
    end if
    !>  Make space for rhs
    call lev%ulevel%factory%create_single(this%rhs, level_index,   lev%lev_shape)

  end subroutine imex_bisicles_initialize

  !>  Subroutine to deallocate sweeper
  subroutine imex_bisicles_destroy(this, pf,level_index)
    class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
    type(pf_pfasst_t),  target,  intent(inout) :: pf
    integer,              intent(in)    :: level_index

    type(pf_level_t), pointer  :: lev        !  Current level
    lev => pf%levels(level_index)   !  Assign level pointer

    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)

    call lev%ulevel%factory%destroy_single(this%rhs)

  end subroutine imex_bisicles_destroy


  !> Subroutine to compute  Picard integral of function values
  subroutine imex_bisicles_integrate(this,pf,level_index, qSDC, fSDC, dt, fintSDC, flags)
    class(pf_imex_sweeper_bisicles_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    class(pf_encap_t), intent(in   ) :: qSDC(:)      !!  Solution values
    class(pf_encap_t), intent(in   ) :: fSDC(:, :)   !!  RHS Function values
    real(pfdp),        intent(in   ) :: dt           !!  Time step
    class(pf_encap_t), intent(inout) :: fintSDC(:)   !!  Integral from t_n to t_m
    integer, optional, intent(in   ) :: flags

    integer :: n, m
    type(pf_level_t), pointer :: lev        !  Current level
    lev => pf%levels(level_index)   !  Assign level pointer
    
    do n = 1, lev%nnodes-1
       call fintSDC(n)%setval(0.0_pfdp)
       do m = 1, lev%nnodes
          if (this%explicit) &
               !print *,'imex_sweeper in interpolate -------'
               call fintSDC(n)%axpy(dt*lev%sdcmats%qmat(n,m), fSDC(m,1))
          if (this%implicit) &
               call fintSDC(n)%axpy(dt*lev%sdcmats%qmat(n,m), fSDC(m,2))
       end do
    end do


!    if (this%explicit) call pf_apply_mat(fintSDC, dt, lev%sdcmats%Qmat, fSDC(:,1), .false.)    
!    if (this%implicit) call pf_apply_mat(fintSDC, dt, lev%sdcmats%Qmat, fSDC(:,2), .false.)    
  end subroutine imex_bisicles_integrate


  !> Subroutine to compute  Residual
  subroutine imex_bisicles_residual(this, pf, level_index, dt, flags)
    class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: dt           !!  Time step
    integer, intent(in), optional   :: flags

    integer :: m
    type(pf_level_t), pointer :: lev        !  Current level
    lev => pf%levels(level_index)   !  Assign level pointer
    !print *, 'begin integrate in residual..... lev%nnodes-1 dt ',lev%nnodes-1, dt
    !print *, '   lev%I(2) before integrate '
    !call lev%I(2)%eprint()
    call imex_bisicles_integrate(this,pf,level_index, lev%Q, lev%F, dt, lev%I, flags)
    !print *, '   lev%I(2) after integrate '
    !call lev%I(2)%eprint()
    !> add tau if it exists
    if (lev%index < pf%state%finest_level) then    
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), flags)
       end do
    end if
    do m = 1, lev%nnodes-1    
       !print *,'imex_sweeper_bisicles 44444 '  
       call lev%R(m)%copy(lev%I(m))
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
       !print *,'imex_sweeper_bisicles 55555 '
       if (present(flags)) then
          if (flags .eq. 0) then
             call lev%R(m)%axpy(1.0_pfdp, lev%q0)
          else
             call lev%R(m)%axpy(1.0_pfdp, lev%Q(1))
          end if
          !print *,'imex_sweeper_bisicles 66666 '
       else
          call lev%R(m)%axpy(1.0_pfdp, lev%Q(1))
          !print *,'imex_sweeper_bisicles 77777 '
       end if
    end do
    
!    call pf_generic_residual(this, pf, level_index, dt)
  end subroutine imex_bisicles_residual


  subroutine imex_bisicles_spreadq0(this, pf,level_index, t0, flags, step)
    class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: t0
    integer, optional,   intent(in)    :: flags, step

    call pf_generic_spreadq0(this, pf,level_index, t0)
  end subroutine imex_bisicles_spreadq0

  subroutine imex_bisicles_compute_dt(this,pf,level_index,  t0, dt,flags)
    class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
    type(pf_pfasst_t), target, intent(inout) :: pf
    integer,              intent(in)    :: level_index
    real(pfdp),        intent(in   ) :: t0
    real(pfdp),        intent(inout) :: dt
    integer, optional,   intent(in)    :: flags

    type(pf_level_t),    pointer :: lev
    lev => pf%levels(level_index)   !!  Assign level pointer
    !  Do nothing now
    return
  end subroutine imex_bisicles_compute_dt

  !> Subroutine to evaluate function value at node m
  subroutine imex_bisicles_evaluate(this, pf,level_index, t, m, flags, step)
    class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: t    !!  Time at which to evaluate
    integer,           intent(in   ) :: m    !!  Node at which to evaluate
    integer, intent(in), optional   :: flags, step

    type(pf_level_t), pointer :: lev        !  Current level
    type(c_ptr)      :: c_AmrIceHolderPtr
    lev => pf%levels(level_index)   !  Assign level pointer

    c_AmrIceHolderPtr=pf%cptr_AmrIceHolder
    !print *,'pf_imex_bisicles_sweeper.f90 1111 c_AmrIceHolderPtr ', c_AmrIceHolderPtr

    if (this%explicit) then
       !print *,'imex_sweeper_bisicles 33333 m ', m
       call pf_start_timer(pf,T_FEVAL,level_index)
       !call pf%levels(level_index)%Q(2)%eprint()
       call this%f_eval(lev%Q(m), t, level_index, lev%F(m,1),c_AmrIceHolderPtr)
       call pf_stop_timer(pf,T_FEVAL,level_index)
    end if
    
    if (this%implicit) then
       call pf_start_timer(pf,T_FEVAL,level_index)       
       call this%f_eval(lev%Q(m), t, level_index, lev%F(m,2),c_AmrIceHolderPtr)
       call pf_stop_timer(pf,T_FEVAL,level_index)
    end if
    
  end subroutine imex_bisicles_evaluate

  !> Subroutine to evaluate the function values at all nodes
  subroutine imex_bisicles_evaluate_all(this, pf,level_index, t, flags, step)
    class(pf_imex_sweeper_bisicles_t),  intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf    !!  PFASST structure
    integer,           intent(in)    :: level_index  !!  level on which to initialize
    real(pfdp),        intent(in   ) :: t(:)  !!  Array of times at each node
    integer, intent(in), optional   :: flags, step
    
    call pf_generic_evaluate_all(this, pf,level_index, t)
  end subroutine imex_bisicles_evaluate_all

end module pf_mod_imex_sweeper_bisicles

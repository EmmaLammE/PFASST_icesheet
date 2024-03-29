!!  Restriction operators
!
! This file is part of LIBPFASST.
!
module pf_mod_restrict
  !!  Module to restrict solutions between pfasst levels and create the FAS tau correction
  use pf_mod_dtype
  use pf_mod_timer
  use pf_mod_hooks
  use pf_mod_utils
  implicit none
contains



  subroutine restrict_time_space_fas(pf, t0, dt, level_index, flags, mystep)
    !! Restrict (in time and space) fine level to coarse and set coarse level FAS correction.
    !!
    !! The coarse function values are re-evaluated after restriction.
    !! Note that even if the number of variables and nodes is the same,
    !! we should still compute the FAS correction since the function
    !! evaluations may be different.
    type(pf_pfasst_t), intent(inout),target :: pf
    real(pfdp),        intent(in)    :: t0            !!  time at beginning of step
    real(pfdp),        intent(in)    :: dt            !!  time step
    integer,           intent(in)    :: level_index   !! defines which level to restrict
    integer, optional, intent(in)    :: flags, mystep    

    !>  Local variables
    class(pf_level_t), pointer :: c_lev    
    class(pf_level_t), pointer :: f_lev

    integer    :: m, step,ierr

    real(pfdp), allocatable :: c_times(:)  !!  Simulation time at coarse nodes  
    real(pfdp), allocatable :: f_times(:)  !!  Simulation time at fine nodes
    
    
    f_lev => pf%levels(level_index);
    c_lev => pf%levels(level_index-1)
    
    step = pf%state%step+1
    if(present(mystep)) step = mystep
    
    call call_hooks(pf, level_index, PF_PRE_RESTRICT_ALL)
    call pf_start_timer(pf, T_RESTRICT, level_index)
    

    allocate(c_times(c_lev%nnodes),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
    allocate(f_times(f_lev%nnodes),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    

    !> restrict q's and recompute f's
    c_times = t0 + dt*c_lev%nodes
    f_times = t0 + dt*f_lev%nodes

    call restrict_ts(f_lev, c_lev, f_lev%Q, c_lev%Q, f_times, flags)

    !>  Recompute the functions
    !call pf%levels(level_index-1)%Q(2)%eprint()
    call c_lev%ulevel%sweeper%evaluate_all(pf,level_index-1, c_times, flags=flags, step=step)

    !>  Compute  FAS correction
    do m = 1, c_lev%nnodes-1
       call c_lev%tauQ(m)%setval(0.0_pfdp, flags)
    end do
    if (pf%state%iter >= pf%taui0)  then
      !print *,'pf_restrict 22222 '
       ! compute '0 to node' integral on the coarse level
      call c_lev%ulevel%sweeper%integrate(pf,level_index-1,  c_lev%Q, &
           c_lev%F, dt, c_lev%I, flags)
       ! compute '0 to node' integral on the fine level
      call f_lev%ulevel%sweeper%integrate(pf,level_index, f_lev%Q, &
        f_lev%F, dt, f_lev%I, flags)
       !  put tau in on fine level
      if (level_index < pf%state%finest_level) then
          do m = 1, f_lev%nnodes-1
             call f_lev%I(m)%axpy(1.0_pfdp, f_lev%tauQ(m), flags)
             if (level_index .eq. 2) then
                !print *,'     pf_restrict put tau in on fine level f_lev%I(m) '
                !call f_lev%I(m)%eprint()
             end if
          end do
       end if
       !  Subtract coarse integral
       do m = 1, c_lev%nnodes-1
          call c_lev%tauQ(m)%axpy(-1.0_pfdp, c_lev%I(m), flags)
       end do
       !print *,'     pf_restrict put tau in on coarse level  '
       !call c_lev%I(1)%eprint()

       ! restrict '0 to node' integral on the fine level  in time and space
       call restrict_ts_integral(f_lev, c_lev, f_lev%I, c_lev%I, f_times, flags)
       !print *,'     pf_restrict after putting f_lev%I on c_lev%I f_lev%I(m) '
                !call f_lev%I(1)%eprint()
                !call c_lev%I(1)%eprint()

       ! Add fine restriction of fine integral (stored on coarse)
       do m = 1, c_lev%nnodes-1
          !print *, '         pf_restrict add fine restrict to fine integral, tauQ(m), I(m)'
          !call c_lev%tauQ(m)%eprint()
          !call c_lev%I(m)%eprint()
          call c_lev%tauQ(m)%axpy(1.0_pfdp, c_lev%I(m), flags)
          !call c_lev%tauQ(m)%eprint()
       end do
       !print *,'     pf_restrict put tau in on coarse level  '
       !call c_lev%tauQ(1)%eprint()
!!$       if (pf%use_Sform) then
!!$          do m = c_lev%nnodes-1,2,-1
!!$!             call c_lev%tauQ(m)%axpy(-1.0_pfdp, c_lev%tauQ(m-1), flags)
!!$          end do
!!$       end if
       
    end if

    call pf_stop_timer(pf, T_RESTRICT,level_index)
    call call_hooks(pf, level_index, PF_POST_RESTRICT_ALL)

    deallocate(c_times)
    deallocate(f_times)

  end subroutine restrict_time_space_fas


  subroutine restrict_ts(f_lev, c_lev, f_encap_array, c_encap_array, f_time, flags)

    !! Restrict (in time and space) f_sol_array  to c_sol_array
    !! This version is for point values (either functions or solutions)
    
    class(pf_level_t),  intent(inout) :: f_lev   !!   pointer to fine level
    class(pf_level_t),  intent(inout) :: c_lev   !!   pointer to coarse level
    class(pf_encap_t),  intent(inout) :: f_encap_array(:)   !! array of fine level data to be restricted
    class(pf_encap_t),  intent(inout) :: c_encap_array(:)   !! array of coarse level data to be computed
    real(pfdp),         intent(in) :: f_time(:)             !! time at the fine nodes
    integer, optional, intent(in)    :: flags    

    integer :: m,j
    integer :: f_nnodes,c_nnodes


    f_nnodes = f_lev%nnodes ! num of SDC nodes on fine level (can think of num of time steps)
    c_nnodes = c_lev%nnodes ! num of SDC nodes on coarse level (can think of num of time steps)
   !  print *, 'EL DEBUG --        restrict ts, f_nnodes:', f_nnodes, 'c_nnodes:',c_nnodes
    !!  Create a temp array for the spatial restriction
    if (f_lev%restrict_workspace_allocated   .eqv. .false.) then      
       call c_lev%ulevel%factory%create_array(f_lev%f_encap_array_c, f_nnodes, c_lev%index, c_lev%lev_shape)
       f_lev%restrict_workspace_allocated  = .true.
    end if
 
    !  spatial restriction
    do m = 1, f_nnodes
       !print *, '                    in restrict_ts, before spatial restrict, f_encap_array(m), m',m
       !call f_encap_array(m)%eprint()
       call f_lev%ulevel%restrict(f_lev, c_lev, f_encap_array(m), f_lev%f_encap_array_c(m), f_time(m), flags)
       !print *,'                                      after spatial restrict,f_lev%f_encap_array_c(m)'
       !call f_lev%f_encap_array_c(m)%eprint()
    end do

    ! temporal restriction, rmat is calculated in pf_time_interpolation_matrix in pf_pfasst.f90
    if (present(flags)) then
       if ((flags .eq. 0) .or. (flags .eq. 1)) &
            call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev%rmat, f_lev%f_encap_array_c, .true., flags)
       if ((flags .eq. 0) .or. (flags .eq. 2)) &
            call pf_apply_mat_backward(c_encap_array, 1.0_pfdp, f_lev%rmat, f_lev%f_encap_array_c, .true., flags=2)
    else
       call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev%rmat, f_lev%f_encap_array_c, .true.)
    end if

  end subroutine restrict_ts

  subroutine restrict_ts_integral(f_lev, c_lev, f_encap_array, c_encap_array,f_time, flags)

    !! Restrict (in time and space) f_sol_array  to c_sol_array
    !! This version is for integrals
    
    class(pf_level_t),  intent(inout) :: f_lev   !!   pointer to fine level
    class(pf_level_t),  intent(inout) :: c_lev   !!   pointer to coarse level
    class(pf_encap_t),  intent(inout) :: f_encap_array(:)   !! array of fine level data to be restricted
    class(pf_encap_t),  intent(inout) :: c_encap_array(:)   !! array of coarse level data to be computed
    real(pfdp),         intent(in) :: f_time(:)             !! time at the fine nodes
    integer, optional, intent(in)    :: flags    

    class(pf_encap_t), allocatable :: f_encap_array_c(:)  !!  fine solution restricted in space only
    integer :: m,j
    integer :: f_nnodes,c_nnodes


    f_nnodes = f_lev%nnodes
    c_nnodes = c_lev%nnodes

    !!  Create a temp array for the spatial restriction
    if (f_lev%restrict_workspace_allocated   .eqv. .false.) then      
       call c_lev%ulevel%factory%create_array(f_lev%f_encap_array_c, f_nnodes, c_lev%index, c_lev%lev_shape)
       f_lev%restrict_workspace_allocated  = .true.
    end if

    !  spatial restriction, this restrict function is at user level in bisicle folder
    do m = 1, f_nnodes-1
       call f_lev%ulevel%restrict(f_lev, c_lev, f_encap_array(m), f_lev%f_encap_array_c(m), f_time(m), flags)
       !print *, '                    in restrict_ts_spatial after copy f_ecanp_array, f_lev%f_encap_array_c(m)'
       !call f_lev%f_encap_array_c(m)%eprint()
    end do
    
    ! temporal restriction
    ! when restricting '0 to node' integral terms, skip the first entry since it is zero
    if (present(flags)) then
       if ((flags .eq. 0) .or. (flags .eq. 1)) &
            call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev%rmat(2:,2:), f_lev%f_encap_array_c, .true., flags=1)
       if ((flags .eq. 0) .or. (flags .eq. 2)) &
            call pf_apply_mat_backward(c_encap_array, 1.0_pfdp, f_lev%rmat(2:,2:), f_lev%f_encap_array_c, .true., flags=2)
    else
       call pf_apply_mat(c_encap_array, 1.0_pfdp, f_lev%rmat(2:,2:), f_lev%f_encap_array_c, .true.)
    end if

  end subroutine restrict_ts_integral


end module pf_mod_restrict

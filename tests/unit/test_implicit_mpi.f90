program test_implicit_mpi
    use mpi
    use flow_data
    use mpi_wrapper
    use mod_ImplicitState
    use mod_ImplicitMicro, only: SolveMicroLUSGS
    use mod_ImplicitDiagnostics, only: SelectGlobalAcceptedStateAdmissibilityMaxima, &
        CollectImplicitTimingMaxima
    use mod_SteadyImplicit, only: SelectGlobalImplicitFailure
    use TestSupport
    implicit none
    integer :: ierr_local, rank_local, size_local
    type(ImplicitFailureRecord) :: local_failure,global_failure
    real(KREAL), allocatable :: delta_w_local(:,:,:)
    real(KREAL), allocatable :: delta_h_local(:,:,:), delta_b_local(:,:,:)
    real(KREAL) :: local_metric_values(5),global_metric_values(5)
    integer(KINT) :: local_metric_i(5),local_metric_j(5)
    integer(KINT) :: global_metric_i(5),global_metric_j(5),global_metric_owner(5)
    real(KREAL) :: timing_global(8)
    real(KREAL), target :: h_state(2),b_state(2)
    real(KREAL) :: expected_h(2),expected_b(2)
    logical :: solve_ok

    call MPI_Init(ierr_local)
    call MPI_Comm_rank(MPI_COMM_WORLD,rank_local,ierr_local)
    call MPI_Comm_size(MPI_COMM_WORLD,size_local,ierr_local)
    call AssertTrue(size_local==2,'implicit MPI test requires two ranks')

    myID=rank_local; nx=2; ny=1; Nx_all=4; Ny_all=1; Nvel=3
    ID_Y_minus=MPI_PROC_NULL; ID_Y_plus=MPI_PROC_NULL
    if (rank_local==0) then
        ID_X_minus=MPI_PROC_NULL; ID_X_plus=1
    else
        ID_X_minus=0; ID_X_plus=MPI_PROC_NULL
    endif

    allocate(delta_w_local(4,1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(delta_h_local(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(delta_b_local(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
    delta_w_local=-1.0_KREAL; delta_h_local=-1.0_KREAL; delta_b_local=-1.0_KREAL
    delta_w_local(:,1:nx,1:ny)=real(rank_local+1,KREAL)
    delta_h_local(:,1:nx,1:ny)=real(rank_local+1,KREAL)
    delta_b_local(:,1:nx,1:ny)=real(rank_local+1,KREAL)

    call MPIExchangeDeltaMacro(delta_w_local)
    call MPIExchangeDeltaMicro(delta_h_local,delta_b_local)
    if (rank_local==0) then
        call AssertTrue(all(delta_w_local(:,nx+1,1)==2.0_KREAL),'macro right halo')
        call AssertTrue(all(delta_h_local(:,nx+1,1)==2.0_KREAL),'micro right halo')
    else
        call AssertTrue(all(delta_w_local(:,0,1)==1.0_KREAL),'macro left halo')
        call AssertTrue(all(delta_b_local(:,0,1)==1.0_KREAL),'micro left halo')
    endif

    implicit_pr_seconds=real(rank_local+1,KREAL)
    implicit_lusgs_seconds=real(rank_local+1,KREAL)
    implicit_exchange_seconds=real(rank_local+1,KREAL)
    implicit_evolution_seconds=real(rank_local+1,KREAL)
    implicit_micro_forward_seconds=real(rank_local+1,KREAL)
    implicit_micro_backward_seconds=real(rank_local+1,KREAL)
    implicit_micro_forward_backward_exchange_seconds=real(rank_local+1,KREAL)
    call CollectImplicitTimingMaxima(real(rank_local+1,KREAL),timing_global)
    call AssertTrue(all(timing_global==2.0_KREAL), &
        'implicit timing output selects maximum rank time')

    if (rank_local==0) then
        local_metric_values=[1.0_KREAL,9.0_KREAL,6.0_KREAL,8.0_KREAL,5.0_KREAL]
        local_metric_i=[1,3,3,4,2]
        local_metric_j=[1,1,1,1,1]
    else
        local_metric_values=[7.0_KREAL,4.0_KREAL,6.0_KREAL,2.0_KREAL,10.0_KREAL]
        local_metric_i=[2,4,2,1,3]
        local_metric_j=[1,1,2,2,1]
    endif
    call SelectGlobalAcceptedStateAdmissibilityMaxima(local_metric_values,local_metric_i,local_metric_j, &
        global_metric_values,global_metric_i,global_metric_j,global_metric_owner)
    call AssertTrue(all(global_metric_values==[7.0_KREAL,9.0_KREAL,6.0_KREAL,8.0_KREAL,10.0_KREAL]), &
        'accepted-state diagnostic global maxima')
    call AssertTrue(all(global_metric_i==[2,3,3,4,3]) .and. &
        all(global_metric_j==[1,1,1,1,1]) .and. &
        all(global_metric_owner==[1,0,0,0,1]), &
        'accepted-state diagnostic owner and location transfer')

    call ResetImplicitFailureRecord(local_failure)
    local_failure%valid=.true.
    local_failure%reason=IMPLICIT_FAILURE_H_NEGATIVE_FRACTION
    local_failure%owner_rank=rank_local
    local_failure%local_i=1
    local_failure%local_j=1
    if (rank_local==0) then
        local_failure%global_i=4
        local_failure%global_j=1
        local_failure%candidate_fraction=4.0e-4_KREAL
    else
        local_failure%global_i=2
        local_failure%global_j=1
        local_failure%candidate_fraction=2.0e-4_KREAL
    endif
    local_failure%velocity_k=2
    local_failure%velocity_u=real(rank_local,KREAL)
    local_failure%velocity_v=-real(rank_local,KREAL)
    local_failure%ncfl_attempt=3.90625_KREAL
    local_failure%alpha_tested=2.0_KREAL**(-19)
    local_failure%alpha_next=2.0_KREAL**(-20)
    local_failure%micro_breakdown_valid=.true.
    local_failure%micro_current=10.0_KREAL+real(rank_local,KREAL)
    local_failure%micro_delta=11.0_KREAL+real(rank_local,KREAL)
    local_failure%micro_candidate=12.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_value=13.0_KREAL+real(rank_local,KREAL)
    local_failure%stream_rhs=14.0_KREAL+real(rank_local,KREAL)
    local_failure%collision_rhs=15.0_KREAL+real(rank_local,KREAL)
    local_failure%neighbor_rhs=16.0_KREAL+real(rank_local,KREAL)
    local_failure%dt_pseudo=17.0_KREAL+real(rank_local,KREAL)
    local_failure%tau_value=18.0_KREAL+real(rank_local,KREAL)
    local_failure%diagonal_total=19.0_KREAL+real(rank_local,KREAL)
    local_failure%rhs_total=20.0_KREAL+real(rank_local,KREAL)
    local_failure%post_sweep_delta=21.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_breakdown_valid=.false.
    local_failure%equilibrium_maxwell_value=30.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_shakhov_correction=31.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_preprojection_value=32.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_projection_correction=33.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_final_value=34.0_KREAL+real(rank_local,KREAL)
    local_failure%q_frozen_x=35.0_KREAL+real(rank_local,KREAL)
    local_failure%q_frozen_y=36.0_KREAL+real(rank_local,KREAL)
    local_failure%prim_tilde_rho=37.0_KREAL+real(rank_local,KREAL)
    local_failure%prim_tilde_u=38.0_KREAL+real(rank_local,KREAL)
    local_failure%prim_tilde_v=39.0_KREAL+real(rank_local,KREAL)
    local_failure%prim_tilde_temperature=40.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_maxwell_negative_fraction=41.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_preprojection_negative_fraction=42.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_final_negative_fraction=43.0_KREAL+real(rank_local,KREAL)
    local_failure%equilibrium_breakdown_stage= &
        IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH
    local_failure%reconstruction_max_abs_h=44.0_KREAL+real(rank_local,KREAL)
    local_failure%reconstruction_max_abs_b=45.0_KREAL+real(rank_local,KREAL)
    local_failure%reconstruction_max_rel_h=46.0_KREAL+real(rank_local,KREAL)
    local_failure%reconstruction_max_rel_b=47.0_KREAL+real(rank_local,KREAL)
    local_failure%reconstruction_max_abs_h_k=5+rank_local
    local_failure%reconstruction_max_abs_b_k=6+rank_local
    local_failure%reconstructed_h_at_failure_k=48.0_KREAL+real(rank_local,KREAL)
    local_failure%stored_h_at_failure_k=49.0_KREAL+real(rank_local,KREAL)
    local_failure%reconstructed_b_at_failure_k=50.0_KREAL+real(rank_local,KREAL)
    local_failure%stored_b_at_failure_k=51.0_KREAL+real(rank_local,KREAL)
    call SelectGlobalImplicitFailure(local_failure,global_failure)
    call AssertTrue(global_failure%valid,'global failure record valid')
    call AssertTrue(global_failure%owner_rank==1,'lowest global cell owner selected')
    call AssertTrue(global_failure%global_i==2 .and. global_failure%global_j==1, &
        'lowest global cell selected')
    call AssertClose(global_failure%candidate_fraction,2.0e-4_KREAL, &
        epsilon(1.0_KREAL),'selected failure payload broadcast')
    call AssertTrue(global_failure%micro_breakdown_valid,'micro breakdown validity broadcast')
    call AssertClose(global_failure%micro_current,11.0_KREAL,epsilon(1.0_KREAL), &
        'micro current payload broadcast')
    call AssertClose(global_failure%micro_delta,12.0_KREAL,epsilon(1.0_KREAL), &
        'micro delta payload broadcast')
    call AssertClose(global_failure%micro_candidate,13.0_KREAL,epsilon(1.0_KREAL), &
        'micro candidate payload broadcast')
    call AssertClose(global_failure%equilibrium_value,14.0_KREAL,epsilon(1.0_KREAL), &
        'equilibrium payload broadcast')
    call AssertClose(global_failure%stream_rhs,15.0_KREAL,epsilon(1.0_KREAL), &
        'stream rhs payload broadcast')
    call AssertClose(global_failure%collision_rhs,16.0_KREAL,epsilon(1.0_KREAL), &
        'collision rhs payload broadcast')
    call AssertClose(global_failure%neighbor_rhs,17.0_KREAL,epsilon(1.0_KREAL), &
        'neighbor rhs payload broadcast')
    call AssertClose(global_failure%dt_pseudo,18.0_KREAL,epsilon(1.0_KREAL), &
        'pseudo time step payload broadcast')
    call AssertClose(global_failure%tau_value,19.0_KREAL,epsilon(1.0_KREAL), &
        'tau payload broadcast')
    call AssertClose(global_failure%diagonal_total,20.0_KREAL,epsilon(1.0_KREAL), &
        'diagonal payload broadcast')
    call AssertClose(global_failure%rhs_total,21.0_KREAL,epsilon(1.0_KREAL), &
        'rhs total payload broadcast')
    call AssertClose(global_failure%post_sweep_delta,22.0_KREAL,epsilon(1.0_KREAL), &
        'post-sweep delta payload broadcast')
    call AssertTrue(.not.global_failure%equilibrium_breakdown_valid .and. &
        global_failure%equilibrium_breakdown_stage== &
        IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH, &
        'equilibrium mismatch stage and invalidity broadcast')
    call AssertClose(global_failure%equilibrium_maxwell_value,31.0_KREAL,epsilon(1.0_KREAL), &
        'equilibrium Maxwell payload broadcast')
    call AssertClose(global_failure%equilibrium_shakhov_correction,32.0_KREAL,epsilon(1.0_KREAL), &
        'equilibrium Shakhov payload broadcast')
    call AssertClose(global_failure%equilibrium_preprojection_value,33.0_KREAL,epsilon(1.0_KREAL), &
        'equilibrium preprojection payload broadcast')
    call AssertClose(global_failure%equilibrium_projection_correction,34.0_KREAL,epsilon(1.0_KREAL), &
        'equilibrium projection payload broadcast')
    call AssertClose(global_failure%equilibrium_final_value,35.0_KREAL,epsilon(1.0_KREAL), &
        'equilibrium final payload broadcast')
    call AssertClose(global_failure%q_frozen_x,36.0_KREAL,epsilon(1.0_KREAL), &
        'frozen qx payload broadcast')
    call AssertClose(global_failure%q_frozen_y,37.0_KREAL,epsilon(1.0_KREAL), &
        'frozen qy payload broadcast')
    call AssertClose(global_failure%prim_tilde_rho,38.0_KREAL,epsilon(1.0_KREAL), &
        'predicted rho payload broadcast')
    call AssertClose(global_failure%prim_tilde_u,39.0_KREAL,epsilon(1.0_KREAL), &
        'predicted u payload broadcast')
    call AssertClose(global_failure%prim_tilde_v,40.0_KREAL,epsilon(1.0_KREAL), &
        'predicted v payload broadcast')
    call AssertClose(global_failure%prim_tilde_temperature,41.0_KREAL,epsilon(1.0_KREAL), &
        'predicted temperature payload broadcast')
    call AssertClose(global_failure%equilibrium_maxwell_negative_fraction,42.0_KREAL, &
        epsilon(1.0_KREAL),'Maxwell negative fraction payload broadcast')
    call AssertClose(global_failure%equilibrium_preprojection_negative_fraction,43.0_KREAL, &
        epsilon(1.0_KREAL),'preprojection negative fraction payload broadcast')
    call AssertClose(global_failure%equilibrium_final_negative_fraction,44.0_KREAL, &
        epsilon(1.0_KREAL),'final negative fraction payload broadcast')
    call AssertClose(global_failure%reconstruction_max_abs_h,45.0_KREAL, &
        epsilon(1.0_KREAL),'h reconstruction absolute error payload broadcast')
    call AssertClose(global_failure%reconstruction_max_abs_b,46.0_KREAL, &
        epsilon(1.0_KREAL),'b reconstruction absolute error payload broadcast')
    call AssertClose(global_failure%reconstruction_max_rel_h,47.0_KREAL, &
        epsilon(1.0_KREAL),'h reconstruction relative error payload broadcast')
    call AssertClose(global_failure%reconstruction_max_rel_b,48.0_KREAL, &
        epsilon(1.0_KREAL),'b reconstruction relative error payload broadcast')
    call AssertTrue(global_failure%reconstruction_max_abs_h_k==6 .and. &
        global_failure%reconstruction_max_abs_b_k==7, &
        'reconstruction error locations payload broadcast')
    call AssertClose(global_failure%reconstructed_h_at_failure_k,49.0_KREAL, &
        epsilon(1.0_KREAL),'reconstructed h payload broadcast')
    call AssertClose(global_failure%stored_h_at_failure_k,50.0_KREAL, &
        epsilon(1.0_KREAL),'stored h payload broadcast')
    call AssertClose(global_failure%reconstructed_b_at_failure_k,51.0_KREAL, &
        epsilon(1.0_KREAL),'reconstructed b payload broadcast')
    call AssertClose(global_failure%stored_b_at_failure_k,52.0_KREAL, &
        epsilon(1.0_KREAL),'stored b payload broadcast')

    nx=1; ny=1; Nx_all=2; Ny_all=1; Nu=2; Nv=1; Nvel=2
    CFL=0.5_KREAL; Implicit_nCFL=10.0_KREAL; Implicit_MicroSweeps=1
    allocate(uVel(Nvel),vVel(Nvel),wVel(Nvel))
    allocate(ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    uVel=[1.0_KREAL,-1.0_KREAL]; vVel=0.0_KREAL; wVel=0.5_KREAL
    h_state=0.0_KREAL; b_state=0.0_KREAL
    ctr(1,1)%h_vec=>h_state; ctr(1,1)%b_vec=>b_state
    ctr(1,1)%area=1.0_KREAL; ctr(1,1)%dt_local=0.1_KREAL
    vface%length=1.0_KREAL; vface%cosx=1.0_KREAL; vface%cosy=0.0_KREAL
    hface%length=0.0_KREAL; hface%cosx=0.0_KREAL; hface%cosy=1.0_KREAL
    BC_Type=BCconst_freestream
    ID_Y_minus=MPI_PROC_NULL; ID_Y_plus=MPI_PROC_NULL
    if (rank_local==0) then
        ID_X_minus=MPI_PROC_NULL; ID_X_plus=1
        expected_h=[0.4_KREAL,0.56_KREAL]
        expected_b=[0.8_KREAL,1.12_KREAL]
    else
        ID_X_minus=0; ID_X_plus=MPI_PROC_NULL
        expected_h=[0.56_KREAL,0.4_KREAL]
        expected_b=[1.12_KREAL,0.8_KREAL]
    endif
    call AllocateImplicitState()
    implicit_ncfl_active=Implicit_nCFL
    delta_h=0.0_KREAL; delta_b=0.0_KREAL
    stream_h=1.0_KREAL; stream_b=2.0_KREAL
    equilibrium_h=0.0_KREAL; equilibrium_b=0.0_KREAL; tau_tilde=1.0_KREAL
    implicit_micro_exchange_calls=0
    call SolveMicroLUSGS(solve_ok)
    call AssertTrue(solve_ok .and. implicit_micro_exchange_calls==1, &
        'two-rank one-sweep micro solve uses one exchange')
    call AssertTrue(maxval(abs(delta_h(:,1,1)-expected_h))<1.0e-14_KREAL .and. &
        maxval(abs(delta_b(:,1,1)-expected_b))<1.0e-14_KREAL, &
        'two-rank forward-backward dependency matches redundant-exchange result')
    call ResetImplicitCorrections()
    implicit_micro_exchange_calls=0
    Implicit_MicroSweeps=2
    call SolveMicroLUSGS(solve_ok)
    call AssertTrue(solve_ok .and. implicit_micro_exchange_calls==3, &
        'two-rank two-sweep micro solve uses three exchanges')

    call FinishTests('implicit_mpi')
    call MPI_Finalize(ierr_local)
end program test_implicit_mpi

program test_implicit_micro
    use mpi
    use flow_data
    use Tools
    use mod_ImplicitState
    use mod_ImplicitBoundary, only: BuildStructuredMirrorMaps, ComputeDiffuseWallFaceDelta, &
        ComputeMaxwellWallFaceDelta
    use mod_ImplicitMicro
    use mod_ImplicitDiagnostics, only: AcceptedStateAdmissibilityDiagnostics, &
        CollectAcceptedStateAdmissibilityDiagnostics, InitializeImplicitDiagnostics
    use TestSupport
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer :: ierr_local
    logical :: ok
    type(ImplicitFailureRecord) :: failure
    real(KREAL), target :: h_store(4),b_store(4)
    real(KREAL) :: prim_test(4),delta_value,expected,norm_h_before,norm_b_before,norm_h_after,norm_b_after
    real(KREAL) :: maxwell_h(4),maxwell_b(4),h_plus(4),b_plus(4)
    real(KREAL) :: h_preprojection(4),b_preprojection(4),prim_tilde_check(4)
    real(KREAL) :: un,face_length,neighbor_h,neighbor_b,rho_delta
    real(KREAL) :: diffuse_face_h(4),diffuse_face_b(4)
    real(KREAL) :: frozen_source_h(4),frozen_source_b(4),frozen_face_h(4),frozen_face_b(4)
    real(KREAL) :: accommodations(3)
    real(KREAL) :: wall_normal_x(4),wall_normal_y(4)
    real(KREAL) :: current_value,stored_delta,equilibrium_value,stream_rhs,collision_rhs
    real(KREAL) :: neighbor_rhs,dt_pseudo,tau_value,diagonal_total,rhs_total,post_sweep_delta
    real(KREAL) :: delta_h_before(4),delta_b_before(4)
    real(KREAL) :: w_tilde_before(4),q_frozen_before(2),equilibrium_h_before(4),equilibrium_b_before(4)
    real(KREAL) :: h_store_before(4),b_store_before(4),tau_tilde_before
    real(KREAL) :: collector_convars_before(4),collector_h_before(4),collector_b_before(4)
    real(KREAL) :: equilibrium_maxwell_value,equilibrium_shakhov_correction
    real(KREAL) :: equilibrium_preprojection_value,equilibrium_projection_correction
    real(KREAL) :: equilibrium_final_value,q_frozen_x,q_frozen_y
    real(KREAL) :: prim_tilde_rho,prim_tilde_u,prim_tilde_v,prim_tilde_temperature
    real(KREAL) :: equilibrium_maxwell_negative_fraction
    real(KREAL) :: equilibrium_preprojection_negative_fraction
    real(KREAL) :: equilibrium_final_negative_fraction
    real(KREAL) :: reconstruction_max_abs_h,reconstruction_max_abs_b
    real(KREAL) :: reconstruction_max_rel_h,reconstruction_max_rel_b
    real(KREAL) :: reconstructed_h_at_failure_k,stored_h_at_failure_k
    real(KREAL) :: reconstructed_b_at_failure_k,stored_b_at_failure_k
    integer(KINT) :: reconstruction_max_abs_h_k,reconstruction_max_abs_b_k
    integer(KINT) :: equilibrium_breakdown_stage
    type(AcceptedStateAdmissibilityDiagnostics) :: accepted_diagnostic
    real(KREAL) :: explicit_q(2),explicit_equilibrium_h(4),explicit_equilibrium_b(4)
    real(KREAL) :: explicit_expected(5)
    real(KREAL) :: micro_timing_before(5)
    integer(KINT) :: accommodation_index,k_test,face_test,wall_side(4)
    type(FrozenMicroWallState) :: frozen_wall_state

    call MPI_Init(ierr_local)
    call MPI_Comm_rank(MPI_COMM_WORLD,myID,ierr_local)
    failure%micro_breakdown_valid=.true.
    failure%equilibrium_breakdown_valid=.true.
    failure%micro_current=1.0_KREAL; failure%micro_delta=2.0_KREAL
    failure%micro_candidate=3.0_KREAL; failure%equilibrium_value=4.0_KREAL
    failure%stream_rhs=5.0_KREAL; failure%collision_rhs=6.0_KREAL
    failure%neighbor_rhs=7.0_KREAL; failure%dt_pseudo=8.0_KREAL
    failure%tau_value=9.0_KREAL; failure%diagonal_total=10.0_KREAL
    failure%rhs_total=11.0_KREAL; failure%post_sweep_delta=12.0_KREAL
    failure%equilibrium_maxwell_value=13.0_KREAL
    failure%equilibrium_shakhov_correction=14.0_KREAL
    failure%equilibrium_preprojection_value=15.0_KREAL
    failure%equilibrium_projection_correction=16.0_KREAL
    failure%equilibrium_final_value=17.0_KREAL
    failure%q_frozen_x=18.0_KREAL; failure%q_frozen_y=19.0_KREAL
    failure%prim_tilde_rho=20.0_KREAL; failure%prim_tilde_u=21.0_KREAL
    failure%prim_tilde_v=22.0_KREAL; failure%prim_tilde_temperature=23.0_KREAL
    failure%equilibrium_maxwell_negative_fraction=24.0_KREAL
    failure%equilibrium_preprojection_negative_fraction=25.0_KREAL
    failure%equilibrium_final_negative_fraction=26.0_KREAL
    failure%equilibrium_breakdown_stage=99
    failure%reconstruction_max_abs_h=27.0_KREAL
    failure%reconstruction_max_abs_b=28.0_KREAL
    failure%reconstruction_max_rel_h=29.0_KREAL
    failure%reconstruction_max_rel_b=30.0_KREAL
    failure%reconstruction_max_abs_h_k=31
    failure%reconstruction_max_abs_b_k=32
    failure%reconstructed_h_at_failure_k=33.0_KREAL
    failure%stored_h_at_failure_k=34.0_KREAL
    failure%reconstructed_b_at_failure_k=35.0_KREAL
    failure%stored_b_at_failure_k=36.0_KREAL
    call ResetImplicitFailureRecord(failure)
    call AssertTrue(.not.failure%micro_breakdown_valid .and. .not.failure%equilibrium_breakdown_valid .and. &
        all([failure%micro_current, &
        failure%micro_delta,failure%micro_candidate,failure%equilibrium_value,failure%stream_rhs, &
        failure%collision_rhs,failure%neighbor_rhs,failure%dt_pseudo,failure%tau_value, &
        failure%diagonal_total,failure%rhs_total,failure%post_sweep_delta, &
        failure%equilibrium_maxwell_value,failure%equilibrium_shakhov_correction, &
        failure%equilibrium_preprojection_value,failure%equilibrium_projection_correction, &
        failure%equilibrium_final_value,failure%q_frozen_x,failure%q_frozen_y, &
        failure%prim_tilde_rho,failure%prim_tilde_u,failure%prim_tilde_v, &
        failure%prim_tilde_temperature,failure%equilibrium_maxwell_negative_fraction, &
        failure%equilibrium_preprojection_negative_fraction, &
        failure%equilibrium_final_negative_fraction,failure%reconstruction_max_abs_h, &
        failure%reconstruction_max_abs_b,failure%reconstruction_max_rel_h, &
        failure%reconstruction_max_rel_b,failure%reconstructed_h_at_failure_k, &
        failure%stored_h_at_failure_k,failure%reconstructed_b_at_failure_k, &
        failure%stored_b_at_failure_k]<=0.0_KREAL) .and. &
        failure%equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE .and. &
        failure%reconstruction_max_abs_h_k==0 .and. failure%reconstruction_max_abs_b_k==0, &
        'micro and equilibrium breakdown reset clears all fields')
    implicit_failure=failure
    implicit_failure%equilibrium_breakdown_valid=.true.
    implicit_failure%equilibrium_final_value=1.0_KREAL
    call ResetImplicitFailureRecord()
    call AssertTrue(.not.implicit_failure%equilibrium_breakdown_valid .and. &
        implicit_failure%equilibrium_final_value==0.0_KREAL, &
        'global equilibrium breakdown reset clears validity and values')
    nx=1; ny=1; Nx_all=1; Ny_all=1; Nu=2; Nv=2; Nvel=4; ck=1; gamma=5.0_KREAL/3.0_KREAL
    CFL=0.5_KREAL; Implicit_nCFL=10.0_KREAL; Implicit_MicroSweeps=1
    flag_vislaw=PowerLaw; mu_ref=0.05_KREAL; omega=0.5_KREAL
    allocate(uVel(4),vVel(4),wVel(4),ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    uVel=[-1.0_KREAL,1.0_KREAL,-1.0_KREAL,1.0_KREAL]
    vVel=[-1.0_KREAL,-1.0_KREAL,1.0_KREAL,1.0_KREAL]; wVel=0.25_KREAL
    h_store=0.2_KREAL; b_store=0.1_KREAL
    ctr(1,1)%h_vec=>h_store; ctr(1,1)%b_vec=>b_store
    ctr(1,1)%area=1.0_KREAL; ctr(1,1)%dt_local=0.1_KREAL
    vface%length=0.0_KREAL; vface%cosx=1.0_KREAL; vface%cosy=0.0_KREAL
    hface%length=0.0_KREAL; hface%cosx=0.0_KREAL; hface%cosy=1.0_KREAL
    ID_X_minus=MPI_PROC_NULL; ID_X_plus=MPI_PROC_NULL
    ID_Y_minus=MPI_PROC_NULL; ID_Y_plus=MPI_PROC_NULL; BC_Type=BCconst_freestream
    call AllocateImplicitState(); implicit_ncfl_active=Implicit_nCFL
    implicit_micro_forward_seconds=1.0_KREAL
    implicit_micro_backward_seconds=2.0_KREAL
    implicit_micro_initial_exchange_seconds=3.0_KREAL
    implicit_micro_forward_backward_exchange_seconds=4.0_KREAL
    implicit_micro_post_backward_exchange_seconds=5.0_KREAL
    implicit_micro_exchange_calls=6
    call InitializeImplicitDiagnostics()
    call AssertTrue(all([implicit_micro_forward_seconds,implicit_micro_backward_seconds, &
        implicit_micro_initial_exchange_seconds,implicit_micro_forward_backward_exchange_seconds, &
        implicit_micro_post_backward_exchange_seconds]==0.0_KREAL) .and. &
        implicit_micro_exchange_calls==0,'micro LU-SGS timing counters initialize to zero')
    call BuildStructuredMirrorMaps()
    BG_geometry=BGconst_general
    Prim_wall=[1.0_KREAL,0.0_KREAL,0.0_KREAL,1.0_KREAL]
    BC_Type(3)=BCconst_wall_isothermal
    FluxWall_type=FluxWall_Maxwell
    ACC=1.0_KREAL
    delta_h(:,1,1)=[0.4_KREAL,-0.2_KREAL,0.1_KREAL,0.3_KREAL]
    delta_b(:,1,1)=[0.2_KREAL,-0.1_KREAL,0.05_KREAL,0.15_KREAL]
    prim_test=[1.0_KREAL,0.2_KREAL,-0.1_KREAL,1.0_KREAL]
    ctr(1,1)%conVars=GetConserved(prim_test)
    w_tilde(:,1,1)=GetConserved(prim_test)
    call SolveCollisionOnlyIncrement(0.5_KREAL,0.25_KREAL,2.0_KREAL,3.0_KREAL,delta_value)
    expected=((3.0_KREAL-2.0_KREAL)/0.25_KREAL)/(1.0_KREAL/0.5_KREAL+1.0_KREAL/0.25_KREAL)
    call AssertClose(delta_value,expected,1.0e-14_KREAL,'collision-only implicit increment')
    Pr=1.0_KREAL
    call DiscreteMaxwell_Vec(maxwell_h,maxwell_b,prim_test)
    call ShakhovPart_Vec(maxwell_h,maxwell_b,[0.3_KREAL,-0.2_KREAL],prim_test,h_plus,b_plus)
    call AssertTrue(maxval(abs(h_plus))<1.0e-13_KREAL,'Pr=1 BGK h limit')
    call AssertTrue(maxval(abs(b_plus))<1.0e-13_KREAL,'Pr=1 BGK b limit')
    call BuildPredictedEquilibrium(ok)
    call AssertTrue(ok .and. all(ieee_is_finite(equilibrium_h)),'predicted equilibrium finite')
    call AssertTrue(ok .and. all(ieee_is_finite(equilibrium_b)),'predicted b equilibrium finite')
    prim_tilde_check=GetPrimary(w_tilde(:,1,1))
    call DiscreteMaxwell_Vec(maxwell_h,maxwell_b,prim_tilde_check)
    call ShakhovPart_Vec(maxwell_h,maxwell_b,q_frozen(:,1,1),prim_tilde_check,h_plus,b_plus)
    h_preprojection=maxwell_h+h_plus
    b_preprojection=maxwell_b+b_plus
    delta_h_before=delta_h(:,1,1); delta_b_before=delta_b(:,1,1)
    w_tilde_before=w_tilde(:,1,1); q_frozen_before=q_frozen(:,1,1)
    equilibrium_h_before=equilibrium_h(:,1,1); equilibrium_b_before=equilibrium_b(:,1,1)
    h_store_before=h_store; b_store_before=b_store; tau_tilde_before=tau_tilde(1,1)
    call EvaluatePredictedEquilibriumBreakdown(1,1,1,.true., &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
        prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction, &
        reconstruction_max_abs_h,reconstruction_max_abs_b,reconstruction_max_rel_h, &
        reconstruction_max_rel_b,reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k,reconstructed_b_at_failure_k, &
        stored_b_at_failure_k,equilibrium_breakdown_stage,ok)
    call AssertTrue(ok,'Pr=1 h predicted-equilibrium breakdown succeeds')
    call AssertTrue(equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS, &
        'Pr=1 h exact reconstruction stage')
    call AssertClose(equilibrium_shakhov_correction,0.0_KREAL,100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 h Shakhov correction vanishes')
    call AssertClose(equilibrium_preprojection_value,equilibrium_maxwell_value+ &
        equilibrium_shakhov_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 h preprojection sum identity')
    call AssertClose(equilibrium_final_value,equilibrium_preprojection_value+ &
        equilibrium_projection_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 h final sum identity')
    call AssertClose(equilibrium_final_value,equilibrium_h_before(1), &
        100.0_KREAL*epsilon(1.0_KREAL),'Pr=1 h final saved point value')
    call AssertClose(equilibrium_maxwell_negative_fraction,WeightedNegativeFraction(maxwell_h,wVel), &
        100.0_KREAL*epsilon(1.0_KREAL),'Pr=1 h Maxwellian negative fraction')
    call AssertClose(equilibrium_preprojection_negative_fraction, &
        WeightedNegativeFraction(h_preprojection,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 h preprojection negative fraction')
    call AssertClose(equilibrium_final_negative_fraction, &
        WeightedNegativeFraction(equilibrium_h_before,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 h final negative fraction')
    call EvaluatePredictedEquilibriumBreakdown(1,1,2,.false., &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
        prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction, &
        reconstruction_max_abs_h,reconstruction_max_abs_b,reconstruction_max_rel_h, &
        reconstruction_max_rel_b,reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k,reconstructed_b_at_failure_k, &
        stored_b_at_failure_k,equilibrium_breakdown_stage,ok)
    call AssertTrue(ok,'Pr=1 b predicted-equilibrium breakdown succeeds')
    call AssertTrue(equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS, &
        'Pr=1 b exact reconstruction stage')
    call AssertClose(equilibrium_shakhov_correction,0.0_KREAL,100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 b Shakhov correction vanishes')
    call AssertClose(equilibrium_preprojection_value,equilibrium_maxwell_value+ &
        equilibrium_shakhov_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 b preprojection sum identity')
    call AssertClose(equilibrium_final_value,equilibrium_preprojection_value+ &
        equilibrium_projection_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 b final sum identity')
    call AssertClose(equilibrium_final_value,equilibrium_b_before(2), &
        100.0_KREAL*epsilon(1.0_KREAL),'Pr=1 b final saved point value')
    call AssertClose(equilibrium_maxwell_negative_fraction,WeightedNegativeFraction(maxwell_b,wVel), &
        100.0_KREAL*epsilon(1.0_KREAL),'Pr=1 b Maxwellian negative fraction')
    call AssertClose(equilibrium_preprojection_negative_fraction, &
        WeightedNegativeFraction(b_preprojection,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 b preprojection negative fraction')
    call AssertClose(equilibrium_final_negative_fraction, &
        WeightedNegativeFraction(equilibrium_b_before,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'Pr=1 b final negative fraction')
    call AssertTrue(maxval(abs(w_tilde(:,1,1)-w_tilde_before))<=0.0_KREAL .and. &
        maxval(abs(q_frozen(:,1,1)-q_frozen_before))<=0.0_KREAL .and. &
        maxval(abs(equilibrium_h(:,1,1)-equilibrium_h_before))<=0.0_KREAL .and. &
        maxval(abs(equilibrium_b(:,1,1)-equilibrium_b_before))<=0.0_KREAL .and. &
        maxval(abs(delta_h(:,1,1)-delta_h_before))<=0.0_KREAL .and. &
        maxval(abs(delta_b(:,1,1)-delta_b_before))<=0.0_KREAL .and. &
        maxval(abs(h_store-h_store_before))<=0.0_KREAL .and. &
        maxval(abs(b_store-b_store_before))<=0.0_KREAL .and. &
        tau_tilde(1,1)==tau_tilde_before,'Pr=1 equilibrium breakdown has no side effects')
    Pr=0.6667_KREAL
    explicit_q=GetHeatFlux_Vec(h_store,b_store,prim_test)
    call DiscreteMaxwell_Vec(maxwell_h,maxwell_b,prim_test)
    call ShakhovPart_Vec(maxwell_h,maxwell_b,explicit_q,prim_test,h_plus,b_plus)
    explicit_equilibrium_h=maxwell_h+h_plus
    explicit_equilibrium_b=maxwell_b+b_plus
    explicit_expected=[WeightedNegativeFraction(h_store,wVel), &
        WeightedNegativeFraction(b_store,wVel), &
        WeightedNegativeFraction(explicit_equilibrium_h,wVel), &
        WeightedNegativeFraction(explicit_equilibrium_b,wVel), &
        sqrt(explicit_q(1)**2+explicit_q(2)**2)]
    collector_convars_before=ctr(1,1)%conVars
    collector_h_before=h_store
    collector_b_before=b_store
    call CollectAcceptedStateAdmissibilityDiagnostics(accepted_diagnostic)
    call AssertTrue(all(ctr(1,1)%conVars==collector_convars_before), &
        'accepted-state collector does not mutate conVars')
    call AssertTrue(all(h_store==collector_h_before), &
        'accepted-state collector does not mutate h_store')
    call AssertTrue(all(b_store==collector_b_before), &
        'accepted-state collector does not mutate b_store')
    call AssertClose(accepted_diagnostic%max_value(1),explicit_expected(1), &
        100.0_KREAL*epsilon(1.0_KREAL),'accepted-state current h negative fraction')
    call AssertClose(accepted_diagnostic%max_value(2),explicit_expected(2), &
        100.0_KREAL*epsilon(1.0_KREAL),'accepted-state current b negative fraction')
    call AssertClose(accepted_diagnostic%max_value(3),explicit_expected(3), &
        100.0_KREAL*epsilon(1.0_KREAL),'accepted-state Shakhov h negative fraction')
    call AssertClose(accepted_diagnostic%max_value(4),explicit_expected(4), &
        100.0_KREAL*epsilon(1.0_KREAL),'accepted-state Shakhov b negative fraction')
    call AssertClose(accepted_diagnostic%max_value(5),explicit_expected(5), &
        100.0_KREAL*epsilon(1.0_KREAL),'accepted-state heat flux magnitude')
    call AssertTrue(all(accepted_diagnostic%global_i==1) .and. &
        all(accepted_diagnostic%global_j==1) .and. &
        all(accepted_diagnostic%owner_rank==0),'accepted-state diagnostic owner and location')
    call BuildPredictedEquilibrium(ok)
    call AssertTrue(ok,'non-BGK predicted equilibrium finite')
    prim_tilde_check=GetPrimary(w_tilde(:,1,1))
    call DiscreteMaxwell_Vec(maxwell_h,maxwell_b,prim_tilde_check)
    call ShakhovPart_Vec(maxwell_h,maxwell_b,q_frozen(:,1,1),prim_tilde_check,h_plus,b_plus)
    h_preprojection=maxwell_h+h_plus
    b_preprojection=maxwell_b+b_plus
    BC_Type=BCconst_freestream
    stream_h=0.0_KREAL; stream_b=0.0_KREAL
    delta_h(:,1,1)=[0.4_KREAL,-0.2_KREAL,0.1_KREAL,0.3_KREAL]
    delta_b(:,1,1)=[0.2_KREAL,-0.1_KREAL,0.05_KREAL,0.15_KREAL]
    delta_h_before=delta_h(:,1,1); delta_b_before=delta_b(:,1,1)
    w_tilde_before=w_tilde(:,1,1); q_frozen_before=q_frozen(:,1,1)
    equilibrium_h_before=equilibrium_h(:,1,1); equilibrium_b_before=equilibrium_b(:,1,1)
    h_store_before=h_store; b_store_before=b_store; tau_tilde_before=tau_tilde(1,1)
    call EvaluatePredictedEquilibriumBreakdown(1,1,1,.true., &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
        prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction, &
        reconstruction_max_abs_h,reconstruction_max_abs_b,reconstruction_max_rel_h, &
        reconstruction_max_rel_b,reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k,reconstructed_b_at_failure_k, &
        stored_b_at_failure_k,equilibrium_breakdown_stage,ok)
    call AssertTrue(ok,'h predicted-equilibrium breakdown succeeds')
    call AssertTrue(equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS, &
        'h exact reconstruction stage')
    call AssertClose(equilibrium_preprojection_value,equilibrium_maxwell_value+ &
        equilibrium_shakhov_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'h preprojection sum identity')
    call AssertClose(equilibrium_maxwell_value,maxwell_h(1),100.0_KREAL*epsilon(1.0_KREAL), &
        'h Maxwellian point value')
    call AssertClose(equilibrium_shakhov_correction,h_plus(1),100.0_KREAL*epsilon(1.0_KREAL), &
        'h Shakhov point value')
    call AssertClose(equilibrium_preprojection_value,h_preprojection(1), &
        100.0_KREAL*epsilon(1.0_KREAL),'h preprojection point value')
    call AssertClose(equilibrium_projection_correction,equilibrium_h(1,1,1)-h_preprojection(1), &
        100.0_KREAL*epsilon(1.0_KREAL),'h projection correction identity')
    call AssertClose(equilibrium_final_value,equilibrium_h(1,1,1), &
        100.0_KREAL*epsilon(1.0_KREAL),'h final point value')
    call AssertClose(equilibrium_maxwell_negative_fraction,WeightedNegativeFraction(maxwell_h,wVel), &
        100.0_KREAL*epsilon(1.0_KREAL),'h Maxwellian negative fraction')
    call AssertClose(equilibrium_preprojection_negative_fraction, &
        WeightedNegativeFraction(h_preprojection,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'h preprojection negative fraction')
    call AssertClose(equilibrium_final_negative_fraction, &
        WeightedNegativeFraction(equilibrium_h(:,1,1),wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'h final negative fraction')
    call AssertClose(q_frozen_x,q_frozen(1,1,1),100.0_KREAL*epsilon(1.0_KREAL),'h frozen qx')
    call AssertClose(q_frozen_y,q_frozen(2,1,1),100.0_KREAL*epsilon(1.0_KREAL),'h frozen qy')
    call AssertClose(prim_tilde_rho,prim_tilde_check(1),100.0_KREAL*epsilon(1.0_KREAL), &
        'h predicted density')
    call AssertClose(prim_tilde_u,prim_tilde_check(2),100.0_KREAL*epsilon(1.0_KREAL), &
        'h predicted x velocity')
    call AssertClose(prim_tilde_v,prim_tilde_check(3),100.0_KREAL*epsilon(1.0_KREAL), &
        'h predicted y velocity')
    call AssertClose(prim_tilde_temperature,1.0_KREAL/prim_tilde_check(4), &
        100.0_KREAL*epsilon(1.0_KREAL),'h predicted temperature')
    call EvaluatePredictedEquilibriumBreakdown(1,1,2,.false., &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
        prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction, &
        reconstruction_max_abs_h,reconstruction_max_abs_b,reconstruction_max_rel_h, &
        reconstruction_max_rel_b,reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k,reconstructed_b_at_failure_k, &
        stored_b_at_failure_k,equilibrium_breakdown_stage,ok)
    call AssertTrue(ok,'b predicted-equilibrium breakdown succeeds')
    call AssertTrue(equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS, &
        'b exact reconstruction stage')
    call AssertClose(equilibrium_preprojection_value,equilibrium_maxwell_value+ &
        equilibrium_shakhov_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'b preprojection sum identity')
    call AssertClose(equilibrium_maxwell_value,maxwell_b(2),100.0_KREAL*epsilon(1.0_KREAL), &
        'b Maxwellian point value')
    call AssertClose(equilibrium_shakhov_correction,b_plus(2),100.0_KREAL*epsilon(1.0_KREAL), &
        'b Shakhov point value')
    call AssertClose(equilibrium_preprojection_value,b_preprojection(2), &
        100.0_KREAL*epsilon(1.0_KREAL),'b preprojection point value')
    call AssertClose(equilibrium_projection_correction,equilibrium_b(2,1,1)-b_preprojection(2), &
        100.0_KREAL*epsilon(1.0_KREAL),'b projection correction identity')
    call AssertClose(equilibrium_final_value,equilibrium_b(2,1,1), &
        100.0_KREAL*epsilon(1.0_KREAL),'b final point value')
    call AssertClose(equilibrium_maxwell_negative_fraction,WeightedNegativeFraction(maxwell_b,wVel), &
        100.0_KREAL*epsilon(1.0_KREAL),'b Maxwellian negative fraction')
    call AssertClose(equilibrium_preprojection_negative_fraction, &
        WeightedNegativeFraction(b_preprojection,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'b preprojection negative fraction')
    call AssertClose(equilibrium_final_negative_fraction, &
        WeightedNegativeFraction(equilibrium_b(:,1,1),wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'b final negative fraction')
    call AssertTrue(maxval(abs(w_tilde(:,1,1)-w_tilde_before))<=0.0_KREAL .and. &
        maxval(abs(q_frozen(:,1,1)-q_frozen_before))<=0.0_KREAL .and. &
        maxval(abs(equilibrium_h(:,1,1)-equilibrium_h_before))<=0.0_KREAL .and. &
        maxval(abs(equilibrium_b(:,1,1)-equilibrium_b_before))<=0.0_KREAL .and. &
        maxval(abs(delta_h(:,1,1)-delta_h_before))<=0.0_KREAL .and. &
        maxval(abs(delta_b(:,1,1)-delta_b_before))<=0.0_KREAL .and. &
        maxval(abs(h_store-h_store_before))<=0.0_KREAL .and. &
        maxval(abs(b_store-b_store_before))<=0.0_KREAL .and. &
        tau_tilde(1,1)==tau_tilde_before,'equilibrium breakdown has no side effects')
    equilibrium_h(1,1,1)=equilibrium_h_before(1)+spacing(equilibrium_h_before(1))
    call EvaluatePredictedEquilibriumBreakdown(1,1,1,.true., &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
        prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction, &
        reconstruction_max_abs_h,reconstruction_max_abs_b,reconstruction_max_rel_h, &
        reconstruction_max_rel_b,reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k,reconstructed_b_at_failure_k, &
        stored_b_at_failure_k,equilibrium_breakdown_stage,ok)
    call AssertTrue(.not.ok .and. equilibrium_breakdown_stage== &
        IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH, &
        'saved h perturbation reports reconstruction mismatch')
    call AssertTrue(reconstruction_max_abs_h>0.0_KREAL .and. reconstruction_max_rel_h>0.0_KREAL .and. &
        reconstruction_max_abs_h_k==1,'h reconstruction mismatch magnitude and location')
    call AssertClose(reconstruction_max_rel_h,reconstruction_max_abs_h/max( &
        abs(stored_h_at_failure_k),tiny(1.0_KREAL)),100.0_KREAL*epsilon(1.0_KREAL), &
        'h reconstruction mismatch relative error')
    call AssertClose(reconstruction_max_abs_b,0.0_KREAL,100.0_KREAL*epsilon(1.0_KREAL), &
        'unperturbed b reconstruction absolute error')
    call AssertClose(reconstruction_max_rel_b,0.0_KREAL,100.0_KREAL*epsilon(1.0_KREAL), &
        'unperturbed b reconstruction relative error')
    call AssertTrue(reconstruction_max_abs_b_k==1,'unperturbed b reconstruction location')
    call AssertClose(reconstruction_max_abs_h,abs(reconstructed_h_at_failure_k- &
        stored_h_at_failure_k),100.0_KREAL*epsilon(1.0_KREAL), &
        'h reconstruction mismatch absolute error')
    call AssertClose(reconstructed_h_at_failure_k,equilibrium_final_value, &
        100.0_KREAL*epsilon(1.0_KREAL),'h reconstructed final value retained')
    call AssertClose(stored_h_at_failure_k,equilibrium_h(1,1,1), &
        100.0_KREAL*epsilon(1.0_KREAL),'h stored perturbed value retained')
    call AssertClose(reconstructed_b_at_failure_k,equilibrium_b_before(1), &
        100.0_KREAL*epsilon(1.0_KREAL),'unperturbed b reconstructed value retained')
    call AssertClose(stored_b_at_failure_k,equilibrium_b_before(1), &
        100.0_KREAL*epsilon(1.0_KREAL),'unperturbed b stored value retained')
    call AssertClose(equilibrium_preprojection_value,equilibrium_maxwell_value+ &
        equilibrium_shakhov_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'mismatch h preprojection sum identity')
    call AssertClose(equilibrium_final_value,equilibrium_preprojection_value+ &
        equilibrium_projection_correction,100.0_KREAL*epsilon(1.0_KREAL), &
        'mismatch h final sum identity')
    call AssertClose(equilibrium_maxwell_negative_fraction,WeightedNegativeFraction(maxwell_h,wVel), &
        100.0_KREAL*epsilon(1.0_KREAL),'mismatch h Maxwellian fraction retained')
    call AssertClose(equilibrium_preprojection_negative_fraction, &
        WeightedNegativeFraction(h_preprojection,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'mismatch h preprojection fraction retained')
    call AssertClose(equilibrium_final_negative_fraction, &
        WeightedNegativeFraction(equilibrium_h_before,wVel),100.0_KREAL*epsilon(1.0_KREAL), &
        'mismatch h final fraction retained')
    equilibrium_h(:,1,1)=equilibrium_h_before
    call EvaluatePredictedEquilibriumBreakdown(1,1,1,.true., &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
        prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction, &
        reconstruction_max_abs_h,reconstruction_max_abs_b,reconstruction_max_rel_h, &
        reconstruction_max_rel_b,reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k,reconstructed_b_at_failure_k, &
        stored_b_at_failure_k,equilibrium_breakdown_stage,ok)
    call AssertTrue(ok .and. equilibrium_breakdown_stage== &
        IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS,'restored h exact reconstruction succeeds')
    call EvaluateMicroPointBreakdown(1,1,1,.true.,current_value,stored_delta,equilibrium_value, &
        stream_rhs,collision_rhs,neighbor_rhs,dt_pseudo,tau_value,diagonal_total,rhs_total, &
        post_sweep_delta,ok)
    call AssertTrue(ok,'h micro breakdown succeeds')
    call AssertClose(current_value,h_store(1),1.0e-14_KREAL,'h breakdown current value')
    call AssertClose(stored_delta,delta_h_before(1),1.0e-14_KREAL,'h breakdown stored delta')
    call AssertClose(equilibrium_value,equilibrium_h(1,1,1),1.0e-14_KREAL,'h breakdown equilibrium')
    call AssertClose(neighbor_rhs,0.0_KREAL,1.0e-14_KREAL,'zero-face h neighbor rhs')
    call AssertClose(dt_pseudo,ctr(1,1)%dt_local*implicit_ncfl_active/CFL,1.0e-14_KREAL, &
        'micro breakdown pseudo time step')
    call AssertClose(diagonal_total,1.0_KREAL/dt_pseudo+1.0_KREAL/tau_value,1.0e-14_KREAL, &
        'zero-face h diagonal')
    call AssertClose(rhs_total,stream_rhs+collision_rhs+neighbor_rhs,1.0e-14_KREAL, &
        'h breakdown total rhs')
    call AssertClose(post_sweep_delta,rhs_total/diagonal_total,1.0e-14_KREAL, &
        'h breakdown post-sweep delta')
    call EvaluateMicroPointBreakdown(1,1,2,.false.,current_value,stored_delta,equilibrium_value, &
        stream_rhs,collision_rhs,neighbor_rhs,dt_pseudo,tau_value,diagonal_total,rhs_total, &
        post_sweep_delta,ok)
    call AssertTrue(ok,'b micro breakdown succeeds')
    call AssertClose(current_value,b_store(2),1.0e-14_KREAL,'b breakdown current value')
    call AssertClose(stored_delta,delta_b_before(2),1.0e-14_KREAL,'b breakdown stored delta')
    call AssertClose(equilibrium_value,equilibrium_b(2,1,1),1.0e-14_KREAL,'b breakdown equilibrium')
    call AssertClose(neighbor_rhs,0.0_KREAL,1.0e-14_KREAL,'zero-face b neighbor rhs')
    call AssertClose(diagonal_total,1.0_KREAL/dt_pseudo+1.0_KREAL/tau_value,1.0e-14_KREAL, &
        'zero-face b diagonal')
    call AssertClose(rhs_total,stream_rhs+collision_rhs+neighbor_rhs,1.0e-14_KREAL, &
        'b breakdown total rhs')
    call AssertClose(post_sweep_delta,rhs_total/diagonal_total,1.0e-14_KREAL, &
        'b breakdown post-sweep delta')
    call AssertTrue(maxval(abs(delta_h(:,1,1)-delta_h_before))<=0.0_KREAL .and. &
        maxval(abs(delta_b(:,1,1)-delta_b_before))<=0.0_KREAL, &
        'micro breakdown does not mutate stored deltas')
    call ComputeMicroDefectNorm(norm_h_before,norm_b_before)
    call SolveMicroLUSGS(ok,freeze_wall_closure=.true.)
    call ComputeMicroDefectNorm(norm_h_after,norm_b_after)
    call AssertTrue(ok .and. norm_h_after<norm_h_before .and. norm_b_after<norm_b_before, &
                    'collision-only defects decrease')
    call AssertTrue(implicit_micro_exchange_calls==1, &
        'one micro sweep records only the forward-backward exchange')
    call AssertTrue(all([implicit_micro_forward_seconds,implicit_micro_backward_seconds, &
        implicit_micro_initial_exchange_seconds,implicit_micro_forward_backward_exchange_seconds, &
        implicit_micro_post_backward_exchange_seconds]>=0.0_KREAL), &
        'micro LU-SGS timing components are nonnegative')
    micro_timing_before=[implicit_micro_forward_seconds,implicit_micro_backward_seconds, &
        implicit_micro_initial_exchange_seconds,implicit_micro_forward_backward_exchange_seconds, &
        implicit_micro_post_backward_exchange_seconds]
    call SolveMicroLUSGS(ok,freeze_wall_closure=.true.)
    call AssertTrue(ok .and. implicit_micro_exchange_calls==2, &
        'micro LU-SGS exchange call counter accumulates')
    call AssertTrue(all([implicit_micro_forward_seconds,implicit_micro_backward_seconds, &
        implicit_micro_initial_exchange_seconds,implicit_micro_forward_backward_exchange_seconds, &
        implicit_micro_post_backward_exchange_seconds]>=micro_timing_before), &
        'micro LU-SGS timing components accumulate')
    call ResetImplicitCorrections()
    implicit_micro_exchange_calls=0
    Implicit_MicroSweeps=2
    call SolveMicroLUSGS(ok,freeze_wall_closure=.true.)
    call AssertTrue(ok .and. implicit_micro_exchange_calls==3, &
        'two micro sweeps record forward-backward and inter-sweep exchanges')
    Implicit_MicroSweeps=1
    BC_Type(3)=BCconst_wall_isothermal
    ACC=0.0_KREAL
    call GetMicroFaceValue(1,1,3,3,un,face_length,neighbor_h,neighbor_b,ok)
    call AssertTrue(ok .and. un<0.0_KREAL,'specular bottom-wall outgoing velocity accepted')
    call AssertClose(neighbor_h,delta_h(mirror_v_index(3),1,1),1.0e-14_KREAL, &
        'micro face dispatches specular h correction')
    call AssertClose(neighbor_b,delta_b(mirror_v_index(3),1,1),1.0e-14_KREAL, &
        'micro face dispatches specular b correction')

    ACC=1.0_KREAL
    call ComputeDiffuseWallFaceDelta(Prim_wall,0.0_KREAL,-1.0_KREAL,Lside, &
        delta_h(:,1,1),delta_b(:,1,1),diffuse_face_h,diffuse_face_b,rho_delta,ok)
    call GetMicroFaceValue(1,1,3,3,un,face_length,neighbor_h,neighbor_b,ok)
    call AssertClose(neighbor_h,diffuse_face_h(3),1.0e-14_KREAL, &
        'micro face dispatches diffuse h endpoint')
    call AssertClose(neighbor_b,diffuse_face_b(3),1.0e-14_KREAL, &
        'micro face dispatches diffuse b endpoint')

    accommodations=[0.0_KREAL,0.6_KREAL,1.0_KREAL]
    wall_normal_x=[-1.0_KREAL,1.0_KREAL,0.0_KREAL,0.0_KREAL]
    wall_normal_y=[0.0_KREAL,0.0_KREAL,-1.0_KREAL,1.0_KREAL]
    wall_side=[Lside,Rside,Lside,Rside]
    frozen_source_h=[0.7_KREAL,-0.3_KREAL,0.2_KREAL,0.5_KREAL]
    frozen_source_b=[0.4_KREAL,-0.1_KREAL,0.3_KREAL,0.6_KREAL]
    do accommodation_index=1,size(accommodations)
        ACC=accommodations(accommodation_index)
        do face_test=1,4
            BC_Type=BCconst_freestream
            BC_Type(face_test)=BCconst_wall_isothermal
            delta_h(:,1,1)=frozen_source_h
            delta_b(:,1,1)=frozen_source_b
            call ComputeMaxwellWallFaceDelta(Prim_wall,wall_normal_x(face_test), &
                wall_normal_y(face_test),wall_side(face_test),ACC,frozen_source_h, &
                frozen_source_b,frozen_face_h,frozen_face_b,rho_delta,ok)
            call AssertTrue(ok,'frozen Maxwell wall oracle succeeds')
            call PrepareFrozenMicroWallState(1,1,frozen_wall_state,ok)
            call AssertTrue(ok .and. frozen_wall_state%build_calls==1, &
                'one frozen Maxwell wall vector is built per cell sweep')
            delta_h(:,1,1)=frozen_source_h+[0.11_KREAL,0.12_KREAL,0.13_KREAL,0.14_KREAL]
            delta_b(:,1,1)=frozen_source_b+[0.21_KREAL,0.22_KREAL,0.23_KREAL,0.24_KREAL]
            do k_test=1,Nvel
                call GetMicroFaceValue(1,1,k_test,face_test,un,face_length,neighbor_h, &
                    neighbor_b,ok,frozen_wall_state)
                call AssertTrue(ok,'frozen Maxwell wall scalar lookup succeeds')
                call AssertClose(neighbor_h,frozen_face_h(k_test),1.0e-14_KREAL, &
                    'frozen Maxwell h lookup matches full-vector oracle')
                call AssertClose(neighbor_b,frozen_face_b(k_test),1.0e-14_KREAL, &
                    'frozen Maxwell b lookup matches full-vector oracle')
            enddo
        enddo
    enddo
    FluxWall_type=FluxWall_Kinetic
    BC_Type=BCconst_freestream
    BC_Type(3)=BCconst_wall_isothermal
    delta_h(:,1,1)=frozen_source_h
    delta_b(:,1,1)=frozen_source_b
    call ComputeDiffuseWallFaceDelta(Prim_wall,0.0_KREAL,-1.0_KREAL,Lside, &
        frozen_source_h,frozen_source_b,frozen_face_h,frozen_face_b,rho_delta,ok)
    call AssertTrue(ok,'frozen diffuse wall oracle succeeds')
    call PrepareFrozenMicroWallState(1,1,frozen_wall_state,ok)
    call AssertTrue(ok .and. frozen_wall_state%build_calls==1, &
        'one frozen diffuse wall vector is built per cell sweep')
    delta_h(:,1,1)=frozen_source_h+0.5_KREAL
    delta_b(:,1,1)=frozen_source_b+0.5_KREAL
    do k_test=1,Nvel
        call GetMicroFaceValue(1,1,k_test,3,un,face_length,neighbor_h,neighbor_b,ok, &
            frozen_wall_state)
        call AssertClose(neighbor_h,frozen_face_h(k_test),1.0e-14_KREAL, &
            'frozen diffuse h lookup matches full-vector oracle')
        call AssertClose(neighbor_b,frozen_face_b(k_test),1.0e-14_KREAL, &
            'frozen diffuse b lookup matches full-vector oracle')
    enddo
    FluxWall_type=FluxWall_Maxwell
    ACC=1.0_KREAL
    hface(1,1)%cosx=0.6_KREAL
    hface(1,1)%cosy=0.8_KREAL
    delta_h(:,1,1)=frozen_source_h
    delta_b(:,1,1)=frozen_source_b
    call PrepareFrozenMicroWallState(1,1,frozen_wall_state,ok)
    call AssertTrue(ok .and. frozen_wall_state%build_calls==0 .and. &
        frozen_wall_state%fallback(3),'slanted wall retains full-vector fallback')
    call ComputeMaxwellWallFaceDelta(Prim_wall,-0.6_KREAL,-0.8_KREAL,Lside,ACC, &
        delta_h(:,1,1),delta_b(:,1,1),frozen_face_h,frozen_face_b,rho_delta,ok)
    call GetMicroFaceValue(1,1,3,3,un,face_length,neighbor_h,neighbor_b,ok,frozen_wall_state)
    call AssertClose(neighbor_h,frozen_face_h(3),1.0e-14_KREAL, &
        'slanted wall fallback preserves h closure')
    call AssertClose(neighbor_b,frozen_face_b(3),1.0e-14_KREAL, &
        'slanted wall fallback preserves b closure')
    call FinishTests('implicit_micro')
    call MPI_Finalize(ierr_local)
end program test_implicit_micro

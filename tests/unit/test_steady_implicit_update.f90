program test_steady_implicit_update
    use mpi
    use flow_data
    use Tools, only: GetConserved, DiscreteMaxwell_Vec, RelativeOrAbsoluteResidual, &
        IsNearZeroResidualScale
    use mod_ImplicitState
    use mod_SteadyImplicit
    use mod_ImplicitDiagnostics, only: InitializeImplicitDiagnostics, UpdateImplicitResidualDiagnostics
    use TestSupport
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer :: ierr_local
    integer(KINT) :: m
    logical :: ok,accepted
    type(ImplicitFailureRecord) :: diagnostic
    real(KREAL), target :: h_store(4),b_store(4)
    real(KREAL) :: prim_test(4),alpha,kinetic_perturbation(4)
    real(KREAL) :: h_accepted(4),b_accepted(4),w_accepted(4),expected_macro(4)
    real(KREAL) :: candidate_h(4),candidate_b(4)
    real(KREAL) :: u_original(4),v_original(4),u_huge(4),v_huge(4)
    real(KREAL) :: min_alpha_original

    call MPI_Init(ierr_local)
    call MPI_Comm_rank(MPI_COMM_WORLD,myID,ierr_local)
    nx=1; ny=1; Nx_all=1; Ny_all=1; Nu=2; Nv=2; Nvel=4; ck=1; gamma=5.0_KREAL/3.0_KREAL
    CFL=0.5_KREAL; Implicit_nCFL=10.0_KREAL; Implicit_MacroSweeps=2; Implicit_MicroSweeps=1
    Implicit_MinAlpha=1.0e-6_KREAL; Implicit_MaxRetries=2; Implicit_MomentTol=1.0e-12_KREAL
    Implicit_ResidualTol=1.0e-8_KREAL
    flag_vislaw=PowerLaw; mu_ref=0.05_KREAL; omega=0.5_KREAL; Pr=1.0_KREAL
    allocate(uVel(4),vVel(4),wVel(4),ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    allocate(res_macro_flux_store(4,nx,ny),res_h_flux_store(Nvel,nx,ny),res_b_flux_store(Nvel,nx,ny))
    uVel=[-1.0_KREAL,1.0_KREAL,-1.0_KREAL,1.0_KREAL]
    vVel=[-1.0_KREAL,-1.0_KREAL,1.0_KREAL,1.0_KREAL]; wVel=0.25_KREAL

    candidate_h=1.0_KREAL; candidate_b=1.0_KREAL
    candidate_h(2)=ieee_value(0.0_KREAL,ieee_quiet_nan)
    call ClassifyMicroCandidate([1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL], &
        [1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL],candidate_h,candidate_b, &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(diagnostic%reason==IMPLICIT_FAILURE_H_NONFINITE .and. &
        diagnostic%velocity_k==2,'h nonfinite classified')

    candidate_h=1.0_KREAL; candidate_b=1.0_KREAL
    candidate_b(3)=ieee_value(0.0_KREAL,ieee_quiet_nan)
    call ClassifyMicroCandidate([1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL], &
        [1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL],candidate_h,candidate_b, &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(diagnostic%reason==IMPLICIT_FAILURE_B_NONFINITE .and. &
        diagnostic%velocity_k==3,'b nonfinite classified')

    candidate_h=[-0.01_KREAL,-0.02_KREAL,1.0_KREAL,1.0_KREAL]
    candidate_b=1.0_KREAL
    call ClassifyMicroCandidate([1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL], &
        [1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL],candidate_h,candidate_b, &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(.not.diagnostic%valid,'finite physical h negativity is diagnostic only')
    call AssertTrue(WeightedNegativeFraction(candidate_h,wVel)>1.0e-4_KREAL, &
        'h diagnostic exceeds former threshold')

    candidate_h=1.0_KREAL
    candidate_b=[-0.01_KREAL,-0.03_KREAL,1.0_KREAL,1.0_KREAL]
    call ClassifyMicroCandidate([1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL], &
        [1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL],candidate_h,candidate_b, &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(.not.diagnostic%valid,'finite physical b negativity is diagnostic only')
    call AssertTrue(WeightedNegativeFraction(candidate_b,wVel)>1.0e-4_KREAL, &
        'b diagnostic exceeds former threshold')

    candidate_h=0.0_KREAL; candidate_b=0.0_KREAL
    call ClassifyMicroCandidate([1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL], &
        [1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL],candidate_h,candidate_b, &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(diagnostic%reason==IMPLICIT_FAILURE_MOMENT_NONPHYSICAL .and. &
        diagnostic%velocity_k==0,'nonphysical moment classified')

    candidate_h=1.0_KREAL; candidate_b=-3.0_KREAL
    call ClassifyMicroCandidate(candidate_h,candidate_h,candidate_h,candidate_b, &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(diagnostic%reason==IMPLICIT_FAILURE_MOMENT_NONPHYSICAL, &
        'nonpositive internal energy classified')

    candidate_h=1.0_KREAL; candidate_b=1.0_KREAL
    u_huge=huge(1.0_KREAL); v_huge=0.0_KREAL
    call ClassifyMicroCandidate([1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL], &
        [1.0_KREAL,1.0_KREAL,1.0_KREAL,1.0_KREAL],candidate_h,candidate_b, &
        wVel,u_huge,v_huge,diagnostic)
    call AssertTrue(all(ieee_is_finite(candidate_h)) .and. all(ieee_is_finite(candidate_b)) .and. &
        .not.all(ieee_is_finite(diagnostic%candidate_moment)) .and. &
        diagnostic%reason==IMPLICIT_FAILURE_MOMENT_NONPHYSICAL, &
        'finite candidate with nonfinite recovered moment is rejected')

    prim_test=[1.0_KREAL,0.0_KREAL,0.0_KREAL,0.5_KREAL]
    call DiscreteMaxwell_Vec(h_store,b_store,prim_test)
    ctr(1,1)%h_vec=>h_store; ctr(1,1)%b_vec=>b_store
    ctr(1,1)%conVars=GetConserved(prim_test)
    ctr(1,1)%area=1.0_KREAL; ctr(1,1)%dt_local=0.1_KREAL
    ctr(1,1)%length=[1.0_KREAL,1.0_KREAL]; ctr(1,1)%x=0.5_KREAL; ctr(1,1)%y=0.5_KREAL
    vface%length=0.0_KREAL; vface%cosx=1.0_KREAL; vface%cosy=0.0_KREAL
    hface%length=0.0_KREAL; hface%cosx=0.0_KREAL; hface%cosy=1.0_KREAL
    ID_X_minus=MPI_PROC_NULL; ID_X_plus=MPI_PROC_NULL
    ID_Y_minus=MPI_PROC_NULL; ID_Y_plus=MPI_PROC_NULL; BC_Type=BCconst_freestream
    call InitializeSteadyImplicit(ok)
    call AssertTrue(ok,'implicit initialization')
    call AssertTrue(all(h_store>0.0_KREAL) .and. all(b_store>0.0_KREAL), &
        'legal Maxwellian remains the accepted state')

    call AssertClose(RelativeOrAbsoluteResidual(4.0e-13_KREAL,1.28e-8_KREAL,1.0e-6_KREAL), &
        4.0e-13_KREAL,epsilon(1.0_KREAL),'absolute floor wins after near-equilibrium convergence')
    call AssertClose(RelativeOrAbsoluteResidual(5.0e-9_KREAL,1.0e-2_KREAL,1.0e-6_KREAL), &
        5.0e-7_KREAL,epsilon(1.0_KREAL),'ordinary nonzero residual retains relative criterion')
    call AssertTrue(IsNearZeroResidualScale(1.0e-10_KREAL,1.0_KREAL), &
        'vanishing conserved component selects absolute residual')
    call AssertTrue(.not.IsNearZeroResidualScale(1.0e-4_KREAL,1.0_KREAL), &
        'ordinary conserved component retains relative residual')

    rate_w(:,1,1)=[2.0_KREAL,3.0_KREAL,4.0_KREAL,5.0_KREAL]
    stream_h(:,1,1)=0.0_KREAL
    stream_b(:,1,1)=0.0_KREAL
    call UpdateImplicitResidualDiagnostics()
    expected_macro=sqrt(real(Nx_all*Ny_all,KREAL)*(ctr(1,1)%area*rate_w(:,1,1))**2)/ &
        (ctr(1,1)%area*abs(ctr(1,1)%conVars+ctr(1,1)%dt_local*rate_w(:,1,1))+SMV)
    do m=1,4
        call AssertClose(implicit_macro_residual(m),expected_macro(m), &
            100.0_KREAL*epsilon(1.0_KREAL)*max(1.0_KREAL,expected_macro(m)), &
            'component implicit macro residual')
    enddo

    min_alpha_original=Implicit_MinAlpha
    Implicit_MinAlpha=0.75_KREAL
    h_accepted=h_store; b_accepted=b_store; w_accepted=ctr(1,1)%conVars
    delta_h=0.0_KREAL; delta_b=0.0_KREAL
    delta_b(1,1,1)=-1.1_KREAL*b_store(1)
    call AssertTrue(WeightedNegativeFraction(b_store+delta_b(:,1,1),wVel)>1.0e-4_KREAL, &
        'line-search candidate exceeds former b threshold')
    call ClassifyMicroCandidate(h_store,b_store,h_store,b_store+delta_b(:,1,1), &
        wVel,uVel,vVel,diagnostic)
    call AssertTrue(.not.diagnostic%valid,'line-search candidate has physical moments')
    call FindAdmissibleMicroAlpha(alpha,ok,.false.)
    call AssertTrue(ok .and. alpha==1.0_KREAL,'negative Shakhov-compatible candidate accepts alpha one')
    call AssertTrue(.not.implicit_failure%valid,'accepted trial leaves no stale failure record')
    call AssertTrue(maxval(abs(h_store-h_accepted))==0.0_KREAL .and. &
        maxval(abs(b_store-b_accepted))==0.0_KREAL .and. &
        maxval(abs(ctr(1,1)%conVars-w_accepted))==0.0_KREAL, &
        'alpha search does not mutate accepted state')

    call FindAdmissibleMicroAlpha(alpha,ok)
    call AssertTrue(.not.ok .and. alpha==0.5_KREAL .and. implicit_failure%valid .and. &
        implicit_failure%reason==IMPLICIT_FAILURE_B_NEGATIVE_FRACTION, &
        'default alpha search preserves the negative fraction gate')
    call AssertTrue(maxval(abs(h_store-h_accepted))==0.0_KREAL .and. &
        maxval(abs(b_store-b_accepted))==0.0_KREAL .and. &
        maxval(abs(ctr(1,1)%conVars-w_accepted))==0.0_KREAL, &
        'rejected negative-fraction trial does not mutate accepted state')

    u_original=uVel; v_original=vVel
    uVel=u_huge; vVel=v_huge
    h_accepted=h_store; b_accepted=b_store; w_accepted=ctr(1,1)%conVars
    delta_h=0.0_KREAL; delta_b=0.0_KREAL
    call FindAdmissibleMicroAlpha(alpha,ok)
    call AssertTrue(.not.ok .and. alpha==0.5_KREAL .and. implicit_failure%valid .and. &
        implicit_failure%reason==IMPLICIT_FAILURE_MOMENT_NONPHYSICAL, &
        'nonphysical moment causes terminal alpha-search failure')
    call AssertTrue(maxval(abs(h_store-h_accepted))==0.0_KREAL .and. &
        maxval(abs(b_store-b_accepted))==0.0_KREAL .and. &
        maxval(abs(ctr(1,1)%conVars-w_accepted))==0.0_KREAL, &
        'nonphysical moment trial does not mutate accepted state')
    uVel=u_original; vVel=v_original
    Implicit_MinAlpha=min_alpha_original

    res_macro_flux_store=0.0_KREAL; res_h_flux_store=0.0_KREAL; res_b_flux_store=0.0_KREAL
    call ApplySteadyImplicitUpdate(accepted)
    call AssertTrue(accepted,'zero-residual update accepted')
    call AssertTrue(maxval(abs(delta_w))==0.0_KREAL,'zero macro correction')
    call AssertTrue(maxval(abs(delta_h))<1.0e-13_KREAL,'equilibrium zero h correction')
    call AssertTrue(maxval(abs(delta_b))<1.0e-13_KREAL,'equilibrium zero b correction')
    call AssertTrue(implicit_moment_error<=Implicit_MomentTol,'accepted moment consistency')
    call AssertTrue(ImplicitConverged(),'zero residual satisfies implicit stopping criterion')

    kinetic_perturbation=1.0e-3_KREAL*[1.0_KREAL,-1.0_KREAL,-1.0_KREAL,1.0_KREAL]
    h_store=h_store+kinetic_perturbation
    call InitializeImplicitDiagnostics()
    call AssertTrue(implicit_retry_total==0,'diagnostic initialization resets cumulative retries')
    call ApplySteadyImplicitUpdate(accepted)
    call AssertTrue(accepted,'non-Maxwellian collision update accepted')
    call AssertTrue(.not.ImplicitConverged(), &
        'finite collision update does not satisfy physical steady residual')
    call FinishTests('steady_implicit_update')
    call MPI_Finalize(ierr_local)
end program test_steady_implicit_update

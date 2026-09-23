program test_unsteady_micro
    use mpi
    use Precision
    use flow_data
    use const_para, only: PowerLaw
    use InitializeFlow, only: InitVelocityGaussHermite, FlattenVelocitySpace
    use control_para, only: TimeIntegration_method, Unsteady_Theta, Unsteady_MaxInner, &
        FluxWall_type
    use flow_para, only: Ma, Pr, Re, flag_vislaw, omega, mu_ref
    use mpi_para, only: ID_X_minus, ID_X_plus, ID_Y_minus, ID_Y_plus
    use Tools, only: DiscreteMaxwell_Vec
    use mod_ImplicitCommon, only: MomentConserved, ProjectToConservedMoments
    use mod_UnsteadyImplicit, only: MicroPhysicalResidual, MicroUnsteadyDiagonal, CollisionForState
    use mod_ImplicitState
    use mod_ImplicitBoundary, only: BuildStructuredMirrorMaps, ComputeMaxwellWallFaceDelta
    use mod_ImplicitMicro, only: SolveMicroLUSGS, PrepareFrozenMicroWallState, GetMicroFaceValue, &
        FrozenMicroWallState
    use mod_ImplicitDiagnostics, only: InitializeImplicitDiagnostics
    use TestSupport
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer :: ierr_local,k
    logical :: ok
    real(KREAL) :: residual
    real(KREAL), target :: h_store(4),b_store(4)
    real(KREAL) :: physical_dt,flux_weight(1,1),un,face_length,neighbor_h,neighbor_b
    real(KREAL) :: rho_delta,reference_h(4),reference_b(4)
    real(KREAL) :: seed_h(4),seed_b(4),changed_h(4),changed_b(4)
    real(KREAL) :: first_face_h(4),first_face_b(4),second_face_h(4),second_face_b(4)
    real(KREAL) :: collision_target(4),collision_moment(4),collision_prim(4),collision_scale
    real(KREAL), allocatable :: collision_reference_h(:),collision_reference_b(:)
    real(KREAL), allocatable :: collision_h(:),collision_b(:),collision_input_h(:),collision_input_b(:)
    type(FrozenMicroWallState) :: first_wall_state,second_wall_state

    call MPI_Init(ierr_local)
    residual=MicroPhysicalResidual(1.0_KREAL,1.1_KREAL,0.2_KREAL,0.4_KREAL, &
        0.3_KREAL,0.5_KREAL,0.1_KREAL,0.2_KREAL,0.4_KREAL,0.5_KREAL,0.25_KREAL,0.5_KREAL)
    call AssertClose(residual,0.65_KREAL,1.0e-14_KREAL,'Zhu micro residual algebra')
    call AssertClose(MicroUnsteadyDiagonal(0.5_KREAL,0.5_KREAL,0.25_KREAL,0.2_KREAL,4.0_KREAL), &
        5.5_KREAL,1.0e-14_KREAL,'micro physical collision and upwind diagonal')

    !GH16 collision-moment contract: a projected old state must have a
    !discrete Maxwell/Shakhov collision term with all four conserved moments zero.
    Nu=16; Nv=16; Nvel=Nu*Nv; CK=1; Ma=0.1_KREAL; Re=1.0_KREAL; Pr=2.0_KREAL/3.0_KREAL
    gamma=5.0_KREAL/3.0_KREAL
    flag_vislaw=PowerLaw; omega=0.5_KREAL; mu_ref=sqrt(gamma/2.0_KREAL)*Ma/Re
    allocate(uSpace(Nvel),vSpace(Nvel),weight(Nvel))
    allocate(uVel(Nvel),vVel(Nvel),wVel(Nvel))
    call InitVelocityGaussHermite(Nu,Nv)
    call FlattenVelocitySpace()
    allocate(collision_reference_h(Nvel),collision_reference_b(Nvel),collision_input_h(Nvel), &
        collision_input_b(Nvel),collision_h(Nvel),collision_b(Nvel))
    collision_prim=[1.0_KREAL,Ma*sqrt(gamma/2.0_KREAL),0.0_KREAL,1.0_KREAL]
    call DiscreteMaxwell_Vec(collision_reference_h,collision_reference_b,collision_prim)
    call MomentConserved(collision_reference_h,collision_reference_b,uVel,vVel,wVel,collision_target)
    collision_input_h=collision_reference_h*(1.0_KREAL+2.0e-3_KREAL*(1.0_KREAL+0.1_KREAL*uVel+0.05_KREAL*vVel))
    collision_input_b=collision_reference_b*(1.0_KREAL+3.0e-3_KREAL*(1.0_KREAL-0.03_KREAL*uVel+0.02_KREAL*vVel))
    call ProjectToConservedMoments(collision_input_h,collision_input_b,collision_target, &
        uVel,vVel,wVel,collision_reference_h,collision_reference_b,ok)
    call AssertTrue(ok,'GH16 collision input moment projection succeeds')
    call MomentConserved(collision_input_h,collision_input_b,uVel,vVel,wVel,collision_moment)
    collision_scale=max(1.0_KREAL,maxval(abs(collision_target)))
    call AssertTrue(maxval(abs(collision_moment-collision_target))<=1.0e-12_KREAL*collision_scale, &
        'GH16 collision input is discretely compatible')
    call CollisionForState(collision_target,collision_input_h,collision_input_b,collision_h,collision_b,ok)
    call AssertTrue(ok,'GH16 captured collision state is finite')
    call MomentConserved(collision_h,collision_b,uVel,vVel,wVel,collision_moment)
    write(*,'(A,1P,4E16.8)') 'captured collision moments=',collision_moment
    call AssertTrue(maxval(abs(collision_moment))<=1.0e-12_KREAL*collision_scale, &
        'GH16 captured collision moments conserve rho rhoU rhoV energy')
    deallocate(collision_reference_h,collision_reference_b,collision_input_h,collision_input_b, &
        collision_h,collision_b)
    deallocate(uSpace,vSpace,weight,uVel,vVel,wVel)

    nx=1; ny=1; Nx_all=1; Ny_all=1; Nu=2; Nv=2; Nvel=4; ck=1
    CFL=0.5_KREAL; Implicit_nCFL=10.0_KREAL; Implicit_MicroSweeps=1
    Unsteady_Theta=0.5_KREAL; Unsteady_MaxInner=80
    TimeIntegration_method='ImplicitUnsteady'
    allocate(uVel(Nvel),vVel(Nvel),wVel(Nvel))
    allocate(ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    uVel=[-1.0_KREAL,1.0_KREAL,-1.0_KREAL,1.0_KREAL]
    vVel=[-1.0_KREAL,-1.0_KREAL,1.0_KREAL,1.0_KREAL]
    wVel=0.25_KREAL
    h_store=0.2_KREAL; b_store=0.1_KREAL
    ctr(1,1)%h_vec=>h_store; ctr(1,1)%b_vec=>b_store
    ctr(1,1)%area=1.0_KREAL; ctr(1,1)%dt_local=0.1_KREAL
    vface%length=0.0_KREAL; vface%cosx=1.0_KREAL; vface%cosy=0.0_KREAL
    hface%length=0.0_KREAL; hface%cosx=0.0_KREAL; hface%cosy=1.0_KREAL
    hface(1,1)%length=1.0_KREAL
    ID_X_minus=MPI_PROC_NULL; ID_X_plus=MPI_PROC_NULL
    ID_Y_minus=MPI_PROC_NULL; ID_Y_plus=MPI_PROC_NULL
    BC_Type=BCconst_freestream; BC_Type(3)=BCconst_wall_isothermal
    FluxWall_type=FluxWall_Maxwell; ACC=0.6_KREAL
    Prim_wall=[1.0_KREAL,0.0_KREAL,0.0_KREAL,1.0_KREAL]
    physical_dt=0.25_KREAL; flux_weight=0.5_KREAL
    seed_h=[0.7_KREAL,-0.3_KREAL,0.2_KREAL,0.5_KREAL]
    seed_b=[0.4_KREAL,-0.1_KREAL,0.3_KREAL,0.6_KREAL]
    changed_h=seed_h+[0.31_KREAL,0.27_KREAL,0.23_KREAL,0.19_KREAL]
    changed_b=seed_b+[0.17_KREAL,0.13_KREAL,0.29_KREAL,0.11_KREAL]
    call AllocateImplicitState()
    call BuildStructuredMirrorMaps()
    tau_tilde(1,1)=1.0_KREAL
    stream_h(:,1,1)=[0.11_KREAL,-0.07_KREAL,0.05_KREAL,-0.03_KREAL]
    stream_b(:,1,1)=[0.09_KREAL,-0.06_KREAL,0.04_KREAL,-0.02_KREAL]
    equilibrium_h(:,1,1)=0.8_KREAL
    equilibrium_b(:,1,1)=0.4_KREAL
    call InitializeImplicitDiagnostics()

    delta_h(:,1,1)=seed_h; delta_b(:,1,1)=seed_b
    call PrepareFrozenMicroWallState(1,1,first_wall_state,ok)
    call AssertTrue(ok .and. first_wall_state%build_calls==1, &
        'unsteady fixture builds one frozen ACC=0.6 wall vector')
    call ComputeMaxwellWallFaceDelta(Prim_wall,0.0_KREAL,-1.0_KREAL,Lside,ACC,seed_h,seed_b, &
        reference_h,reference_b,rho_delta,ok)
    call AssertTrue(ok,'unsteady frozen wall reference succeeds')
    do k=1,Nvel
        call GetMicroFaceValue(1,1,k,3,un,face_length,neighbor_h,neighbor_b,ok,first_wall_state)
        call AssertTrue(ok .and. ieee_is_finite(neighbor_h) .and. ieee_is_finite(neighbor_b), &
            'unsteady frozen wall lookup is finite')
        call AssertClose(neighbor_h,reference_h(k),1.0e-14_KREAL, &
            'unsteady lookup uses one frozen h face vector')
        call AssertClose(neighbor_b,reference_b(k),1.0e-14_KREAL, &
            'unsteady lookup uses one frozen b face vector')
    enddo

    implicit_micro_exchange_calls=0
    call SolveMicroLUSGS(ok,physical_dt,flux_weight,Unsteady_Theta, &
        freeze_wall_closure=.true.)
    call AssertTrue(ok .and. implicit_micro_exchange_calls==1, &
        'unsteady physical-time micro sweep succeeds with one halo exchange')
    call AssertTrue(all(ieee_is_finite(delta_h(:,1,1))) .and. all(ieee_is_finite(delta_b(:,1,1))), &
        'unsteady physical-time corrections are finite')
    first_face_h=first_wall_state%delta_h_face(:,3)
    first_face_b=first_wall_state%delta_b_face(:,3)

    delta_h(:,1,1)=changed_h; delta_b(:,1,1)=changed_b
    call PrepareFrozenMicroWallState(1,1,second_wall_state,ok)
    call AssertTrue(ok .and. second_wall_state%build_calls==1, &
        'second unsteady iteration rebuilds the frozen wall vector')
    second_face_h=second_wall_state%delta_h_face(:,3)
    second_face_b=second_wall_state%delta_b_face(:,3)
    call AssertTrue(maxval(abs(second_face_h-first_face_h))>1.0e-8_KREAL .and. &
        maxval(abs(second_face_b-first_face_b))>1.0e-8_KREAL, &
        'changed corrections do not reuse the preceding frozen wall vector')
    implicit_micro_exchange_calls=0
    call SolveMicroLUSGS(ok,physical_dt,flux_weight,Unsteady_Theta, &
        freeze_wall_closure=.true.)
    call AssertTrue(ok .and. implicit_micro_exchange_calls==1, &
        'second unsteady physical-time micro call rebuilds and succeeds')

    call FinishTests('unsteady_micro')
    call MPI_Finalize(ierr_local)
end program test_unsteady_micro

module mod_S2NativeAdapter
    use mod_S2CartesianMPI
    use iso_fortran_env, only: int64
    use Tools, only: GetPrimary
    use mod_S2NativeRestart, only: ReadNativeMetadata,CommitNativeRestart
    use mod_CorrectedWF2DEquilibrium, only: analytic_maxwellian_2d
    use mod_S2MPI2IO, only: load_input,legacy_restart_path=>restart_path
    use mod_S2MPI1Layout, only: slab_lo,slab_hi,communication_clock,communication_seconds
    use mod_S2PhysicalGas, only: gas_tau_min,gas_tau_max,gas_tau_calls,s2_tau
    use mod_VelocityTables, only: VelocityMode,VelocityCount
    use mod_TimeIntegrator, only: RecordAcceptedMacroUpdateRate
    use Precision, only: KREAL,KINT
    use flow_data
    use flow_para, only: ck,Pr,S2PhysicalGas,BC_Type,S2W1_wall_mode,S2W1_time_centered,S2W1D_face_time_consistent,S2W2R_shared_predictor, &
        S2W1_wall_ux,S2W1_wall_theta,S2CartesianMesh,S2Diagnostics,S2CheckpointImport
    use control_para, only: dt_global,TimeStep_method,TimeIntegration_method,PeriodicComparisonProfile, &
        S2D1_initial_mode,MAX_ITER,iter,ComparisonFaceTauMin,ComparisonFaceTauMax, &
        ComparisonCellTauMin,ComparisonCellTauMax,ComparisonFaceTauCalls,ComparisonCellTauCalls
    use Flow_DIS, only: time_dis,Iflag_DIS,Iflag_Fringe
    use mpi_para, only: myID,NP_X,NP_Y,NP_XY,Block_XY,ID_Y_minus,ID_Y_plus
    use const_para, only: BCconst_periodic,BCconst_wall_isothermal,Time_global_fixed
    use mod_S2D1, only: IsS2CN2,IsS2Steady,compute_flux,build_targets,advance_step,init_velocity, &
        write_state,heat,make_mke,s2w2r_correct,s2w1_compute_wall_flux,s2w1d_instantaneous_flux, &
        w2_reconstruction_s,w2_wall_s,w2_wall_calls,w2_closure_s
    use mod_CorrectedWF2DClosure, only: moment_2d
    use mod_S2SteadyTermination, only: SteadyGateObserve,SteadyNeedsEvaluation
    use mod_ImplicitState, only: rate_w,stream_h,stream_b,s2_steady_rhs
    use mpi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: S2NativeInitializeState,S2NativeEvolution,S2NativeApplyUpdate,S2NativeMarkResidual
    public :: S2NativeWriteOutput,S2NativeIsReady,S2NativeRecordLowOrderStep
    public :: S2NativeAssembleSteadyResidual,S2NativeInvalidateTargets,S2NativeSaveCurrent
    public :: S2NativeSteadyConverged,S2NativeSteadyNeedsEvaluation,S2NativeDiagnosticEnabled
    ! Representation-only aliases: physical cells 1:nx,1:ny; no halo/cache rebuild.
    public :: S2NativePackState,S2NativePublishState,S2NativeReadMetadata,S2NativeCommitRestart,S2NativeOutputHeat,S2NativeRunSummary
    interface S2NativePackState
        module procedure pack_state
    end interface
    interface S2NativePublishState
        module procedure publish_state
    end interface
    public :: S2NativePhysicalCollision
    interface S2NativePublishWallDiagnostics
        module procedure PublishWallDiagnostics
    end interface
    public :: S2NativePublishWallDiagnostics
    public :: S2NativeSteadyMaxNorm

    real(KREAL),parameter::pi_d=acos(-1.0_KREAL)
    real(KREAL),parameter::rho0=1.0_KREAL,u0=0.6_KREAL,v0=0.3_KREAL,t0=1.0_KREAL
    real(KREAL),parameter::amp=0.01_KREAL,gate=2.0e-12_KREAL,tau0=0.1_KREAL
    real(KREAL),allocatable,save::cw(:,:,:),ch(:,:,:),cb(:,:,:)
    real(KREAL),allocatable,save::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:)
    real(KREAL),allocatable,save::eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
    real(KREAL),allocatable,save::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:)
    real(KREAL),allocatable,save::lower_w(:,:),lower_h(:,:),lower_b(:,:)
    real(KREAL),allocatable,save::leading_xh(:,:,:),leading_xb(:,:,:),leading_yh(:,:,:),leading_yb(:,:,:)
    real(KREAL),allocatable,save::time_xh(:,:,:),time_xb(:,:,:),time_yh(:,:,:),time_yb(:,:,:)
    real(KREAL),allocatable,save::w2r_l0h(:,:),w2r_l0b(:,:),w2r_lkh(:,:),w2r_lkb(:,:)
    real(KREAL),save::wall_balance(4)=0.0_KREAL,cumulative_bottom(4)=0.0_KREAL,cumulative_top(4)=0.0_KREAL
    integer,save::wall_unit=-1,wall_segment_unit=-1
    integer,save::w1t_cost_unit=-1
    real(KREAL),save::w1t_predictor_step_s=0.0_KREAL,w1t_predictor_wall_s=0.0_KREAL
    real(KREAL),save::w1t_official_face_s=0.0_KREAL
    real(KREAL),save::face_error=0.0_KREAL,total0(4)=0.0_KREAL
    real(KREAL),save::layout_seconds=0.0_KREAL,face_seconds=0.0_KREAL,map_seconds=0.0_KREAL
    real(KREAL),save::collision_seconds=0.0_KREAL,diagnostic_seconds=0.0_KREAL,io_seconds=0.0_KREAL
    integer,save::limited_count=0,signed_count=0,init_calls=0,evolution_calls=0,residual_calls=0
    integer,save::update_calls=0,publish_calls=0,output_calls=0,steps_unit=-1
    logical,save::initialized=.false.,evolution_ready=.false.,atomic_test_passed=.false.
    logical,save::targets_ready=.false.,steps_open=.false.,force_global_reject=.false.
    logical,save::o1_bench=.false.
    logical,save::w2_cost_enabled=.false.
    real(KREAL),save::w2_off_reconstruction_s=0.0_KREAL,w2_off_wall_s=0.0_KREAL
    integer,save::w2_off_wall_calls=0,w2_cost_unit=-1
    real(KREAL),save::o1_wall_start=0.0_KREAL
    real(KREAL),save::i1_scale(6)=1.0_KREAL,i1_initial_norm(6)=0.0_KREAL
    real(KREAL),save::i1_last_norm(6)=0.0_KREAL,i1_start_wall=0.0_KREAL
    integer,save::i1_evaluations=0,i1_consecutive=0,i1_unit=-1
    integer,save::i1_last_state_iter=-1
    logical,save::i1_diag_enabled=.false.,i1_below=.false.
    logical,save::w1d_diagnostic=.false.
    integer,save::w1d_audit_unit=-1,w1d_mid_limited=0,w1d_mid_signed=0
    real(KREAL),save::w1d_prediction_diag(12)=0.0_KREAL
    real(KREAL),save::w1d_instant_s=0.0_KREAL
    integer,save::w1d_cost_unit=-1,w1d_instant_limited=0,w1d_instant_signed=0

    integer,save::gas_unit=-1
    logical,save::w2r_metrics=.false.
    integer,save::w2r_unit=-1
    real(KREAL),save::w2r_predict_s=0,w2r_instant_s=0,w2r_apply_s=0,w2r_off_closure_s=0,w2r_immutable_error=0
    real(KREAL),save::w2r_diag(6)=0
    integer,save::w2r_active=0
contains
    subroutine native_build_targets(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,error,ok)
        real(KREAL),intent(in)::w(:,:,:),u(:),v(:),weight(:)
        real(KREAL),intent(out)::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:),eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
        real(KREAL),intent(out)::error
        logical,intent(out)::ok
        if(cart_native)then
            call CartesianTargets(w,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,error,ok)
        else
            call build_targets(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,error,ok)
        endif
    end subroutine
    function w2r_snapshot() result(data)
        real(KREAL),allocatable::data(:),pw(:,:,:),ph(:,:,:),pb(:,:,:)
        allocate(pw(4,nx,ny),ph(nx,ny,Nvel),pb(nx,ny,Nvel));call pack_state(pw,ph,pb)
        data=[reshape(pw,[size(pw)]),reshape(ph,[size(ph)]),reshape(pb,[size(pb)]), &
            reshape(cw,[size(cw)]),reshape(ch,[size(ch)]),reshape(cb,[size(cb)]), &
            reshape(mh,[size(mh)]),reshape(mb,[size(mb)]),reshape(kxh,[size(kxh)]),reshape(kxb,[size(kxb)]), &
            reshape(kyh,[size(kyh)]),reshape(kyb,[size(kyb)]),reshape(eah,[size(eah)]),reshape(eab,[size(eab)]), &
            reshape(kappa,[size(kappa)]),wall_balance,total0,real([update_calls,publish_calls],KREAL), &
            merge(1.0_KREAL,0.0_KREAL,targets_ready),time_dis]
    end function w2r_snapshot

    logical function package_mode()
        package_mode=IsS2CN2().or.IsS2Steady().or. &
            trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau'
    end function package_mode

    logical function S2NativeIsReady()
        S2NativeIsReady=initialized
    end function S2NativeIsReady

    logical function S2NativeDiagnosticEnabled()
        S2NativeDiagnosticEnabled=i1_diag_enabled
    end function S2NativeDiagnosticEnabled

    logical function S2NativeSteadyConverged()
        S2NativeSteadyConverged=IsS2Steady().and.i1_consecutive>=5
    end function S2NativeSteadyConverged

    logical function S2NativeSteadyNeedsEvaluation()
        S2NativeSteadyNeedsEvaluation=IsS2Steady().and.SteadyNeedsEvaluation(iter,i1_last_state_iter)
    end function S2NativeSteadyNeedsEvaluation

    real(KREAL) function S2NativeSteadyMaxNorm()
        S2NativeSteadyMaxNorm=maxval(i1_last_norm)
    end function S2NativeSteadyMaxNorm

    subroutine abort_native(message,code)
        character(len=*),intent(in)::message
        integer,intent(in)::code
        integer::ierr
        write(*,'(A,I0,A)')'S2CN2 FAILURE rank=',myID,': '//trim(message)
        if(initialized.and.IsS2CN2())call SaveRejectedAcceptedState(message)
        call MPI_Abort(MPI_COMM_WORLD,code,ierr)
    end subroutine abort_native

    subroutine pack_state(w,h,b)
        real(KREAL),intent(out)::w(:,:,:),h(:,:,:),b(:,:,:)
        integer::i,j
        do j=1,ny;do i=1,nx
            w(:,i,j)=ctr(i,j)%conVars
            h(i,j,:)=ctr(i,j)%h_vec
            b(i,j,:)=ctr(i,j)%b_vec
        enddo;enddo
    end subroutine pack_state

    subroutine publish_state(w,h,b,accept)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:)
        logical,intent(in)::accept
        integer::i,j
        if(.not.accept)then;evolution_ready=.false.;return;endif
        do j=1,ny;do i=1,nx
            ctr(i,j)%conVars=w(:,i,j)
            ctr(i,j)%h_vec=h(i,j,:)
            ctr(i,j)%b_vec=b(i,j,:)
        enddo;enddo
        publish_calls=publish_calls+1;evolution_ready=.false.
    end subroutine publish_state

    subroutine check_roundtrip_and_atomicity()
        real(KREAL),allocatable::w(:,:,:),h(:,:,:),b(:,:,:),w0(:,:,:),h0(:,:,:),b0(:,:,:)
        real(KREAL)::err
        allocate(w(4,nx,ny),h(nx,ny,Nvel),b(nx,ny,Nvel),w0(4,nx,ny),h0(nx,ny,Nvel),b0(nx,ny,Nvel))
        call pack_state(w,h,b);w0=w;h0=h;b0=b
        call publish_state(w,h,b,.true.);call pack_state(w,h,b)
        err=max(maxval(abs(w-w0)),maxval(abs(h-h0)),maxval(abs(b-b0)))
        if(err/=0.0_KREAL)call abort_native('pack/unpack round trip is not bitwise identical',411)
        w=w0+1.0_KREAL;h=h0+1.0_KREAL;b=b0+1.0_KREAL
        call publish_state(w,h,b,.false.);call pack_state(w,h,b)
        err=max(maxval(abs(w-w0)),maxval(abs(h-h0)),maxval(abs(b-b0)))
        if(err/=0.0_KREAL.or.evolution_ready)call abort_native('failed proposal changed accepted state/cache',412)
        atomic_test_passed=.true.;publish_calls=0
        deallocate(w,h,b,w0,h0,b0)
    end subroutine check_roundtrip_and_atomicity

    subroutine S2NativeInitializeState(native_state,cold_start)
        logical,intent(in),optional::native_state,cold_start
        real(KREAL)::dx,dy,sincx,sincy,x,y,rho,err,cfl,tstart,tend,tol
        real(KREAL)::tu(Nvel),tv(Nvel),tw(Nvel),raw_node(Nu),raw_weight(Nu),momerr
        integer::i,j,m,env_status,ru,rx,ry,rv
        character(len=8)::reject_env
        character(len=1024)::restart_path
        real(KREAL)::restart_time,ruv(Nvel),rvv(Nvel),rwt(Nvel)
        logical::ok,rect_enabled,from_ctr,cold,resumed
        integer(int64)::imported_steps
        if(.not.package_mode())return
        call cpu_time(tstart);init_calls=init_calls+1
        if(IsS2CN2())then
            call CartesianSetup()
        elseif(NP_XY/=1)then
            call abort_native('non-S2CN2 research profile is single rank',401)
        endif
        reject_env=''
        call get_environment_variable('S2_RECT_WALL',reject_env,status=env_status)
        rect_enabled=S2CartesianMesh.or.(env_status==0.and.trim(reject_env)=='1')
        if(Nvel/=VelocityCount()**2.or.Nu/=VelocityCount().or.Nv/=VelocityCount().or.(ck/=1.and.S2PhysicalGas==0)) &
            call abort_native('velocity table count/CK mismatch',402)
        if(S2W1_wall_mode==1)then
            if(any(BC_Type(1:2)/=BCconst_periodic).or.any(BC_Type(3:4)/=BCconst_wall_isothermal)) &
                call abort_native('S2W1 boundary mode mismatch',403)
        elseif(S2W1_wall_mode==2)then
            if(any(BC_Type/=BCconst_wall_isothermal))call abort_native('cavity requires four isothermal walls',403)
            if(.not.cart_native)call abort_native('cavity requires native Cartesian adapter',403)
        else
            if(any(BC_Type/=BCconst_periodic))call abort_native('all four boundaries must be periodic',403)
        endif
        if((TimeStep_method/=Time_global_fixed.and..not.UnifiedPhysicalOperators).or..not.ieee_is_finite(dt_global).or.dt_global<=0.0_KREAL) &
            call abort_native('requires positive fixed dt',404)
        if(Iflag_DIS/=0.or.Iflag_Fringe/=0)call abort_native('DIS and fringe must be disabled',405)
        dx=ctr(1,1)%length(1);dy=ctr(1,1)%length(2)
        tol=1.0e-14_KREAL
        if(.not.UnifiedPhysicalOperators.and.abs(dx-2.0_KREAL*pi_d/real(Nx_all,KREAL))>tol) &
            call abort_native('requires uniform x Cartesian mesh',406)
        if(S2W1_wall_mode>=1)then
            if(abs(dy-1.0_KREAL/real(Ny_all,KREAL))>tol) &
                call abort_native('S2W1 requires uniform y in [0,1]',406)
        else
            if(.not.UnifiedPhysicalOperators.and.abs(dy-2.0_KREAL*pi_d/real(Ny_all,KREAL))>tol)call abort_native('requires uniform [0,2pi]^2 Cartesian mesh',406)
        endif
        do j=1,ny+1;do i=1,nx+1
            if(abs(gridx(i,j)-real(i-1+Block_XY(myID)%Origin_I-1,KREAL)*dx)>tol.or. &
               abs(gridy(i,j)-real(j-1+Block_XY(myID)%Origin_J-1,KREAL)*dy)>tol) &
                call abort_native('Cartesian coordinates or periodic seam mismatch',406)
        enddo;enddo
        do j=1,ny;do i=1,nx
            if(abs(ctr(i,j)%length(1)-dx)>tol.or.abs(ctr(i,j)%length(2)-dy)>tol.or. &
                abs(ctr(i,j)%area-dx*dy)>tol.or.abs(ctr(i,j)%x-(real(i+Block_XY(myID)%Origin_I-1,KREAL)-0.5_KREAL)*dx)>tol.or. &
                abs(ctr(i,j)%y-(real(j+Block_XY(myID)%Origin_J-1,KREAL)-0.5_KREAL)*dy)>tol) &
                call abort_native('nonuniform/non-Cartesian cell detected',407)
        enddo;enddo
        do j=1,ny;do i=1,nx+1
            if(abs(vface(i,j)%length-dy)>tol.or.abs(vface(i,j)%cosx-1.0_KREAL)>tol.or. &
               abs(vface(i,j)%cosy)>tol)call abort_native('invalid vertical face geometry',407)
        enddo;enddo
        do j=1,ny+1;do i=1,nx
            if(abs(hface(i,j)%length-dx)>tol.or.abs(hface(i,j)%cosx)>tol.or. &
               abs(hface(i,j)%cosy-1.0_KREAL)>tol)call abort_native('invalid horizontal face geometry',407)
        enddo;enddo
        if(len_trim(StructuredVelocityFile)>0)then
            ! InitVelocityFromStructuredFile already validates and broadcasts the
            ! selected quadrature. All methods use that same authoritative table.
            ! Do not regenerate exp/sqrt weights and demand cross-compiler identity.
            tu=uSpace;tv=vSpace;tw=weight
        else
            ! Retain the frozen generated-table contract when no file is selected.
            call init_velocity(tu,tv,tw,raw_node,raw_weight)
        endif
        if(.not.all(ieee_is_finite(uVel)).or..not.all(ieee_is_finite(vVel)).or. &
            .not.all(ieee_is_finite(wVel)).or..not.all(wVel>0.0_KREAL)) &
            call abort_native('invalid native velocity table',408)
        err=max(maxval(abs(tu-uVel)),maxval(abs(tv-vVel)),maxval(abs(tw-wVel)))
        if(err/=0.0_KREAL)call abort_native('native velocity arrays differ from selected source',408)
        cfl=dt_global*(maxval(abs(uVel))/dx+maxval(abs(vVel))/dy)
        if(cfl>0.4_KREAL)call abort_native('transport CFL exceeds 0.4',409)
        allocate(cw(4,nx,ny),ch(nx,ny,Nvel),cb(nx,ny,Nvel),mh(nx,ny,Nvel),mb(nx,ny,Nvel), &
            kxh(nx,ny,Nvel),kxb(nx,ny,Nvel),kyh(nx,ny,Nvel),kyb(nx,ny,Nvel), &
            eah(nx,ny,Nvel),eab(nx,ny,Nvel),kappa(2,2,nx,ny),fxw(4,nx,ny),fyw(4,nx,ny), &
            fxh(nx,ny,Nvel),fxb(nx,ny,Nvel),fyh(nx,ny,Nvel),fyb(nx,ny,Nvel))
        from_ctr=.false.;cold=.false.
        if(present(native_state))from_ctr=native_state
        if(present(cold_start))cold=cold_start
        if(from_ctr)call pack_state(cw,ch,cb)
        sincx=sin(0.5_KREAL*dx)/(0.5_KREAL*dx);sincy=sin(0.5_KREAL*dy)/(0.5_KREAL*dy)
        do j=1,ny;y=ctr(1,j)%y
            do i=1,nx;x=ctr(i,j)%x
                if(.not.from_ctr)then
                select case(S2D1_initial_mode)
                case(0);rho=rho0+amp*sincx*sincy*sin(x+y)
                case(3);rho=rho0
                case(4,5);rho=rho0
                case default;call abort_native('native package supports initial modes 0 and 3 only',410)
                end select
                if(S2W1_wall_mode>=1)then
                    x=0.0_KREAL
                    if(S2D1_initial_mode==5)x=0.1_KREAL*(y-0.5_KREAL)
                    cw(:,i,j)=[rho,rho*x,0.0_KREAL,rho*(0.5_KREAL*x*x+0.5_KREAL*real(ck+2,KREAL))]
                else
                    cw(:,i,j)=[rho,rho*u0,rho*v0,rho*(0.5_KREAL*(u0*u0+v0*v0)+0.5_KREAL*real(ck+2,KREAL)*t0)]
                endif
                endif
                if(.not.from_ctr.or.cold)then
                call make_mke(cw(:,i,j),uVel,vVel,wVel,mh(i,j,:),mb(i,j,:),kxh(i,j,:),kxb(i,j,:), &
                    kyh(i,j,:),kyb(i,j,:),eah(i,j,:),eab(i,j,:),kappa(:,:,i,j),ok,momerr)
                if(.not.ok)call abort_native('discrete conservative initialization failed',413)
                ch(i,j,:)=mh(i,j,:);cb(i,j,:)=mb(i,j,:)
                endif
            enddo
        enddo
        if(len_trim(S2CheckpointImport)>0)then
            if(.not.cold)call abort_native('native restart and old checkpoint import are mutually exclusive',456)
            if(NP_XY/=1)call abort_native('legacy checkpoint import currently requires one rank; then native repartition',456)
            slab_lo=1;slab_hi=ny
            legacy_restart_path=S2CheckpointImport
            call load_input(cw,ch,cb,uVel,vVel,wVel,time_dis,imported_steps,total0,wall_balance, &
                cumulative_bottom,cumulative_top,resumed)
            if(.not.resumed)call abort_native('old checkpoint did not resume',457)
            if(imported_steps>=int(MAX_ITER,int64))call abort_native('checkpoint step must be below target MAX_ITER',460)
            iter=int(imported_steps)
        endif
        call native_build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,err,ok)
        if(.not.ok)call abort_native('initial volume target construction failed',414)
        targets_ready=.true.
        call publish_state(cw,ch,cb,.true.);publish_calls=0
        do j=1,ny;do i=1,nx
            ctr(i,j)%dt_local=dt_global;ctr(i,j)%conVars0=ctr(i,j)%conVars
            ctr(i,j)%h0_vec=ctr(i,j)%h_vec;ctr(i,j)%b0_vec=ctr(i,j)%b_vec
            ctr(i,j)%res_conVars0=0.0_KREAL;ctr(i,j)%res_h0_vec=0.0_KREAL;ctr(i,j)%res_b0_vec=0.0_KREAL
        enddo;enddo
        if(len_trim(S2CheckpointImport)==0)then
            total0=0.0_KREAL
            do i=1,4;total0(i)=sum(cw(i,:,:));enddo
        endif
        if(cart_native)then
            total0=GlobalSum4(total0)
            call CartesianLedger(total0,wall_balance,cumulative_bottom,cumulative_top)
        endif
        call check_roundtrip_and_atomicity()
        restart_path=''
        call get_environment_variable('S2DO1_RESTART_PATH',restart_path,status=env_status)
        if(.not.from_ctr.and.env_status==0.and.len_trim(restart_path)>0)then
            open(newunit=ru,file=trim(restart_path),status='old',access='stream',form='unformatted',action='read')
            read(ru)rx,ry,rv,restart_time
            if(rx/=nx.or.ry/=ny.or.rv/=Nvel)call abort_native('O1 restart dimensions disagree',431)
            read(ru)cw;read(ru)ch;read(ru)cb;read(ru)ruv;read(ru)rvv;read(ru)rwt
            close(ru)
            if(max(maxval(abs(ruv-uVel)),maxval(abs(rvv-vVel)),maxval(abs(rwt-wVel)))/=0.0_KREAL) &
                call abort_native('O1 restart velocity table disagrees',432)
            call native_build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,err,ok)
            if(.not.ok)call abort_native('O1 restart target construction failed',433)
            call publish_state(cw,ch,cb,.true.);publish_calls=0
            do j=1,ny;do i=1,nx
                ctr(i,j)%conVars0=ctr(i,j)%conVars
                ctr(i,j)%h0_vec=ctr(i,j)%h_vec;ctr(i,j)%b0_vec=ctr(i,j)%b_vec
            enddo;enddo
            total0=0.0_KREAL
            do i=1,4;total0(i)=sum(cw(i,:,:));enddo
            open(newunit=ru,file='o1_restart_origin.txt',status='replace',action='write')
            write(ru,'(A)')trim(restart_path)
            write(ru,'(ES24.16)')restart_time
            close(ru)
        endif
        reject_env=''
        call get_environment_variable('S2DO1_BENCH',reject_env,status=env_status)
        o1_bench=env_status==0.and.trim(reject_env)=='1'
        reject_env=''
        call get_environment_variable('S2W2_COST',reject_env,status=env_status)
        w2_cost_enabled=env_status==0.and.trim(reject_env)=='1'
        if(w2_cost_enabled)then
            open(newunit=w2_cost_unit,file='s2w2_step_cost.csv',status='replace',action='write')
            write(w2_cost_unit,'(A)')'step,official_micro_reconstruction_cpu_s,official_wall_cpu_s,wall_node_calls'
        endif
        call get_environment_variable('S2W2R_METRICS',reject_env,status=env_status)
        w2r_metrics=env_status==0.and.trim(reject_env)=='1'
        if(w2r_metrics)then
            open(newunit=w2r_unit,file='s2w2r_step_cost.csv',status='replace',action='write')
            write(w2r_unit,'(A)')'step,enabled,native_wall_s,old_reconstruction_cpu_s,old_closure_cpu_s,'// &
                'old_wall_cpu_s,predict_cpu_s,instant_reconstruction_cpu_s,apply_cpu_s,cn_cpu_s,'// &
                'immutable_maxabs,pred_moment_error,pred_h_min,pred_b_min,pred_limited,pred_signed,wall_mass,'// &
                'accepted_calls,halfstep_calls,old_recon_calls,star_recon_calls,old_target_reuses,new_target_builds'
        endif
        call get_environment_variable('S2W1D_DIAGNOSTIC' ,reject_env,status=env_status)
        w1d_diagnostic=env_status==0.and.trim(reject_env)=='1'
        if(w1d_diagnostic)then
            if(S2W1_wall_mode/=1)call abort_native('W1D diagnostics require flat diffuse walls',450)
            open(newunit=w1d_audit_unit,file='w1d_predictor_audit.csv',status='replace',action='write')
            write(w1d_audit_unit,'(A)')'step,immutable_max_difference,half_limited,half_signed,'// &
                'mid_wall_limited,mid_wall_signed,half_target_moment,half_cn_residual'
        endif
        reject_env=''
        call get_environment_variable('S2I1_DIAGNOSTIC',reject_env,status=env_status)
        i1_diag_enabled=IsS2Steady().or.(env_status==0.and.trim(reject_env)=='1')
        if(i1_diag_enabled)then
            do m=1,4
                i1_scale(m)=max(1.0_KREAL,maxval(abs(cw(m,:,:))))
            enddo
            i1_scale(5)=maxval(abs(ch));i1_scale(6)=maxval(abs(cb))
            if(minval(i1_scale)<=0.0_KREAL)call abort_native('I1 initial scale is zero',446)
            i1_evaluations=0;i1_consecutive=0;i1_below=.false.;i1_last_state_iter=-1
            open(newunit=i1_unit,file='s2i1_true_residual.csv',status='replace',action='write')
            write(i1_unit,'(A)')'accepted_iter,eval,wall_s,W1_L2,W2_L2,W3_L2,W4_L2,h_L2,b_L2,'// &
                'W1_Linf,W2_Linf,W3_Linf,W4_Linf,h_Linf,b_Linf,max_relative,consecutive,unique_state_eval'
            i1_start_wall=MPI_Wtime()
        endif
        if((S2Diagnostics.or..not.IsS2CN2()).and..not.cart_native)then
        call DumpLocalAccepted(cw,ch,cb,time_dis,'initial')
        open(newunit=steps_unit,file='native_steps.csv',status='replace',action='write')
        if(S2PhysicalGas==1)then
            open(newunit=gas_unit,file='physical_gas.csv',status='replace',action='write')
            write(gas_unit,'(A)')'step,body_tau_min,face_tau_min,wall_tau_min,body_tau_max,face_tau_max,wall_tau_max,body_calls,face_calls,wall_calls'
        endif
        steps_open=.true.
        write(steps_unit,'(A)')'step,time,face_moment,target_moment,cn_residual,macro_micro,limited,signed,'// &
            'dt_tau,cfl,h_min,h_max,b_min,b_max,global_drift,layout_s,face_s,map_s,collision_s,diagnostic_s'
        endif
        if(S2W1_wall_mode>=1)then
            allocate(lower_w(4,nx),lower_h(nx,Nvel),lower_b(nx,Nvel))
            if((S2Diagnostics.or..not.IsS2CN2()).and..not.cart_native)then
            open(newunit=wall_unit,file='wall_balance.csv',status='replace',action='write')
            write(wall_unit,'(A)')'step,time,delta_mass,delta_mx,delta_my,delta_energy,'// &
                'bottom_mass,bottom_mx,bottom_my,bottom_energy,top_mass,top_mx,top_my,top_energy,'// &
                'balance_mass,balance_mx,balance_my,balance_energy,max_relative_error'
            open(newunit=wall_segment_unit,file='wall_segment.csv',status='replace',action='write')
            write(wall_segment_unit,'(A)')'step,i,side,mass,mx,my,energy,shear,heat'
            endif
            if(S2W1_time_centered==1)then
                open(newunit=w1t_cost_unit,file='s2w1t_step_cost.csv',status='replace',action='write')
                write(w1t_cost_unit,'(A)')'step,full_step_wall_s,official_face_s,predictor_total_s,'// &
                    'predictor_wall_only_s,official_collision_s,official_face_calls,predictor_half_calls,'// &
                    'predictor_wall_calls'
            endif
            if(S2W1D_face_time_consistent==1)then
                allocate(leading_xh(nx,ny,Nvel),leading_xb(nx,ny,Nvel), &
                    leading_yh(nx,ny,Nvel),leading_yb(nx,ny,Nvel),time_xh(nx,ny,Nvel), &
                    time_xb(nx,ny,Nvel),time_yh(nx,ny,Nvel),time_yb(nx,ny,Nvel))
                open(newunit=w1d_cost_unit,file='s2w1d_step_cost.csv',status='replace',action='write')
                write(w1d_cost_unit,'(A)')'step,instantaneous_face_s,instantaneous_face_calls,'// &
                    'instantaneous_limited,instantaneous_signed'
            endif
        endif
        reject_env=''
        call get_environment_variable('S2D3_TEST_GLOBAL_REJECT',reject_env,status=env_status)
        force_global_reject=IsS2CN2().and.env_status==0.and.trim(reject_env)=='1'
        if(S2W2R_shared_predictor==1)then
            w2r_active=1
            allocate(leading_xh(nx,ny,Nvel),leading_xb(nx,ny,Nvel),leading_yh(nx,ny,Nvel),leading_yb(nx,ny,Nvel), &
                time_xh(nx,ny,Nvel),time_xb(nx,ny,Nvel),time_yh(nx,ny,Nvel),time_yb(nx,ny,Nvel), &
                w2r_l0h(nx,Nvel),w2r_l0b(nx,Nvel),w2r_lkh(nx,Nvel),w2r_lkb(nx,Nvel))
            if(.not.w2r_metrics)call abort_native('W2R requires S2W2R_METRICS=1 for research admissibility record',455)
        endif
        if(cart_native.and.S2Diagnostics)call DumpLocalAccepted(cw,ch,cb,time_dis,'initial')
        initialized=.true.;evolution_ready=.false.;call cpu_time(tend);layout_seconds=layout_seconds+tend-tstart
    end subroutine S2NativeInitializeState


    subroutine M1DumpAudit(phase)
        integer,intent(in)::phase
        integer::st,fu,i,j
        character(len=8)::flag
        character(len=64)::path
        flag='';call get_environment_variable('M1_AUDIT',flag,status=st)
        if(st/=0.or.trim(flag)/='1'.or.update_calls/=0)return
        if(phase==0)then
            open(newunit=fu,file='m1_yfaces.bin',status='replace',access='stream',form='unformatted')
            do j=1,ny+1;do i=1,nx
                write(fu)hface(i,j)%L_flux/hface(i,j)%length, &
                    hface(i,j)%L_flux_h_vec/hface(i,j)%length,hface(i,j)%L_flux_b_vec/hface(i,j)%length
            enddo;enddo
            close(fu)
        endif
        write(path,'("m1_targets_",I1,".bin")')phase
        open(newunit=fu,file=trim(path),status='replace',access='stream',form='unformatted')
        write(fu)mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa
        close(fu)
    end subroutine M1DumpAudit

    subroutine S2NativeEvolution()
        real(KREAL)::dx,dy,t0c,t1c
        real(KREAL),allocatable::w1d_before(:),w1d_after(:),w2r_before(:),w2r_after(:)
        integer::i,j
        logical::ok
        if(.not.(IsS2CN2().or.IsS2Steady()))return
        if(.not.initialized.or.(evolution_ready.and..not.IsHighOrderImplicit()))call abort_native('invalid Evolution adapter state',421)
        if(cart_native)call CartesianLedger(total0,wall_balance,cumulative_bottom,cumulative_top)
        o1_wall_start=MPI_Wtime()
        call cpu_time(t0c);call pack_state(cw,ch,cb);call cpu_time(t1c);layout_seconds=layout_seconds+t1c-t0c
        if(IsHighOrderImplicit())targets_ready=.false.
        if(.not.targets_ready)then
            call native_build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,face_error,ok)
            if(.not.ok)call abort_native('accepted-state target rebuild failed',420)
            targets_ready=.true.
        endif
        dx=ctr(1,1)%length(1);dy=ctr(1,1)%length(2)
        call cpu_time(t0c)
        if(cart_native)then
            call CartesianEvolution(cw,ch,cb,fxw,fyw,fxh,fxb,fyh,fyb,lower_w,lower_h,lower_b,face_error,limited_count,signed_count,ok)
            evolution_calls=evolution_calls+1;evolution_ready=.true.
            if(S2W1_wall_mode==2)call PublishWallDiagnostics()
            return
        endif
        if(S2W2R_shared_predictor==1)then
            call compute_flux(cw,ch,cb,uVel,vVel,wVel,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb, &
                face_error,limited_count,signed_count,ok,lower_w,lower_h,lower_b, &
                leading_xh,leading_xb,leading_yh,leading_yb,time_xh,time_xb,time_yh,time_yb, &
                w2r_l0h,w2r_l0b,w2r_lkh,w2r_lkb)
        elseif(S2W1D_face_time_consistent==1)then
            call compute_flux(cw,ch,cb,uVel,vVel,wVel,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb, &
                face_error,limited_count,signed_count,ok,lower_w,lower_h,lower_b, &
                leading_xh,leading_xb,leading_yh,leading_yb,time_xh,time_xb,time_yh,time_yb)
        elseif(S2W1_wall_mode>=1)then
            call compute_flux(cw,ch,cb,uVel,vVel,wVel,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb, &
                face_error,limited_count,signed_count,ok,lower_w,lower_h,lower_b)
        else
            call compute_flux(cw,ch,cb,uVel,vVel,wVel,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb, &
                face_error,limited_count,signed_count,ok)
        endif
        w2r_predict_s=0;w2r_instant_s=0;w2r_apply_s=0;w2r_diag=0;w2r_immutable_error=0
        w2r_off_closure_s=w2_closure_s
        w2_off_reconstruction_s=w2_reconstruction_s;w2_off_wall_s=w2_wall_s;w2_off_wall_calls=w2_wall_calls
        call cpu_time(t1c);w1t_official_face_s=t1c-t0c;face_seconds=face_seconds+w1t_official_face_s
        if(.not.ok)then
            if(IsS2Steady())call write_state('s2i1_face_rejected_last_accepted.bin', &
                cw,ch,cb,uVel,vVel,wVel,0.0_KREAL)
            call abort_native('face kernel rejected state',422)
        endif
        if(w2r_metrics)w2r_before=w2r_snapshot()
        if(S2W2R_shared_predictor==1)then
            call s2w2r_correct(cw,ch,cb,mh,mb,uVel,vVel,wVel,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb, &
                lower_w,lower_h,lower_b,leading_xh,leading_xb,leading_yh,leading_yb, &
                time_xh,time_xb,time_yh,time_yb,w2r_l0h,w2r_l0b,w2r_lkh,w2r_lkb, &
                w2r_diag,w2r_predict_s,w2r_instant_s,w2r_apply_s,ok)
            if(.not.ok)call abort_native('W2R explicit shared predictor or corrected face rejected',456)
        endif
        if(w2r_metrics)then
            w2r_after=w2r_snapshot();w2r_immutable_error=maxval(abs(w2r_after-w2r_before))
            if(w2r_immutable_error/=0.0_KREAL)call abort_native('W2R predictor changed accepted state/cache/ledger',454)
        endif
        if(S2W1_time_centered==1)then
            if(w1d_diagnostic)w1d_before=w1d_immutable_snapshot()
            call S2W1PredictMidpointWall(dx,dy,ok)
            if(.not.ok)then
                call write_state('native_predictor_rejected_accepted.bin',cw,ch,cb,uVel,vVel,wVel,time_dis)
                call abort_native('isolated midpoint wall predictor rejected state',449)
            endif
            if(w1d_diagnostic)then
                w1d_after=w1d_immutable_snapshot()
                if(any(w1d_after/=w1d_before))call abort_native('predictor altered formal input/interior faces',451)
                write(w1d_audit_unit,'(I0,7(",",ES24.16))')evolution_calls+1, &
                    maxval(abs(w1d_after-w1d_before)),w1d_prediction_diag(5:6), &
                    real(w1d_mid_limited,KREAL),real(w1d_mid_signed,KREAL),w1d_prediction_diag(2:3)
                flush(w1d_audit_unit)
            endif
        endif
        call cpu_time(t0c)
        do j=1,ny
            do i=1,nx-1
                vface(i+1,j)%L_flux=fxw(:,i,j)*vface(i+1,j)%length
                vface(i+1,j)%R_flux=vface(i+1,j)%L_flux
                vface(i+1,j)%L_flux_h_vec=fxh(i,j,:)*vface(i+1,j)%length
                vface(i+1,j)%R_flux_h_vec=vface(i+1,j)%L_flux_h_vec
                vface(i+1,j)%L_flux_b_vec=fxb(i,j,:)*vface(i+1,j)%length
                vface(i+1,j)%R_flux_b_vec=vface(i+1,j)%L_flux_b_vec
            enddo
            vface(1,j)%L_flux=fxw(:,nx,j)*vface(1,j)%length;vface(1,j)%R_flux=vface(1,j)%L_flux
            vface(nx+1,j)%L_flux=vface(1,j)%L_flux;vface(nx+1,j)%R_flux=vface(1,j)%L_flux
            vface(1,j)%L_flux_h_vec=fxh(nx,j,:)*vface(1,j)%length
            vface(1,j)%R_flux_h_vec=vface(1,j)%L_flux_h_vec
            vface(nx+1,j)%L_flux_h_vec=vface(1,j)%L_flux_h_vec;vface(nx+1,j)%R_flux_h_vec=vface(1,j)%L_flux_h_vec
            vface(1,j)%L_flux_b_vec=fxb(nx,j,:)*vface(1,j)%length
            vface(1,j)%R_flux_b_vec=vface(1,j)%L_flux_b_vec
            vface(nx+1,j)%L_flux_b_vec=vface(1,j)%L_flux_b_vec;vface(nx+1,j)%R_flux_b_vec=vface(1,j)%L_flux_b_vec
        enddo
        do i=1,nx
            do j=1,ny-1
                hface(i,j+1)%L_flux=fyw(:,i,j)*hface(i,j+1)%length
                hface(i,j+1)%R_flux=hface(i,j+1)%L_flux
                hface(i,j+1)%L_flux_h_vec=fyh(i,j,:)*hface(i,j+1)%length
                hface(i,j+1)%R_flux_h_vec=hface(i,j+1)%L_flux_h_vec
                hface(i,j+1)%L_flux_b_vec=fyb(i,j,:)*hface(i,j+1)%length
                hface(i,j+1)%R_flux_b_vec=hface(i,j+1)%L_flux_b_vec
            enddo
            if(S2W1_wall_mode>=1)then
                hface(i,1)%L_flux=lower_w(:,i)*hface(i,1)%length
                hface(i,ny+1)%L_flux=fyw(:,i,ny)*hface(i,ny+1)%length
                hface(i,1)%L_flux_h_vec=lower_h(i,:)*hface(i,1)%length
            else
                hface(i,1)%L_flux=fyw(:,i,ny)*hface(i,1)%length
                hface(i,ny+1)%L_flux=hface(i,1)%L_flux
                hface(i,1)%L_flux_h_vec=fyh(i,ny,:)*hface(i,1)%length
            endif
            hface(i,1)%R_flux=hface(i,1)%L_flux
            hface(i,ny+1)%R_flux=hface(i,ny+1)%L_flux
            hface(i,1)%R_flux_h_vec=hface(i,1)%L_flux_h_vec
            hface(i,ny+1)%L_flux_h_vec=fyh(i,ny,:)*hface(i,ny+1)%length
            hface(i,ny+1)%R_flux_h_vec=hface(i,ny+1)%L_flux_h_vec
            if(S2W1_wall_mode>=1)then
                hface(i,1)%L_flux_b_vec=lower_b(i,:)*hface(i,1)%length
            else
                hface(i,1)%L_flux_b_vec=fyb(i,ny,:)*hface(i,1)%length
            endif
            hface(i,1)%R_flux_b_vec=hface(i,1)%L_flux_b_vec
            hface(i,ny+1)%L_flux_b_vec=fyb(i,ny,:)*hface(i,ny+1)%length
            hface(i,ny+1)%R_flux_b_vec=hface(i,ny+1)%L_flux_b_vec
        enddo
        call cpu_time(t1c);map_seconds=map_seconds+t1c-t0c
        call M1DumpAudit(0)
        evolution_calls=evolution_calls+1;evolution_ready=.true.
    end subroutine S2NativeEvolution

    ! An isolated W1 half step predicts the accepted state at t^n+dt/2.
    ! Its transport uses both interior UGKS faces and the original diffuse wall;
    ! advance_step applies the actual CN collision. No ctr, time, accepted target,
    ! face ledger, or native residual cache is changed by this temporary step.
    subroutine S2W1PredictMidpointWall(dx,dy,ok)
        real(KREAL),intent(in)::dx,dy
        logical,intent(out)::ok
        real(KREAL),allocatable::pw(:,:,:),ph(:,:,:),pb(:,:,:),sw(:,:,:),sh(:,:,:),sb(:,:,:)
        real(KREAL),allocatable::pmh(:,:,:),pmb(:,:,:),pkxh(:,:,:),pkxb(:,:,:),pkyh(:,:,:),pkyb(:,:,:)
        real(KREAL),allocatable::peah(:,:,:),peab(:,:,:),pkap(:,:,:,:)
        real(KREAL),allocatable::pfxw(:,:,:),pfyw(:,:,:),pfxh(:,:,:),pfxb(:,:,:),pfyh(:,:,:),pfyb(:,:,:)
        real(KREAL),allocatable::plw(:,:),plh(:,:),plb(:,:),mlw(:,:),mlh(:,:),mlb(:,:)
        real(KREAL),allocatable::muw(:,:),muh(:,:),mub(:,:)
        real(KREAL),allocatable::mfxh(:,:,:),mfxb(:,:,:),mfyh(:,:,:),mfyb(:,:,:)
        real(KREAL)::half,face_err,diag(12),flux_s,collision_s,tstart,twall
        integer::i,j,im,limited,signed,mid_limited,mid_signed
        logical::lok
        ok=.false.;half=0.5_KREAL*dt_global;tstart=MPI_Wtime()
        allocate(pw(4,nx,ny),ph(nx,ny,Nvel),pb(nx,ny,Nvel),sw(4,nx,ny),sh(nx,ny,Nvel),sb(nx,ny,Nvel))
        allocate(pmh(nx,ny,Nvel),pmb(nx,ny,Nvel),pkxh(nx,ny,Nvel),pkxb(nx,ny,Nvel), &
            pkyh(nx,ny,Nvel),pkyb(nx,ny,Nvel),peah(nx,ny,Nvel),peab(nx,ny,Nvel),pkap(2,2,nx,ny))
        allocate(pfxw(4,nx,ny),pfyw(4,nx,ny),pfxh(nx,ny,Nvel),pfxb(nx,ny,Nvel), &
            pfyh(nx,ny,Nvel),pfyb(nx,ny,Nvel),plw(4,nx),plh(nx,Nvel),plb(nx,Nvel))
        allocate(mlw(4,nx),mlh(nx,Nvel),mlb(nx,Nvel),muw(4,nx),muh(nx,Nvel),mub(nx,Nvel))
        pw=cw;ph=ch;pb=cb;pmh=mh;pmb=mb;pkxh=kxh;pkxb=kxb;pkyh=kyh;pkyb=kyb
        peah=eah;peab=eab;pkap=kappa
        call compute_flux(cw,ch,cb,uVel,vVel,wVel,dx,dy,half,pfxw,pfyw,pfxh,pfxb,pfyh,pfyb, &
            face_err,limited,signed,lok,plw,plh,plb)
        if(.not.lok)return
        do j=1,ny
            do i=1,nx
                im=modulo(i-2,nx)+1
                sw(:,i,j)=cw(:,i,j)-(pfxw(:,i,j)-pfxw(:,im,j))/dx
                sh(i,j,:)=ch(i,j,:)-(pfxh(i,j,:)-pfxh(im,j,:))/dx
                sb(i,j,:)=cb(i,j,:)-(pfxb(i,j,:)-pfxb(im,j,:))/dx
                if(j==1)then
                    sw(:,i,j)=sw(:,i,j)-(pfyw(:,i,j)-plw(:,i))/dy
                    sh(i,j,:)=sh(i,j,:)-(pfyh(i,j,:)-plh(i,:))/dy
                    sb(i,j,:)=sb(i,j,:)-(pfyb(i,j,:)-plb(i,:))/dy
                else
                    sw(:,i,j)=sw(:,i,j)-(pfyw(:,i,j)-pfyw(:,i,j-1))/dy
                    sh(i,j,:)=sh(i,j,:)-(pfyh(i,j,:)-pfyh(i,j-1,:))/dy
                    sb(i,j,:)=sb(i,j,:)-(pfyb(i,j,:)-pfyb(i,j-1,:))/dy
                endif
            enddo
        enddo
        call advance_step(pw,ph,pb,uVel,vVel,wVel,dx,dy,half,pmh,pmb,pkxh,pkxb,pkyh,pkyb, &
            peah,peab,pkap,diag,flux_s,collision_s,lok,sw,sh,sb,face_err,limited,signed)
        if(.not.lok)return
        if(w1d_diagnostic)w1d_prediction_diag=diag
        twall=MPI_Wtime()
        call s2w1_compute_wall_flux(ph,pb,uVel,vVel,wVel,dx,dy,dt_global, &
            mlw,mlh,mlb,muw,muh,mub,mid_limited,mid_signed,lok)
        w1t_predictor_wall_s=MPI_Wtime()-twall
        if(.not.lok)return
        if(w1d_diagnostic)then
            w1d_mid_limited=mid_limited;w1d_mid_signed=mid_signed
        endif
        lower_w=mlw;lower_h=mlh;lower_b=mlb
        fyw(:,:,ny)=muw;fyh(:,ny,:)=muh;fyb(:,ny,:)=mub
        if(S2W1D_face_time_consistent==1)then
            allocate(mfxh(nx,ny,Nvel),mfxb(nx,ny,Nvel),mfyh(nx,ny,Nvel),mfyb(nx,ny,Nvel))
            twall=MPI_Wtime()
            call s2w1d_instantaneous_flux(ph,pb,uVel,vVel,dx,dy,mfxh,mfxb,mfyh,mfyb, &
                w1d_instant_limited,w1d_instant_signed,lok)
            if(.not.lok)return
            fxh=fxh+dt_global*(mfxh-leading_xh)-0.5_KREAL*dt_global**2*time_xh
            fxb=fxb+dt_global*(mfxb-leading_xb)-0.5_KREAL*dt_global**2*time_xb
            fyh(:,1:ny-1,:)=fyh(:,1:ny-1,:)+dt_global*(mfyh(:,1:ny-1,:)-leading_yh(:,1:ny-1,:)) &
                -0.5_KREAL*dt_global**2*time_yh(:,1:ny-1,:)
            fyb(:,1:ny-1,:)=fyb(:,1:ny-1,:)+dt_global*(mfyb(:,1:ny-1,:)-leading_yb(:,1:ny-1,:)) &
                -0.5_KREAL*dt_global**2*time_yb(:,1:ny-1,:)
            do j=1,ny;do i=1,nx
                call moment_2d(fxh(i,j,:),fxb(i,j,:),uVel,vVel,wVel,fxw(:,i,j),lok)
                if(.not.lok)return
                if(j<ny)then
                    call moment_2d(fyh(i,j,:),fyb(i,j,:),uVel,vVel,wVel,fyw(:,i,j),lok)
                    if(.not.lok)return
                endif
            enddo;enddo
            w1d_instant_s=MPI_Wtime()-twall
        endif
        w1t_predictor_step_s=MPI_Wtime()-tstart
        ok=.true.
    end subroutine S2W1PredictMidpointWall

    ! Research-only snapshots of data subsequently consumed by the formal step.
    ! Boundary slots are excluded because replacing them is the intended W1T action.
    function w1d_immutable_snapshot() result(data)
        real(KREAL),allocatable::data(:)
        data=[reshape(cw,[size(cw)]),reshape(ch,[size(ch)]),reshape(cb,[size(cb)]), &
            reshape(mh,[size(mh)]),reshape(mb,[size(mb)]),reshape(kxh,[size(kxh)]), &
            reshape(kxb,[size(kxb)]),reshape(kyh,[size(kyh)]),reshape(kyb,[size(kyb)]), &
            reshape(eah,[size(eah)]),reshape(eab,[size(eab)]),reshape(kappa,[size(kappa)])]
        if(S2W1D_face_time_consistent==0)data=[data, &
            reshape(fxw,[size(fxw)]),reshape(fxh,[size(fxh)]),reshape(fxb,[size(fxb)]), &
            reshape(fyw(:,:,1:ny-1),[4*nx*(ny-1)]), &
            reshape(fyh(:,1:ny-1,:),[nx*(ny-1)*Nvel]),reshape(fyb(:,1:ny-1,:),[nx*(ny-1)*Nvel])]
        if(S2W1D_face_time_consistent==1)data=[data, &
            reshape(leading_xh,[size(leading_xh)]),reshape(leading_xb,[size(leading_xb)]), &
            reshape(leading_yh,[size(leading_yh)]),reshape(leading_yb,[size(leading_yb)]), &
            reshape(time_xh,[size(time_xh)]),reshape(time_xb,[size(time_xb)]), &
            reshape(time_yh,[size(time_yh)]),reshape(time_yb,[size(time_yb)])]
    end function w1d_immutable_snapshot

    subroutine w1d_dump_pre(step,wn,hs,bs)
        integer,intent(in)::step
        real(KREAL),intent(in)::wn(:,:,:),hs(:,:,:),bs(:,:,:)
        character(len=64)::path
        integer::du
        write(path,'("w1d_step",I3.3,"_pre.bin")')step
        open(newunit=du,file=trim(path),status='replace',access='stream',form='unformatted')
        write(du)int(nx,KINT),int(ny,KINT),int(Nvel,KINT),int(S2W1_time_centered,KINT)
        write(du)dt_global,ctr(1,1)%length(1),ctr(1,1)%length(2),time_dis
        write(du)ctr(1,1)%area,hface(1,1)%length,vface(1,1)%length
        write(du)cw,ch,cb,mh,mb,eah,eab
        write(du)fxw,fyw,fxh,fxb,fyh,fyb,lower_w,lower_h,lower_b
        write(du)res_macro_flux_store,res_h_flux_store,res_b_flux_store,wn,hs,bs
        close(du)
    end subroutine w1d_dump_pre

    subroutine w1d_dump_post(step,w,h,b,diag)
        integer,intent(in)::step
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),diag(12)
        character(len=64)::path
        integer::du
        write(path,'("w1d_step",I3.3,"_post.bin")')step
        open(newunit=du,file=trim(path),status='replace',access='stream',form='unformatted')
        write(du)int(nx,KINT),int(ny,KINT),int(Nvel,KINT),int(S2W1_time_centered,KINT)
        write(du)dt_global,time_dis+dt_global,diag
        write(du)w,h,b,mh,mb,eah,eab
        close(du)
    end subroutine w1d_dump_post

    subroutine S2NativePhysicalCollision(eqh,eqb,colh,colb)
        real(KREAL),intent(out)::eqh(:,:,:),eqb(:,:,:),colh(:,:,:),colb(:,:,:)
        real(KREAL)::q(2),tau,mom(4)
        integer::i,j
        logical::ok
        if(.not.targets_ready.or..not.evolution_ready)call abort_native('stale S2 collision target',441)
        do j=1,ny;do i=1,nx
            call heat(cw(:,i,j),ch(i,j,:)-mh(i,j,:),cb(i,j,:)-mb(i,j,:),uVel,vVel,wVel,q,ok)
            if(.not.ok)call abort_native('S2 physical heat closure',443)
            tau=s2_tau(cw(:,i,j))
            eqh(:,i,j)=mh(i,j,:)+kxh(i,j,:)*q(1)+kyh(i,j,:)*q(2)
            eqb(:,i,j)=mb(i,j,:)+kxb(i,j,:)*q(1)+kyb(i,j,:)*q(2)
            colh(:,i,j)=(eqh(:,i,j)-ch(i,j,:))/tau
            colb(:,i,j)=(eqb(:,i,j)-cb(i,j,:))/tau
            mom=[sum(wVel*colh(:,i,j)),sum(wVel*uVel*colh(:,i,j)),sum(wVel*vVel*colh(:,i,j)), &
                0.5_KREAL*sum(wVel*((uVel*uVel+vVel*vVel)*colh(:,i,j)+colb(:,i,j)))]
            ! S2's established gate is on equilibrium conserved moments. A collision
            ! rate divides this error by tau; check the original dimensionless quantity.
            if(tau*maxval(abs(mom))/max(1.0_KREAL,maxval(abs(cw(:,i,j))))>gate)then
                write(*,'(A,I0,2I6,10ES24.16)')'collision diagnostic rank/cell/tau/mom/hgap/bgap:',myID,i,j,tau,mom, &
                    maxval(abs(eqh(:,i,j)-ch(i,j,:))),maxval(abs(eqb(:,i,j)-cb(i,j,:)))
                call abort_native('S2 physical collision moment gate',444)
            endif
        enddo;enddo
    end subroutine

    subroutine S2NativeMarkResidual()
        if(.not.(IsS2CN2().or.IsS2Steady()))return
        if(.not.evolution_ready)call abort_native('residual assembled without fresh Evolution cache',423)
        residual_calls=residual_calls+1
    end subroutine S2NativeMarkResidual

    subroutine S2NativeAssembleSteadyResidual()
        real(KREAL)::q(2),mom(4),crh(Nvel),crb(Nvel),rh(Nvel),rb(Nvel),rw(4),invdt,scale
        real(KREAL)::sums(6),linf(6),norm(6),relative(6),area_total,weight_norm
        real(KREAL)::qdump(2,nx,ny)
        integer::i,j,k,m,env_status,fu
        character(len=8)::diag_env
        logical::ok,counted,valid
        if(.not.(IsS2Steady().or.i1_diag_enabled))return
        if(.not.evolution_ready.or.residual_calls/=evolution_calls) &
            call abort_native('S2 steady residual lacks fresh high-order face',440)
        if(.not.targets_ready)call abort_native('S2 steady body target is stale',441)
        if(maxval(abs(res_macro_store-res_macro_flux_store))/=0.0_KREAL) &
            call abort_native('unsupported source in S2 steady residual',442)
        invdt=1.0_KREAL/dt_global
        sums=0.0_KREAL;linf=0.0_KREAL;area_total=0.0_KREAL
        weight_norm=sum(wVel)
        do j=1,ny;do i=1,nx
            rw=res_macro_flux_store(:,i,j)*invdt
            if(IsS2Steady())rate_w(:,i,j)=rw
            call heat(cw(:,i,j),ch(i,j,:)-mh(i,j,:),cb(i,j,:)-mb(i,j,:),uVel,vVel,wVel,q,ok)
            if(.not.ok)call abort_native('S2 steady heat closure failed',443)
            qdump(:,i,j)=q
            crh=(mh(i,j,:)+kxh(i,j,:)*q(1)+kyh(i,j,:)*q(2)-ch(i,j,:))/tau0
            crb=(mb(i,j,:)+kxb(i,j,:)*q(1)+kyb(i,j,:)*q(2)-cb(i,j,:))/tau0
            mom(1)=sum(wVel*crh)
            mom(2)=sum(wVel*uVel*crh)
            mom(3)=sum(wVel*vVel*crh)
            mom(4)=0.5_KREAL*sum(wVel*((uVel*uVel+vVel*vVel)*crh+crb))
            scale=max(1.0_KREAL,maxval(abs(cw(:,i,j))))
            if(maxval(abs(mom))/scale>gate)call abort_native('S2 steady collision moment gate',444)
            do k=1,Nvel
                rh(k)=res_h_flux_store(k,i,j)*invdt+crh(k)
                rb(k)=res_b_flux_store(k,i,j)*invdt+crb(k)
            enddo
            if(IsS2Steady())then
                stream_h(:,i,j)=rh;stream_b(:,i,j)=rb
            endif
            do m=1,4
                sums(m)=sums(m)+ctr(i,j)%area*(rw(m)/i1_scale(m))**2
                linf(m)=max(linf(m),abs(rw(m))/i1_scale(m))
            enddo
            sums(5)=sums(5)+ctr(i,j)%area*sum(wVel*(rh/i1_scale(5))**2)/weight_norm
            sums(6)=sums(6)+ctr(i,j)%area*sum(wVel*(rb/i1_scale(6))**2)/weight_norm
            linf(5)=max(linf(5),maxval(abs(rh))/i1_scale(5))
            linf(6)=max(linf(6),maxval(abs(rb))/i1_scale(6))
            area_total=area_total+ctr(i,j)%area
        enddo;enddo
        norm=sqrt(sums/area_total)
        if(.not.all(ieee_is_finite(norm)).or..not.all(ieee_is_finite(linf))) &
            call abort_native('nonfinite S2 steady residual',445)
        if(IsS2Steady())then
            if(.not.all(ieee_is_finite(rate_w)).or..not.all(ieee_is_finite(stream_h)).or. &
               .not.all(ieee_is_finite(stream_b)))call abort_native('nonfinite S2 steady RHS',447)
            s2_steady_rhs=.true.
        endif
        i1_evaluations=i1_evaluations+1
        if(IsS2Steady().and.i1_evaluations==1)then
            diag_env=''
            call get_environment_variable('S2I1_R1_DIAG',diag_env,status=env_status)
            if(env_status==0.and.trim(diag_env)=='1')then
                open(newunit=fu,file='s2i1_rhs_arrays.bin',status='replace',access='stream',form='unformatted')
                write(fu)nx,ny,Nvel
                write(fu)ctr(1,1)%length(1),ctr(1,1)%length(2),dt_global,tau0
                write(fu)cw,ch,cb,mh,mb,kxh,kxb,kyh,kyb,qdump
                write(fu)fxw,fyw,fxh,fxb,fyh,fyb
                write(fu)res_macro_flux_store,res_h_flux_store,res_b_flux_store
                write(fu)rate_w,stream_h,stream_b,uVel,vVel,wVel
                close(fu)
            endif
        endif
        if(i1_evaluations==1)i1_initial_norm=norm
        i1_last_norm=norm;relative=0.0_KREAL
        do m=1,6
            if(i1_initial_norm(m)>100.0_KREAL*epsilon(1.0_KREAL)) &
                relative(m)=norm(m)/i1_initial_norm(m)
        enddo
        i1_below=all(norm<=1.0e-8_KREAL)
        if(IsS2Steady())then
            call SteadyGateObserve(iter,i1_below,i1_last_state_iter,i1_consecutive,counted,valid)
            if(.not.valid)call abort_native('S2 steady accepted-state number regressed',449)
        else
            counted=.false.
        endif
        write(i1_unit,'(I0,",",I0,",",ES24.16,12(",",ES24.16),",",ES24.16,",",I0,",",L1)') &
            iter,i1_evaluations,MPI_Wtime()-i1_start_wall,norm,linf,maxval(relative),i1_consecutive,counted
        flush(i1_unit)
    end subroutine S2NativeAssembleSteadyResidual

    subroutine S2NativeInvalidateTargets()
        if(.not.IsS2Steady())return
        targets_ready=.false.;evolution_ready=.false.
    end subroutine S2NativeInvalidateTargets

    subroutine S2NativeSaveCurrent(path)
        character(len=*),intent(in)::path
        if(.not.IsS2Steady())return
        call pack_state(cw,ch,cb)
        call write_state(path,cw,ch,cb,uVel,vVel,wVel,time_dis)
    end subroutine S2NativeSaveCurrent

    subroutine S2NativeApplyUpdate(accepted)
        logical,intent(out)::accepted
        real(KREAL),allocatable::w(:,:,:),h(:,:,:),b(:,:,:),wn(:,:,:),hs(:,:,:),bs(:,:,:)
        real(KREAL)::diag(12),dummy_flux,coll,t0c,t1c,total(4),drift,dx,dy,diag_step
        real(KREAL)::bottom(4),top(4),balance_next(4),delta(4),world_drift,xnet(4)
        integer::i,j,k,step,o1unit,side,mpi_e
        logical::ok
        accepted=.false.;if(.not.IsS2CN2())return
        if(.not.evolution_ready.or.residual_calls/=evolution_calls) &
            call abort_native('update lacks one native residual assembly',424)
        allocate(w(4,nx,ny),h(nx,ny,Nvel),b(nx,ny,Nvel),wn(4,nx,ny),hs(nx,ny,Nvel),bs(nx,ny,Nvel))
        call pack_state(w,h,b)
        if(max(maxval(abs(w-cw)),maxval(abs(h-ch)),maxval(abs(b-cb)))/=0.0_KREAL) &
            call abort_native('authoritative ctr changed between Evolution and update',425)
        if(maxval(abs(res_macro_store-res_macro_flux_store))/=0.0_KREAL) &
            call abort_native('unsupported source entered native residual',426)
        do j=1,ny;do i=1,nx
            wn(:,i,j)=w(:,i,j)+res_macro_flux_store(:,i,j)
            do k=1,Nvel
                hs(i,j,k)=h(i,j,k)+res_h_flux_store(k,i,j)
                bs(i,j,k)=b(i,j,k)+res_b_flux_store(k,i,j)
            enddo
        enddo;enddo
        dx=ctr(1,1)%length(1);dy=ctr(1,1)%length(2)
        if(cart_native)then;dx=mesh_dx;dy=mesh_dy;endif
        if(w1d_diagnostic)call w1d_dump_pre(update_calls+1,wn,hs,bs)
        call cpu_time(t0c)
        if(cart_native)then
            call CartesianAdvance(w,h,b,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,diag,dummy_flux,coll,ok,wn,hs,bs,face_error,limited_count,signed_count)
        else
        call advance_step(w,h,b,uVel,vVel,wVel,dx,dy,dt_global,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa, &
            diag,dummy_flux,coll,ok,wn,hs,bs,face_error,limited_count,signed_count)
        endif
        call cpu_time(t1c);collision_seconds=collision_seconds+t1c-t0c
        if(.not.ok)then
            if(S2Diagnostics.and..not.cart_native)call write_state('native_failed_stream.bin',wn,hs,bs,uVel,vVel,wVel,time_dis)
            if(S2Diagnostics.and..not.cart_native)call write_state('native_failed_proposal.bin',w,h,b,uVel,vVel,wVel,time_dis+dt_global)
            targets_ready=.false.
            call publish_state(w,h,b,.false.)
            call SaveRejectedAcceptedState('collision_or_target',.true.)
            deallocate(w,h,b,wn,hs,bs);return
        endif
        call cpu_time(t0c);total=0.0_KREAL
        do i=1,4;total(i)=sum(w(i,:,:));enddo
        if(cart_native)total=GlobalSum4(total)
        if(S2W1_wall_mode>=1)then
            bottom=0.0_KREAL;top=0.0_KREAL
            do i=1,nx
                if(ID_Y_minus==MPI_PROC_NULL)bottom=bottom+lower_w(:,i)/dy
                if(ID_Y_plus==MPI_PROC_NULL)top=top+fyw(:,i,ny)/dy
            enddo
            if(cart_native)then
                bottom=GlobalSum4(bottom);top=GlobalSum4(top)
            endif
            xnet=0.0_KREAL
            if(S2W1_wall_mode==2)then
                do j=1,ny
                    if(ID_X_minus==MPI_PROC_NULL)xnet=xnet-vface(1,j)%L_flux/ctr(1,j)%area
                    if(ID_X_plus==MPI_PROC_NULL)xnet=xnet+vface(nx+1,j)%R_flux/ctr(nx,j)%area
                enddo
                if(cart_native)xnet=GlobalSum4(xnet)
            endif
            balance_next=wall_balance+top-bottom+xnet
            drift=maxval(abs(total+balance_next-total0)/max(1.0_KREAL,abs(total0)))
        else
            drift=maxval(abs(total-total0)/max(1.0_KREAL,abs(total0)))
        endif
        if(force_global_reject)drift=max(drift,2.0_KREAL*gate)
        if(cart_native)then
            communication_clock=MPI_Wtime()
            call MPI_Allreduce(drift,world_drift,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_e)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            drift=world_drift
        endif
        if(drift>gate)then
            if(S2Diagnostics.and..not.cart_native)call write_state('native_failed_proposal.bin',w,h,b,uVel,vVel,wVel,time_dis+dt_global)
            targets_ready=.false.
            call publish_state(w,h,b,.false.)
            if(myID==0)then
            open(unit=93,file='native_rejection.csv',status='replace',action='write')
            write(93,'(A)')'reason,drift,accepted_calls,publish_calls,time,targets_ready,ctr_bitwise'
            call pack_state(wn,hs,bs)
            write(93,'(A,ES24.16,",",I0,",",I0,",",ES24.16,",",L1,",",L1)') &
                'global_conservation,',drift,update_calls,publish_calls,time_dis,targets_ready, &
                max(maxval(abs(wn-cw)),maxval(abs(hs-ch)),maxval(abs(bs-cb)))==0.0_KREAL
            close(93)
            endif
            call SaveRejectedAcceptedState('global_conservation',.true.)
            deallocate(w,h,b,wn,hs,bs);return
        endif
        if(cart_native.and.S2Diagnostics)call DumpLocalAccepted(w,h,b,time_dis+dt_global,'accepted')
        call M1DumpAudit(1)
        call publish_state(w,h,b,.true.);update_calls=update_calls+1;accepted=.true.
        if(IsS2CN2())call RecordAcceptedMacroUpdateRate(cw,dt_global)
        step=update_calls
        if(S2PhysicalGas==1.and.steps_open)then
            write(gas_unit,'(I0,6(",",ES24.16),3(",",I0))')step,gas_tau_min,gas_tau_max,gas_tau_calls
            flush(gas_unit)
        endif
        if(w1d_diagnostic)call w1d_dump_post(step,w,h,b,diag)
        if(S2W1_wall_mode>=1)then
            wall_balance=balance_next;delta=total-total0
            cumulative_bottom=cumulative_bottom+bottom;cumulative_top=cumulative_top+top
            call PublishWallDiagnostics()
            if(steps_open)then
            write(wall_unit,'(I0,18(",",ES24.16))')step,time_dis+dt_global,delta,bottom,top,wall_balance,drift
            flush(wall_unit)
            do i=1,nx
                do side=1,2
                    if(side==1)then
                        bottom=lower_w(:,i)
                    else
                        bottom=fyw(:,i,ny)
                    endif
                    write(wall_segment_unit,'(3(I0,","),6(ES24.16,:,","))')step,i,side,bottom, &
                        bottom(2)/dt_global,(bottom(4)-S2W1_wall_ux(side)*bottom(2))/dt_global
                enddo
            enddo
            flush(wall_segment_unit)
            endif
            if(S2W1_time_centered==1)then
                write(w1t_cost_unit,'(I0,5(",",ES24.16),3(",",I0))')step, &
                    MPI_Wtime()-o1_wall_start,w1t_official_face_s,w1t_predictor_step_s, &
                    w1t_predictor_wall_s,coll,1,1,1
                flush(w1t_cost_unit)
            endif
            if(S2W1D_face_time_consistent==1)then
                write(w1d_cost_unit,'(I0,",",ES24.16,3(",",I0))')step,w1d_instant_s,1, &
                    w1d_instant_limited,w1d_instant_signed
                flush(w1d_cost_unit)
            endif
        endif
        if(w2r_metrics)then
            write(w2r_unit,'(I0,",",I0,15(",",ES24.16),6(",",I0))')step,w2r_active, &
                MPI_Wtime()-o1_wall_start,w2_off_reconstruction_s,w2r_off_closure_s,w2_off_wall_s, &
                w2r_predict_s,w2r_instant_s,w2r_apply_s,coll,w2r_immutable_error,w2r_diag, &
                1,0,1,w2r_active,w2r_active,1
            flush(w2r_unit)
        endif
        if(w2_cost_enabled)then
            write(w2_cost_unit,'(I0,2(",",ES24.16),",",I0)')step, &
                w2_off_reconstruction_s,w2_off_wall_s,w2_off_wall_calls
            flush(w2_cost_unit)
        endif
        call cpu_time(t1c);diag_step=t1c-t0c;diagnostic_seconds=diagnostic_seconds+diag_step
        if(steps_open)then
        write(steps_unit,'(I0,19(",",ES16.8))')step,time_dis+dt_global,diag,drift,layout_seconds, &
            face_seconds,map_seconds,coll,diag_step
        flush(steps_unit)
        if(step==1.and..not.o1_bench) &
            call write_state('native_state_prefix.bin',w,h,b,uVel,vVel,wVel,time_dis+dt_global)
        endif
        deallocate(w,h,b,wn,hs,bs)
        if(o1_bench)then
            open(newunit=o1unit,file='o1_native_step_wall.csv',status='replace',action='write')
            write(o1unit,'(A)')'full_native_step_wall_s'
            write(o1unit,'(ES24.16)')MPI_Wtime()-o1_wall_start
            close(o1unit)
        endif
    end subroutine S2NativeApplyUpdate

    subroutine S2NativeRecordLowOrderStep()
        real(KREAL)::total(4),drift,moment(4),err
        integer::i,j,k
        if(trim(adjustl(PeriodicComparisonProfile))/='BGK_ConstTau')return
        total=0.0_KREAL;err=0.0_KREAL
        do j=1,ny;do i=1,nx
            total=total+ctr(i,j)%conVars
            moment=0.0_KREAL
            do k=1,Nvel
                moment=moment+wVel(k)*[ctr(i,j)%h_vec(k),uVel(k)*ctr(i,j)%h_vec(k), &
                    vVel(k)*ctr(i,j)%h_vec(k), &
                    0.5_KREAL*((uVel(k)**2+vVel(k)**2)*ctr(i,j)%h_vec(k)+ctr(i,j)%b_vec(k))]
            enddo
            err=max(err,maxval(abs(moment-ctr(i,j)%conVars)))
        enddo;enddo
        drift=maxval(abs(total-total0)/max(1.0_KREAL,abs(total0)))
        update_calls=update_calls+1;publish_calls=publish_calls+1
        if(steps_open)then
        write(steps_unit,'(I0,19(",",ES16.8))')update_calls,time_dis+dt_global, &
                0.0_KREAL,0.0_KREAL,0.0_KREAL,err,0.0_KREAL,0.0_KREAL,dt_global/tau0, &
                dt_global*(maxval(abs(uVel))/ctr(1,1)%length(1)+maxval(abs(vVel))/ctr(1,1)%length(2)), &
                minval(cell_h_store(:,1:nx,1:ny)),maxval(cell_h_store(:,1:nx,1:ny)), &
                minval(cell_b_store(:,1:nx,1:ny)),maxval(cell_b_store(:,1:nx,1:ny)),drift, &
                layout_seconds,face_seconds,map_seconds,collision_seconds,diagnostic_seconds
            flush(steps_unit)
        endif
    end subroutine S2NativeRecordLowOrderStep

    subroutine S2NativeReadMetadata(validate_only)
        logical,intent(in)::validate_only
        call ReadNativeMetadata(validate_only,total0,wall_balance,cumulative_bottom,cumulative_top)
        if(.not.validate_only)then
            targets_ready=.false.;evolution_ready=.false.
            if(S2Diagnostics)then
                call pack_state(cw,ch,cb)
                call DumpLocalAccepted(cw,ch,cb,time_dis,'initial')
            endif
            ! InitFlow performs native physical and MPI halo reconstruction afterwards.
        endif
    end subroutine

    subroutine S2NativeCommitRestart()
        call CommitNativeRestart(total0,wall_balance,cumulative_bottom,cumulative_top)
    end subroutine

    function S2NativeOutputHeat(i,j) result(q)
        integer,intent(in)::i,j
        real(KREAL)::q(2),ah(Nvel),ab(Nvel)
        logical::ok
        call analytic_maxwellian_2d(ctr(i,j)%conVars,ck,uVel,vVel,ah,ab,ok)
        if(.not.ok)call abort_native('output analytic target failed',428)
        call heat(ctr(i,j)%conVars,ctr(i,j)%h_vec-ah,ctr(i,j)%b_vec-ab,uVel,vVel,wVel,q,ok)
        if(.not.ok)call abort_native('output heat evaluation failed',428)
    end function

    subroutine PublishWallDiagnostics()
        integer::side,i,j
        real(KREAL)::hh(Nvel),bb(Nvel),ww(4),prim(4),rate(4),measure
        logical::ok
        do side=1,2
            if(side==1.and.ID_Y_minus/=MPI_PROC_NULL)cycle
            if(side==2.and.ID_Y_plus/=MPI_PROC_NULL)cycle
            j=1;if(side==2)j=ny+1
            do i=1,nx
                measure=dt_global*hface(i,j)%length
                ! Stored face flux includes dt and face length; v is nonzero in supported tables.
                hh=hface(i,j)%L_flux_h_vec/(measure*vVel)
                bb=hface(i,j)%L_flux_b_vec/(measure*vVel)
                call moment_2d(hh,bb,uVel,vVel,wVel,ww,ok)
                if(.not.ok)call abort_native('wall output moment failed',459)
                prim=GetPrimary(ww)
                hface(i,j)%FlowVars=[prim(1:3),1.0_KREAL/prim(4),0.5_KREAL*prim(1)/prim(4)]
                rate=hface(i,j)%L_flux/measure
                ! Native frame n=+y, t=-x. Heat removes moving-wall mechanical work.
                hface(i,j)%Aerodynamics=[rate(3),-rate(2),rate(4)-S2W1_wall_ux(side)*rate(2)]
            enddo
        enddo
        if(S2W1_wall_mode==2)then
            do side=1,2
                if(side==1.and.ID_X_minus/=MPI_PROC_NULL)cycle
                if(side==2.and.ID_X_plus/=MPI_PROC_NULL)cycle
                i=1;if(side==2)i=nx+1
                do j=1,ny
                    measure=dt_global*vface(i,j)%length
                    hh=vface(i,j)%L_flux_h_vec/(measure*uVel)
                    bb=vface(i,j)%L_flux_b_vec/(measure*uVel)
                    call moment_2d(hh,bb,uVel,vVel,wVel,ww,ok)
                    if(.not.ok)call abort_native('vertical wall output moment failed',459)
                    prim=GetPrimary(ww)
                    vface(i,j)%FlowVars=[prim(1:3),1.0_KREAL/prim(4),0.5_KREAL*prim(1)/prim(4)]
                    rate=vface(i,j)%L_flux/measure
                    vface(i,j)%Aerodynamics=[rate(2),rate(3),rate(4)]
                enddo
            enddo
        endif
    end subroutine

    subroutine SaveRejectedAcceptedState(reason,collective)
        character(len=*),intent(in)::reason
        logical,intent(in),optional::collective
        integer::fu,e
        call pack_state(cw,ch,cb)
        if(cart_native)then
            call DumpLocalAccepted(cw,ch,cb,time_dis,'last_accepted')
        else
            call write_state('last_accepted.bin',cw,ch,cb,uVel,vVel,wVel,time_dis)
        endif
        if(myID==0)then
        open(newunit=fu,file='last_accepted_ledger.txt',status='replace')
        write(fu,'(A)')reason
        write(fu,*)iter,time_dis
        write(fu,*)total0,wall_balance,cumulative_bottom,cumulative_top
        write(fu,*)'targets_ready=',targets_ready,' evolution_ready=',evolution_ready,' update/publish=',update_calls,publish_calls
        close(fu)
        endif
        if(present(collective))then
            if(collective)call MPI_Barrier(MPI_COMM_WORLD,e)
        endif
    end subroutine

    subroutine S2NativeRunSummary()
        integer::fu,i,j
        real(KREAL)::area,total(4)
        area=ctr(1,1)%area;total=0.0_KREAL
        do j=1,ny;do i=1,nx
            total=total+ctr(i,j)%conVars
        enddo;enddo
        if(cart_native)total=GlobalSum4(total)
        if(myID/=0)return
        open(newunit=fu,file='Simulation_Log.txt',status='old',position='append')
        write(fu,'(A,I0,A,ES24.16)')'Accepted final step=',iter,' physical time=',time_dis
        write(fu,'(A,4ES24.16)')'Accepted macro residual=',res
        if(S2W1_wall_mode>=1)then
            write(fu,'(A,4ES24.16)')'Integrated bottom +y flux (mass,mx,my,E)=',area*cumulative_bottom
            write(fu,'(A,4ES24.16)')'Integrated top +y flux (mass,mx,my,E)=',area*cumulative_top
            write(fu,'(A,4ES24.16)')'Integrated outward net wall flux=',area*wall_balance
            write(fu,'(A,4ES24.16)')'Volume change plus net wall flux=',area*(total-total0+wall_balance)
        endif
        close(fu)
    end subroutine

    subroutine S2NativeWriteOutput()
        real(KREAL)::rho,ux,vy,temp,p,q(2),x,y,t0c,t1c
        integer::i,j,unit
        logical::ok
        if(.not.package_mode())return
        if(cart_native)then
            call pack_state(cw,ch,cb)
            call DumpLocalAccepted(cw,ch,cb,time_dis,'final')
            return
        endif
        call cpu_time(t0c);call pack_state(cw,ch,cb)
        if(IsS2Steady())then
            call native_build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,rho,ok)
            if(.not.ok)call abort_native('S2 steady final target construction failed',448)
            targets_ready=.true.
        elseif(.not.IsS2CN2())then
            call native_build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,rho,ok)
            if(.not.ok)call abort_native('low-order output target construction failed',428)
        endif
        call write_state('native_state_final.bin',cw,ch,cb,uVel,vVel,wVel,time_dis)
        open(newunit=unit,file='native_field.csv',status='replace',action='write')
        write(unit,'(A)')'x,y,rho,U,V,p,T,qx,qy,time'
        do j=1,ny;do i=1,nx
            rho=cw(1,i,j);ux=cw(2,i,j)/rho;vy=cw(3,i,j)/rho
            p=(2.0_KREAL/real(ck+2,KREAL))*(cw(4,i,j)-0.5_KREAL*rho*(ux*ux+vy*vy));temp=p/rho
            call heat(cw(:,i,j),ch(i,j,:)-eah(i,j,:),cb(i,j,:)-eab(i,j,:),uVel,vVel,wVel,q,ok)
            if(.not.ok)call abort_native('output heat evaluation failed',428)
            x=ctr(i,j)%x;y=ctr(i,j)%y
            write(unit,'(10(ES24.16,:,","))')x,y,rho,ux,vy,p,temp,q(1),q(2),time_dis
        enddo;enddo
        close(unit);output_calls=output_calls+1
        if(steps_open.and.iter>=MAX_ITER)then
            close(steps_unit);steps_unit=-1;steps_open=.false.
        endif
        if(S2W1_wall_mode>=1.and.iter>=MAX_ITER.and.(S2Diagnostics.or..not.IsS2CN2()))then
            close(wall_unit);close(wall_segment_unit)
            if(S2W1_time_centered==1)close(w1t_cost_unit)
            if(S2W1D_face_time_consistent==1)close(w1d_cost_unit)
        endif
        open(newunit=unit,file='native_call_counts.csv',status='replace',action='write')
        write(unit,'(A)')'init,evolution,residual,update,publish,output,atomic_failure_passed'
        write(unit,'(6(I0,","),L1)')init_calls,evolution_calls,residual_calls,update_calls,publish_calls,output_calls, &
            atomic_test_passed
        close(unit)
        if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau')then
            open(newunit=unit,file='native_tau.csv',status='replace',action='write')
            write(unit,'(A)')'location,min_tau,max_tau,calls'
            write(unit,'(A,",",ES24.16,",",ES24.16,",",I0)') &
                'face',ComparisonFaceTauMin,ComparisonFaceTauMax,ComparisonFaceTauCalls
            write(unit,'(A,",",ES24.16,",",ES24.16,",",I0)') &
                'cell',ComparisonCellTauMin,ComparisonCellTauMax,ComparisonCellTauCalls
            close(unit)
        endif
        call cpu_time(t1c);io_seconds=io_seconds+t1c-t0c
        open(newunit=unit,file='native_timing.csv',status='replace',action='write')
        write(unit,'(A)')'layout_s,face_s,face_map_s,collision_s,diagnostic_s,io_s'
        write(unit,'(6(ES24.16,:,","))')layout_seconds,face_seconds,map_seconds,collision_seconds, &
            diagnostic_seconds,io_seconds
        close(unit)
    end subroutine S2NativeWriteOutput
end module mod_S2NativeAdapter

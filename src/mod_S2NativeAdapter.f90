module mod_S2NativeAdapter
    use Precision, only: KREAL,KINT
    use flow_data
    use flow_para, only: ck,Pr,BC_Type
    use control_para, only: dt_global,TimeStep_method,TimeIntegration_method,PeriodicComparisonProfile, &
        S2D1_initial_mode,MAX_ITER,iter,ComparisonFaceTauMin,ComparisonFaceTauMax, &
        ComparisonCellTauMin,ComparisonCellTauMax,ComparisonFaceTauCalls,ComparisonCellTauCalls
    use Flow_DIS, only: time_dis,Iflag_DIS,Iflag_Fringe
    use mpi_para, only: myID,NP_X,NP_Y,NP_XY
    use const_para, only: BCconst_periodic,Time_global_fixed
    use mod_S2D1, only: IsS2CN2,IsS2Steady,compute_flux,build_targets,advance_step,init_velocity,write_state,heat,make_mke
    use mod_S2SteadyTermination, only: SteadyGateObserve,SteadyNeedsEvaluation
    use mod_ImplicitState, only: rate_w,stream_h,stream_b,s2_steady_rhs
    use mpi, only: MPI_COMM_WORLD,MPI_Abort,MPI_Wtime
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: S2NativeInitializeState,S2NativeEvolution,S2NativeApplyUpdate,S2NativeMarkResidual
    public :: S2NativeWriteOutput,S2NativeIsReady,S2NativeRecordLowOrderStep
    public :: S2NativeAssembleSteadyResidual,S2NativeInvalidateTargets,S2NativeSaveCurrent
    public :: S2NativeSteadyConverged,S2NativeSteadyNeedsEvaluation,S2NativeDiagnosticEnabled
    public :: S2NativeSteadyMaxNorm

    real(KREAL),parameter::pi_d=acos(-1.0_KREAL)
    real(KREAL),parameter::rho0=1.0_KREAL,u0=0.6_KREAL,v0=0.3_KREAL,t0=1.0_KREAL
    real(KREAL),parameter::amp=0.01_KREAL,gate=2.0e-12_KREAL,tau0=0.1_KREAL
    real(KREAL),allocatable,save::cw(:,:,:),ch(:,:,:),cb(:,:,:)
    real(KREAL),allocatable,save::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:)
    real(KREAL),allocatable,save::eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
    real(KREAL),allocatable,save::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:)
    real(KREAL),save::face_error=0.0_KREAL,total0(4)=0.0_KREAL
    real(KREAL),save::layout_seconds=0.0_KREAL,face_seconds=0.0_KREAL,map_seconds=0.0_KREAL
    real(KREAL),save::collision_seconds=0.0_KREAL,diagnostic_seconds=0.0_KREAL,io_seconds=0.0_KREAL
    integer,save::limited_count=0,signed_count=0,init_calls=0,evolution_calls=0,residual_calls=0
    integer,save::update_calls=0,publish_calls=0,output_calls=0,steps_unit=-1
    logical,save::initialized=.false.,evolution_ready=.false.,atomic_test_passed=.false.
    logical,save::targets_ready=.false.,steps_open=.false.,force_global_reject=.false.
    logical,save::o1_bench=.false.
    real(KREAL),save::o1_wall_start=0.0_KREAL
    real(KREAL),save::i1_scale(6)=1.0_KREAL,i1_initial_norm(6)=0.0_KREAL
    real(KREAL),save::i1_last_norm(6)=0.0_KREAL,i1_start_wall=0.0_KREAL
    integer,save::i1_evaluations=0,i1_consecutive=0,i1_unit=-1
    integer,save::i1_last_state_iter=-1
    logical,save::i1_diag_enabled=.false.,i1_below=.false.

contains

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
        if(myID==0)write(*,'(A)')'S2CN2 FAILURE: '//trim(message)
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

    subroutine S2NativeInitializeState()
        real(KREAL)::dx,dy,sincx,sincy,x,y,rho,err,cfl,tstart,tend,tol
        real(KREAL)::tu(Nvel),tv(Nvel),tw(Nvel),raw_node(16),raw_weight(16),momerr
        integer::i,j,m,env_status,ru,rx,ry,rv
        character(len=8)::reject_env
        character(len=1024)::restart_path
        real(KREAL)::restart_time,ruv(Nvel),rvv(Nvel),rwt(Nvel)
        logical::ok
        if(.not.package_mode())return
        call cpu_time(tstart);init_calls=init_calls+1
        if(NP_X/=1.or.NP_Y/=1.or.NP_XY/=1)call abort_native('only one MPI rank is supported',401)
        if(nx/=ny.or.all(nx/=[8,16,32,64]).or.Nvel/=256.or.ck/=1) &
            call abort_native('requires N8/N16/N32/N64, GH16xGH16, CK=1',402)
        if(any(BC_Type/=BCconst_periodic))call abort_native('all four boundaries must be periodic',403)
        if(TimeStep_method/=Time_global_fixed.or..not.ieee_is_finite(dt_global).or.dt_global<=0.0_KREAL) &
            call abort_native('requires positive fixed dt',404)
        if(Iflag_DIS/=0.or.Iflag_Fringe/=0)call abort_native('DIS and fringe must be disabled',405)
        dx=ctr(1,1)%length(1);dy=ctr(1,1)%length(2)
        tol=1.0e-14_KREAL
        if(abs(dx-2.0_KREAL*pi_d/real(nx,KREAL))>tol.or.abs(dy-dx)>tol) &
            call abort_native('requires uniform [0,2pi]^2 Cartesian mesh',406)
        do j=1,ny+1;do i=1,nx+1
            if(abs(gridx(i,j)-real(i-1,KREAL)*dx)>tol.or. &
               abs(gridy(i,j)-real(j-1,KREAL)*dy)>tol) &
                call abort_native('Cartesian coordinates or periodic seam mismatch',406)
        enddo;enddo
        do j=1,ny;do i=1,nx
            if(abs(ctr(i,j)%length(1)-dx)>tol.or.abs(ctr(i,j)%length(2)-dy)>tol.or. &
                abs(ctr(i,j)%area-dx*dy)>tol.or.abs(ctr(i,j)%x-(real(i,KREAL)-0.5_KREAL)*dx)>tol.or. &
                abs(ctr(i,j)%y-(real(j,KREAL)-0.5_KREAL)*dy)>tol) &
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
        call init_velocity(tu,tv,tw,raw_node,raw_weight)
        err=max(maxval(abs(tu-uVel)),maxval(abs(tv-vVel)),maxval(abs(tw-wVel)))
        if(err/=0.0_KREAL)call abort_native('native velocity table differs from frozen R1 table',408)
        cfl=dt_global*(maxval(abs(uVel))/dx+maxval(abs(vVel))/dy)
        if(cfl>0.4_KREAL)call abort_native('transport CFL exceeds 0.4',409)
        allocate(cw(4,nx,ny),ch(nx,ny,Nvel),cb(nx,ny,Nvel),mh(nx,ny,Nvel),mb(nx,ny,Nvel), &
            kxh(nx,ny,Nvel),kxb(nx,ny,Nvel),kyh(nx,ny,Nvel),kyb(nx,ny,Nvel), &
            eah(nx,ny,Nvel),eab(nx,ny,Nvel),kappa(2,2,nx,ny),fxw(4,nx,ny),fyw(4,nx,ny), &
            fxh(nx,ny,Nvel),fxb(nx,ny,Nvel),fyh(nx,ny,Nvel),fyb(nx,ny,Nvel))
        sincx=sin(0.5_KREAL*dx)/(0.5_KREAL*dx);sincy=sin(0.5_KREAL*dy)/(0.5_KREAL*dy)
        do j=1,ny;y=ctr(1,j)%y
            do i=1,nx;x=ctr(i,j)%x
                select case(S2D1_initial_mode)
                case(0);rho=rho0+amp*sincx*sincy*sin(x+y)
                case(3);rho=rho0
                case default;call abort_native('native package supports initial modes 0 and 3 only',410)
                end select
                cw(:,i,j)=[rho,rho*u0,rho*v0,rho*(0.5_KREAL*(u0*u0+v0*v0)+1.5_KREAL*t0)]
                call make_mke(cw(:,i,j),uVel,vVel,wVel,mh(i,j,:),mb(i,j,:),kxh(i,j,:),kxb(i,j,:), &
                    kyh(i,j,:),kyb(i,j,:),eah(i,j,:),eab(i,j,:),kappa(:,:,i,j),ok,momerr)
                if(.not.ok)call abort_native('discrete conservative initialization failed',413)
                ch(i,j,:)=mh(i,j,:);cb(i,j,:)=mb(i,j,:)
            enddo
        enddo
        call build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,err,ok)
        if(.not.ok)call abort_native('initial volume target construction failed',414)
        targets_ready=.true.
        call publish_state(cw,ch,cb,.true.);publish_calls=0
        do j=1,ny;do i=1,nx
            ctr(i,j)%dt_local=dt_global;ctr(i,j)%conVars0=ctr(i,j)%conVars
            ctr(i,j)%h0_vec=ctr(i,j)%h_vec;ctr(i,j)%b0_vec=ctr(i,j)%b_vec
            ctr(i,j)%res_conVars0=0.0_KREAL;ctr(i,j)%res_h0_vec=0.0_KREAL;ctr(i,j)%res_b0_vec=0.0_KREAL
        enddo;enddo
        total0=0.0_KREAL
        do i=1,4;total0(i)=sum(cw(i,:,:));enddo
        call check_roundtrip_and_atomicity()
        restart_path=''
        call get_environment_variable('S2DO1_RESTART_PATH',restart_path,status=env_status)
        if(env_status==0.and.len_trim(restart_path)>0)then
            open(newunit=ru,file=trim(restart_path),status='old',access='stream',form='unformatted',action='read')
            read(ru)rx,ry,rv,restart_time
            if(rx/=nx.or.ry/=ny.or.rv/=Nvel)call abort_native('O1 restart dimensions disagree',431)
            read(ru)cw;read(ru)ch;read(ru)cb;read(ru)ruv;read(ru)rvv;read(ru)rwt
            close(ru)
            if(max(maxval(abs(ruv-uVel)),maxval(abs(rvv-vVel)),maxval(abs(rwt-wVel)))/=0.0_KREAL) &
                call abort_native('O1 restart velocity table disagrees',432)
            call build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,err,ok)
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
        call write_state('native_state_initial.bin',cw,ch,cb,uVel,vVel,wVel,0.0_KREAL)
        open(newunit=steps_unit,file='native_steps.csv',status='replace',action='write')
        steps_open=.true.
        write(steps_unit,'(A)')'step,time,face_moment,target_moment,cn_residual,macro_micro,limited,signed,'// &
            'dt_tau,cfl,h_min,h_max,b_min,b_max,global_drift,layout_s,face_s,map_s,collision_s,diagnostic_s'
        reject_env=''
        call get_environment_variable('S2D3_TEST_GLOBAL_REJECT',reject_env,status=env_status)
        force_global_reject=IsS2CN2().and.env_status==0.and.trim(reject_env)=='1'
        initialized=.true.;evolution_ready=.false.;call cpu_time(tend);layout_seconds=layout_seconds+tend-tstart
    end subroutine S2NativeInitializeState

    subroutine S2NativeEvolution()
        real(KREAL)::dx,dy,t0c,t1c
        integer::i,j
        logical::ok
        if(.not.(IsS2CN2().or.IsS2Steady()))return
        if(.not.initialized.or.evolution_ready)call abort_native('invalid Evolution adapter state',421)
        o1_wall_start=MPI_Wtime()
        call cpu_time(t0c);call pack_state(cw,ch,cb);call cpu_time(t1c);layout_seconds=layout_seconds+t1c-t0c
        if(.not.targets_ready)then
            call build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,face_error,ok)
            if(.not.ok)call abort_native('accepted-state target rebuild failed',420)
            targets_ready=.true.
        endif
        dx=ctr(1,1)%length(1);dy=ctr(1,1)%length(2)
        call cpu_time(t0c)
        call compute_flux(cw,ch,cb,uVel,vVel,wVel,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb, &
            face_error,limited_count,signed_count,ok)
        call cpu_time(t1c);face_seconds=face_seconds+t1c-t0c
        if(.not.ok)then
            if(IsS2Steady())call write_state('s2i1_face_rejected_last_accepted.bin', &
                cw,ch,cb,uVel,vVel,wVel,0.0_KREAL)
            call abort_native('face kernel rejected state',422)
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
            hface(i,1)%L_flux=fyw(:,i,ny)*hface(i,1)%length;hface(i,1)%R_flux=hface(i,1)%L_flux
            hface(i,ny+1)%L_flux=hface(i,1)%L_flux;hface(i,ny+1)%R_flux=hface(i,1)%L_flux
            hface(i,1)%L_flux_h_vec=fyh(i,ny,:)*hface(i,1)%length
            hface(i,1)%R_flux_h_vec=hface(i,1)%L_flux_h_vec
            hface(i,ny+1)%L_flux_h_vec=hface(i,1)%L_flux_h_vec;hface(i,ny+1)%R_flux_h_vec=hface(i,1)%L_flux_h_vec
            hface(i,1)%L_flux_b_vec=fyb(i,ny,:)*hface(i,1)%length
            hface(i,1)%R_flux_b_vec=hface(i,1)%L_flux_b_vec
            hface(i,ny+1)%L_flux_b_vec=hface(i,1)%L_flux_b_vec;hface(i,ny+1)%R_flux_b_vec=hface(i,1)%L_flux_b_vec
        enddo
        call cpu_time(t1c);map_seconds=map_seconds+t1c-t0c
        evolution_calls=evolution_calls+1;evolution_ready=.true.
    end subroutine S2NativeEvolution

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
        integer::i,j,k,step,o1unit
        logical::ok
        accepted=.false.;if(.not.IsS2CN2())return
        if(.not.evolution_ready.or.residual_calls/=evolution_calls)call abort_native('update lacks one native residual assembly',424)
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
        call cpu_time(t0c)
        call advance_step(w,h,b,uVel,vVel,wVel,dx,dy,dt_global,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa, &
            diag,dummy_flux,coll,ok,wn,hs,bs,face_error,limited_count,signed_count)
        call cpu_time(t1c);collision_seconds=collision_seconds+t1c-t0c
        if(.not.ok)then
            call write_state('native_failed_stream.bin',wn,hs,bs,uVel,vVel,wVel,time_dis)
            call write_state('native_failed_proposal.bin',w,h,b,uVel,vVel,wVel,time_dis+dt_global)
            targets_ready=.false.
            call publish_state(w,h,b,.false.);deallocate(w,h,b,wn,hs,bs);return
        endif
        call cpu_time(t0c);total=0.0_KREAL
        do i=1,4;total(i)=sum(w(i,:,:));enddo
        drift=maxval(abs(total-total0)/max(1.0_KREAL,abs(total0)))
        if(force_global_reject)drift=max(drift,2.0_KREAL*gate)
        if(drift>gate)then
            call write_state('native_failed_proposal.bin',w,h,b,uVel,vVel,wVel,time_dis+dt_global)
            targets_ready=.false.
            call publish_state(w,h,b,.false.)
            open(unit=93,file='native_rejection.csv',status='replace',action='write')
            write(93,'(A)')'reason,drift,accepted_calls,publish_calls,time,targets_ready,ctr_bitwise'
            call pack_state(wn,hs,bs)
            write(93,'(A,ES24.16,",",I0,",",I0,",",ES24.16,",",L1,",",L1)') &
                'global_conservation,',drift,update_calls,publish_calls,time_dis,targets_ready, &
                max(maxval(abs(wn-cw)),maxval(abs(hs-ch)),maxval(abs(bs-cb)))==0.0_KREAL
            close(93)
            deallocate(w,h,b,wn,hs,bs);return
        endif
        call publish_state(w,h,b,.true.);update_calls=update_calls+1;accepted=.true.
        step=update_calls
        call cpu_time(t1c);diag_step=t1c-t0c;diagnostic_seconds=diagnostic_seconds+diag_step
        write(steps_unit,'(I0,19(",",ES16.8))')step,time_dis+dt_global,diag,drift,layout_seconds, &
            face_seconds,map_seconds,coll,diag_step
        flush(steps_unit)
        if(step==1.and..not.o1_bench) &
            call write_state('native_state_prefix.bin',w,h,b,uVel,vVel,wVel,time_dis+dt_global)
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

    subroutine S2NativeWriteOutput()
        real(KREAL)::rho,ux,vy,temp,p,q(2),x,y,t0c,t1c
        integer::i,j,unit
        logical::ok
        if(.not.package_mode())return
        call cpu_time(t0c);call pack_state(cw,ch,cb)
        if(IsS2Steady())then
            call build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,rho,ok)
            if(.not.ok)call abort_native('S2 steady final target construction failed',448)
            targets_ready=.true.
        elseif(.not.IsS2CN2())then
            call build_targets(cw,uVel,vVel,wVel,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,rho,ok)
            if(.not.ok)call abort_native('low-order output target construction failed',428)
        endif
        call write_state('native_state_final.bin',cw,ch,cb,uVel,vVel,wVel,time_dis)
        open(newunit=unit,file='native_field.csv',status='replace',action='write')
        write(unit,'(A)')'x,y,rho,U,V,p,T,qx,qy,time'
        do j=1,ny;do i=1,nx
            rho=cw(1,i,j);ux=cw(2,i,j)/rho;vy=cw(3,i,j)/rho
            p=(2.0_KREAL/3.0_KREAL)*(cw(4,i,j)-0.5_KREAL*rho*(ux*ux+vy*vy));temp=p/rho
            call heat(cw(:,i,j),ch(i,j,:)-eah(i,j,:),cb(i,j,:)-eab(i,j,:),uVel,vVel,wVel,q,ok)
            if(.not.ok)call abort_native('output heat evaluation failed',428)
            x=ctr(i,j)%x;y=ctr(i,j)%y
            write(unit,'(10(ES24.16,:,","))')x,y,rho,ux,vy,p,temp,q(1),q(2),time_dis
        enddo;enddo
        close(unit);output_calls=output_calls+1
        if(steps_open.and.iter>=MAX_ITER)then
            close(steps_unit);steps_unit=-1;steps_open=.false.
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

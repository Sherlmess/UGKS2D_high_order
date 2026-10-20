module ImportPara
    use mod_S2MPI1Layout, only: slab_enabled
    use mod_VelocityTables, only: VelocityMode,VelocityCount,VelocityMax
    use flow_para
    use control_para
    use tools
    use Flow_DIS
    use mpi
    use, intrinsic :: ieee_arithmetic
    implicit none

    type MethodName
        character(len=30) :: name
    end type MethodName

contains

    !--------------------------------------------------
    !>Initialize parameters
    !--------------------------------------------------
    subroutine InitPara()
        implicit none
        logical :: ext1

        block
            character(len=8)::mpi1_flag
            call get_environment_variable("S2_MPI1",mpi1_flag)
            slab_enabled=trim(mpi1_flag)=="1"
        end block
        call Set_default_parameters

        if (myid == 0) then
            inquire(file="Config.ini", exist=ext1)
            if (ext1) then
                call ReadParameters
            else
                print *, "Can not find 'Config.ini', stop !"
                call MPI_Abort(MPI_COMM_WORLD,472,ierr)
            endif
        endif

        call bcast_para
        call ValidateSupportedOptions
        call set_const_para
        if (myid == 0) call write_log

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine InitPara

    !--------------------------------------------------
    !>Initialize parameters (set default value)
    !--------------------------------------------------
    subroutine Set_default_parameters()
        Ma = 0.0
        Re = 0.0
        AOA = 0.d0
        Tinf = 0.d0
        Tw = 0.d0

        ck = 1
        Pr = 2.0_KREAL / 3.0_KREAL
        flag_vislaw = PowerLaw
        omega = 0.5
        Tcm = 0.4044

        CFL = 0.1
        MAX_ITER = 10000
        EPS = 1.0E-7
        Nstep_show = 10
        Nstep_write = 10000
        Nstep_restart=100
        Nstep_audit=0
        Nstep_monitor=1
        Quadrature_method = Quad_NewtonCotes4
        TimeStep_method = Time_global_varying
        Iflag_InitFlow = Init_Freestream
        Iflag_WriteSurf = 0
        Iflag_WriteRestart = 0
        TimeIntegration_method = 'Explicit'
        PeriodicComparisonProfile = 'Off'
        ComparisonFaceTauMin=huge(1.0_KREAL)
        ComparisonFaceTauMax=0.0_KREAL
        ComparisonCellTauMin=huge(1.0_KREAL)
        ComparisonCellTauMax=0.0_KREAL
        ComparisonFaceTauCalls=0
        ComparisonCellTauCalls=0
        S2D1_initial_mode = 0
        S2PhysicalGas = 0
        S2_mu_ref = 0.1_KREAL
        S2W1_wall_mode = 0
        S2W2R_shared_predictor = 0
        S2W1_time_centered = 0
        S2W1D_face_time_consistent = 0
        S2W1_wall_ux = 0.0_KREAL
        S2W1_wall_theta = -huge(1.0_KREAL) ! omitted per-side override
        Implicit_nCFL = 100.0_KREAL
        Implicit_MacroSweeps = 2
        Implicit_MicroSweeps = 1
        Implicit_ResidualTol = 1.0e-8_KREAL
        Implicit_MomentTol = 1.0e-12_KREAL
        Implicit_MinAlpha = 1.0e-6_KREAL
        Implicit_MaxRetries = 8
        Implicit_MaxIter = 100000
        Unsteady_Theta = 0.5_KREAL
        Unsteady_InnerTol = 1.0e-8_KREAL
        Unsteady_MaxInner = 50
        simTime = 0.0
        iter = 0
        time_dis=0.0_KREAL
        dt_global = 1.0E-6
        res(4) = 0.0

        Nu = 28
        Nv = 28
        U_MIN = -4.0
        U_MAX = +4.0
        V_MIN = -4.0
        V_MAX = +4.0

        BC_Type = 3
        BG_geometry = BGconst_general
        FluxWall_type = FluxWall_Kinetic
        ACC = 1.d0

        call date_and_time(DATE=date, TIME=time)
        fileName = '_' // date

        Iflag_DIS = 0
        DISBC_Type = 0
        Iflag_Fringe = 0
    end subroutine Set_default_parameters

    !--------------------------------------------------
    !>Read parameters
    !--------------------------------------------------
    subroutine ReadParameters()
        call ReadBaseParameters
        call ReadDISParameters
    end subroutine ReadParameters

    subroutine ReadBaseParameters()
        namelist /config_para/ Ma, Re, AOA, Tinf, Tw, &
                          ck, Pr, flag_vislaw, omega, Tcm, &
                          dt_global, CFL, Nstep_show, Nstep_write, Nstep_restart, Nstep_audit, Nstep_monitor, EPS, MAX_ITER, &
                          Quadrature_method, TimeStep_method, TimeIntegration_method, PeriodicComparisonProfile, &
                          Implicit_nCFL, Implicit_MacroSweeps, Implicit_MicroSweeps, &
                          Implicit_ResidualTol, Implicit_MomentTol, Implicit_MinAlpha, &
                          Implicit_MaxRetries, Implicit_MaxIter, UnifiedPhysicalOperators, &
                          CouetteMuRef,StructuredVelocityFile,Implicit_ClosedMassConstraint, &
                          Implicit_ZhuMovingWall,Implicit_InitHaloResync,ClosedCouetteMassRescale, &
                          Unsteady_Theta, Unsteady_InnerTol, Unsteady_MaxInner, &
                          Iflag_InitFlow, Iflag_WriteSurf, Iflag_WriteRestart, S2D1_initial_mode, &
                          S2W1_wall_mode,S2W1_time_centered,S2W1D_face_time_consistent,S2W2R_shared_predictor, &
                          S2W1_wall_ux,S2W1_wall_theta,S2PhysicalGas,S2_mu_ref, &
                          S2WallTNative,S2VelocityTable,S2VelocityScale,S2CartesianMesh,S2Diagnostics,S2CheckpointImport, &
                          Nx_all, Ny_all, Nu, Nv, U_MIN, U_MAX, V_MIN, V_MAX, NP_X, NP_Y
        namelist /boundary_para/ BC_Type, BG_geometry, FluxWall_type, ACC, WallUx, CouetteUseW2BWall

        call CheckDuplicateS2Walls()
        open(1, file="Config.ini")
        read(1, nml=config_para,iostat=ierr)
        if(ierr/=0)call AbortInvalidImplicitOption('Config.ini config namelist parse failure')
        read(1, nml=boundary_para,iostat=ierr)
        if(ierr/=0)call AbortInvalidImplicitOption('Config.ini boundary namelist parse failure')
        close(1)
        if(UnifiedPhysicalOperators)then
            S2_mu_ref=CouetteMuRef
            S2W1_wall_ux=WallUx(3:4)
        endif
        if(Nstep_show<1.or.Nstep_write<1.or.Nstep_restart<0.or.Nstep_audit<0.or.Nstep_monitor<1) &
            call AbortInvalidImplicitOption('invalid output frequency in Config.ini')
        call ResolveS2WallTemperatures()
    end subroutine ReadBaseParameters

    subroutine CheckDuplicateS2Walls()
        character(len=4096)::line
        character(len=:),allocatable::text,key
        real(KREAL)::seen(2),values(2),tn
        integer::fu,ios,p,a,z,v,q,side,last_z,last_kind
        logical::in_config
        namelist /wall_check/ values
        text='';in_config=.false.
        open(newunit=fu,file='Config.ini',status='old')
        do
            read(fu,'(A)',iostat=ios)line
            if(ios/=0)exit
            q=0
            do p=1,len_trim(line)
                if(q==0.and.(line(p:p)=="'".or.line(p:p)=='"'))then
                    q=iachar(line(p:p))
                elseif(q/=0.and.iachar(line(p:p))==q)then
                    q=0
                elseif(q==0.and.line(p:p)=='!')then
                    line=line(:p-1);exit
                endif
            enddo
            if(index(LowerS2(line),'&config_para')>0)in_config=.true.
            if(in_config)text=text//' '//trim(line)
            if(in_config.and.trim(adjustl(line))=='/')exit
        enddo
        close(fu)
        seen=-huge(1.0_KREAL);p=1;last_z=0;last_kind=0;q=0
        do while(p<=len(text))
            if(q==0.and.(text(p:p)=="'".or.text(p:p)=='"'))then
                q=iachar(text(p:p));p=p+1;cycle
            elseif(q/=0)then
                if(iachar(text(p:p))==q)q=0
                p=p+1;cycle
            endif
            if(text(p:p)=='/')exit
            if(index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_',text(p:p))==0)then
                p=p+1;cycle
            endif
            a=p
            do while(p<=len(text))
                if(index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789',text(p:p))==0)exit
                p=p+1
            enddo
            z=p-1;key=LowerS2(text(a:z))
            do while(p<=len(text))
                if(text(p:p)/=' '.and.iachar(text(p:p))/=9)exit
                p=p+1
            enddo
            if(p<=len(text))then
                if(text(p:p)=='(')then
                    v=index(text(p:),')');if(v==0)exit
                    p=p+v
                    do while(p<=len(text))
                        if(text(p:p)/=' ')exit
                        p=p+1
                    enddo
                endif
            endif
            if(p>len(text))exit
            if(text(p:p)/='=')cycle
            if(last_kind/=0)call CheckAssignment(text(last_z+1:a-1),last_kind)
            last_z=z;last_kind=0
            if(key=='s2walltnative')last_kind=1
            if(key=='s2w1_wall_theta')last_kind=2
            p=p+1
        enddo
        if(last_kind/=0)call CheckAssignment(text(last_z+1:p-1),last_kind)
    contains
        subroutine CheckAssignment(tail,which)
            character(len=*),intent(in)::tail
            integer,intent(in)::which
            character(len=:),allocatable::nml
            values=-huge(1.0_KREAL)
            nml='&wall_check values'//tail//' /'
            read(nml,nml=wall_check,iostat=ios)
            if(ios/=0)call AbortInvalidImplicitOption('cannot parse S2 wall temperature assignment')
            do side=1,2
                if(values(side)==-huge(1.0_KREAL))cycle
                tn=values(side);if(which==2)tn=2.0_KREAL*tn
                if(seen(side)/=-huge(1.0_KREAL))then
                    if(abs(tn-seen(side))>8.0_KREAL*epsilon(tn)*max(abs(tn),abs(seen(side)))) &
                        call AbortInvalidImplicitOption('repeated same-side S2 wall temperature assignment disagrees')
                endif
                seen(side)=tn
            enddo
        end subroutine
    end subroutine
    function LowerS2(text) result(lower)
        character(len=*),intent(in)::text
        character(len=len(text))::lower
        integer::i,n
        lower=text
        do i=1,len(text)
            n=iachar(text(i:i))
            if(n>=iachar('A').and.n<=iachar('Z'))lower(i:i)=achar(n+32)
        enddo
    end function

    subroutine ResolveS2WallTemperatures()
        integer :: side
        logical :: has_theta,has_native
        real(KREAL) :: tn
        do side=1,2
            has_theta=S2W1_wall_theta(side)/=-huge(1.0_KREAL)
            has_native=S2WallTNative(side)/=-huge(1.0_KREAL)
            if(has_theta.and.has_native)then
                if(abs(S2WallTNative(side)-2.0_KREAL*S2W1_wall_theta(side))> &
                    8.0_KREAL*epsilon(tn)*max(abs(S2WallTNative(side)),abs(2.0_KREAL*S2W1_wall_theta(side)))) &
                    call AbortInvalidImplicitOption('same-side S2WallTNative and S2W1_wall_theta disagree')
            endif
            tn=Tw
            if(has_theta)tn=2.0_KREAL*S2W1_wall_theta(side)
            if(has_native)tn=S2WallTNative(side)
            if(S2W1_wall_mode==1)then
                if(.not.ieee_is_finite(tn).or.tn<=0.0_KREAL) &
                    call AbortInvalidImplicitOption('effective S2 wall T_native must be positive and finite')
            elseif(.not.has_theta.and..not.has_native)then
                tn=2.0_KREAL ! unused legacy periodic defaults
            endif
            S2WallTNative(side)=tn;S2W1_wall_theta(side)=0.5_KREAL*tn
            if(S2W1_wall_mode==1)write(*,'(A,I0,A,3ES24.16)') &
                'S2 wall side ',side,' T_native, theta, lambda: ',tn,0.5_KREAL*tn,1.0_KREAL/tn
        enddo
    end subroutine ResolveS2WallTemperatures

    subroutine ReadDISParameters()
        integer(KINT) :: i

        open(1, file="Config.ini")
        call FindMarker(1, '$DIS_para')
        read(1, *) Iflag_DIS
        read(1, *) Iflag_Fringe
        read(1, *) DISBC_Type

        if (DISBC_Type == DISBC_Inflow_TS) then
            call FindMarker(1, '#Eigv_TS')
            read(1, *) Num_dis_ts
            do i = 1, Num_dis_ts
                read(1, *) wave_ts(i)%name, wave_ts(i)%omer, wave_ts(i)%omei, &
                           wave_ts(i)%alfr, wave_ts(i)%alfi, wave_ts(i)%amp
            enddo
        elseif (DISBC_Type == DISBC_Free_Wave) then
            call FindMarker(1, '#Eigv_free')
            read(1, *) Num_dis_free
            do i = 1, Num_dis_free
                read(1, *) wave_free(i)%name, wave_free(i)%omer, &
                           wave_free(i)%theta, wave_free(i)%amp
            enddo
        endif

        if (Iflag_Fringe == 1) then
            call FindMarker(1, '#Eigv_fringe')
            read(1, *) Iters_fringe, rfringe(1:4)
        endif

        close(1)
    end subroutine ReadDISParameters

    subroutine FindMarker(unit_id, marker)
        integer(KINT), intent(in)    :: unit_id
        character(len=*), intent(in) :: marker
        character(len=64)            :: line

        do
            read(unit_id, '(a)') line
            if (trim(adjustl(line)) == marker) exit
        enddo
    end subroutine FindMarker

    !--------------------------------------------------
    !>Validate options retained in the simplified branch
    !--------------------------------------------------
    pure logical function IsValidMaxwellAccommodation(value)
        real(KREAL), intent(in) :: value

        IsValidMaxwellAccommodation=ieee_is_finite(value) .and. &
            value>=0.0_KREAL .and. value<=1.0_KREAL
    end function IsValidMaxwellAccommodation

    pure logical function IsSupportedImplicitMaxwellAccommodation(method,value)
        character(len=*), intent(in) :: method
        real(KREAL), intent(in) :: value

        IsSupportedImplicitMaxwellAccommodation=.false.
        if (.not.IsValidMaxwellAccommodation(value)) return
        select case(trim(adjustl(method)))
        case('ImplicitSteady','ImplicitUnsteady')
            IsSupportedImplicitMaxwellAccommodation=.true.
        case default
            IsSupportedImplicitMaxwellAccommodation=.true.
        end select
    end function IsSupportedImplicitMaxwellAccommodation

    subroutine ValidateSupportedOptions()
        integer(KINT) :: boundary_index
        integer::rect_status
        character(len=8)::rect_flag
        logical::rect_enabled
        if(IsHighOrderImplicit().and..not.UnifiedPhysicalOperators) &
            call AbortInvalidImplicitOption('HighOrderImplicit requires UnifiedPhysicalOperators')
        if(IsHighOrderImplicit())then
            if(.not.ieee_is_finite(Implicit_nCFL).or.Implicit_nCFL<=0.0_KREAL.or. &
                .not.ieee_is_finite(Implicit_ResidualTol).or.Implicit_ResidualTol<=0.0_KREAL.or. &
                .not.ieee_is_finite(Implicit_MomentTol).or.Implicit_MomentTol<=0.0_KREAL.or. &
                .not.ieee_is_finite(Implicit_MinAlpha).or.Implicit_MinAlpha<=0.0_KREAL.or. &
                Implicit_MinAlpha>1.0_KREAL.or.Implicit_MaxIter<1) &
                call AbortInvalidImplicitOption('invalid high-order implicit tolerance, CFL or iteration limit')
        endif
        if(UnifiedPhysicalOperators)then
            if(TimeStep_method/=Time_global_varying) &
                call AbortInvalidImplicitOption('shared physical operators require global physical CFL')
            if(.not.ieee_is_finite(CouetteMuRef).or.CouetteMuRef<=0) &
                call AbortInvalidImplicitOption('shared operators require positive viscosity')
            if(IsHighOrderImplicit().and.(Implicit_MacroSweeps/=2.or.Implicit_MicroSweeps/=1.or.Implicit_MaxRetries/=0)) &
                call AbortInvalidImplicitOption('HighOrderImplicit requires PR2, one micro pair, no retry')
        endif
        if(ClosedCouetteMassRescale)then
            if(.not.IsImplicitSteady().or.any(BC_Type/=[10,10,2,2]).or.BG_geometry/=0.or. &
                .not.CouetteUseW2BWall.or.FluxWall_type/=1.or.ACC/=1.or.Iflag_DIS/=0.or. &
                Iflag_Fringe/=0.or.DISBC_Type/=0.or.Implicit_ClosedMassConstraint) &
                call AbortInvalidImplicitOption('mass rescale only for closed periodic Couette implicit')
        endif
        if(S2VelocityScale/=-1.0_KREAL)then
            if(.not.ieee_is_finite(S2VelocityScale).or.S2VelocityScale<=0.0_KREAL) &
                call AbortInvalidImplicitOption('S2VelocityScale must be positive and finite')
            if(VelocityMode()=='full16'.and.S2VelocityScale/=1.0_KREAL) &
                call AbortInvalidImplicitOption('full16 keeps its original scale 1 velocity table')
        endif
        if(slab_enabled)then
            if(NP_X/=1.or.all(NP_Y/=[1,2,4]).or.Nx_all/=8.or.all(Ny_all/=[32,64]).or. &
                (trim(TimeIntegration_method)/='S2CN2'.and..not.IsHighOrderImplicit()).or.FluxWall_type/=1.or.S2W1_wall_mode/=1.or. &
                VelocityMode()/='half32'.or.S2PhysicalGas/=1.or.ck/=1.or. &
                abs(Pr-2.0_KREAL/3.0_KREAL)>epsilon(Pr).or.trim(PeriodicComparisonProfile)/='Off')call AbortInvalidImplicitOption('MPI1 scope violation')
        endif
        if(S2PhysicalGas/=0.and.S2PhysicalGas/=1)call AbortInvalidImplicitOption('S2PhysicalGas must be 0 or 1')
        if(S2PhysicalGas==1.and.(trim(TimeIntegration_method)=='S2CN2'.or.IsHighOrderImplicit()))then
            if((trim(TimeIntegration_method)/='S2CN2'.and..not.IsHighOrderImplicit()).or.(S2W1_wall_mode==1.and.FluxWall_type/=1).or. &
                S2W2R_shared_predictor/=0.or.S2W1_time_centered/=0.or.S2W1D_face_time_consistent/=0.or. &
                all(ck/=[1,3]))call AbortInvalidImplicitOption('physical gas requires W2 type1 single-rank CK1/3')
            if(.not.ieee_is_finite(Pr).or.Pr<=0.0_KREAL.or.Pr>1.0_KREAL.or. &
                .not.ieee_is_finite(S2_mu_ref).or.S2_mu_ref<=0.0_KREAL.or. &
                .not.ieee_is_finite(omega).or..not.ieee_is_finite(Tcm).or.Tcm<0.0_KREAL.or. &
                all(flag_vislaw/=[PowerLaw,Sutherland]))call AbortInvalidImplicitOption('invalid physical transport parameters')
        endif
        if(VelocityMode()/='full16'.and.(trim(TimeIntegration_method)=='S2CN2'.or.IsHighOrderImplicit().or..not.UnifiedPhysicalOperators))then
            if(Nu/=VelocityCount().or.Nv/=VelocityCount().or.(S2PhysicalGas==0.and.(ck/=1.or.Pr/=1.0_KREAL)).or. &
                Quadrature_method/=Quad_GaussHermite.or.FluxWall_type/=1.or.(S2W1_wall_mode==1.and.any(BC_Type/=[10,10,2,2])).or. &
                Iflag_DIS/=0.or.Iflag_Fringe/=0.or.(trim(TimeIntegration_method)/='S2CN2'.and..not.IsHighOrderImplicit()).or. &
                S2W2R_shared_predictor/=0.or.S2W1_time_centered/=0.or.S2W1D_face_time_consistent/=0) &
                call AbortInvalidImplicitOption('new velocity tables require approved single-rank flat BGK range')
        endif


        if(trim(adjustl(PeriodicComparisonProfile))/='Off'.and. &
           trim(adjustl(PeriodicComparisonProfile))/='BGK_ConstTau') &
            call AbortInvalidImplicitOption('PeriodicComparisonProfile is unsupported')
        if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau') then
            if(trim(adjustl(TimeIntegration_method))/='Explicit') &
                call AbortInvalidImplicitOption('BGK_ConstTau requires native Explicit/VL/CN')
            if(NP_X/=1.or.NP_Y/=1) call AbortInvalidImplicitOption('BGK_ConstTau requires one rank')
            if(Nx_all/=Ny_all.or.all(Nx_all/=[8,16,32])) &
                call AbortInvalidImplicitOption('BGK_ConstTau requires N8/N16/N32 square mesh')
            if(Quadrature_method/=Quad_GaussHermite.or.Nu/=16.or.Nv/=16.or.ck/=1.or. &
               abs(Pr-1.0_KREAL)>epsilon(Pr)) &
                call AbortInvalidImplicitOption('BGK_ConstTau requires BGK CK1 frozen GH16xGH16')
            if((TimeStep_method/=Time_global_fixed.and..not.UnifiedPhysicalOperators).or..not.ieee_is_finite(dt_global).or.dt_global<=0.0_KREAL) &
                call AbortInvalidImplicitOption('BGK_ConstTau requires positive fixed dt')
            if(any(BC_Type/=BCconst_periodic).or.Iflag_DIS/=0.or.Iflag_Fringe/=0) &
                call AbortInvalidImplicitOption('BGK_ConstTau requires double periodic and DIS/fringe off')
            if(S2D1_initial_mode/=0.and.S2D1_initial_mode/=3) &
                call AbortInvalidImplicitOption('BGK_ConstTau initial mode unsupported')
            if(MAX_ITER<1) call AbortInvalidImplicitOption('BGK_ConstTau requires MAX_ITER>=1')
        endif

        select case (trim(adjustl(TimeIntegration_method)))
        case ('Explicit', 'ImplicitSteady', 'ImplicitUnsteady', 'S2D1', 'S2D1VL', 'S2CN2', 'S2Steady','HighOrderImplicit')
        case default
            call AbortInvalidImplicitOption('TimeIntegration_method is unsupported')
        end select
        if(S2W1_wall_mode/=0.and.(trim(adjustl(TimeIntegration_method))/='S2CN2'.and..not.IsHighOrderImplicit()).and..not.UnifiedPhysicalOperators) &
            call AbortInvalidImplicitOption('S2W1 wall mode requires S2CN2')
        if(S2W1_time_centered<0.or.S2W1_time_centered>1) &
            call AbortInvalidImplicitOption('S2W1_time_centered must be 0 or 1')
        if(S2W1_time_centered==1.and.S2W1_wall_mode/=1) &
            call AbortInvalidImplicitOption('S2W1 midpoint time option requires the flat wall')
        if(S2W1D_face_time_consistent<0.or.S2W1D_face_time_consistent>1) &
            call AbortInvalidImplicitOption('S2W1D_face_time_consistent must be 0 or 1')
        if(S2W1D_face_time_consistent==1.and.(S2W1_wall_mode/=1.or.S2W1_time_centered/=1)) &
            call AbortInvalidImplicitOption('S2W1D face correction requires W1T midpoint wall')

        if(S2W2R_shared_predictor<0.or.S2W2R_shared_predictor>1) &
            call AbortInvalidImplicitOption('S2W2R_shared_predictor must be 0 or 1')
        if(S2W2R_shared_predictor==1)then
            if((trim(adjustl(TimeIntegration_method))/='S2CN2'.and..not.IsHighOrderImplicit()).or.S2W1_wall_mode/=1.or.FluxWall_type/=1.or. &
                ck/=1.or.Pr/=1.0_KREAL.or.S2W1_time_centered/=0.or.S2W1D_face_time_consistent/=0.or. &
                NP_X/=1.or.NP_Y/=1)call AbortInvalidImplicitOption('W2R supports only frozen flat type1 BGK setup')
        endif
        select case (Quadrature_method)
        case (Quad_NewtonCotes4, Quad_GaussHermite, Quad_Unstruct)
        case default
            write(*, '(a,I0)') 'Error: Quadrature_method must be 1, 2, or 3. Current value=', Quadrature_method
            stop
        end select

        select case (FluxWall_type)
        case (FluxWall_Kinetic, FluxWall_Multiscale, FluxWall_Maxwell)
        case default
            write(*, '(a,I0)') 'Error: FluxWall_type must be 0, 1, or 2. Current value=', FluxWall_type
            stop
        end select

        if (FluxWall_type==FluxWall_Maxwell) then
            if (.not.IsValidMaxwellAccommodation(ACC)) &
                call AbortInvalidImplicitOption('FluxWall_Maxwell requires finite ACC in [0,1]')
            if (.not.IsSupportedImplicitMaxwellAccommodation(TimeIntegration_method,ACC)) &
                call AbortInvalidImplicitOption('FluxWall_Maxwell is unsupported for the selected time integration method')
        endif

        if((trim(adjustl(TimeIntegration_method))=='S2CN2'.or.IsHighOrderImplicit()).or.IsS2Steady()) then
            if(NP_X<1.or.NP_Y<1)call AbortInvalidImplicitOption('S2CN2 requires positive NP_X/NP_Y')
            if(IsS2Steady())then
            rect_flag='';call get_environment_variable('S2_RECT_WALL',rect_flag,status=rect_status)
            rect_enabled=S2CartesianMesh.or.(rect_status==0.and.trim(rect_flag)=='1')
            if(Nx_all/=Ny_all)then
                if(.not.rect_enabled.or.all(Nx_all/=[8,16]).or.all(Ny_all/=[16,32,64]).or. &
                    (trim(adjustl(TimeIntegration_method))/='S2CN2'.and..not.IsHighOrderImplicit()).or.S2W1_wall_mode/=1.or. &
                    FluxWall_type/=1.or.S2W2R_shared_predictor/=0.or.S2W1_time_centered/=0.or. &
                    S2W1D_face_time_consistent/=0.or.(S2PhysicalGas==0.and.abs(Pr-1.0_KREAL)>epsilon(Pr))) &
                    call AbortInvalidImplicitOption('rectangle requires enabled W2 type1 Nx8/16 Ny16/32/64 BGK')
            elseif(all(Nx_all/=[8,16,32,64]))then
                call AbortInvalidImplicitOption('S2CN2 requires N8/N16/N32/N64 square mesh')
            endif
            endif
            if(Quadrature_method/=Quad_GaussHermite.or.Nu/=VelocityCount().or.Nv/=VelocityCount()) &
                call AbortInvalidImplicitOption('S2CN2 velocity count disagrees with selected table')
            if(ck/=1.and.S2PhysicalGas==0) call AbortInvalidImplicitOption('S2CN2 legacy requires CK=1')
            if(S2PhysicalGas==0.and.abs(Pr-1.0_KREAL)>epsilon(Pr).and.abs(Pr-2.0_KREAL/3.0_KREAL)>epsilon(Pr)) &
                call AbortInvalidImplicitOption('S2CN2 supports only Pr=1 or Pr=2/3 in this package')
            if((TimeStep_method/=Time_global_fixed.and..not.UnifiedPhysicalOperators).or..not.ieee_is_finite(dt_global).or.dt_global<=0.0_KREAL) &
                call AbortInvalidImplicitOption('S2CN2 requires positive fixed dt_global')
            if(S2W1_wall_mode==0)then
                if(any(BC_Type/=BCconst_periodic)) call AbortInvalidImplicitOption('S2CN2 requires four periodic boundaries')
            elseif(S2W1_wall_mode==1.and.(trim(adjustl(TimeIntegration_method))=='S2CN2'.or.IsHighOrderImplicit()))then
                if(any(BC_Type(1:2)/=BCconst_periodic).or.any(BC_Type(3:4)/=BCconst_wall_isothermal)) &
                    call AbortInvalidImplicitOption('S2W1 requires x periodic and two y isothermal walls')
                if((FluxWall_type/=FluxWall_Kinetic.and.FluxWall_type/=FluxWall_Multiscale).or. &
                    (S2PhysicalGas==0.and.abs(Pr-1.0_KREAL)>epsilon(Pr)).or. &
                    any(S2W1_wall_theta<=0.0_KREAL).or..not.all(ieee_is_finite(S2W1_wall_theta)).or. &
                    .not.all(ieee_is_finite(S2W1_wall_ux))) &
                    call AbortInvalidImplicitOption('S2 flat walls require kinetic/multiscale BGK and finite wall state')
                if(FluxWall_type==FluxWall_Multiscale.and. &
                    (S2W1_time_centered/=0.or.S2W1D_face_time_consistent/=0)) &
                    call AbortInvalidImplicitOption('S2W2 type1 requires centered=0 and face_time_consistent=0')
            else
                call AbortInvalidImplicitOption('S2W1 wall mode is unsupported for this method')
            endif
            if(Iflag_DIS/=0.or.Iflag_Fringe/=0) call AbortInvalidImplicitOption('S2CN2 does not support DIS/fringe')
            if(.not.IsS2Steady())then
                ! Native S2 cold/restart initialization is selected by Iflag_InitFlow.
                continue
            elseif(S2W1_wall_mode==1)then
                if(S2D1_initial_mode/=4.and.S2D1_initial_mode/=5) &
                    call AbortInvalidImplicitOption('S2W1 requires initial mode 4 or 5')
            else
                if(S2D1_initial_mode/=0.and.S2D1_initial_mode/=3) &
                    call AbortInvalidImplicitOption('S2CN2 supports initial modes 0 and 3 only')
            endif
            if(MAX_ITER<1) call AbortInvalidImplicitOption('S2CN2 requires MAX_ITER>=1')
            if(IsS2Steady())then
                if(Nx_all/=8.or.abs(Pr-1.0_KREAL)>epsilon(Pr)) &
                    call AbortInvalidImplicitOption('S2Steady I1 requires N8 BGK')
                if(.not.ieee_is_finite(CFL).or.CFL<=0.0_KREAL.or.CFL>1.0_KREAL) &
                    call AbortInvalidImplicitOption('S2Steady requires valid transport CFL')
                if(.not.ieee_is_finite(Implicit_nCFL).or.Implicit_nCFL<=0.0_KREAL.or. &
                   Implicit_MacroSweeps/=2.or.Implicit_MicroSweeps/=1) &
                    call AbortInvalidImplicitOption('S2Steady I1 requires nCFL and 2/1 sweeps')
                if(Implicit_MaxIter<1.or.Implicit_MaxIter>1000) &
                    call AbortInvalidImplicitOption('S2Steady I1 iteration limit out of range')
            endif
            return
        endif

        if (IsImplicitUnsteady()) then
            if (Quadrature_method /= Quad_NewtonCotes4 .and. Quadrature_method /= Quad_GaussHermite) &
                call AbortInvalidImplicitOption('ImplicitUnsteady supports Quadrature_method 1 or 2 only')
            if (TimeStep_method /= Time_global_fixed) &
                call AbortInvalidImplicitOption('ImplicitUnsteady requires TimeStep_method=Time_global_fixed')
            if (.not.ieee_is_finite(dt_global) .or. dt_global<=0.0_KREAL) &
                call AbortInvalidImplicitOption('ImplicitUnsteady requires positive dt_global')
            if (.not.ieee_is_finite(CFL) .or. CFL<=0.0_KREAL .or. CFL>1.0_KREAL) &
                call AbortInvalidImplicitOption('ImplicitUnsteady flux CFL must lie in (0,1]')
            if (.not.ieee_is_finite(Unsteady_Theta) .or. Unsteady_Theta<0.5_KREAL .or. Unsteady_Theta>1.0_KREAL) &
                call AbortInvalidImplicitOption('Unsteady_Theta must lie in [0.5,1]')
            if (.not.ieee_is_finite(Unsteady_InnerTol) .or. Unsteady_InnerTol<=0.0_KREAL .or. Unsteady_MaxInner<1) &
                call AbortInvalidImplicitOption('unsteady inner tolerance and iteration limit must be positive')
            return
        endif

        if (.not.IsImplicitSteady()) return

        if (Quadrature_method /= Quad_NewtonCotes4 .and. Quadrature_method /= Quad_GaussHermite) &
            call AbortInvalidImplicitOption('ImplicitSteady supports Quadrature_method 1 or 2 only')
        if (Iflag_DIS /= 0) call AbortInvalidImplicitOption('Iflag_DIS must equal 0')
        if (Iflag_Fringe /= 0) call AbortInvalidImplicitOption('Iflag_Fringe must equal 0')
        if (.not.ieee_is_finite(CFL) .or. CFL <= 0.0_KREAL .or. CFL > 1.0_KREAL) &
            call AbortInvalidImplicitOption('implicit physical CFL must lie in (0,1]')
        if (.not.ieee_is_finite(Implicit_nCFL) .or. Implicit_nCFL <= 0.0_KREAL) &
            call AbortInvalidImplicitOption('Implicit_nCFL must be positive')
        if (TimeStep_method==Time_global_fixed) &
            call AbortInvalidImplicitOption('fixed global dt is incompatible with implicit pCFL/nCFL separation')
        if (Implicit_MacroSweeps < 1 .or. Implicit_MicroSweeps < 1) &
            call AbortInvalidImplicitOption('implicit sweep counts must be positive')
        if (.not.ieee_is_finite(Implicit_MinAlpha) .or. &
            Implicit_MinAlpha <= 0.0_KREAL .or. Implicit_MinAlpha > 1.0_KREAL) &
            call AbortInvalidImplicitOption('Implicit_MinAlpha must lie in (0,1]')
        if (.not.ieee_is_finite(Implicit_ResidualTol) .or. .not.ieee_is_finite(Implicit_MomentTol) .or. &
            Implicit_ResidualTol <= 0.0_KREAL .or. Implicit_MomentTol <= 0.0_KREAL) &
            call AbortInvalidImplicitOption('implicit tolerances must be positive')
        if (Implicit_MaxRetries < 0 .or. Implicit_MaxIter < 1) &
            call AbortInvalidImplicitOption('implicit retry and iteration limits are invalid')

        do boundary_index = 1, size(BC_Type)
            select case (BC_Type(boundary_index))
            case (BCconst_wall_isothermal, BCconst_symmetry, BCconst_freestream, BCconst_extrapolation,BCconst_periodic)
            case default
                call AbortInvalidImplicitOption('boundary types must be 2, 3, 4, or 6')
            end select
        enddo

    end subroutine ValidateSupportedOptions

    subroutine AbortInvalidImplicitOption(message)
        character(len=*), intent(in) :: message

        write(*,'(A,I0,2A)') 'Invalid configuration on rank ', myID, ': ', trim(message)
        call MPI_Abort(MPI_COMM_WORLD, 200, ierr)
    end subroutine AbortInvalidImplicitOption

    !--------------------------------------------------
    !>Root (0) processor bcast all parameters
    !--------------------------------------------------
    subroutine bcast_para
        implicit none
        integer(KINT) :: Ipara(200), ierr
        real(KREAL)   :: Rpara(200)
        real(KREAL)   :: Ipara_real(200)

        Ipara = 0
        Rpara = 0.d0
        Ipara_real = 0.d0

        call PackIntegerParameters(Ipara)
        call PackRealParameters(Rpara)

        Ipara_real = real(Ipara, KREAL)
        call MPI_bcast(Rpara, 200, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_bcast(Ipara_real, 200, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(TimeIntegration_method, len(TimeIntegration_method), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(PeriodicComparisonProfile, len(PeriodicComparisonProfile), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        Ipara = int(Ipara_real, KINT)

        call MPI_Bcast(CouetteMuRef,1,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(CouetteUseW2BWall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(Implicit_ClosedMassConstraint,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(Implicit_ZhuMovingWall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(Implicit_InitHaloResync,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ClosedCouetteMassRescale,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(WallUx,4,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(StructuredVelocityFile,len(StructuredVelocityFile),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(UnifiedPhysicalOperators,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(Nstep_restart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(Nstep_audit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(Nstep_monitor,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call UnpackIntegerParameters(Ipara)
        call UnpackRealParameters(Rpara)
        call MPI_Bcast(S2VelocityTable,len(S2VelocityTable),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(S2VelocityScale,1,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(S2WallTNative,2,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(S2CartesianMesh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(S2Diagnostics,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(S2CheckpointImport,len(S2CheckpointImport),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    end subroutine bcast_para

    subroutine PackIntegerParameters(Ipara)
        integer(KINT), intent(out) :: Ipara(200)
        integer(KINT)              :: i, is

        Ipara(1) = ck
        Ipara(2) = flag_vislaw
        Ipara(3) = Nstep_show
        Ipara(4) = Nstep_write
        Ipara(5) = MAX_ITER
        Ipara(7) = Quadrature_method
        Ipara(8) = TimeStep_method
        Ipara(10) = Iflag_InitFlow
        Ipara(11) = Iflag_WriteSurf
        Ipara(12) = Iflag_WriteRestart
        Ipara(13) = Nx_all
        Ipara(14) = Ny_all
        Ipara(15) = Nu
        Ipara(16) = Nv
        Ipara(17) = NP_X
        Ipara(18) = NP_Y
        Ipara(19:22) = BC_Type(1:4)
        Ipara(23) = BG_geometry
        Ipara(24) = FluxWall_type
        Ipara(25) = Implicit_MacroSweeps
        Ipara(26) = Implicit_MicroSweeps
        Ipara(27) = Implicit_MaxRetries
        Ipara(28) = Implicit_MaxIter
        Ipara(29) = Unsteady_MaxInner
        Ipara(30) = S2D1_initial_mode
        Ipara(31) = S2W1_wall_mode
        Ipara(32) = S2W1_time_centered
        Ipara(33) = S2W1D_face_time_consistent
        Ipara(34) = S2W2R_shared_predictor
        Ipara(35) = S2PhysicalGas
        Ipara(40) = Iflag_DIS
        Ipara(41) = DISBC_Type
        Ipara(42) = Iflag_Fringe

        if (Iflag_DIS == 1) then
            is = 43
            if (DISBC_Type == DISBC_Inflow_TS) then
                Ipara(is) = Num_dis_ts
                is = is + 1
            elseif (DISBC_Type == DISBC_Free_Wave) then
                Ipara(is) = Num_dis_free
                do i = 1, Num_dis_free
                    Ipara(is + i) = EncodeFreeWaveName(wave_free(i)%name)
                enddo
                is = is + Num_dis_free + 1
            endif
            if (Iflag_Fringe == 1) Ipara(is) = Iters_fringe
        endif
    end subroutine PackIntegerParameters

    subroutine PackRealParameters(Rpara)
        real(KREAL), intent(out) :: Rpara(200)
        integer(KINT)            :: i, is

        Rpara(1) = Ma
        Rpara(2) = Re
        Rpara(3) = AOA
        Rpara(4) = Tinf
        Rpara(5) = Tw
        Rpara(6) = Pr
        Rpara(7) = omega
        Rpara(8) = Tcm
        Rpara(9) = dt_global
        Rpara(10) = CFL
        Rpara(11) = EPS
        Rpara(12) = U_MIN
        Rpara(13) = U_MAX
        Rpara(14) = V_MIN
        Rpara(15) = V_MAX
        Rpara(16) = ACC
        Rpara(17) = Implicit_nCFL
        Rpara(18) = Implicit_ResidualTol
        Rpara(19) = Implicit_MomentTol
        Rpara(20) = Implicit_MinAlpha
        Rpara(21) = Unsteady_Theta
        Rpara(22) = Unsteady_InnerTol
        Rpara(23:24) = S2W1_wall_ux
        Rpara(25:26) = S2W1_wall_theta
        Rpara(27) = S2_mu_ref

        if (Iflag_DIS == 1) then
            is = 40
            if (DISBC_Type == DISBC_Inflow_TS) then
                do i = 1, Num_dis_ts
                    Rpara(is + (i - 1) * 5) = wave_ts(i)%omer
                    Rpara(is + 1 + (i - 1) * 5) = wave_ts(i)%omei
                    Rpara(is + 2 + (i - 1) * 5) = wave_ts(i)%alfr
                    Rpara(is + 3 + (i - 1) * 5) = wave_ts(i)%alfi
                    Rpara(is + 4 + (i - 1) * 5) = wave_ts(i)%amp
                enddo
                is = is + 5 * Num_dis_ts
            elseif (DISBC_Type == DISBC_Free_Wave) then
                do i = 1, Num_dis_free
                    Rpara(is + (i - 1) * 3) = wave_free(i)%omer
                    Rpara(is + 1 + (i - 1) * 3) = wave_free(i)%theta
                    Rpara(is + 2 + (i - 1) * 3) = wave_free(i)%amp
                enddo
                is = is + 3 * Num_dis_free
            endif
            if (Iflag_Fringe == 1) Rpara(is:is + 3) = rfringe(1:4)
        endif
    end subroutine PackRealParameters

    subroutine UnpackIntegerParameters(Ipara)
        integer(KINT), intent(in) :: Ipara(200)
        integer(KINT)             :: i, is

        ck = Ipara(1)
        flag_vislaw = Ipara(2)
        Nstep_show = Ipara(3)
        Nstep_write = Ipara(4)
        MAX_ITER = Ipara(5)
        Quadrature_method = Ipara(7)
        TimeStep_method = Ipara(8)
        Iflag_InitFlow = Ipara(10)
        Iflag_WriteSurf = Ipara(11)
        Iflag_WriteRestart = Ipara(12)
        Nx_all = Ipara(13)
        Ny_all = Ipara(14)
        Nu = Ipara(15)
        Nv = Ipara(16)
        NP_X = Ipara(17)
        NP_Y = Ipara(18)
        BC_Type(1:4) = Ipara(19:22)
        BG_geometry = Ipara(23)
        FluxWall_type = Ipara(24)
        Implicit_MacroSweeps = Ipara(25)
        Implicit_MicroSweeps = Ipara(26)
        Implicit_MaxRetries = Ipara(27)
        Implicit_MaxIter = Ipara(28)
        Unsteady_MaxInner = Ipara(29)
        S2D1_initial_mode = Ipara(30)
        S2W1_wall_mode = Ipara(31)
        S2W1_time_centered = Ipara(32)
        S2W1D_face_time_consistent = Ipara(33)
        S2W2R_shared_predictor = Ipara(34)
        S2PhysicalGas = Ipara(35)
        Iflag_DIS = Ipara(40)
        DISBC_Type = Ipara(41)
        Iflag_Fringe = Ipara(42)

        if (Iflag_DIS == 1) then
            is = 43
            if (DISBC_Type == DISBC_Inflow_TS) then
                Num_dis_ts = Ipara(is)
                wave_ts(:)%name = 'TS'
                is = is + 1
            elseif (DISBC_Type == DISBC_Free_Wave) then
                Num_dis_free = Ipara(is)
                do i = 1, Num_dis_free
                    wave_free(i)%name = DecodeFreeWaveName(Ipara(is + i))
                enddo
                is = is + Num_dis_free + 1
            endif
            if (Iflag_Fringe == 1) Iters_fringe = Ipara(is)
        endif
    end subroutine UnpackIntegerParameters

    subroutine UnpackRealParameters(Rpara)
        real(KREAL), intent(in) :: Rpara(200)
        integer(KINT)           :: i, is

        Ma = Rpara(1)
        Re = Rpara(2)
        AOA = Rpara(3)
        Tinf = Rpara(4)
        Tw = Rpara(5)
        Pr = Rpara(6)
        omega = Rpara(7)
        Tcm = Rpara(8)
        dt_global = Rpara(9)
        CFL = Rpara(10)
        EPS = Rpara(11)
        U_MIN = Rpara(12)
        U_MAX = Rpara(13)
        V_MIN = Rpara(14)
        V_MAX = Rpara(15)
        ACC = Rpara(16)
        Implicit_nCFL = Rpara(17)
        Implicit_ResidualTol = Rpara(18)
        Implicit_MomentTol = Rpara(19)
        Implicit_MinAlpha = Rpara(20)
        Unsteady_Theta = Rpara(21)
        Unsteady_InnerTol = Rpara(22)
        S2W1_wall_ux = Rpara(23:24)
        S2W1_wall_theta = Rpara(25:26)
        S2_mu_ref = Rpara(27)

        if (Iflag_DIS == 1) then
            is = 40
            if (DISBC_Type == DISBC_Inflow_TS) then
                do i = 1, Num_dis_ts
                    wave_ts(i)%omer = Rpara(is + (i - 1) * 5)
                    wave_ts(i)%omei = Rpara(is + 1 + (i - 1) * 5)
                    wave_ts(i)%alfr = Rpara(is + 2 + (i - 1) * 5)
                    wave_ts(i)%alfi = Rpara(is + 3 + (i - 1) * 5)
                    wave_ts(i)%amp = Rpara(is + 4 + (i - 1) * 5)
                enddo
                is = is + 5 * Num_dis_ts
            elseif (DISBC_Type == DISBC_Free_Wave) then
                do i = 1, Num_dis_free
                    wave_free(i)%omer = Rpara(is + (i - 1) * 3)
                    wave_free(i)%theta = Rpara(is + 1 + (i - 1) * 3)
                    wave_free(i)%amp = Rpara(is + 2 + (i - 1) * 3)
                enddo
                is = is + 3 * Num_dis_free
            endif
            if (Iflag_Fringe == 1) rfringe(1:4) = Rpara(is:is + 3)
        endif
    end subroutine UnpackRealParameters

    integer(KINT) function EncodeFreeWaveName(name)
        character(len=*), intent(in) :: name

        select case (trim(adjustl(name)))
        case ('F')
            EncodeFreeWaveName = 1
        case ('S')
            EncodeFreeWaveName = 2
        case ('E')
            EncodeFreeWaveName = 3
        case ('V')
            EncodeFreeWaveName = 4
        case default
            EncodeFreeWaveName = 0
        end select
    end function EncodeFreeWaveName

    character(len=2) function DecodeFreeWaveName(code)
        integer(KINT), intent(in) :: code

        select case (code)
        case (1)
            DecodeFreeWaveName = 'F'
        case (2)
            DecodeFreeWaveName = 'S'
        case (3)
            DecodeFreeWaveName = 'E'
        case (4)
            DecodeFreeWaveName = 'V'
        case default
            DecodeFreeWaveName = '?'
        end select
    end function DecodeFreeWaveName

    !--------------------------------------------------
    !>Calculate constant parameters
    !--------------------------------------------------
    subroutine set_const_para
        implicit none
        real(KREAL) :: Vinf, uf, vf

        gamma = real(CK + 4, KREAL) / real(CK + 2, KREAL)

        mu_ref = sqrt(gamma / 2.d0) * Ma / Re
        Kn = sqrt(pi * gamma / 2.d0) * Ma / Re
        if(S2PhysicalGas==1)then
            mu_ref=S2_mu_ref
            Kn=sqrt(pi)*mu_ref
        endif


        if(CouetteMuRef>0)then
            mu_ref=CouetteMuRef;Kn=sqrt(PI)*mu_ref
        endif
        Vinf = Ma * sqrt(gamma / 2.d0)
        uf = Vinf * cosd(AOA)
        vf = Vinf * sind(AOA)
        Prim_inf = [1.d0, uf, vf, 1.d0 / Tinf]

        Prim_wall = [1.d0, 0.d0, 0.d0, 1.d0 / Tw]

        if(S2W2R_shared_predictor<0.or.S2W2R_shared_predictor>1) &
            call AbortInvalidImplicitOption('S2W2R_shared_predictor must be 0 or 1')
        if(S2W2R_shared_predictor==1)then
            if((trim(adjustl(TimeIntegration_method))/='S2CN2'.and..not.IsHighOrderImplicit()).or.S2W1_wall_mode/=1.or.FluxWall_type/=1.or. &
                ck/=1.or.Pr/=1.0_KREAL.or.S2W1_time_centered/=0.or.S2W1D_face_time_consistent/=0.or. &
                NP_X/=1.or.NP_Y/=1)call AbortInvalidImplicitOption('W2R supports only frozen flat type1 BGK setup')
        endif
        select case (Quadrature_method)
        case (Quad_NewtonCotes4)
            Nu = (Nu / 4) * 4 + 1
            Nv = (Nv / 4) * 4 + 1
        case (Quad_GaussHermite)
            if(len_trim(StructuredVelocityFile)>0)then
                continue
            elseif(VelocityMode()/='full16')then
                if(Nu/=VelocityCount().or.Nv/=VelocityCount())call AbortInvalidImplicitOption('selected table count mismatch')
                U_MIN=-VelocityMax();U_MAX=-U_MIN;V_MIN=U_MIN;V_MAX=U_MAX
            elseif (Nu >= 28) then
                Nu = 28
                Nv = 28
                U_MIN = -0.5392407922630E+01
                U_MAX = +0.5392407922630E+01
                V_MIN = U_MIN
                V_MAX = U_MAX
            elseif (Nu >= 16) then
                Nu = 16
                Nv = 16
                U_MIN = -0.3686007162724397E+1
                U_MAX = 0.3686007162724397E+1
                V_MIN = U_MIN
                V_MAX = U_MAX
            elseif (Nu >= 8) then
                Nu = 8
                Nv = 8
                U_MIN = -0.2262664477010362E+1
                U_MAX = 0.2262664477010362E+1
                V_MIN = U_MIN
                V_MAX = U_MAX
            endif
        case (Quad_Unstruct)
            if (myid == 0) then
                open(10, file='microProperties2d.dat', status='old', action='read')
                read(10, *)
                read(10, *) Nu
                close(10)
            endif

            call MPI_BCAST(Nu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            Nv = 1
            if (myid == 0) print *, "Unstructured Velocity Grid detected: Nu=", Nu, " Nv=", Nv
        case default
            print *, 'Error: the quadrature method is illegal: Quadrature_method =', Quadrature_method
            stop
        end select
    end subroutine set_const_para

    !--------------------------------------------------
    !>Write parameters (log file)
    !--------------------------------------------------
    subroutine write_log
        implicit none
        integer :: log_unit
        type(MethodName) :: TimeStep_name(0:20), Quadrature_name(0:20)
        type(MethodName) :: Geometry_name(0:20), FluxBC_name(0:20), Boundary_name(0:100)
        type(MethodName) :: Initflow_name(0:20), VisLaw_name(0:20), DISBC_name(0:20), Option_name(0:1)

        call InitializeMethodNames(TimeStep_name, Quadrature_name, &
                                   Geometry_name, FluxBC_name, Boundary_name, Initflow_name, &
                                   VisLaw_name, DISBC_name, Option_name)

        call WriteParameterReport(6, .true., TimeStep_name, Quadrature_name, &
                                  Geometry_name, FluxBC_name, Boundary_name, Initflow_name, &
                                  VisLaw_name, DISBC_name, Option_name)

        log_unit = 99
        open(log_unit, file="Simulation_Log.txt", action="write")
        call WriteParameterReport(log_unit, .false., TimeStep_name, Quadrature_name, &
                                  Geometry_name, FluxBC_name, Boundary_name, Initflow_name, &
                                  VisLaw_name, DISBC_name, Option_name)
        close(log_unit)
    end subroutine write_log

    subroutine InitializeMethodNames(TimeStep_name, Quadrature_name, &
                                     Geometry_name, FluxBC_name, Boundary_name, Initflow_name, &
                                     VisLaw_name, DISBC_name, Option_name)
        type(MethodName), intent(out) :: TimeStep_name(0:20), Quadrature_name(0:20)
        type(MethodName), intent(out) :: Geometry_name(0:20), FluxBC_name(0:20), Boundary_name(0:100)
        type(MethodName), intent(out) :: Initflow_name(0:20), VisLaw_name(0:20), DISBC_name(0:20), Option_name(0:1)

        TimeStep_name(Time_global_fixed)%name = "Time_global_fixed"
        TimeStep_name(Time_global_varying)%name = "Time_global_varying"
        TimeStep_name(Time_local_varying)%name = "Time_local_varying"
        Quadrature_name(Quad_NewtonCotes4)%name = "Quad_NewtonCotes4"
        Quadrature_name(Quad_GaussHermite)%name = "Quad_GaussHermite"
        Quadrature_name(Quad_Unstruct)%name = "Quad_Unstruct"
        Geometry_name(BGconst_General)%name = "BGconst_General"
        Geometry_name(BGconst_FlatPlate)%name = "BGconst_FlatPlate"
        FluxBC_name(FluxWall_Kinetic)%name = "FluxWall_Kinetic"
        FluxBC_name(FluxWall_Multiscale)%name = "FluxWall_Multiscale"
        FluxBC_name(FluxWall_Maxwell)%name = "FluxWall_Maxwell"
        Boundary_name(BCconst_wall_isothermal)%name = "Wall_isothermal"
        Boundary_name(BCconst_wall_insulate)%name = "Wall_insulate"
        Boundary_name(BCconst_freestream)%name = "Freestream"
        Boundary_name(BCconst_symmetry)%name = "Symmetry"
        Boundary_name(BCconst_inflow)%name = "Inflow"
        Boundary_name(BCconst_extrapolation)%name = "Extrapolation"
        Boundary_name(BCconst_periodic)%name = "Periodic"
        Initflow_name(Init_Freestream)%name = "Init_Freestream"
        Initflow_name(Init_Macroflow)%name = "Init_Macroflow"
        Initflow_name(Init_Restart)%name = "Init_Restart"
        VisLaw_name(PowerLaw)%name = "PowerLaw"
        VisLaw_name(Sutherland)%name = "Sutherland"
        DISBC_name(DISBC_None)%name = "None"
        DISBC_name(DISBC_Inflow_TS)%name = "Inflow_TS"
        DISBC_name(DISBC_Free_wave)%name = "Free_wave"
        Option_name(0)%name = "No"
        Option_name(1)%name = "Yes"
    end subroutine InitializeMethodNames

    subroutine WriteParameterReport(unit_id, is_screen, TimeStep_name, Quadrature_name, &
                                    Geometry_name, FluxBC_name, Boundary_name, Initflow_name, &
                                    VisLaw_name, DISBC_name, Option_name)
        integer, intent(in)          :: unit_id
        logical, intent(in)          :: is_screen
        type(MethodName), intent(in) :: TimeStep_name(0:20), Quadrature_name(0:20)
        type(MethodName), intent(in) :: Geometry_name(0:20), FluxBC_name(0:20), Boundary_name(0:100)
        type(MethodName), intent(in) :: Initflow_name(0:20), VisLaw_name(0:20), DISBC_name(0:20), Option_name(0:1)
        integer :: i

        if (is_screen) then
            write(unit_id, '(a)') "-------------------------------Parameters check----------------------------------"
        else
            write(unit_id, *)
            write(unit_id, '(a)') "======================UGKS2D parallel computation (Ver 2.0)======================="
            write(unit_id, '(a)') "--------------------------------Control parameters--------------------------------"
            write(unit_id, '(2(a15,a10))') "Date=", date, "Time=", time
        endif

        write(unit_id, '(2x,3(a15,f15.8))') "Ma=", Ma, "Re=", Re, "Kn=", Kn
        write(unit_id, '(2x,3(a15,f15.8))') "AOA=", AOA, "Tinf=", Tinf, "Tw=", Tw
        write(unit_id, '(2x,3(a15,f15.8))') "ck=", real(ck), "gamma=", gamma, "Pr=", Pr
        write(unit_id, '(2x,2(a15,f15.8))') "omega=", omega, "Tcm=", Tcm
        write(unit_id, '(2x,2(a15,I15))') "Nstep_show=", Nstep_show, "Nstep_write=", Nstep_write
        write(unit_id, '(2x,3(a15,e15.8))') "CFL=", CFL, "EPS=", EPS, "dt_global=", dt_global
        write(unit_id, '(2x,4(a15,I15))') "Nx_all=", Nx_all, "Ny_all=", Ny_all, "Nu=", Nu, "Nv=", Nv
        write(unit_id, '(2x,2(a15,I15))') "NP_X=", NP_X, "NP_Y=", NP_Y
        write(unit_id, '(2x,4(a15,f15.8))') "U_MIN=", U_MIN, "U_MAX=", U_MAX, "V_MIN=", U_MIN, "V_MAX=", V_MAX
        write(unit_id, '(2x,2a30)') "flag_vislaw  =  ", VisLaw_name(flag_vislaw)
        write(unit_id, '(2x,2a30)') "TimeStep_method  =  ", TimeStep_name(TimeStep_method)
        write(unit_id, '(2x,2a30)') "TimeIntegration_method  =  ", trim(TimeIntegration_method)
        write(unit_id, '(2x,a30,1p,e15.8)') "Implicit_nCFL  =  ", Implicit_nCFL
        write(unit_id, '(2x,2(a30,I8))') "Implicit_MacroSweeps  =  ", Implicit_MacroSweeps, &
                                        "Implicit_MicroSweeps  =  ", Implicit_MicroSweeps
        write(unit_id, '(2x,3(a30,1p,e15.8))') "Implicit_ResidualTol  =  ", Implicit_ResidualTol, &
                                                  "Implicit_MomentTol  =  ", Implicit_MomentTol, &
                                                  "Implicit_MinAlpha  =  ", Implicit_MinAlpha
        write(unit_id, '(2x,2(a30,I8))') "Implicit_MaxRetries  =  ", Implicit_MaxRetries, &
                                        "Implicit_MaxIter  =  ", Implicit_MaxIter
        write(unit_id, '(2x,2(a30,1p,e15.8),a30,I8)') "Unsteady_Theta  =  ", Unsteady_Theta, &
                                        "Unsteady_InnerTol  =  ", Unsteady_InnerTol, &
                                        "Unsteady_MaxInner  =  ", Unsteady_MaxInner
        write(unit_id, '(2x,2a30)') "Reconstruction_method  =  ", "Reconst_Vanleer"
        write(unit_id, '(2x,2a30)') "Quadrature_method  =  ", Quadrature_name(Quadrature_method)
        write(unit_id, '(2x,2a30)') "Iflag_InitFlow  =  ", Initflow_name(Iflag_InitFlow)
        write(unit_id, '(a)') "-------------------------------Boundary information------------------------------"
        write(unit_id, '(2x,2a30)') "BG_geometry  =  ", Geometry_name(BG_geometry)
        write(unit_id, '(2x,2a30)') "FluxWall_type  =  ", FluxBC_name(FluxWall_type)
        if (is_screen) then
            write(unit_id, '(2x,a30,f15.8)') "Acc coeff.  =  ", ACC
        else
            write(unit_id, '(2x,a30,f12.8)') "Acc coeff.  =  ", ACC
        endif
        write(unit_id, '(2x,2a30)') "Left B.C.  =  ", Boundary_name(BC_Type(1))
        write(unit_id, '(2x,2a30)') "Right B.C.  =  ", Boundary_name(BC_Type(2))
        write(unit_id, '(2x,2a30)') "Bottom B.C.  =  ", Boundary_name(BC_Type(3))
        write(unit_id, '(2x,2a30)') "Upper B.C.  =  ", Boundary_name(BC_Type(4))
        write(unit_id, '(a)') "-------------------------------Disturbance information------------------------------"
        if (is_screen) then
            write(unit_id, '(2x,2a30)') "Iflag_unsteady  =  ", Option_name(Iflag_DIS)
        else
            write(unit_id, '(2x,2a30)') "If_DIS_simulation  =  ", Option_name(Iflag_DIS)
        endif
        write(unit_id, '(2x,2a30)') "Iflag_fringe  =  ", Option_name(Iflag_Fringe)
        write(unit_id, '(2x,2a30)') "DISBC_Type  =  ", DISBC_name(DISBC_Type)

        if (DISBC_Type == DISBC_Inflow_TS) then
            write(unit_id, '(2x,a)') "#Eigv_TS"
            do i = 1, Num_dis_ts
                write(unit_id, '(2x,a8,5f12.8)') wave_ts(i)%name, wave_ts(i)%omer, wave_ts(i)%omei, &
                                                 wave_ts(i)%alfr, wave_ts(i)%alfi, wave_ts(i)%amp
            enddo
        elseif (DISBC_Type == DISBC_Free_Wave) then
            write(unit_id, '(2x,a)') "#Eigv_free"
            do i = 1, Num_dis_free
                write(unit_id, '(2x,a8,3f12.8)') wave_free(i)%name, wave_free(i)%omer, wave_free(i)%theta, wave_free(i)%amp
            enddo
        endif

        if (Iflag_Fringe == 1) then
            write(unit_id, '(2x,a)') "#Eigv_fringe"
            write(unit_id, '(2x,I8,4f12.6)') Iters_fringe, rfringe(1:4)
        endif

        if (is_screen) then
            write(unit_id, '(a)') "--------------------------------------------------------------------------------"
        else
            write(unit_id, '(a)') "================================================================================="
        endif
        if(S2W1_wall_mode==1)then
            write(unit_id,'(A,2ES24.16)')'S2 effective wall T_native: ',S2WallTNative
            write(unit_id,'(A,2ES24.16)')'S2 effective wall theta: ',S2W1_wall_theta
            write(unit_id,'(A,2ES24.16)')'S2 effective wall lambda: ',1.0_KREAL/S2WallTNative
        endif
    end subroutine WriteParameterReport

end module ImportPara

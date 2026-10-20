module mod_R5Observe
    use mod_ClosedCouetteMass
    use mod_ImplicitDiagnostics, only: UpdateImplicitResidualDiagnostics, ImplicitConverged, GetImplicitReferences, AbortSteadyImplicit, InheritImplicitReferences
    use Boundary
    use mod_S2D1,only:IsS2CN2
    use mod_Interpolation
    use mod_Evolution
    use mod_Residual
    use mpi_wrapper
    implicit none
    real(KREAL) :: native_science(22),last_view_science(22)
    real(KREAL) :: profile_seconds=0.0_KREAL,diag_seconds=0.0_KREAL,checkpoint_seconds=0.0_KREAL,csv_seconds=0.0_KREAL
    real(KREAL) :: cost_prev(9)=0.0_KREAL,control_begin=0.0_KREAL,control_limit=230.0_KREAL
    integer :: old_forward=0,old_backward=0,origin_iteration=100,verified_streak=0
    integer :: evaluation_state=-1,last_checkpoint=-1,profile_state=100
    logical :: verify_active=.false.,references_written=.false.,last_view_pass=.false.
    character(len=64) :: r5_stop_reason='running'
    logical :: r6_preflight=.false.
    integer :: reference_epoch_state=0
    logical :: monitor_initialized=.false.,light_initialized=.false.,cost_initialized=.false.
    logical :: profile_force=.false.
    logical :: wall_output_initialized=.false.
    integer :: last_monitor=-1
    real(KREAL) :: monitor_residual(8),monitor_flux(16),monitor_totals(5),monitor_metrics(5),monitor_h,monitor_b

contains
    function SciencePack() result(v)
        real(KREAL) :: v(22)
        v=[implicit_macro_residual,implicit_macro_relative,implicit_macro_absolute, &
            implicit_h_norm,implicit_b_norm,implicit_h_relative,implicit_b_relative,implicit_h_absolute,implicit_b_absolute, &
            implicit_moment_error,implicit_mass_defect,r5_raw_h,r5_raw_b]
    end function
    subroutine RestoreScience(v)
        real(KREAL),intent(in) :: v(22)
        implicit_macro_residual=v(1:4);implicit_macro_relative=v(5:8);implicit_macro_absolute=v(9:12)
        implicit_h_norm=v(13);implicit_b_norm=v(14);implicit_h_relative=v(15);implicit_b_relative=v(16)
        implicit_h_absolute=v(17);implicit_b_absolute=v(18);implicit_moment_error=v(19);implicit_mass_defect=v(20)
        r5_raw_h=v(21);r5_raw_b=v(22)
    end subroutine
    subroutine PhysicalDiagnostic(residual,flux,rebuild)
        real(KREAL),intent(out) :: residual(8),flux(16)
        logical,intent(in) :: rebuild
        real(KREAL) :: start,profile_start,fl(4,4),fg(4,4),saved_science(22)
        real(KREAL) :: refh(Nvel),refb(Nvel),dh(Nvel),db(Nvel)
        integer :: i,j
        logical :: ok
        type(CellCenter),allocatable :: view_ctr(:,:)
        type(CellInterface),allocatable :: view_vface(:,:),view_hface(:,:)
        real(KREAL),allocatable :: view_cell_h_store(:,:,:)
        real(KREAL),allocatable :: view_cell_b_store(:,:,:)
        real(KREAL),allocatable :: view_cell_sh_store(:,:,:,:)
        real(KREAL),allocatable :: view_cell_sb_store(:,:,:,:)
        real(KREAL),allocatable :: view_v_l_flux_h_store(:,:,:)
        real(KREAL),allocatable :: view_v_l_flux_b_store(:,:,:)
        real(KREAL),allocatable :: view_v_r_flux_h_store(:,:,:)
        real(KREAL),allocatable :: view_v_r_flux_b_store(:,:,:)
        real(KREAL),allocatable :: view_h_l_flux_h_store(:,:,:)
        real(KREAL),allocatable :: view_h_l_flux_b_store(:,:,:)
        real(KREAL),allocatable :: view_h_r_flux_h_store(:,:,:)
        real(KREAL),allocatable :: view_h_r_flux_b_store(:,:,:)
        real(KREAL),allocatable :: view_res_macro_flux_store(:,:,:)
        real(KREAL),allocatable :: view_res_macro_store(:,:,:)
        real(KREAL),allocatable :: view_res_h_flux_store(:,:,:)
        real(KREAL),allocatable :: view_res_b_flux_store(:,:,:)
        real(KREAL),allocatable :: view_rate_w(:,:,:)
        real(KREAL),allocatable :: view_stream_h(:,:,:)
        real(KREAL),allocatable :: view_stream_b(:,:,:)
        start=MPI_Wtime()
        saved_science=SciencePack();profile_start=profile_seconds
        if (rebuild) then
            view_ctr=ctr;view_vface=vface;view_hface=hface
            view_cell_h_store=cell_h_store
            view_cell_b_store=cell_b_store
            view_cell_sh_store=cell_sh_store
            view_cell_sb_store=cell_sb_store
            view_v_l_flux_h_store=v_l_flux_h_store
            view_v_l_flux_b_store=v_l_flux_b_store
            view_v_r_flux_h_store=v_r_flux_h_store
            view_v_r_flux_b_store=v_r_flux_b_store
            view_h_l_flux_h_store=h_l_flux_h_store
            view_h_l_flux_b_store=h_l_flux_b_store
            view_h_r_flux_h_store=h_r_flux_h_store
            view_h_r_flux_b_store=h_r_flux_b_store
            view_res_macro_flux_store=res_macro_flux_store
            view_res_macro_store=res_macro_store
            view_res_h_flux_store=res_h_flux_store
            view_res_b_flux_store=res_b_flux_store
            view_rate_w=rate_w
            view_stream_h=stream_h
            view_stream_b=stream_b
            call MPIUpdate_CellConVars()
            call MPIUpdate_CellVDF()
            call MPIUpdate_Conner_conVars()
            call MPIUpdate_Conner_VDF()
            call Set_Ghost_PhysBoundary()
            if(.not.IsS2CN2())call Interpolation()
            call Evolution()
            call ComputeResidual()
            call AssembleSteadyRateResiduals()
        endif
        if (rebuild) call UpdateImplicitResidualDiagnostics()
        residual=[implicit_macro_absolute,r5_raw_h,r5_raw_b,implicit_h_absolute,implicit_b_absolute]
        last_view_science=SciencePack()
        fl=0.0_KREAL
        do j=1,ny;do i=1,nx
            if (i==1 .and. IP_X==0) fl(:,1)=fl(:,1)-vface(i,j)%L_flux/ctr(i,j)%dt_local
            if (i==nx .and. IP_X==NP_X-1) fl(:,2)=fl(:,2)+vface(i+1,j)%R_flux/ctr(i,j)%dt_local
            if (j==1 .and. IP_Y==0) fl(:,3)=fl(:,3)-hface(i,j)%L_flux/ctr(i,j)%dt_local
            if (j==ny .and. IP_Y==NP_Y-1) fl(:,4)=fl(:,4)+hface(i,j+1)%R_flux/ctr(i,j)%dt_local
        enddo;enddo
        call MPI_Allreduce(fl,fg,16,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        flux=reshape(fg,[16])
        if (.not.all(ieee_is_finite(residual)) .or. .not.all(ieee_is_finite(flux))) &
            call MassAbort('nonfinite original physical diagnostic')
        if(mod(profile_state,Nstep_monitor)==0.or.profile_force.or. &
            (Nstep_audit>0.and.mod(profile_state,max(1,Nstep_audit))==0))call WriteProfilesAndWalls()
        if (rebuild) then
            call RestoreScience(saved_science)
            ctr=view_ctr;vface=view_vface;hface=view_hface
            cell_h_store=view_cell_h_store
            cell_b_store=view_cell_b_store
            cell_sh_store=view_cell_sh_store
            cell_sb_store=view_cell_sb_store
            v_l_flux_h_store=view_v_l_flux_h_store
            v_l_flux_b_store=view_v_l_flux_b_store
            v_r_flux_h_store=view_v_r_flux_h_store
            v_r_flux_b_store=view_v_r_flux_b_store
            h_l_flux_h_store=view_h_l_flux_h_store
            h_l_flux_b_store=view_h_l_flux_b_store
            h_r_flux_h_store=view_h_r_flux_h_store
            h_r_flux_b_store=view_h_r_flux_b_store
            res_macro_flux_store=view_res_macro_flux_store
            res_macro_store=view_res_macro_store
            res_h_flux_store=view_res_h_flux_store
            res_b_flux_store=view_res_b_flux_store
            rate_w=view_rate_w
            stream_h=view_stream_h
            stream_b=view_stream_b
        endif
        if (rebuild) diag_seconds=diag_seconds+MPI_Wtime()-start-(profile_seconds-profile_start)
    end subroutine

    subroutine R5FiniteFluxGuard()
        logical :: local,global
        integer :: i,j
        if (.not.ClosedCouetteMassRescale) return
        local=.true.
        do j=1,ny;do i=1,nx+1
            local=local .and. all(ieee_is_finite(vface(i,j)%L_flux)) .and. all(ieee_is_finite(vface(i,j)%R_flux))
        enddo;enddo
        do j=1,ny+1;do i=1,nx
            local=local .and. all(ieee_is_finite(hface(i,j)%L_flux)) .and. all(ieee_is_finite(hface(i,j)%R_flux))
        enddo;enddo
        call MPI_Allreduce(local,global,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global) then
            r5_stop_reason='nonfinite_original_flux'
            call WriteR5Checkpoint(.true.)
            call SaveRejectedProposal()
            call AbortSteadyImplicit('R5 nonfinite flux; last accepted physical state and proposal caches saved')
        endif
    end subroutine

    subroutine StartR5Control()
        character(len=80) :: value
        integer :: status,unit,ios,k
        real(KREAL) :: h,b,rtol,mtol
        logical :: exists
        origin_iteration=restart_step
        reference_epoch_state=origin_iteration
        if (.not.ClosedCouetteMassRescale) then
            r5_force_true_updates=.false.
            return
        endif
        if (origin_iteration<0) call MassAbort('invalid checkpoint origin iteration')
        r5_force_true_updates=ClosedCouetteMassRescale .and. IsImplicitSteady()
        inquire(file='residual_reference.in',exist=exists)
        if (exists) then
            open(newunit=unit,file='residual_reference.in',status='old',iostat=ios)
            if (ios/=0) call MassAbort('cannot open residual reference inheritance file')
            read(unit,*,iostat=ios) reference_epoch_state,h,b,rtol,mtol
            close(unit)
            if (ios/=0) call MassAbort('cannot read residual reference inheritance file')
            if (rtol/=Implicit_ResidualTol .or. mtol/=Implicit_MomentTol) &
                call MassAbort('reference inheritance convergence thresholds mismatch')
            if (reference_epoch_state>origin_iteration) call MassAbort('reference epoch after checkpoint')
            call InheritImplicitReferences(h,b)
        else if (Iflag_InitFlow==2) then
            call MassAbort('Couette restart requires residual_reference.in; refusing reference reset')
        endif
        call get_environment_variable('R6_MODE',value,status=status)
        if (status==0) r6_preflight=trim(value)=='preflight'
        call get_environment_variable('R6_INITIAL_STREAK',value,status=status)
        if (status==0) then
            read(value,*) verified_streak
            if (verified_streak<0 .or. verified_streak>3) call MassAbort('invalid inherited streak')
        endif
        if (myID==0) then
            open(newunit=unit,file='velocity_loaded.dat',status='replace')
            write(unit,'(i0)') Nvel
            do k=1,Nvel
                write(unit,'(*(es25.16e3,:,","))') uVel(k),vVel(k),wVel(k)
            enddo
            close(unit)
        endif
        control_begin=MPI_Wtime()
        call get_environment_variable('R5_WALL_LIMIT',value,status=status)
        if (status==0) read(value,*) control_limit
    end subroutine

    logical function R5ShouldStop()
        logical :: request,global_request
        R5ShouldStop=.false.
        if (.not.r5_force_true_updates) return
        request=MPI_Wtime()-control_begin>=control_limit
        if (myID==0) then
            inquire(file='STOP_REQUEST',exist=global_request)
            request=request .or. global_request
        endif
        call MPI_Allreduce(request,global_request,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
        R5ShouldStop=global_request
        if (global_request) r5_stop_reason='resource_budget_controlled_stop'
    end function

    subroutine BeginR5Step()
        old_forward=r4_forward_scans;old_backward=r4_backward_scans
    end subroutine

    subroutine AfterNativeResidual()
        real(KREAL) :: h,b,start
        integer :: unit
        if (.not.ClosedCouetteMassRescale) return
        native_science=SciencePack()
        if (.not.references_written) then
            call GetImplicitReferences(h,b)
            if (myID==0) then
                open(newunit=unit,file='reference_epoch.csv',status='replace')
                write(unit,'(a)') 'origin_state,reference_h,reference_b,residual_tolerance,moment_tolerance'
                write(unit,'(i0,",")',advance='no') reference_epoch_state
                write(unit,'(*(es25.16e3,:,","))') h,b,Implicit_ResidualTol,Implicit_MomentTol
                close(unit)
            endif
            references_written=.true.
        endif
        if (ImplicitConverged()) verify_active=.true.
        ! R6 records each post-update accepted state; the initial state was checked separately.
    end subroutine

    subroutine CheckR5Policy(ok)
        logical,intent(inout) :: ok
        logical :: local,global
        local=ok .and. implicit_alpha==1.0_KREAL .and. r5_macro_alpha==1.0_KREAL .and. &
            implicit_ncfl_active==1000.0_KREAL .and. implicit_retry_count==0 .and. &
            r4_forward_scans-old_forward==1 .and. r4_backward_scans-old_backward==1
        call MPI_Allreduce(local,global,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global
    end subroutine

    subroutine RecordAcceptedLight()
        integer :: unit,k
        real(KREAL) :: start
        if(mod(origin_iteration+iter+1,Nstep_monitor)/=0)return
        start=MPI_Wtime()
        if (myID==0) then
            open(newunit=unit,file='light_steps.csv',status='unknown',position='append')
            if (.not.light_initialized) then
                write(unit,'(a)',advance='no') 'local_update,cumulative_update,native_residual_state,M0,c_minus1,Mstar_W,Mstar_h,Mnew_W,Mnew_h'
                write(unit,'(a)',advance='no') ',alpha,macro_alpha,nCFL,retries,forward_scans,backward_scans,moment_star,moment_new,min_rho,min_internal_energy,finite,original_gate,scale_gate'
                write(unit,'(a)',advance='no') ',dM,dPx,dPy,dE,ledger_M,ledger_Px,ledger_Py,ledger_E'
                do k=1,22
                    write(unit,'(a,i0)',advance='no') ',native_',k
                enddo
                write(unit,*)
            endif
            write(unit,'(3(i0,","))',advance='no') iter+1,origin_iteration+iter+1,origin_iteration+iter
            write(unit,'(*(es25.16e3,:,","))') closed_mass_target,mass_scale-1.0_KREAL, &
                mass_before(1),mass_before(5),mass_after(1),mass_after(5), &
                implicit_alpha,r5_macro_alpha,implicit_ncfl_active,real(implicit_retry_count,KREAL), &
                real(r4_forward_scans-old_forward,KREAL),real(r4_backward_scans-old_backward,KREAL), &
                mass_before_metrics(1),mass_after_metrics(1),mass_after_metrics(4:5), &
                1.0_KREAL,1.0_KREAL,1.0_KREAL,scale_increment,correction_ledger,native_science
            close(unit)
            light_initialized=.true.
        endif
        csv_seconds=csv_seconds+MPI_Wtime()-start
    end subroutine

    subroutine AfterR5Accepted()
        if (.not.ClosedCouetteMassRescale) return
        ! Exactly one distinct state follows each genuine PR(2) + one forward/backward sweep.
        call EvaluateSparseState('accepted_convergence_verification',.true.)
        if (last_view_pass) then
            verified_streak=verified_streak+1
        else
            verified_streak=0
        endif
        call WriteR5Checkpoint(.false.)
    end subroutine

    subroutine EvaluateSparseState(reason,rebuild)
        character(len=*),intent(in) :: reason
        logical,intent(in) :: rebuild
        real(KREAL) :: residual(8),flux(16),totals(5),metrics(5),start,h,b
        integer :: k,unit
        logical :: ok
        if (evaluation_state==origin_iteration+iter) return
        profile_state=origin_iteration+iter
        start=MPI_Wtime()
        call MassStatistics(totals,metrics,ok)
        if (.not.ok) call MassAbort('R5 sparse accepted-state gate failed')
        call PhysicalDiagnostic(residual,flux,rebuild)
        last_view_pass=all(last_view_science(1:4)<=Implicit_ResidualTol) .and. &
            all(last_view_science(13:14)<=Implicit_ResidualTol) .and. metrics(1)<=Implicit_MomentTol
        call GetImplicitReferences(h,b)
        monitor_residual=residual;monitor_flux=flux;monitor_totals=totals;monitor_metrics=metrics
        monitor_h=h;monitor_b=b
        call WriteMonitorRows(reason,.false.)
        evaluation_state=profile_state
        ! PhysicalDiagnostic accounts full evaluation including the temporary cache view.
        csv_seconds=csv_seconds+MPI_Wtime()-start ! Inclusive sparse dispatch (separately documented).
    end subroutine

    subroutine WriteMonitorRows(reason,force)
        character(len=*),intent(in)::reason
        logical,intent(in)::force
        integer::unit,k
        if(last_monitor==profile_state)return
        if(.not.force.and.mod(profile_state,Nstep_monitor)/=0)return
        if (myID==0) then
            open(newunit=unit,file='sparse_residual.csv',status='unknown',position='append')
            if (.not.monitor_initialized) then
                write(unit,'(a)',advance='no') 'local_state,cumulative_state,source,original_gate,moment,M_W,M_h'
                do k=1,8
                    write(unit,'(a,i0)',advance='no') ',R',k
                enddo
                do k=1,22
                    write(unit,'(a,i0)',advance='no') ',science_',k
                enddo
                do k=1,16
                    write(unit,'(a,i0)',advance='no') ',flux_',k
                enddo
                write(unit,'(a)') ',reference_h,reference_b,c_minus_1,mass_rescale_enabled'
            endif
            write(unit,'(2(i0,","),a,",")',advance='no') iter,profile_state,trim(reason)
            write(unit,'(*(es25.16e3,:,","))') merge(1.0_KREAL,0.0_KREAL,last_view_pass),monitor_metrics(1),monitor_totals(1),monitor_totals(5), &
                monitor_residual,last_view_science,monitor_flux,monitor_h,monitor_b, &
                merge(mass_scale-1.0_KREAL,0.0_KREAL,ClosedCouetteMassRescale), &
                merge(1.0_KREAL,0.0_KREAL,ClosedCouetteMassRescale)
            close(unit)
        endif
        if (myID==0) then
            open(newunit=unit,file='consecutive_states.csv',status='unknown',position='append')
            if (.not.monitor_initialized) write(unit,'(a)') 'cumulative_state,pass,streak_after_evaluation'
            write(unit,'(i0,",",l1,",",i0)') profile_state,last_view_pass, &
                merge(merge(verified_streak,verified_streak+1,r6_preflight.or.force),0,last_view_pass)
            close(unit)
        endif
        last_monitor=profile_state;monitor_initialized=.true.
    end subroutine WriteMonitorRows

    subroutine WriteProfilesAndWalls()
        real(KREAL) :: prim(4),q(2),stress,fw(4),shear,heat,begin_io
        real(KREAL) :: local(10,2),global(10,2)
        integer :: i,j,unit,wall
        character(len=96) :: name
        begin_io=MPI_Wtime()
        if(Nstep_audit>0.and.(mod(profile_state,max(1,Nstep_audit))==0.or.profile_force))then
        write(name,'(a,i6.6,a,i0,a)') 'profile_',profile_state,'_rank',myID,'.csv'
        open(newunit=unit,file=trim(name),status='replace')
        write(unit,'(a)') 'x,y,rho,U,V,T_native,qx,qy,tau,Pxy'
        do j=1,ny;do i=1,nx
            prim=GetPrimary(ctr(i,j)%conVars)
            q=GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim)
            stress=sum(wVel*(uVel-prim(2))*(vVel-prim(3))*ctr(i,j)%h_vec)
            write(unit,'(*(es25.16e3,:,","))') ctr(i,j)%x,ctr(i,j)%y,prim(1:3),1.0_KREAL/prim(4),q,GetTau(prim),stress
        enddo;enddo
        close(unit)
        endif
        if(mod(profile_state,Nstep_monitor)/=0.and..not.profile_force)return
        local=0.0_KREAL
        do i=1,nx
            if (IP_Y==0) then
                fw=-hface(i,1)%L_flux/ctr(i,1)%dt_local
                local(1:4,1)=local(1:4,1)+fw
                local(5,1)=local(5,1)+hface(i,1)%length
                local(6,1)=local(6,1)+fw(4)-WallUx(3)*fw(2)
                local(7:9,1)=local(7:9,1)+hface(i,1)%Aerodynamics*hface(i,1)%length
            endif
            if (IP_Y==NP_Y-1) then
                fw=hface(i,ny+1)%R_flux/ctr(i,ny)%dt_local
                local(1:4,2)=local(1:4,2)+fw
                local(5,2)=local(5,2)+hface(i,ny+1)%length
                local(6,2)=local(6,2)+fw(4)-WallUx(4)*fw(2)
                local(7:9,2)=local(7:9,2)+hface(i,ny+1)%Aerodynamics*hface(i,ny+1)%length
            endif
        enddo
        call MPI_Allreduce(local,global,20,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        if (myID==0) then
            open(newunit=unit,file='wall_history.csv',status='unknown',position='append')
            if (.not.wall_output_initialized) write(unit,'(a)') 'cumulative_state,wall,out_mass,out_Px,out_Py,out_E,length,out_heat,shear_per_length,heat_per_length,native_pnn,native_pnt,native_qn'
            do wall=1,2
                write(unit,'(2(i0,","))',advance='no') profile_state,wall
                write(unit,'(*(es25.16e3,:,","))') global(1:6,wall),global(2,wall)/global(5,wall), &
                    global(6,wall)/global(5,wall),global(7:9,wall)/global(5,wall)
            enddo
            close(unit)
        endif
        wall_output_initialized=.true.
        profile_seconds=profile_seconds+MPI_Wtime()-begin_io
    end subroutine

    subroutine WriteR5Checkpoint(force)
        logical,intent(in) :: force
        integer :: unit,i,j,ios,version,cumulative
        real(KREAL) :: start,mass_ref_h,mass_ref_b
        character(len=96) :: name
        if (.not.ClosedCouetteMassRescale) return
        cumulative=origin_iteration+iter
        if (Nstep_audit==0) return
        if (.not.force .and. mod(cumulative,max(1,Nstep_audit))/=0) return
        if (last_checkpoint==cumulative) return
        start=MPI_Wtime()
        write(name,'(a,i6.6,a,i0,a)') 'state_',cumulative,'_rank',myID,'.bin'
        open(newunit=unit,file=trim(name),access='stream',form='unformatted',status='replace',iostat=ios)
        if (ios/=0) call MassAbort('R5 physical checkpoint open failed')
        write(unit) nx,ny,Nvel,cumulative,myID,closed_mass_target,correction_ledger
        do j=1,ny;do i=1,nx
            write(unit) ctr(i,j)%x,ctr(i,j)%y,ctr(i,j)%area,ctr(i,j)%conVars, &
                ctr(i,j)%h_vec,ctr(i,j)%b_vec,GetTau(GetPrimary(ctr(i,j)%conVars))
        enddo;enddo
        close(unit)
        write(name,'(a,i6.6,a,i0,a)') 'restart_',cumulative,'_rank',myID,'.bin'
        open(newunit=unit,file=trim(name),access='stream',form='unformatted',status='replace',iostat=ios)
        if (ios/=0) call MassAbort('R5 complete checkpoint open failed')
        version=4
        write(unit) version,nx,ny,Nvel,LAP,cumulative,closed_mass_target,correction_ledger,RestartIdentity()
        do j=1-LAP,ny+LAP;do i=1-LAP,nx+LAP
            write(unit) ctr(i,j)%conVars,ctr(i,j)%h_vec,ctr(i,j)%b_vec
        enddo;enddo
        close(unit)
        call GetImplicitReferences(mass_ref_h,mass_ref_b)
        if (myID==0) then
            write(name,'(a,i6.6,a)') 'reference_',cumulative,'.in'
            open(newunit=unit,file=trim(name),status='replace')
            write(unit,'(i0,1x,4(es25.16e3,1x))') reference_epoch_state,mass_ref_h,mass_ref_b,Implicit_ResidualTol,Implicit_MomentTol
            close(unit)
        endif
        last_checkpoint=cumulative
        checkpoint_seconds=checkpoint_seconds+MPI_Wtime()-start
    end subroutine

    subroutine SaveRejectedProposal()
        integer :: unit,i,j
        character(len=64) :: name
        write(name,'(a,i0,a)') 'rejected_proposal_rank',myID,'.bin'
        open(newunit=unit,file=trim(name),access='stream',form='unformatted',status='replace')
        write(unit) origin_iteration+iter+1,closed_mass_target,correction_ledger,implicit_alpha,r5_macro_alpha,implicit_ncfl_active
        do j=1,ny;do i=1,nx
            write(unit) ctr(i,j)%conVars,ctr(i,j)%h_vec,ctr(i,j)%b_vec,delta_w(:,i,j),delta_h(:,i,j),delta_b(:,i,j),tau_tilde(i,j)
        enddo;enddo
        close(unit)
    end subroutine

    subroutine FinalR5State()
        integer :: unit
        if (.not.ClosedCouetteMassRescale) return
        profile_force=.true.
        call EvaluateSparseState('final_accepted_state',.true.)
        call WriteMonitorRows('final_accepted_state',.true.)
        profile_force=.false.
        call WriteR5Checkpoint(.true.)
        if (myID==0) then
            open(newunit=unit,file='R5_status.txt',status='replace')
            write(unit,'(a)') 'local_steps cumulative_state verified_streak final_original_gate stop_reason'
            write(unit,'(3(i0,1x),l1,1x,a)') iter,origin_iteration+iter,verified_streak,last_view_pass,trim(r5_stop_reason)
            close(unit)
        endif
    end subroutine

    subroutine RecordR5Cost(begin,output_seconds)
        real(KREAL),intent(in) :: begin,output_seconds
        real(KREAL) :: cumulative(9),parts(9),global(9),loop,loopglobal
        integer :: unit
        cumulative=[implicit_pr_seconds,implicit_lusgs_seconds,mass_scale_seconds,mass_refresh_seconds, &
            diag_seconds,checkpoint_seconds,profile_seconds,output_seconds,csv_seconds]
        parts=cumulative-cost_prev;cost_prev=cumulative
        call MPI_Allreduce(parts,global,9,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        loop=MPI_Wtime()-begin
        call MPI_Allreduce(loop,loopglobal,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        if (myID==0.and.mod(origin_iteration+iter,Nstep_monitor)==0) then
            open(newunit=unit,file='cost_steps.csv',position='append',status='unknown')
            if (.not.cost_initialized) write(unit,'(a)') 'local_update,cumulative_update,PR,micro_single_sweep,mass_reduce_scale,dependent_refresh,extra_residual_diag,checkpoint_io,profile_wall_io,native_io,observer_dispatch_inclusive,full_loop'
            write(unit,'(2(i0,","))',advance='no') iter,origin_iteration+iter
            write(unit,'(*(es25.16e3,:,","))') global,loopglobal
            close(unit)
            cost_initialized=.true.
        endif
    end subroutine

    subroutine RecordTerminalIO(seconds)
        real(KREAL),intent(in) :: seconds
        real(KREAL) :: global
        integer :: unit,cumulative
        character(len=80) :: name
        if(.not.ClosedCouetteMassRescale)return
        cumulative=origin_iteration+iter
        ! Only the output epoch header changes; W/h/b, M0 and ledger bytes remain untouched.
        if (ClosedCouetteMassRescale) then
        write(name,'(a,i0,a)') 'closed_mass_rank',myID,'.bin'
        open(newunit=unit,file=trim(name),access='stream',form='unformatted',status='old',action='readwrite')
        write(unit,pos=21) cumulative
        close(unit)
        endif
        call MPI_Allreduce(seconds,global,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        if (myID==0) then
            open(newunit=unit,file='terminal_io.csv',status='replace')
            write(unit,'(a)') 'terminal_native_io_seconds'
            write(unit,'(es25.16e3)') global
            close(unit)
        endif
    end subroutine
end module mod_R5Observe

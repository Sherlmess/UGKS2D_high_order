module mod_UnsteadyImplicit
    use flow_data
    use Flow_DIS
    use tools, only: GetPrimary, GetSoundSpeed, GetTau, GetHeatFlux_Vec, &
        DiscreteMaxwell_Vec, ShakhovPart_Vec
    use mod_ImplicitState
    use mod_ImplicitCommon, only: ProjectToConservedMoments
    use mod_ImplicitMacro, only: SolveMacroPR2, BuildPredictedMacroState
    use mod_ImplicitMicro, only: BuildPredictedEquilibrium, SolveMicroLUSGS, &
        RecoverMacroFromKineticState, EnforceStreamingMomentCompatibility, &
        EnforceFreestreamGhostMomentCompatibility
    use mod_SteadyImplicit, only: FindAdmissibleMicroAlpha
    use mod_Residual, only: AssembleSteadyRateResiduals, ComputeResidual
    use mod_Interpolation, only: Interpolation
    use mod_Evolution, only: Evolution
    use Boundary, only: Set_Ghost_PhysBoundary
    use mpi_wrapper, only: MPIUpdate_CellConVars, MPIUpdate_CellVDF, &
        MPIUpdate_Conner_conVars, MPIUpdate_Conner_VDF
    use mpi
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(KREAL), allocatable :: old_w(:,:,:),old_h(:,:,:),old_b(:,:,:)
    real(KREAL), allocatable :: old_macro_flux(:,:,:),old_stream_h(:,:,:),old_stream_b(:,:,:)
    real(KREAL), allocatable :: old_collision_h(:,:,:),old_collision_b(:,:,:)
    real(KREAL), allocatable :: flux_weight(:,:),old_fringe(:,:)
    real(KREAL) :: unsteady_macro_residual(4)=0.0_KREAL
    real(KREAL) :: unsteady_h_residual=0.0_KREAL,unsteady_b_residual=0.0_KREAL
    integer(KINT) :: unsteady_inner_iterations=0
contains
    pure real(KREAL) function CappedFluxTime(cfl_interval,physical_interval)
        real(KREAL), intent(in) :: cfl_interval,physical_interval
        CappedFluxTime=min(cfl_interval,physical_interval)
    end function CappedFluxTime

    pure real(KREAL) function ModifiedFluxWeight(dt_s,dt_physical,theta)
        real(KREAL), intent(in) :: dt_s,dt_physical,theta
        ModifiedFluxWeight=theta*(1.0_KREAL-dt_s/dt_physical)
    end function ModifiedFluxWeight

    pure real(KREAL) function MacroPhysicalResidual(old_state,current_state,old_flux,current_flux, &
        base_rate,old_source,current_source,dt_physical,flux_weight,theta)
        real(KREAL), intent(in) :: old_state,current_state,old_flux,current_flux,base_rate
        real(KREAL), intent(in) :: old_source,current_source,dt_physical,flux_weight,theta
        MacroPhysicalResidual=(old_state-current_state)/dt_physical+ &
            (1.0_KREAL-flux_weight)*old_flux+flux_weight*current_flux-base_rate+ &
            (1.0_KREAL-theta)*old_source+theta*current_source
    end function MacroPhysicalResidual

    pure real(KREAL) function MicroPhysicalResidual(old_state,current_state,old_stream,current_stream, &
        old_collision,current_collision,base_rate,old_source,current_source,dt_physical,flux_weight,theta)
        real(KREAL), intent(in) :: old_state,current_state,old_stream,current_stream
        real(KREAL), intent(in) :: old_collision,current_collision,base_rate
        real(KREAL), intent(in) :: old_source,current_source,dt_physical,flux_weight,theta
        MicroPhysicalResidual=(old_state-current_state)/dt_physical+ &
            (1.0_KREAL-flux_weight)*old_stream+flux_weight*current_stream+ &
            (1.0_KREAL-theta)*old_collision+theta*current_collision-base_rate+ &
            (1.0_KREAL-theta)*old_source+theta*current_source
    end function MicroPhysicalResidual

    pure real(KREAL) function MacroUnsteadyDiagonal(dt_physical,flux_weight,sum_gamma,cell_area)
        real(KREAL), intent(in) :: dt_physical,flux_weight,sum_gamma,cell_area
        MacroUnsteadyDiagonal=1.0_KREAL/dt_physical+0.5_KREAL*flux_weight*sum_gamma/cell_area
    end function MacroUnsteadyDiagonal

    pure real(KREAL) function MicroUnsteadyDiagonal(dt_physical,collision_weight,flux_weight,tau_value, &
        upwind_sum)
        real(KREAL), intent(in) :: dt_physical,collision_weight,flux_weight,tau_value,upwind_sum
        MicroUnsteadyDiagonal=1.0_KREAL/dt_physical+collision_weight/tau_value+flux_weight*upwind_sum
    end function MicroUnsteadyDiagonal

    subroutine ComputeUnsteadyFluxTimeScales()
        integer(KINT) :: i,j
        real(KREAL) :: prim(4),sos,cfl_interval

        ctr(:,:)%dt_local=dt_global
        do j=1,ny
            do i=1,nx
                prim=GetPrimary(ctr(i,j)%conVars)
                sos=GetSoundSpeed(prim)
                prim(2)=max(U_MAX,abs(prim(2))+sos)
                prim(3)=max(V_MAX,abs(prim(3))+sos)
                cfl_interval=CFL*ctr(i,j)%area/&
                    (ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3))
                ctr(i,j)%dt_local=CappedFluxTime(cfl_interval,dt_global)
            enddo
        enddo
    end subroutine ComputeUnsteadyFluxTimeScales

    subroutine AllocateUnsteadyImplicitState()
        if (.not.allocated(rate_w)) call AllocateImplicitState()
        if (allocated(old_w)) return
        allocate(old_w(4,1:nx,1:ny),old_h(Nvel,1:nx,1:ny),old_b(Nvel,1:nx,1:ny))
        allocate(old_macro_flux(4,1:nx,1:ny))
        allocate(old_stream_h(Nvel,1:nx,1:ny),old_stream_b(Nvel,1:nx,1:ny))
        allocate(old_collision_h(Nvel,1:nx,1:ny),old_collision_b(Nvel,1:nx,1:ny))
        allocate(flux_weight(1:nx,1:ny),old_fringe(1:nx,1:ny))
        old_w=0.0_KREAL; old_h=0.0_KREAL; old_b=0.0_KREAL
        old_macro_flux=0.0_KREAL; old_stream_h=0.0_KREAL; old_stream_b=0.0_KREAL
        old_collision_h=0.0_KREAL; old_collision_b=0.0_KREAL
        flux_weight=0.0_KREAL; old_fringe=0.0_KREAL
    end subroutine AllocateUnsteadyImplicitState

    subroutine FillFringeAtStep(step_index,field)
        integer(KINT), intent(in) :: step_index
        real(KREAL), intent(out) :: field(:,:)
        real(KREAL), parameter :: lambda_max=1.0_KREAL
        real(KREAL) :: fs,flength,fc,fend,fe,ramp_length,aa,bb,shape,xx
        integer(KINT) :: i

        field=0.0_KREAL
        if (Iflag_DIS/=1 .or. Iflag_Fringe/=1) return
        fend=rfringe(4); flength=rfringe(2); fc=rfringe(3)
        fs=rfringe(1)+fc*real(step_index-Iters_fringe,KREAL)
        fe=fs+flength
        if (fe>fend) then
            fe=fend; fs=fe-flength
        endif
        ramp_length=2.0_KREAL*flength/3.0_KREAL
        fe=fs+ramp_length
        do i=1,nx
            xx=ctr(i,1)%x
            if (xx<=fs) then
                field(i,:)=0.0_KREAL
            elseif (xx<=fe) then
                aa=(xx-fs)/ramp_length; bb=1.0_KREAL
                shape=0.5_KREAL*(tanh(7.0_KREAL*aa-(7.0_KREAL+bb)/2.0_KREAL)+1.0_KREAL)
                field(i,:)=lambda_max*shape
            else
                field(i,:)=lambda_max
            endif
        enddo
    end subroutine FillFringeAtStep

    subroutine ProjectBaseResidualMoments(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        real(KREAL) :: work_h(Nvel),work_b(Nvel),reference_h(Nvel),reference_b(Nvel)
        logical :: point_ok,local_ok,global_ok

        local_ok=.true.
        do j=1,ny; do i=1,nx
            work_h=ctr(i,j)%res_h0_vec; work_b=ctr(i,j)%res_b0_vec
            reference_h=ctr(i,j)%h_vec; reference_b=ctr(i,j)%b_vec
            call ProjectToConservedMoments(work_h,work_b,ctr(i,j)%res_conVars0, &
                uVel,vVel,wVel,reference_h,reference_b,point_ok)
            if (point_ok) then
                ctr(i,j)%res_h0_vec=work_h; ctr(i,j)%res_b0_vec=work_b
            endif
            local_ok=local_ok .and. point_ok
        enddo; enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    end subroutine ProjectBaseResidualMoments

    subroutine CaptureOldStateAndOperator(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        logical :: local_ok,global_ok

        call AllocateUnsteadyImplicitState()
        local_ok=.true.
        do j=1,ny; do i=1,nx
            old_w(:,i,j)=ctr(i,j)%conVars
            old_h(:,i,j)=ctr(i,j)%h_vec; old_b(:,i,j)=ctr(i,j)%b_vec
            old_macro_flux(:,i,j)=rate_w(:,i,j)
            old_stream_h(:,i,j)=stream_h(:,i,j); old_stream_b(:,i,j)=stream_b(:,i,j)
            flux_weight(i,j)=ModifiedFluxWeight(ctr(i,j)%dt_local,dt_global,Unsteady_Theta)
            local_ok=local_ok .and. flux_weight(i,j)>=0.0_KREAL .and. &
                flux_weight(i,j)<=Unsteady_Theta .and. ieee_is_finite(flux_weight(i,j))
            call CollisionForState(old_w(:,i,j),old_h(:,i,j),old_b(:,i,j), &
                old_collision_h(:,i,j),old_collision_b(:,i,j),ok)
            local_ok=local_ok .and. ok
        enddo; enddo
        old_fringe=Xlamd
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    end subroutine CaptureOldStateAndOperator

    subroutine CollisionForState(w,h,b,collision_h,collision_b,ok)
        real(KREAL), intent(in) :: w(4),h(:),b(:)
        real(KREAL), intent(out) :: collision_h(:),collision_b(:)
        logical, intent(out) :: ok
        real(KREAL) :: prim(4),tau_value,qf(2)
        real(KREAL) :: equilibrium_h_local(size(h)),equilibrium_b_local(size(b))
        real(KREAL) :: plus_h(size(h)),plus_b(size(b))
        real(KREAL) :: reference_h(size(h)),reference_b(size(b))
        logical :: projection_ok

        prim=GetPrimary(w); tau_value=GetTau(prim)
        call DiscreteMaxwell_Vec(equilibrium_h_local,equilibrium_b_local,prim)
        qf=GetHeatFlux_Vec(h,b,prim)
        call ShakhovPart_Vec(equilibrium_h_local,equilibrium_b_local,qf,prim,plus_h,plus_b)
        equilibrium_h_local=equilibrium_h_local+plus_h
        equilibrium_b_local=equilibrium_b_local+plus_b
        reference_h=equilibrium_h_local; reference_b=equilibrium_b_local
        call ProjectToConservedMoments(equilibrium_h_local,equilibrium_b_local,w, &
            uVel,vVel,wVel,reference_h,reference_b,projection_ok)
        collision_h=(equilibrium_h_local-h)/tau_value
        collision_b=(equilibrium_b_local-b)/tau_value
        ok=tau_value>0.0_KREAL .and. projection_ok .and. &
            all(ieee_is_finite(collision_h)) .and. &
            all(ieee_is_finite(collision_b))
    end subroutine CollisionForState

    subroutine AssembleUnsteadyMacroResidual()
        integer(KINT) :: i,j,m
        real(KREAL) :: old_source,current_source

        do j=1,ny; do i=1,nx; do m=1,4
            old_source=old_fringe(i,j)*(ctr(i,j)%conVars0(m)-old_w(m,i,j))
            current_source=Xlamd(i,j)*(ctr(i,j)%conVars0(m)-ctr(i,j)%conVars(m))
            rate_w(m,i,j)=MacroPhysicalResidual(old_w(m,i,j),ctr(i,j)%conVars(m), &
                old_macro_flux(m,i,j),rate_w(m,i,j),ctr(i,j)%res_conVars0(m), &
                old_source,current_source,dt_global,flux_weight(i,j),Unsteady_Theta)
        enddo; enddo; enddo
    end subroutine AssembleUnsteadyMacroResidual

    subroutine AssembleUnsteadyMicroResidual()
        integer(KINT) :: i,j,k
        real(KREAL) :: current_stream_h,current_stream_b,current_collision_h,current_collision_b
        real(KREAL) :: old_source_h,old_source_b,current_source_h,current_source_b

        do j=1,ny; do i=1,nx; do k=1,Nvel
            current_stream_h=stream_h(k,i,j); current_stream_b=stream_b(k,i,j)
            current_collision_h=(equilibrium_h(k,i,j)-ctr(i,j)%h_vec(k))/tau_tilde(i,j)
            current_collision_b=(equilibrium_b(k,i,j)-ctr(i,j)%b_vec(k))/tau_tilde(i,j)
            old_source_h=old_fringe(i,j)*(ctr(i,j)%h0_vec(k)-old_h(k,i,j))
            old_source_b=old_fringe(i,j)*(ctr(i,j)%b0_vec(k)-old_b(k,i,j))
            current_source_h=Xlamd(i,j)*(ctr(i,j)%h0_vec(k)-ctr(i,j)%h_vec(k))
            current_source_b=Xlamd(i,j)*(ctr(i,j)%b0_vec(k)-ctr(i,j)%b_vec(k))
            stream_h(k,i,j)=MicroPhysicalResidual(old_h(k,i,j),ctr(i,j)%h_vec(k), &
                old_stream_h(k,i,j),current_stream_h,old_collision_h(k,i,j),current_collision_h, &
                ctr(i,j)%res_h0_vec(k),old_source_h,current_source_h,dt_global,flux_weight(i,j),Unsteady_Theta)
            stream_b(k,i,j)=MicroPhysicalResidual(old_b(k,i,j),ctr(i,j)%b_vec(k), &
                old_stream_b(k,i,j),current_stream_b,old_collision_b(k,i,j),current_collision_b, &
                ctr(i,j)%res_b0_vec(k),old_source_b,current_source_b,dt_global,flux_weight(i,j),Unsteady_Theta)
        enddo; enddo; enddo
    end subroutine AssembleUnsteadyMicroResidual

    subroutine UpdateUnsteadyResidualNorms()
        integer(KINT) :: i,j,m,k
        real(KREAL) :: local_macro(4),global_macro(4),local_micro(2),global_micro(2)

        local_macro=0.0_KREAL; local_micro=0.0_KREAL
        do j=1,ny; do i=1,nx
            do m=1,4
                local_macro(m)=max(local_macro(m),abs(dt_global*rate_w(m,i,j))/ &
                    max(abs(old_w(m,i,j)),1.0_KREAL))
            enddo
            do k=1,Nvel
                local_micro(1)=max(local_micro(1),abs(dt_global*stream_h(k,i,j))/ &
                    max(abs(old_h(k,i,j)),1.0e-12_KREAL))
                local_micro(2)=max(local_micro(2),abs(dt_global*stream_b(k,i,j))/ &
                    max(abs(old_b(k,i,j)),1.0e-12_KREAL))
            enddo
        enddo; enddo
        call MPI_Allreduce(local_macro,global_macro,4,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(local_micro,global_micro,2,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        unsteady_macro_residual=global_macro
        unsteady_h_residual=global_micro(1); unsteady_b_residual=global_micro(2)
        res=unsteady_macro_residual
    end subroutine UpdateUnsteadyResidualNorms

    logical function UnsteadyInnerConverged()
        UnsteadyInnerConverged=max(maxval(unsteady_macro_residual), &
            unsteady_h_residual,unsteady_b_residual)<=Unsteady_InnerTol
    end function UnsteadyInnerConverged

    subroutine ApplyUnsteadyImplicitStep(accepted)
        logical, intent(out) :: accepted
        integer(KINT) :: inner,i,j
        real(KREAL) :: saved_time,macro_alpha,micro_alpha
        logical :: ok
        character(len=64) :: failure_stage

        accepted=.false.; saved_time=time_dis; failure_stage='initialization'
        call AssembleSteadyRateResiduals()
        failure_stage='old streaming moment projection'
        call EnforceStreamingMomentCompatibility(ok)
        if (.not.ok) return
        failure_stage='old state capture'
        call CaptureOldStateAndOperator(ok)
        if (.not.ok) return
        call FillFringeAtStep(iter+1,Xlamd)
        time_dis=saved_time+dt_global

        do inner=0,Unsteady_MaxInner
            failure_stage='target boundary moment projection'
            call Set_Ghost_PhysBoundary()
            call EnforceFreestreamGhostMomentCompatibility(ok)
            if (.not.ok) exit
            call Interpolation()
            call Evolution()
            call ComputeResidual()
            call AssembleSteadyRateResiduals()
            failure_stage='current streaming moment projection'
            call EnforceStreamingMomentCompatibility(ok)
            if (.not.ok) exit
            call AssembleUnsteadyMacroResidual()
            call ResetImplicitCorrections()
            failure_stage='macro PR2'
            call SolveMacroPR2(ok,dt_global,flux_weight)
            if (.not.ok) exit
            failure_stage='macro predictor'
            call BuildPredictedMacroState(macro_alpha,ok)
            if (.not.ok) exit
            failure_stage='predicted equilibrium'
            call BuildPredictedEquilibrium(ok)
            if (.not.ok) exit
            call AssembleUnsteadyMicroResidual()
            call UpdateUnsteadyResidualNorms()
            unsteady_inner_iterations=inner
            if (UnsteadyInnerConverged()) then
                accepted=.true.
                exit
            endif
            if (inner==Unsteady_MaxInner) exit
            failure_stage='micro LU-SGS'
            call SolveMicroLUSGS(ok,dt_global,flux_weight,Unsteady_Theta, &
                freeze_wall_closure=.true.)
            if (.not.ok) exit
            failure_stage='micro admissibility'
            call FindAdmissibleMicroAlpha(micro_alpha,ok)
            if (.not.ok) exit
            do j=1,ny; do i=1,nx
                ctr(i,j)%h_vec=ctr(i,j)%h_vec+micro_alpha*delta_h(:,i,j)
                ctr(i,j)%b_vec=ctr(i,j)%b_vec+micro_alpha*delta_b(:,i,j)
            enddo; enddo
            failure_stage='macro recovery'
            call RecoverMacroFromKineticState(ok)
            if (.not.ok) exit
            call MPIUpdate_CellConVars(); call MPIUpdate_CellVDF()
            call MPIUpdate_Conner_conVars(); call MPIUpdate_Conner_VDF()
        enddo

        time_dis=saved_time
        if (.not.accepted) then
            if (myID==0) write(*,'(A,A,A,I0,A,1P,6E12.4)') 'Unsteady inner failure stage=', &
                trim(failure_stage),' inner=',unsteady_inner_iterations,' residuals=', &
                unsteady_macro_residual,unsteady_h_residual,unsteady_b_residual
            do j=1,ny; do i=1,nx
                ctr(i,j)%conVars=old_w(:,i,j)
                ctr(i,j)%h_vec=old_h(:,i,j); ctr(i,j)%b_vec=old_b(:,i,j)
            enddo; enddo
        endif
    end subroutine ApplyUnsteadyImplicitStep

    subroutine WriteUnsteadyImplicitResidual()
        if (myID/=0) return
        write(*,'(A,I0,A,I0,A,1P,6E12.4)') 'Unsteady step=',iter,' inner=',unsteady_inner_iterations, &
            ' Rrho/RrhoU/RrhoV/RE/Rh/Rb=',unsteady_macro_residual,unsteady_h_residual,unsteady_b_residual
    end subroutine WriteUnsteadyImplicitResidual
end module mod_UnsteadyImplicit

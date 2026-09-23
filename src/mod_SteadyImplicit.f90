module mod_SteadyImplicit
    use flow_data
    use mod_ImplicitState
    use mod_ImplicitCommon, only: MomentConserved, IsPhysicalState
    use mod_ImplicitBoundary, only: BuildStructuredMirrorMaps, ValidateImplicitSymmetryFaces
    use mod_ImplicitMacro, only: SolveMacroPR2, BuildPredictedMacroState
    use mod_ImplicitMicro, only: BuildPredictedEquilibrium, SolveMicroLUSGS, RecoverMacroFromKineticState, &
        EnforceStreamingMomentCompatibility, EnforceAcceptedMomentCompatibility, &
        EnforceFreestreamGhostMomentCompatibility, EvaluateMicroPointBreakdown, &
        EvaluatePredictedEquilibriumBreakdown, WeightedNegativeFractionCore=>WeightedNegativeFraction
    use mod_ImplicitDiagnostics
    use mod_Residual, only: AssembleSteadyRateResiduals
    use mpi
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    public :: InitializeSteadyImplicit, ApplySteadyImplicitUpdate, FindAdmissibleMicroAlpha
    public :: ImplicitConverged, WriteImplicitResidual, AbortSteadyImplicit
    public :: EnforceImplicitBoundaryMomentCompatibility
    public :: WeightedNegativeFraction
    public :: ClassifyMicroCandidate, SelectGlobalImplicitFailure
    integer(KINT), parameter :: IMPLICIT_FAILURE_INTEGER_PAYLOAD_SIZE=12
    integer(KINT), parameter :: IMPLICIT_FAILURE_REAL_PAYLOAD_SIZE=47
contains
    subroutine InitializeSteadyImplicit(ok)
        logical, intent(out), optional :: ok
        logical :: initialized
        call AllocateImplicitState()
        call BuildStructuredMirrorMaps()
        call ValidateImplicitSymmetryFaces()
        call EnforceAcceptedMomentCompatibility(initialized)
        call InitializeImplicitDiagnostics()
        if (present(ok)) ok=initialized
        if (.not.initialized) call AbortSteadyImplicit('initial accepted-state moment projection failed')
    end subroutine InitializeSteadyImplicit

    subroutine EnforceImplicitBoundaryMomentCompatibility(ok)
        logical, intent(out), optional :: ok
        logical :: compatible

        call EnforceFreestreamGhostMomentCompatibility(compatible)
        if (present(ok)) ok=compatible
        if (.not.compatible) call AbortSteadyImplicit('freestream ghost moment projection failed')
    end subroutine EnforceImplicitBoundaryMomentCompatibility

    subroutine FindAdmissibleMicroAlpha(alpha,ok,enforce_negative_fraction)
        real(KREAL), intent(out) :: alpha
        logical, intent(out) :: ok
        logical, intent(in), optional :: enforce_negative_fraction
        real(KREAL) :: h_candidate(Nvel),b_candidate(Nvel)
        real(KREAL) :: alpha_next
        logical :: local_ok,global_ok,enforce_negative_fraction_active
        type(ImplicitFailureRecord) :: cell_failure,local_failure
        integer(KINT) :: i,j,global_i,global_j,cell_index,record_index

        enforce_negative_fraction_active=.true.
        if (present(enforce_negative_fraction)) then
            enforce_negative_fraction_active=enforce_negative_fraction
        endif
        alpha=1.0_KREAL; ok=.false.
        call ResetImplicitFailureRecord()
        do
            local_ok=.true.
            call ResetImplicitFailureRecord(local_failure)
            record_index=huge(1_KINT)
            do j=1,ny; do i=1,nx
                h_candidate=ctr(i,j)%h_vec+alpha*delta_h(:,i,j)
                b_candidate=ctr(i,j)%b_vec+alpha*delta_b(:,i,j)
                call ClassifyMicroCandidate(ctr(i,j)%h_vec,ctr(i,j)%b_vec, &
                    h_candidate,b_candidate,wVel,uVel,vVel,cell_failure, &
                    enforce_negative_fraction_active)
                if (cell_failure%valid) then
                    local_ok=.false.
                    if (allocated(Block_XY)) then
                        global_i=Block_XY(IP_XY)%Origin_I+i-1
                        global_j=Block_XY(IP_XY)%Origin_J+j-1
                    else
                        global_i=i; global_j=j
                    endif
                    cell_index=(global_j-1)*Nx_all+global_i
                    if (cell_index<record_index) then
                        local_failure=cell_failure
                        local_failure%owner_rank=myID
                        local_failure%local_i=i; local_failure%local_j=j
                        local_failure%global_i=global_i; local_failure%global_j=global_j
                        local_failure%ncfl_attempt=implicit_ncfl_active
                        local_failure%alpha_tested=alpha
                        record_index=cell_index
                    endif
                endif
            enddo; enddo
            call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
            if (global_ok) then
                call ResetImplicitFailureRecord()
                ok=.true.; return
            endif
            alpha_next=0.5_KREAL*alpha
            if (alpha_next<Implicit_MinAlpha) then
                if (local_failure%valid) local_failure%alpha_next=alpha_next
                call SelectGlobalImplicitFailure(local_failure,implicit_failure)
                if (implicit_failure%valid .and. implicit_failure%owner_rank==myID) then
                    call PopulateMicroFailureBreakdown(local_failure)
                endif
                call SelectGlobalImplicitFailure(local_failure,implicit_failure)
                alpha=alpha_next
                return
            endif
            alpha=alpha_next
        enddo
    end subroutine FindAdmissibleMicroAlpha

    subroutine PopulateMicroFailureBreakdown(record)
        type(ImplicitFailureRecord), intent(inout) :: record
        logical :: breakdown_ok,equilibrium_breakdown_ok

        record%micro_breakdown_valid=.false.
        record%equilibrium_breakdown_valid=.false.
        record%equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE
        record%equilibrium_maxwell_value=0.0_KREAL
        record%equilibrium_shakhov_correction=0.0_KREAL
        record%equilibrium_preprojection_value=0.0_KREAL
        record%equilibrium_projection_correction=0.0_KREAL
        record%equilibrium_final_value=0.0_KREAL
        record%q_frozen_x=0.0_KREAL; record%q_frozen_y=0.0_KREAL
        record%prim_tilde_rho=0.0_KREAL; record%prim_tilde_u=0.0_KREAL
        record%prim_tilde_v=0.0_KREAL; record%prim_tilde_temperature=0.0_KREAL
        record%equilibrium_maxwell_negative_fraction=0.0_KREAL
        record%equilibrium_preprojection_negative_fraction=0.0_KREAL
        record%equilibrium_final_negative_fraction=0.0_KREAL
        record%reconstruction_max_abs_h=0.0_KREAL
        record%reconstruction_max_abs_b=0.0_KREAL
        record%reconstruction_max_rel_h=0.0_KREAL
        record%reconstruction_max_rel_b=0.0_KREAL
        record%reconstruction_max_abs_h_k=0
        record%reconstruction_max_abs_b_k=0
        record%reconstructed_h_at_failure_k=0.0_KREAL
        record%stored_h_at_failure_k=0.0_KREAL
        record%reconstructed_b_at_failure_k=0.0_KREAL
        record%stored_b_at_failure_k=0.0_KREAL
        if (.not.record%valid .or. record%velocity_k<=0) return
        if (record%reason==IMPLICIT_FAILURE_H_NEGATIVE_FRACTION) then
            call EvaluateMicroPointBreakdown(record%local_i,record%local_j,record%velocity_k,.true., &
                record%micro_current,record%micro_delta,record%equilibrium_value,record%stream_rhs, &
                record%collision_rhs,record%neighbor_rhs,record%dt_pseudo,record%tau_value, &
                record%diagonal_total,record%rhs_total,record%post_sweep_delta,breakdown_ok)
        elseif (record%reason==IMPLICIT_FAILURE_B_NEGATIVE_FRACTION) then
            call EvaluateMicroPointBreakdown(record%local_i,record%local_j,record%velocity_k,.false., &
                record%micro_current,record%micro_delta,record%equilibrium_value,record%stream_rhs, &
                record%collision_rhs,record%neighbor_rhs,record%dt_pseudo,record%tau_value, &
                record%diagonal_total,record%rhs_total,record%post_sweep_delta,breakdown_ok)
        else
            return
        endif
        record%micro_candidate=record%micro_current+record%alpha_tested*record%micro_delta
        if (breakdown_ok) record%micro_breakdown_valid=.true.
        call EvaluatePredictedEquilibriumBreakdown(record%local_i,record%local_j,record%velocity_k, &
            record%reason==IMPLICIT_FAILURE_H_NEGATIVE_FRACTION, &
            record%equilibrium_maxwell_value,record%equilibrium_shakhov_correction, &
            record%equilibrium_preprojection_value,record%equilibrium_projection_correction, &
            record%equilibrium_final_value,record%q_frozen_x,record%q_frozen_y, &
            record%prim_tilde_rho,record%prim_tilde_u,record%prim_tilde_v, &
            record%prim_tilde_temperature,record%equilibrium_maxwell_negative_fraction, &
            record%equilibrium_preprojection_negative_fraction, &
            record%equilibrium_final_negative_fraction,record%reconstruction_max_abs_h, &
            record%reconstruction_max_abs_b,record%reconstruction_max_rel_h, &
            record%reconstruction_max_rel_b,record%reconstruction_max_abs_h_k, &
            record%reconstruction_max_abs_b_k,record%reconstructed_h_at_failure_k, &
            record%stored_h_at_failure_k,record%reconstructed_b_at_failure_k, &
            record%stored_b_at_failure_k,record%equilibrium_breakdown_stage, &
            equilibrium_breakdown_ok)
        if (equilibrium_breakdown_ok .and. breakdown_ok .and. &
            record%equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS) then
            record%equilibrium_breakdown_valid=record%equilibrium_final_value==record%equilibrium_value
        endif
    end subroutine PopulateMicroFailureBreakdown

    pure real(KREAL) function WeightedNegativeFraction(distribution,weights) result(fraction)
        real(KREAL), intent(in) :: distribution(:),weights(:)
        fraction=WeightedNegativeFractionCore(distribution,weights)
    end function WeightedNegativeFraction

    pure subroutine EvaluateDistributionCandidate(current,candidate,weights,current_fraction, &
        candidate_fraction,allowed_fraction,allowed)
        real(KREAL), intent(in) :: current(:),candidate(:),weights(:)
        real(KREAL), intent(out) :: current_fraction,candidate_fraction,allowed_fraction
        logical, intent(out) :: allowed
        current_fraction=WeightedNegativeFraction(current,weights)
        candidate_fraction=WeightedNegativeFraction(candidate,weights)
        allowed_fraction=max(current_fraction,1.0e-4_KREAL)
        allowed=candidate_fraction<=allowed_fraction
    end subroutine EvaluateDistributionCandidate

    subroutine ClassifyMicroCandidate(current_h,current_b,candidate_h,candidate_b,weights,u,v,record, &
        enforce_negative_fraction)
        real(KREAL), intent(in) :: current_h(:),current_b(:),candidate_h(:),candidate_b(:)
        real(KREAL), intent(in) :: weights(:),u(:),v(:)
        type(ImplicitFailureRecord), intent(out) :: record
        logical, intent(in), optional :: enforce_negative_fraction
        logical :: allowed,enforce_negative_fraction_active
        integer(KINT) :: k
        real(KREAL) :: current_fraction,candidate_fraction,allowed_fraction
        real(KREAL) :: contribution(size(weights))

        enforce_negative_fraction_active=.false.
        if (present(enforce_negative_fraction)) then
            enforce_negative_fraction_active=enforce_negative_fraction
        endif
        call ResetImplicitFailureRecord(record)
        do k=1,size(candidate_h)
            if (.not.ieee_is_finite(candidate_h(k))) then
                record%valid=.true.; record%reason=IMPLICIT_FAILURE_H_NONFINITE
                record%velocity_k=k; record%velocity_u=u(k); record%velocity_v=v(k)
                return
            endif
        enddo
        do k=1,size(candidate_b)
            if (.not.ieee_is_finite(candidate_b(k))) then
                record%valid=.true.; record%reason=IMPLICIT_FAILURE_B_NONFINITE
                record%velocity_k=k; record%velocity_u=u(k); record%velocity_v=v(k)
                return
            endif
        enddo

        call MomentConserved(candidate_h,candidate_b,u,v,weights,record%candidate_moment)
        if (all(ieee_is_finite(record%candidate_moment)) .and. &
            record%candidate_moment(1)>0.0_KREAL) then
            record%candidate_internal_energy=record%candidate_moment(4)-0.5_KREAL* &
                (record%candidate_moment(2)**2+record%candidate_moment(3)**2)/ &
                record%candidate_moment(1)
        endif

        if (enforce_negative_fraction_active) then
            call EvaluateDistributionCandidate(current_h,candidate_h,weights,current_fraction, &
                candidate_fraction,allowed_fraction,allowed)
            if (.not.allowed) then
                contribution=weights*max(-candidate_h,0.0_KREAL)
                k=maxloc(contribution,dim=1)
                record%valid=.true.; record%reason=IMPLICIT_FAILURE_H_NEGATIVE_FRACTION
                record%velocity_k=k; record%velocity_u=u(k); record%velocity_v=v(k)
                record%current_fraction=current_fraction
                record%candidate_fraction=candidate_fraction
                record%allowed_fraction=allowed_fraction
                return
            endif

            call EvaluateDistributionCandidate(current_b,candidate_b,weights,current_fraction, &
                candidate_fraction,allowed_fraction,allowed)
            if (.not.allowed) then
                contribution=weights*max(-candidate_b,0.0_KREAL)
                k=maxloc(contribution,dim=1)
                record%valid=.true.; record%reason=IMPLICIT_FAILURE_B_NEGATIVE_FRACTION
                record%velocity_k=k; record%velocity_u=u(k); record%velocity_v=v(k)
                record%current_fraction=current_fraction
                record%candidate_fraction=candidate_fraction
                record%allowed_fraction=allowed_fraction
                return
            endif
        endif

        if (.not.IsPhysicalState(record%candidate_moment)) then
            record%valid=.true.; record%reason=IMPLICIT_FAILURE_MOMENT_NONPHYSICAL
        endif
    end subroutine ClassifyMicroCandidate

    subroutine SelectGlobalImplicitFailure(local_record,global_record)
        type(ImplicitFailureRecord), intent(in) :: local_record
        type(ImplicitFailureRecord), intent(out) :: global_record
        integer(KINT) :: local_index,global_index,local_owner,owner
        integer(KINT) :: integer_payload(IMPLICIT_FAILURE_INTEGER_PAYLOAD_SIZE)
        real(KREAL) :: real_payload(IMPLICIT_FAILURE_REAL_PAYLOAD_SIZE)

        local_index=huge(1_KINT)
        if (local_record%valid) local_index=(local_record%global_j-1)*Nx_all+local_record%global_i
        call MPI_Allreduce(local_index,global_index,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

        local_owner=huge(1_KINT)
        if (local_record%valid .and. local_index==global_index) local_owner=myID
        call MPI_Allreduce(local_owner,owner,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

        call ResetImplicitFailureRecord(global_record)
        if (global_index==huge(1_KINT) .or. owner==huge(1_KINT)) return

        integer_payload=0; real_payload=0.0_KREAL
        if (myID==owner) then
            integer_payload=[local_record%reason,owner,local_record%local_i, &
                local_record%local_j,local_record%global_i,local_record%global_j, &
                local_record%velocity_k,merge(1_KINT,0_KINT,local_record%micro_breakdown_valid), &
                merge(1_KINT,0_KINT,local_record%equilibrium_breakdown_valid), &
                local_record%equilibrium_breakdown_stage, &
                local_record%reconstruction_max_abs_h_k,local_record%reconstruction_max_abs_b_k]
            real_payload=[local_record%velocity_u,local_record%velocity_v, &
                local_record%ncfl_attempt,local_record%alpha_tested,local_record%alpha_next, &
                local_record%current_fraction,local_record%candidate_fraction, &
                local_record%allowed_fraction,local_record%candidate_moment, &
                local_record%candidate_internal_energy,local_record%micro_current, &
                local_record%micro_delta,local_record%micro_candidate,local_record%equilibrium_value, &
                local_record%stream_rhs,local_record%collision_rhs,local_record%neighbor_rhs, &
                local_record%dt_pseudo,local_record%tau_value,local_record%diagonal_total, &
                local_record%rhs_total,local_record%post_sweep_delta, &
                local_record%equilibrium_maxwell_value,local_record%equilibrium_shakhov_correction, &
                local_record%equilibrium_preprojection_value, &
                local_record%equilibrium_projection_correction,local_record%equilibrium_final_value, &
                local_record%q_frozen_x,local_record%q_frozen_y, &
                local_record%prim_tilde_rho,local_record%prim_tilde_u,local_record%prim_tilde_v, &
                local_record%prim_tilde_temperature, &
                local_record%equilibrium_maxwell_negative_fraction, &
                local_record%equilibrium_preprojection_negative_fraction, &
                local_record%equilibrium_final_negative_fraction, &
                local_record%reconstruction_max_abs_h,local_record%reconstruction_max_abs_b, &
                local_record%reconstruction_max_rel_h,local_record%reconstruction_max_rel_b, &
                local_record%reconstructed_h_at_failure_k,local_record%stored_h_at_failure_k, &
                local_record%reconstructed_b_at_failure_k,local_record%stored_b_at_failure_k]
        endif
        call MPI_Bcast(integer_payload,IMPLICIT_FAILURE_INTEGER_PAYLOAD_SIZE,MPI_INTEGER,owner, &
            MPI_COMM_WORLD,ierr)
        call MPI_Bcast(real_payload,IMPLICIT_FAILURE_REAL_PAYLOAD_SIZE,MPI_KREAL,owner, &
            MPI_COMM_WORLD,ierr)

        global_record%valid=.true.
        global_record%reason=integer_payload(1)
        global_record%owner_rank=integer_payload(2)
        global_record%local_i=integer_payload(3); global_record%local_j=integer_payload(4)
        global_record%global_i=integer_payload(5); global_record%global_j=integer_payload(6)
        global_record%velocity_k=integer_payload(7)
        global_record%micro_breakdown_valid=integer_payload(8)/=0
        global_record%equilibrium_breakdown_valid=integer_payload(9)/=0
        global_record%equilibrium_breakdown_stage=integer_payload(10)
        global_record%reconstruction_max_abs_h_k=integer_payload(11)
        global_record%reconstruction_max_abs_b_k=integer_payload(12)
        global_record%velocity_u=real_payload(1); global_record%velocity_v=real_payload(2)
        global_record%ncfl_attempt=real_payload(3)
        global_record%alpha_tested=real_payload(4); global_record%alpha_next=real_payload(5)
        global_record%current_fraction=real_payload(6)
        global_record%candidate_fraction=real_payload(7)
        global_record%allowed_fraction=real_payload(8)
        global_record%candidate_moment=real_payload(9:12)
        global_record%candidate_internal_energy=real_payload(13)
        global_record%micro_current=real_payload(14)
        global_record%micro_delta=real_payload(15)
        global_record%micro_candidate=real_payload(16)
        global_record%equilibrium_value=real_payload(17)
        global_record%stream_rhs=real_payload(18)
        global_record%collision_rhs=real_payload(19)
        global_record%neighbor_rhs=real_payload(20)
        global_record%dt_pseudo=real_payload(21)
        global_record%tau_value=real_payload(22)
        global_record%diagonal_total=real_payload(23)
        global_record%rhs_total=real_payload(24)
        global_record%post_sweep_delta=real_payload(25)
        global_record%equilibrium_maxwell_value=real_payload(26)
        global_record%equilibrium_shakhov_correction=real_payload(27)
        global_record%equilibrium_preprojection_value=real_payload(28)
        global_record%equilibrium_projection_correction=real_payload(29)
        global_record%equilibrium_final_value=real_payload(30)
        global_record%q_frozen_x=real_payload(31); global_record%q_frozen_y=real_payload(32)
        global_record%prim_tilde_rho=real_payload(33)
        global_record%prim_tilde_u=real_payload(34); global_record%prim_tilde_v=real_payload(35)
        global_record%prim_tilde_temperature=real_payload(36)
        global_record%equilibrium_maxwell_negative_fraction=real_payload(37)
        global_record%equilibrium_preprojection_negative_fraction=real_payload(38)
        global_record%equilibrium_final_negative_fraction=real_payload(39)
        global_record%reconstruction_max_abs_h=real_payload(40)
        global_record%reconstruction_max_abs_b=real_payload(41)
        global_record%reconstruction_max_rel_h=real_payload(42)
        global_record%reconstruction_max_rel_b=real_payload(43)
        global_record%reconstructed_h_at_failure_k=real_payload(44)
        global_record%stored_h_at_failure_k=real_payload(45)
        global_record%reconstructed_b_at_failure_k=real_payload(46)
        global_record%stored_b_at_failure_k=real_payload(47)
    end subroutine SelectGlobalImplicitFailure

    subroutine ApplySteadyImplicitUpdate(accepted)
        logical, intent(out) :: accepted
        logical :: ok
        integer(KINT) :: retry,i,j
        real(KREAL) :: macro_alpha,start_total,start_stage,moment_local,moment_global,moment_value(4)

        accepted=.false.; start_total=MPI_Wtime()
        call AssembleSteadyRateResiduals()
        call EnforceStreamingMomentCompatibility(ok)
        if (.not.ok) then
            implicit_total_seconds=implicit_total_seconds+MPI_Wtime()-start_total
            return
        endif
        call UpdateImplicitResidualDiagnostics()
        if (ImplicitConverged()) then
            call ResetImplicitCorrections()
            accepted=.true.
            implicit_total_seconds=implicit_total_seconds+MPI_Wtime()-start_total
            return
        endif
        implicit_ncfl_active=Implicit_nCFL
        do retry=0,Implicit_MaxRetries
            call ResetImplicitCorrections()
            implicit_retry_count=retry
            start_stage=MPI_Wtime(); call SolveMacroPR2(ok)
            implicit_pr_seconds=implicit_pr_seconds+MPI_Wtime()-start_stage
            if (ok) call BuildPredictedMacroState(macro_alpha,ok)
            if (ok) call BuildPredictedEquilibrium(ok)
            if (ok) then
                start_stage=MPI_Wtime(); call SolveMicroLUSGS(ok,freeze_wall_closure=.true.)
                implicit_lusgs_seconds=implicit_lusgs_seconds+MPI_Wtime()-start_stage
            endif
            if (ok) call FindAdmissibleMicroAlpha(implicit_alpha,ok,.false.)
            if (ok) exit
            implicit_ncfl_active=0.5_KREAL*implicit_ncfl_active
        enddo
        implicit_retry_total=implicit_retry_total+implicit_retry_count
        if (.not.ok) then
            implicit_total_seconds=implicit_total_seconds+MPI_Wtime()-start_total
            return
        endif
        do j=1,ny; do i=1,nx
            ctr(i,j)%h_vec=ctr(i,j)%h_vec+implicit_alpha*delta_h(:,i,j)
            ctr(i,j)%b_vec=ctr(i,j)%b_vec+implicit_alpha*delta_b(:,i,j)
        enddo; enddo
        call RecoverMacroFromKineticState(ok)
        if (.not.ok) then
            implicit_total_seconds=implicit_total_seconds+MPI_Wtime()-start_total
            return
        endif
        moment_local=0.0_KREAL
        do j=1,ny; do i=1,nx
            call MomentConserved(ctr(i,j)%h_vec,ctr(i,j)%b_vec,uVel,vVel,wVel,moment_value)
            moment_local=max(moment_local,maxval(abs(ctr(i,j)%conVars-moment_value)/ &
                max(abs(ctr(i,j)%conVars),1.0_KREAL)))
        enddo; enddo
        call MPI_Allreduce(moment_local,moment_global,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        implicit_moment_error=moment_global
        ok=implicit_moment_error<=Implicit_MomentTol
        accepted=ok
        implicit_total_seconds=implicit_total_seconds+MPI_Wtime()-start_total
    end subroutine ApplySteadyImplicitUpdate
end module mod_SteadyImplicit

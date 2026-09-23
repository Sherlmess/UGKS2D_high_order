module mod_ImplicitMicro
    use flow_data
    use Tools, only: GetPrimary, GetTau, GetHeatFlux_Vec, DiscreteMaxwell_Vec, ShakhovPart_Vec
    use mod_ImplicitState
    use mod_ImplicitCommon, only: PseudoTimeStep, MomentConserved, IsPhysicalState, ProjectToConservedMoments
    use mod_ImplicitBoundary, only: GetSymmetryMicroValue, ComputeDiffuseWallFaceDelta, &
        ComputeMaxwellWallFaceDelta, ResolveImplicitBoundaryType, IsAxisAlignedSymmetryNormal
    use mpi_wrapper, only: MPIExchangeDeltaMicro
    use mpi
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    !Zhu-style wall closure frozen from the cell correction at one sweep stage.
    type, public :: FrozenMicroWallState
        logical :: active(4)=.false.
        logical :: fallback(4)=.false.
        integer(KINT) :: build_calls=0
        real(KREAL), allocatable :: delta_h_face(:,:),delta_b_face(:,:)
    end type FrozenMicroWallState
    public :: BuildPredictedEquilibrium, SolveMicroLUSGS, ComputeMicroDefectNorm
    public :: EvaluateMicroHDefect, EvaluateMicroBDefect, RecoverMacroFromKineticState
    public :: SolveCollisionOnlyIncrement, EnforceStreamingMomentCompatibility
    public :: EnforceAcceptedMomentCompatibility
    public :: EnforceFreestreamGhostMomentCompatibility
    public :: GetMicroFaceValue
    public :: EvaluateMicroPointBreakdown
    public :: EvaluatePredictedEquilibriumBreakdown, WeightedNegativeFraction
    public :: PrepareFrozenMicroWallState
contains
    pure real(KREAL) function WeightedNegativeFraction(distribution,weights) result(fraction)
        real(KREAL), intent(in) :: distribution(:),weights(:)

        fraction=sum(weights*max(-distribution,0.0_KREAL))/ &
            max(sum(weights*abs(distribution)),tiny(1.0_KREAL))
    end function WeightedNegativeFraction

    subroutine SolveCollisionOnlyIncrement(dt_pseudo,tau_value,current,equilibrium,delta_value)
        real(KREAL), intent(in) :: dt_pseudo,tau_value,current,equilibrium
        real(KREAL), intent(out) :: delta_value
        delta_value=((equilibrium-current)/tau_value)/(1.0_KREAL/dt_pseudo+1.0_KREAL/tau_value)
    end subroutine SolveCollisionOnlyIncrement

    subroutine BuildPredictedEquilibrium(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        real(KREAL) :: prim_old(4),prim_tilde(4),h_plus(Nvel),b_plus(Nvel)
        real(KREAL) :: reference_h(Nvel),reference_b(Nvel)
        real(KREAL) :: old_equilibrium_h(Nvel),old_equilibrium_b(Nvel)
        real(KREAL) :: old_reference_h(Nvel),old_reference_b(Nvel)
        logical :: projection_ok
        logical :: local_ok,global_ok

        local_ok=.true.
        do j=1,ny
            do i=1,nx
                prim_old=GetPrimary(ctr(i,j)%conVars)
                q_frozen(:,i,j)=GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim_old)
                prim_tilde=GetPrimary(w_tilde(:,i,j))
                if(s2_steady_rhs)then
                    tau_tilde(i,j)=0.1_KREAL
                else
                    tau_tilde(i,j)=GetTau(prim_tilde)
                endif
                call DiscreteMaxwell_Vec(equilibrium_h(:,i,j),equilibrium_b(:,i,j),prim_tilde)
                call ShakhovPart_Vec(equilibrium_h(:,i,j),equilibrium_b(:,i,j),q_frozen(:,i,j), &
                                     prim_tilde,h_plus,b_plus)
                equilibrium_h(:,i,j)=equilibrium_h(:,i,j)+h_plus
                equilibrium_b(:,i,j)=equilibrium_b(:,i,j)+b_plus
                reference_h=equilibrium_h(:,i,j); reference_b=equilibrium_b(:,i,j)
                call ProjectToConservedMoments(equilibrium_h(:,i,j),equilibrium_b(:,i,j), &
                    w_tilde(:,i,j),uVel,vVel,wVel,reference_h,reference_b,projection_ok)
                if(s2_steady_rhs.and.projection_ok)then
                    call DiscreteMaxwell_Vec(old_equilibrium_h,old_equilibrium_b,prim_old)
                    call ShakhovPart_Vec(old_equilibrium_h,old_equilibrium_b,q_frozen(:,i,j), &
                        prim_old,h_plus,b_plus)
                    old_equilibrium_h=old_equilibrium_h+h_plus
                    old_equilibrium_b=old_equilibrium_b+b_plus
                    old_reference_h=old_equilibrium_h;old_reference_b=old_equilibrium_b
                    call ProjectToConservedMoments(old_equilibrium_h,old_equilibrium_b, &
                        ctr(i,j)%conVars,uVel,vVel,wVel,old_reference_h,old_reference_b,projection_ok)
                    if(projection_ok)then
                        equilibrium_increment_h(:,i,j)=equilibrium_h(:,i,j)-old_equilibrium_h
                        equilibrium_increment_b(:,i,j)=equilibrium_b(:,i,j)-old_equilibrium_b
                    endif
                endif
                local_ok=local_ok .and. tau_tilde(i,j)>0.0_KREAL .and. &
                    projection_ok .and. &
                    all(ieee_is_finite(equilibrium_h(:,i,j))) .and. &
                    all(ieee_is_finite(equilibrium_b(:,i,j)))
            enddo
        enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    end subroutine BuildPredictedEquilibrium

    subroutine EvaluatePredictedEquilibriumBreakdown(i,j,k,is_h, &
        equilibrium_maxwell_value,equilibrium_shakhov_correction, &
        equilibrium_preprojection_value,equilibrium_projection_correction, &
        equilibrium_final_value,q_frozen_x,q_frozen_y, &
        prim_tilde_rho,prim_tilde_u,prim_tilde_v,prim_tilde_temperature, &
        equilibrium_maxwell_negative_fraction, &
        equilibrium_preprojection_negative_fraction, &
        equilibrium_final_negative_fraction,reconstruction_max_abs_h, &
        reconstruction_max_abs_b,reconstruction_max_rel_h,reconstruction_max_rel_b, &
        reconstruction_max_abs_h_k,reconstruction_max_abs_b_k, &
        reconstructed_h_at_failure_k,stored_h_at_failure_k, &
        reconstructed_b_at_failure_k,stored_b_at_failure_k, &
        equilibrium_breakdown_stage,ok)
        integer(KINT), intent(in) :: i,j,k
        logical, intent(in) :: is_h
        real(KREAL), intent(out) :: equilibrium_maxwell_value
        real(KREAL), intent(out) :: equilibrium_shakhov_correction
        real(KREAL), intent(out) :: equilibrium_preprojection_value
        real(KREAL), intent(out) :: equilibrium_projection_correction
        real(KREAL), intent(out) :: equilibrium_final_value
        real(KREAL), intent(out) :: q_frozen_x,q_frozen_y
        real(KREAL), intent(out) :: prim_tilde_rho,prim_tilde_u,prim_tilde_v
        real(KREAL), intent(out) :: prim_tilde_temperature
        real(KREAL), intent(out) :: equilibrium_maxwell_negative_fraction
        real(KREAL), intent(out) :: equilibrium_preprojection_negative_fraction
        real(KREAL), intent(out) :: equilibrium_final_negative_fraction
        real(KREAL), intent(out) :: reconstruction_max_abs_h,reconstruction_max_abs_b
        real(KREAL), intent(out) :: reconstruction_max_rel_h,reconstruction_max_rel_b
        integer(KINT), intent(out) :: reconstruction_max_abs_h_k,reconstruction_max_abs_b_k
        real(KREAL), intent(out) :: reconstructed_h_at_failure_k,stored_h_at_failure_k
        real(KREAL), intent(out) :: reconstructed_b_at_failure_k,stored_b_at_failure_k
        integer(KINT), intent(out) :: equilibrium_breakdown_stage
        logical, intent(out) :: ok
        real(KREAL) :: prim_tilde(4)
        real(KREAL), allocatable :: maxwell_h(:),maxwell_b(:)
        real(KREAL), allocatable :: shakhov_h(:),shakhov_b(:)
        real(KREAL), allocatable :: pre_h(:),pre_b(:)
        real(KREAL), allocatable :: projected_h(:),projected_b(:)
        real(KREAL), allocatable :: reference_h(:),reference_b(:)
        logical :: projection_ok,exact_reconstruction

        equilibrium_maxwell_value=0.0_KREAL
        equilibrium_shakhov_correction=0.0_KREAL
        equilibrium_preprojection_value=0.0_KREAL
        equilibrium_projection_correction=0.0_KREAL
        equilibrium_final_value=0.0_KREAL
        q_frozen_x=0.0_KREAL; q_frozen_y=0.0_KREAL
        prim_tilde_rho=0.0_KREAL; prim_tilde_u=0.0_KREAL
        prim_tilde_v=0.0_KREAL; prim_tilde_temperature=0.0_KREAL
        equilibrium_maxwell_negative_fraction=0.0_KREAL
        equilibrium_preprojection_negative_fraction=0.0_KREAL
        equilibrium_final_negative_fraction=0.0_KREAL
        reconstruction_max_abs_h=0.0_KREAL; reconstruction_max_abs_b=0.0_KREAL
        reconstruction_max_rel_h=0.0_KREAL; reconstruction_max_rel_b=0.0_KREAL
        reconstruction_max_abs_h_k=0; reconstruction_max_abs_b_k=0
        reconstructed_h_at_failure_k=0.0_KREAL; stored_h_at_failure_k=0.0_KREAL
        reconstructed_b_at_failure_k=0.0_KREAL; stored_b_at_failure_k=0.0_KREAL
        equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE
        ok=.false.

        if (i<1 .or. i>nx .or. j<1 .or. j>ny .or. k<1 .or. k>Nvel) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_INDEX
            return
        endif
        if (Nvel<=0 .or. nx<=0 .or. ny<=0) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_UNAVAILABLE_STATE
            return
        endif
        if (.not.allocated(w_tilde) .or. .not.allocated(q_frozen) .or. &
            .not.allocated(equilibrium_h) .or. .not.allocated(equilibrium_b) .or. &
            .not.allocated(uVel) .or. .not.allocated(vVel) .or. .not.allocated(wVel)) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_UNAVAILABLE_STATE
            return
        endif
        if (size(w_tilde,1)/=4 .or. size(w_tilde,2)/=nx .or. size(w_tilde,3)/=ny) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_SHAPE
            return
        endif
        if (size(q_frozen,1)/=2 .or. size(q_frozen,2)/=nx .or. size(q_frozen,3)/=ny) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_SHAPE
            return
        endif
        if (size(equilibrium_h,1)/=Nvel .or. size(equilibrium_h,2)/=nx .or. &
            size(equilibrium_h,3)/=ny .or. size(equilibrium_b,1)/=Nvel .or. &
            size(equilibrium_b,2)/=nx .or. size(equilibrium_b,3)/=ny) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_SHAPE
            return
        endif
        if (size(uVel)/=Nvel .or. size(vVel)/=Nvel .or. size(wVel)/=Nvel) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_SHAPE
            return
        endif
        if (.not.all(ieee_is_finite(w_tilde(:,i,j)))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_INPUT
            return
        endif
        if (.not.IsPhysicalState(w_tilde(:,i,j))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONPHYSICAL_W_TILDE
            return
        endif
        if (.not.all(ieee_is_finite(q_frozen(:,i,j))) .or. &
            .not.all(ieee_is_finite(uVel)) .or. .not.all(ieee_is_finite(vVel)) .or. &
            .not.all(ieee_is_finite(wVel))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_INPUT
            return
        endif
        if (.not.all(ieee_is_finite(equilibrium_h(:,i,j))) .or. &
            .not.all(ieee_is_finite(equilibrium_b(:,i,j)))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_INPUT
            return
        endif

        prim_tilde=GetPrimary(w_tilde(:,i,j))
        if (.not.all(ieee_is_finite(prim_tilde))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_INPUT
            return
        endif
        allocate(maxwell_h(Nvel),maxwell_b(Nvel),shakhov_h(Nvel),shakhov_b(Nvel))
        allocate(pre_h(Nvel),pre_b(Nvel),projected_h(Nvel),projected_b(Nvel))
        allocate(reference_h(Nvel),reference_b(Nvel))
        call DiscreteMaxwell_Vec(maxwell_h,maxwell_b,prim_tilde)
        call ShakhovPart_Vec(maxwell_h,maxwell_b,q_frozen(:,i,j),prim_tilde,shakhov_h,shakhov_b)
        if (.not.all(ieee_is_finite(maxwell_h)) .or. .not.all(ieee_is_finite(maxwell_b)) .or. &
            .not.all(ieee_is_finite(shakhov_h)) .or. .not.all(ieee_is_finite(shakhov_b))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_MAXWELL_OR_SHAKHOV
            return
        endif
        pre_h=maxwell_h+shakhov_h
        pre_b=maxwell_b+shakhov_b
        if (.not.all(ieee_is_finite(pre_h)) .or. .not.all(ieee_is_finite(pre_b))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_MAXWELL_OR_SHAKHOV
            return
        endif
        reference_h=pre_h; reference_b=pre_b
        projected_h=pre_h; projected_b=pre_b
        call ProjectToConservedMoments(projected_h,projected_b,w_tilde(:,i,j), &
            uVel,vVel,wVel,reference_h,reference_b,projection_ok)
        if (.not.projection_ok .or. .not.all(ieee_is_finite(projected_h)) .or. &
            .not.all(ieee_is_finite(projected_b))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_PROJECTION_FAILED
            return
        endif

        reconstruction_max_abs_h=maxval(abs(projected_h-equilibrium_h(:,i,j)))
        reconstruction_max_abs_b=maxval(abs(projected_b-equilibrium_b(:,i,j)))
        reconstruction_max_rel_h=maxval(abs(projected_h-equilibrium_h(:,i,j))/ &
            max(abs(equilibrium_h(:,i,j)),tiny(1.0_KREAL)))
        reconstruction_max_rel_b=maxval(abs(projected_b-equilibrium_b(:,i,j))/ &
            max(abs(equilibrium_b(:,i,j)),tiny(1.0_KREAL)))
        reconstruction_max_abs_h_k=maxloc(abs(projected_h-equilibrium_h(:,i,j)),dim=1)
        reconstruction_max_abs_b_k=maxloc(abs(projected_b-equilibrium_b(:,i,j)),dim=1)

        q_frozen_x=q_frozen(1,i,j); q_frozen_y=q_frozen(2,i,j)
        prim_tilde_rho=prim_tilde(1); prim_tilde_u=prim_tilde(2)
        prim_tilde_v=prim_tilde(3)
        prim_tilde_temperature=1.0_KREAL/prim_tilde(4)
        reconstructed_h_at_failure_k=projected_h(k)
        stored_h_at_failure_k=equilibrium_h(k,i,j)
        reconstructed_b_at_failure_k=projected_b(k)
        stored_b_at_failure_k=equilibrium_b(k,i,j)
        if (is_h) then
            equilibrium_maxwell_value=maxwell_h(k)
            equilibrium_shakhov_correction=shakhov_h(k)
            equilibrium_preprojection_value=pre_h(k)
            equilibrium_projection_correction=projected_h(k)-pre_h(k)
            equilibrium_final_value=projected_h(k)
            equilibrium_maxwell_negative_fraction=WeightedNegativeFraction(maxwell_h,wVel)
            equilibrium_preprojection_negative_fraction=WeightedNegativeFraction(pre_h,wVel)
            equilibrium_final_negative_fraction=WeightedNegativeFraction(projected_h,wVel)
        else
            equilibrium_maxwell_value=maxwell_b(k)
            equilibrium_shakhov_correction=shakhov_b(k)
            equilibrium_preprojection_value=pre_b(k)
            equilibrium_projection_correction=projected_b(k)-pre_b(k)
            equilibrium_final_value=projected_b(k)
            equilibrium_maxwell_negative_fraction=WeightedNegativeFraction(maxwell_b,wVel)
            equilibrium_preprojection_negative_fraction=WeightedNegativeFraction(pre_b,wVel)
            equilibrium_final_negative_fraction=WeightedNegativeFraction(projected_b,wVel)
        endif
        if (.not.all(ieee_is_finite([equilibrium_maxwell_value,equilibrium_shakhov_correction, &
            equilibrium_preprojection_value,equilibrium_projection_correction, &
            equilibrium_final_value,q_frozen_x,q_frozen_y,prim_tilde_rho,prim_tilde_u, &
            prim_tilde_v,prim_tilde_temperature,equilibrium_maxwell_negative_fraction, &
            equilibrium_preprojection_negative_fraction,equilibrium_final_negative_fraction]))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_MAXWELL_OR_SHAKHOV
            return
        endif
        if (.not.all(ieee_is_finite([reconstruction_max_abs_h,reconstruction_max_abs_b, &
            reconstruction_max_rel_h,reconstruction_max_rel_b,reconstructed_h_at_failure_k, &
            stored_h_at_failure_k,reconstructed_b_at_failure_k,stored_b_at_failure_k]))) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_PROJECTION_FAILED
            return
        endif
        exact_reconstruction=all(projected_h==equilibrium_h(:,i,j)) .and. &
            all(projected_b==equilibrium_b(:,i,j))
        if (exact_reconstruction) then
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS
            ok=.true.
        else
            equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH
        endif
    end subroutine EvaluatePredictedEquilibriumBreakdown

    subroutine EnforceStreamingMomentCompatibility(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        logical :: projection_ok,local_ok,global_ok

        local_ok=.true.
        do j=1,ny; do i=1,nx
            call ProjectToConservedMoments(stream_h(:,i,j),stream_b(:,i,j),rate_w(:,i,j), &
                uVel,vVel,wVel,ctr(i,j)%h_vec,ctr(i,j)%b_vec,projection_ok)
            local_ok=local_ok .and. projection_ok
        enddo; enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    end subroutine EnforceStreamingMomentCompatibility

    subroutine EnforceAcceptedMomentCompatibility(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        real(KREAL) :: reference_h(Nvel),reference_b(Nvel)
        logical :: projection_ok,local_ok,global_ok

        local_ok=.true.
        do j=1,ny; do i=1,nx
            reference_h=ctr(i,j)%h_vec; reference_b=ctr(i,j)%b_vec
            call ProjectToConservedMoments(ctr(i,j)%h_vec,ctr(i,j)%b_vec,ctr(i,j)%conVars, &
                uVel,vVel,wVel,reference_h,reference_b,projection_ok)
            local_ok=local_ok .and. projection_ok
        enddo; enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    end subroutine EnforceAcceptedMomentCompatibility

    subroutine EnforceFreestreamGhostMomentCompatibility(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        logical :: projection_ok,local_ok,global_ok

        local_ok=.true.
        if (IP_X==0 .and. BC_Type(1)==BCconst_freestream) then
            do j=1,ny; do i=1-GHOST,0
                call ProjectGhostCell(i,j,projection_ok); local_ok=local_ok .and. projection_ok
            enddo; enddo
        endif
        if (IP_X==NP_X-1 .and. BC_Type(2)==BCconst_freestream) then
            do j=1,ny; do i=nx+1,nx+GHOST
                call ProjectGhostCell(i,j,projection_ok); local_ok=local_ok .and. projection_ok
            enddo; enddo
        endif
        if (IP_Y==0 .and. BC_Type(3)==BCconst_freestream) then
            do j=1-GHOST,0; do i=1,nx
                call ProjectGhostCell(i,j,projection_ok); local_ok=local_ok .and. projection_ok
            enddo; enddo
        endif
        if (IP_Y==NP_Y-1 .and. BC_Type(4)==BCconst_freestream) then
            do j=ny+1,ny+GHOST; do i=1,nx
                call ProjectGhostCell(i,j,projection_ok); local_ok=local_ok .and. projection_ok
            enddo; enddo
        endif
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    contains
        subroutine ProjectGhostCell(cell_i,cell_j,cell_ok)
            integer(KINT), intent(in) :: cell_i,cell_j
            logical, intent(out) :: cell_ok
            real(KREAL) :: reference_h(Nvel),reference_b(Nvel)

            reference_h=ctr(cell_i,cell_j)%h_vec
            reference_b=ctr(cell_i,cell_j)%b_vec
            call ProjectToConservedMoments(ctr(cell_i,cell_j)%h_vec,ctr(cell_i,cell_j)%b_vec, &
                ctr(cell_i,cell_j)%conVars,uVel,vVel,wVel,reference_h,reference_b,cell_ok)
        end subroutine ProjectGhostCell
    end subroutine EnforceFreestreamGhostMomentCompatibility

    subroutine SolveMicroLUSGS(ok,physical_dt,flux_weight,collision_weight,freeze_wall_closure)
        logical, intent(out) :: ok
        real(KREAL), intent(in), optional :: physical_dt,collision_weight
        real(KREAL), intent(in), optional :: flux_weight(:,:)
        logical, intent(in), optional :: freeze_wall_closure
        integer(KINT) :: sweep,i,j,k
        real(KREAL) :: start_stage
        logical :: point_ok,local_ok,global_ok,freeze_wall_closure_active
        type(FrozenMicroWallState) :: wall_state

        ok=.false.
        !ImplicitUnsteady keeps its historical per-point wall closure unless explicitly enabled.
        freeze_wall_closure_active=.false.
        if (present(freeze_wall_closure)) freeze_wall_closure_active=freeze_wall_closure
        do sweep=1,Implicit_MicroSweeps
            local_ok=.true.
            start_stage=MPI_Wtime()
            do j=1,ny
                do i=1,nx
                    if (freeze_wall_closure_active) then
                        call PrepareFrozenMicroWallState(i,j,wall_state,point_ok)
                        local_ok=local_ok .and. point_ok
                        if (.not.point_ok) cycle
                        do k=1,Nvel
                            call RelaxMicroPoint(i,j,k,point_ok,physical_dt,flux_weight, &
                                collision_weight,wall_state)
                            local_ok=local_ok .and. point_ok
                        enddo
                    else
                        do k=1,Nvel
                            call RelaxMicroPoint(i,j,k,point_ok,physical_dt,flux_weight,collision_weight)
                            local_ok=local_ok .and. point_ok
                        enddo
                    endif
                enddo
            enddo
            implicit_micro_forward_seconds=implicit_micro_forward_seconds+MPI_Wtime()-start_stage
            start_stage=MPI_Wtime()
            call MPIExchangeDeltaMicro(delta_h,delta_b)
            implicit_micro_forward_backward_exchange_seconds= &
                implicit_micro_forward_backward_exchange_seconds+MPI_Wtime()-start_stage
            implicit_micro_exchange_calls=implicit_micro_exchange_calls+1
            if(s2_steady_rhs)call RecordR1MicroForward(sweep)
            start_stage=MPI_Wtime()
            do j=ny,1,-1
                do i=nx,1,-1
                    if (freeze_wall_closure_active) then
                        call PrepareFrozenMicroWallState(i,j,wall_state,point_ok)
                        local_ok=local_ok .and. point_ok
                        if (.not.point_ok) cycle
                        do k=Nvel,1,-1
                            call RelaxMicroPoint(i,j,k,point_ok,physical_dt,flux_weight, &
                                collision_weight,wall_state)
                            local_ok=local_ok .and. point_ok
                        enddo
                    else
                        do k=Nvel,1,-1
                            call RelaxMicroPoint(i,j,k,point_ok,physical_dt,flux_weight,collision_weight)
                            local_ok=local_ok .and. point_ok
                        enddo
                    endif
                enddo
            enddo
            implicit_micro_backward_seconds=implicit_micro_backward_seconds+MPI_Wtime()-start_stage
            call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
            if (.not.global_ok) return
            if (sweep<Implicit_MicroSweeps) then
                start_stage=MPI_Wtime()
                call MPIExchangeDeltaMicro(delta_h,delta_b)
                implicit_micro_post_backward_exchange_seconds= &
                    implicit_micro_post_backward_exchange_seconds+MPI_Wtime()-start_stage
                implicit_micro_exchange_calls=implicit_micro_exchange_calls+1
            endif
        enddo
        ok=.true.
    end subroutine SolveMicroLUSGS

    subroutine RecordR1MicroForward(sweep)
        integer,intent(in)::sweep
        real(KREAL)::mom(4),tot(4),hdef,bdef
        integer::i,j,env_status,fu
        character(len=8)::env_value
        logical,save::opened=.false.
        env_value=''
        call get_environment_variable('S2I1_R1_DIAG',env_value,status=env_status)
        if(env_status/=0.or.trim(env_value)/='1')return
        tot=0.0_KREAL
        do j=1,ny;do i=1,nx
            call MomentConserved(delta_h(:,i,j),delta_b(:,i,j),uVel,vVel,wVel,mom)
            tot=tot+ctr(i,j)%area*mom
        enddo;enddo
        call ComputeMicroDefectNorm(hdef,bdef)
        if(myID==0)then
            if(.not.opened)then
                open(newunit=fu,file='s2i1_micro_forward.csv',status='replace',action='write')
                write(fu,'(A)')'iter,retry,sweep,dMf_sum1,dMf_sum2,dMf_sum3,dMf_sum4,h_linear_defect,b_linear_defect'
                close(fu);opened=.true.
            endif
            open(newunit=fu,file='s2i1_micro_forward.csv',status='old',position='append',action='write')
            write(fu,'(I0,",",I0,",",I0,6(",",ES24.16))')iter,implicit_retry_count,sweep,tot,hdef,bdef
            close(fu)
        endif
    end subroutine RecordR1MicroForward

    subroutine RelaxMicroPoint(i,j,k,ok,physical_dt,flux_weight,collision_weight,wall_state)
        integer(KINT), intent(in) :: i,j,k
        logical, intent(out) :: ok
        real(KREAL), intent(in), optional :: physical_dt,collision_weight
        real(KREAL), intent(in), optional :: flux_weight(:,:)
        type(FrozenMicroWallState), intent(in), optional :: wall_state
        integer(KINT) :: face_number
        real(KREAL) :: rhs_h,rhs_b,diagonal,dt_pseudo,un,face_length,neighbor_h,neighbor_b
        real(KREAL) :: stream_coefficient,collision_coefficient
        logical :: face_ok

        if (present(physical_dt) .and. present(flux_weight) .and. present(collision_weight)) then
            stream_coefficient=flux_weight(i,j)
            collision_coefficient=collision_weight
            rhs_h=stream_h(k,i,j)
            rhs_b=stream_b(k,i,j)
            diagonal=1.0_KREAL/physical_dt+collision_coefficient/tau_tilde(i,j)
        else
            stream_coefficient=1.0_KREAL
            collision_coefficient=1.0_KREAL
            dt_pseudo=PseudoTimeStep(ctr(i,j)%dt_local,implicit_ncfl_active,CFL)
            if(s2_steady_rhs)then
                if(s2_steady_coupling_on)then
                    rhs_h=stream_h(k,i,j)+equilibrium_increment_h(k,i,j)/tau_tilde(i,j)
                    rhs_b=stream_b(k,i,j)+equilibrium_increment_b(k,i,j)/tau_tilde(i,j)
                else
                    rhs_h=stream_h(k,i,j)
                    rhs_b=stream_b(k,i,j)
                endif
            else
                rhs_h=stream_h(k,i,j)+(equilibrium_h(k,i,j)-ctr(i,j)%h_vec(k))/tau_tilde(i,j)
                rhs_b=stream_b(k,i,j)+(equilibrium_b(k,i,j)-ctr(i,j)%b_vec(k))/tau_tilde(i,j)
            endif
            diagonal=1.0_KREAL/dt_pseudo+1.0_KREAL/tau_tilde(i,j)
        endif
        do face_number=1,4
            call GetMicroFaceValue(i,j,k,face_number,un,face_length,neighbor_h,neighbor_b, &
                face_ok,wall_state)
            if (.not.face_ok) then; ok=.false.; return; endif
            diagonal=diagonal+stream_coefficient*max(un,0.0_KREAL)*face_length/ctr(i,j)%area
            rhs_h=rhs_h-stream_coefficient*min(un,0.0_KREAL)*face_length*neighbor_h/ctr(i,j)%area
            rhs_b=rhs_b-stream_coefficient*min(un,0.0_KREAL)*face_length*neighbor_b/ctr(i,j)%area
        enddo
        ok=diagonal>0.0_KREAL .and. ieee_is_finite(diagonal) .and. &
           ieee_is_finite(rhs_h) .and. ieee_is_finite(rhs_b)
        if (ok) then
            delta_h(k,i,j)=rhs_h/diagonal
            delta_b(k,i,j)=rhs_b/diagonal
        endif
    end subroutine RelaxMicroPoint

    subroutine EvaluateMicroPointBreakdown(i,j,k,is_h,current_value,stored_delta, &
        equilibrium_value,stream_rhs,collision_rhs,neighbor_rhs,dt_pseudo,tau_value, &
        diagonal_total,rhs_total,post_sweep_delta,ok)
        integer(KINT), intent(in) :: i,j,k
        logical, intent(in) :: is_h
        real(KREAL), intent(out) :: current_value,stored_delta,equilibrium_value,stream_rhs
        real(KREAL), intent(out) :: collision_rhs,neighbor_rhs,dt_pseudo,tau_value
        real(KREAL), intent(out) :: diagonal_total,rhs_total,post_sweep_delta
        logical, intent(out) :: ok
        integer(KINT) :: face_number
        real(KREAL) :: un,face_length,neighbor_h,neighbor_b,neighbor_delta,area
        logical :: face_ok

        current_value=0.0_KREAL; stored_delta=0.0_KREAL
        equilibrium_value=0.0_KREAL; stream_rhs=0.0_KREAL
        collision_rhs=0.0_KREAL; neighbor_rhs=0.0_KREAL
        dt_pseudo=0.0_KREAL; tau_value=0.0_KREAL
        diagonal_total=0.0_KREAL; rhs_total=0.0_KREAL
        post_sweep_delta=0.0_KREAL; ok=.false.

        if (i<1 .or. i>nx .or. j<1 .or. j>ny .or. k<1 .or. k>Nvel) return
        if (.not.allocated(ctr) .or. .not.allocated(vface) .or. .not.allocated(hface)) return
        if (.not.allocated(delta_h) .or. .not.allocated(delta_b) .or. &
            .not.allocated(stream_h) .or. .not.allocated(stream_b) .or. &
            .not.allocated(tau_tilde) .or. .not.allocated(equilibrium_h) .or. &
            .not.allocated(equilibrium_b)) return
        if (.not.allocated(uVel) .or. .not.allocated(vVel)) return
        if (.not.associated(ctr(i,j)%h_vec) .or. .not.associated(ctr(i,j)%b_vec)) return
        if (.not.ieee_is_finite(ctr(i,j)%dt_local) .or. ctr(i,j)%dt_local<=0.0_KREAL) return
        if (.not.ieee_is_finite(implicit_ncfl_active) .or. implicit_ncfl_active<=0.0_KREAL) return
        if (.not.ieee_is_finite(CFL) .or. CFL<=0.0_KREAL) return

        if (is_h) then
            current_value=ctr(i,j)%h_vec(k)
            stored_delta=delta_h(k,i,j)
            equilibrium_value=equilibrium_h(k,i,j)
            stream_rhs=stream_h(k,i,j)
        else
            current_value=ctr(i,j)%b_vec(k)
            stored_delta=delta_b(k,i,j)
            equilibrium_value=equilibrium_b(k,i,j)
            stream_rhs=stream_b(k,i,j)
        endif
        dt_pseudo=PseudoTimeStep(ctr(i,j)%dt_local,implicit_ncfl_active,CFL)
        tau_value=tau_tilde(i,j)
        if (.not.all(ieee_is_finite([current_value,stored_delta,equilibrium_value,stream_rhs, &
            dt_pseudo,tau_value])) .or. dt_pseudo<=0.0_KREAL .or. tau_value<=0.0_KREAL) return

        collision_rhs=(equilibrium_value-current_value)/tau_value
        area=ctr(i,j)%area
        if (.not.ieee_is_finite(area) .or. area<=0.0_KREAL) return
        neighbor_rhs=0.0_KREAL
        diagonal_total=1.0_KREAL/dt_pseudo+1.0_KREAL/tau_value
        do face_number=1,4
            call GetMicroFaceValue(i,j,k,face_number,un,face_length,neighbor_h,neighbor_b,face_ok)
            if (.not.face_ok) return
            if (.not.all(ieee_is_finite([un,face_length,neighbor_h,neighbor_b]))) return
            if (is_h) then
                neighbor_delta=neighbor_h
            else
                neighbor_delta=neighbor_b
            endif
            neighbor_rhs=neighbor_rhs-min(un,0.0_KREAL)*face_length*neighbor_delta/area
            diagonal_total=diagonal_total+max(un,0.0_KREAL)*face_length/area
        enddo
        rhs_total=stream_rhs+collision_rhs+neighbor_rhs
        if (.not.ieee_is_finite(diagonal_total) .or. diagonal_total<=0.0_KREAL) return
        post_sweep_delta=rhs_total/diagonal_total
        ok=all(ieee_is_finite([current_value,stored_delta,equilibrium_value,stream_rhs, &
            collision_rhs,neighbor_rhs,dt_pseudo,tau_value,diagonal_total,rhs_total, &
            post_sweep_delta]))
    end subroutine EvaluateMicroPointBreakdown

    real(KREAL) function EvaluateMicroHDefect(i,j,k) result(defect)
        integer(KINT), intent(in) :: i,j,k
        defect=EvaluateMicroDefectComponent(i,j,k,.true.)
    end function EvaluateMicroHDefect

    real(KREAL) function EvaluateMicroBDefect(i,j,k) result(defect)
        integer(KINT), intent(in) :: i,j,k
        defect=EvaluateMicroDefectComponent(i,j,k,.false.)
    end function EvaluateMicroBDefect

    real(KREAL) function EvaluateMicroDefectComponent(i,j,k,is_h) result(defect)
        integer(KINT), intent(in) :: i,j,k
        logical, intent(in) :: is_h
        integer(KINT) :: face_number
        real(KREAL) :: rhs,lhs,dt_pseudo,un,face_length,neighbor_h,neighbor_b,current_delta
        logical :: face_ok

        dt_pseudo=PseudoTimeStep(ctr(i,j)%dt_local,implicit_ncfl_active,CFL)
        if (is_h) then
            rhs=stream_h(k,i,j)+(equilibrium_h(k,i,j)-ctr(i,j)%h_vec(k))/tau_tilde(i,j)
            if(s2_steady_rhs)then
                rhs=stream_h(k,i,j)
                if(s2_steady_coupling_on)rhs=rhs+equilibrium_increment_h(k,i,j)/tau_tilde(i,j)
            endif
            current_delta=delta_h(k,i,j)
        else
            rhs=stream_b(k,i,j)+(equilibrium_b(k,i,j)-ctr(i,j)%b_vec(k))/tau_tilde(i,j)
            if(s2_steady_rhs)then
                rhs=stream_b(k,i,j)
                if(s2_steady_coupling_on)rhs=rhs+equilibrium_increment_b(k,i,j)/tau_tilde(i,j)
            endif
            current_delta=delta_b(k,i,j)
        endif
        lhs=current_delta*(1.0_KREAL/dt_pseudo+1.0_KREAL/tau_tilde(i,j))
        do face_number=1,4
            call GetMicroFaceValue(i,j,k,face_number,un,face_length,neighbor_h,neighbor_b,face_ok)
            if (.not.face_ok) then; defect=huge(1.0_KREAL); return; endif
            if (is_h) then
                lhs=lhs+face_length*(max(un,0.0_KREAL)*current_delta+min(un,0.0_KREAL)*neighbor_h)/ctr(i,j)%area
            else
                lhs=lhs+face_length*(max(un,0.0_KREAL)*current_delta+min(un,0.0_KREAL)*neighbor_b)/ctr(i,j)%area
            endif
        enddo
        defect=rhs-lhs
    end function EvaluateMicroDefectComponent

    subroutine ComputeMicroDefectNorm(norm_h,norm_b)
        real(KREAL), intent(out) :: norm_h,norm_b
        integer(KINT) :: i,j,k
        real(KREAL) :: local_sums(2),global_sums(2),dh,db
        local_sums=0.0_KREAL
        do j=1,ny; do i=1,nx; do k=1,Nvel
            dh=EvaluateMicroHDefect(i,j,k); db=EvaluateMicroBDefect(i,j,k)
            local_sums(1)=local_sums(1)+wVel(k)*dh*dh
            local_sums(2)=local_sums(2)+wVel(k)*db*db
        enddo; enddo; enddo
        call MPI_Allreduce(local_sums,global_sums,2,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        norm_h=sqrt(global_sums(1)); norm_b=sqrt(global_sums(2))
    end subroutine ComputeMicroDefectNorm

    subroutine RecoverMacroFromKineticState(ok)
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        real(KREAL) :: recovered(4)
        logical :: local_ok,global_ok
        local_ok=.true.
        do j=1,ny; do i=1,nx
            call MomentConserved(ctr(i,j)%h_vec,ctr(i,j)%b_vec,uVel,vVel,wVel,recovered)
            local_ok=local_ok .and. IsPhysicalState(recovered)
            w_tilde(:,i,j)=recovered
        enddo; enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
        if (.not.ok) return
        do j=1,ny; do i=1,nx
            ctr(i,j)%conVars=w_tilde(:,i,j)
        enddo; enddo
    end subroutine RecoverMacroFromKineticState

    subroutine PrepareFrozenMicroWallState(i,j,state,ok)
        integer(KINT), intent(in) :: i,j
        type(FrozenMicroWallState), intent(inout) :: state
        logical, intent(out) :: ok
        integer(KINT) :: face_number,btype,side
        real(KREAL) :: nx_out,ny_out,face_length,rho_delta
        logical :: physical,face_ok

        ok=.false.
        state%active=.false.
        state%fallback=.false.
        state%build_calls=0
        if (i<1 .or. i>nx .or. j<1 .or. j>ny) return
        if (.not.allocated(state%delta_h_face)) then
            allocate(state%delta_h_face(Nvel,4),state%delta_b_face(Nvel,4))
        elseif (size(state%delta_h_face,1)/=Nvel) then
            deallocate(state%delta_h_face,state%delta_b_face)
            allocate(state%delta_h_face(Nvel,4),state%delta_b_face(Nvel,4))
        endif
        do face_number=1,4
            physical=.false.; btype=0; side=Rside
            select case(face_number)
            case(1)
                nx_out=-vface(i,j)%cosx; ny_out=-vface(i,j)%cosy
                face_length=vface(i,j)%length
                physical=(i==1 .and. ID_X_minus==MPI_PROC_NULL)
                btype=BC_Type(1); side=Lside
                if (physical) btype=ResolveImplicitBoundaryType(i,j,IDIRC,Lside,btype)
            case(2)
                nx_out=vface(i+1,j)%cosx; ny_out=vface(i+1,j)%cosy
                face_length=vface(i+1,j)%length
                physical=(i==nx .and. ID_X_plus==MPI_PROC_NULL)
                btype=BC_Type(2); side=Rside
                if (physical) btype=ResolveImplicitBoundaryType(i,j,IDIRC,Rside,btype)
            case(3)
                nx_out=-hface(i,j)%cosx; ny_out=-hface(i,j)%cosy
                face_length=hface(i,j)%length
                physical=(j==1 .and. ID_Y_minus==MPI_PROC_NULL)
                btype=BC_Type(3); side=Lside
                if (physical) btype=ResolveImplicitBoundaryType(i,j,JDIRC,Lside,btype)
            case default
                nx_out=hface(i,j+1)%cosx; ny_out=hface(i,j+1)%cosy
                face_length=hface(i,j+1)%length
                physical=(j==ny .and. ID_Y_plus==MPI_PROC_NULL)
                btype=BC_Type(4); side=Rside
                if (physical) btype=ResolveImplicitBoundaryType(i,j,JDIRC,Rside,btype)
            end select
            if (.not.physical .or. btype/=BCconst_wall_isothermal) cycle
            if (.not.ieee_is_finite(face_length)) return
            if (.not.IsAxisAlignedSymmetryNormal(nx_out,ny_out)) then
                state%fallback(face_number)=.true.
                cycle
            endif
            !Build the complete wall correction once; all velocity points in this cell sweep reuse it.
            if (FluxWall_type==FluxWall_Maxwell) then
                call ComputeMaxwellWallFaceDelta(Prim_wall,nx_out,ny_out,side,ACC, &
                    delta_h(:,i,j),delta_b(:,i,j),state%delta_h_face(:,face_number), &
                    state%delta_b_face(:,face_number),rho_delta,face_ok)
            else
                call ComputeDiffuseWallFaceDelta(Prim_wall,nx_out,ny_out,side, &
                    delta_h(:,i,j),delta_b(:,i,j),state%delta_h_face(:,face_number), &
                    state%delta_b_face(:,face_number),rho_delta,face_ok)
            endif
            if (.not.face_ok) return
            state%active(face_number)=.true.
            state%build_calls=state%build_calls+1
        enddo
        ok=.true.
    end subroutine PrepareFrozenMicroWallState

    subroutine GetMicroFaceValue(i,j,k,face_number,un,face_length,neighbor_h,neighbor_b,ok,wall_state)
        integer(KINT), intent(in) :: i,j,k,face_number
        real(KREAL), intent(out) :: un,face_length,neighbor_h,neighbor_b
        logical, intent(out) :: ok
        type(FrozenMicroWallState), intent(in), optional :: wall_state
        integer(KINT) :: ni,nj,btype,side
        logical :: physical
        real(KREAL) :: nx_out,ny_out,dh_boundary(Nvel),db_boundary(Nvel),rho_delta

        ni=i; nj=j; physical=.false.; btype=0; side=Rside
        select case(face_number)
        case(1)
            ni=i-1; nx_out=-vface(i,j)%cosx; ny_out=-vface(i,j)%cosy; face_length=vface(i,j)%length
            physical=(i==1 .and. ID_X_minus==MPI_PROC_NULL); btype=BC_Type(1); side=Lside
        case(2)
            ni=i+1; nx_out=vface(i+1,j)%cosx; ny_out=vface(i+1,j)%cosy; face_length=vface(i+1,j)%length
            physical=(i==nx .and. ID_X_plus==MPI_PROC_NULL); btype=BC_Type(2); side=Rside
        case(3)
            nj=j-1; nx_out=-hface(i,j)%cosx; ny_out=-hface(i,j)%cosy; face_length=hface(i,j)%length
            physical=(j==1 .and. ID_Y_minus==MPI_PROC_NULL); btype=BC_Type(3); side=Lside
        case default
            nj=j+1; nx_out=hface(i,j+1)%cosx; ny_out=hface(i,j+1)%cosy; face_length=hface(i,j+1)%length
            physical=(j==ny .and. ID_Y_plus==MPI_PROC_NULL); btype=BC_Type(4); side=Rside
        end select
        un=uVel(k)*nx_out+vVel(k)*ny_out; ok=.true.
        if (.not.physical) then
            neighbor_h=delta_h(k,ni,nj); neighbor_b=delta_b(k,ni,nj); return
        endif
        select case(face_number)
        case(1); btype=ResolveImplicitBoundaryType(i,j,IDIRC,Lside,btype)
        case(2); btype=ResolveImplicitBoundaryType(i,j,IDIRC,Rside,btype)
        case(3); btype=ResolveImplicitBoundaryType(i,j,JDIRC,Lside,btype)
        case default; btype=ResolveImplicitBoundaryType(i,j,JDIRC,Rside,btype)
        end select
        select case(btype)
        case(BCconst_freestream,BCconst_inflow)
            neighbor_h=0.0_KREAL; neighbor_b=0.0_KREAL
        case(BCconst_extrapolation)
            neighbor_h=delta_h(k,i,j); neighbor_b=delta_b(k,i,j)
        case(BCconst_symmetry)
            call GetSymmetryMicroValue(delta_h(:,i,j),nx_out,ny_out,k,neighbor_h,ok)
            if (ok) call GetSymmetryMicroValue(delta_b(:,i,j),nx_out,ny_out,k,neighbor_b,ok)
        case(BCconst_wall_isothermal)
            if (present(wall_state)) then
                if (wall_state%active(face_number)) then
                    if (.not.allocated(wall_state%delta_h_face) .or. &
                        .not.allocated(wall_state%delta_b_face) .or. &
                        size(wall_state%delta_h_face,1)/=Nvel .or. &
                        size(wall_state%delta_b_face,1)/=Nvel) then
                        ok=.false.
                        return
                    endif
                    neighbor_h=wall_state%delta_h_face(k,face_number)
                    neighbor_b=wall_state%delta_b_face(k,face_number)
                    return
                elseif (.not.wall_state%fallback(face_number)) then
                    ok=.false.
                    return
                endif
            endif
            if (FluxWall_type==FluxWall_Maxwell) then
                call ComputeMaxwellWallFaceDelta(Prim_wall,nx_out,ny_out,side,ACC, &
                    delta_h(:,i,j),delta_b(:,i,j),dh_boundary,db_boundary,rho_delta,ok)
            else
                call ComputeDiffuseWallFaceDelta(Prim_wall,nx_out,ny_out,side, &
                    delta_h(:,i,j),delta_b(:,i,j),dh_boundary,db_boundary,rho_delta,ok)
            endif
            if (ok) then
                neighbor_h=dh_boundary(k)
                neighbor_b=db_boundary(k)
            endif
        case default
            ok=.false.
        end select
    end subroutine GetMicroFaceValue
end module mod_ImplicitMicro

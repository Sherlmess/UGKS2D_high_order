module mod_ImplicitState
    use flow_data
    implicit none
    integer(KINT), parameter :: IMPLICIT_FAILURE_NONE=0
    integer(KINT), parameter :: IMPLICIT_FAILURE_H_NONFINITE=1
    integer(KINT), parameter :: IMPLICIT_FAILURE_B_NONFINITE=2
    integer(KINT), parameter :: IMPLICIT_FAILURE_H_NEGATIVE_FRACTION=3
    integer(KINT), parameter :: IMPLICIT_FAILURE_B_NEGATIVE_FRACTION=4
    integer(KINT), parameter :: IMPLICIT_FAILURE_MOMENT_NONPHYSICAL=5
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE=0
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS=1
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_INDEX=2
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_UNAVAILABLE_STATE=3
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_SHAPE=4
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONPHYSICAL_W_TILDE=5
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_INPUT=6
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_MAXWELL_OR_SHAKHOV=7
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_PROJECTION_FAILED=8
    integer(KINT), parameter :: IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH=9

    type :: ImplicitFailureRecord
        logical :: valid=.false.
        integer(KINT) :: reason=IMPLICIT_FAILURE_NONE
        integer(KINT) :: owner_rank=-1
        integer(KINT) :: local_i=0,local_j=0
        integer(KINT) :: global_i=0,global_j=0
        integer(KINT) :: velocity_k=0
        real(KREAL) :: velocity_u=0.0_KREAL,velocity_v=0.0_KREAL
        real(KREAL) :: ncfl_attempt=0.0_KREAL
        real(KREAL) :: alpha_tested=0.0_KREAL,alpha_next=0.0_KREAL
        real(KREAL) :: current_fraction=0.0_KREAL
        real(KREAL) :: candidate_fraction=0.0_KREAL
        real(KREAL) :: allowed_fraction=0.0_KREAL
        real(KREAL) :: candidate_moment(4)=0.0_KREAL
        real(KREAL) :: candidate_internal_energy=0.0_KREAL
        logical :: micro_breakdown_valid=.false.
        logical :: equilibrium_breakdown_valid=.false.
        integer(KINT) :: equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE
        real(KREAL) :: micro_current=0.0_KREAL
        real(KREAL) :: micro_delta=0.0_KREAL
        real(KREAL) :: micro_candidate=0.0_KREAL
        real(KREAL) :: equilibrium_value=0.0_KREAL
        real(KREAL) :: equilibrium_maxwell_value=0.0_KREAL
        real(KREAL) :: equilibrium_shakhov_correction=0.0_KREAL
        real(KREAL) :: equilibrium_preprojection_value=0.0_KREAL
        real(KREAL) :: equilibrium_projection_correction=0.0_KREAL
        real(KREAL) :: equilibrium_final_value=0.0_KREAL
        real(KREAL) :: q_frozen_x=0.0_KREAL,q_frozen_y=0.0_KREAL
        real(KREAL) :: prim_tilde_rho=0.0_KREAL,prim_tilde_u=0.0_KREAL
        real(KREAL) :: prim_tilde_v=0.0_KREAL,prim_tilde_temperature=0.0_KREAL
        real(KREAL) :: equilibrium_maxwell_negative_fraction=0.0_KREAL
        real(KREAL) :: equilibrium_preprojection_negative_fraction=0.0_KREAL
        real(KREAL) :: equilibrium_final_negative_fraction=0.0_KREAL
        real(KREAL) :: reconstruction_max_abs_h=0.0_KREAL
        real(KREAL) :: reconstruction_max_abs_b=0.0_KREAL
        real(KREAL) :: reconstruction_max_rel_h=0.0_KREAL
        real(KREAL) :: reconstruction_max_rel_b=0.0_KREAL
        integer(KINT) :: reconstruction_max_abs_h_k=0
        integer(KINT) :: reconstruction_max_abs_b_k=0
        real(KREAL) :: reconstructed_h_at_failure_k=0.0_KREAL
        real(KREAL) :: stored_h_at_failure_k=0.0_KREAL
        real(KREAL) :: reconstructed_b_at_failure_k=0.0_KREAL
        real(KREAL) :: stored_b_at_failure_k=0.0_KREAL
        real(KREAL) :: stream_rhs=0.0_KREAL
        real(KREAL) :: collision_rhs=0.0_KREAL
        real(KREAL) :: neighbor_rhs=0.0_KREAL
        real(KREAL) :: dt_pseudo=0.0_KREAL
        real(KREAL) :: tau_value=0.0_KREAL
        real(KREAL) :: diagonal_total=0.0_KREAL
        real(KREAL) :: rhs_total=0.0_KREAL
        real(KREAL) :: post_sweep_delta=0.0_KREAL
    end type ImplicitFailureRecord

    type(ImplicitFailureRecord) :: implicit_failure
    real(KREAL), allocatable :: rate_w(:,:,:)
    real(KREAL), allocatable :: stream_h(:,:,:), stream_b(:,:,:)
    real(KREAL), allocatable :: delta_w(:,:,:), delta_h(:,:,:), delta_b(:,:,:)
    real(KREAL), allocatable :: w_tilde(:,:,:), tau_tilde(:,:), q_frozen(:,:,:)
    real(KREAL), allocatable :: equilibrium_h(:,:,:), equilibrium_b(:,:,:)
    real(KREAL), allocatable :: equilibrium_increment_h(:,:,:), equilibrium_increment_b(:,:,:)
    logical :: s2_steady_rhs=.false.
    logical :: s2_steady_coupling_on=.false.
    integer(KINT), allocatable :: mirror_u_index(:), mirror_v_index(:)
    logical :: mirror_u_valid=.false., mirror_v_valid=.false.
    real(KREAL) :: implicit_ncfl_active=0.0_KREAL, implicit_alpha=1.0_KREAL
    real(KREAL) :: implicit_macro_residual(4)=0.0_KREAL
    real(KREAL) :: implicit_macro_relative(4)=0.0_KREAL
    real(KREAL) :: implicit_macro_absolute(4)=0.0_KREAL
    real(KREAL) :: implicit_h_norm=0.0_KREAL, implicit_b_norm=0.0_KREAL
    real(KREAL) :: implicit_h_relative=0.0_KREAL, implicit_b_relative=0.0_KREAL
    real(KREAL) :: implicit_h_absolute=0.0_KREAL, implicit_b_absolute=0.0_KREAL
    real(KREAL) :: implicit_moment_error=0.0_KREAL, implicit_mass_defect=0.0_KREAL
    real(KREAL) :: implicit_pr_seconds=0.0_KREAL, implicit_lusgs_seconds=0.0_KREAL
    real(KREAL) :: implicit_exchange_seconds=0.0_KREAL
    real(KREAL) :: implicit_micro_forward_seconds=0.0_KREAL
    real(KREAL) :: implicit_micro_backward_seconds=0.0_KREAL
    real(KREAL) :: implicit_micro_initial_exchange_seconds=0.0_KREAL
    real(KREAL) :: implicit_micro_forward_backward_exchange_seconds=0.0_KREAL
    real(KREAL) :: implicit_micro_post_backward_exchange_seconds=0.0_KREAL
    real(KREAL) :: implicit_evolution_seconds=0.0_KREAL, implicit_total_seconds=0.0_KREAL
    integer(KINT) :: implicit_retry_count=0, implicit_retry_total=0
    integer(KINT) :: implicit_micro_exchange_calls=0
contains
    subroutine AllocateImplicitState()
        if (allocated(rate_w)) return
        allocate(rate_w(4,1:nx,1:ny), stream_h(Nvel,1:nx,1:ny), stream_b(Nvel,1:nx,1:ny))
        allocate(delta_w(4,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(delta_h(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(delta_b(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(w_tilde(4,1:nx,1:ny), tau_tilde(1:nx,1:ny), q_frozen(2,1:nx,1:ny))
        allocate(equilibrium_h(Nvel,1:nx,1:ny), equilibrium_b(Nvel,1:nx,1:ny))
        allocate(equilibrium_increment_h(Nvel,1:nx,1:ny),equilibrium_increment_b(Nvel,1:nx,1:ny))
        if (.not.allocated(mirror_u_index)) allocate(mirror_u_index(Nvel), mirror_v_index(Nvel))
        call ResetImplicitCorrections()
        rate_w=0.0_KREAL; stream_h=0.0_KREAL; stream_b=0.0_KREAL
        w_tilde=0.0_KREAL; tau_tilde=0.0_KREAL; q_frozen=0.0_KREAL
        equilibrium_h=0.0_KREAL; equilibrium_b=0.0_KREAL
        equilibrium_increment_h=0.0_KREAL;equilibrium_increment_b=0.0_KREAL
        implicit_ncfl_active=Implicit_nCFL
    end subroutine AllocateImplicitState

    subroutine ResetImplicitCorrections()
        delta_w=0.0_KREAL; delta_h=0.0_KREAL; delta_b=0.0_KREAL
        implicit_alpha=1.0_KREAL
        call ResetImplicitFailureRecord()
    end subroutine ResetImplicitCorrections

    subroutine ResetImplicitFailureRecord(record)
        type(ImplicitFailureRecord), intent(out), optional :: record
        if (present(record)) then
            record=ImplicitFailureRecord()
            record%micro_breakdown_valid=.false.
            record%equilibrium_breakdown_valid=.false.
            record%equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE
            record%micro_current=0.0_KREAL; record%micro_delta=0.0_KREAL
            record%micro_candidate=0.0_KREAL; record%equilibrium_value=0.0_KREAL
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
            record%stream_rhs=0.0_KREAL; record%collision_rhs=0.0_KREAL
            record%neighbor_rhs=0.0_KREAL; record%dt_pseudo=0.0_KREAL
            record%tau_value=0.0_KREAL; record%diagonal_total=0.0_KREAL
            record%rhs_total=0.0_KREAL; record%post_sweep_delta=0.0_KREAL
        else
            implicit_failure=ImplicitFailureRecord()
            implicit_failure%micro_breakdown_valid=.false.
            implicit_failure%equilibrium_breakdown_valid=.false.
            implicit_failure%equilibrium_breakdown_stage=IMPLICIT_EQUILIBRIUM_BREAKDOWN_NOT_APPLICABLE
            implicit_failure%micro_current=0.0_KREAL; implicit_failure%micro_delta=0.0_KREAL
            implicit_failure%micro_candidate=0.0_KREAL; implicit_failure%equilibrium_value=0.0_KREAL
            implicit_failure%equilibrium_maxwell_value=0.0_KREAL
            implicit_failure%equilibrium_shakhov_correction=0.0_KREAL
            implicit_failure%equilibrium_preprojection_value=0.0_KREAL
            implicit_failure%equilibrium_projection_correction=0.0_KREAL
            implicit_failure%equilibrium_final_value=0.0_KREAL
            implicit_failure%q_frozen_x=0.0_KREAL; implicit_failure%q_frozen_y=0.0_KREAL
            implicit_failure%prim_tilde_rho=0.0_KREAL; implicit_failure%prim_tilde_u=0.0_KREAL
            implicit_failure%prim_tilde_v=0.0_KREAL; implicit_failure%prim_tilde_temperature=0.0_KREAL
            implicit_failure%equilibrium_maxwell_negative_fraction=0.0_KREAL
            implicit_failure%equilibrium_preprojection_negative_fraction=0.0_KREAL
            implicit_failure%equilibrium_final_negative_fraction=0.0_KREAL
            implicit_failure%reconstruction_max_abs_h=0.0_KREAL
            implicit_failure%reconstruction_max_abs_b=0.0_KREAL
            implicit_failure%reconstruction_max_rel_h=0.0_KREAL
            implicit_failure%reconstruction_max_rel_b=0.0_KREAL
            implicit_failure%reconstruction_max_abs_h_k=0
            implicit_failure%reconstruction_max_abs_b_k=0
            implicit_failure%reconstructed_h_at_failure_k=0.0_KREAL
            implicit_failure%stored_h_at_failure_k=0.0_KREAL
            implicit_failure%reconstructed_b_at_failure_k=0.0_KREAL
            implicit_failure%stored_b_at_failure_k=0.0_KREAL
            implicit_failure%stream_rhs=0.0_KREAL; implicit_failure%collision_rhs=0.0_KREAL
            implicit_failure%neighbor_rhs=0.0_KREAL; implicit_failure%dt_pseudo=0.0_KREAL
            implicit_failure%tau_value=0.0_KREAL; implicit_failure%diagonal_total=0.0_KREAL
            implicit_failure%rhs_total=0.0_KREAL; implicit_failure%post_sweep_delta=0.0_KREAL
        endif
    end subroutine ResetImplicitFailureRecord

    pure function ImplicitFailureReasonLabel(reason) result(label)
        integer(KINT), intent(in) :: reason
        character(len=32) :: label
        select case(reason)
        case(IMPLICIT_FAILURE_H_NONFINITE)
            label='h_nonfinite'
        case(IMPLICIT_FAILURE_B_NONFINITE)
            label='b_nonfinite'
        case(IMPLICIT_FAILURE_H_NEGATIVE_FRACTION)
            label='h_negative_fraction'
        case(IMPLICIT_FAILURE_B_NEGATIVE_FRACTION)
            label='b_negative_fraction'
        case(IMPLICIT_FAILURE_MOMENT_NONPHYSICAL)
            label='moment_nonphysical'
        case default
            label='none'
        end select
    end function ImplicitFailureReasonLabel

    pure function ImplicitEquilibriumBreakdownStageLabel(stage) result(label)
        integer(KINT), intent(in) :: stage
        character(len=40) :: label
        select case(stage)
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS)
            label='success'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_INDEX)
            label='invalid_index'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_UNAVAILABLE_STATE)
            label='unavailable_state'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_INVALID_SHAPE)
            label='invalid_shape'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONPHYSICAL_W_TILDE)
            label='nonphysical_w_tilde'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_INPUT)
            label='nonfinite_input'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_NONFINITE_MAXWELL_OR_SHAKHOV)
            label='nonfinite_maxwell_or_shakhov'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_PROJECTION_FAILED)
            label='projection_failed'
        case(IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH)
            label='stored_reconstruction_mismatch'
        case default
            label='not_applicable'
        end select
    end function ImplicitEquilibriumBreakdownStageLabel

    subroutine ReleaseImplicitState()
        if (allocated(rate_w)) then
            deallocate(rate_w, stream_h, stream_b, delta_w, delta_h, delta_b)
            deallocate(w_tilde, tau_tilde, q_frozen, equilibrium_h, equilibrium_b)
        endif
        if (allocated(mirror_u_index)) deallocate(mirror_u_index, mirror_v_index)
        mirror_u_valid=.false.; mirror_v_valid=.false.
    end subroutine ReleaseImplicitState
end module mod_ImplicitState

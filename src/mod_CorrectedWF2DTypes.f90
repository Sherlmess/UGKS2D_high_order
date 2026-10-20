module mod_CorrectedWF2DTypes
    use Precision, only: KINT, KREAL
    implicit none
    private

    public :: invalidate_corrected_wf_workspace

    type, public :: GasThermo2D
        integer(KINT) :: ck=-1
    end type GasThermo2D

    type, abstract, public :: CollisionPointModel2D
        type(GasThermo2D) :: thermo
    contains
        procedure(point_target_value_ifc), deferred :: point_target_value
        procedure(point_target_direction_ifc), deferred :: point_target_direction
    end type CollisionPointModel2D

    type, public :: CorrectedWF2DDiagnostics
        logical :: success=.false.
        logical :: snapshot_hit=.false.
        integer(KINT) :: failure_code=0
        integer(KINT) :: stage=0
        character(len=160) :: phase=''
        character(len=256) :: reason=''
        real(KREAL) :: cfl=0.0_KREAL
        real(KREAL) :: dt_over_tau=0.0_KREAL
        real(KREAL) :: moment_error=0.0_KREAL
    end type CorrectedWF2DDiagnostics

    type, public :: CorrectedWF2DWorkspace
        logical :: snapshot_valid=.false.
        integer(KINT) :: primal_capture_count=0
        integer(KINT) :: direction_hit_count=0
        integer(KINT) :: invalidation_count=0
        integer(KINT) :: reject_count=0
        integer(KINT) :: generation=0
        integer(KINT) :: nx=0,ny=0,nvel=0,ck=-1
        real(KREAL) :: dx=0.0_KREAL,dy=0.0_KREAL,physical_dt=0.0_KREAL,tau0=0.0_KREAL
        real(KREAL), allocatable :: w_key(:,:,:),h_key(:,:,:),b_key(:,:,:)
        real(KREAL), allocatable :: u_key(:),v_key(:),weight_key(:)
    end type CorrectedWF2DWorkspace

    abstract interface
        subroutine point_target_value_ifc(self,w_point,h_point,b_point,uvel,vvel,weight, &
            g_h,g_b,tau,ok)
            import :: CollisionPointModel2D,KREAL
            class(CollisionPointModel2D), intent(in) :: self
            real(KREAL), intent(in) :: w_point(4),h_point(:),b_point(:)
            real(KREAL), intent(in) :: uvel(:),vvel(:),weight(:)
            real(KREAL), intent(out) :: g_h(:),g_b(:),tau
            logical, intent(out) :: ok
        end subroutine point_target_value_ifc

        subroutine point_target_direction_ifc(self,w_point,h_point,b_point,w_direction, &
            h_direction,b_direction,uvel,vvel,weight,g_h,g_b,tau,g_h_direction, &
            g_b_direction,tau_direction,ok)
            import :: CollisionPointModel2D,KREAL
            class(CollisionPointModel2D), intent(in) :: self
            real(KREAL), intent(in) :: w_point(4),h_point(:),b_point(:)
            real(KREAL), intent(in) :: w_direction(4),h_direction(:),b_direction(:)
            real(KREAL), intent(in) :: uvel(:),vvel(:),weight(:)
            real(KREAL), intent(out) :: g_h(:),g_b(:),tau
            real(KREAL), intent(out) :: g_h_direction(:),g_b_direction(:),tau_direction
            logical, intent(out) :: ok
        end subroutine point_target_direction_ifc
    end interface

contains

    subroutine invalidate_corrected_wf_workspace(workspace,rejected)
        type(CorrectedWF2DWorkspace), intent(inout) :: workspace
        logical, intent(in), optional :: rejected
        workspace%snapshot_valid=.false.
        workspace%invalidation_count=workspace%invalidation_count+1
        if (present(rejected)) then
            if (rejected) workspace%reject_count=workspace%reject_count+1
        end if
    end subroutine invalidate_corrected_wf_workspace

end module mod_CorrectedWF2DTypes



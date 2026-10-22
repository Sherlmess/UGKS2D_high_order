module mod_ImplicitMacro
    use flow_data
    use Tools, only: GetPrimary, GetTau
    use mod_ImplicitState, only: rate_w, delta_w, w_tilde, implicit_ncfl_active,macro_pr_passes
    use mod_ImplicitCommon, only: IsPhysicalState, EulerFluxNormal, PseudoTimeStep
    use mod_ImplicitBoundary, only: GetMacroBoundaryDelta, ResolveImplicitBoundaryType, ZhuMovingWallMap
    use mpi_wrapper, only: MPIExchangeDeltaMacro
    use mpi
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    public :: SolveMacroPR2, ComputeMacroDefectNorm, BuildPredictedMacroState
    public :: DumpCornerMacroDefect
    public :: ImplicitSpectralRadius
contains
    real(KREAL) function ImplicitSpectralRadius(i,j,dirc,side) result(gamma_face)
        integer(KINT), intent(in) :: i,j,dirc,side
        integer(KINT) :: ni,nj,boundary_type
        logical :: physical_boundary
        real(KREAL) :: nx_out,ny_out,face_length,projected_distance
        real(KREAL) :: wi(4),wj(4),primi(4),primj(4),wall_prim(4)
        integer(KINT) :: face_number
        real(KREAL) :: rhoi,rhoj,rho_avg,u_avg,v_avg,p_i,p_j,p_avg
        real(KREAL) :: tau_i,tau_j,mu_avg,sound_speed

        call GetFaceGeometry(i,j,dirc,side,ni,nj,nx_out,ny_out,face_length, &
                             boundary_type,physical_boundary)
        wi=ctr(i,j)%conVars
        if (physical_boundary) then
            wj=wi
            if (Implicit_ZhuMovingWall .and. boundary_type==BCconst_wall_isothermal) then
                face_number=2*(dirc-1)+1
                if (side==Rside) face_number=face_number+1
                wall_prim=WallPrimary(face_number)
                wj=ZhuMovingWallMap(wi,wall_prim(2),wall_prim(3))
            endif
            projected_distance=ctr(i,j)%length(dirc)
        else
            wj=ctr(ni,nj)%conVars
            projected_distance=abs((ctr(ni,nj)%x-ctr(i,j)%x)*nx_out+ &
                                   (ctr(ni,nj)%y-ctr(i,j)%y)*ny_out)
        endif
        primi=GetPrimary(wi); primj=GetPrimary(wj)
        rhoi=primi(1); rhoj=primj(1); rho_avg=0.5_KREAL*(rhoi+rhoj)
        u_avg=0.5_KREAL*(primi(2)+primj(2))
        v_avg=0.5_KREAL*(primi(3)+primj(3))
        p_i=0.5_KREAL*primi(1)/primi(4)
        p_j=0.5_KREAL*primj(1)/primj(4)
        p_avg=0.5_KREAL*(p_i+p_j)
        tau_i=GetTau(primi); tau_j=GetTau(primj)
        mu_avg=0.5_KREAL*(tau_i*p_i+tau_j*p_j)
        sound_speed=sqrt(gamma*p_avg/rho_avg)
        gamma_face=abs(u_avg*nx_out+v_avg*ny_out)+sound_speed+ &
            2.0_KREAL*mu_avg/(rho_avg*max(projected_distance, &
            1000.0_KREAL*epsilon(1.0_KREAL)))
    end function ImplicitSpectralRadius

    subroutine SolveMacroPR2(ok,physical_dt,flux_weight)
        logical, intent(out) :: ok
        real(KREAL), intent(in), optional :: physical_dt
        real(KREAL), intent(in), optional :: flux_weight(:,:)
        integer(KINT) :: sweep,i,j
        logical :: point_ok,local_ok,global_ok

        ok=.false.
        call MPIExchangeDeltaMacro(delta_w)
        do sweep=1,Implicit_MacroSweeps
            macro_pr_passes=macro_pr_passes+1
            local_ok=.true.
            do j=1,ny
                do i=1,nx
                    call RelaxMacroPoint(i,j,point_ok,physical_dt,flux_weight)
                    local_ok=local_ok .and. point_ok
                enddo
            enddo
            do j=ny,1,-1
                do i=nx,1,-1
                    call RelaxMacroPoint(i,j,point_ok,physical_dt,flux_weight)
                    local_ok=local_ok .and. point_ok
                enddo
            enddo
            call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
            if (.not.global_ok) return
            call MPIExchangeDeltaMacro(delta_w)
        enddo
        local_ok=all(ieee_is_finite(delta_w(:,1:nx,1:ny)))
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
    end subroutine SolveMacroPR2

    subroutine RelaxMacroPoint(i,j,ok,physical_dt,flux_weight)
        integer(KINT), intent(in) :: i,j
        logical, intent(out) :: ok
        real(KREAL), intent(in), optional :: physical_dt
        real(KREAL), intent(in), optional :: flux_weight(:,:)
        integer(KINT) :: face_number
        real(KREAL) :: delta_t_sum(4),delta_neighbor(4),sum_delta_t(4),sum_gamma_neighbor(4)
        real(KREAL) :: gamma_face,face_length,sum_gamma,dt_pseudo,diagonal,rhs_local(4),coefficient
        logical :: face_ok

        sum_delta_t=0.0_KREAL; sum_gamma_neighbor=0.0_KREAL; sum_gamma=0.0_KREAL
        do face_number=1,4
            call GetMacroFaceTerms(i,j,face_number,delta_t_sum,gamma_face, &
                                   delta_neighbor,face_length,face_ok)
            if (.not.face_ok) then
                ok=.false.; return
            endif
            sum_delta_t=sum_delta_t+face_length*delta_t_sum
            sum_gamma_neighbor=sum_gamma_neighbor+face_length*gamma_face*delta_neighbor
            sum_gamma=sum_gamma+face_length*gamma_face
        enddo
        if (present(physical_dt) .and. present(flux_weight)) then
            coefficient=flux_weight(i,j)
            diagonal=1.0_KREAL/physical_dt+0.5_KREAL*coefficient*sum_gamma/ctr(i,j)%area
        else
            coefficient=1.0_KREAL
            dt_pseudo=PseudoTimeStep(ctr(i,j)%dt_local,implicit_ncfl_active,CFL)
            diagonal=1.0_KREAL/dt_pseudo+0.5_KREAL*sum_gamma/ctr(i,j)%area
        endif
        rhs_local=rate_w(:,i,j)-0.5_KREAL*coefficient*sum_delta_t/ctr(i,j)%area+ &
                  0.5_KREAL*coefficient*sum_gamma_neighbor/ctr(i,j)%area
        ok=ieee_is_finite(diagonal) .and. diagonal>0.0_KREAL .and. all(ieee_is_finite(rhs_local))
        if (ok) delta_w(:,i,j)=rhs_local/diagonal
    end subroutine RelaxMacroPoint

    real(KREAL) function ComputeMacroDefectNorm() result(defect_norm)
        integer(KINT) :: i,j,face_number
        real(KREAL) :: delta_t_sum(4),delta_neighbor(4),approximate_lhs(4),defect(4)
        real(KREAL) :: gamma_face,face_length,dt_pseudo,local_sum,global_sum
        logical :: face_ok,local_ok,global_ok

        local_sum=0.0_KREAL; local_ok=.true.
        do j=1,ny
            do i=1,nx
                dt_pseudo=PseudoTimeStep(ctr(i,j)%dt_local,implicit_ncfl_active,CFL)
                approximate_lhs=delta_w(:,i,j)/dt_pseudo
                do face_number=1,4
                    call GetMacroFaceTerms(i,j,face_number,delta_t_sum,gamma_face, &
                                           delta_neighbor,face_length,face_ok)
                    if (.not.face_ok) then
                        local_ok=.false.
                        cycle
                    endif
                    approximate_lhs=approximate_lhs+0.5_KREAL*face_length* &
                        (delta_t_sum+gamma_face*(delta_w(:,i,j)-delta_neighbor))/ctr(i,j)%area
                enddo
                defect=rate_w(:,i,j)-approximate_lhs
                local_sum=local_sum+sum(defect*defect)
            enddo
        enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(local_sum,global_sum,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        if (global_ok) then
            defect_norm=sqrt(global_sum)
        else
            defect_norm=huge(1.0_KREAL)
        endif
    end function ComputeMacroDefectNorm
    subroutine DumpCornerMacroDefect()
        integer(KINT) :: i,j,face_number,fu
        character(len=80)::filename
        real(KREAL) :: delta_t_sum(4),delta_neighbor(4),approximate_lhs(4),defect(4)
        real(KREAL) :: gamma_face,face_length,dt_pseudo,local_sum,global_sum
        logical :: face_ok,local_ok,global_ok

        write(filename,'(a,i0,a)')'macro_defect_rank',myID,'.bin'
        open(newunit=fu,file=trim(filename),status='new',access='stream',form='unformatted')
        write(fu)nx,ny,Block_XY(myID)%Origin_I,Block_XY(myID)%Origin_J
        local_sum=0.0_KREAL; local_ok=.true.
        do j=1,ny
            do i=1,nx
                dt_pseudo=PseudoTimeStep(ctr(i,j)%dt_local,implicit_ncfl_active,CFL)
                approximate_lhs=delta_w(:,i,j)/dt_pseudo
                do face_number=1,4
                    call GetMacroFaceTerms(i,j,face_number,delta_t_sum,gamma_face, &
                                           delta_neighbor,face_length,face_ok)
                    if (.not.face_ok) then
                        local_ok=.false.
                        cycle
                    endif
                    approximate_lhs=approximate_lhs+0.5_KREAL*face_length* &
                        (delta_t_sum+gamma_face*(delta_w(:,i,j)-delta_neighbor))/ctr(i,j)%area
                enddo
                defect=rate_w(:,i,j)-approximate_lhs
                write(fu)dt_pseudo,rate_w(:,i,j),delta_w(:,i,j),w_tilde(:,i,j),defect
                if(.not.local_ok)error stop 'diagnostic macro face failure' 
            enddo
        enddo
        close(fu)
    end subroutine DumpCornerMacroDefect

    subroutine BuildPredictedMacroState(alpha,ok)
        real(KREAL), intent(out) :: alpha
        logical, intent(out) :: ok
        integer(KINT) :: i,j
        logical :: local_ok,global_ok

        alpha=1.0_KREAL; ok=.false.
        do
            local_ok=.true.
            do j=1,ny
                do i=1,nx
                    if (.not.IsPhysicalState(ctr(i,j)%conVars+alpha*delta_w(:,i,j))) local_ok=.false.
                enddo
            enddo
            call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
            if (global_ok) exit
            alpha=0.5_KREAL*alpha
            if (alpha<Implicit_MinAlpha) return
        enddo
        delta_w=alpha*delta_w
        do j=1,ny
            do i=1,nx
                w_tilde(:,i,j)=ctr(i,j)%conVars+delta_w(:,i,j)
            enddo
        enddo
        ok=.true.
    end subroutine BuildPredictedMacroState

    subroutine GetMacroFaceTerms(i,j,face_number,delta_t_sum,gamma_face, &
                                 delta_neighbor,face_length,ok)
        integer(KINT), intent(in) :: i,j,face_number
        real(KREAL), intent(out) :: delta_t_sum(4),gamma_face,delta_neighbor(4),face_length
        logical, intent(out) :: ok
        integer(KINT) :: dirc,side,ni,nj,boundary_type
        logical :: physical_boundary
        real(KREAL) :: nx_out,ny_out,wi(4),wj(4),flux_old(4),flux_new(4),delta_t_i(4),delta_t_j(4),wall_prim(4)

        call FaceNumberToDirection(face_number,dirc,side)
        call GetFaceGeometry(i,j,dirc,side,ni,nj,nx_out,ny_out,face_length, &
                             boundary_type,physical_boundary)
        wi=ctr(i,j)%conVars
        if (physical_boundary) then
            wj=wi
            if (Implicit_ZhuMovingWall .and. boundary_type==BCconst_wall_isothermal) then
                wall_prim=WallPrimary(face_number)
                wj=ZhuMovingWallMap(wi,wall_prim(2),wall_prim(3))
            endif
            call GetMacroBoundaryDelta(boundary_type,nx_out,ny_out,delta_w(:,i,j),delta_neighbor,ok,face_number)
            if (.not.ok) return
        else
            wj=ctr(ni,nj)%conVars
            delta_neighbor=delta_w(:,ni,nj)
        endif
        call EulerFluxNormal(wi+delta_w(:,i,j),nx_out,ny_out,flux_new)
        call EulerFluxNormal(wi,nx_out,ny_out,flux_old)
        delta_t_i=flux_new-flux_old
        call EulerFluxNormal(wj+delta_neighbor,nx_out,ny_out,flux_new)
        call EulerFluxNormal(wj,nx_out,ny_out,flux_old)
        delta_t_j=flux_new-flux_old
        gamma_face=ImplicitSpectralRadius(i,j,dirc,side)
        delta_t_sum=delta_t_i+delta_t_j
        ok=all(ieee_is_finite(delta_t_sum)) .and. ieee_is_finite(gamma_face) .and. gamma_face>0.0_KREAL
    end subroutine GetMacroFaceTerms

    subroutine FaceNumberToDirection(face_number,dirc,side)
        integer(KINT), intent(in) :: face_number
        integer(KINT), intent(out) :: dirc,side

        select case(face_number)
        case(1); dirc=IDIRC; side=Lside
        case(2); dirc=IDIRC; side=Rside
        case(3); dirc=JDIRC; side=Lside
        case default; dirc=JDIRC; side=Rside
        end select
    end subroutine FaceNumberToDirection

    subroutine GetFaceGeometry(i,j,dirc,side,ni,nj,nx_out,ny_out,face_length, &
                               boundary_type,physical_boundary)
        integer(KINT), intent(in) :: i,j,dirc,side
        integer(KINT), intent(out) :: ni,nj,boundary_type
        real(KREAL), intent(out) :: nx_out,ny_out,face_length
        logical, intent(out) :: physical_boundary

        ni=i; nj=j; physical_boundary=.false.; boundary_type=0
        if (dirc==IDIRC) then
            ni=i+side
            if (side==Lside) then
                nx_out=-vface(i,j)%cosx; ny_out=-vface(i,j)%cosy; face_length=vface(i,j)%length
                physical_boundary=(i==1 .and. ID_X_minus==MPI_PROC_NULL); boundary_type=BC_Type(1)
            else
                nx_out=vface(i+1,j)%cosx; ny_out=vface(i+1,j)%cosy; face_length=vface(i+1,j)%length
                physical_boundary=(i==nx .and. ID_X_plus==MPI_PROC_NULL); boundary_type=BC_Type(2)
            endif
        else
            nj=j+side
            if (side==Lside) then
                nx_out=-hface(i,j)%cosx; ny_out=-hface(i,j)%cosy; face_length=hface(i,j)%length
                physical_boundary=(j==1 .and. ID_Y_minus==MPI_PROC_NULL); boundary_type=BC_Type(3)
            else
                nx_out=hface(i,j+1)%cosx; ny_out=hface(i,j+1)%cosy; face_length=hface(i,j+1)%length
                physical_boundary=(j==ny .and. ID_Y_plus==MPI_PROC_NULL); boundary_type=BC_Type(4)
            endif
        endif
        if (physical_boundary) boundary_type=ResolveImplicitBoundaryType(i,j,dirc,side,boundary_type)
    end subroutine GetFaceGeometry
end module mod_ImplicitMacro

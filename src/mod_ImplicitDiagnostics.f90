module mod_ImplicitDiagnostics
    use flow_data
    use Tools, only: GetPrimary, GetTau, GetHeatFlux_Vec, DiscreteMaxwell_Vec, ShakhovPart_Vec, &
        RelativeOrAbsoluteResidual, IsNearZeroResidualScale
    use mod_ImplicitState
    use mod_ImplicitMicro, only: WeightedNegativeFraction
    use mod_ImplicitCommon, only: IsPhysicalState, ProjectToConservedMoments
    use mpi
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    private
    public :: InitializeImplicitDiagnostics, UpdateImplicitResidualDiagnostics
    public :: ImplicitConverged, WriteImplicitResidual, AbortSteadyImplicit
    public :: CollectImplicitTimingMaxima
    public :: AcceptedStateAdmissibilityDiagnostics, CollectAcceptedStateAdmissibilityDiagnostics
    public :: SelectGlobalAcceptedStateAdmissibilityMaxima, WriteAcceptedStateAdmissibilityDiagnostics
    integer(KINT), parameter, public :: ADMISSIBILITY_METRIC_COUNT=5

    type :: AcceptedStateAdmissibilityDiagnostics
        real(KREAL) :: max_value(ADMISSIBILITY_METRIC_COUNT)=0.0_KREAL
        integer(KINT) :: global_i(ADMISSIBILITY_METRIC_COUNT)=0
        integer(KINT) :: global_j(ADMISSIBILITY_METRIC_COUNT)=0
        integer(KINT) :: owner_rank(ADMISSIBILITY_METRIC_COUNT)=-1
    end type AcceptedStateAdmissibilityDiagnostics

    real(KREAL), save :: reference_h=0.0_KREAL, reference_b=0.0_KREAL
    logical, save :: references_initialized=.false.
    logical, save :: residual_output_initialized=.false.
contains
    subroutine InitializeImplicitDiagnostics()
        reference_h=0.0_KREAL; reference_b=0.0_KREAL
        references_initialized=.false.
        implicit_macro_residual=0.0_KREAL
        implicit_macro_relative=0.0_KREAL; implicit_macro_absolute=0.0_KREAL
        implicit_h_norm=0.0_KREAL; implicit_b_norm=0.0_KREAL
        implicit_h_relative=0.0_KREAL; implicit_b_relative=0.0_KREAL
        implicit_h_absolute=0.0_KREAL; implicit_b_absolute=0.0_KREAL
        implicit_moment_error=0.0_KREAL; implicit_mass_defect=0.0_KREAL
        implicit_pr_seconds=0.0_KREAL; implicit_lusgs_seconds=0.0_KREAL
        implicit_exchange_seconds=0.0_KREAL; implicit_evolution_seconds=0.0_KREAL
        implicit_micro_forward_seconds=0.0_KREAL
        implicit_micro_backward_seconds=0.0_KREAL
        implicit_micro_initial_exchange_seconds=0.0_KREAL
        implicit_micro_forward_backward_exchange_seconds=0.0_KREAL
        implicit_micro_post_backward_exchange_seconds=0.0_KREAL
        implicit_micro_exchange_calls=0
        implicit_total_seconds=0.0_KREAL
        implicit_retry_count=0; implicit_retry_total=0
        residual_output_initialized=.false.
    end subroutine InitializeImplicitDiagnostics

    subroutine CollectAcceptedStateAdmissibilityDiagnostics(diagnostics)
        type(AcceptedStateAdmissibilityDiagnostics), intent(out) :: diagnostics
        integer(KINT) :: i,j,m,cell_global_i,cell_global_j
        real(KREAL) :: local_values(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT) :: local_global_i(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT) :: local_global_j(ADMISSIBILITY_METRIC_COUNT)
        real(KREAL), allocatable :: maxwell_h(:),maxwell_b(:),shakhov_h(:),shakhov_b(:)
        real(KREAL), allocatable :: equilibrium_h(:),equilibrium_b(:)
        real(KREAL) :: cell_values(ADMISSIBILITY_METRIC_COUNT)
        logical :: cell_ok

        diagnostics%max_value=0.0_KREAL
        diagnostics%global_i=0
        diagnostics%global_j=0
        diagnostics%owner_rank=-1
        local_values=0.0_KREAL
        local_global_i=0
        local_global_j=0

        if (nx<=0 .or. ny<=0 .or. Nvel<=0 .or. .not.allocated(ctr) .or. &
            .not.allocated(wVel) .or. .not.allocated(uVel) .or. .not.allocated(vVel)) then
            call SelectGlobalAcceptedStateAdmissibilityMaxima(local_values,local_global_i,local_global_j, &
                diagnostics%max_value,diagnostics%global_i,diagnostics%global_j,diagnostics%owner_rank)
            return
        endif
        if (lbound(ctr,1)>1 .or. ubound(ctr,1)<nx .or. lbound(ctr,2)>1 .or. ubound(ctr,2)<ny .or. &
            size(wVel)/=Nvel .or. size(uVel)/=Nvel .or. size(vVel)/=Nvel) then
            call SelectGlobalAcceptedStateAdmissibilityMaxima(local_values,local_global_i,local_global_j, &
                diagnostics%max_value,diagnostics%global_i,diagnostics%global_j,diagnostics%owner_rank)
            return
        endif

        allocate(maxwell_h(Nvel),maxwell_b(Nvel),shakhov_h(Nvel),shakhov_b(Nvel))
        allocate(equilibrium_h(Nvel),equilibrium_b(Nvel))
        do j=1,ny
            do i=1,nx
                if (.not.associated(ctr(i,j)%h_vec) .or. .not.associated(ctr(i,j)%b_vec)) cycle
                call EvaluateAcceptedStateAdmissibilityCell(ctr(i,j)%h_vec,ctr(i,j)%b_vec,ctr(i,j)%conVars, &
                    maxwell_h,maxwell_b,shakhov_h,shakhov_b,equilibrium_h,equilibrium_b,cell_values,cell_ok)
                if (.not.cell_ok) cycle
                if (allocated(Block_XY)) then
                    if (IP_XY>=lbound(Block_XY,1) .and. IP_XY<=ubound(Block_XY,1)) then
                        cell_global_i=Block_XY(IP_XY)%Origin_I+i-1
                        cell_global_j=Block_XY(IP_XY)%Origin_J+j-1
                    else
                        cell_global_i=i
                        cell_global_j=j
                    endif
                else
                    cell_global_i=i
                    cell_global_j=j
                endif
                do m=1,ADMISSIBILITY_METRIC_COUNT
                    if (local_global_i(m)==0 .or. cell_values(m)>local_values(m)) then
                        local_values(m)=cell_values(m)
                        local_global_i(m)=cell_global_i
                        local_global_j(m)=cell_global_j
                    endif
                enddo
            enddo
        enddo
        deallocate(maxwell_h,maxwell_b,shakhov_h,shakhov_b,equilibrium_h,equilibrium_b)
        call SelectGlobalAcceptedStateAdmissibilityMaxima(local_values,local_global_i,local_global_j, &
            diagnostics%max_value,diagnostics%global_i,diagnostics%global_j,diagnostics%owner_rank)
    end subroutine CollectAcceptedStateAdmissibilityDiagnostics

    subroutine EvaluateAcceptedStateAdmissibilityCell(h_vec,b_vec,con_vars,maxwell_h,maxwell_b, &
        shakhov_h,shakhov_b,equilibrium_h,equilibrium_b,metric_values,ok)
        real(KREAL), intent(in) :: h_vec(:),b_vec(:),con_vars(4)
        real(KREAL), intent(out) :: maxwell_h(:),maxwell_b(:),shakhov_h(:),shakhov_b(:)
        real(KREAL), intent(out) :: equilibrium_h(:),equilibrium_b(:)
        real(KREAL), intent(out) :: metric_values(ADMISSIBILITY_METRIC_COUNT)
        logical, intent(out) :: ok
        real(KREAL) :: prim(4),qf(2)

        metric_values=0.0_KREAL
        ok=.false.
        if (size(h_vec)/=Nvel .or. size(b_vec)/=Nvel .or. size(maxwell_h)/=Nvel .or. &
            size(maxwell_b)/=Nvel .or. size(shakhov_h)/=Nvel .or. size(shakhov_b)/=Nvel .or. &
            size(equilibrium_h)/=Nvel .or. size(equilibrium_b)/=Nvel) return
        if (.not.IsPhysicalState(con_vars) .or. .not.all(ieee_is_finite(h_vec)) .or. &
            .not.all(ieee_is_finite(b_vec))) return
        prim=GetPrimary(con_vars)
        qf=GetHeatFlux_Vec(h_vec,b_vec,prim)
        if (.not.all(ieee_is_finite(prim)) .or. .not.all(ieee_is_finite(qf))) return
        call DiscreteMaxwell_Vec(maxwell_h,maxwell_b,prim)
        call ShakhovPart_Vec(maxwell_h,maxwell_b,qf,prim,shakhov_h,shakhov_b)
        equilibrium_h=maxwell_h+shakhov_h
        equilibrium_b=maxwell_b+shakhov_b
        if (.not.all(ieee_is_finite(maxwell_h)) .or. .not.all(ieee_is_finite(maxwell_b)) .or. &
            .not.all(ieee_is_finite(shakhov_h)) .or. .not.all(ieee_is_finite(shakhov_b)) .or. &
            .not.all(ieee_is_finite(equilibrium_h)) .or. .not.all(ieee_is_finite(equilibrium_b))) return
        metric_values(1)=WeightedNegativeFraction(h_vec,wVel)
        metric_values(2)=WeightedNegativeFraction(b_vec,wVel)
        metric_values(3)=WeightedNegativeFraction(equilibrium_h,wVel)
        metric_values(4)=WeightedNegativeFraction(equilibrium_b,wVel)
        metric_values(5)=sqrt(qf(1)*qf(1)+qf(2)*qf(2))
        ok=all(ieee_is_finite(metric_values))
    end subroutine EvaluateAcceptedStateAdmissibilityCell

    subroutine SelectGlobalAcceptedStateAdmissibilityMaxima(local_values,local_global_i,local_global_j, &
        global_values,global_i,global_j,owner_rank)
        real(KREAL), intent(in) :: local_values(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT), intent(in) :: local_global_i(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT), intent(in) :: local_global_j(ADMISSIBILITY_METRIC_COUNT)
        real(KREAL), intent(out) :: global_values(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT), intent(out) :: global_i(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT), intent(out) :: global_j(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT), intent(out) :: owner_rank(ADMISSIBILITY_METRIC_COUNT)
        integer(KINT) :: m,candidate_index,selected_index,candidate_owner,selected_owner

        call MPI_Allreduce(local_values,global_values,ADMISSIBILITY_METRIC_COUNT, &
            MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        do m=1,ADMISSIBILITY_METRIC_COUNT
            candidate_index=huge(1_KINT)
            if (Nx_all>0 .and. Ny_all>0 .and. local_global_i(m)>=1 .and. &
                local_global_i(m)<=Nx_all .and. local_global_j(m)>=1 .and. &
                local_global_j(m)<=Ny_all .and. local_values(m)==global_values(m)) then
                candidate_index=(local_global_j(m)-1)*Nx_all+local_global_i(m)
            endif
            call MPI_Allreduce(candidate_index,selected_index,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
            if (selected_index/=huge(1_KINT) .and. selected_index>0) then
                global_i(m)=mod(selected_index-1,Nx_all)+1
                global_j(m)=(selected_index-1)/Nx_all+1
                candidate_owner=huge(1_KINT)
                if (candidate_index==selected_index) candidate_owner=myID
                call MPI_Allreduce(candidate_owner,selected_owner,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
                if (selected_owner/=huge(1_KINT)) then
                    owner_rank(m)=selected_owner
                else
                    owner_rank(m)=-1
                endif
            else
                global_i(m)=0
                global_j(m)=0
                owner_rank(m)=-1
            endif
        enddo
    end subroutine SelectGlobalAcceptedStateAdmissibilityMaxima

    subroutine WriteAcceptedStateAdmissibilityDiagnostics(mode_label)
        character(len=*), intent(in) :: mode_label
        type(AcceptedStateAdmissibilityDiagnostics) :: diagnostics

        call CollectAcceptedStateAdmissibilityDiagnostics(diagnostics)
        if (myID/=0) return
        write(output_unit,'(A)') trim(mode_label)//' admissibility diagnostics:'
        write(output_unit,'(1P,A,E16.8,A,I0,A,I0,A,I0)') &
            'current_h_negative_fraction_max=',diagnostics%max_value(1), &
            ' global_i=',diagnostics%global_i(1),' global_j=',diagnostics%global_j(1), &
            ' owner_rank=',diagnostics%owner_rank(1)
        write(output_unit,'(1P,A,E16.8,A,I0,A,I0,A,I0)') &
            'current_b_negative_fraction_max=',diagnostics%max_value(2), &
            ' global_i=',diagnostics%global_i(2),' global_j=',diagnostics%global_j(2), &
            ' owner_rank=',diagnostics%owner_rank(2)
        write(output_unit,'(1P,A,E16.8,A,I0,A,I0,A,I0)') &
            'shakhov_equilibrium_h_negative_fraction_max=',diagnostics%max_value(3), &
            ' global_i=',diagnostics%global_i(3),' global_j=',diagnostics%global_j(3), &
            ' owner_rank=',diagnostics%owner_rank(3)
        write(output_unit,'(1P,A,E16.8,A,I0,A,I0,A,I0)') &
            'shakhov_equilibrium_b_negative_fraction_max=',diagnostics%max_value(4), &
            ' global_i=',diagnostics%global_i(4),' global_j=',diagnostics%global_j(4), &
            ' owner_rank=',diagnostics%owner_rank(4)
        write(output_unit,'(1P,A,E16.8,A,I0,A,I0,A,I0)') &
            'heat_flux_magnitude_max=',diagnostics%max_value(5), &
            ' global_i=',diagnostics%global_i(5),' global_j=',diagnostics%global_j(5), &
            ' owner_rank=',diagnostics%owner_rank(5)
        flush(output_unit)
    end subroutine WriteAcceptedStateAdmissibilityDiagnostics

    subroutine UpdateImplicitResidualDiagnostics()
        integer(KINT) :: i,j,k
        real(KREAL) :: macro_num_local(4),macro_num_global(4)
        real(KREAL) :: macro_den_local(4),macro_den_global(4)
        real(KREAL) :: macro_weight_local,macro_weight_global
        real(KREAL) :: micro_sums_local(4),micro_sums_global(4)
        real(KREAL) :: dh,db,mass_local,mass_global,h_raw,b_raw
        real(KREAL) :: prim(4),tau_value,q_current(2),ordinary_macro_scale
        real(KREAL) :: near_zero_micro_reference
        real(KREAL) :: physical_equilibrium_h(Nvel),physical_equilibrium_b(Nvel)
        real(KREAL) :: h_plus(Nvel),b_plus(Nvel),equilibrium_reference_h(Nvel),equilibrium_reference_b(Nvel)
        logical :: projection_ok,local_ok,global_ok

        macro_num_local=0.0_KREAL; macro_den_local=0.0_KREAL; macro_weight_local=0.0_KREAL
        micro_sums_local=0.0_KREAL; mass_local=0.0_KREAL; local_ok=.true.
        do j=1,ny; do i=1,nx
            macro_num_local=macro_num_local+(ctr(i,j)%area*rate_w(:,i,j))**2
            macro_den_local=macro_den_local+ctr(i,j)%area* &
                abs(ctr(i,j)%conVars+ctr(i,j)%dt_local*rate_w(:,i,j))
            macro_weight_local=macro_weight_local+ctr(i,j)%area**2
            prim=GetPrimary(ctr(i,j)%conVars)
            tau_value=GetTau(prim)
            q_current=GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim)
            call DiscreteMaxwell_Vec(physical_equilibrium_h,physical_equilibrium_b,prim)
            call ShakhovPart_Vec(physical_equilibrium_h,physical_equilibrium_b,q_current,prim,h_plus,b_plus)
            physical_equilibrium_h=physical_equilibrium_h+h_plus
            physical_equilibrium_b=physical_equilibrium_b+b_plus
            equilibrium_reference_h=physical_equilibrium_h
            equilibrium_reference_b=physical_equilibrium_b
            call ProjectToConservedMoments(physical_equilibrium_h,physical_equilibrium_b, &
                ctr(i,j)%conVars,uVel,vVel,wVel,equilibrium_reference_h,equilibrium_reference_b,projection_ok)
            local_ok=local_ok .and. projection_ok .and. tau_value>0.0_KREAL .and. &
                all(ieee_is_finite(physical_equilibrium_h)) .and. &
                all(ieee_is_finite(physical_equilibrium_b))
            do k=1,Nvel
                dh=stream_h(k,i,j)+(physical_equilibrium_h(k)-ctr(i,j)%h_vec(k))/tau_value
                db=stream_b(k,i,j)+(physical_equilibrium_b(k)-ctr(i,j)%b_vec(k))/tau_value
                micro_sums_local(1)=micro_sums_local(1)+ctr(i,j)%area*wVel(k)*dh*dh
                micro_sums_local(2)=micro_sums_local(2)+ctr(i,j)%area*wVel(k)*db*db
            enddo
            micro_sums_local(3)=micro_sums_local(3)+ctr(i,j)%area*sum(wVel* &
                max(ctr(i,j)%h_vec**2,physical_equilibrium_h**2))
            micro_sums_local(4)=micro_sums_local(4)+ctr(i,j)%area*sum(wVel* &
                max(ctr(i,j)%b_vec**2,physical_equilibrium_b**2))
            mass_local=mass_local+rate_w(1,i,j)*ctr(i,j)%area
        enddo; enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global_ok) call AbortSteadyImplicit('physical microscopic residual evaluation failed')
        call MPI_Allreduce(macro_num_local,macro_num_global,4,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(macro_den_local,macro_den_global,4,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(macro_weight_local,macro_weight_global,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(micro_sums_local,micro_sums_global,4,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(mass_local,mass_global,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        implicit_macro_relative=sqrt(real(Nx_all*Ny_all,KREAL)*macro_num_global)/(macro_den_global+SMV)
        implicit_macro_absolute=sqrt(macro_num_global/max(macro_weight_global,SMV))
        ordinary_macro_scale=max(macro_den_global(1),macro_den_global(4),SMV)
        do k=1,4
            if (IsNearZeroResidualScale(macro_den_global(k),ordinary_macro_scale)) then
                implicit_macro_residual(k)=implicit_macro_absolute(k)
            else
                implicit_macro_residual(k)=implicit_macro_relative(k)
            endif
        enddo
        h_raw=sqrt(micro_sums_global(1)/max(micro_sums_global(3),SMV))
        b_raw=sqrt(micro_sums_global(2)/max(micro_sums_global(4),SMV))
        implicit_mass_defect=abs(mass_global)
        if (.not.references_initialized) then
            reference_h=h_raw
            reference_b=b_raw
            references_initialized=.true.
        endif
        implicit_h_absolute=h_raw; implicit_b_absolute=b_raw
        if (reference_h>tiny(1.0_KREAL)) then
            implicit_h_relative=h_raw/reference_h
        else
            implicit_h_relative=h_raw
        endif
        if (reference_b>tiny(1.0_KREAL)) then
            implicit_b_relative=b_raw/reference_b
        else
            implicit_b_relative=b_raw
        endif
        near_zero_micro_reference=max(100.0_KREAL*Implicit_ResidualTol,1000.0_KREAL*epsilon(1.0_KREAL))
        implicit_h_norm=RelativeOrAbsoluteResidual(h_raw,reference_h,near_zero_micro_reference)
        implicit_b_norm=RelativeOrAbsoluteResidual(b_raw,reference_b,near_zero_micro_reference)
    end subroutine UpdateImplicitResidualDiagnostics

    logical function ImplicitConverged()
        ImplicitConverged=all(implicit_macro_residual<=Implicit_ResidualTol) .and. &
            implicit_h_norm<=Implicit_ResidualTol .and. &
            implicit_b_norm<=Implicit_ResidualTol .and. &
            implicit_moment_error<=Implicit_MomentTol
    end function ImplicitConverged

    subroutine CollectImplicitTimingMaxima(wall_time,timing_global)
        real(KREAL), intent(in) :: wall_time
        real(KREAL), intent(out) :: timing_global(8)
        real(KREAL) :: timing_local(8)

        timing_local=[wall_time,implicit_pr_seconds,implicit_lusgs_seconds, &
            implicit_exchange_seconds,implicit_evolution_seconds,implicit_micro_forward_seconds, &
            implicit_micro_backward_seconds,implicit_micro_forward_backward_exchange_seconds]
        call MPI_Allreduce(timing_local,timing_global,8,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    end subroutine CollectImplicitTimingMaxima

    subroutine WriteImplicitResidual(wall_time)
        real(KREAL), intent(in) :: wall_time
        real(KREAL) :: dt_local_min,dt_global_min,timing_global(8)
        integer(KINT) :: i,j

        dt_local_min=huge(1.0_KREAL)
        do j=1,ny; do i=1,nx
            dt_local_min=min(dt_local_min,ctr(i,j)%dt_local)
        enddo; enddo
        call MPI_Allreduce(dt_local_min,dt_global_min,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
        call CollectImplicitTimingMaxima(wall_time,timing_global)
        if (myID==0) then
            if (.not.residual_output_initialized) then
                open(104,file='ResidualImplicit.plt',status='replace',action='write')
                write(104,'(A)') 'VARIABLES = iter dtPhysicalMin nCFL alpha rhoResidual rhoUResidual ' // &
                    'rhoVResidual energyResidual hResidual bResidual momentError massDefect retries wallTime ' // &
                    'prSeconds lusgsSeconds correctionExchangeSeconds evolutionSeconds retryTotal ' // &
                    'rhoRelative rhoURelative rhoVRelative energyRelative rhoAbsolute rhoUAbsolute ' // &
                    'rhoVAbsolute energyAbsolute hRelative bRelative hAbsolute bAbsolute ' // &
                    'microForwardSeconds microBackwardSeconds microForwardBackwardExchangeSeconds'
                residual_output_initialized=.true.
            else
                open(104,file='ResidualImplicit.plt',status='unknown',action='write',position='append')
            endif
            write(104,'(I12,1P,11E18.9,I8,5E18.9,I12,15E18.9)') iter,dt_global_min,implicit_ncfl_active,implicit_alpha, &
                implicit_macro_residual,implicit_h_norm,implicit_b_norm,implicit_moment_error,implicit_mass_defect, &
                implicit_retry_count,timing_global(1:5),implicit_retry_total,implicit_macro_relative, &
                implicit_macro_absolute,implicit_h_relative,implicit_b_relative, &
                implicit_h_absolute,implicit_b_absolute,timing_global(6:8)
            close(104)
            write(output_unit,'(A,I0,1P,11(A,E12.4),2(A,I0))') &
                'Implicit iter=',iter,' nCFL=',implicit_ncfl_active,' alpha=',implicit_alpha, &
                ' Rrho=',implicit_macro_residual(1),' RrhoU=',implicit_macro_residual(2), &
                ' RrhoV=',implicit_macro_residual(3),' RE=',implicit_macro_residual(4), &
                ' Rh=',implicit_h_norm,' Rb=',implicit_b_norm, &
                ' moment=',implicit_moment_error,' mass=',implicit_mass_defect,' wall=',wall_time, &
                ' retries=',implicit_retry_count,' retryTotal=',implicit_retry_total
            flush(output_unit)
        endif
    end subroutine WriteImplicitResidual

    subroutine AbortSteadyImplicit(stage)
        character(len=*), intent(in) :: stage
        character(len=32) :: reason_label

        if (implicit_failure%valid) then
            if (myID==0) then
                reason_label=ImplicitFailureReasonLabel(implicit_failure%reason)
                write(output_unit,'(A,I0,A,A,A,I0)') &
                    'Steady implicit failure rank=',myID,' stage=',trim(stage),' iter=',iter
                write(output_unit,'(A,I0,A,I0,A,I0,A,I0,A,I0)') &
                    'owner_rank=',implicit_failure%owner_rank, &
                    ' global_i=',implicit_failure%global_i, &
                    ' global_j=',implicit_failure%global_j, &
                    ' local_i=',implicit_failure%local_i, &
                    ' local_j=',implicit_failure%local_j
                write(output_unit,'(A,A,A,I0)') &
                    'reason=',trim(reason_label),' velocity_k=',implicit_failure%velocity_k
                if (implicit_failure%velocity_k>0) then
                    write(output_unit,'(1P,2(A,E16.8))') &
                        'u=',implicit_failure%velocity_u,' v=',implicit_failure%velocity_v
                else
                    write(output_unit,'(A)') 'u=not_applicable v=not_applicable'
                endif
                write(output_unit,'(1P,4(A,E16.8))') &
                    'nCFL_attempt=',implicit_failure%ncfl_attempt, &
                    ' alpha_tested=',implicit_failure%alpha_tested, &
                    ' alpha_next=',implicit_failure%alpha_next, &
                    ' minAlpha=',Implicit_MinAlpha
                if (implicit_failure%reason==IMPLICIT_FAILURE_H_NEGATIVE_FRACTION .or. &
                    implicit_failure%reason==IMPLICIT_FAILURE_B_NEGATIVE_FRACTION) then
                    write(output_unit,'(1P,3(A,E16.8))') &
                        'current_fraction=',implicit_failure%current_fraction, &
                        ' candidate_fraction=',implicit_failure%candidate_fraction, &
                        ' allowed_fraction=',implicit_failure%allowed_fraction
                else
                    write(output_unit,'(A)') &
                        'current_fraction=not_applicable candidate_fraction=not_applicable allowed_fraction=not_applicable'
                endif
                if (implicit_failure%reason==IMPLICIT_FAILURE_H_NONFINITE .or. &
                    implicit_failure%reason==IMPLICIT_FAILURE_B_NONFINITE) then
                    write(output_unit,'(A)') &
                        'candidate_rho=not_applicable candidate_rhoU=not_applicable ' // &
                        'candidate_rhoV=not_applicable candidate_E=not_applicable'
                    write(output_unit,'(A)') 'candidate_internal_energy=not_applicable'
                else
                    write(output_unit,'(1P,4(A,E16.8))') &
                        'candidate_rho=',implicit_failure%candidate_moment(1), &
                        ' candidate_rhoU=',implicit_failure%candidate_moment(2), &
                        ' candidate_rhoV=',implicit_failure%candidate_moment(3), &
                        ' candidate_E=',implicit_failure%candidate_moment(4)
                    if (all(ieee_is_finite(implicit_failure%candidate_moment)) .and. &
                        implicit_failure%candidate_moment(1)>0.0_KREAL .and. &
                        ieee_is_finite(implicit_failure%candidate_internal_energy)) then
                        write(output_unit,'(1P,A,E16.8)') 'candidate_internal_energy=', &
                            implicit_failure%candidate_internal_energy
                    else
                        write(output_unit,'(A)') 'candidate_internal_energy=not_applicable'
                    endif
                endif
                if (implicit_failure%micro_breakdown_valid) then
                    if (implicit_failure%reason==IMPLICIT_FAILURE_H_NEGATIVE_FRACTION) then
                        write(output_unit,'(A)') 'micro_component=h'
                    else
                        write(output_unit,'(A)') 'micro_component=b'
                    endif
                    write(output_unit,'(1P,A,E16.8)') 'micro_current=',implicit_failure%micro_current
                    write(output_unit,'(1P,A,E16.8)') 'micro_delta=',implicit_failure%micro_delta
                    write(output_unit,'(1P,A,E16.8)') 'micro_candidate=',implicit_failure%micro_candidate
                    write(output_unit,'(1P,A,E16.8)') 'equilibrium_value=',implicit_failure%equilibrium_value
                    write(output_unit,'(1P,A,E16.8)') 'stream_rhs=',implicit_failure%stream_rhs
                    write(output_unit,'(1P,A,E16.8)') 'collision_rhs=',implicit_failure%collision_rhs
                    write(output_unit,'(1P,A,E16.8)') 'neighbor_rhs=',implicit_failure%neighbor_rhs
                    write(output_unit,'(1P,A,E16.8)') 'rhs_total=',implicit_failure%rhs_total
                    write(output_unit,'(1P,A,E16.8)') 'dt_pseudo=',implicit_failure%dt_pseudo
                    write(output_unit,'(1P,A,E16.8)') 'tau_value=',implicit_failure%tau_value
                    write(output_unit,'(1P,A,E16.8)') 'diagonal_total=',implicit_failure%diagonal_total
                    write(output_unit,'(1P,A,E16.8)') 'post_sweep_delta=',implicit_failure%post_sweep_delta
                else
                    write(output_unit,'(A)') 'micro_breakdown=not_applicable'
                endif
                write(output_unit,'(A,A)') 'equilibrium_breakdown_stage=', &
                    trim(ImplicitEquilibriumBreakdownStageLabel(implicit_failure%equilibrium_breakdown_stage))
                if (implicit_failure%equilibrium_breakdown_stage==IMPLICIT_EQUILIBRIUM_BREAKDOWN_SUCCESS .or. &
                    implicit_failure%equilibrium_breakdown_stage== &
                    IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH) then
                    if (implicit_failure%reason==IMPLICIT_FAILURE_H_NEGATIVE_FRACTION) then
                        write(output_unit,'(A)') 'equilibrium_component=h'
                    else
                        write(output_unit,'(A)') 'equilibrium_component=b'
                    endif
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_maxwell_value=',implicit_failure%equilibrium_maxwell_value
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_shakhov_correction=',implicit_failure%equilibrium_shakhov_correction
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_preprojection_value=',implicit_failure%equilibrium_preprojection_value
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_projection_correction=',implicit_failure%equilibrium_projection_correction
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_final_value=',implicit_failure%equilibrium_final_value
                    write(output_unit,'(1P,A,E16.8)') 'q_frozen_x=',implicit_failure%q_frozen_x
                    write(output_unit,'(1P,A,E16.8)') 'q_frozen_y=',implicit_failure%q_frozen_y
                    write(output_unit,'(1P,A,E16.8)') 'prim_tilde_rho=',implicit_failure%prim_tilde_rho
                    write(output_unit,'(1P,A,E16.8)') 'prim_tilde_u=',implicit_failure%prim_tilde_u
                    write(output_unit,'(1P,A,E16.8)') 'prim_tilde_v=',implicit_failure%prim_tilde_v
                    write(output_unit,'(1P,A,E16.8)') &
                        'prim_tilde_temperature=',implicit_failure%prim_tilde_temperature
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_maxwell_negative_fraction=', &
                        implicit_failure%equilibrium_maxwell_negative_fraction
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_preprojection_negative_fraction=', &
                        implicit_failure%equilibrium_preprojection_negative_fraction
                    write(output_unit,'(1P,A,E16.8)') &
                        'equilibrium_final_negative_fraction=', &
                        implicit_failure%equilibrium_final_negative_fraction
                    if (implicit_failure%equilibrium_breakdown_stage== &
                        IMPLICIT_EQUILIBRIUM_BREAKDOWN_STORED_RECONSTRUCTION_MISMATCH) then
                        write(output_unit,'(1P,A,E16.8)') 'reconstruction_max_abs_h=', &
                            implicit_failure%reconstruction_max_abs_h
                        write(output_unit,'(1P,A,E16.8)') 'reconstruction_max_abs_b=', &
                            implicit_failure%reconstruction_max_abs_b
                        write(output_unit,'(1P,A,E16.8)') 'reconstruction_max_rel_h=', &
                            implicit_failure%reconstruction_max_rel_h
                        write(output_unit,'(1P,A,E16.8)') 'reconstruction_max_rel_b=', &
                            implicit_failure%reconstruction_max_rel_b
                        write(output_unit,'(A,I0)') 'reconstruction_max_abs_h_k=', &
                            implicit_failure%reconstruction_max_abs_h_k
                        write(output_unit,'(A,I0)') 'reconstruction_max_abs_b_k=', &
                            implicit_failure%reconstruction_max_abs_b_k
                        write(output_unit,'(1P,A,E16.8)') 'reconstructed_h_at_failure_k=', &
                            implicit_failure%reconstructed_h_at_failure_k
                        write(output_unit,'(1P,A,E16.8)') 'stored_h_at_failure_k=', &
                            implicit_failure%stored_h_at_failure_k
                        write(output_unit,'(1P,A,E16.8)') 'reconstructed_b_at_failure_k=', &
                            implicit_failure%reconstructed_b_at_failure_k
                        write(output_unit,'(1P,A,E16.8)') 'stored_b_at_failure_k=', &
                            implicit_failure%stored_b_at_failure_k
                    endif
                else
                    write(output_unit,'(A)') 'equilibrium_breakdown=not_applicable'
                endif
                flush(output_unit)
            endif
        else
            write(*,'(A,I0,A,A,A,I0,A,I0,A,I0,A,I0,1P,3(A,E16.8))') &
                'Steady implicit failure rank=',myID,' stage=',trim(stage),' iter=',iter, &
                ' cell_i=',0,' cell_j=',0,' velocity=',0, &
                ' nCFL=',implicit_ncfl_active,' alpha=',implicit_alpha, &
                ' minAlpha=',Implicit_MinAlpha
        endif
        call MPI_Abort(MPI_COMM_WORLD,203,ierr)
    end subroutine AbortSteadyImplicit
end module mod_ImplicitDiagnostics

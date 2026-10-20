module Solver
    use mod_ClosedCouetteMass,only:mass_scale
    use mod_S2MPI1Layout,only:communication_seconds,communication_clock
    use flow_data
    use mod_W3BObserver, only: W3B_Init, W3B_Observe
    use Flow_DIS
    use mpi
    use mpi_wrapper
    use tools
    use mod_Interpolation
    use mod_Evolution
    use mod_Residual
    use mod_TimeIntegrator
    use mod_SteadyImplicit
    use mod_ImplicitDiagnostics, only: WriteAcceptedStateAdmissibilityDiagnostics
    use mod_UnsteadyImplicit, only: ComputeUnsteadyFluxTimeScales, ApplyUnsteadyImplicitStep, &
        WriteUnsteadyImplicitResidual, ProjectBaseResidualMoments, unsteady_macro_residual, &
        unsteady_h_residual, unsteady_b_residual
    use mod_ImplicitState, only: implicit_evolution_seconds,implicit_macro_residual,implicit_h_norm, &
        implicit_b_norm,implicit_moment_error,r5_force_true_updates
    use Boundary
    use InitializeFlow, only: Init_fringe_x
    use mod_R5Observe,only:StartR5Control,WriteR5Checkpoint,R5ShouldStop,FinalR5State, &
        AfterR5Accepted,RecordR5Cost,verified_streak,r5_stop_reason,r6_preflight,EvaluateSparseState
    use mod_NativeWiringChecks,only:DumpPhysicalOperators,RecordWiringState
    use Output
    use mod_S2D1, only: IsS2CN2,IsS2Steady
    use mod_S2NativeAdapter, only: S2NativeApplyUpdate,S2NativeRecordLowOrderStep
    use mod_S2NativeAdapter, only: S2NativeSteadyConverged,S2NativeSteadyMaxNorm,S2NativeSaveCurrent, &
        S2NativeDiagnosticEnabled,S2NativeInvalidateTargets,S2NativeSteadyNeedsEvaluation
    use mod_S2SteadyIteration, only: S2NativeApplySteadyUpdate
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(KREAL), save :: s2_loop_start=0.0_KREAL,s2_loop_total=0.0_KREAL,s2_comm_start=0.0_KREAL,s2_comm_loop=0.0_KREAL
    real(KREAL), save                               :: time_step_last = 0.0_KREAL
    real(KREAL), save                               :: time_boundary_last = 0.0_KREAL
    real(KREAL), save                               :: time_interp_last = 0.0_KREAL
    real(KREAL), save                               :: time_evolution_last = 0.0_KREAL
    real(KREAL), save                               :: time_update_last = 0.0_KREAL
    real(KREAL), save                               :: time_mpi_last = 0.0_KREAL
    real(KREAL), save                               :: time_output_last = 0.0_KREAL
    real(KREAL), save                               :: time_step_acc = 0.0_KREAL
    real(KREAL), save                               :: time_boundary_acc = 0.0_KREAL
    real(KREAL), save                               :: time_interp_acc = 0.0_KREAL
    real(KREAL), save                               :: time_evolution_acc = 0.0_KREAL
    real(KREAL), save                               :: time_update_acc = 0.0_KREAL
    real(KREAL), save                               :: time_mpi_acc = 0.0_KREAL
    real(KREAL), save                               :: time_output_acc = 0.0_KREAL
    logical, save :: stop_requested_by_user=.false.
    logical, save                                   :: residual_output_initialized = .false.
    logical, save                                   :: timing_output_initialized = .false.
    
contains
    subroutine Solve_UGKS
        implicit none
        real(KREAL)                         :: time_now, time_pre_comm, time_end_comm
        real(KREAL)                         :: time_iter_start, time_block_start, time_block_end
        logical                             :: update_accepted,stop_now,converged_now

        s2_loop_start=MPI_Wtime()
        call W3B_Init()
        residual_output_initialized=.false.
        timing_output_initialized=.false.
        !Compute Res0 for the initial steady flow
        call Compute_InitRes0
        
        !call cpu_time(time_pre)
        if(.not.IsS2CN2().or.IsHighOrderImplicit())time_dis = 0.d0
        call Reset_Timing()
        s2_comm_start=communication_seconds

        call RecordWiringState()
        call StartR5Control()
        call WriteR5Checkpoint(.true.)
        if(r6_preflight)then
            call TimeStep()
            call Set_Ghost_PhysBoundary()
            if(.not.IsS2CN2())call Interpolation()
            call Evolution();call ComputeResidual()
            call EvaluateSparseState('restart_readonly_verification',.true.)
            r5_stop_reason='readonly_no_update';call FinalR5State();return
        endif
        !Time iteration
        do while(.true.)
            if(R5ShouldStop())then
                call FinalR5State();exit
            endif
            time_iter_start = MPI_Wtime()

            !Calculate the time step
            time_block_start = MPI_Wtime()
            call TimeStep() 
            if(Iflag_DIS==1 .and. Iflag_Fringe==1) call Init_fringe_x()
            time_block_end = MPI_Wtime()
            time_step_last = time_block_end-time_block_start
            time_step_acc = time_step_acc+time_step_last
            
            !Set ghost value for physical boundary (except for isothermal wall)
            time_block_start = MPI_Wtime()
            call Set_Ghost_PhysBoundary() 
            if(IsImplicitSteady() .or. IsImplicitUnsteady()) call EnforceImplicitBoundaryMomentCompatibility()
            time_block_end = MPI_Wtime()
            time_boundary_last = time_block_end-time_block_start
            time_boundary_acc = time_boundary_acc+time_boundary_last
            
            !Calculate the slope of distribution function
            time_block_start = MPI_Wtime()
            if(.not.IsS2CN2())call Interpolation() 
            time_block_end = MPI_Wtime()
            time_interp_last = time_block_end-time_block_start
            time_interp_acc = time_interp_acc+time_interp_last
            
            !Calculate flux across the interfaces
            time_block_start = MPI_Wtime()
            call Evolution()
            call CheckFiniteFlux('Solve_UGKS loop')
            call ComputeResidual()
            call DumpPhysicalOperators()
            block
                character(len=8)::value
                integer::status
                call get_environment_variable('NATIVE_OPERATOR_ONLY',value,status=status)
                if(status==0.and.trim(value)=='1')exit
            end block
            time_block_end = MPI_Wtime()
            time_evolution_last = time_block_end-time_block_start
            time_evolution_acc = time_evolution_acc+time_evolution_last
            if(IsImplicitSteady()) implicit_evolution_seconds=time_evolution_acc

            ! The current accepted state has just received its unique true-residual
            ! evaluation.  A fifth passing state terminates before any new proposal.
            if(IsS2Steady().and.S2NativeSteadyConverged())then
                call WriteConvergenceStatus(.true.)
                exit
            endif
            
            !Dispatch the accepted-state update while sharing the physical UGKS residual caches.
            time_block_start = MPI_Wtime()
            if(IsHighOrderImplicit())then
                call ApplySteadyImplicitUpdate(update_accepted)
                if(.not.update_accepted)then
                    r5_stop_reason='scientific_rejection'
                    call WriteR5Checkpoint(.true.);call Output_result()
                    call AbortSteadyImplicit('high-order implicit rejected; accepted restart saved')
                endif
            elseif(IsS2CN2()) then
                call S2NativeApplyUpdate(update_accepted)
                if(.not.update_accepted) call AbortSteadyImplicit('S2CN2 physical step was rejected')
            elseif(IsS2Steady()) then
                call S2NativeApplySteadyUpdate(update_accepted)
                if(.not.update_accepted)then
                    call S2NativeSaveCurrent('s2i1_last_accepted.bin')
                    call AbortSteadyImplicit('S2 steady correction rejected')
                endif
            elseif(IsImplicitSteady()) then
                call ApplySteadyImplicitUpdate(update_accepted)
                if(.not.update_accepted) call AbortSteadyImplicit('unaccepted update returned to Solve_UGKS')
            elseif(IsImplicitUnsteady()) then
                call ApplyUnsteadyImplicitStep(update_accepted)
                if(.not.update_accepted) call AbortSteadyImplicit('unsteady physical step failed to converge')
            else
                call ApplyExplicitUpdate()
                if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau')call S2NativeRecordLowOrderStep()
                update_accepted=.true.
            endif
            time_block_end = MPI_Wtime()
            time_update_last = time_block_end-time_block_start
            time_update_acc = time_update_acc+time_update_last
            
            !Update the overlaped region using MPI_Sendrecv
            time_pre_comm = MPI_Wtime( )
            call MPIUpdate_CellConVars()
            call MPIUpdate_CellVDF()
            call MPIUpdate_Conner_conVars()
            call MPIUpdate_Conner_VDF()
            
            !Synchronize different processor 
            communication_clock=MPI_Wtime()
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            
            iter = iter+1
            call RecordWiringState()
            if(((IsS2CN2().and..not.IsHighOrderImplicit()) .or. trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau' .or. &
                IsImplicitUnsteady() .or. Iflag_DIS==1) .and. update_accepted) &
                call AdvancePhysicalTime(dt_global)
            time_end_comm = MPI_Wtime( )
            time_mpi_last = time_end_comm-time_pre_comm
            time_mpi_acc = time_mpi_acc+time_mpi_last
            time_now = MPI_Wtime( )
            !call cpu_time(time_now)
            dcommtime = time_mpi_last
            dsimtime = time_now-time_iter_start
            simtime = simtime+dsimtime
            call W3B_Observe(iter,simtime)
            call AfterR5Accepted()
            call WriteR5Checkpoint(.false.)
    
            !Output the residuals
            if(mod(iter,Nstep_show)==0) then
                if(IsS2Steady())then
                    call WriteAcceptedStateAdmissibilityDiagnostics('S2Steady')
                elseif(IsImplicitSteady()) then
                    call WriteImplicitResidual(simtime,merge(mass_scale-1.0_KREAL,0.0_KREAL,ClosedCouetteMassRescale))
                    call WriteAcceptedStateAdmissibilityDiagnostics('ImplicitSteady')
                elseif(IsImplicitUnsteady()) then
                    call WriteUnsteadyImplicitResidual()
                else
                    call Write_residuals
                    if(.not.IsS2CN2().or.S2Diagnostics)call WriteAcceptedStateAdmissibilityDiagnostics('Explicit')
                endif
            end if
    
            !Determine whether this accepted state is terminal before periodic flow output.
            if(IsS2Steady())then
                converged_now=S2NativeSteadyConverged()
                stop_now=converged_now.or.iter>=Implicit_MaxIter
            elseif(IsImplicitSteady()) then
                converged_now=ImplicitConverged()
                if(r5_force_true_updates)converged_now=verified_streak>=5
                stop_now=converged_now .or. iter>=Implicit_MaxIter
                if(converged_now)r5_stop_reason='five_distinct_accepted_states_passed'
                if(.not.converged_now.and.iter>=Implicit_MaxIter)r5_stop_reason='iteration_limit'
            elseif(IsS2CN2() .or. trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau' .or. &
                IsImplicitUnsteady() .or. Iflag_DIS==1) then
                converged_now=.false.
                stop_now=NonsteadyStopRequested()
            else
                converged_now=all(res<EPS)
                stop_now=converged_now .or. iter>=MAX_ITER
            endif

            !Output non-terminal periodic flow; the main program owns the terminal output.
            if (mod(iter,Nstep_write)==0 .and. .not.stop_now) then
                !call Write_flow
                time_block_start = MPI_Wtime()
                call Output_result(.false.)
                time_block_end = MPI_Wtime()
                time_output_last = time_block_end-time_block_start
                time_output_acc = time_output_acc+time_output_last
            end if
            if(Nstep_restart>0.and..not.stop_now)then
                if(mod(iter,max(1,Nstep_restart))==0)then
                    time_block_start=MPI_Wtime()
                    call Output_restart_only()
                    time_output_acc=time_output_acc+MPI_Wtime()-time_block_start
                endif
            endif
            if(stop_now)call FinalR5State()
            call RecordR5Cost(time_iter_start,time_output_acc)
            !Write terminal residual/status; main_UGKS2D writes the terminal field once.
            if(stop_now) then
                if((IsS2Steady().and.S2NativeSteadyNeedsEvaluation()).or. &
                    (IsS2CN2().and.S2NativeDiagnosticEnabled()))then
                    ! The loop residual preceded the last accepted update.
                    ! Evaluate the actual terminal state once for both methods.
                    if(IsS2Steady())then
                        ! The accepted final state has new ghosts and slopes.
                        ! Rebuild the full spatial entry before its one terminal
                        ! residual evaluation; a stale prior-state face is invalid.
                        call S2NativeInvalidateTargets()
                        call TimeStep()
                        call Set_Ghost_PhysBoundary()
                        if(IsImplicitSteady())call EnforceImplicitBoundaryMomentCompatibility()
                        if(.not.IsS2CN2())call Interpolation()
                    endif
                    call Evolution()
                    call ComputeResidual()
            call DumpPhysicalOperators()
            block
                character(len=8)::value
                integer::status
                call get_environment_variable('NATIVE_OPERATOR_ONLY',value,status=status)
                if(status==0.and.trim(value)=='1')exit
            end block
                    if(IsS2Steady())converged_now=S2NativeSteadyConverged()
                endif
                if(mod(iter,Nstep_show)/=0) then
                    if(IsS2Steady())then
                        continue
                    elseif(IsImplicitSteady()) then
                        call WriteImplicitResidual(simtime,merge(mass_scale-1.0_KREAL,0.0_KREAL,ClosedCouetteMassRescale))
                    elseif(IsImplicitUnsteady()) then
                        call WriteUnsteadyImplicitResidual()
                    else
                        call Write_residuals
                    endif
                endif
                call WriteConvergenceStatus(converged_now)
                exit
            endif
            
        end do
      
        s2_loop_total=MPI_Wtime()-s2_loop_start
        s2_comm_loop=communication_seconds-s2_comm_start
    end subroutine Solve_UGKS

    subroutine AdvancePhysicalTime(accepted_dt)
        real(KREAL), intent(in) :: accepted_dt
        time_dis=time_dis+accepted_dt
    end subroutine AdvancePhysicalTime

    logical function NonsteadyStopRequested()
        integer::local_stop,global_stop,local_request,global_request
        logical::request
        inquire(file='STOP_REQUEST',exist=request)
        local_request=merge(1,0,request)
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(local_request,global_request,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        stop_requested_by_user=global_request/=0
        local_stop=merge(1,0,request.or.iter>=MAX_ITER)
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(local_stop,global_stop,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        NonsteadyStopRequested=global_stop/=0
    end function NonsteadyStopRequested

    subroutine WriteConvergenceStatus(converged)
        logical, intent(in) :: converged
        real(KREAL) :: maximum_residual

        if (myID/=0) return
        if (IsS2Steady()) then
            maximum_residual=S2NativeSteadyMaxNorm()
        elseif (IsImplicitSteady()) then
            maximum_residual=max(maxval(implicit_macro_residual),implicit_h_norm,implicit_b_norm,implicit_moment_error)
        elseif(IsImplicitUnsteady()) then
            maximum_residual=max(maxval(unsteady_macro_residual),unsteady_h_residual,unsteady_b_residual)
        else
            maximum_residual=maxval(res)
        endif
        open(unit=105,file='ConvergenceStatus.txt',status='replace',action='write')
        write(105,'(A)') 'mode converged iteration maxResidual'
        write(105,'(A,1X,L1,1X,I0,1X,1P,E18.9)') trim(adjustl(TimeIntegration_method)), &
            converged,iter,maximum_residual
        close(105)
        if(stop_requested_by_user)then
            write(*,'(A,I0)')'Controlled STOP_REQUEST at accepted step ',iter
            return
        endif
        if (.not.converged) write(*,'(A,A,A,I0,A,1P,E16.8)') &
            'WARNING: ',trim(adjustl(TimeIntegration_method)),' reached iteration limit at ',iter, &
            ' with max residual ',maximum_residual
    end subroutine WriteConvergenceStatus


    !--------------------------------------------------
    !>Reset accumulated profiling timers
    !--------------------------------------------------
    subroutine Reset_Timing()
        time_step_last = 0.0_KREAL
        time_boundary_last = 0.0_KREAL
        time_interp_last = 0.0_KREAL
        time_evolution_last = 0.0_KREAL
        time_update_last = 0.0_KREAL
        time_mpi_last = 0.0_KREAL
        time_output_last = 0.0_KREAL
        time_step_acc = 0.0_KREAL
        time_boundary_acc = 0.0_KREAL
        time_interp_acc = 0.0_KREAL
        time_evolution_acc = 0.0_KREAL
        time_update_acc = 0.0_KREAL
        time_mpi_acc = 0.0_KREAL
        time_output_acc = 0.0_KREAL
    end subroutine Reset_Timing


    !--------------------------------------------------
    !>Calculate time step
    !--------------------------------------------------
    subroutine TimeStep()
        real(KREAL)                                     :: prim(4) !Primary variables
        real(KREAL)                                     :: sos !Speed of sound
        real(KREAL)                                     :: tMin !Min dt allowed
        real(KREAL)                                     :: tbuf_send(1),tbuf_recv(1)
        integer(KINT)                                   :: i,j
        
        if(IsImplicitUnsteady()) then
            call ComputeUnsteadyFluxTimeScales()
            return
        elseif(TimeStep_Method==Time_global_fixed) then
            !dt=dt_global: dt is the default value or reading from config.ini
            ctr(:,:)%dt_local = dt_global
        else
            do j=1,ny
                do i=1,nx
                    !Convert conservative variables to primary variables
                    prim = GetPrimary(ctr(i,j)%conVars)

                    !Get sound speed
                    sos = GetSoundSpeed(prim)
                    !------------- For large Knudsen number ------------
                    prim(2) = max(U_MAX,abs(prim(2))+sos)
                    prim(3) = max(V_MAX,abs(prim(3))+sos)

                    !------------- For large Reynolds number -----------
                    ! prim(2) = abs(prim(2))+sos
                    ! prim(3) = abs(prim(3))+sos

                    !Maximum 1/dt allowed
                    !tMax = max(tMax,(ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3))/ctr(i,j)%area)
                
                    !Calculate local time step
                    ctr(i,j)%dt_local = CFL*ctr(i,j)%area/(ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3)) !Record local time step
                end do
            end do
        endif
        
        !Time step
        if(TimeStep_Method==Time_global_varying) then
            tMin = minval(ctr(1:nx,1:ny)%dt_local)       !minimum of dt(1:nx, 1:ny)
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            tbuf_send(1) = tMin
            call MPI_AllReduce(tbuf_send,tbuf_recv,1,MPI_KREAL, MPI_MIN, MPI_COMM_WORLD,ierr)
            dt_global = tbuf_recv(1)
            ctr(1-Ghost:nx+Ghost,1-Ghost:ny+Ghost)%dt_local = dt_global     !the timestep of the ghost cell and physical cell is dt_global
        endif
  
    end subroutine TimeStep


    !--------------------------------------------------
    !>Write residuals to screen and file
    !--------------------------------------------------
    subroutine Write_residuals
        if(myID==0) then
            !Write screen
            write(*,"(A35,I15,3E15.7)") "iter,simTime,dsimtime,dcommtime:",iter,simTime,dsimTime,dcommtime
            write(*,"(A35,E15.7)") "dtPhysical:",dt_global
            write(*,"(A35,4E15.7)") "res:",res
            call Write_timing_profile
        
            !Write file
            if(.not.residual_output_initialized) then
                open(unit=101,file='Residuals.plt',status="replace",action="write") !Open residual file
                write(101,'(a)') "VARIABLES = iter, resRho, resRhoU, resRhoV, resE, simtime, c_minus_1, mass_rescale_enabled" !Global L2-normalized residuals
                write(101,"(I15,5E15.7,2ES25.16E3)") iter,res,simTime,0.0_KREAL,0.0_KREAL
                close(101)

                open(unit=102,file='Residuals_detail.plt',status="replace",action="write")
                write(102,'(a)') "VARIABLES = iter, resRho, resRhoU, resRhoV, resE, " // &
                    "numRho, numRhoU, numRhoV, numE, denRho, denRhoU, denRhoV, denE, " // &
                    "rmsRho, rmsRhoU, rmsRhoV, rmsE, maxRho, maxRhoU, maxRhoV, maxE, " // &
                    "relativeRho, relativeRhoU, relativeRhoV, relativeE, " // &
                    "absoluteRho, absoluteRhoU, absoluteRhoV, absoluteE, simtime, c_minus_1, mass_rescale_enabled"
                write(102,"(I15,29E15.7,2ES25.16E3)") iter,res,res_num,res_den,res_rms,res_max,res_relative,res_absolute,simTime,0.0_KREAL,0.0_KREAL
                close(102)
                residual_output_initialized=.true.
            else
                open(unit=101,file='Residuals.plt',status="old",action="write",position="append") !Open residual file
                write(101,"(I15,5E15.7,2ES25.16E3)") iter,res,simTime,0.0_KREAL,0.0_KREAL
                close(101)

                open(unit=102,file='Residuals_detail.plt',status="old",action="write",position="append")
                write(102,"(I15,29E15.7,2ES25.16E3)") iter,res,res_num,res_den,res_rms,res_max,res_relative,res_absolute,simTime,0.0_KREAL,0.0_KREAL
                close(102)
            endif   
        endif
    end subroutine Write_residuals


    !--------------------------------------------------
    !>Write profiling timers to screen and file
    !--------------------------------------------------
    subroutine Write_timing_profile
        real(KREAL)                                     :: total_acc,total_last
        real(KREAL)                                     :: pct_step,pct_boundary,pct_interp,pct_evolution
        real(KREAL)                                     :: pct_update,pct_mpi,pct_output
        total_acc = time_step_acc+time_boundary_acc+time_interp_acc+time_evolution_acc+&
                    time_update_acc+time_mpi_acc+time_output_acc
        total_last = time_step_last+time_boundary_last+time_interp_last+time_evolution_last+&
                     time_update_last+time_mpi_last+time_output_last

        if(total_acc>0.0_KREAL) then
            pct_step = 100.0_KREAL*time_step_acc/total_acc
            pct_boundary = 100.0_KREAL*time_boundary_acc/total_acc
            pct_interp = 100.0_KREAL*time_interp_acc/total_acc
            pct_evolution = 100.0_KREAL*time_evolution_acc/total_acc
            pct_update = 100.0_KREAL*time_update_acc/total_acc
            pct_mpi = 100.0_KREAL*time_mpi_acc/total_acc
            pct_output = 100.0_KREAL*time_output_acc/total_acc
        else
            pct_step = 0.0_KREAL
            pct_boundary = 0.0_KREAL
            pct_interp = 0.0_KREAL
            pct_evolution = 0.0_KREAL
            pct_update = 0.0_KREAL
            pct_mpi = 0.0_KREAL
            pct_output = 0.0_KREAL
        end if

        write(*,"(A35,8E12.4)") "timing last total/parts:",total_last,time_step_last,time_boundary_last,&
             time_interp_last,time_evolution_last,time_update_last,time_mpi_last,time_output_last
        write(*,"(A35,7F9.2)") "timing cumulative pct:",pct_step,pct_boundary,pct_interp,&
             pct_evolution,pct_update,pct_mpi,pct_output

        if(.not.timing_output_initialized) then
            open(unit=103,file='Timing_detail.plt',status="replace",action="write")
            write(103,'(a)') "VARIABLES = iter, step_last, boundary_last, interp_last, evolution_last, update_last, mpi_last, output_last, total_last, step_pct, boundary_pct, interp_pct, evolution_pct, update_pct, mpi_pct, output_pct"
            write(103,"(I15,15E15.7)") iter,time_step_last,time_boundary_last,time_interp_last,&
                time_evolution_last,time_update_last,time_mpi_last,time_output_last,total_last,&
                pct_step,pct_boundary,pct_interp,pct_evolution,pct_update,pct_mpi,pct_output
            close(103)
            timing_output_initialized=.true.
        else
            open(unit=103,file='Timing_detail.plt',status="old",action="write",position="append")
            write(103,"(I15,15E15.7)") iter,time_step_last,time_boundary_last,time_interp_last,&
                time_evolution_last,time_update_last,time_mpi_last,time_output_last,total_last,&
                pct_step,pct_boundary,pct_interp,pct_evolution,pct_update,pct_mpi,pct_output
            close(103)
        endif
    end subroutine Write_timing_profile
  
    
    !--------------------------------------------------
    !>Compute residuals of the initial flow
    !--------------------------------------------------
    subroutine Compute_InitRes0
        real(KREAL), allocatable, dimension(:)          :: H_old_vec,B_old_vec                              !Pure 1D old-time equilibrium arrays
        real(KREAL), allocatable, dimension(:)          :: H_plus_vec,B_plus_vec                            !Pure 1D Shakhov correction arrays
        real(KREAL)                                     :: conVars_old(4),prim_old(4)                       !Conversative and primary variables at t^n
        real(KREAL)                                     :: tau_old                                          !Collision time at t^n
        real(KREAL)                                     :: qf(2)
        real(KREAL)                                     :: dt
        real(KREAL)                                     :: implicit_inv                                     !Precomputed collision scalar for Nvel residual loop
        integer(KINT)                                   :: i,j,k                                            !k is flattened velocity-space index
        logical                                         :: base_ok

        if(Iflag_DIS==0) then               !Steady base flow simulation: Res0 is set 0
            do j=1,ny
                do i=1,nx
                    ctr(i,j)%res_conVars0 = 0.d0
                    !Clear initial VDF residuals through pure 1D velocity-space views.
                    ctr(i,j)%res_h0_vec = 0.d0
                    ctr(i,j)%res_b0_vec = 0.d0
                enddo
            enddo
            
        elseif(Iflag_DIS==1) then           !Unsteady DIS simulation: Res0 is needed
            Iflag_DIS = 0   !!!!Compute steady res.
            
            !Allocate pure 1D Maxwellian/Shakhov arrays for initial residual storage.
            allocate(H_old_vec(Nvel))
            allocate(B_old_vec(Nvel))
            allocate(H_plus_vec(Nvel))
            allocate(B_plus_vec(Nvel))
        
            !Calculate the time step
            call TimeStep() 
            
            !Set ghost value for physical boundary (except for isothermal wall)
            call Set_Ghost_PhysBoundary() 
            
            !Calculate the slope of distribution function
            if(.not.IsS2CN2())call Interpolation() 
            
            !Calculate flux across the interfaces
            call Evolution()
            call CheckFiniteFlux('Compute_InitRes0')
            call ComputeResidual()


            !Compute initial res. (Update)
            do j=1,ny
                do i=1,nx
                    !Get local time step
                    dt = ctr(i,j)%dt_local
                
                    !--------------------------------------------------
                    !##Residuals of conVars at initial state (notice: res. contain dt)
                    !--------------------------------------------------
                    !Store conVars^n and calculate H^n,B^n,\tau^n
                    conVars_old = ctr(i,j)%conVars !Store conVars^n
                    prim_old = GetPrimary(conVars_old) !Convert to primary variables
                    call DiscreteMaxwell_Vec(H_old_vec,B_old_vec,prim_old) !Calculate Maxwellian
                    tau_old = GetTau(prim_old) !Calculate collision time \tau^n

                    !Unsteady implicit equations consume an operator rate; the legacy explicit path
                    !keeps its step-integrated base-flow increment.
                    if(IsImplicitUnsteady()) then
                        ctr(i,j)%res_conVars0 = res_macro_flux_store(1:4,i,j)/dt
                    else
                        ctr(i,j)%res_conVars0 = res_macro_flux_store(1:4,i,j)
                    endif

                
                    !--------------------------------------------------
                    !##Residuals of VDF at initial state (notice res. contain dt)
                    !--------------------------------------------------
                    !Shakhov part
                    !Calculate heat flux at t=t^n directly from flattened cell-center VDF views.
                    qf = GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim_old)

                    !h^+ = H+H^+ at t=t^n
                    call ShakhovPart_Vec(H_old_vec,B_old_vec,qf,prim_old,H_plus_vec,B_plus_vec) !H^+ and B^+
                    !Add Shakhov correction through 1D views for initial residual storage.
                    do k=1,Nvel
                        H_old_vec(k) = H_old_vec(k)+H_plus_vec(k)
                        B_old_vec(k) = B_old_vec(k)+B_plus_vec(k)
                    enddo

                    !Res. of VDF
                    !Store initial residual VDFs from the shared micro face-flux operator cache.
                    if(IsImplicitUnsteady()) then
                        do k=1,Nvel
                            ctr(i,j)%res_h0_vec(k) = res_h_flux_store(k,i,j)/dt+ &
                                (H_old_vec(k)-ctr(i,j)%h_vec(k))/tau_old
                            ctr(i,j)%res_b0_vec(k) = res_b_flux_store(k,i,j)/dt+ &
                                (B_old_vec(k)-ctr(i,j)%b_vec(k))/tau_old
                        enddo
                    else
                        implicit_inv = 1.0_KREAL/(1.0_KREAL+0.5_KREAL*dt/tau_old)
                        do k=1,Nvel
                            ctr(i,j)%res_h0_vec(k) = (res_h_flux_store(k,i,j)+&
                                dt*(H_old_vec(k)-ctr(i,j)%h_vec(k))/tau_old)*implicit_inv
                            ctr(i,j)%res_b0_vec(k) = (res_b_flux_store(k,i,j)+&
                                dt*(B_old_vec(k)-ctr(i,j)%b_vec(k))/tau_old)*implicit_inv
                        enddo
                    endif
                end do
            end do
            
            if(bc_type(3)==BCconst_wall_isothermal) then    !bottom boundary
                do i=1,nx
                    hface(i,1)%FlowVars0(1:5)       = hface(i,1)%FlowVars(1:5)
                    hface(i,1)%Aerodynamics0(1:3)   = hface(i,1)%Aerodynamics(1:3)
                end do
            endif
        
            !Deallocate arrays
            deallocate(H_old_vec)
            deallocate(B_old_vec)
            deallocate(H_plus_vec)
            deallocate(B_plus_vec)
            
            Iflag_DIS = 1   !!!!when steady res. is computed, set Iflag_DIS=1 to continue unsteady simulation
            if(IsImplicitUnsteady()) then
                call ProjectBaseResidualMoments(base_ok)
                if(.not.base_ok) call AbortSteadyImplicit('base residual moment projection failed')
            endif
        endif

        do j=1,ny
            do i=1,nx
                if(.not.ieee_is_finite(ctr(i,j)%area) .or. ctr(i,j)%area<=0._KREAL) then
                    print '(A,3I6,1P,E15.7)', 'Invalid cell area detected (rank,i,j,area):', myID,i,j,ctr(i,j)%area
                    call MPI_Abort(MPI_COMM_WORLD, 9, ierr)
                end if
                if(.not.all(ieee_is_finite(ctr(i,j)%res_conVars0))) then
                    print '(A,3I6,1P,4E15.7)', 'Non-finite res_conVars0 (rank,i,j,res):', myID,i,j,ctr(i,j)%res_conVars0
                    call MPI_Abort(MPI_COMM_WORLD, 10, ierr)
                end if
                if(.not.all(ieee_is_finite(ctr(i,j)%res_h0_vec))) then
                    print '(A,3I6)', 'Non-finite res_h0 (rank,i,j):', myID,i,j
                    call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
                end if
                if(.not.all(ieee_is_finite(ctr(i,j)%res_b0_vec))) then
                    print '(A,3I6)', 'Non-finite res_b0 (rank,i,j):', myID,i,j
                    call MPI_Abort(MPI_COMM_WORLD, 12, ierr)
                end if
            end do
        end do
        
    end subroutine Compute_InitRes0

    subroutine CheckFiniteFlux(stage)
        character(len=*),intent(in) :: stage
        integer(KINT) :: i,j

        do j=1,ny
            do i=1,nx+1
                if(.not.all(ieee_is_finite(vface(i,j)%L_flux))) then
                    print *, 'Non-finite vface L_flux, stage=', trim(stage), ' rank=', myID, ' i=', i, ' j=', j
                    print '(1P,4E15.7)', vface(i,j)%L_flux
                    call MPI_Abort(MPI_COMM_WORLD, 13, ierr)
                end if
                if(.not.all(ieee_is_finite(vface(i,j)%R_flux))) then
                    print *, 'Non-finite vface R_flux, stage=', trim(stage), ' rank=', myID, ' i=', i, ' j=', j
                    print '(1P,4E15.7)', vface(i,j)%R_flux
                    call MPI_Abort(MPI_COMM_WORLD, 14, ierr)
                end if
            end do
        end do

        do j=1,ny+1
            do i=1,nx
                if(.not.all(ieee_is_finite(hface(i,j)%L_flux))) then
                    print *, 'Non-finite hface L_flux, stage=', trim(stage), ' rank=', myID, ' i=', i, ' j=', j
                    print '(1P,4E15.7)', hface(i,j)%L_flux
                    call MPI_Abort(MPI_COMM_WORLD, 15, ierr)
                end if
                if(.not.all(ieee_is_finite(hface(i,j)%R_flux))) then
                    print *, 'Non-finite hface R_flux, stage=', trim(stage), ' rank=', myID, ' i=', i, ' j=', j
                    print '(1P,4E15.7)', hface(i,j)%R_flux
                    call MPI_Abort(MPI_COMM_WORLD, 16, ierr)
                end if
            end do
        end do

    end subroutine CheckFiniteFlux

end module Solver

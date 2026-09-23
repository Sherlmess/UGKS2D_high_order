module mod_TimeIntegrator
    use flow_data
    use Flow_DIS
    use mpi
    use tools
    use, intrinsic :: ieee_arithmetic
    implicit none

contains

    !--------------------------------------------------
    !>Apply the current explicit time update to cell averaged values and VDFs.
    !--------------------------------------------------
    subroutine ApplyExplicitUpdate()
        real(KREAL), allocatable, dimension(:)          :: H_old_vec,B_old_vec,H_vec,B_vec    !Pure 1D Maxwellian arrays
        real(KREAL), allocatable, dimension(:)          :: H_plus_vec,B_plus_vec              !Pure 1D Shakhov correction arrays
        real(KREAL)                                     :: conVars_old(4),prim_old(4),prim(4) !Conservative and primary variables at t^n and t^{n+1}
        real(KREAL)                                     :: tau_old,tau                       !Collision time at t^n and t^{n+1}
        real(KREAL)                                     :: sumRes(4),sumAvg(4),sumRes_all(4),sumAvg_all(4)
        real(KREAL)                                     :: areaSquared,areaSquared_all
        real(KREAL)                                     :: maxCell_local(4),maxCell_all(4)
        real(KREAL)                                     :: dt,ordinary_residual_scale
        real(KREAL)                                     :: dbg_area
        integer(KINT)                                   :: i,j

        !Allocate pure 1D Maxwellian/Shakhov arrays for the update loop.
        allocate(H_old_vec(Nvel))
        allocate(B_old_vec(Nvel))
        allocate(H_vec(Nvel))
        allocate(B_vec(Nvel))
        allocate(H_plus_vec(Nvel))
        allocate(B_plus_vec(Nvel))

        !Reset residual accumulators for this explicit update.
        res = 0.0
        sumRes = 0.0
        sumAvg = 0.0
        areaSquared = 0.0_KREAL
        maxCell_local = 0.0_KREAL

        do j=1,ny
            do i=1,nx
                !Validate geometry and local physical timestep before using them in update denominators.
                dt = ctr(i,j)%dt_local
                dbg_area = ctr(i,j)%area

                if (.not.ieee_is_finite(dt) .or. dt<=0._KREAL) then
                    print '(A,3I6,1P,2E15.7)', 'Invalid dt_local detected (rank,i,j,dt,area):', myID,i,j,dt,dbg_area
                    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
                end if

                if (dbg_area<=0._KREAL .or. .not.ieee_is_finite(dbg_area)) then
                    print '(A,3I6,1P,2E15.7)', 'Invalid cell area detected (rank,i,j,area,dt):', myID,i,j,dbg_area,dt
                    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
                end if

                !Store W^n and build the old-time equilibrium needed by the trapezoidal collision update.
                conVars_old = ctr(i,j)%conVars
                prim_old = GetPrimary(conVars_old)
                call DiscreteMaxwell_Vec(H_old_vec,B_old_vec,prim_old)
                tau_old = GetTau(prim_old)
                if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau')then
                    ComparisonCellTauMin=min(ComparisonCellTauMin,tau_old)
                    ComparisonCellTauMax=max(ComparisonCellTauMax,tau_old)
                    ComparisonCellTauCalls=ComparisonCellTauCalls+1
                endif
                if (.not.ieee_is_finite(tau_old) .or. tau_old<=0._KREAL) then
                    print '(A,3I6,1P,4E15.7)', 'Invalid tau_old (rank,i,j,tau_old,prim_old):', myID,i,j,tau_old,prim_old
                    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
                end if

                !Apply the explicit macro RHS before constructing the new-time equilibrium.
                call ApplyMacroExplicitUpdate(i,j)

                !Build the new-time equilibrium from W^{n+1}; the micro update uses it immediately below.
                prim = GetPrimary(ctr(i,j)%conVars)
                call DiscreteMaxwell_Vec(H_vec,B_vec,prim)
                tau = GetTau(prim)
                if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau')then
                    ComparisonCellTauMin=min(ComparisonCellTauMin,tau)
                    ComparisonCellTauMax=max(ComparisonCellTauMax,tau)
                    ComparisonCellTauCalls=ComparisonCellTauCalls+1
                endif
                if (.not.ieee_is_finite(tau) .or. tau<=0._KREAL) then
                    print '(A,3I6,1P,5E15.7)', 'Invalid tau (rank,i,j,tau,prim):', myID,i,j,tau,prim
                    call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
                end if

                !Accumulate convergence diagnostics from the macro update only.
                call AccumulateMacroResidualNorm(i,j,dt,conVars_old,sumRes,sumAvg,maxCell_local)
                areaSquared=areaSquared+ctr(i,j)%area**2

                !Apply the explicit micro update after macro variables define the new-time equilibrium.
                call ApplyMicroExplicitUpdate(i,j,dt,prim_old,prim,tau_old,tau,&
                    H_old_vec,B_old_vec,H_vec,B_vec,H_plus_vec,B_plus_vec)
            end do
        end do

        !Calculate global residual diagnostics for output and convergence checks.
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_AllReduce(sumRes(1:4),sumRes_all(1:4),4,MPI_KREAL, MPI_SUM, MPI_COMM_WORLD,ierr)
        call MPI_AllReduce(sumAvg(1:4),sumAvg_all(1:4),4,MPI_KREAL, MPI_SUM, MPI_COMM_WORLD,ierr)
        call MPI_AllReduce(areaSquared,areaSquared_all,1,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_AllReduce(maxCell_local(1:4),maxCell_all(1:4),4,MPI_KREAL, MPI_MAX, MPI_COMM_WORLD,ierr)

        res_num = sqrt(Nx_all*Ny_all*sumRes_all)
        res_den = sumAvg_all+SMV
        res_rms = sqrt(sumRes_all/(Nx_all*Ny_all+SMV))
        res_max = maxCell_all
        res_relative=res_num/res_den
        res_absolute=sqrt(sumRes_all/max(areaSquared_all,SMV))
        ordinary_residual_scale=max(sumAvg_all(1),sumAvg_all(4),SMV)
        do i=1,4
            if (IsNearZeroResidualScale(sumAvg_all(i),ordinary_residual_scale)) then
                res(i)=res_absolute(i)
            else
                res(i)=res_relative(i)
            endif
        enddo

        if(.not.all(ieee_is_finite(res))) then
            print '(A,I6,1P,4E15.7)', 'Non-finite residual detected (rank,res):', myID,res
            call MPI_Abort(MPI_COMM_WORLD, 8, ierr)
        endif

        deallocate(H_old_vec)
        deallocate(B_old_vec)
        deallocate(H_vec)
        deallocate(B_vec)
        deallocate(H_plus_vec)
        deallocate(B_plus_vec)
    end subroutine ApplyExplicitUpdate


    !--------------------------------------------------
    !>Apply the already assembled macro explicit RHS to one cell.
    !--------------------------------------------------
    subroutine ApplyMacroExplicitUpdate(i,j)
        integer(KINT), intent(in)                       :: i,j
        real(KREAL)                                     :: macro_rhs(4)

        !Update conVars through the assembled macro residual/update RHS.
        macro_rhs = res_macro_store(1:4,i,j)
        if(.not.all(ieee_is_finite(macro_rhs))) then
            print *, 'Non-finite macro residual RHS detected at (i,j)=', i,j,' rank=',myID
            print '(1P,4E15.7)', macro_rhs
            call MPI_Abort(MPI_COMM_WORLD, 17, ierr)
        end if
        ctr(i,j)%conVars = ctr(i,j)%conVars+macro_rhs

        if (.not.all(ieee_is_finite(ctr(i,j)%conVars))) then
            print '(A,3I6,1P,4E15.7)', 'Non-finite conservative vars (rank,i,j,conVars):', myID,i,j,ctr(i,j)%conVars
            call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
        end if
    end subroutine ApplyMacroExplicitUpdate


    !--------------------------------------------------
    !>Accumulate macro residual norms after the macro explicit update.
    !--------------------------------------------------
    subroutine AccumulateMacroResidualNorm(i,j,dt,conVars_old,sumRes,sumAvg,maxCell_local)
        integer(KINT), intent(in)                       :: i,j
        real(KREAL), intent(in)                         :: dt
        real(KREAL), intent(in)                         :: conVars_old(4)
        real(KREAL), intent(inout)                      :: sumRes(4),sumAvg(4),maxCell_local(4)

        !Keep the legacy residual definition: macro conservative change divided by local dt.
        ctr(i,j)%res_cell = abs((conVars_old-ctr(i,j)%conVars)*ctr(i,j)%area/dt)
        maxCell_local = max(maxCell_local, ctr(i,j)%res_cell)
        sumRes = sumRes+((conVars_old-ctr(i,j)%conVars)*ctr(i,j)%area/dt)**2
        if (.not.all(ieee_is_finite(sumRes))) then
            print '(A,3I6,1P,12E15.7)', 'Non-finite sumRes (rank,i,j,conVars_old,conVars,area,dt,sumRes):', myID,i,j,&
                  conVars_old,ctr(i,j)%conVars,ctr(i,j)%area,dt,sumRes
            call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
        end if
        if (.not.all(ieee_is_finite(ctr(i,j)%res_cell))) then
            print '(A,3I6,1P,8E15.7)', 'Non-finite cell residual (rank,i,j,res_cell,area,dt):', myID,i,j,&
                  ctr(i,j)%res_cell,ctr(i,j)%area,dt
            call MPI_Abort(MPI_COMM_WORLD, 20, ierr)
        end if

        !sumAvg is the normalization denominator used by Write_residuals.
        sumAvg = sumAvg+abs(ctr(i,j)%conVars)*ctr(i,j)%area
        if (.not.all(ieee_is_finite(sumAvg))) then
            print '(A,3I6,1P,9E15.7)', 'Non-finite sumAvg (rank,i,j,conVars,area,sumAvg):', myID,i,j,&
                  ctr(i,j)%conVars,ctr(i,j)%area,sumAvg
            call MPI_Abort(MPI_COMM_WORLD, 7, ierr)
        end if
    end subroutine AccumulateMacroResidualNorm


    !--------------------------------------------------
    !>Apply the explicit VDF update after macro variables define the new equilibrium.
    !--------------------------------------------------
    subroutine ApplyMicroExplicitUpdate(i,j,dt,prim_old,prim,tau_old,tau,&
        H_old_vec,B_old_vec,H_vec,B_vec,H_plus_vec,B_plus_vec)
        integer(KINT), intent(in)                       :: i,j
        real(KREAL), intent(in)                         :: dt,tau_old,tau
        real(KREAL), intent(in)                         :: prim_old(4),prim(4)
        real(KREAL), intent(inout), dimension(:)        :: H_old_vec,B_old_vec,H_vec,B_vec
        real(KREAL), intent(inout), dimension(:)        :: H_plus_vec,B_plus_vec
        real(KREAL)                                     :: qf(2)
        real(KREAL)                                     :: h_old_val,b_old_val
        real(KREAL)                                     :: h_flux_rhs,b_flux_rhs
        real(KREAL)                                     :: implicit_inv,relax_rate
        integer(KINT)                                   :: k

        !Calculate heat flux at t=t^n directly from flattened cell-center VDF views.
        qf = GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim_old)

        !Add the Shakhov correction to the old-time equilibrium arrays.
        call ShakhovPart_Vec(H_old_vec,B_old_vec,qf,prim_old,H_plus_vec,B_plus_vec)
        do k=1,Nvel
            H_old_vec(k) = H_old_vec(k)+H_plus_vec(k)
            B_old_vec(k) = B_old_vec(k)+B_plus_vec(k)
        enddo

        !Add the Shakhov correction to the new-time equilibrium arrays.
        call ShakhovPart_Vec(H_vec,B_vec,qf,prim,H_plus_vec,B_plus_vec)
        do k=1,Nvel
            H_vec(k) = H_vec(k)+H_plus_vec(k)
            B_vec(k) = B_vec(k)+B_plus_vec(k)
        enddo

        !Update micro distribution values through 1D views.
        implicit_inv = 1.0_KREAL/(1.0_KREAL+0.5_KREAL*dt/tau)
        relax_rate = dt*Xlamd(i,j)
        do k=1,Nvel
            h_old_val = ctr(i,j)%h_vec(k)
            b_old_val = ctr(i,j)%b_vec(k)
            !Read the cached streaming operator so this updater only handles time/collision/source algebra.
            h_flux_rhs = res_h_flux_store(k,i,j)
            b_flux_rhs = res_b_flux_store(k,i,j)
            ctr(i,j)%h_vec(k) = (h_old_val+&
                h_flux_rhs+&
                0.5_KREAL*dt*(H_vec(k)/tau+(H_old_vec(k)-h_old_val)/tau_old))*implicit_inv-&
                ctr(i,j)%res_h0_vec(k)+relax_rate*(ctr(i,j)%h0_vec(k)-h_old_val)
            ctr(i,j)%b_vec(k) = (b_old_val+&
                b_flux_rhs+&
                0.5_KREAL*dt*(B_vec(k)/tau+(B_old_vec(k)-b_old_val)/tau_old))*implicit_inv-&
                ctr(i,j)%res_b0_vec(k)+relax_rate*(ctr(i,j)%b0_vec(k)-b_old_val)
        end do
    end subroutine ApplyMicroExplicitUpdate

end module mod_TimeIntegrator

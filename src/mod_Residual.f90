module mod_Residual
    use flow_data
    use Flow_DIS
    use mpi
    use mod_S2D1, only: IsS2CN2,IsS2Steady
    use mod_S2NativeAdapter, only: S2NativeMarkResidual,S2NativeAssembleSteadyResidual,S2NativeDiagnosticEnabled
    implicit none

contains

    !--------------------------------------------------
    !>Compute the current residual/operator caches from the face fluxes built by Evolution().
    !--------------------------------------------------
    subroutine ComputeResidual()
        !This is the single R(U) entry point for explicit, high-order, and future implicit paths.
        call AssembleMacroResidual()
        call AssembleMicroFluxResidual()
        if(IsS2CN2().or.IsS2Steady()) call S2NativeMarkResidual()
        if(IsS2Steady().or.S2NativeDiagnosticEnabled()) call S2NativeAssembleSteadyResidual()
    end subroutine ComputeResidual


    !--------------------------------------------------
    !>Assemble macro residual/update RHS from face fluxes and current source-like terms.
    !--------------------------------------------------
    subroutine AssembleMacroResidual()
        real(KREAL)                                     :: flux_contrib(4),res_contrib(4),relax_contrib(4)
        integer(KINT)                                   :: i,j

        do j=1,ny
            do i=1,nx
                !Face flux contribution is the physical UGKS operator evaluated from Evolution().
                flux_contrib = (vface(i,j)%L_flux-vface(i+1,j)%R_flux+&
                                hface(i,j)%L_flux-hface(i,j+1)%R_flux)/ctr(i,j)%area
                res_macro_flux_store(1:4,i,j) = flux_contrib

                !Disturbance simulations subtract the stored base-flow residual.
                res_contrib = ctr(i,j)%res_conVars0

                !Fringe forcing is kept in the assembled RHS so the explicit update remains unchanged.
                if (abs(Xlamd(i,j))<=epsilon(1.0_KREAL)) then
                    relax_contrib = 0.0_KREAL
                else
                    relax_contrib = ctr(i,j)%dt_local*Xlamd(i,j)*(ctr(i,j)%conVars0-ctr(i,j)%conVars)
                endif

                res_macro_store(1:4,i,j) = flux_contrib-res_contrib+relax_contrib
            enddo
        enddo
    end subroutine AssembleMacroResidual


    !--------------------------------------------------
    !>Assemble micro face-flux contributions from h/b interface fluxes.
    !--------------------------------------------------
    subroutine AssembleMicroFluxResidual()
        real(KREAL)                                     :: area_inv
        integer(KINT)                                   :: i,j,k

        do j=1,ny
            do i=1,nx
                !Only the streaming/face-flux part is cached here; collision and sources stay in ApplyExplicitUpdate().
                area_inv = 1.0_KREAL/ctr(i,j)%area
                do k=1,Nvel
                    res_h_flux_store(k,i,j) = (vface(i,j)%L_flux_h_vec(k)-vface(i+1,j)%R_flux_h_vec(k)+&
                        hface(i,j)%L_flux_h_vec(k)-hface(i,j+1)%R_flux_h_vec(k))*area_inv
                    res_b_flux_store(k,i,j) = (vface(i,j)%L_flux_b_vec(k)-vface(i+1,j)%R_flux_b_vec(k)+&
                        hface(i,j)%L_flux_b_vec(k)-hface(i,j+1)%R_flux_b_vec(k))*area_inv
                enddo
            enddo
        enddo
    end subroutine AssembleMicroFluxResidual

    !--------------------------------------------------
    !>Convert physical UGKS step increments into steady residual rates.
    !--------------------------------------------------
    subroutine AssembleSteadyRateResiduals()
        use mod_ImplicitState, only: rate_w, stream_h, stream_b
        use, intrinsic :: ieee_arithmetic
        integer(KINT) :: i, j
        real(KREAL) :: dt_inv

        do j=1,ny
            do i=1,nx
                if (.not.ieee_is_finite(ctr(i,j)%dt_local) .or. ctr(i,j)%dt_local<=0.0_KREAL) then
                    write(*,'(A,3I8,1P,E16.8)') 'Invalid physical dt in residual-rate conversion:', &
                                                myID, i, j, ctr(i,j)%dt_local
                    call MPI_Abort(MPI_COMM_WORLD, 201, ierr)
                endif
                dt_inv=1.0_KREAL/ctr(i,j)%dt_local
                rate_w(:,i,j)=res_macro_flux_store(:,i,j)*dt_inv
                stream_h(:,i,j)=res_h_flux_store(:,i,j)*dt_inv
                stream_b(:,i,j)=res_b_flux_store(:,i,j)*dt_inv
            enddo
        enddo
    end subroutine AssembleSteadyRateResiduals

end module mod_Residual

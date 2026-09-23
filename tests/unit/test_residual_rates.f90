program test_residual_rates
    use mpi
    use flow_data
    use mod_ImplicitState
    use mod_Residual
    use mod_TimeIntegrator, only: ApplyMacroExplicitUpdate
    use TestSupport
    implicit none
    integer :: ierr_local,i,j
    real(KREAL) :: w_old(4)

    call MPI_Init(ierr_local)
    nx=1; ny=1; Nvel=2
    allocate(ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    allocate(res_macro_flux_store(4,nx,ny), res_h_flux_store(Nvel,nx,ny), &
             res_b_flux_store(Nvel,nx,ny),res_macro_store(4,nx,ny),Xlamd(nx,ny))
    do j=1,ny
        do i=1,nx+1
            allocate(vface(i,j)%L_flux_h_vec(Nvel),vface(i,j)%L_flux_b_vec(Nvel))
            allocate(vface(i,j)%R_flux_h_vec(Nvel),vface(i,j)%R_flux_b_vec(Nvel))
            vface(i,j)%L_flux_h_vec=0.0_KREAL; vface(i,j)%L_flux_b_vec=0.0_KREAL
            vface(i,j)%R_flux_h_vec=0.0_KREAL; vface(i,j)%R_flux_b_vec=0.0_KREAL
        enddo
    enddo
    do j=1,ny+1
        do i=1,nx
            allocate(hface(i,j)%L_flux_h_vec(Nvel),hface(i,j)%L_flux_b_vec(Nvel))
            allocate(hface(i,j)%R_flux_h_vec(Nvel),hface(i,j)%R_flux_b_vec(Nvel))
            hface(i,j)%L_flux_h_vec=0.0_KREAL; hface(i,j)%L_flux_b_vec=0.0_KREAL
            hface(i,j)%R_flux_h_vec=0.0_KREAL; hface(i,j)%R_flux_b_vec=0.0_KREAL
        enddo
    enddo
    call AllocateImplicitState()
    ctr(1,1)%dt_local=0.25_KREAL; ctr(1,1)%area=2.0_KREAL
    ctr(1,1)%conVars=[1.0_KREAL,0.0_KREAL,0.0_KREAL,2.0_KREAL]
    ctr(1,1)%res_conVars0=0.0_KREAL; Xlamd=0.0_KREAL
    vface(1,1)%L_flux=[1.0_KREAL,2.0_KREAL,3.0_KREAL,4.0_KREAL]
    vface(2,1)%R_flux=[0.2_KREAL,0.4_KREAL,0.6_KREAL,0.8_KREAL]
    hface(1,1)%L_flux=[0.3_KREAL,0.1_KREAL,-0.1_KREAL,0.5_KREAL]
    hface(1,2)%R_flux=[0.1_KREAL,0.2_KREAL,0.3_KREAL,0.4_KREAL]
    vface(1,1)%L_flux_h_vec=[0.4_KREAL,-0.2_KREAL]
    vface(2,1)%R_flux_h_vec=[0.1_KREAL,0.05_KREAL]
    hface(1,1)%L_flux_h_vec=[0.2_KREAL,0.3_KREAL]
    hface(1,2)%R_flux_h_vec=[0.05_KREAL,0.1_KREAL]
    vface(1,1)%L_flux_b_vec=[0.6_KREAL,0.2_KREAL]
    vface(2,1)%R_flux_b_vec=[0.2_KREAL,0.1_KREAL]
    hface(1,1)%L_flux_b_vec=[0.1_KREAL,0.4_KREAL]
    hface(1,2)%R_flux_b_vec=[0.3_KREAL,0.2_KREAL]
    call ComputeResidual()
    call AssembleSteadyRateResiduals()
    call AssertTrue(maxval(abs(rate_w(:,1,1)-[2.0_KREAL,3.0_KREAL,4.0_KREAL,6.6_KREAL])) &
                    <1.0e-14_KREAL, 'macro rate')
    call AssertTrue(maxval(abs(stream_h(:,1,1)-[0.9_KREAL,-0.1_KREAL]))<1.0e-14_KREAL, 'h rate')
    call AssertTrue(maxval(abs(stream_b(:,1,1)-[0.4_KREAL,0.6_KREAL]))<1.0e-14_KREAL, 'b rate')
    call AssertTrue(maxval(abs(res_h_flux_store(:,1,1)-ctr(1,1)%dt_local*stream_h(:,1,1))) &
                    <1.0e-14_KREAL,'h physical increment identity')
    call AssertTrue(maxval(abs(res_b_flux_store(:,1,1)-ctr(1,1)%dt_local*stream_b(:,1,1))) &
                    <1.0e-14_KREAL,'b physical increment identity')
    w_old=ctr(1,1)%conVars
    call ApplyMacroExplicitUpdate(1,1)
    call AssertTrue(maxval(abs((ctr(1,1)%conVars-w_old)-ctr(1,1)%dt_local*rate_w(:,1,1))) &
                    <1.0e-14_KREAL,'explicit macro increment identity')
    call FinishTests('residual_rates')
    call MPI_Finalize(ierr_local)
end program test_residual_rates

program test_implicit_macro
    use mpi
    use flow_data
    use Tools, only: GetConserved
    use mod_ImplicitState
    use mod_ImplicitMacro
    use TestSupport
    implicit none
    integer :: ierr_local, i, j
    logical :: ok
    real(KREAL) :: gamma_face, defect_before, defect_after, alpha, w_uniform(4)

    call MPI_Init(ierr_local)
    call MPI_Comm_rank(MPI_COMM_WORLD,myID,ierr_local)
    nx=2; ny=1; Nvel=1
    gamma=1.4_KREAL; CFL=0.5_KREAL; Implicit_nCFL=10.0_KREAL
    Implicit_MacroSweeps=2; Implicit_MinAlpha=1.0e-6_KREAL
    flag_vislaw=PowerLaw; mu_ref=0.01_KREAL; omega=0.5_KREAL
    BC_Type=BCconst_freestream
    ID_X_minus=MPI_PROC_NULL; ID_X_plus=MPI_PROC_NULL
    ID_Y_minus=MPI_PROC_NULL; ID_Y_plus=MPI_PROC_NULL
    allocate(ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    w_uniform=GetConserved([1.0_KREAL,0.2_KREAL,0.0_KREAL,1.0_KREAL])
    do j=1-LAP,ny+LAP
        do i=1-LAP,nx+LAP
            ctr(i,j)%conVars=w_uniform
            ctr(i,j)%x=real(i,KREAL)-0.5_KREAL
            ctr(i,j)%y=0.5_KREAL
            ctr(i,j)%area=1.0_KREAL
            ctr(i,j)%length=[1.0_KREAL,1.0_KREAL]
            ctr(i,j)%dt_local=0.1_KREAL
        enddo
    enddo
    vface%length=1.0_KREAL; vface%cosx=1.0_KREAL; vface%cosy=0.0_KREAL
    hface%length=0.0_KREAL; hface%cosx=0.0_KREAL; hface%cosy=1.0_KREAL
    call AllocateImplicitState()
    implicit_ncfl_active=Implicit_nCFL
    rate_w=0.0_KREAL
    rate_w(1,1,1)=0.01_KREAL
    rate_w(1,2,1)=-0.01_KREAL

    gamma_face=ImplicitSpectralRadius(1,1,IDIRC,Rside)
    call AssertTrue(gamma_face>0.0_KREAL,'positive spectral radius')
    defect_before=ComputeMacroDefectNorm()
    call SolveMacroPR2(ok)
    defect_after=ComputeMacroDefectNorm()
    call AssertTrue(ok,'macro solve accepted')
    call AssertTrue(defect_after<defect_before,'PR2 reduces macro defect')
    call BuildPredictedMacroState(alpha,ok)
    call AssertTrue(ok .and. alpha>0.0_KREAL,'physical macro predictor')
    call FinishTests('implicit_macro')
    call MPI_Finalize(ierr_local)
end program test_implicit_macro

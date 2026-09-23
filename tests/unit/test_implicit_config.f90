program test_implicit_config
    use mpi
    use control_para
    use ImportPara, only: Set_default_parameters, IsValidMaxwellAccommodation, &
        IsSupportedImplicitMaxwellAccommodation
    use TestSupport
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer :: ierr_local
    real(KREAL) :: nan_value

    call MPI_Init(ierr_local)
    call Set_default_parameters()
    call AssertTrue(trim(TimeIntegration_method)=='Explicit', 'default mode is Explicit')
    call AssertClose(Implicit_nCFL, 100.0_KREAL, 0.0_KREAL, 'default nCFL')
    call AssertTrue(Implicit_MacroSweeps==2, 'default macro sweep pairs')
    call AssertTrue(Implicit_MicroSweeps==1, 'default micro sweep pairs')
    call AssertClose(Unsteady_Theta,0.5_KREAL,0.0_KREAL,'default unsteady theta')
    call AssertTrue(Unsteady_MaxInner==50,'default unsteady inner limit')
    call AssertTrue(.not.IsImplicitSteady(), 'explicit mode query')
    TimeIntegration_method = 'ImplicitSteady'
    call AssertTrue(IsImplicitSteady(), 'implicit mode query')
    TimeIntegration_method = 'ImplicitUnsteady'
    call AssertTrue(IsImplicitUnsteady(), 'unsteady implicit mode query')
    nan_value=ieee_value(0.0_KREAL,ieee_quiet_nan)
    call AssertTrue(IsValidMaxwellAccommodation(0.0_KREAL), 'ACC=0 is valid')
    call AssertTrue(IsValidMaxwellAccommodation(0.5_KREAL), 'intermediate ACC is valid')
    call AssertTrue(IsValidMaxwellAccommodation(1.0_KREAL), 'ACC=1 is valid')
    call AssertTrue(.not.IsValidMaxwellAccommodation(-epsilon(1.0_KREAL)), 'negative ACC rejected')
    call AssertTrue(.not.IsValidMaxwellAccommodation(1.0_KREAL+epsilon(1.0_KREAL)), 'ACC above one rejected')
    call AssertTrue(.not.IsValidMaxwellAccommodation(nan_value), 'nonfinite ACC rejected')
    call AssertTrue(IsSupportedImplicitMaxwellAccommodation('ImplicitSteady',0.0_KREAL), &
        'steady implicit accepts specular Maxwell wall')
    call AssertTrue(IsSupportedImplicitMaxwellAccommodation('ImplicitSteady',0.5_KREAL), &
        'steady implicit accepts partial Maxwell wall')
    call AssertTrue(IsSupportedImplicitMaxwellAccommodation('ImplicitUnsteady',1.0_KREAL), &
        'unsteady implicit accepts diffuse Maxwell wall')
    call AssertTrue(IsSupportedImplicitMaxwellAccommodation('ImplicitUnsteady',0.6_KREAL), &
        'ImplicitUnsteady accepts partial Maxwell accommodation')
    call AssertTrue(.not.IsSupportedImplicitMaxwellAccommodation('ImplicitUnsteady',-0.1_KREAL), &
        'ImplicitUnsteady rejects ACC below zero')
    call AssertTrue(.not.IsSupportedImplicitMaxwellAccommodation('ImplicitUnsteady',1.1_KREAL), &
        'ImplicitUnsteady rejects ACC above one')
    call FinishTests('implicit_config')
    call MPI_Finalize(ierr_local)
end program test_implicit_config

program test_unsteady_time_scales
    use Precision
    use mod_UnsteadyImplicit, only: CappedFluxTime, ModifiedFluxWeight
    use TestSupport
    implicit none

    call AssertClose(CappedFluxTime(0.2_KREAL,0.5_KREAL),0.2_KREAL,0.0_KREAL, &
        'local CFL interval retained below physical dt')
    call AssertClose(CappedFluxTime(0.8_KREAL,0.5_KREAL),0.5_KREAL,0.0_KREAL, &
        'flux interval capped by physical dt')
    call AssertClose(ModifiedFluxWeight(0.5_KREAL,0.5_KREAL,0.5_KREAL),0.0_KREAL,1.0e-15_KREAL, &
        'Zhu flux weight vanishes in explicit limit')
    call AssertClose(ModifiedFluxWeight(0.1_KREAL,0.5_KREAL,0.5_KREAL),0.4_KREAL,1.0e-15_KREAL, &
        'Zhu modified Crank-Nicolson weight')
    call AssertTrue(ModifiedFluxWeight(0.1_KREAL,0.5_KREAL,0.5_KREAL)>=0.0_KREAL .and. &
                    ModifiedFluxWeight(0.1_KREAL,0.5_KREAL,0.5_KREAL)<=0.5_KREAL, &
        'modified flux weight stays bounded')
    call FinishTests('unsteady_time_scales')
end program test_unsteady_time_scales

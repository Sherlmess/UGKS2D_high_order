program test_unsteady_macro
    use Precision
    use mod_UnsteadyImplicit, only: MacroPhysicalResidual, MacroUnsteadyDiagonal
    use TestSupport
    implicit none
    real(KREAL) :: residual

    residual=MacroPhysicalResidual(1.0_KREAL,1.2_KREAL,0.3_KREAL,0.5_KREAL, &
        0.1_KREAL,0.2_KREAL,0.4_KREAL,0.5_KREAL,0.25_KREAL,0.5_KREAL)
    call AssertClose(residual,0.15_KREAL,1.0e-14_KREAL,'Zhu macro residual algebra')
    call AssertClose(MacroUnsteadyDiagonal(0.5_KREAL,0.25_KREAL,8.0_KREAL,2.0_KREAL), &
        2.5_KREAL,1.0e-14_KREAL,'macro physical-time Rusanov diagonal')
    call FinishTests('unsteady_macro')
end program test_unsteady_macro

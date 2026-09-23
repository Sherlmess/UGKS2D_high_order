program test_unsteady_baseline
    use Precision
    use control_para
    use Flow_DIS
    use Solver, only: AdvancePhysicalTime, NonsteadyStopRequested
    use TestSupport
    implicit none

    time_dis=0.0_KREAL
    call AdvancePhysicalTime(0.1_KREAL)
    call AdvancePhysicalTime(0.2_KREAL)
    call AssertClose(time_dis,0.3_KREAL,1.0e-14_KREAL,'physical time accumulates unequal steps')

    Iflag_DIS=1; iter=2; MAX_ITER=3; res=0.0_KREAL
    call AssertTrue(.not.NonsteadyStopRequested(),'nonsteady run ignores steady residual convergence')
    iter=3
    call AssertTrue(NonsteadyStopRequested(),'nonsteady run stops at configured physical step count')

    print '(A)', 'PASS: unsteady_baseline'
end program test_unsteady_baseline

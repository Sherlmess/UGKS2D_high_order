program test_termination
    use mod_S2SteadyTermination
    implicit none
    integer::last,streak,i
    logical::counted,valid
    last=-1;streak=0
    call SteadyGateObserve(0,.false.,last,streak,counted,valid)
    if(.not.valid.or..not.counted.or.last/=0.or.streak/=0)error stop 1
    do i=1,4
        call SteadyGateObserve(i,.true.,last,streak,counted,valid)
        if(.not.valid.or..not.counted.or.streak/=i)error stop 2
    enddo
    ! A repeated residual evaluation of state 4 cannot create state 5.
    call SteadyGateObserve(4,.true.,last,streak,counted,valid)
    if(.not.valid.or.counted.or.streak/=4)error stop 3
    ! A rejected proposal does not increment the accepted-state number.
    call SteadyGateObserve(4,.true.,last,streak,counted,valid)
    if(counted.or.streak/=4)error stop 4
    ! Four consecutive passes followed by a rebound reset the window.
    call SteadyGateObserve(5,.false.,last,streak,counted,valid)
    if(.not.counted.or.streak/=0)error stop 5
    do i=6,10
        call SteadyGateObserve(i,.true.,last,streak,counted,valid)
    enddo
    if(streak/=5.or.last/=10)error stop 6
    ! An iteration limit immediately after update 11 requires evaluation.
    if(.not.SteadyNeedsEvaluation(11,last))error stop 7
    call SteadyGateObserve(11,.true.,last,streak,counted,valid)
    if(.not.counted.or.SteadyNeedsEvaluation(11,last))error stop 8
    ! State-number regression is rejected, and the tracked state is unchanged.
    call SteadyGateObserve(10,.true.,last,streak,counted,valid)
    if(valid.or.counted.or.last/=11.or.streak/=6)error stop 9
    print '(A)','TERMINATION_CONTROL_PASS'
end program test_termination

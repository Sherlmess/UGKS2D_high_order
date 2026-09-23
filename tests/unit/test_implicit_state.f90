program test_implicit_state
    use mpi
    use flow_data
    use mod_ImplicitState
    use TestSupport
    implicit none
    integer :: ierr_local

    call MPI_Init(ierr_local)
    nx=3; ny=2; Nvel=5
    call AllocateImplicitState()
    call AssertTrue(allocated(rate_w), 'rate_w allocated')
    call AssertTrue(all(shape(rate_w)==[4,3,2]), 'rate_w shape')
    call AssertTrue(lbound(delta_w,2)==1-LAP .and. ubound(delta_w,2)==nx+LAP, 'delta_w x halo')
    call AssertTrue(all(shape(delta_h)==[Nvel,nx+2*LAP,ny+2*LAP]), 'delta_h shape')
    delta_w=1.0_KREAL; delta_h=1.0_KREAL; delta_b=1.0_KREAL
    call ResetImplicitCorrections()
    call AssertTrue(maxval(abs(delta_w))==0.0_KREAL, 'delta_w reset')
    call AssertTrue(maxval(abs(delta_h))==0.0_KREAL, 'delta_h reset')
    call ReleaseImplicitState()
    call AssertTrue(.not.allocated(delta_h), 'workspace released')
    call FinishTests('implicit_state')
    call MPI_Finalize(ierr_local)
end program test_implicit_state

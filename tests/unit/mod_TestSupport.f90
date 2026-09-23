module TestSupport
    use Precision
    implicit none
    integer(KINT) :: failures = 0
contains
    subroutine AssertTrue(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not.condition) then
            failures = failures + 1
            write(*,'(A)') 'FAIL: '//trim(message)
        endif
    end subroutine AssertTrue

    subroutine AssertClose(actual, expected, tolerance, message)
        real(KREAL), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message
        call AssertTrue(abs(actual-expected)<=tolerance, message)
    end subroutine AssertClose

    subroutine FinishTests(name)
        character(len=*), intent(in) :: name
        if (failures/=0) then
            write(*,'(A,I0)') trim(name)//' failures=', failures
            error stop 1
        endif
        write(*,'(A)') 'PASS: '//trim(name)
    end subroutine FinishTests
end module TestSupport

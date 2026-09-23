module mod_S2SteadyTermination
    implicit none
    private
    public :: SteadyGateObserve,SteadyNeedsEvaluation
contains
    pure subroutine SteadyGateObserve(accepted_iter,below,last_iter,consecutive,counted,valid)
        integer,intent(in)::accepted_iter
        logical,intent(in)::below
        integer,intent(inout)::last_iter,consecutive
        logical,intent(out)::counted,valid
        valid=accepted_iter>=last_iter
        counted=.false.
        if(.not.valid.or.accepted_iter==last_iter)return
        last_iter=accepted_iter
        counted=.true.
        if(below)then
            consecutive=consecutive+1
        else
            consecutive=0
        endif
    end subroutine SteadyGateObserve

    pure logical function SteadyNeedsEvaluation(accepted_iter,last_iter)
        integer,intent(in)::accepted_iter,last_iter
        SteadyNeedsEvaluation=accepted_iter/=last_iter
    end function SteadyNeedsEvaluation
end module mod_S2SteadyTermination

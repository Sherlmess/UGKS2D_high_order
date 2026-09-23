! Value/first-direction scalar kernels are ported from the accepted A1 source
! tree at SHA 7d48095697801ad3fb896d19cad557bc1db99fba
! (production source commit 902d7aa330943d8c17d13b5c696ef8ebde3353a4).
module WENO5AOScalar
    use Precision, only: KREAL
    implicit none

    integer, parameter :: rk = KREAL

    type, public :: WENO5AOCellValueCache
        real(rk) :: final_weight(4) = 0.0_rk
        real(rk) :: beta(4) = 0.0_rk
        real(rk) :: coefficient(0:4,4) = 0.0_rk
        logical :: valid = .false.
    end type WENO5AOCellValueCache

    type, public :: WENO5AOCellFirstDirectionCache
        real(rk) :: final_weight(4) = 0.0_rk
        real(rk) :: final_weight_direction(4) = 0.0_rk
        real(rk) :: coefficient(0:4,4,2) = 0.0_rk
        logical :: valid = .false.
    end type WENO5AOCellFirstDirectionCache

contains

    pure integer function pidx(i,n)
        integer, intent(in) :: i,n
        pidx = modulo(i-1,n)+1
    end function pidx

    pure real(rk) function weno5ao_cell_point(q,i,xi) result(qpoint)
        real(rk), intent(in) :: q(:),xi
        integer, intent(in) :: i
        type(WENO5AOCellValueCache) :: cache
        call prepare_weno5ao_cell_value(q,i,cache)
        call evaluate_weno5ao_cell_value(cache,xi,qpoint)
    end function weno5ao_cell_point

    pure real(rk) function reconstruct_weno5ao_cell_point(q,i,xi) result(qpoint)
        real(rk), intent(in) :: q(:),xi
        integer, intent(in) :: i
        qpoint=weno5ao_cell_point(q,i,xi)
    end function reconstruct_weno5ao_cell_point

    pure subroutine prepare_weno5ao_cell_value(q,i,cache)
        real(rk), intent(in) :: q(:)
        integer, intent(in) :: i
        type(WENO5AOCellValueCache), intent(out) :: cache
        integer :: n
        integer :: stencil0(3),stencil1(3),stencil2(3),stencil3(5)
        real(rk) :: local(5),beta(4),tau,alpha(4),omega(4)
        real(rk), parameter :: eps=1.0e-8_rk
        real(rk), parameter :: gamma_hi=0.85_rk,gamma_lo=0.85_rk
        real(rk), parameter :: gamma(4)=[(1.0_rk-gamma_hi)*(1.0_rk-gamma_lo)*0.5_rk, &
            (1.0_rk-gamma_hi)*gamma_lo,(1.0_rk-gamma_hi)*(1.0_rk-gamma_lo)*0.5_rk,gamma_hi]

        cache=WENO5AOCellValueCache()
        n=size(q)
        if (n<5 .or. i<1 .or. i>n) return
        local=[q(pidx(i-2,n)),q(pidx(i-1,n)),q(i),q(pidx(i+1,n)),q(pidx(i+2,n))]
        call weno5ao_smoothness(local,beta)
        cache%beta=beta
        tau=(abs(beta(4)-beta(1))+abs(beta(4)-beta(2))+abs(beta(4)-beta(3)))/3.0_rk
        alpha=gamma*(1.0_rk+(tau/(beta+eps))**2)
        omega=alpha/sum(alpha)
        cache%final_weight(4)=omega(4)/gamma(4)
        cache%final_weight(1:3)=omega(1:3)-cache%final_weight(4)*gamma(1:3)

        stencil0=[-2,-1,0]
        stencil1=[-1,0,1]
        stencil2=[0,1,2]
        stencil3=[-2,-1,0,1,2]
        call cell_average_polynomial_coefficients(q,i,stencil0,cache%coefficient(:,1))
        call cell_average_polynomial_coefficients(q,i,stencil1,cache%coefficient(:,2))
        call cell_average_polynomial_coefficients(q,i,stencil2,cache%coefficient(:,3))
        call cell_average_polynomial_coefficients(q,i,stencil3,cache%coefficient(:,4))
        cache%valid=.true.
    end subroutine prepare_weno5ao_cell_value

    pure subroutine evaluate_weno5ao_cell_value(cache,xi,qpoint)
        type(WENO5AOCellValueCache), intent(in) :: cache
        real(rk), intent(in) :: xi
        real(rk), intent(out) :: qpoint
        real(rk) :: power(0:4),candidate_value(4)
        integer :: degree,candidate
        qpoint=0.0_rk
        if (.not.cache%valid) return
        power(0)=1.0_rk
        do degree=1,4
            power(degree)=power(degree-1)*xi
        end do
        do candidate=1,4
            candidate_value(candidate)=sum(cache%coefficient(:,candidate)*power)
        end do
        qpoint=sum(cache%final_weight*candidate_value)
    end subroutine evaluate_weno5ao_cell_value

    pure subroutine prepare_weno5ao_cell_first_direction(q,q_direction,i,cache,differentiable)
        real(rk), intent(in) :: q(:),q_direction(:)
        integer, intent(in) :: i
        type(WENO5AOCellFirstDirectionCache), intent(out) :: cache
        logical, intent(out) :: differentiable
        integer :: n,candidate
        integer :: stencil0(3),stencil1(3),stencil2(3),stencil3(5)
        real(rk) :: local(5),local_direction(5),beta(4),beta_direction(4),denominator(4)
        real(rk) :: delta(3),delta_direction(3),branch_sign,branch_tolerance,path_scale
        real(rk) :: tau,tau_direction,ratio(4),ratio_direction(4)
        real(rk) :: alpha(4),alpha_direction(4),sum_alpha,sum_alpha_direction
        real(rk) :: omega(4),omega_direction(4)
        real(rk), parameter :: eps=1.0e-8_rk
        real(rk), parameter :: gamma_hi=0.85_rk,gamma_lo=0.85_rk
        real(rk), parameter :: gamma(4)=[(1.0_rk-gamma_hi)*(1.0_rk-gamma_lo)*0.5_rk, &
            (1.0_rk-gamma_hi)*gamma_lo,(1.0_rk-gamma_hi)*(1.0_rk-gamma_lo)*0.5_rk,gamma_hi]

        cache=WENO5AOCellFirstDirectionCache()
        differentiable=size(q)>=5 .and. size(q_direction)==size(q) .and. i>=1 .and. i<=size(q)
        if (.not.differentiable) return
        n=size(q)
        local=[q(pidx(i-2,n)),q(pidx(i-1,n)),q(i),q(pidx(i+1,n)),q(pidx(i+2,n))]
        local_direction=[q_direction(pidx(i-2,n)),q_direction(pidx(i-1,n)),q_direction(i), &
            q_direction(pidx(i+1,n)),q_direction(pidx(i+2,n))]
        call weno5ao_smoothness_first_direction(local,local_direction,beta,beta_direction)
        delta=beta(4)-beta(1:3)
        delta_direction=beta_direction(4)-beta_direction(1:3)
        path_scale=max(maxval(abs(local))+maxval(abs(local_direction)),sqrt(tiny(1.0_rk)))
        branch_tolerance=1024.0_rk*epsilon(1.0_rk)*path_scale**2
        tau=0.0_rk
        tau_direction=0.0_rk
        do candidate=1,3
            if (abs(delta(candidate))>tiny(1.0_rk)) then
                branch_sign=sign(1.0_rk,delta(candidate))
                tau=tau+abs(delta(candidate))
                tau_direction=tau_direction+branch_sign*delta_direction(candidate)
            else if (abs(delta_direction(candidate))>branch_tolerance) then
                differentiable=.false.
                return
            end if
        end do
        tau=tau/3.0_rk
        tau_direction=tau_direction/3.0_rk
        denominator=beta+eps
        ratio=tau/denominator
        ratio_direction=tau_direction/denominator-tau*beta_direction/denominator**2
        alpha=gamma*(1.0_rk+ratio**2)
        alpha_direction=2.0_rk*gamma*ratio*ratio_direction
        sum_alpha=sum(alpha)
        sum_alpha_direction=sum(alpha_direction)
        omega=alpha/sum_alpha
        omega_direction=alpha_direction/sum_alpha-alpha*sum_alpha_direction/sum_alpha**2
        cache%final_weight(4)=omega(4)/gamma(4)
        cache%final_weight_direction(4)=omega_direction(4)/gamma(4)
        cache%final_weight(1:3)=omega(1:3)-cache%final_weight(4)*gamma(1:3)
        cache%final_weight_direction(1:3)=omega_direction(1:3)-cache%final_weight_direction(4)*gamma(1:3)
        cache%final_weight(4)=1.0_rk-sum(cache%final_weight(1:3))
        cache%final_weight_direction(4)=-sum(cache%final_weight_direction(1:3))

        stencil0=[-2,-1,0]
        stencil1=[-1,0,1]
        stencil2=[0,1,2]
        stencil3=[-2,-1,0,1,2]
        call cell_average_polynomial_coefficients(q,i,stencil0,cache%coefficient(:,1,1))
        call cell_average_polynomial_coefficients(q_direction,i,stencil0,cache%coefficient(:,1,2))
        call cell_average_polynomial_coefficients(q,i,stencil1,cache%coefficient(:,2,1))
        call cell_average_polynomial_coefficients(q_direction,i,stencil1,cache%coefficient(:,2,2))
        call cell_average_polynomial_coefficients(q,i,stencil2,cache%coefficient(:,3,1))
        call cell_average_polynomial_coefficients(q_direction,i,stencil2,cache%coefficient(:,3,2))
        call cell_average_polynomial_coefficients(q,i,stencil3,cache%coefficient(:,4,1))
        call cell_average_polynomial_coefficients(q_direction,i,stencil3,cache%coefficient(:,4,2))
        cache%valid=.true.
    end subroutine prepare_weno5ao_cell_first_direction

    pure subroutine evaluate_weno5ao_cell_first_direction(cache,xi,qpoint,qpoint_direction)
        type(WENO5AOCellFirstDirectionCache), intent(in) :: cache
        real(rk), intent(in) :: xi
        real(rk), intent(out) :: qpoint,qpoint_direction
        real(rk) :: power(0:4),p(4),p_direction(4)
        integer :: degree,candidate
        qpoint=0.0_rk
        qpoint_direction=0.0_rk
        if (.not.cache%valid) return
        power(0)=1.0_rk
        do degree=1,4
            power(degree)=power(degree-1)*xi
        end do
        do candidate=1,4
            p(candidate)=sum(cache%coefficient(:,candidate,1)*power)
            p_direction(candidate)=sum(cache%coefficient(:,candidate,2)*power)
        end do
        qpoint=sum(cache%final_weight*p)
        qpoint_direction=sum(cache%final_weight_direction*p+cache%final_weight*p_direction)
    end subroutine evaluate_weno5ao_cell_first_direction

    pure subroutine weno5ao_smoothness_first_direction(local,local_direction,beta,beta_direction)
        real(rk), intent(in) :: local(5),local_direction(5)
        real(rk), intent(out) :: beta(4),beta_direction(4)
        real(rk) :: beta_direction_quadratic(4),beta_sum(4)
        call weno5ao_smoothness(local,beta)
        call weno5ao_smoothness(local_direction,beta_direction_quadratic)
        call weno5ao_smoothness(local+local_direction,beta_sum)
        beta_direction=beta_sum-beta-beta_direction_quadratic
    end subroutine weno5ao_smoothness_first_direction

    pure subroutine weno5ao_smoothness(local,beta)
        real(rk), intent(in) :: local(5)
        real(rk), intent(out) :: beta(4)
        real(rk) :: qm2,qm1,qi,qp1,qp2
        qm2=local(1); qm1=local(2); qi=local(3); qp1=local(4); qp2=local(5)
        beta(1)=(13.0_rk/12.0_rk)*(qm2-2.0_rk*qm1+qi)**2+0.25_rk*(qm2-4.0_rk*qm1+3.0_rk*qi)**2
        beta(2)=(13.0_rk/12.0_rk)*(qm1-2.0_rk*qi+qp1)**2+0.25_rk*(qm1-qp1)**2
        beta(3)=(13.0_rk/12.0_rk)*(qi-2.0_rk*qp1+qp2)**2+0.25_rk*(3.0_rk*qi-4.0_rk*qp1+qp2)**2
        beta(4)=whole_stencil_beta(qm2,qm1,qi,qp1,qp2)
    end subroutine weno5ao_smoothness

    pure real(rk) function whole_stencil_beta(qm2,qm1,qi,qp1,qp2) result(beta)
        real(rk), intent(in) :: qm2,qm1,qi,qp1,qp2
        beta=(231153.0_rk*qi*qi+104963.0_rk*qm1*qm1+6908.0_rk*qm2*qm2-38947.0_rk*qm2*qp1+ &
            104963.0_rk*qp1*qp1+qm1*(-51001.0_rk*qm2+179098.0_rk*qp1-38947.0_rk*qp2)- &
            3.0_rk*qi*(99692.0_rk*qm1-22641.0_rk*qm2+99692.0_rk*qp1-22641.0_rk*qp2)+ &
            8209.0_rk*qm2*qp2-51001.0_rk*qp1*qp2+6908.0_rk*qp2*qp2)/5040.0_rk
    end function whole_stencil_beta

    pure subroutine cell_average_polynomial_coefficients(q,i,offsets,coefficient)
        real(rk), intent(in) :: q(:)
        integer, intent(in) :: i,offsets(:)
        real(rk), intent(out) :: coefficient(0:)
        real(rk) :: matrix(size(offsets),size(offsets)),rhs(size(offsets)),solution(size(offsets))
        integer :: row,degree,n
        n=size(q)
        coefficient=0.0_rk
        do row=1,size(offsets)
            rhs(row)=q(pidx(i+offsets(row),n))
            do degree=0,size(offsets)-1
                matrix(row,degree+1)=((real(offsets(row),rk)+0.5_rk)**(degree+1)- &
                    (real(offsets(row),rk)-0.5_rk)**(degree+1))/real(degree+1,rk)
            end do
        end do
        call solve_small_linear_system(matrix,rhs,solution)
        coefficient(0:size(offsets)-1)=solution
    end subroutine cell_average_polynomial_coefficients

    pure subroutine solve_small_linear_system(matrix,rhs,solution)
        real(rk), intent(in) :: matrix(:,:),rhs(:)
        real(rk), intent(out) :: solution(:)
        real(rk) :: work(size(rhs),size(rhs)),work_rhs(size(rhs)),factor,pivot_value
        integer :: i,pivot,row,n
        n=size(rhs)
        work=matrix
        work_rhs=rhs
        do pivot=1,n-1
            row=pivot-1+maxloc(abs(work(pivot:n,pivot)),dim=1)
            if (row/=pivot) then
                call swap_small_system_rows(work,pivot,row)
                pivot_value=work_rhs(pivot); work_rhs(pivot)=work_rhs(row); work_rhs(row)=pivot_value
            end if
            do i=pivot+1,n
                factor=work(i,pivot)/work(pivot,pivot)
                work(i,pivot:n)=work(i,pivot:n)-factor*work(pivot,pivot:n)
                work_rhs(i)=work_rhs(i)-factor*work_rhs(pivot)
            end do
        end do
        solution=0.0_rk
        do i=n,1,-1
            solution(i)=(work_rhs(i)-sum(work(i,i+1:n)*solution(i+1:n)))/work(i,i)
        end do
    end subroutine solve_small_linear_system

    pure subroutine swap_small_system_rows(matrix,first,second)
        real(rk), intent(inout) :: matrix(:,:)
        integer, intent(in) :: first,second
        real(rk) :: temporary(size(matrix,2))
        temporary=matrix(first,:)
        matrix(first,:)=matrix(second,:)
        matrix(second,:)=temporary
    end subroutine swap_small_system_rows

end module WENO5AOScalar



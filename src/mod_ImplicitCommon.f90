module mod_ImplicitCommon
    use Precision, only: KREAL
    use flow_para, only: gamma
    use, intrinsic :: ieee_arithmetic
    implicit none
    private
    public :: IsPhysicalState, EulerFluxNormal, MomentConserved, PseudoTimeStep
    public :: ProjectToConservedMoments
contains
    pure logical function IsPhysicalState(w)
        real(KREAL), intent(in) :: w(4)
        real(KREAL) :: internal_energy

        if (.not.all(ieee_is_finite(w)) .or. w(1)<=0.0_KREAL) then
            IsPhysicalState=.false.
            return
        endif
        internal_energy=w(4)-0.5_KREAL*(w(2)*w(2)+w(3)*w(3))/w(1)
        IsPhysicalState=ieee_is_finite(internal_energy) .and. internal_energy>0.0_KREAL
    end function IsPhysicalState

    subroutine EulerFluxNormal(w,nx_face,ny_face,flux)
        real(KREAL), intent(in) :: w(4),nx_face,ny_face
        real(KREAL), intent(out) :: flux(4)
        real(KREAL) :: u,v,un,pressure

        u=w(2)/w(1); v=w(3)/w(1); un=u*nx_face+v*ny_face
        pressure=(gamma-1.0_KREAL)*(w(4)-0.5_KREAL*(w(2)*w(2)+w(3)*w(3))/w(1))
        flux=[w(1)*un,w(2)*un+pressure*nx_face,w(3)*un+pressure*ny_face,(w(4)+pressure)*un]
    end subroutine EulerFluxNormal

    subroutine MomentConserved(h,b,u,v,weight_local,w)
        real(KREAL), intent(in) :: h(:),b(:),u(:),v(:),weight_local(:)
        real(KREAL), intent(out) :: w(4)

        w(1)=sum(weight_local*h)
        w(2)=sum(weight_local*u*h)
        w(3)=sum(weight_local*v*h)
        w(4)=0.5_KREAL*sum(weight_local*((u*u+v*v)*h+b))
    end subroutine MomentConserved

    subroutine ProjectToConservedMoments(h,b,target,u,v,weight_local,reference_h,reference_b,ok)
        real(KREAL), intent(inout) :: h(:),b(:)
        real(KREAL), intent(in) :: target(4),u(:),v(:),weight_local(:),reference_h(:),reference_b(:)
        logical, intent(out) :: ok
        real(KREAL) :: gram(4,4),rhs(4),lambda(4),moment(4),kernel_h(4),kernel_b(4)
        real(KREAL) :: scale_h,scale_b,weight_h,weight_b,residual_scale
        integer :: k,a,c,n

        n=size(h); ok=.false.
        if (size(b)/=n .or. size(u)/=n .or. size(v)/=n .or. size(weight_local)/=n .or. &
            size(reference_h)/=n .or. size(reference_b)/=n) return
        if (.not.all(ieee_is_finite(target)) .or. .not.all(ieee_is_finite(h)) .or. &
            .not.all(ieee_is_finite(b))) return
        call MomentConserved(h,b,u,v,weight_local,moment)
        rhs=target-moment; gram=0.0_KREAL
        scale_h=max(maxval(abs(reference_h)),1.0_KREAL)
        scale_b=max(maxval(abs(reference_b)),1.0_KREAL)
        kernel_b=[0.0_KREAL,0.0_KREAL,0.0_KREAL,0.5_KREAL]
        do k=1,n
            kernel_h=[1.0_KREAL,u(k),v(k),0.5_KREAL*(u(k)*u(k)+v(k)*v(k))]
            weight_h=max(abs(reference_h(k)),epsilon(1.0_KREAL)*scale_h)
            weight_b=max(abs(reference_b(k)),epsilon(1.0_KREAL)*scale_b)
            do a=1,4; do c=1,4
                gram(a,c)=gram(a,c)+weight_local(k)*(weight_h*kernel_h(a)*kernel_h(c)+ &
                    weight_b*kernel_b(a)*kernel_b(c))
            enddo; enddo
        enddo
        call SolveDense4(gram,rhs,lambda,ok)
        if (.not.ok) return
        do k=1,n
            kernel_h=[1.0_KREAL,u(k),v(k),0.5_KREAL*(u(k)*u(k)+v(k)*v(k))]
            weight_h=max(abs(reference_h(k)),epsilon(1.0_KREAL)*scale_h)
            weight_b=max(abs(reference_b(k)),epsilon(1.0_KREAL)*scale_b)
            h(k)=h(k)+weight_h*dot_product(kernel_h,lambda)
            b(k)=b(k)+0.5_KREAL*weight_b*lambda(4)
        enddo
        call MomentConserved(h,b,u,v,weight_local,moment)
        residual_scale=max(maxval(abs(target)),1.0_KREAL)
        ok=all(ieee_is_finite(h)) .and. all(ieee_is_finite(b)) .and. &
            maxval(abs(moment-target))<=1.0e-10_KREAL*residual_scale
    end subroutine ProjectToConservedMoments

    subroutine SolveDense4(matrix,rhs,solution,ok)
        real(KREAL), intent(in) :: matrix(4,4),rhs(4)
        real(KREAL), intent(out) :: solution(4)
        logical, intent(out) :: ok
        real(KREAL) :: augmented(4,5),row_swap(5),factor,pivot_scale
        integer :: column,pivot,row

        augmented(:,1:4)=matrix; augmented(:,5)=rhs
        pivot_scale=max(maxval(abs(matrix)),tiny(1.0_KREAL)); ok=.false.
        do column=1,4
            pivot=column
            do row=column+1,4
                if (abs(augmented(row,column))>abs(augmented(pivot,column))) pivot=row
            enddo
            if (abs(augmented(pivot,column))<=1000.0_KREAL*epsilon(1.0_KREAL)*pivot_scale) return
            if (pivot/=column) then
                row_swap=augmented(column,:)
                augmented(column,:)=augmented(pivot,:)
                augmented(pivot,:)=row_swap
            endif
            augmented(column,:)=augmented(column,:)/augmented(column,column)
            do row=1,4
                if (row==column) cycle
                factor=augmented(row,column)
                augmented(row,:)=augmented(row,:)-factor*augmented(column,:)
            enddo
        enddo
        solution=augmented(:,5)
        ok=all(ieee_is_finite(solution))
    end subroutine SolveDense4

    pure real(KREAL) function PseudoTimeStep(dt_phys,ncfl,pcfl)
        real(KREAL), intent(in) :: dt_phys,ncfl,pcfl

        PseudoTimeStep=dt_phys*ncfl/pcfl
    end function PseudoTimeStep
end module mod_ImplicitCommon

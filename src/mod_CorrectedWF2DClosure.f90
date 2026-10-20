module mod_CorrectedWF2DClosure
    use Precision, only: KINT, KREAL
    use mod_CorrectedWF2DEquilibrium, only: primitive_from_conserved_2d, &
        analytic_maxwellian_2d, analytic_maxwellian_direction_2d
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: moment_2d, discrete_equilibrium_q_2d, dq_2d
    public :: projector_2d
    public :: ClosureFirstOrderContext,closure_first_prepare,closure_first_direction
    public :: closure_first_project,closure_first_matches,closure_first_copy

    type :: ClosureFirstOrderContext
        private
        logical :: valid=.false.
        integer(KINT) :: ck=-1
        real(KREAL) :: w(4)=0.0_KREAL,theta(4)=0.0_KREAL,jac(4,4)=0.0_KREAL
        real(KREAL),allocatable :: u(:),v(:),weight(:),qh(:),qb(:),qht(:,:),qbt(:,:)
    end type ClosureFirstOrderContext

contains

    subroutine moment_2d(h,b,uvel,vvel,weight,m,ok)
        real(KREAL), intent(in) :: h(:),b(:),uvel(:),vvel(:),weight(:)
        real(KREAL), intent(out) :: m(4)
        logical, intent(out) :: ok
        integer :: k,n

        m=0.0_KREAL; ok=.false.
        n=size(h)
        if (size(b)/=n .or. size(uvel)/=n .or. size(vvel)/=n .or. size(weight)/=n) return
        if (.not.all(ieee_is_finite(h)) .or. .not.all(ieee_is_finite(b)) .or. &
            .not.all(ieee_is_finite(uvel)) .or. .not.all(ieee_is_finite(vvel)) .or. &
            .not.all(ieee_is_finite(weight)) .or. any(weight<=0.0_KREAL)) return
        do k=1,n
            m(1)=m(1)+weight(k)*h(k)
            m(2)=m(2)+weight(k)*uvel(k)*h(k)
            m(3)=m(3)+weight(k)*vvel(k)*h(k)
            m(4)=m(4)+0.5_KREAL*weight(k)*((uvel(k)**2+vvel(k)**2)*h(k)+b(k))
        end do
        ok=all(ieee_is_finite(m))
    end subroutine moment_2d

    ! One quadrature point owns one short-lived context.  Evaluation takes no
    ! external W/velocity table, so a different state cannot silently use it.
    subroutine closure_first_invalidate(ctx)
        type(ClosureFirstOrderContext),intent(inout)::ctx
        ctx%valid=.false.;ctx%ck=-1;ctx%w=0.0_KREAL;ctx%theta=0.0_KREAL;ctx%jac=0.0_KREAL
        if(allocated(ctx%u))deallocate(ctx%u)
        if(allocated(ctx%v))deallocate(ctx%v)
        if(allocated(ctx%weight))deallocate(ctx%weight)
        if(allocated(ctx%qh))deallocate(ctx%qh)
        if(allocated(ctx%qb))deallocate(ctx%qb)
        if(allocated(ctx%qht))deallocate(ctx%qht)
        if(allocated(ctx%qbt))deallocate(ctx%qbt)
    end subroutine closure_first_invalidate

    subroutine closure_first_prepare(ctx,w,ck,uvel,vvel,weight,ok)
        type(ClosureFirstOrderContext),intent(inout)::ctx
        real(KREAL),intent(in)::w(4),uvel(:),vvel(:),weight(:)
        integer(KINT),intent(in)::ck
        logical,intent(out)::ok
        integer::n
        call closure_first_invalidate(ctx);ok=.false.;n=size(uvel)
        if(n<1.or.ck<0.or.size(vvel)/=n.or.size(weight)/=n)return
        if(.not.all(ieee_is_finite(uvel)).or..not.all(ieee_is_finite(vvel)).or. &
            .not.all(ieee_is_finite(weight)).or.any(weight<=0.0_KREAL))return
        allocate(ctx%u(n),ctx%v(n),ctx%weight(n),ctx%qh(n),ctx%qb(n),ctx%qht(n,4),ctx%qbt(n,4))
        ctx%w=w;ctx%ck=ck;ctx%u=uvel;ctx%v=vvel;ctx%weight=weight
        call solve_q_first(w,ck,ctx%u,ctx%v,ctx%weight,ctx%qh,ctx%qb,ctx%qht,ctx%qbt,ctx%theta,ctx%jac,ok)
        if(.not.ok)then
            call closure_first_invalidate(ctx)
            return
        endif
        ctx%valid=.true.
    end subroutine closure_first_prepare

    subroutine closure_first_matches(ctx,w,ck,uvel,vvel,weight,matched)
        type(ClosureFirstOrderContext),intent(in)::ctx
        real(KREAL),intent(in)::w(4),uvel(:),vvel(:),weight(:)
        integer(KINT),intent(in)::ck
        logical,intent(out)::matched
        integer::n
        matched=.false.
        if(.not.ctx%valid)return
        n=size(ctx%u)
        if(size(uvel)/=n.or.size(vvel)/=n.or.size(weight)/=n.or.ck/=ctx%ck)return
        matched=all(w==ctx%w).and.all(uvel==ctx%u).and.all(vvel==ctx%v).and.all(weight==ctx%weight)
    end subroutine closure_first_matches

    subroutine closure_first_copy(ctx,qh,qb,theta,jac,qht,qbt,ok)
        type(ClosureFirstOrderContext),intent(in)::ctx
        real(KREAL),intent(out)::qh(:),qb(:),theta(4),jac(4,4),qht(:,:),qbt(:,:)
        logical,intent(out)::ok
        qh=0.0_KREAL;qb=0.0_KREAL;theta=0.0_KREAL;jac=0.0_KREAL;qht=0.0_KREAL;qbt=0.0_KREAL;ok=.false.
        if(.not.ctx%valid)return
        if(size(qh)/=size(ctx%qh).or.size(qb)/=size(ctx%qb).or. &
            any(shape(qht)/=shape(ctx%qht)).or.any(shape(qbt)/=shape(ctx%qbt)))return
        qh=ctx%qh;qb=ctx%qb;theta=ctx%theta;jac=ctx%jac;qht=ctx%qht;qbt=ctx%qbt;ok=.true.
    end subroutine closure_first_copy

    subroutine closure_first_direction(ctx,m_direction,dq_h,dq_b,ok)
        type(ClosureFirstOrderContext),intent(in)::ctx
        real(KREAL),intent(in)::m_direction(4)
        real(KREAL),intent(out)::dq_h(:),dq_b(:)
        logical,intent(out)::ok
        real(KREAL)::theta_direction(4)
        integer::a
        dq_h=0.0_KREAL;dq_b=0.0_KREAL;ok=.false.
        if(.not.ctx%valid)return
        if(size(dq_h)/=size(ctx%qh).or.size(dq_b)/=size(ctx%qb))return
        call solve_dense_4x4(ctx%jac,m_direction,theta_direction,ok)
        if(.not.ok)return
        do a=1,4
            dq_h=dq_h+ctx%qht(:,a)*theta_direction(a)
            dq_b=dq_b+ctx%qbt(:,a)*theta_direction(a)
        enddo
        ok=all(ieee_is_finite(dq_h)).and.all(ieee_is_finite(dq_b))
    end subroutine closure_first_direction

    subroutine closure_first_project(ctx,z_h,z_b,p_h,p_b,ok)
        type(ClosureFirstOrderContext),intent(in)::ctx
        real(KREAL),intent(in)::z_h(:),z_b(:)
        real(KREAL),intent(out)::p_h(:),p_b(:)
        logical,intent(out)::ok
        real(KREAL)::m(4)
        real(KREAL),allocatable::ch(:),cb(:)
        logical::moment_ok
        p_h=0.0_KREAL;p_b=0.0_KREAL;ok=.false.
        if(.not.ctx%valid)return
        if(size(z_h)/=size(ctx%qh).or.size(z_b)/=size(ctx%qb).or. &
            size(p_h)/=size(z_h).or.size(p_b)/=size(z_h))return
        call moment_2d(z_h,z_b,ctx%u,ctx%v,ctx%weight,m,moment_ok)
        if(.not.moment_ok)return
        allocate(ch(size(z_h)),cb(size(z_h)))
        call closure_first_direction(ctx,m,ch,cb,ok)
        if(ok)then
            p_h=z_h-ch;p_b=z_b-cb
            ok=all(ieee_is_finite(p_h)).and.all(ieee_is_finite(p_b))
        endif
        deallocate(ch,cb)
    end subroutine closure_first_project

    subroutine discrete_equilibrium_q_2d(w,ck,uvel,vvel,weight,q_h,q_b,theta,jac,ok)
        real(KREAL), intent(in) :: w(4),uvel(:),vvel(:),weight(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: q_h(:),q_b(:),theta(4),jac(4,4)
        logical, intent(out) :: ok
        real(KREAL), allocatable :: qht(:,:),qbt(:,:)

        allocate(qht(size(uvel),4),qbt(size(uvel),4))
        call solve_q_first(w,ck,uvel,vvel,weight,q_h,q_b,qht,qbt,theta,jac,ok)
        deallocate(qht,qbt)
    end subroutine discrete_equilibrium_q_2d

    subroutine dq_2d(w,m_direction,ck,uvel,vvel,weight,dq_h,dq_b,ok)
        real(KREAL), intent(in) :: w(4),m_direction(4),uvel(:),vvel(:),weight(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: dq_h(:),dq_b(:)
        logical, intent(out) :: ok
        integer :: n,a
        real(KREAL) :: theta(4),jac(4,4),theta_direction(4)
        real(KREAL), allocatable :: qh(:),qb(:),qh_t(:,:),qb_t(:,:)

        n=size(uvel); dq_h=0.0_KREAL; dq_b=0.0_KREAL; ok=.false.
        if (size(dq_h)/=n .or. size(dq_b)/=n .or. size(vvel)/=n .or. size(weight)/=n) return
        allocate(qh(n),qb(n),qh_t(n,4),qb_t(n,4))
        call solve_q_first(w,ck,uvel,vvel,weight,qh,qb,qh_t,qb_t,theta,jac,ok)
        if (.not.ok) then
            deallocate(qh,qb,qh_t,qb_t)
            return
        end if
        call solve_dense_4x4(jac,m_direction,theta_direction,ok)
        if (.not.ok) then
            deallocate(qh,qb,qh_t,qb_t)
            return
        end if
        do a=1,4
            dq_h=dq_h+qh_t(:,a)*theta_direction(a)
            dq_b=dq_b+qb_t(:,a)*theta_direction(a)
        end do
        ok=all(ieee_is_finite(dq_h)) .and. all(ieee_is_finite(dq_b))
        deallocate(qh,qb,qh_t,qb_t)
    end subroutine dq_2d

    subroutine projector_2d(w,z_h,z_b,ck,uvel,vvel,weight,p_h,p_b,ok)
        real(KREAL), intent(in) :: w(4),z_h(:),z_b(:),uvel(:),vvel(:),weight(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: p_h(:),p_b(:)
        logical, intent(out) :: ok
        real(KREAL) :: m(4)
        real(KREAL), allocatable :: correction_h(:),correction_b(:)
        logical :: moment_ok

        p_h=0.0_KREAL; p_b=0.0_KREAL; ok=.false.
        if (size(z_b)/=size(z_h) .or. size(p_h)/=size(z_h) .or. size(p_b)/=size(z_h)) return
        call moment_2d(z_h,z_b,uvel,vvel,weight,m,moment_ok)
        if (.not.moment_ok) return
        allocate(correction_h(size(z_h)),correction_b(size(z_h)))
        call dq_2d(w,m,ck,uvel,vvel,weight,correction_h,correction_b,ok)
        if (ok) then
            p_h=z_h-correction_h
            p_b=z_b-correction_b
            ok=all(ieee_is_finite(p_h)) .and. all(ieee_is_finite(p_b))
        end if
        deallocate(correction_h,correction_b)
    end subroutine projector_2d




    ! First-order preparation keeps the legacy Newton steps and acceptance rule,
    ! but never forms the velocity-by-4-by-4 Hessian used only by tangent APIs.
    subroutine solve_q_first(w,ck,uvel,vvel,weight,qh,qb,qh_theta,qb_theta,theta,jac,ok)
        real(KREAL), intent(in) :: w(4),uvel(:),vvel(:),weight(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: qh(:),qb(:),qh_theta(:,:),qb_theta(:,:),theta(4),jac(4,4)
        logical, intent(out) :: ok
        real(KREAL) :: rho,u,v,temp,residual(4),new_theta(4),new_jac(4,4),new_residual(4)
        real(KREAL) :: step(4),norm_old,norm_new,alpha
        real(KREAL), allocatable :: new_qh(:),new_qb(:),new_qht(:,:),new_qbt(:,:)
        logical :: primitive_ok,solve_ok,accepted
        integer :: iter
        real(KREAL), parameter :: tolerance=5.0e-13_KREAL,backtracking_floor=2.0_KREAL**(-20)

        qh=0.0_KREAL; qb=0.0_KREAL; qh_theta=0.0_KREAL; qb_theta=0.0_KREAL
        theta=0.0_KREAL; jac=0.0_KREAL; ok=.false.
        if (ck<0 .or. size(qb)/=size(qh) .or. size(qh_theta,1)/=size(qh) .or. size(qh_theta,2)/=4 .or. &
            size(qb_theta,1)/=size(qh) .or. size(qb_theta,2)/=4) return
        call primitive_from_conserved_2d(w,ck,rho,u,v,temp,primitive_ok)
        if (.not.primitive_ok) return
        theta=[rho,u,v,temp]
        call q_state_first(theta,ck,uvel,vvel,weight,qh,qb,qh_theta,qb_theta,jac,solve_ok)
        if (.not.solve_ok) return
        call moment_2d(qh,qb,uvel,vvel,weight,residual,solve_ok)
        if (.not.solve_ok) return
        residual=residual-w
        norm_old=scaled_residual(residual,w)
        do iter=1,12
            if (norm_old<=tolerance) exit
            call solve_dense_4x4(jac,-residual,step,solve_ok)
            if (.not.solve_ok) then
                ok=.false.
                return
            end if
            alpha=1.0_KREAL
            accepted=.false.
            do while (alpha>=backtracking_floor)
                new_theta=theta+alpha*step
                if (new_theta(1)>0.0_KREAL .and. new_theta(4)>0.0_KREAL .and. &
                    all(ieee_is_finite(new_theta))) then
                    allocate(new_qh(size(qh)),new_qb(size(qh)),new_qht(size(qh),4),new_qbt(size(qh),4))
                    call q_state_first(new_theta,ck,uvel,vvel,weight,new_qh,new_qb,new_qht,new_qbt,new_jac,solve_ok)
                    if (solve_ok) then
                        call moment_2d(new_qh,new_qb,uvel,vvel,weight,new_residual,solve_ok)
                        if (solve_ok) then
                            new_residual=new_residual-w
                            norm_new=scaled_residual(new_residual,w)
                            if (norm_new<norm_old .or. norm_new<=tolerance) accepted=.true.
                        end if
                    end if
                    if (accepted) then
                        theta=new_theta; qh=new_qh; qb=new_qb; qh_theta=new_qht; qb_theta=new_qbt
                        jac=new_jac; residual=new_residual; norm_old=norm_new
                    end if
                    deallocate(new_qh,new_qb,new_qht,new_qbt)
                    if (accepted) exit
                end if
                alpha=0.5_KREAL*alpha
            end do
            if (.not.accepted) then
                ok=.false.
                return
            end if
        end do
        ok=norm_old<=tolerance .and. all(ieee_is_finite(qh)) .and. all(ieee_is_finite(qb)) .and. &
            all(qh>0.0_KREAL) .and. all(qb>0.0_KREAL)
    end subroutine solve_q_first

    subroutine q_state_first(theta,ck,uvel,vvel,weight,qh,qb,qh_theta,qb_theta,jac,ok)
        real(KREAL), intent(in) :: theta(4),uvel(:),vvel(:),weight(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: qh(:),qb(:),qh_theta(:,:),qb_theta(:,:),jac(4,4)
        logical, intent(out) :: ok
        real(KREAL) :: logd(4),logdb(4),du,dv,r,theta_width,c2,m(4)
        integer :: k,a
        logical :: moment_ok

        qh=0.0_KREAL; qb=0.0_KREAL; qh_theta=0.0_KREAL; qb_theta=0.0_KREAL
        jac=0.0_KREAL; ok=.false.
        if (theta(1)<=0.0_KREAL .or. theta(4)<=0.0_KREAL .or. size(qb)/=size(qh) .or. &
            size(vvel)/=size(qh) .or. size(weight)/=size(qh)) return
        r=theta(1); theta_width=theta(4)
        do k=1,size(qh)
            du=uvel(k)-theta(2); dv=vvel(k)-theta(3); c2=du*du+dv*dv
            qh(k)=r/(acos(-1.0_KREAL)*theta_width)*exp(-c2/theta_width)
            logd=[1.0_KREAL/r,2.0_KREAL*du/theta_width,2.0_KREAL*dv/theta_width, &
                -1.0_KREAL/theta_width+c2/theta_width**2]
            do a=1,4
                qh_theta(k,a)=qh(k)*logd(a)
            end do
            if (ck>0) then
                qb(k)=0.5_KREAL*real(ck,KREAL)*theta_width*qh(k)
                logdb=logd; logdb(4)=logdb(4)+1.0_KREAL/theta_width
                do a=1,4
                    qb_theta(k,a)=qb(k)*logdb(a)
                end do
            end if
        end do
        do a=1,4
            call moment_2d(qh_theta(:,a),qb_theta(:,a),uvel,vvel,weight,m,moment_ok)
            if (.not.moment_ok) return
            jac(:,a)=m
        end do
        ok=all(ieee_is_finite(qh)) .and. all(ieee_is_finite(qb)) .and. all(ieee_is_finite(jac))
    end subroutine q_state_first




    real(KREAL) function scaled_residual(residual,target)
        real(KREAL), intent(in) :: residual(4),target(4)
        scaled_residual=maxval(abs(residual)/max(1.0_KREAL,abs(target)))
    end function scaled_residual

    subroutine solve_dense_4x4(matrix,rhs,solution,ok)
        real(KREAL), intent(in) :: matrix(4,4),rhs(4)
        real(KREAL), intent(out) :: solution(4)
        logical, intent(out) :: ok
        real(KREAL) :: a(4,4),b(4),factor,pivot_value,tmp
        integer :: p,row,i,j

        a=matrix; b=rhs; solution=0.0_KREAL; ok=.false.
        do p=1,3
            row=p-1+maxloc(abs(a(p:4,p)),dim=1)
            if (abs(a(row,p))<=100.0_KREAL*tiny(1.0_KREAL) .or. .not.ieee_is_finite(a(row,p))) return
            if (row/=p) then
                do j=1,4
                    tmp=a(p,j); a(p,j)=a(row,j); a(row,j)=tmp
                end do
                tmp=b(p); b(p)=b(row); b(row)=tmp
            end if
            do i=p+1,4
                factor=a(i,p)/a(p,p)
                a(i,p:4)=a(i,p:4)-factor*a(p,p:4)
                b(i)=b(i)-factor*b(p)
            end do
        end do
        if (abs(a(4,4))<=100.0_KREAL*tiny(1.0_KREAL) .or. .not.ieee_is_finite(a(4,4))) return
        do i=4,1,-1
            pivot_value=b(i)
            if (i<4) pivot_value=pivot_value-sum(a(i,i+1:4)*solution(i+1:4))
            solution(i)=pivot_value/a(i,i)
        end do
        ok=all(ieee_is_finite(solution))
    end subroutine solve_dense_4x4

end module mod_CorrectedWF2DClosure

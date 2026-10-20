module mod_CorrectedWF2DEquilibrium
    use Precision, only: KINT, KREAL
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: primitive_from_conserved_2d, analytic_maxwellian_2d
    public :: analytic_maxwellian_direction_2d


contains

    ! Public keyword temp retained: it means T_native=1/lambda, not theta.
    subroutine primitive_from_conserved_2d(w,ck,rho,u,v,temp,ok)
        real(KREAL), intent(in) :: w(4)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: rho,u,v,temp
        logical, intent(out) :: ok
        real(KREAL) :: kinetic,einternal

        rho=0.0_KREAL; u=0.0_KREAL; v=0.0_KREAL; temp=0.0_KREAL; ok=.false.
        if (ck<0 .or. .not.all(ieee_is_finite(w))) return
        rho=w(1)
        if (.not.ieee_is_finite(rho) .or. rho<=0.0_KREAL) return
        u=w(2)/rho
        v=w(3)/rho
        kinetic=0.5_KREAL*(u*u+v*v)
        einternal=w(4)/rho-kinetic
        temp=4.0_KREAL*einternal/real(ck+2,KREAL)
        ok=ieee_is_finite(u) .and. ieee_is_finite(v) .and. ieee_is_finite(temp) .and. temp>0.0_KREAL
    end subroutine primitive_from_conserved_2d

    subroutine analytic_maxwellian_2d(w,ck,uvel,vvel,mh,mb,ok)
        real(KREAL), intent(in) :: w(4),uvel(:),vvel(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: mh(:),mb(:)
        logical, intent(out) :: ok
        real(KREAL) :: rho,u,v,temp,lambda,c2,pi_value
        integer :: k
        logical :: primitive_ok

        mh=0.0_KREAL; mb=0.0_KREAL; ok=.false.
        if (size(mh)/=size(uvel) .or. size(mb)/=size(uvel) .or. size(vvel)/=size(uvel)) return
        call primitive_from_conserved_2d(w,ck,rho,u,v,temp,primitive_ok)
        if (.not.primitive_ok) return
        lambda=1.0_KREAL/temp
        pi_value=acos(-1.0_KREAL)
        do k=1,size(uvel)
            c2=(uvel(k)-u)**2+(vvel(k)-v)**2
            mh(k)=rho*lambda/pi_value*exp(-lambda*c2)
            mb(k)=0.5_KREAL*real(ck,KREAL)*temp*mh(k)
        end do
        ok=all(ieee_is_finite(mh)) .and. all(ieee_is_finite(mb))
    end subroutine analytic_maxwellian_2d

    subroutine analytic_maxwellian_direction_2d(w,w_direction,ck,uvel,vvel,mh,mb,mh_direction, &
        mb_direction,ok)
        real(KREAL), intent(in) :: w(4),w_direction(4),uvel(:),vvel(:)
        integer(KINT), intent(in) :: ck
        real(KREAL), intent(out) :: mh(:),mb(:),mh_direction(:),mb_direction(:)
        logical, intent(out) :: ok
        real(KREAL) :: rho,u,v,temp,lambda,du,dv,drho,dtemp,c2,log_direction
        real(KREAL) :: gamma_factor,pi_value
        integer :: k
        logical :: primitive_ok

        mh=0.0_KREAL; mb=0.0_KREAL; mh_direction=0.0_KREAL; mb_direction=0.0_KREAL; ok=.false.
        if (size(mh)/=size(uvel) .or. size(mb)/=size(uvel) .or. size(mh_direction)/=size(uvel) .or. &
            size(mb_direction)/=size(uvel) .or. size(vvel)/=size(uvel)) return
        call primitive_from_conserved_2d(w,ck,rho,u,v,temp,primitive_ok)
        if (.not.primitive_ok .or. .not.all(ieee_is_finite(w_direction))) return
        drho=w_direction(1)
        du=(w_direction(2)-u*drho)/rho
        dv=(w_direction(3)-v*drho)/rho
        gamma_factor=4.0_KREAL/real(ck+2,KREAL)
        dtemp=gamma_factor*(w_direction(4)/rho-(w(4)/rho**2)*drho-u*du-v*dv)
        lambda=1.0_KREAL/temp
        pi_value=acos(-1.0_KREAL)
        do k=1,size(uvel)
            c2=(uvel(k)-u)**2+(vvel(k)-v)**2
            mh(k)=rho*lambda/pi_value*exp(-lambda*c2)
            mb(k)=0.5_KREAL*real(ck,KREAL)*temp*mh(k)
            log_direction=drho/rho + 2.0_KREAL*((uvel(k)-u)*du+(vvel(k)-v)*dv)/temp + &
                (c2/temp-1.0_KREAL)*dtemp/temp
            mh_direction(k)=mh(k)*log_direction
            mb_direction(k)=mb(k)*(dtemp/temp+log_direction)
        end do
        ok=all(ieee_is_finite(mh)) .and. all(ieee_is_finite(mb)) .and. &
            all(ieee_is_finite(mh_direction)) .and. all(ieee_is_finite(mb_direction))
    end subroutine analytic_maxwellian_direction_2d





end module mod_CorrectedWF2DEquilibrium



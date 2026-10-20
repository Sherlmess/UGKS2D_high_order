module mod_CorrectedWF2DEquilibrium
    use Precision, only: KINT, KREAL
    use mod_CorrectedWF2DTypes, only: CollisionPointModel2D
    use WENO5AO2D, only: GQ3_NODE, GQ3_WEIGHT, reconstruct_weno5ao_tensor_point, &
        reconstruct_weno5ao_tensor_point_direction
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: primitive_from_conserved_2d, analytic_maxwellian_2d
    public :: analytic_maxwellian_direction_2d
    public :: constant_bgk_cell_average_target_2d
    public :: constant_bgk_cell_average_target_direction_2d

    type, extends(CollisionPointModel2D), public :: ConstantBGKModel2D
        real(KREAL) :: tau0=-1.0_KREAL
    contains
        procedure :: point_target_value => constant_bgk_point_target_value
        procedure :: point_target_direction => constant_bgk_point_target_direction
    end type ConstantBGKModel2D

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

    subroutine constant_bgk_point_target_value(self,w_point,h_point,b_point,uvel,vvel,weight,g_h,g_b,tau,ok)
        class(ConstantBGKModel2D), intent(in) :: self
        real(KREAL), intent(in) :: w_point(4),h_point(:),b_point(:),uvel(:),vvel(:),weight(:)
        real(KREAL), intent(out) :: g_h(:),g_b(:),tau
        logical, intent(out) :: ok

        g_h=0.0_KREAL; g_b=0.0_KREAL; tau=0.0_KREAL; ok=.false.
        if (size(h_point)/=size(uvel) .or. size(b_point)/=size(uvel) .or. size(vvel)/=size(uvel) .or. &
            size(weight)/=size(uvel)) return
        if (.not.all(ieee_is_finite(h_point)) .or. .not.all(ieee_is_finite(b_point)) .or. &
            .not.all(ieee_is_finite(uvel)) .or. .not.all(ieee_is_finite(vvel)) .or. &
            .not.all(ieee_is_finite(weight)) .or. any(weight<=0.0_KREAL)) return
        if (.not.ieee_is_finite(self%tau0) .or. self%tau0<=0.0_KREAL) return
        call analytic_maxwellian_2d(w_point,self%thermo%ck,uvel,vvel,g_h,g_b,ok)
        tau=self%tau0
        ok=ok .and. ieee_is_finite(tau)
    end subroutine constant_bgk_point_target_value

    subroutine constant_bgk_point_target_direction(self,w_point,h_point,b_point,w_direction, &
        h_direction,b_direction,uvel,vvel,weight,g_h,g_b,tau,g_h_direction,g_b_direction,tau_direction,ok)
        class(ConstantBGKModel2D), intent(in) :: self
        real(KREAL), intent(in) :: w_point(4),h_point(:),b_point(:),w_direction(4),h_direction(:),b_direction(:)
        real(KREAL), intent(in) :: uvel(:),vvel(:),weight(:)
        real(KREAL), intent(out) :: g_h(:),g_b(:),tau,g_h_direction(:),g_b_direction(:),tau_direction
        logical, intent(out) :: ok

        g_h=0.0_KREAL; g_b=0.0_KREAL; tau=0.0_KREAL
        g_h_direction=0.0_KREAL; g_b_direction=0.0_KREAL; tau_direction=0.0_KREAL; ok=.false.
        if (size(h_point)/=size(uvel) .or. size(b_point)/=size(uvel) .or. size(vvel)/=size(uvel) .or. &
            size(weight)/=size(uvel) .or. size(h_direction)/=size(uvel) .or. size(b_direction)/=size(uvel)) return
        if (.not.all(ieee_is_finite(h_point)) .or. .not.all(ieee_is_finite(b_point)) .or. &
            .not.all(ieee_is_finite(h_direction)) .or. .not.all(ieee_is_finite(b_direction)) .or. &
            .not.all(ieee_is_finite(uvel)) .or. .not.all(ieee_is_finite(vvel)) .or. &
            .not.all(ieee_is_finite(weight)) .or. any(weight<=0.0_KREAL)) return
        if (.not.ieee_is_finite(self%tau0) .or. self%tau0<=0.0_KREAL) return
        call analytic_maxwellian_direction_2d(w_point,w_direction,self%thermo%ck,uvel,vvel,g_h,g_b, &
            g_h_direction,g_b_direction,ok)
        tau=self%tau0
        tau_direction=0.0_KREAL
        ok=ok .and. ieee_is_finite(tau)
    end subroutine constant_bgk_point_target_direction

    subroutine constant_bgk_cell_average_target_2d(w,model,uvel,vvel,weight,gbar_h,gbar_b,ok)
        real(KREAL), intent(in) :: w(:,:,:),uvel(:),vvel(:),weight(:)
        type(ConstantBGKModel2D), intent(in) :: model
        real(KREAL), intent(inout) :: gbar_h(:,:,:),gbar_b(:,:,:)
        logical, intent(out) :: ok
        integer :: nx,ny,nvel,i,j,a,b,m
        real(KREAL) :: what(4),what_avg(4),wc(4),mh(size(uvel)),mb(size(uvel))
        real(KREAL), allocatable :: trial_h(:,:,:),trial_b(:,:,:)
        logical :: local_ok

        ok=.false.
        nx=size(w,2); ny=size(w,3); nvel=size(uvel)
        if (size(w,1)/=4 .or. nx<5 .or. ny<5 .or. size(vvel)/=nvel .or. size(weight)/=nvel .or. &
            size(gbar_h,1)/=nvel .or. size(gbar_h,2)/=nx .or. size(gbar_h,3)/=ny .or. &
            size(gbar_b,1)/=nvel .or. size(gbar_b,2)/=nx .or. size(gbar_b,3)/=ny) return
        if (.not.all(ieee_is_finite(w)) .or. .not.all(ieee_is_finite(uvel)) .or. &
            .not.all(ieee_is_finite(vvel)) .or. .not.all(ieee_is_finite(weight)) .or. any(weight<=0.0_KREAL)) return
        if (model%thermo%ck<0 .or. model%tau0<=0.0_KREAL .or. .not.ieee_is_finite(model%tau0)) return
        allocate(trial_h(nvel,nx,ny),trial_b(nvel,nx,ny))
        trial_h=0.0_KREAL; trial_b=0.0_KREAL

        do j=1,ny
            do i=1,nx
                what_avg=0.0_KREAL
                do b=1,3
                    do a=1,3
                        do m=1,4
                            call reconstruct_weno5ao_tensor_point(w(m,:,:),i,j,0.5_KREAL*GQ3_NODE(a), &
                                0.5_KREAL*GQ3_NODE(b),what(m),local_ok)
                            if (.not.local_ok) return
                            what_avg(m)=what_avg(m)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*what(m)
                        end do
                    end do
                end do
                do b=1,3
                    do a=1,3
                        wc=0.0_KREAL
                        do m=1,4
                            call reconstruct_weno5ao_tensor_point(w(m,:,:),i,j,0.5_KREAL*GQ3_NODE(a), &
                                0.5_KREAL*GQ3_NODE(b),what(m),local_ok)
                            if (.not.local_ok) return
                            wc(m)=what(m)+w(m,i,j)-what_avg(m)
                        end do
                        call analytic_maxwellian_2d(wc,model%thermo%ck,uvel,vvel,mh,mb,local_ok)
                        if (.not.local_ok) return
                        trial_h(:,i,j)=trial_h(:,i,j)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*mh
                        trial_b(:,i,j)=trial_b(:,i,j)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*mb
                    end do
                end do
            end do
        end do
        ok=all(ieee_is_finite(trial_h)) .and. all(ieee_is_finite(trial_b))
        if (ok) then
            gbar_h=trial_h; gbar_b=trial_b
        end if
    end subroutine constant_bgk_cell_average_target_2d

    subroutine constant_bgk_cell_average_target_direction_2d(w,w_direction,model,uvel,vvel,weight, &
        gbar_h,gbar_b,gbar_h_direction,gbar_b_direction,ok)
        real(KREAL), intent(in) :: w(:,:,:),w_direction(:,:,:),uvel(:),vvel(:),weight(:)
        type(ConstantBGKModel2D), intent(in) :: model
        real(KREAL), intent(inout) :: gbar_h(:,:,:),gbar_b(:,:,:),gbar_h_direction(:,:,:),gbar_b_direction(:,:,:)
        logical, intent(out) :: ok
        integer :: nx,ny,nvel,i,j,a,b,m
        real(KREAL) :: what(4),what_direction(4),what_avg(4),what_avg_direction(4),wc(4),wc_direction(4)
        real(KREAL) :: mh(size(uvel)),mb(size(uvel)),mhd(size(uvel)),mbd(size(uvel))
        real(KREAL), allocatable :: trial_h(:,:,:),trial_b(:,:,:),trial_hd(:,:,:),trial_bd(:,:,:)
        logical :: local_ok

        ok=.false.
        nx=size(w,2); ny=size(w,3); nvel=size(uvel)
        if (size(w,1)/=4 .or. size(w_direction,1)/=4 .or. size(w_direction,2)/=nx .or. &
            size(w_direction,3)/=ny .or. nx<5 .or. ny<5 .or. size(vvel)/=nvel .or. size(weight)/=nvel .or. &
            size(gbar_h,1)/=nvel .or. size(gbar_h,2)/=nx .or. size(gbar_h,3)/=ny .or. &
            size(gbar_b,1)/=nvel .or. size(gbar_b,2)/=nx .or. size(gbar_b,3)/=ny .or. &
            size(gbar_h_direction,1)/=nvel .or. size(gbar_h_direction,2)/=nx .or. size(gbar_h_direction,3)/=ny .or. &
            size(gbar_b_direction,1)/=nvel .or. size(gbar_b_direction,2)/=nx .or. size(gbar_b_direction,3)/=ny) return
        if (.not.all(ieee_is_finite(w)) .or. .not.all(ieee_is_finite(w_direction)) .or. &
            .not.all(ieee_is_finite(uvel)) .or. .not.all(ieee_is_finite(vvel)) .or. &
            .not.all(ieee_is_finite(weight)) .or. any(weight<=0.0_KREAL)) return
        if (model%thermo%ck<0 .or. model%tau0<=0.0_KREAL .or. .not.ieee_is_finite(model%tau0)) return
        allocate(trial_h(nvel,nx,ny),trial_b(nvel,nx,ny),trial_hd(nvel,nx,ny),trial_bd(nvel,nx,ny))
        trial_h=0.0_KREAL; trial_b=0.0_KREAL; trial_hd=0.0_KREAL; trial_bd=0.0_KREAL

        do j=1,ny
            do i=1,nx
                what_avg=0.0_KREAL; what_avg_direction=0.0_KREAL
                do b=1,3
                    do a=1,3
                        do m=1,4
                            call reconstruct_weno5ao_tensor_point_direction(w(m,:,:),w_direction(m,:,:),i,j, &
                                0.5_KREAL*GQ3_NODE(a),0.5_KREAL*GQ3_NODE(b),what(m),what_direction(m),local_ok)
                            if (.not.local_ok) return
                            what_avg(m)=what_avg(m)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*what(m)
                            what_avg_direction(m)=what_avg_direction(m)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*what_direction(m)
                        end do
                    end do
                end do
                do b=1,3
                    do a=1,3
                        wc=0.0_KREAL; wc_direction=0.0_KREAL
                        do m=1,4
                            call reconstruct_weno5ao_tensor_point_direction(w(m,:,:),w_direction(m,:,:),i,j, &
                                0.5_KREAL*GQ3_NODE(a),0.5_KREAL*GQ3_NODE(b),what(m),what_direction(m),local_ok)
                            if (.not.local_ok) return
                            wc(m)=what(m)+w(m,i,j)-what_avg(m)
                            wc_direction(m)=what_direction(m)+w_direction(m,i,j)-what_avg_direction(m)
                        end do
                        call analytic_maxwellian_direction_2d(wc,wc_direction,model%thermo%ck,uvel,vvel, &
                            mh,mb,mhd,mbd,local_ok)
                        if (.not.local_ok) return
                        trial_h(:,i,j)=trial_h(:,i,j)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*mh
                        trial_b(:,i,j)=trial_b(:,i,j)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*mb
                        trial_hd(:,i,j)=trial_hd(:,i,j)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*mhd
                        trial_bd(:,i,j)=trial_bd(:,i,j)+0.25_KREAL*GQ3_WEIGHT(a)*GQ3_WEIGHT(b)*mbd
                    end do
                end do
            end do
        end do
        ok=all(ieee_is_finite(trial_h)) .and. all(ieee_is_finite(trial_b)) .and. &
            all(ieee_is_finite(trial_hd)) .and. all(ieee_is_finite(trial_bd))
        if (ok) then
            gbar_h=trial_h; gbar_b=trial_b; gbar_h_direction=trial_hd; gbar_b_direction=trial_bd
        end if
    end subroutine constant_bgk_cell_average_target_direction_2d

end module mod_CorrectedWF2DEquilibrium



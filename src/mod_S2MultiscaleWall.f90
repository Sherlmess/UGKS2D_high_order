! Flat Cartesian BGK migration from implicitUGKS 753607a, MultiscaleFluxBoundary.
! See MAPPING_AND_FORMULAS.md for explicit missing factors, orientation, and
! analytic-vs-discrete macro/reflection changes. All fluxes include dt, not length.
module mod_S2MultiscaleWall
    use mod_S2MultiscaleWallLegacy, only: legacy_wall_node=>s2w2_wall_node
    use Precision, only: KINT,KREAL
    use flow_para, only: ck,Pr,S2PhysicalGas
    use mod_S2PhysicalGas, only: s2_tau
    use mod_CorrectedWF2DClosure, only: discrete_equilibrium_q_2d
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: s2w2_wall_node,s2w2_time_coefficients,s2w2_instantaneous_node
    integer,parameter::MNUM=6,MTUM=5
    real(KREAL),parameter::PI=4.0_KREAL*atan(1.0_KREAL)
    real(KREAL),parameter::gate=2.0e-12_KREAL
contains
    subroutine s2w2_time_coefficients(dt,tau,mt)
        real(KREAL),intent(in)::dt,tau
        real(KREAL),intent(out)::mt(5)
        real(KREAL)::x,t,ex
        integer::n
        x=dt/tau
        if(x<0.1_KREAL)then
            mt=0.0_KREAL;t=1.0_KREAL
            do n=1,24
                t=t*(-x)/real(n,KREAL)
                mt(4)=mt(4)-tau*t
                if(n>=2)then
                    mt(1)=mt(1)+tau*t
                    mt(5)=mt(5)+tau*tau*real(n-1,KREAL)*t
                endif
                if(n>=3)then
                    mt(2)=mt(2)+tau*tau*real(n-2,KREAL)*t
                    mt(3)=mt(3)-tau*tau*t
                endif
            enddo
        else
            ex=exp(-x);mt(4)=tau*(1.0_KREAL-ex)
            mt(5)=-tau*dt*ex+tau*mt(4);mt(1)=dt-mt(4)
            mt(2)=-tau*mt(1)+mt(5);mt(3)=0.5_KREAL*dt**2-tau*mt(1)
        endif
    end subroutine

    subroutine s2w2_wall_node(side,dt,dn,wcell,hin,bin,hn,bn,ht,bt,u,v,weight,ux,theta, &
        fh,fb,fw,ok,incidence_out,rho_out,mt_out,f0_h,f0_b,k_h,k_b)
        integer,intent(in)::side
        real(KREAL),intent(in)::dt,dn,wcell(4),hin(:),bin(:),hn(:),bn(:),ht(:),bt(:)
        real(KREAL),intent(in)::u(:),v(:),weight(:),ux,theta
        real(KREAL),intent(out)::fh(:),fb(:),fw(4)
        logical,intent(out)::ok
        real(KREAL),intent(out),optional::incidence_out,rho_out,mt_out(5),f0_h(:),f0_b(:),k_h(:),k_b(:)
        real(KREAL)::a0,a1,den,r0,r1,qh(size(u)),qb(size(u))
        real(KREAL)::tau_wall,heat_local(2),corrh(size(u)),corrb(size(u)),coeff,qdot,shapeh
        real(KREAL)::prim_g(4),wg(4),wc(4),sw(4),sw_t(4),an(4),at(4),aTtime(4)
        real(KREAL)::Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2),MuL(0:MNUM),MuR(0:MNUM),mt(5)
        real(KREAL)::mh(size(u)),mb(size(u)),wall(size(u)),wallb(size(u)),ginh(size(u)),ginb(size(u))
        real(KREAL)::vn,vt,v2,lambda,p,incidence,reflection,rhow,Han,Hatan,Hatime,Ban,Batan,Batime
        integer::k
        if(S2PhysicalGas==0)then
            call legacy_wall_node(side,dt,dn,wcell,hin,bin,hn,bn,ht,bt,u,v,weight,ux,theta, &
                fh,fb,fw,ok,incidence_out,rho_out,mt_out,f0_h,f0_b,k_h,k_b)
            return
        endif
        ok=.false.;fh=0.0_KREAL;fb=0.0_KREAL;fw=0.0_KREAL
        if(abs(side)/=1.or.dt<=0.0_KREAL.or.dn<=0.0_KREAL.or.theta<=0.0_KREAL)return
        if(.not.all(ieee_is_finite(wcell)).or.wcell(1)<=0.0_KREAL)return
        p=(2.0_KREAL/real(ck+2,KREAL))*(wcell(4)-0.5_KREAL* &
            (wcell(2)**2+wcell(3)**2)/wcell(1))
        if(p<=0.0_KREAL.or..not.ieee_is_finite(p))return
        lambda=0.5_KREAL/theta
        prim_g=[p/theta,0.0_KREAL,-ux,lambda]
        ! Preserve original pressure-matched wall equilibrium and center-wall
        ! normal closure; corrected derivative is along the same +y as vn.
        wc=[wcell(1),wcell(3),-wcell(2),wcell(4)]
        wg=[prim_g(1),0.0_KREAL,-prim_g(1)*ux, &
            prim_g(1)*(real(ck+2,KREAL)/(4.0_KREAL*lambda)+0.5_KREAL*ux*ux)]
        sw=(wg-wc)/(real(side,KREAL)*0.5_KREAL*dn)
        an=MicroSlope(prim_g,sw)
        sw_t=0.0_KREAL
        do k=1,size(u)
            vn=v(k);vt=-u(k);v2=vn*vn+vt*vt
            ! ht/bt arrive as Cartesian +x derivatives. Local tangent is -x.
            sw_t=sw_t+weight(k)*[-ht(k),-vn*ht(k),-vt*ht(k), &
                -0.5_KREAL*(v2*ht(k)+bt(k))]
        enddo
        ! Use the same discrete conservative Q as the S2 volume target.
        ! Both reduced distributions must share its fitted thermal width.
        call wall_equilibrium(prim_g(1),ux,theta,u,v,weight,mh,mb,ok)
        if(.not.ok)return
        wall=mh/prim_g(1);wallb=mb/prim_g(1)
        at=MicroSlope(prim_g,sw_t)
        call CalcMoment(prim_g,Mu,Mv,Mxi,MuL,MuR)
        sw=-prim_g(1)*(Moment_auvxi(an,Mu,Mv,Mxi,1,0)+Moment_auvxi(at,Mu,Mv,Mxi,0,1))
        aTtime=MicroSlope(prim_g,sw)
        tau_wall=s2_tau([wg(1),-wg(3),wg(2),wg(4)],3)
        if(S2PhysicalGas==1)then
            call s2w2_time_coefficients(dt,tau_wall,mt)
        else
            call s2w2_time_coefficients(dt,0.1_KREAL,mt)
        endif
        corrh=0.0_KREAL;corrb=0.0_KREAL
        if(S2PhysicalGas==1.and.Pr/=1.0_KREAL)then
            ! Same local velocity frame and old reconstructed incident heat as native type1.
            heat_local=0.0_KREAL
            do k=1,size(u)
                vn=v(k);vt=-u(k)+ux;v2=vn*vn+vt*vt
                heat_local=heat_local+0.5_KREAL*weight(k)*[vn,vt]*(v2*hin(k)+bin(k))
            enddo
            coeff=0.8_KREAL*(1.0_KREAL-Pr)*lambda*lambda/prim_g(1)
            do k=1,size(u)
                vn=v(k);vt=-u(k)+ux;v2=vn*vn+vt*vt
                qdot=vn*heat_local(1)+vt*heat_local(2);shapeh=2.0_KREAL*lambda*v2-4.0_KREAL
                corrh(k)=coeff*qdot*shapeh*mh(k)
                corrb(k)=coeff*qdot*(shapeh+2.0_KREAL/real(ck,KREAL))*mb(k)
            enddo
        endif
        incidence=0.0_KREAL;reflection=0.0_KREAL
        do k=1,size(u)
            vn=v(k);vt=-u(k);v2=vn*vn+vt*vt
            Han=(an(1)+an(2)*vn+an(3)*vt)*mh(k)+0.5_KREAL*an(4)*(v2*mh(k)+mb(k))
            Hatan=(at(1)+at(2)*vn+at(3)*vt)*mh(k)+0.5_KREAL*at(4)*(v2*mh(k)+mb(k))
            Hatime=(aTtime(1)+aTtime(2)*vn+aTtime(3)*vt)*mh(k) &
                +0.5_KREAL*aTtime(4)*(v2*mh(k)+mb(k))
            Ban=(an(1)+an(2)*vn+an(3)*vt)*mb(k)+0.5_KREAL*an(4)*(v2*mb(k)+Mxi(2)*mh(k))
            Batan=(at(1)+at(2)*vn+at(3)*vt)*mb(k)+0.5_KREAL*at(4)*(v2*mb(k)+Mxi(2)*mh(k))
            Batime=(aTtime(1)+aTtime(2)*vn+aTtime(3)*vt)*mb(k) &
                +0.5_KREAL*aTtime(4)*(v2*mb(k)+Mxi(2)*mh(k))
            ginh(k)=mt(1)*mh(k)+mt(2)*(vn*Han+vt*Hatan)+mt(3)*Hatime &
                +mt(4)*hin(k)-mt(5)*(vn*hn(k)+u(k)*ht(k))
            ginb(k)=mt(1)*mb(k)+mt(2)*(vn*Ban+vt*Batan)+mt(3)*Batime &
                +mt(4)*bin(k)-mt(5)*(vn*bn(k)+u(k)*bt(k))
            if(S2PhysicalGas==1)then
                ginh(k)=ginh(k)+mt(1)*corrh(k)
                ginb(k)=ginb(k)+mt(1)*corrb(k)
            endif
            if(real(side,KREAL)*vn>0.0_KREAL)then
                incidence=incidence+weight(k)*vn*ginh(k)
            else
                reflection=reflection+dt*weight(k)*vn*wall(k)
            endif
        enddo
        if(.not.ieee_is_finite(reflection).or.abs(reflection)<=tiny(reflection))return
        rhow=-incidence/reflection
        if(.not.ieee_is_finite(rhow).or.rhow<=0.0_KREAL)return
        do k=1,size(u)
            if(real(side,KREAL)*v(k)>0.0_KREAL)then
                fh(k)=v(k)*ginh(k);fb(k)=v(k)*ginb(k)
            else
                fh(k)=dt*v(k)*rhow*wall(k);fb(k)=dt*v(k)*rhow*wallb(k)
            endif
            fw=fw+weight(k)*[fh(k),u(k)*fh(k),v(k)*fh(k), &
                0.5_KREAL*((u(k)**2+v(k)**2)*fh(k)+fb(k))]
        enddo
        if(present(f0_h))then
            if(.not.(present(f0_b).and.present(k_h).and.present(k_b)))return
            qh=(mh+corrh-hin)/tau_wall-v*hn-u*ht;qb=(mb+corrb-bin)/tau_wall-v*bn-u*bt
            a0=0;a1=0;den=0
            do k=1,size(u)
                if(side*v(k)>0.0_KREAL)then
                    a0=a0+weight(k)*v(k)*hin(k);a1=a1+weight(k)*v(k)*qh(k)
                else
                    den=den+weight(k)*v(k)*wall(k)
                endif
            enddo
            r0=-a0/den;r1=-a1/den
            if(.not.all(ieee_is_finite([r0,r1])).or.r0<=0.0_KREAL)return
            do k=1,size(u)
                if(side*v(k)>0.0_KREAL)then
                    f0_h(k)=v(k)*hin(k);f0_b(k)=v(k)*bin(k)
                    k_h(k)=v(k)*qh(k);k_b(k)=v(k)*qb(k)
                else
                    f0_h(k)=v(k)*r0*wall(k);f0_b(k)=v(k)*r0*wallb(k)
                    k_h(k)=v(k)*r1*wall(k);k_b(k)=v(k)*r1*wallb(k)
                endif
            enddo
            if(abs(sum(weight*f0_h))>gate.or.abs(sum(weight*k_h))>gate)return
        endif
        if(present(incidence_out))incidence_out=incidence
        if(present(rho_out))rho_out=rhow
        if(present(mt_out))mt_out=mt
        ok=all(ieee_is_finite(fh)).and.all(ieee_is_finite(fb)).and.all(ieee_is_finite(fw)) &
            .and.abs(fw(1))<=gate*max(1.0_KREAL,abs(incidence))
    end subroutine s2w2_wall_node

    ! Exact dt->0 diffuse limit, not an integrated flux divided by dt.
    subroutine s2w2_instantaneous_node(side,hin,bin,u,v,weight,ux,theta,fh,fb,fw,ok)
        integer,intent(in)::side
        real(KREAL),intent(in)::hin(:),bin(:),u(:),v(:),weight(:),ux,theta
        real(KREAL),intent(out)::fh(:),fb(:),fw(4)
        logical,intent(out)::ok
        real(KREAL)::wall(size(u)),wallb(size(u)),lambda,inc,den,rho
        integer::k
        ok=.false.;lambda=0.5_KREAL/theta;inc=0;den=0;fw=0
        if(S2PhysicalGas==0)then
            wall=lambda/PI*exp(-lambda*(v*v+(u-ux)**2))
            wallb=real(ck,KREAL)*theta*wall
        else
            call wall_equilibrium(1.0_KREAL,ux,theta,u,v,weight,wall,wallb,ok)
            if(.not.ok)return
        endif
        do k=1,size(u)
            if(side*v(k)>0.0_KREAL)then
                inc=inc+weight(k)*v(k)*hin(k)
            else
                den=den+weight(k)*v(k)*wall(k)
            endif
        enddo
        rho=-inc/den
        if(.not.ieee_is_finite(rho).or.rho<=0.0_KREAL)return
        do k=1,size(u)
            if(side*v(k)>0.0_KREAL)then
                fh(k)=v(k)*hin(k);fb(k)=v(k)*bin(k)
            else
                fh(k)=v(k)*rho*wall(k)
                if(S2PhysicalGas==0)then
                    ! Preserve legacy floating-point evaluation order as well.
                    fb(k)=real(ck,KREAL)*theta*fh(k)
                else
                    fb(k)=v(k)*rho*wallb(k)
                endif
            endif
            fw=fw+weight(k)*[fh(k),u(k)*fh(k),v(k)*fh(k), &
                0.5_KREAL*((u(k)**2+v(k)**2)*fh(k)+fb(k))]
        enddo
        ok=all(ieee_is_finite(fh)).and.all(ieee_is_finite(fb)).and.abs(fw(1))<=gate*max(1.0_KREAL,abs(inc))
    end subroutine s2w2_instantaneous_node

    ! Diffuse wall shape with the existing S2 four-moment closure. No flux
    ! projection: rhow still closes the actual incident half-space mass flux.
    subroutine wall_equilibrium(rho,ux,theta,u,v,weight,mh,mb,ok)
        real(KREAL),intent(in)::rho,ux,theta,u(:),v(:),weight(:)
        real(KREAL),intent(out)::mh(:),mb(:)
        logical,intent(out)::ok
        real(KREAL)::w(4),fit(4),jac(4,4)
        w=[rho,rho*ux,0.0_KREAL,rho*(real(ck+2,KREAL)*theta/2.0_KREAL+ux*ux/2.0_KREAL)]
        call discrete_equilibrium_q_2d(w,ck,u,v,weight,mh,mb,fit,jac,ok)
    end subroutine wall_equilibrium

    function MicroSlope(prim,sw)
        real(KREAL), intent(in)                         :: prim(4),sw(4)
        real(KREAL)                                     :: MicroSlope(4)

        MicroSlope(4) = 4.0*prim(4)**2/((CK+2)*prim(1))*(2.0*sw(4)-2.0*prim(2)*sw(2)-2.0*prim(3)*sw(3) &
            +sw(1)*(prim(2)**2+prim(3)**2-0.5*(CK+2)/prim(4)))
        MicroSlope(3) = 2.0*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*MicroSlope(4)
        MicroSlope(2) = 2.0*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*MicroSlope(4)
        MicroSlope(1) = sw(1)/prim(1)-prim(2)*MicroSlope(2)-prim(3)*MicroSlope(3)-0.5*(prim(2)**2+prim(3)**2 &
            +0.5*(CK+2)/prim(4))*MicroSlope(4)
    end function MicroSlope
    
    
    !--------------------------------------------------
    !>calculate moments of velocity and \xi
    !>@param[in] prim       :primary variables
    !>@param[out] Mu,Mv     :<u^n>,<v^m>
    !>@param[out] Mxi       :<\xi^2n>
    !>@param[out] MuL,MuR   :<u^n>_{>0},<u^n>_{<0}
    !--------------------------------------------------
    subroutine CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL), intent(out)                        :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM)
        real(KREAL), intent(out)                        :: Mv(0:MTUM)
        real(KREAL), intent(out)                        :: Mxi(0:2)
        integer :: i

        !Moments of normal velocity
        MuL(0) = 0.5*erfc(-sqrt(prim(4))*prim(2))
        MuL(1) = prim(2)*MuL(0)+0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))
        MuR(0) = 0.5*erfc(sqrt(prim(4))*prim(2))
        MuR(1) = prim(2)*MuR(0)-0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))

        do i=2,MNUM
            MuL(i) = prim(2)*MuL(i-1)+0.5*(i-1)*MuL(i-2)/prim(4)
            MuR(i) = prim(2)*MuR(i-1)+0.5*(i-1)*MuR(i-2)/prim(4)
        end do

        Mu = MuL+MuR

        !Moments of tangential velocity
        Mv(0) = 1.0
        Mv(1) = prim(3)

        do i=2,MTUM
            Mv(i) = prim(3)*Mv(i-1)+0.5*(i-1)*Mv(i-2)/prim(4)
        end do

        !Moments of \xi
        Mxi(0) = 1.0 !<\xi^0>
        Mxi(1) = 0.5*CK/prim(4) !<\xi^2>
        Mxi(2) = CK*(CK+2.0)/(4.0*prim(4)**2) !<\xi^4>
    end subroutine CalcMoment

    
    !--------------------------------------------------
    !>Calculate <u^\alpha*v^\beta*\xi^\delta*\psi>
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^delta>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@param[in] delta      :exponential index of \xi
    !>@return    Moment_uvxi :moment of <u^\alpha*v^\beta*\xi^\delta*\psi>
    !--------------------------------------------------
    function Moment_uvxi(Mu,Mv,Mxi,alpha,beta,delta)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,beta,delta
        real(KREAL)                                     :: Moment_uvxi(4)

        Moment_uvxi(1) = Mu(alpha)*Mv(beta)*Mxi(delta/2)
        Moment_uvxi(2) = Mu(alpha+1)*Mv(beta)*Mxi(delta/2)
        Moment_uvxi(3) = Mu(alpha)*Mv(beta+1)*Mxi(delta/2)
        Moment_uvxi(4) = 0.5*(Mu(alpha+2)*Mv(beta)*Mxi(delta/2)+Mu(alpha)*Mv(beta+2)*Mxi(delta/2) &
            +Mu(alpha)*Mv(beta)*Mxi((delta+2)/2))
    end function Moment_uvxi

    
    !--------------------------------------------------
    !>Calculate <a*u^\alpha*v^\beta>
    !>@param[in] a          :micro slope of Maxwellian
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^l>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@return    Moment_auvxi  :moment of <a*u^\alpha*v^\beta>
    !--------------------------------------------------
    function Moment_auvxi(a,Mu,Mv,Mxi,alpha,beta)
        real(KREAL), intent(in)                         :: a(4)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,beta
        real(KREAL)                                     :: Moment_auvxi(4)

        Moment_auvxi = a(1)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+0,0)+&
                    a(2)*Moment_uvxi(Mu,Mv,Mxi,alpha+1,beta+0,0)+&
                    a(3)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+1,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+2,beta+0,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+2,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+0,2)
    end function Moment_auvxi
end module mod_S2MultiscaleWall

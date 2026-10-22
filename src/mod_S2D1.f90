module mod_S2D1
    use mod_S2MPI1Layout
    use mod_S2PhysicalGas, only: s2_tau
    use mod_VelocityTables, only: VelocityMode,FillVelocity
    use Precision, only: KINT,KREAL
    use flow_para, only: Nx_all,Ny_all,Nu,Nv,ck,Pr,BC_Type,S2W1_wall_mode, &
        S2W1_wall_ux,S2W1_wall_theta,FluxWall_type,S2PhysicalGas
    use control_para, only: dt_global,MAX_ITER,TimeIntegration_method,Nstep_show,S2D1_initial_mode,IsS2Steady,IsHighOrderImplicit
    use mpi_para, only: myID
    use const_para, only: PI,BCconst_periodic
    use mod_S2MultiscaleWall, only: s2w2_wall_node,s2w2_instantaneous_node
    use mod_S2D1Reconstruction, only: s2d1_face3,s2d1_side3,s2d1_cell9,s2d1_cell_faces9, &
        FirstLayerPlane,prepare_first_layer_plane,s2d1_cell_faces9_cached,s2d1_cell_faces9_wall, &
        VelocityTileWorkspace,VT_VELOCITIES,prepare_velocity_tile,reconstruct_velocity_tile_cell
    use mod_CorrectedWF2DClosure, only: moment_2d,discrete_equilibrium_q_2d,dq_2d,projector_2d, &
        ClosureFirstOrderContext,closure_first_prepare,closure_first_copy,closure_first_project, &
        closure_first_direction,closure_first_matches
    use mod_CorrectedWF2DEquilibrium, only: analytic_maxwellian_2d
    use mpi, only: MPI_COMM_WORLD,MPI_Comm_size,MPI_Abort
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: IsS2D1,IsS2CN2,IsS2Steady,RunS2D1
    public :: compute_flux,build_targets,advance_step,init_velocity,write_state,heat,make_mke
    public :: s2w1_wall_node
    public :: s2w1_compute_wall_flux
    public :: s2w1d_instantaneous_flux
    public :: s2w2r_correct
    public :: reconstruct_cell_pair ! research-only access to unchanged wall input constructor
    public :: w2_reconstruction_s,w2_wall_s,w2_wall_calls,w2_closure_s

    real(KREAL), parameter :: tau0=0.1_KREAL,pi_d=acos(-1.0_KREAL)
    real(KREAL), parameter :: rho0=1.0_KREAL,u0=0.6_KREAL,v0=0.3_KREAL,t0=1.0_KREAL
    real(KREAL), parameter :: amp=0.01_KREAL,gate=2.0e-12_KREAL
    real(KREAL), parameter :: micro_margin=64.0_KREAL*epsilon(1.0_KREAL)
    real(KREAL), parameter :: gw(3)=[5.0_KREAL/18.0_KREAL,4.0_KREAL/9.0_KREAL,5.0_KREAL/18.0_KREAL]

    real(KREAL),save::w2_reconstruction_s=0.0_KREAL,w2_wall_s=0.0_KREAL
    real(KREAL),save::w2_closure_s=0.0_KREAL
    integer,save::w2r_research_calls=0
    integer,save::w2_wall_calls=0

contains

    logical function IsS2D1()
        character(len=24) :: method
        method=trim(adjustl(TimeIntegration_method))
        IsS2D1=method=='S2D1'.or.method=='S2D1VL'
    end function IsS2D1

    logical function IsS2CN2()
        IsS2CN2=trim(adjustl(TimeIntegration_method))=='S2CN2'.or.IsHighOrderImplicit()
    end function IsS2CN2

    integer pure function ip(i,n)
        integer,intent(in)::i,n
        ip=modulo(i-1,n)+1
    end function ip

    subroutine abort_s2(message,code)
        character(len=*),intent(in)::message
        integer,intent(in)::code
        integer::ierr
        if(myID==0)write(*,'(A)')'S2D1 FAILURE: '//trim(message)
        call MPI_Abort(MPI_COMM_WORLD,code,ierr)
    end subroutine abort_s2

    ! theta_gas=p/rho=T_native/2; this private interface is not native prim.
    subroutine primitive(w,rho,ux,vy,theta_gas,p,ok)
        real(KREAL),intent(in)::w(4)
        real(KREAL),intent(out)::rho,ux,vy,theta_gas,p
        logical,intent(out)::ok
        real(KREAL)::eint
        ok=.false.;rho=w(1)
        if(.not.all(ieee_is_finite(w)).or.rho<=0.0_KREAL)return
        ux=w(2)/rho;vy=w(3)/rho
        eint=w(4)-0.5_KREAL*rho*(ux*ux+vy*vy)
        if(eint<=0.0_KREAL)return
        p=(2.0_KREAL/real(ck+2,KREAL))*eint;theta_gas=p/rho
        ok=all(ieee_is_finite([rho,ux,vy,theta_gas,p])).and.theta_gas>0.0_KREAL
    end subroutine primitive

    subroutine heat(w,h,b,u,v,weight,q,ok)
        real(KREAL),intent(in)::w(4),h(:),b(:),u(:),v(:),weight(:)
        real(KREAL),intent(out)::q(2)
        logical,intent(out)::ok
        real(KREAL)::rho,ux,vy,theta_gas,p,cx,cy,c2
        integer::k
        call primitive(w,rho,ux,vy,theta_gas,p,ok);q=0.0_KREAL
        if(.not.ok)return
        do k=1,size(u)
            cx=u(k)-ux;cy=v(k)-vy;c2=cx*cx+cy*cy
            q(1)=q(1)+0.5_KREAL*weight(k)*cx*(c2*h(k)+b(k))
            q(2)=q(2)+0.5_KREAL*weight(k)*cy*(c2*h(k)+b(k))
        end do
        ok=all(ieee_is_finite(q))
    end subroutine heat

    subroutine make_mke(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,ok,residual)
        real(KREAL),intent(in)::w(4),u(:),v(:),weight(:)
        real(KREAL),intent(out)::mh(:),mb(:),kxh(:),kxb(:),kyh(:),kyb(:),eah(:),eab(:)
        real(KREAL),intent(out)::kappa(2,2),residual
        logical,intent(out)::ok
        type(ClosureFirstOrderContext)::ctx
        call closure_first_prepare(ctx,w,ck,u,v,weight,ok)
        if(.not.ok)return
        call make_mke_prepared(ctx,w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,ok,residual)
    end subroutine make_mke

    subroutine make_mke_prepared(ctx,w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,ok,residual)
        type(ClosureFirstOrderContext),intent(in)::ctx
        real(KREAL),intent(in)::w(4),u(:),v(:),weight(:)
        real(KREAL),intent(out)::mh(:),mb(:),kxh(:),kxb(:),kyh(:),kyb(:),eah(:),eab(:)
        real(KREAL),intent(out)::kappa(2,2),residual
        logical,intent(out)::ok
        real(KREAL)::rho,ux,vy,theta_gas,p,cx,cy,c2,lambda,coeff,shapeh,shapeb
        real(KREAL)::rh(size(u)),rb(size(u)),mom(4),q(2),theta(4),jac(4,4),qht(size(u),4),qbt(size(u),4)
        logical::lok,matched
        integer::k
        ok=.false.;residual=0.0_KREAL;kappa=0.0_KREAL
        call closure_first_matches(ctx,w,ck,u,v,weight,matched);if(.not.matched)return
        call primitive(w,rho,ux,vy,theta_gas,p,lok);if(.not.lok)return
        call analytic_maxwellian_2d(w,ck,u,v,eah,eab,lok);if(.not.lok)return
        call closure_first_copy(ctx,mh,mb,theta,jac,qht,qbt,lok);if(.not.lok)return
        lambda=0.5_KREAL/theta_gas
        do k=1,size(u)
            cx=u(k)-ux;cy=v(k)-vy;c2=cx*cx+cy*cy
            coeff=0.8_KREAL*(1.0_KREAL-Pr)*lambda*lambda/rho
            shapeh=2.0_KREAL*lambda*c2-4.0_KREAL
            shapeb=shapeh+2.0_KREAL/real(ck,KREAL)
            rh(k)=coeff*cx*shapeh*mh(k);rb(k)=coeff*cx*shapeb*mb(k)
        end do
        call closure_first_project(ctx,rh,rb,kxh,kxb,lok);if(.not.lok)return
        do k=1,size(u)
            cx=u(k)-ux;cy=v(k)-vy;c2=cx*cx+cy*cy
            coeff=0.8_KREAL*(1.0_KREAL-Pr)*lambda*lambda/rho
            shapeh=2.0_KREAL*lambda*c2-4.0_KREAL
            shapeb=shapeh+2.0_KREAL/real(ck,KREAL)
            rh(k)=coeff*cy*shapeh*mh(k);rb(k)=coeff*cy*shapeb*mb(k)
        end do
        call closure_first_project(ctx,rh,rb,kyh,kyb,lok);if(.not.lok)return
        call heat(w,kxh,kxb,u,v,weight,q,lok);if(.not.lok)return;kappa(:,1)=q
        call heat(w,kyh,kyb,u,v,weight,q,lok);if(.not.lok)return;kappa(:,2)=q
        call moment_2d(mh,mb,u,v,weight,mom,lok);if(.not.lok)return
        residual=maxval(abs(mom-w))/max(1.0_KREAL,maxval(abs(w)))
        call moment_2d(kxh,kxb,u,v,weight,mom,lok);if(.not.lok)return
        residual=max(residual,maxval(abs(mom)))
        call moment_2d(kyh,kyb,u,v,weight,mom,lok);if(.not.lok)return
        residual=max(residual,maxval(abs(mom)))
        ok=residual<=gate.and.all(ieee_is_finite(kappa))
    end subroutine make_mke_prepared

    subroutine m_direction(w,dw,u,v,weight,dh,db,ok,residual)
        real(KREAL),intent(in)::w(4),dw(4),u(:),v(:),weight(:)
        real(KREAL),intent(out)::dh(:),db(:),residual
        logical,intent(out)::ok
        type(ClosureFirstOrderContext)::ctx
        call closure_first_prepare(ctx,w,ck,u,v,weight,ok)
        if(.not.ok)then;residual=huge(1.0_KREAL);return;endif
        call m_direction_prepared(ctx,w,dw,u,v,weight,dh,db,ok,residual)
    end subroutine m_direction

    subroutine m_direction_prepared(ctx,w,dw,u,v,weight,dh,db,ok,residual)
        type(ClosureFirstOrderContext),intent(in)::ctx
        real(KREAL),intent(in)::w(4),dw(4),u(:),v(:),weight(:)
        real(KREAL),intent(out)::dh(:),db(:),residual
        logical,intent(out)::ok
        real(KREAL)::mom(4)
        logical::lok,matched
        call closure_first_matches(ctx,w,ck,u,v,weight,matched)
        if(.not.matched)then;ok=.false.;residual=huge(1.0_KREAL);return;endif
        call closure_first_direction(ctx,dw,dh,db,lok)
        if(.not.lok)then;ok=.false.;residual=huge(1.0_KREAL);return;endif
        call moment_2d(dh,db,u,v,weight,mom,lok)
        if(.not.lok)then;ok=.false.;residual=huge(1.0_KREAL);return;endif
        residual=maxval(abs(mom-dw))/max(1.0_KREAL,maxval(abs(dw)))
        ok=residual<=gate.and.all(ieee_is_finite(dh)).and.all(ieee_is_finite(db))
    end subroutine m_direction_prepared

    subroutine target_at(w,q,u,v,weight,sh,sb,ok)
        real(KREAL),intent(in)::w(4),q(2),u(:),v(:),weight(:)
        real(KREAL),intent(out)::sh(:),sb(:)
        logical,intent(out)::ok
        real(KREAL)::mh(size(u)),mb(size(u)),kxh(size(u)),kxb(size(u))
        real(KREAL)::kyh(size(u)),kyb(size(u)),eah(size(u)),eab(size(u)),kap(2,2),res
        call make_mke(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kap,ok,res)
        if(ok)then
            sh=mh+kxh*q(1)+kyh*q(2);sb=mb+kxb*q(1)+kyb*q(2)
        endif
    end subroutine target_at

    ! v2.0 full-Van-Leer comparator: analytic Maxwellian and the original,
    ! unprojected Shakhov correction.  This deliberately does not use the
    ! S2 discrete-conservative M/K closure.
    subroutine old_target_at(w,q,u,v,sh,sb,mh,mb,ok)
        real(KREAL),intent(in)::w(4),q(2),u(:),v(:)
        real(KREAL),intent(out)::sh(:),sb(:),mh(:),mb(:)
        logical,intent(out)::ok
        real(KREAL)::rho,ux,vy,theta_gas,p,cx,cy,c2,lambda,coeff,qdot
        integer::k
        call primitive(w,rho,ux,vy,theta_gas,p,ok);if(.not.ok)return
        call analytic_maxwellian_2d(w,ck,u,v,mh,mb,ok);if(.not.ok)return
        lambda=0.5_KREAL/theta_gas
        coeff=0.8_KREAL*(1.0_KREAL-Pr)*lambda*lambda/rho
        do k=1,size(u)
            cx=u(k)-ux;cy=v(k)-vy;c2=cx*cx+cy*cy;qdot=cx*q(1)+cy*q(2)
            sh(k)=mh(k)+coeff*qdot*(2.0_KREAL*lambda*c2-4.0_KREAL)*mh(k)
            sb(k)=mb(k)+coeff*qdot*(2.0_KREAL*lambda*c2-2.0_KREAL)*mb(k)
        enddo
        ok=all(ieee_is_finite(sh)).and.all(ieee_is_finite(sb))
    end subroutine old_target_at

    pure function old_micro_slope(prim,sw) result(a)
        real(KREAL),intent(in)::prim(4),sw(4)
        real(KREAL)::a(4)
        a(4)=4.0_KREAL*prim(4)**2/(3.0_KREAL*prim(1))*(2.0_KREAL*sw(4)- &
            2.0_KREAL*prim(2)*sw(2)-2.0_KREAL*prim(3)*sw(3)+ &
            sw(1)*(prim(2)**2+prim(3)**2-1.5_KREAL/prim(4)))
        a(3)=2.0_KREAL*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*a(4)
        a(2)=2.0_KREAL*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*a(4)
        a(1)=sw(1)/prim(1)-prim(2)*a(2)-prim(3)*a(3)- &
            0.5_KREAL*(prim(2)**2+prim(3)**2+1.5_KREAL/prim(4))*a(4)
    end function old_micro_slope

    subroutine old_calc_moment(prim,mu,mv,mxi,mul,mur)
        real(KREAL),intent(in)::prim(4)
        real(KREAL),intent(out)::mu(0:6),mv(0:5),mxi(0:2),mul(0:6),mur(0:6)
        integer::n
        mu=0.0_KREAL;mv=0.0_KREAL;mxi=0.0_KREAL;mul=0.0_KREAL;mur=0.0_KREAL
        mu(0)=1.0_KREAL;mu(1)=prim(2)
        mv(0)=1.0_KREAL;mv(1)=prim(3)
        do n=2,6;mu(n)=prim(2)*mu(n-1)+0.5_KREAL*real(n-1,KREAL)/prim(4)*mu(n-2);enddo
        do n=2,5;mv(n)=prim(3)*mv(n-1)+0.5_KREAL*real(n-1,KREAL)/prim(4)*mv(n-2);enddo
        mul(0)=0.5_KREAL*erfc(-sqrt(prim(4))*prim(2))
        mul(1)=prim(2)*mul(0)+0.5_KREAL*exp(-prim(4)*prim(2)**2)/sqrt(pi_d*prim(4))
        do n=2,6;mul(n)=prim(2)*mul(n-1)+0.5_KREAL*real(n-1,KREAL)/prim(4)*mul(n-2);enddo
        mur=mu-mul
        mxi(0)=1.0_KREAL;mxi(1)=0.5_KREAL/prim(4);mxi(2)=0.75_KREAL/prim(4)**2
    end subroutine old_calc_moment

    pure function old_moment_uvxi(mu,mv,mxi,aa,bb,dd) result(ans)
        real(KREAL),intent(in)::mu(0:6),mv(0:5),mxi(0:2)
        integer,intent(in)::aa,bb,dd
        real(KREAL)::ans(4)
        ans(1)=mu(aa)*mv(bb)*mxi(dd/2)
        ans(2)=mu(aa+1)*mv(bb)*mxi(dd/2)
        ans(3)=mu(aa)*mv(bb+1)*mxi(dd/2)
        ans(4)=0.5_KREAL*(mu(aa+2)*mv(bb)*mxi(dd/2)+mu(aa)*mv(bb+2)*mxi(dd/2)+ &
            mu(aa)*mv(bb)*mxi((dd+2)/2))
    end function old_moment_uvxi

    pure function old_moment_auvxi(a,mu,mv,mxi,aa,bb) result(ans)
        real(KREAL),intent(in)::a(4),mu(0:6),mv(0:5),mxi(0:2)
        integer,intent(in)::aa,bb
        real(KREAL)::ans(4)
        ans=a(1)*old_moment_uvxi(mu,mv,mxi,aa,bb,0)+ &
            a(2)*old_moment_uvxi(mu,mv,mxi,aa+1,bb,0)+ &
            a(3)*old_moment_uvxi(mu,mv,mxi,aa,bb+1,0)+ &
            0.5_KREAL*a(4)*(old_moment_uvxi(mu,mv,mxi,aa+2,bb,0)+ &
            old_moment_uvxi(mu,mv,mxi,aa,bb+2,0)+old_moment_uvxi(mu,mv,mxi,aa,bb,2))
    end function old_moment_auvxi

    subroutine limit_pair(hmean,bmean,hv,bv,hn,ht,bn,bt,limited,signed)
        real(KREAL),intent(in)::hmean,bmean
        real(KREAL),intent(inout)::hv,bv,hn,ht,bn,bt
        integer,intent(inout)::limited,signed
        real(KREAL)::theta
        if(hmean<=0.0_KREAL.or.bmean<=0.0_KREAL)then
            hv=hmean;bv=bmean;hn=0.0_KREAL;ht=0.0_KREAL;bn=0.0_KREAL;bt=0.0_KREAL
            signed=signed+1;return
        endif
        theta=1.0_KREAL
        if(hv<0.0_KREAL)theta=min(theta,hmean/(hmean-hv))
        if(bv<0.0_KREAL)theta=min(theta,bmean/(bmean-bv))
        if(theta<1.0_KREAL)then
            theta=max(0.0_KREAL,theta)
            hv=hmean+theta*(hv-hmean);bv=bmean+theta*(bv-bmean)
            hn=theta*hn;ht=theta*ht;bn=theta*bn;bt=theta*bt
            limited=limited+1
        endif
    end subroutine limit_pair

    subroutine reconstruct_cell_pair(hfield,bfield,i,j,dx,dy,hface,bface,hnd,htd,bnd,btd, &
        theta,was_signed,ok,cacheh,cacheb)
        real(KREAL),intent(in)::hfield(:,:),bfield(:,:),dx,dy
        integer,intent(in)::i,j
        real(KREAL),intent(out)::hface(4,3),bface(4,3),hnd(4,3),htd(4,3),bnd(4,3),btd(4,3),theta
        logical,intent(out)::was_signed,ok
        type(FirstLayerPlane),intent(in),optional::cacheh,cacheb
        real(KREAL)::hvol(3,3),bvol(3,3),hm,bm,den
        integer::s,a,c
        logical::lok
        hm=hfield(i,j);bm=bfield(i,j);theta=1.0_KREAL;was_signed=.false.;ok=.false.
        if(present(cacheh).neqv.present(cacheb))return
        if(present(cacheh))then
            if(S2W1_wall_mode>=1)then
                call s2d1_cell_faces9_cached(cacheh,i,j,hface,hnd,htd,hvol,lok,hfield,.true.);if(.not.lok)return
                call s2d1_cell_faces9_cached(cacheb,i,j,bface,bnd,btd,bvol,lok,bfield,.true.);if(.not.lok)return
            else
                call s2d1_cell_faces9_cached(cacheh,i,j,hface,hnd,htd,hvol,lok);if(.not.lok)return
                call s2d1_cell_faces9_cached(cacheb,i,j,bface,bnd,btd,bvol,lok);if(.not.lok)return
            endif
        else
            if(S2W1_wall_mode>=1)then
                call s2d1_cell_faces9_wall(hfield,i,j,hface,hnd,htd,hvol,lok);if(.not.lok)return
                call s2d1_cell_faces9_wall(bfield,i,j,bface,bnd,btd,bvol,lok);if(.not.lok)return
            else
                call s2d1_cell_faces9(hfield,i,j,hface,hnd,htd,hvol,lok);if(.not.lok)return
                call s2d1_cell_faces9(bfield,i,j,bface,bnd,btd,bvol,lok);if(.not.lok)return
            endif
        endif
        call limit_cell_outputs(hm,bm,dx,dy,hface,bface,hnd,htd,bnd,btd,hvol,bvol,theta,was_signed,ok)
    end subroutine reconstruct_cell_pair

    subroutine limit_cell_outputs(hm,bm,dx,dy,hface,bface,hnd,htd,bnd,btd,hvol,bvol,theta,was_signed,ok)
        real(KREAL),intent(in)::hm,bm,dx,dy,hvol(3,3),bvol(3,3)
        real(KREAL),intent(inout)::hface(4,3),bface(4,3),hnd(4,3),htd(4,3),bnd(4,3),btd(4,3)
        real(KREAL),intent(out)::theta
        logical,intent(out)::was_signed,ok
        real(KREAL)::den
        integer::s,a,c
        theta=1.0_KREAL;was_signed=.false.;ok=.false.
        if(hm<=0.0_KREAL.or.bm<=0.0_KREAL)then
            theta=0.0_KREAL;was_signed=.true.
        else
            do s=1,4;do a=1,3
                if(hface(s,a)<0.0_KREAL)then;den=hm-hface(s,a);theta=min(theta,(1.0_KREAL-micro_margin)*hm/den);endif
                if(bface(s,a)<0.0_KREAL)then;den=bm-bface(s,a);theta=min(theta,(1.0_KREAL-micro_margin)*bm/den);endif
            enddo;enddo
            do c=1,3;do a=1,3
                if(hvol(a,c)<0.0_KREAL)then;den=hm-hvol(a,c);theta=min(theta,(1.0_KREAL-micro_margin)*hm/den);endif
                if(bvol(a,c)<0.0_KREAL)then;den=bm-bvol(a,c);theta=min(theta,(1.0_KREAL-micro_margin)*bm/den);endif
            enddo;enddo
            theta=max(0.0_KREAL,min(1.0_KREAL,theta))
        endif
        hface=hm+theta*(hface-hm);bface=bm+theta*(bface-bm)
        hnd=theta*hnd;htd=theta*htd;bnd=theta*bnd;btd=theta*btd
        hnd(1:2,:)=hnd(1:2,:)/dx;htd(1:2,:)=htd(1:2,:)/dy
        bnd(1:2,:)=bnd(1:2,:)/dx;btd(1:2,:)=btd(1:2,:)/dy
        hnd(3:4,:)=hnd(3:4,:)/dy;htd(3:4,:)=htd(3:4,:)/dx
        bnd(3:4,:)=bnd(3:4,:)/dy;btd(3:4,:)=btd(3:4,:)/dx
        ok=all(ieee_is_finite(hface)).and.all(ieee_is_finite(bface)).and. &
            all(ieee_is_finite(hnd)).and.all(ieee_is_finite(htd)).and. &
            all(ieee_is_finite(bnd)).and.all(ieee_is_finite(btd)).and.ieee_is_finite(theta)
    end subroutine


    subroutine reconstruct_face(h,b,axis,face,tangent,dx,dy,hl,hr,bl,br,hnl,hnr,htl,htr, &
        bnl,bnr,btl,btr,limited,signed,ok)
        real(KREAL),intent(in)::h(:,:,:),b(:,:,:),dx,dy
        integer,intent(in)::axis,face,tangent
        real(KREAL),intent(out)::hl(:,:),hr(:,:),bl(:,:),br(:,:),hnl(:,:),hnr(:,:),htl(:,:),htr(:,:)
        real(KREAL),intent(out)::bnl(:,:),bnr(:,:),btl(:,:),btr(:,:)
        integer,intent(inout)::limited,signed
        logical,intent(out)::ok
        real(KREAL)::x1(3),x2(3),n1(3),n2(3),t1(3),t2(3),dn,dt
        real(KREAL),parameter::gn(3)=[-0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL),0.0_KREAL, &
            0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL)]
        integer::k,a,n,right_cell
        logical::lok
        n=merge(size(h,1),size(h,2),axis==1);right_cell=ip(face+1,n)
        dn=merge(dx,dy,axis==1);dt=merge(dy,dx,axis==1)
        do k=1,size(h,3)
            if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
                call vl_side3(h(:,:,k),axis,face,tangent,0.5_KREAL,gn,x1,n1,t1)
                call vl_side3(h(:,:,k),axis,right_cell,tangent,-0.5_KREAL,gn,x2,n2,t2)
                lok=.true.
            else
                call s2d1_face3(h(:,:,k),axis,face,tangent,x1,x2,n1,n2,t1,t2,lok)
            endif
            if(.not.lok)then;ok=.false.;return;endif
            hl(k,:)=x1;hr(k,:)=x2;hnl(k,:)=n1/dn;hnr(k,:)=n2/dn;htl(k,:)=t1/dt;htr(k,:)=t2/dt
            if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
                call vl_side3(b(:,:,k),axis,face,tangent,0.5_KREAL,gn,x1,n1,t1)
                call vl_side3(b(:,:,k),axis,right_cell,tangent,-0.5_KREAL,gn,x2,n2,t2)
                lok=.true.
            else
                call s2d1_face3(b(:,:,k),axis,face,tangent,x1,x2,n1,n2,t1,t2,lok)
            endif
            if(.not.lok)then;ok=.false.;return;endif
            bl(k,:)=x1;br(k,:)=x2;bnl(k,:)=n1/dn;bnr(k,:)=n2/dn;btl(k,:)=t1/dt;btr(k,:)=t2/dt
            do a=1,3
                if(axis==1)then
                    call limit_pair(h(face,tangent,k),b(face,tangent,k),hl(k,a),bl(k,a), &
                        hnl(k,a),htl(k,a),bnl(k,a),btl(k,a),limited,signed)
                    call limit_pair(h(right_cell,tangent,k),b(right_cell,tangent,k),hr(k,a),br(k,a), &
                        hnr(k,a),htr(k,a),bnr(k,a),btr(k,a),limited,signed)
                else
                    call limit_pair(h(tangent,face,k),b(tangent,face,k),hl(k,a),bl(k,a), &
                        hnl(k,a),htl(k,a),bnl(k,a),btl(k,a),limited,signed)
                    call limit_pair(h(tangent,right_cell,k),b(tangent,right_cell,k),hr(k,a),br(k,a), &
                        hnr(k,a),htr(k,a),bnr(k,a),btr(k,a),limited,signed)
                endif
            end do
        end do
        ok=.true.
    contains
        real(KREAL) function mm(a,b)
            real(KREAL),intent(in)::a,b
            if(a*b<=0.0_KREAL)then;mm=0.0_KREAL
            else;mm=sign(min(abs(a),abs(b)),a);endif
        end function mm
        subroutine vl_side3(f,ax,c,t,side,eta,val,dn0,dt0)
            real(KREAL),intent(in)::f(:,:),side,eta(3)
            integer,intent(in)::ax,c,t
            real(KREAL),intent(out)::val(3),dn0(3),dt0(3)
            real(KREAL)::sn,st,center
            integer::a,nx,ny
            nx=size(f,1);ny=size(f,2)
            if(ax==1)then
                center=f(c,t);sn=mm(center-f(ip(c-1,nx),t),f(ip(c+1,nx),t)-center)
                st=mm(center-f(c,ip(t-1,ny)),f(c,ip(t+1,ny))-center)
            else
                center=f(t,c);sn=mm(center-f(t,ip(c-1,ny)),f(t,ip(c+1,ny))-center)
                st=mm(center-f(ip(t-1,nx),c),f(ip(t+1,nx),c)-center)
            endif
            do a=1,3;val(a)=center+side*sn+eta(a)*st;dn0(a)=sn;dt0(a)=st;enddo
        end subroutine vl_side3
    end subroutine reconstruct_face

    subroutine ugks_node(axis,dt,hl,hr,bl,br,hnl,hnr,htl,htr,bnl,bnr,btl,btr, &
        macro_nl,macro_nr,macro_t,u,v,weight,fh,fb,fw,face_moment_error,ok,f0_h,f0_b,k_h,k_b)
        integer,intent(in)::axis
        real(KREAL),intent(in)::dt,hl(:),hr(:),bl(:),br(:),hnl(:),hnr(:),htl(:),htr(:)
        real(KREAL),intent(in)::bnl(:),bnr(:),btl(:),btr(:),u(:),v(:),weight(:)
        real(KREAL),intent(in)::macro_nl(4),macro_nr(4),macro_t(4)
        real(KREAL),intent(out)::fh(:),fb(:),fw(4),face_moment_error
        logical,intent(out)::ok
        real(KREAL),intent(out),optional::f0_h(:),f0_b(:),k_h(:),k_b(:)
        real(KREAL)::f0h(size(u)),f0b(size(u)),fnh(size(u)),fnb(size(u)),fth(size(u)),ftb(size(u))
        real(KREAL)::w0(4),q0(2),wtime(4),kap(2,2),closure_res,dres
        real(KREAL)::mh(size(u)),mb(size(u)),kxh(size(u)),kxb(size(u)),kyh(size(u)),kyb(size(u))
        real(KREAL)::eah(size(u)),eab(size(u)),s0h(size(u)),s0b(size(u))
        real(KREAL)::mnlh(size(u)),mnlb(size(u)),mnrh(size(u)),mnrb(size(u))
        real(KREAL)::mth(size(u)),mtb(size(u)),mtimeh(size(u)),mtimeb(size(u)),vn,vt
        real(KREAL)::mt(5),slopeh,slopeb,mom(4),ex,tau_face
        logical::lok
        integer::k
        type(ClosureFirstOrderContext)::closure
        do k=1,size(u)
            vn=merge(u(k),v(k),axis==1);vt=merge(v(k),u(k),axis==1)
            if(vn>=0.0_KREAL)then
                f0h(k)=hl(k);f0b(k)=bl(k);fnh(k)=hnl(k);fnb(k)=bnl(k)
                fth(k)=htl(k);ftb(k)=btl(k)
            else
                f0h(k)=hr(k);f0b(k)=br(k);fnh(k)=hnr(k);fnb(k)=bnr(k)
                fth(k)=htr(k);ftb(k)=btr(k)
            endif
        end do
        call moment_2d(f0h,f0b,u,v,weight,w0,lok);if(.not.lok)then;ok=.false.;return;endif
        call closure_first_prepare(closure,w0,ck,u,v,weight,lok);if(.not.lok)then;ok=.false.;return;endif
        call make_mke_prepared(closure,w0,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kap,lok,closure_res)
        if(.not.lok)then;ok=.false.;return;endif
        call heat(w0,f0h-eah,f0b-eab,u,v,weight,q0,lok);if(.not.lok)then;ok=.false.;return;endif
        s0h=mh+kxh*q0(1)+kyh*q0(2);s0b=mb+kxb*q0(1)+kyb*q0(2)
        call m_direction_prepared(closure,w0,macro_nl,u,v,weight,mnlh,mnlb,lok,dres);if(.not.lok)then;ok=.false.;return;endif
        call m_direction_prepared(closure,w0,macro_nr,u,v,weight,mnrh,mnrb,lok,dres);if(.not.lok)then;ok=.false.;return;endif
        call m_direction_prepared(closure,w0,macro_t,u,v,weight,mth,mtb,lok,dres);if(.not.lok)then;ok=.false.;return;endif
        wtime=0.0_KREAL
        do k=1,size(u)
            vn=merge(u(k),v(k),axis==1);vt=merge(v(k),u(k),axis==1)
            if(vn>=0.0_KREAL)then
                fnh(k)=mnlh(k)
                fnb(k)=mnlb(k)
            else
                fnh(k)=mnrh(k)
                fnb(k)=mnrb(k)
            endif
            fnh(k)=vn*fnh(k)+vt*mth(k)
            fnb(k)=vn*fnb(k)+vt*mtb(k)
        end do
        call moment_2d(fnh,fnb,u,v,weight,wtime,lok);if(.not.lok)then;ok=.false.;return;endif
        wtime=-wtime
        call m_direction_prepared(closure,w0,wtime,u,v,weight,mtimeh,mtimeb,lok,dres);if(.not.lok)then;ok=.false.;return;endif
        if(S2PhysicalGas==1)then
        tau_face=s2_tau(w0,2)
        ex=exp(-dt/tau_face);mt(4)=tau_face*(1.0_KREAL-ex)
        mt(5)=-tau_face*dt*ex+tau_face*mt(4);mt(1)=dt-mt(4)
        mt(2)=-tau_face*mt(1)+mt(5);mt(3)=0.5_KREAL*dt*dt-tau_face*mt(1)
        else
            tau_face=tau0
        ex=exp(-dt/tau0);mt(4)=tau0*(1.0_KREAL-ex)
        mt(5)=-tau0*dt*ex+tau0*mt(4);mt(1)=dt-mt(4)
        mt(2)=-tau0*mt(1)+mt(5);mt(3)=0.5_KREAL*dt*dt-tau0*mt(1)
        endif
        do k=1,size(u)
            vn=merge(u(k),v(k),axis==1);vt=merge(v(k),u(k),axis==1)
            if(vn>=0.0_KREAL)then
                slopeh=vn*mnlh(k)+vt*mth(k);slopeb=vn*mnlb(k)+vt*mtb(k)
            else
                slopeh=vn*mnrh(k)+vt*mth(k);slopeb=vn*mnrb(k)+vt*mtb(k)
            endif
            fh(k)=vn*(mt(1)*s0h(k)+mt(2)*slopeh+mt(3)*mtimeh(k)+mt(4)*f0h(k) &
                -mt(5)*(vn*fnh_from(k,hnl,hnr,u,v,axis)+vt*fth(k)))
            fb(k)=vn*(mt(1)*s0b(k)+mt(2)*slopeb+mt(3)*mtimeb(k)+mt(4)*f0b(k) &
                -mt(5)*(vn*fnb_from(k,bnl,bnr,u,v,axis)+vt*ftb(k)))
            if(present(f0_h))then
                f0_h(k)=vn*f0h(k);f0_b(k)=vn*f0b(k)
                k_h(k)=vn*((s0h(k)-f0h(k))/tau_face-vn*fnh_from(k,hnl,hnr,u,v,axis)-vt*fth(k))
                k_b(k)=vn*((s0b(k)-f0b(k))/tau_face-vn*fnb_from(k,bnl,bnr,u,v,axis)-vt*ftb(k))
            endif
        end do
        call moment_2d(fh,fb,u,v,weight,fw,lok)
        if(.not.lok)then;ok=.false.;return;endif
        call moment_2d(fh,fb,u,v,weight,mom,lok)
        face_moment_error=maxval(abs(mom-fw))/max(1.0_KREAL,maxval(abs(fw)))
        ok=face_moment_error<=gate.and.all(ieee_is_finite(fh)).and.all(ieee_is_finite(fb))
    contains
        real(KREAL) function fnh_from(k,a,b,u,v,axis)
            integer,intent(in)::k,axis
            real(KREAL),intent(in)::a(:),b(:),u(:),v(:)
            real(KREAL)::vel
            vel=merge(u(k),v(k),axis==1);fnh_from=merge(a(k),b(k),vel>=0.0_KREAL)
        end function fnh_from
        real(KREAL) function fnb_from(k,a,b,u,v,axis)
            integer,intent(in)::k,axis
            real(KREAL),intent(in)::a(:),b(:),u(:),v(:)
            real(KREAL)::vel
            vel=merge(u(k),v(k),axis==1);fnb_from=merge(a(k),b(k),vel>=0.0_KREAL)
        end function fnb_from
    end subroutine ugks_node

    subroutine ugks_node_vl(axis,dt,dn,wl,wr,hl,hr,bl,br,hnl,hnr,htl,htr,bnl,bnr,btl,btr, &
        u,v,weight,fh,fb,fw,face_moment_error,ok)
        integer,intent(in)::axis
        real(KREAL),intent(in)::dt,dn,wl(4),wr(4),hl(:),hr(:),bl(:),br(:),hnl(:),hnr(:)
        real(KREAL),intent(in)::htl(:),htr(:),bnl(:),bnr(:),btl(:),btr(:),u(:),v(:),weight(:)
        real(KREAL),intent(out)::fh(:),fb(:),fw(4),face_moment_error
        logical,intent(out)::ok
        real(KREAL)::vn(size(u)),vt(size(u)),f0h(size(u)),f0b(size(u)),snh(size(u)),snb(size(u))
        real(KREAL)::sth(size(u)),stb(size(u)),mh(size(u)),mb(size(u)),th(size(u)),tb(size(u))
        real(KREAL)::w0(4),w0l(4),wll(4),wrr(4),swt(4),q(2),prim(4),al(4),ar(4),as(4),at(4)
        real(KREAL)::mu(0:6),mul(0:6),mur(0:6),mv(0:5),mxi(0:2),mt(5),ex
        real(KREAL)::mau0(4),maul(4),maur(4),mbv(4),maut(4),fplus(4),ffree(4),fslope(4)
        real(KREAL)::c2,haL,haR,hb0,haT,baL,baR,bb0,baT,slopeh,slopeb,mom(4)
        logical::lok
        integer::k
        do k=1,size(u)
            if(axis==1)then;vn(k)=u(k);vt(k)=v(k)
            else;vn(k)=v(k);vt(k)=-u(k);endif
            if(vn(k)>=0.0_KREAL)then
                f0h(k)=hl(k);f0b(k)=bl(k);snh(k)=hnl(k);snb(k)=bnl(k)
                sth(k)=htl(k);stb(k)=btl(k)
            else
                f0h(k)=hr(k);f0b(k)=br(k);snh(k)=hnr(k);snb(k)=bnr(k)
                sth(k)=htr(k);stb(k)=btr(k)
            endif
            if(axis==2)then;sth(k)=-sth(k);stb(k)=-stb(k);endif
        enddo
        call moment_2d(f0h,f0b,vn,vt,weight,w0l,lok);if(.not.lok)then;ok=.false.;return;endif
        if(axis==1)then
            wll=wl;wrr=wr
        else
            wll=[wl(1),wl(3),-wl(2),wl(4)];wrr=[wr(1),wr(3),-wr(2),wr(4)]
        endif
        call primitive(w0l,prim(1),prim(2),prim(3),ex,mt(1),lok)
        if(.not.lok)then;ok=.false.;return;endif
        prim(4)=0.5_KREAL/ex
        al=old_micro_slope(prim,(w0l-wll)/(0.5_KREAL*dn))
        ar=old_micro_slope(prim,(wrr-w0l)/(0.5_KREAL*dn))
        call moment_2d(sth,stb,vn,vt,weight,swt,lok);if(.not.lok)then;ok=.false.;return;endif
        as=old_micro_slope(prim,swt)
        call old_calc_moment(prim,mu,mv,mxi,mul,mur)
        maul=old_moment_auvxi(al,mul,mv,mxi,1,0)
        maur=old_moment_auvxi(ar,mur,mv,mxi,1,0)
        mbv=old_moment_auvxi(as,mu,mv,mxi,0,1)
        at=old_micro_slope(prim,-prim(1)*(maul+maur+mbv))
        call old_target_at(w0l,[0.0_KREAL,0.0_KREAL],vn,vt,mh,mb,th,tb,lok)
        if(.not.lok)then;ok=.false.;return;endif
        call heat(w0l,f0h,f0b,vn,vt,weight,q,lok);if(.not.lok)then;ok=.false.;return;endif
        call old_target_at(w0l,q,vn,vt,th,tb,mh,mb,lok);if(.not.lok)then;ok=.false.;return;endif
        th=th-mh;tb=tb-mb
        ex=exp(-dt/tau0);mt(4)=tau0*(1.0_KREAL-ex)
        mt(5)=-tau0*dt*ex+tau0*mt(4);mt(1)=dt-mt(4)
        mt(2)=-tau0*mt(1)+mt(5);mt(3)=0.5_KREAL*dt*dt-tau0*mt(1)
        mau0=old_moment_uvxi(mu,mv,mxi,1,0,0)
        maul=old_moment_auvxi(al,mul,mv,mxi,2,0)
        maur=old_moment_auvxi(ar,mur,mv,mxi,2,0)
        mbv=old_moment_auvxi(as,mu,mv,mxi,1,1)
        maut=old_moment_auvxi(at,mu,mv,mxi,1,0)
        fw=mt(1)*prim(1)*mau0+mt(2)*prim(1)*(maul+maur+mbv)+mt(3)*prim(1)*maut
        fplus=0.0_KREAL;ffree=0.0_KREAL;fslope=0.0_KREAL
        do k=1,size(u)
            c2=vn(k)**2+vt(k)**2;slopeh=vn(k)*snh(k)+vt(k)*sth(k)
            slopeb=vn(k)*snb(k)+vt(k)*stb(k)
            fplus=fplus+weight(k)*vn(k)*[th(k),vn(k)*th(k),vt(k)*th(k), &
                0.5_KREAL*(c2*th(k)+tb(k))]
            ffree=ffree+weight(k)*vn(k)*[f0h(k),vn(k)*f0h(k),vt(k)*f0h(k), &
                0.5_KREAL*(c2*f0h(k)+f0b(k))]
            fslope=fslope+weight(k)*vn(k)*[slopeh,vn(k)*slopeh,vt(k)*slopeh, &
                0.5_KREAL*(c2*slopeh+slopeb)]
            haL=al(1)*mh(k)+al(2)*vn(k)*mh(k)+al(3)*vt(k)*mh(k)+ &
                0.5_KREAL*al(4)*(c2*mh(k)+mb(k))
            haR=ar(1)*mh(k)+ar(2)*vn(k)*mh(k)+ar(3)*vt(k)*mh(k)+ &
                0.5_KREAL*ar(4)*(c2*mh(k)+mb(k))
            hb0=as(1)*mh(k)+as(2)*vn(k)*mh(k)+as(3)*vt(k)*mh(k)+ &
                0.5_KREAL*as(4)*(c2*mh(k)+mb(k))
            haT=at(1)*mh(k)+at(2)*vn(k)*mh(k)+at(3)*vt(k)*mh(k)+ &
                0.5_KREAL*at(4)*(c2*mh(k)+mb(k))
            baL=al(1)*mb(k)+al(2)*vn(k)*mb(k)+al(3)*vt(k)*mb(k)+ &
                0.5_KREAL*al(4)*(c2*mb(k)+mxi(2)*mh(k))
            baR=ar(1)*mb(k)+ar(2)*vn(k)*mb(k)+ar(3)*vt(k)*mb(k)+ &
                0.5_KREAL*ar(4)*(c2*mb(k)+mxi(2)*mh(k))
            bb0=as(1)*mb(k)+as(2)*vn(k)*mb(k)+as(3)*vt(k)*mb(k)+ &
                0.5_KREAL*as(4)*(c2*mb(k)+mxi(2)*mh(k))
            baT=at(1)*mb(k)+at(2)*vn(k)*mb(k)+at(3)*vt(k)*mb(k)+ &
                0.5_KREAL*at(4)*(c2*mb(k)+mxi(2)*mh(k))
            if(vn(k)>=0.0_KREAL)then
                fh(k)=mt(1)*vn(k)*(mh(k)+th(k))+mt(2)*vn(k)**2*haL+ &
                    mt(2)*vn(k)*vt(k)*hb0+mt(3)*vn(k)*haT+mt(4)*vn(k)*f0h(k)-mt(5)*vn(k)*slopeh
                fb(k)=mt(1)*vn(k)*(mb(k)+tb(k))+mt(2)*vn(k)**2*baL+ &
                    mt(2)*vn(k)*vt(k)*bb0+mt(3)*vn(k)*baT+mt(4)*vn(k)*f0b(k)-mt(5)*vn(k)*slopeb
            else
                fh(k)=mt(1)*vn(k)*(mh(k)+th(k))+mt(2)*vn(k)**2*haR+ &
                    mt(2)*vn(k)*vt(k)*hb0+mt(3)*vn(k)*haT+mt(4)*vn(k)*f0h(k)-mt(5)*vn(k)*slopeh
                fb(k)=mt(1)*vn(k)*(mb(k)+tb(k))+mt(2)*vn(k)**2*baR+ &
                    mt(2)*vn(k)*vt(k)*bb0+mt(3)*vn(k)*baT+mt(4)*vn(k)*f0b(k)-mt(5)*vn(k)*slopeb
            endif
        enddo
        fw=fw+mt(1)*fplus+mt(4)*ffree-mt(5)*fslope
        if(axis==2)fw=[fw(1),-fw(3),fw(2),fw(4)]
        call moment_2d(fh,fb,u,v,weight,mom,lok)
        if(axis==2)then
            ! fh/fb are indexed on the original global velocity table, while their
            ! values were evaluated in the rotated local frame.
            mom=[mom(1),mom(2),mom(3),mom(4)]
        endif
        face_moment_error=maxval(abs(mom-fw))/max(1.0_KREAL,maxval(abs(fw)))
        ok=all(ieee_is_finite(fh)).and.all(ieee_is_finite(fb)).and.all(ieee_is_finite(fw))
    end subroutine ugks_node_vl

    subroutine compute_flux(w,h,b,u,v,weight,dx,dy,dt,fxw,fyw,fxh,fxb,fyh,fyb, &
        face_error,limited,signed,ok,lower_w,lower_h,lower_b, &
        leading_xh,leading_xb,leading_yh,leading_yb,time_xh,time_xb,time_yh,time_yb,lower0h,lower0b,lowerKh,lowerKb)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),weight(:),dx,dy,dt
        real(KREAL),intent(out)::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:)
        real(KREAL),intent(out)::face_error
        integer,intent(out)::limited,signed
        logical,intent(out)::ok
        real(KREAL),intent(out),optional::lower_w(:,:),lower_h(:,:),lower_b(:,:)
        real(KREAL),intent(out),optional::leading_xh(:,:,:),leading_xb(:,:,:),leading_yh(:,:,:),leading_yb(:,:,:)
        real(KREAL),intent(out),optional::time_xh(:,:,:),time_xb(:,:,:),time_yh(:,:,:),time_yb(:,:,:)
        real(KREAL),intent(out),optional::lower0h(:,:),lower0b(:,:),lowerKh(:,:),lowerKb(:,:)
        integer::nx,ny,nvel,i,j,a,m
        real(KREAL)::hl(size(u),3),hr(size(u),3),bl(size(u),3),br(size(u),3)
        real(KREAL)::hnl(size(u),3),hnr(size(u),3),htl(size(u),3),htr(size(u),3)
        real(KREAL)::bnl(size(u),3),bnr(size(u),3),btl(size(u),3),btr(size(u),3)
        real(KREAL)::fh(size(u)),fb(size(u)),fw(4),err
        real(KREAL)::ml(3),mr(3),mnl(3),mnr(3),mtl(3),mtr(3)
        real(KREAL)::macro_nl(4,3),macro_nr(4,3),macro_t(4,3),theta
        real(KREAL),allocatable::ch(:,:,:,:,:),cb(:,:,:,:,:),chnd(:,:,:,:,:),chtd(:,:,:,:,:)
        real(KREAL),allocatable::cbnd(:,:,:,:,:),cbtd(:,:,:,:,:)
        real(KREAL)::thf(4,3),tbf(4,3),thn(4,3),tht(4,3),tbn(4,3),tbt(4,3)
        real(KREAL)::nf0h(size(u)),nf0b(size(u)),nkh(size(u)),nkb(size(u))
        real(KREAL)::w2_start,w2_end
        logical::lok,was_signed,is_s2,want_coeff,want_wallcoeff
        integer::k,il,ih,jl,jh
        integer,allocatable::owner_counts(:,:,:)
        logical::use_owner,recon_ok,use_scope,use_tile,profile_tile
        type(VelocityTileWorkspace)::tile
        integer::k0,nt,kv,fhidx,fbidx,tile_unit
        real(KREAL)::tile_t0,tile_t1,tile_prepare_s,tile_store_s
        character(len=80)::tile_filename
        logical::tile_exists
        character(len=16)::scope_setting
        integer::scope_status
        real(KREAL)::owner_start
        type(FirstLayerPlane)::cacheh,cacheb
        nx=size(h,1);ny=size(h,2);nvel=size(h,3)
        want_coeff=present(leading_xh);want_wallcoeff=present(lower0h)
        if(want_wallcoeff)then
            if(.not.(want_coeff.and.present(lower0b).and.present(lowerKh).and.present(lowerKb)))then
                ok=.false.;return
            endif
            lower0h=0;lower0b=0;lowerKh=0;lowerKb=0
        endif
        if(want_coeff)then
            if(.not.(present(leading_xb).and.present(leading_yh).and.present(leading_yb).and. &
                present(time_xh).and.present(time_xb).and.present(time_yh).and.present(time_yb)))then
                ok=.false.;return
            endif
            leading_xh=0.0_KREAL;leading_xb=0.0_KREAL;leading_yh=0.0_KREAL;leading_yb=0.0_KREAL
            time_xh=0.0_KREAL;time_xb=0.0_KREAL;time_yh=0.0_KREAL;time_yb=0.0_KREAL
        endif
        ! This selects the spatial S2 reconstruction only.  IsS2CN2 retains
        ! its separate meaning for the physical-time CN update dispatcher.
        is_s2=trim(adjustl(TimeIntegration_method))=='S2D1'.or.IsS2CN2().or.IsS2Steady()
        fxw=0.0_KREAL;fyw=0.0_KREAL;fxh=0.0_KREAL;fxb=0.0_KREAL
        fyh=0.0_KREAL;fyb=0.0_KREAL;face_error=0.0_KREAL;limited=0;signed=0
        w2_reconstruction_s=0.0_KREAL;w2_wall_s=0.0_KREAL;w2_wall_calls=0
        call cpu_time(w2_start)
        ! Local cache order: side,node,velocity,x,y. Arithmetic and traversal unchanged.
        if(is_s2)then
            allocate(ch(4,3,nvel,nx,ny),cb(4,3,nvel,nx,ny),chnd(4,3,nvel,nx,ny), &
                chtd(4,3,nvel,nx,ny),cbnd(4,3,nvel,nx,ny),cbtd(4,3,nvel,nx,ny))
            il=xrecon_lo();ih=xrecon_hi(nx);jl=yrecon_lo();jh=recon_hi(ny)
            use_owner=.false.
            use_owner=owner_enabled()
            if(use_owner)then
                il=xlo();ih=xhi(nx);jl=active_lo();jh=active_hi(ny)
                allocate(owner_counts(2,nx,ny));owner_counts=0
            endif
            scope_setting='';call get_environment_variable('S2_FIRST_LAYER_SCOPE',scope_setting,status=scope_status)
            use_scope=use_owner.and.scope_status==0.and.trim(scope_setting)=='1'
            scope_setting='';call get_environment_variable('S2_VELOCITY_TILE',scope_setting,status=scope_status)
            use_tile=use_scope.and.scope_status==0.and.trim(scope_setting)=='1'
            scope_setting='';call get_environment_variable('S2_TILE_PROFILE',scope_setting,status=scope_status)
            profile_tile=scope_status==0.and.trim(scope_setting)=='1'
            tile%profile=profile_tile;tile_prepare_s=0;tile_store_s=0
            scope_setting='';call get_environment_variable('S2_WALL_TILE',scope_setting,status=scope_status)
            tile%wall_batch=use_tile.and.scope_status==0.and.trim(scope_setting)=='1'
            recon_ok=.true.;owner_start=MPI_Wtime()
            if(use_tile)then
                velocity_tiles:do k0=1,nvel,VT_VELOCITIES
                    nt=min(VT_VELOCITIES,nvel-k0+1)
                    if(profile_tile)call cpu_time(tile_t0)
                    call prepare_velocity_tile(tile,h,b,k0,nt,[il,ih,jl,jh],S2W1_wall_mode>=1,lok)
                    if(.not.lok)then;recon_ok=.false.;exit velocity_tiles;endif
                    if(profile_tile)then
                        call cpu_time(tile_t1);tile_prepare_s=tile_prepare_s+tile_t1-tile_t0
                    endif
                    do j=jl,jh;do i=il,ih
                        call reconstruct_velocity_tile_cell(tile,h,b,k0,i,j,S2W1_wall_mode>=1,lok)
                        if(.not.lok)then;recon_ok=.false.;exit velocity_tiles;endif
                        if(profile_tile)call cpu_time(tile_t0)
                        do kv=1,nt
                            k=k0+kv-1;fhidx=2*kv-1;fbidx=2*kv
                            call limit_cell_outputs(h(i,j,k),b(i,j,k),dx,dy, &
                                tile%face(:,:,fhidx),tile%face(:,:,fbidx),tile%dn(:,:,fhidx), &
                                tile%dt(:,:,fhidx),tile%dn(:,:,fbidx),tile%dt(:,:,fbidx), &
                                tile%vol(:,:,fhidx),tile%vol(:,:,fbidx),theta,was_signed,lok)
                            if(.not.lok)then;recon_ok=.false.;exit velocity_tiles;endif
                            ch(:,:,k,i,j)=tile%face(:,:,fhidx);cb(:,:,k,i,j)=tile%face(:,:,fbidx)
                            chnd(:,:,k,i,j)=tile%dn(:,:,fhidx);chtd(:,:,k,i,j)=tile%dt(:,:,fhidx)
                            cbnd(:,:,k,i,j)=tile%dn(:,:,fbidx);cbtd(:,:,k,i,j)=tile%dt(:,:,fbidx)
                            owner_cells=owner_cells+1
                            if(was_signed)owner_counts(2,i,j)=owner_counts(2,i,j)+1
                            if((.not.was_signed).and.theta<1.0_KREAL)owner_counts(1,i,j)=owner_counts(1,i,j)+1
                            if(was_signed)signed=signed+1
                            if((.not.was_signed).and.theta<1.0_KREAL)limited=limited+1
                        enddo
                        if(profile_tile)then
                            call cpu_time(tile_t1);tile_store_s=tile_store_s+tile_t1-tile_t0
                        endif
                    enddo;enddo
                enddo velocity_tiles
            else
            recon_vel: do k=1,nvel
                if(use_scope)then
                    call prepare_first_layer_plane(cacheh,h(:,:,k),lok,S2W1_wall_mode>=1,[il,ih,jl,jh])
                else
                    call prepare_first_layer_plane(cacheh,h(:,:,k),lok,S2W1_wall_mode>=1)
                endif
                if(.not.lok)then;recon_ok=.false.;exit recon_vel;endif
                if(use_scope)then
                    call prepare_first_layer_plane(cacheb,b(:,:,k),lok,S2W1_wall_mode>=1,[il,ih,jl,jh])
                else
                    call prepare_first_layer_plane(cacheb,b(:,:,k),lok,S2W1_wall_mode>=1)
                endif
                if(.not.lok)then;recon_ok=.false.;exit recon_vel;endif
            do j=jl,jh;do i=il,ih
                call reconstruct_cell_pair(h(:,:,k),b(:,:,k),i,j,dx,dy,thf,tbf,thn,tht,tbn,tbt, &
                    theta,was_signed,lok,cacheh,cacheb)
                if(.not.lok)then;recon_ok=.false.;exit recon_vel;endif
                ch(:,:,k,i,j)=thf;cb(:,:,k,i,j)=tbf;chnd(:,:,k,i,j)=thn;chtd(:,:,k,i,j)=tht
                cbnd(:,:,k,i,j)=tbn;cbtd(:,:,k,i,j)=tbt
                owner_cells=owner_cells+1
                if(use_owner)then
                    if(was_signed)owner_counts(2,i,j)=owner_counts(2,i,j)+1
                    if((.not.was_signed).and.theta<1.0_KREAL)owner_counts(1,i,j)=owner_counts(1,i,j)+1
                endif
                if(was_signed)signed=signed+1
                if((.not.was_signed).and.theta<1.0_KREAL)limited=limited+1
            enddo;enddo
            enddo recon_vel
            endif
            if(use_tile.and.profile_tile)then
                write(tile_filename,'("velocity_tile_rank",I4.4,".csv")')myID
                inquire(file=trim(tile_filename),exist=tile_exists)
                open(newunit=tile_unit,file=trim(tile_filename),status='unknown',position='append')
                if(.not.tile_exists)write(tile_unit,'(A)') &
                    'prepare_s,pack_s,kernel_s,unpack_s,limit_store_s,wall_fallback_s,packed_values,stored_values,cells,fallback_cells'
                write(tile_unit,'(6(ES24.16,","),3(I0,","),I0)') &
                    tile_prepare_s,tile%pack_s,tile%kernel_s,tile%unpack_s,tile_store_s,tile%fallback_s, &
                    tile%packed_values,tile%stored_values,tile%cells,tile%fallback_cells
                close(tile_unit)
            endif
            owner_reconstruction_s=owner_reconstruction_s+MPI_Wtime()-owner_start
            if(use_owner)then
                call owner_ready(recon_ok)
                if(.not.recon_ok)then;ok=.false.;return;endif
                call owner_exchange(ch,cb,chnd,chtd,cbnd,cbtd,owner_counts)
                limited=sum(owner_counts(1,xrecon_lo():xrecon_hi(nx),yrecon_lo():recon_hi(ny)))
                signed=sum(owner_counts(2,xrecon_lo():xrecon_hi(nx),yrecon_lo():recon_hi(ny)))
            endif
            if(.not.recon_ok)then;ok=.false.;return;endif
        endif
        call cpu_time(w2_end);w2_reconstruction_s=w2_end-w2_start
        call cpu_time(w2_start)
        do j=yrecon_lo(),active_hi(ny)
            do i=xrecon_lo(),xhi(nx)
                if(.not.(cavity_mode().and.global_i(i)==global_nx(nx)))then
                if(is_s2)then
                    do m=1,4
                        call s2d1_face3(w(m,:,:),1,i,j,ml,mr,mnl,mnr,mtl,mtr,lok,S2W1_wall_mode>=1)
                        if(.not.lok)then;ok=.false.;return;endif
                        macro_nl(m,:)=mnl/dx;macro_nr(m,:)=mnr/dx
                        macro_t(m,:)=0.5_KREAL*(mtl+mtr)/dy
                    enddo
                    hl=transpose(ch(2,:,:,i,j));hr=transpose(ch(1,:,:,ip(i+1,nx),j))
                    bl=transpose(cb(2,:,:,i,j));br=transpose(cb(1,:,:,ip(i+1,nx),j))
                    hnl=transpose(chnd(2,:,:,i,j));hnr=transpose(chnd(1,:,:,ip(i+1,nx),j))
                    htl=transpose(chtd(2,:,:,i,j));htr=transpose(chtd(1,:,:,ip(i+1,nx),j))
                    bnl=transpose(cbnd(2,:,:,i,j));bnr=transpose(cbnd(1,:,:,ip(i+1,nx),j))
                    btl=transpose(cbtd(2,:,:,i,j));btr=transpose(cbtd(1,:,:,ip(i+1,nx),j))
                else
                    call reconstruct_face(h,b,1,i,j,dx,dy,hl,hr,bl,br,hnl,hnr,htl,htr, &
                        bnl,bnr,btl,btr,limited,signed,lok)
                    if(.not.lok)then;ok=.false.;return;endif
                endif
                if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
                    call ugks_node_vl(1,dt,dx,w(:,i,j),w(:,ip(i+1,nx),j),hl(:,2),hr(:,2),bl(:,2),br(:,2), &
                        hnl(:,2),hnr(:,2),htl(:,2),htr(:,2),bnl(:,2),bnr(:,2),btl(:,2),btr(:,2), &
                        u,v,weight,fh,fb,fw,err,lok)
                    if(.not.lok)then;ok=.false.;return;endif
                    fxh(i,j,:)=fh;fxb(i,j,:)=fb;fxw(:,i,j)=fw;face_error=max(face_error,err)
                else
                    do a=1,3
                        if(want_coeff)then
                            call ugks_node(1,dt,hl(:,a),hr(:,a),bl(:,a),br(:,a),hnl(:,a),hnr(:,a), &
                                htl(:,a),htr(:,a),bnl(:,a),bnr(:,a),btl(:,a),btr(:,a), &
                                macro_nl(:,a),macro_nr(:,a),macro_t(:,a),u,v,weight,fh,fb,fw,err,lok, &
                                nf0h,nf0b,nkh,nkb)
                            leading_xh(i,j,:)=leading_xh(i,j,:)+gw(a)*nf0h
                            leading_xb(i,j,:)=leading_xb(i,j,:)+gw(a)*nf0b
                            time_xh(i,j,:)=time_xh(i,j,:)+gw(a)*nkh
                            time_xb(i,j,:)=time_xb(i,j,:)+gw(a)*nkb
                        else
                        call ugks_node(1,dt,hl(:,a),hr(:,a),bl(:,a),br(:,a),hnl(:,a),hnr(:,a), &
                            htl(:,a),htr(:,a),bnl(:,a),bnr(:,a),btl(:,a),btr(:,a), &
                            macro_nl(:,a),macro_nr(:,a),macro_t(:,a),u,v,weight,fh,fb,fw,err,lok)
                        endif
                        if(.not.lok)then;ok=.false.;return;endif
                        fxh(i,j,:)=fxh(i,j,:)+gw(a)*fh;fxb(i,j,:)=fxb(i,j,:)+gw(a)*fb
                        fxw(:,i,j)=fxw(:,i,j)+gw(a)*fw;face_error=max(face_error,err)
                    enddo
                endif
                endif ! interior x face (physical right wall is integrated below)
                if(S2W1_wall_mode>=1.and.global_j(j)==global_ny(ny))cycle
                if(is_s2)then
                    do m=1,4
                        call s2d1_face3(w(m,:,:),2,j,i,ml,mr,mnl,mnr,mtl,mtr,lok,S2W1_wall_mode>=1)
                        if(.not.lok)then;ok=.false.;return;endif
                        macro_nl(m,:)=mnl/dy;macro_nr(m,:)=mnr/dy
                        macro_t(m,:)=0.5_KREAL*(mtl+mtr)/dx
                    enddo
                    hl=transpose(ch(4,:,:,i,j));hr=transpose(ch(3,:,:,i,ip(j+1,ny)))
                    bl=transpose(cb(4,:,:,i,j));br=transpose(cb(3,:,:,i,ip(j+1,ny)))
                    hnl=transpose(chnd(4,:,:,i,j));hnr=transpose(chnd(3,:,:,i,ip(j+1,ny)))
                    htl=transpose(chtd(4,:,:,i,j));htr=transpose(chtd(3,:,:,i,ip(j+1,ny)))
                    bnl=transpose(cbnd(4,:,:,i,j));bnr=transpose(cbnd(3,:,:,i,ip(j+1,ny)))
                    btl=transpose(cbtd(4,:,:,i,j));btr=transpose(cbtd(3,:,:,i,ip(j+1,ny)))
                else
                    call reconstruct_face(h,b,2,j,i,dx,dy,hl,hr,bl,br,hnl,hnr,htl,htr, &
                        bnl,bnr,btl,btr,limited,signed,lok)
                    if(.not.lok)then;ok=.false.;return;endif
                endif
                if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
                    call ugks_node_vl(2,dt,dy,w(:,i,j),w(:,i,ip(j+1,ny)),hl(:,2),hr(:,2),bl(:,2),br(:,2), &
                        hnl(:,2),hnr(:,2),htl(:,2),htr(:,2),bnl(:,2),bnr(:,2),btl(:,2),btr(:,2), &
                        u,v,weight,fh,fb,fw,err,lok)
                    if(.not.lok)then;ok=.false.;return;endif
                    fyh(i,j,:)=fh;fyb(i,j,:)=fb;fyw(:,i,j)=fw;face_error=max(face_error,err)
                else
                    do a=1,3
                        if(want_coeff)then
                            call ugks_node(2,dt,hl(:,a),hr(:,a),bl(:,a),br(:,a),hnl(:,a),hnr(:,a), &
                                htl(:,a),htr(:,a),bnl(:,a),bnr(:,a),btl(:,a),btr(:,a), &
                                macro_nl(:,a),macro_nr(:,a),macro_t(:,a),u,v,weight,fh,fb,fw,err,lok, &
                                nf0h,nf0b,nkh,nkb)
                            leading_yh(i,j,:)=leading_yh(i,j,:)+gw(a)*nf0h
                            leading_yb(i,j,:)=leading_yb(i,j,:)+gw(a)*nf0b
                            time_yh(i,j,:)=time_yh(i,j,:)+gw(a)*nkh
                            time_yb(i,j,:)=time_yb(i,j,:)+gw(a)*nkb
                        else
                        call ugks_node(2,dt,hl(:,a),hr(:,a),bl(:,a),br(:,a),hnl(:,a),hnr(:,a), &
                            htl(:,a),htr(:,a),bnl(:,a),bnr(:,a),btl(:,a),btr(:,a), &
                            macro_nl(:,a),macro_nr(:,a),macro_t(:,a),u,v,weight,fh,fb,fw,err,lok)
                        endif
                        if(.not.lok)then;ok=.false.;return;endif
                        fyh(i,j,:)=fyh(i,j,:)+gw(a)*fh;fyb(i,j,:)=fyb(i,j,:)+gw(a)*fb
                        fyw(:,i,j)=fyw(:,i,j)+gw(a)*fw;face_error=max(face_error,err)
                    enddo
                endif
            enddo
        enddo
        call cpu_time(w2_end);w2_closure_s=w2_end-w2_start
        if(cavity_mode())then
            if(.not.cart_context.or.want_coeff)then;ok=.false.;return;endif
            do j=active_lo(),active_hi(ny)
                do a=1,3
                    if(left_wall())then
                        i=xlo()
                        call s2w2_wall_node(-1,dt,dx,w([1,3,2,4],i,j),ch(1,a,:,i,j),cb(1,a,:,i,j), &
                            chnd(1,a,:,i,j),cbnd(1,a,:,i,j),chtd(1,a,:,i,j),cbtd(1,a,:,i,j), &
                            v,u,weight,0.0_KREAL,S2W1_wall_theta(1),fh,fb,fw,lok)
                        if(.not.lok)then;ok=.false.;return;endif
                        fxh(i-1,j,:)=fxh(i-1,j,:)+gw(a)*fh;fxb(i-1,j,:)=fxb(i-1,j,:)+gw(a)*fb
                        fxw(:,i-1,j)=fxw(:,i-1,j)+gw(a)*fw([1,3,2,4])
                    endif
                    if(right_wall())then
                        i=xhi(nx)
                        call s2w2_wall_node(1,dt,dx,w([1,3,2,4],i,j),ch(2,a,:,i,j),cb(2,a,:,i,j), &
                            chnd(2,a,:,i,j),cbnd(2,a,:,i,j),chtd(2,a,:,i,j),cbtd(2,a,:,i,j), &
                            v,u,weight,0.0_KREAL,S2W1_wall_theta(1),fh,fb,fw,lok)
                        if(.not.lok)then;ok=.false.;return;endif
                        fxh(i,j,:)=fxh(i,j,:)+gw(a)*fh;fxb(i,j,:)=fxb(i,j,:)+gw(a)*fb
                        fxw(:,i,j)=fxw(:,i,j)+gw(a)*fw([1,3,2,4])
                    endif
                enddo
            enddo
        endif
        if(S2W1_wall_mode>=1)then
            if(.not.(present(lower_w).and.present(lower_h).and.present(lower_b)))then
                ok=.false.;return
            endif
            lower_w=0.0_KREAL;lower_h=0.0_KREAL;lower_b=0.0_KREAL
            call cpu_time(w2_start)
            do i=xlo(),xhi(nx)
                do a=1,3
                    if(lower_wall())then
                    if(FluxWall_type==1)then
                        if(want_wallcoeff)then
                        call s2w2_wall_node(-1,dt,dy,w(:,i,active_lo()),ch(3,a,:,i,active_lo()),cb(3,a,:,i,active_lo()), &
                            chnd(3,a,:,i,active_lo()),cbnd(3,a,:,i,active_lo()),chtd(3,a,:,i,active_lo()),cbtd(3,a,:,i,active_lo()), &
                            u,v,weight,S2W1_wall_ux(1),S2W1_wall_theta(1),fh,fb,fw,lok, &
                            f0_h=nf0h,f0_b=nf0b,k_h=nkh,k_b=nkb)
                            lower0h(i,:)=lower0h(i,:)+gw(a)*nf0h
                            lower0b(i,:)=lower0b(i,:)+gw(a)*nf0b
                            lowerKh(i,:)=lowerKh(i,:)+gw(a)*nkh
                            lowerKb(i,:)=lowerKb(i,:)+gw(a)*nkb
                        else
                        call s2w2_wall_node(-1,dt,dy,w(:,i,active_lo()),ch(3,a,:,i,active_lo()),cb(3,a,:,i,active_lo()), &
                            chnd(3,a,:,i,active_lo()),cbnd(3,a,:,i,active_lo()),chtd(3,a,:,i,active_lo()),cbtd(3,a,:,i,active_lo()), &
                            u,v,weight,S2W1_wall_ux(1),S2W1_wall_theta(1),fh,fb,fw,lok)
                        endif
                    else
                    call s2w1_wall_node(-1,dt,ch(3,a,:,i,active_lo()),cb(3,a,:,i,active_lo()),u,v,weight, &
                        S2W1_wall_ux(1),S2W1_wall_theta(1),fh,fb,fw,lok)
                    endif
                    if(.not.lok)then;ok=.false.;return;endif
                    lower_h(i,:)=lower_h(i,:)+gw(a)*fh
                    lower_b(i,:)=lower_b(i,:)+gw(a)*fb
                    lower_w(:,i)=lower_w(:,i)+gw(a)*fw
                    endif
                    if(upper_wall())then
                    if(FluxWall_type==1)then
                        if(want_wallcoeff)then
                        call s2w2_wall_node(1,dt,dy,w(:,i,active_hi(ny)),ch(4,a,:,i,active_hi(ny)),cb(4,a,:,i,active_hi(ny)), &
                            chnd(4,a,:,i,active_hi(ny)),cbnd(4,a,:,i,active_hi(ny)),chtd(4,a,:,i,active_hi(ny)),cbtd(4,a,:,i,active_hi(ny)), &
                            u,v,weight,S2W1_wall_ux(2),S2W1_wall_theta(2),fh,fb,fw,lok, &
                            f0_h=nf0h,f0_b=nf0b,k_h=nkh,k_b=nkb)
                            leading_yh(i,active_hi(ny),:)=leading_yh(i,active_hi(ny),:)+gw(a)*nf0h
                            leading_yb(i,active_hi(ny),:)=leading_yb(i,active_hi(ny),:)+gw(a)*nf0b
                            time_yh(i,active_hi(ny),:)=time_yh(i,active_hi(ny),:)+gw(a)*nkh
                            time_yb(i,active_hi(ny),:)=time_yb(i,active_hi(ny),:)+gw(a)*nkb
                        else
                        call s2w2_wall_node(1,dt,dy,w(:,i,active_hi(ny)),ch(4,a,:,i,active_hi(ny)),cb(4,a,:,i,active_hi(ny)), &
                            chnd(4,a,:,i,active_hi(ny)),cbnd(4,a,:,i,active_hi(ny)),chtd(4,a,:,i,active_hi(ny)),cbtd(4,a,:,i,active_hi(ny)), &
                            u,v,weight,S2W1_wall_ux(2),S2W1_wall_theta(2),fh,fb,fw,lok)
                        endif
                    else
                    call s2w1_wall_node(1,dt,ch(4,a,:,i,active_hi(ny)),cb(4,a,:,i,active_hi(ny)),u,v,weight, &
                        S2W1_wall_ux(2),S2W1_wall_theta(2),fh,fb,fw,lok)
                    endif
                    if(.not.lok)then;ok=.false.;return;endif
                    fyh(i,active_hi(ny),:)=fyh(i,active_hi(ny),:)+gw(a)*fh
                    fyb(i,active_hi(ny),:)=fyb(i,active_hi(ny),:)+gw(a)*fb
                    fyw(:,i,active_hi(ny))=fyw(:,i,active_hi(ny))+gw(a)*fw
                    endif
                enddo
            enddo
            call cpu_time(w2_end);w2_wall_s=w2_end-w2_start;w2_wall_calls=3*nx*(merge(1,0,lower_wall())+merge(1,0,upper_wall()))
        endif
        ok=.true.
    end subroutine compute_flux

    ! Time-integrated free-streaming kinetic wall flux at one tangential GQ3 point.
    subroutine s2w1_wall_node(side,dt,hin,bin,u,v,weight,ux,theta,fh,fb,fw,ok)
        integer,intent(in)::side
        real(KREAL),intent(in)::dt,hin(:),bin(:),u(:),v(:),weight(:),ux,theta
        real(KREAL),intent(out)::fh(:),fb(:),fw(4)
        logical,intent(out)::ok
        real(KREAL)::lambda,incidence,reflection,mh(size(u)),rho_w
        integer::k
        lambda=0.5_KREAL/theta
        incidence=0.0_KREAL;reflection=0.0_KREAL
        do k=1,size(u)
            mh(k)=lambda/pi_d*exp(-lambda*((u(k)-ux)**2+v(k)**2))
            if(real(side,KREAL)*v(k)>0.0_KREAL)then
                incidence=incidence+weight(k)*v(k)*hin(k)
            else
                reflection=reflection+weight(k)*v(k)*mh(k)
            endif
        enddo
        ok=ieee_is_finite(reflection).and.abs(reflection)>tiny(reflection)
        if(.not.ok)return
        rho_w=-incidence/reflection
        ok=ieee_is_finite(rho_w).and.rho_w>0.0_KREAL
        if(.not.ok)return
        do k=1,size(u)
            if(real(side,KREAL)*v(k)>0.0_KREAL)then
                fh(k)=dt*v(k)*hin(k);fb(k)=dt*v(k)*bin(k)
            else
                fh(k)=dt*v(k)*rho_w*mh(k)
                fb(k)=real(ck,KREAL)*theta*fh(k)
            endif
        enddo
        call moment_2d(fh,fb,u,v,weight,fw,ok)
        if(.not.ok)return
        ok=abs(fw(1))<=gate*max(1.0_KREAL,dt*abs(incidence))
    end subroutine s2w1_wall_node

    ! Wall-only reevaluation for a temporary midpoint state.  It deliberately
    ! calls the same cell-local reconstruction, h/b limiter and diffuse kernel
    ! as the full W1 face path, without recomputing any interior UGKS faces.
    subroutine s2w1_compute_wall_flux(h,b,u,v,weight,dx,dy,dt, &
        lower_w,lower_h,lower_b,upper_w,upper_h,upper_b,limited,signed,ok)
        real(KREAL),intent(in)::h(:,:,:),b(:,:,:),u(:),v(:),weight(:),dx,dy,dt
        real(KREAL),intent(out)::lower_w(:,:),lower_h(:,:),lower_b(:,:)
        real(KREAL),intent(out)::upper_w(:,:),upper_h(:,:),upper_b(:,:)
        integer,intent(out)::limited,signed
        logical,intent(out)::ok
        real(KREAL)::lo_h(size(u),3),lo_b(size(u),3),up_h(size(u),3),up_b(size(u),3)
        real(KREAL)::hf(4,3),bf(4,3),hn(4,3),ht(4,3),bn(4,3),bt(4,3)
        real(KREAL)::fh(size(u)),fb(size(u)),fw(4),theta
        integer::nx,ny,nvel,i,k,a
        logical::lok,was_signed
        nx=size(h,1);ny=size(h,2);nvel=size(h,3)
        lower_w=0.0_KREAL;lower_h=0.0_KREAL;lower_b=0.0_KREAL
        upper_w=0.0_KREAL;upper_h=0.0_KREAL;upper_b=0.0_KREAL
        limited=0;signed=0;ok=.false.
        if(S2W1_wall_mode/=1.or.size(u)/=nvel.or.size(b,3)/=nvel)return
        do i=1,nx
            do k=1,nvel
                call reconstruct_cell_pair(h(:,:,k),b(:,:,k),i,1,dx,dy,hf,bf,hn,ht,bn,bt, &
                    theta,was_signed,lok)
                if(.not.lok)return
                lo_h(k,:)=hf(3,:);lo_b(k,:)=bf(3,:)
                if(was_signed)signed=signed+1
                if((.not.was_signed).and.theta<1.0_KREAL)limited=limited+1
                call reconstruct_cell_pair(h(:,:,k),b(:,:,k),i,ny,dx,dy,hf,bf,hn,ht,bn,bt, &
                    theta,was_signed,lok)
                if(.not.lok)return
                up_h(k,:)=hf(4,:);up_b(k,:)=bf(4,:)
                if(was_signed)signed=signed+1
                if((.not.was_signed).and.theta<1.0_KREAL)limited=limited+1
            enddo
            do a=1,3
                call s2w1_wall_node(-1,dt,lo_h(:,a),lo_b(:,a),u,v,weight, &
                    S2W1_wall_ux(1),S2W1_wall_theta(1),fh,fb,fw,lok)
                if(.not.lok)return
                lower_h(i,:)=lower_h(i,:)+gw(a)*fh
                lower_b(i,:)=lower_b(i,:)+gw(a)*fb
                lower_w(:,i)=lower_w(:,i)+gw(a)*fw
                call s2w1_wall_node(1,dt,up_h(:,a),up_b(:,a),u,v,weight, &
                    S2W1_wall_ux(2),S2W1_wall_theta(2),fh,fb,fw,lok)
                if(.not.lok)return
                upper_h(i,:)=upper_h(i,:)+gw(a)*fh
                upper_b(i,:)=upper_b(i,:)+gw(a)*fb
                upper_w(:,i)=upper_w(:,i)+gw(a)*fw
            enddo
        enddo
        ok=.true.
    end subroutine s2w1_compute_wall_flux

    ! The unchanged discrete instantaneous transport flux, without a face collision
    ! closure. Used by the W1D dt^2 coefficient correction at the existing midpoint.
    subroutine s2w1d_instantaneous_flux(h,b,u,v,dx,dy,fxh,fxb,fyh,fyb,limited,signed,ok)
        real(KREAL),intent(in)::h(:,:,:),b(:,:,:),u(:),v(:),dx,dy
        real(KREAL),intent(out)::fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:)
        integer,intent(out)::limited,signed
        logical,intent(out)::ok
        real(KREAL),allocatable::ch(:,:,:,:,:),cb(:,:,:,:,:)
        real(KREAL)::hf(4,3),bf(4,3),hn(4,3),ht(4,3),bn(4,3),bt(4,3),theta
        real(KREAL)::rh(size(u)),rb(size(u))
        type(FirstLayerPlane)::cacheh,cacheb
        integer::nx,ny,nvel,k,i,j,a,ir
        logical::lok,was_signed
        nx=size(h,1);ny=size(h,2);nvel=size(h,3);limited=0;signed=0;ok=.false.
        if(S2W1_wall_mode/=1)return
        allocate(ch(nx,ny,nvel,4,3),cb(nx,ny,nvel,4,3))
        do k=1,nvel
            call prepare_first_layer_plane(cacheh,h(:,:,k),lok,.true.);if(.not.lok)return
            call prepare_first_layer_plane(cacheb,b(:,:,k),lok,.true.);if(.not.lok)return
            do j=1,ny;do i=1,nx
                call reconstruct_cell_pair(h(:,:,k),b(:,:,k),i,j,dx,dy,hf,bf,hn,ht,bn,bt, &
                    theta,was_signed,lok,cacheh,cacheb)
                if(.not.lok)return
                ch(i,j,k,:,:)=hf;cb(i,j,k,:,:)=bf
                if(was_signed)signed=signed+1
                if((.not.was_signed).and.theta<1.0_KREAL)limited=limited+1
            enddo;enddo
        enddo
        fxh=0.0_KREAL;fxb=0.0_KREAL;fyh=0.0_KREAL;fyb=0.0_KREAL
        do j=1,ny;do i=1,nx
            ir=ip(i+1,nx)
            do a=1,3
                rh=merge(ch(i,j,:,2,a),ch(ir,j,:,1,a),u>=0.0_KREAL)
                rb=merge(cb(i,j,:,2,a),cb(ir,j,:,1,a),u>=0.0_KREAL)
                fxh(i,j,:)=fxh(i,j,:)+gw(a)*u*rh
                fxb(i,j,:)=fxb(i,j,:)+gw(a)*u*rb
                if(j<ny)then
                    rh=merge(ch(i,j,:,4,a),ch(i,j+1,:,3,a),v>=0.0_KREAL)
                    rb=merge(cb(i,j,:,4,a),cb(i,j+1,:,3,a),v>=0.0_KREAL)
                    fyh(i,j,:)=fyh(i,j,:)+gw(a)*v*rh
                    fyb(i,j,:)=fyb(i,j,:)+gw(a)*v*rb
                endif
            enddo
        enddo;enddo
        ok=.true.
    end subroutine s2w1d_instantaneous_flux

    ! One additional reconstruction, no equilibrium closure or accepted update.
    subroutine s2w2r_instantaneous(h,b,u,v,weight,dx,dy,fxh,fxb,fyh,fyb,lh,lb,limited,signed,ok)
        real(KREAL),intent(in)::h(:,:,:),b(:,:,:),u(:),v(:),weight(:),dx,dy
        real(KREAL),intent(out)::fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:),lh(:,:),lb(:,:)
        integer,intent(out)::limited,signed
        logical,intent(out)::ok
        real(KREAL),allocatable::ch(:,:,:,:,:),cb(:,:,:,:,:)
        real(KREAL)::hf(4,3),bf(4,3),hn(4,3),ht(4,3),bn(4,3),bt(4,3),theta
        real(KREAL)::rh(size(u)),rb(size(u)),fh(size(u)),fb(size(u)),fw(4)
        type(FirstLayerPlane)::cacheh,cacheb
        integer::nx,ny,nvel,k,i,j,a,ir
        logical::lok,was_signed
        nx=size(h,1);ny=size(h,2);nvel=size(h,3);limited=0;signed=0;ok=.false.
        allocate(ch(nx,ny,nvel,4,3),cb(nx,ny,nvel,4,3))
        do k=1,nvel
            call prepare_first_layer_plane(cacheh,h(:,:,k),lok,.true.);if(.not.lok)return
            call prepare_first_layer_plane(cacheb,b(:,:,k),lok,.true.);if(.not.lok)return
            do j=1,ny;do i=1,nx
                call reconstruct_cell_pair(h(:,:,k),b(:,:,k),i,j,dx,dy,hf,bf,hn,ht,bn,bt, &
                    theta,was_signed,lok,cacheh,cacheb)
                if(.not.lok)return
                ch(i,j,k,:,:)=hf;cb(i,j,k,:,:)=bf
                if(was_signed)signed=signed+1
                if((.not.was_signed).and.theta<1.0_KREAL)limited=limited+1
            enddo;enddo
        enddo
        fxh=0;fxb=0;fyh=0;fyb=0;lh=0;lb=0
        do j=1,ny;do i=1,nx
            ir=ip(i+1,nx)
            do a=1,3
                rh=merge(ch(i,j,:,2,a),ch(ir,j,:,1,a),u>=0.0_KREAL)
                rb=merge(cb(i,j,:,2,a),cb(ir,j,:,1,a),u>=0.0_KREAL)
                fxh(i,j,:)=fxh(i,j,:)+gw(a)*u*rh
                fxb(i,j,:)=fxb(i,j,:)+gw(a)*u*rb
                if(j<ny)then
                    rh=merge(ch(i,j,:,4,a),ch(i,j+1,:,3,a),v>=0.0_KREAL)
                    rb=merge(cb(i,j,:,4,a),cb(i,j+1,:,3,a),v>=0.0_KREAL)
                    fyh(i,j,:)=fyh(i,j,:)+gw(a)*v*rh
                    fyb(i,j,:)=fyb(i,j,:)+gw(a)*v*rb
                endif
            enddo
        enddo;enddo
        do i=1,nx;do a=1,3
            call s2w2_instantaneous_node(-1,ch(i,1,:,3,a),cb(i,1,:,3,a),u,v,weight, &
                S2W1_wall_ux(1),S2W1_wall_theta(1),fh,fb,fw,lok)
            if(.not.lok)return
            lh(i,:)=lh(i,:)+gw(a)*fh;lb(i,:)=lb(i,:)+gw(a)*fb
            call s2w2_instantaneous_node(1,ch(i,ny,:,4,a),cb(i,ny,:,4,a),u,v,weight, &
                S2W1_wall_ux(2),S2W1_wall_theta(2),fh,fb,fw,lok)
            if(.not.lok)return
            fyh(i,ny,:)=fyh(i,ny,:)+gw(a)*fh;fyb(i,ny,:)=fyb(i,ny,:)+gw(a)*fb
        enddo;enddo
        ok=.true.
    end subroutine s2w2r_instantaneous

    subroutine s2w2r_correct(w,h,b,mh,mb,u,v,weight,dx,dy,dt,fxw,fyw,fxh,fxb,fyh,fyb,lw,lh,lb, &
        x0h,x0b,y0h,y0b,xkh,xkb,ykh,ykb,l0h,l0b,lkh,lkb,diag,predict_s,instant_s,apply_s,ok)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),mh(:,:,:),mb(:,:,:),u(:),v(:),weight(:),dx,dy,dt
        real(KREAL),intent(in)::x0h(:,:,:),x0b(:,:,:),y0h(:,:,:),y0b(:,:,:),xkh(:,:,:),xkb(:,:,:),ykh(:,:,:),ykb(:,:,:)
        real(KREAL),intent(in)::l0h(:,:),l0b(:,:),lkh(:,:),lkb(:,:)
        real(KREAL),intent(inout)::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:),lw(:,:),lh(:,:),lb(:,:)
        real(KREAL),intent(out)::diag(6),predict_s,instant_s,apply_s
        logical,intent(out)::ok
        real(KREAL),allocatable::ws(:,:,:),hs(:,:,:),bs(:,:,:),xh(:,:,:),xb(:,:,:),yh(:,:,:),yb(:,:,:),bh(:,:),bb(:,:)
        real(KREAL)::dlh(size(u)),dlb(size(u)),mom(4),rho,ux,vy,theta_gas,p,t0,t1,err,oldmom(4)
        integer::nx,ny,i,j,im,il,is,du,env_status
        character(len=8)::save_env
        character(len=80)::path
        logical::lok
        nx=size(h,1);ny=size(h,2);ok=.false.;diag=0;predict_s=0;instant_s=0;apply_s=0
        allocate(ws(4,nx,ny),hs(nx,ny,size(u)),bs(nx,ny,size(u)),xh(nx,ny,size(u)),xb(nx,ny,size(u)), &
            yh(nx,ny,size(u)),yb(nx,ny,size(u)),bh(nx,size(u)),bb(nx,size(u)))
        ws=w;hs=h;bs=b
        call cpu_time(t0)
        do j=1,ny;do i=1,nx
            im=ip(i-1,nx)
            dlh=-(x0h(i,j,:)-x0h(im,j,:))/dx;dlb=-(x0b(i,j,:)-x0b(im,j,:))/dx
            if(j==1)then
                dlh=dlh-(y0h(i,j,:)-l0h(i,:))/dy;dlb=dlb-(y0b(i,j,:)-l0b(i,:))/dy
            else
                dlh=dlh-(y0h(i,j,:)-y0h(i,j-1,:))/dy;dlb=dlb-(y0b(i,j,:)-y0b(i,j-1,:))/dy
            endif
            call moment_2d(dlh,dlb,u,v,weight,mom,lok);if(.not.lok)return
            ws(:,i,j)=w(:,i,j)+0.5_KREAL*dt*mom
            hs(i,j,:)=h(i,j,:)+0.5_KREAL*dt*(dlh+(mh(i,j,:)-h(i,j,:))/tau0)
            bs(i,j,:)=b(i,j,:)+0.5_KREAL*dt*(dlb+(mb(i,j,:)-b(i,j,:))/tau0)
            call primitive(ws(:,i,j),rho,ux,vy,theta_gas,p,lok)
            if(.not.lok)exit
            call moment_2d(hs(i,j,:),bs(i,j,:),u,v,weight,mom,lok);if(.not.lok)exit
            diag(1)=max(diag(1),maxval(abs(mom-ws(:,i,j)))/max(1.0_KREAL,maxval(abs(ws(:,i,j)))))
        enddo
            if(.not.lok)exit
        enddo
        diag(2)=minval(hs);diag(3)=minval(bs)
        if(.not.lok.or.diag(1)>gate.or..not.all(ieee_is_finite(hs)).or..not.all(ieee_is_finite(bs)))then
            call write_state('w2r_predictor_rejected.bin',ws,hs,bs,u,v,weight,dt/2)
            call write_state('w2r_last_accepted.bin',w,h,b,u,v,weight,0.0_KREAL)
            return
        endif
        call cpu_time(t1);predict_s=t1-t0
        call cpu_time(t0)
        call s2w2r_instantaneous(hs,bs,u,v,weight,dx,dy,xh,xb,yh,yb,bh,bb,il,is,lok)
        diag(4)=real(il,KREAL);diag(5)=real(is,KREAL)
        if(.not.lok)then
            call write_state('w2r_predictor_rejected.bin',ws,hs,bs,u,v,weight,dt/2)
            call write_state('w2r_last_accepted.bin',w,h,b,u,v,weight,0.0_KREAL)
            return
        endif
        call cpu_time(t1);instant_s=t1-t0
        call get_environment_variable('S2W2R_SAVE_PREDICTOR',save_env,status=env_status)
        if(env_status==0.and.trim(save_env)=='1')then
            w2r_research_calls=w2r_research_calls+1
            write(path,'(A,I3.3,A)')'w2r_predictor_step',w2r_research_calls,'.bin'
            call write_state(trim(path),ws,hs,bs,u,v,weight,dt/2)
            write(path,'(A,I3.3,A)')'w2r_coefficients_step',w2r_research_calls,'.bin'
            open(newunit=du,file=trim(path),status='replace',access='stream',form='unformatted')
            write(du)int(nx,KINT),int(ny,KINT),int(size(u),KINT),dt
            write(du)fxh,fxb,fyh,fyb,lh,lb,x0h,x0b,y0h,y0b,l0h,l0b, &
                xkh,xkb,ykh,ykb,lkh,lkb,xh,xb,yh,yb,bh,bb
            close(du)
        endif
        call cpu_time(t0)
        ! Arrays xh/etc become proposed shared integrals before any official flux is replaced.
        xh=fxh+dt*(xh-x0h)-0.5_KREAL*dt*dt*xkh
        xb=fxb+dt*(xb-x0b)-0.5_KREAL*dt*dt*xkb
        yh=fyh+dt*(yh-y0h)-0.5_KREAL*dt*dt*ykh
        yb=fyb+dt*(yb-y0b)-0.5_KREAL*dt*dt*ykb
        bh=lh+dt*(bh-l0h)-0.5_KREAL*dt*dt*lkh
        bb=lb+dt*(bb-l0b)-0.5_KREAL*dt*dt*lkb
        do i=1,nx
            call moment_2d(bh(i,:),bb(i,:),u,v,weight,mom,lok);if(.not.lok)return
            diag(6)=max(diag(6),abs(mom(1)))
            call moment_2d(yh(i,ny,:),yb(i,ny,:),u,v,weight,mom,lok);if(.not.lok)return
            diag(6)=max(diag(6),abs(mom(1)))
        enddo
        if(diag(6)>gate)return
        if(.not.all(ieee_is_finite(xh)).or..not.all(ieee_is_finite(xb)).or. &
            .not.all(ieee_is_finite(yh)).or..not.all(ieee_is_finite(yb)))return
        fxh=xh;fxb=xb;fyh=yh;fyb=yb;lh=bh;lb=bb
        do j=1,ny;do i=1,nx
            call moment_2d(fxh(i,j,:),fxb(i,j,:),u,v,weight,fxw(:,i,j),lok);if(.not.lok)return
            call moment_2d(fyh(i,j,:),fyb(i,j,:),u,v,weight,fyw(:,i,j),lok);if(.not.lok)return
        enddo;enddo
        do i=1,nx
            call moment_2d(lh(i,:),lb(i,:),u,v,weight,lw(:,i),lok);if(.not.lok)return
        enddo
        call cpu_time(t1);apply_s=t1-t0;ok=.true.
    end subroutine s2w2r_correct

    subroutine build_targets(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,target_error,ok)
        real(KREAL),intent(in)::w(:,:,:),u(:),v(:),weight(:)
        real(KREAL),intent(out)::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:)
        real(KREAL),intent(out)::eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
        real(KREAL),intent(out)::target_error
        logical,intent(out)::ok
        integer::nx,ny,i,j,a,c,m
        real(KREAL)::nodes(4,3,3),avg(4),wq(4),tmh(size(u)),tmb(size(u))
        real(KREAL)::tkxh(size(u)),tkxb(size(u)),tkyh(size(u)),tkyb(size(u))
        real(KREAL)::teah(size(u)),teab(size(u)),kap(2,2),err,wt2
        logical::lok
        nx=size(w,2);ny=size(w,3);mh=0.0_KREAL;mb=0.0_KREAL;kxh=0.0_KREAL;kxb=0.0_KREAL
        kyh=0.0_KREAL;kyb=0.0_KREAL;eah=0.0_KREAL;eab=0.0_KREAL;kappa=0.0_KREAL
        target_error=0.0_KREAL
        do j=active_lo(),active_hi(ny)
            do i=xlo(),xhi(nx)
                failure_i=i;failure_j=j
                if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
                    call old_target_at(w(:,i,j),[0.0_KREAL,0.0_KREAL],u,v,tmh,tmb,teah,teab,lok)
                    if(.not.lok)then;ok=.false.;return;endif
                    mh(i,j,:)=tmh;mb(i,j,:)=tmb;eah(i,j,:)=teah;eab(i,j,:)=teab
                    kxh(i,j,:)=0.0_KREAL;kxb(i,j,:)=0.0_KREAL
                    kyh(i,j,:)=0.0_KREAL;kyb(i,j,:)=0.0_KREAL;kappa(:,:,i,j)=0.0_KREAL
                    call moment_2d(tmh,tmb,u,v,weight,wq,lok);if(.not.lok)then;ok=.false.;return;endif
                    err=maxval(abs(wq-w(:,i,j)))/max(1.0_KREAL,maxval(abs(w(:,i,j))))
                    target_error=max(target_error,err)
                else
                    do m=1,4
                        call s2d1_cell9(w(m,:,:),i,j,nodes(m,:,:),lok,S2W1_wall_mode>=1)
                        if(.not.lok)then;ok=.false.;return;endif
                    enddo
                    avg=0.0_KREAL
                    do c=1,3;do a=1,3
                        avg=avg+gw(a)*gw(c)*nodes(:,a,c)
                    enddo;enddo
                    do c=1,3;do a=1,3
                        wq=nodes(:,a,c)+w(:,i,j)-avg;wt2=gw(a)*gw(c)
                        call make_mke(wq,u,v,weight,tmh,tmb,tkxh,tkxb,tkyh,tkyb,teah,teab,kap,lok,err)
                        if(.not.lok)then;ok=.false.;return;endif
                        mh(i,j,:)=mh(i,j,:)+wt2*tmh;mb(i,j,:)=mb(i,j,:)+wt2*tmb
                        kxh(i,j,:)=kxh(i,j,:)+wt2*tkxh;kxb(i,j,:)=kxb(i,j,:)+wt2*tkxb
                        kyh(i,j,:)=kyh(i,j,:)+wt2*tkyh;kyb(i,j,:)=kyb(i,j,:)+wt2*tkyb
                        target_error=max(target_error,err)
                    enddo;enddo
                    call heat(w(:,i,j),kxh(i,j,:),kxb(i,j,:),u,v,weight,kappa(:,1,i,j),lok)
                    if(.not.lok)then;ok=.false.;return;endif
                    call heat(w(:,i,j),kyh(i,j,:),kyb(i,j,:),u,v,weight,kappa(:,2,i,j),lok)
                    if(.not.lok)then;ok=.false.;return;endif
                    call analytic_maxwellian_2d(w(:,i,j),ck,u,v,eah(i,j,:),eab(i,j,:),lok)
                    if(.not.lok)then;ok=.false.;return;endif
                endif
            enddo
        enddo
        ok=.true.
    end subroutine build_targets

    subroutine advance_step(w,h,b,u,v,weight,dx,dy,dt,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa, &
        diag,flux_seconds,collision_seconds,ok,stream_w,stream_h,stream_b,native_face_error, &
        native_limited,native_signed)
        real(KREAL),intent(inout)::w(:,:,:),h(:,:,:),b(:,:,:)
        real(KREAL),intent(in)::u(:),v(:),weight(:),dx,dy,dt
        real(KREAL),intent(inout)::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:)
        real(KREAL),intent(inout)::kyh(:,:,:),kyb(:,:,:),eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
        real(KREAL),intent(out)::diag(12),flux_seconds,collision_seconds
        logical,intent(out)::ok
        real(KREAL),intent(in),optional::stream_w(:,:,:),stream_h(:,:,:),stream_b(:,:,:)
        real(KREAL),intent(in),optional::native_face_error
        integer,intent(in),optional::native_limited,native_signed
        integer::nx,ny,nvel,i,j,k,im,jm,limited,signed
        real(KREAL),allocatable::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:)
        real(KREAL),allocatable::wn(:,:,:),hs(:,:,:),bs(:,:,:),hn(:,:,:),bn(:,:,:)
        real(KREAL),allocatable::mhn(:,:,:),mbn(:,:,:),kxhn(:,:,:),kxbn(:,:,:),kyhn(:,:,:),kybn(:,:,:)
        real(KREAL),allocatable::eahn(:,:,:),eabn(:,:,:),kapn(:,:,:,:)
        real(KREAL)::face_err,target_err,beta,beta_old,rhs(2),rnew(2),rold(2),rhat(2),qtmp(2),a22(2,2),det
        real(KREAL)::rh(size(u)),rb(size(u)),eqh(size(u)),eqb(size(u)),mom(4),scale,cnres,mmerr
        real(KREAL)::t0,t1
        logical::lok
        nx=size(h,1);ny=size(h,2);nvel=size(h,3);diag=0.0_KREAL;ok=.false.
        allocate(fxw(4,nx,ny),fyw(4,nx,ny),fxh(nx,ny,nvel),fxb(nx,ny,nvel), &
            fyh(nx,ny,nvel),fyb(nx,ny,nvel),wn(4,nx,ny),hs(nx,ny,nvel),bs(nx,ny,nvel), &
            hn(nx,ny,nvel),bn(nx,ny,nvel),mhn(nx,ny,nvel),mbn(nx,ny,nvel), &
            kxhn(nx,ny,nvel),kxbn(nx,ny,nvel),kyhn(nx,ny,nvel),kybn(nx,ny,nvel), &
            eahn(nx,ny,nvel),eabn(nx,ny,nvel),kapn(2,2,nx,ny))
        if(present(stream_w).or.present(stream_h).or.present(stream_b).or.present(native_face_error).or. &
            present(native_limited).or.present(native_signed))then
            if(.not.(present(stream_w).and.present(stream_h).and.present(stream_b).and. &
                present(native_face_error).and.present(native_limited).and.present(native_signed)))return
            if(any(shape(stream_w)/=shape(wn)).or.any(shape(stream_h)/=shape(hs)).or. &
                any(shape(stream_b)/=shape(bs)))return
            wn=stream_w;hs=stream_h;bs=stream_b;face_err=native_face_error
            limited=native_limited;signed=native_signed;flux_seconds=0.0_KREAL
        else
            call cpu_time(t0)
            call compute_flux(w,h,b,u,v,weight,dx,dy,dt,fxw,fyw,fxh,fxb,fyh,fyb, &
                face_err,limited,signed,lok)
            call cpu_time(t1);flux_seconds=t1-t0
            if(.not.lok)return
            do j=active_lo(),active_hi(ny);jm=ip(j-1,ny)
                do i=xlo(),xhi(nx)
                failure_i=i;failure_j=j;im=ip(i-1,nx)
                    wn(:,i,j)=w(:,i,j)-(fxw(:,i,j)-fxw(:,im,j))/dx-(fyw(:,i,j)-fyw(:,i,jm))/dy
                    hs(i,j,:)=h(i,j,:)-(fxh(i,j,:)-fxh(im,j,:))/dx-(fyh(i,j,:)-fyh(i,jm,:))/dy
                    bs(i,j,:)=b(i,j,:)-(fxb(i,j,:)-fxb(im,j,:))/dx-(fyb(i,j,:)-fyb(i,jm,:))/dy
                enddo
            enddo
        endif
        call exchange_w(wn)
        hn=h;bn=b
        call cpu_time(t0)
        if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
            mhn=0.0_KREAL;mbn=0.0_KREAL;kxhn=0.0_KREAL;kxbn=0.0_KREAL
            kyhn=0.0_KREAL;kybn=0.0_KREAL;eahn=0.0_KREAL;eabn=0.0_KREAL
            kapn=0.0_KREAL;target_err=0.0_KREAL
        else
            call build_targets(wn,u,v,weight,mhn,mbn,kxhn,kxbn,kyhn,kybn,eahn,eabn,kapn,target_err,lok)
            if(.not.lok)return
        endif
        beta=dt/(2.0_KREAL*tau0);cnres=0.0_KREAL;mmerr=0.0_KREAL
        diag(9)=huge(1.0_KREAL);diag(10)=0.0_KREAL;diag(11)=huge(1.0_KREAL);diag(12)=0.0_KREAL
        do j=active_lo(),active_hi(ny)
            do i=xlo(),xhi(nx)
                failure_i=i;failure_j=j
                beta_old=beta
                if(S2PhysicalGas==1)then
                    beta_old=dt/(2.0_KREAL*s2_tau(w(:,i,j)))
                    beta=dt/(2.0_KREAL*s2_tau(wn(:,i,j)))
                endif
                diag(7)=max(diag(7),2.0_KREAL*beta_old,2.0_KREAL*beta)
                if(trim(adjustl(TimeIntegration_method))=='S2D1VL')then
                    call heat(w(:,i,j),h(i,j,:),b(i,j,:),u,v,weight,rold,lok);if(.not.lok)return
                    call old_target_at(w(:,i,j),rold,u,v,rh,rb,eahn(i,j,:),eabn(i,j,:),lok)
                    if(.not.lok)return
                    call old_target_at(wn(:,i,j),rold,u,v,eqh,eqb,mhn(i,j,:),mbn(i,j,:),lok)
                    if(.not.lok)return
                    eahn(i,j,:)=mhn(i,j,:);eabn(i,j,:)=mbn(i,j,:)
                    hn(i,j,:)=(hs(i,j,:)+beta*(eqh+rh-h(i,j,:)))/(1.0_KREAL+beta)
                    bn(i,j,:)=(bs(i,j,:)+beta*(eqb+rb-b(i,j,:)))/(1.0_KREAL+beta)
                    scale=max(1.0_KREAL,maxval(abs(hs(i,j,:))),maxval(abs(bs(i,j,:))), &
                        maxval(abs(hn(i,j,:))),maxval(abs(bn(i,j,:))))
                    cnres=max(cnres,max(maxval(abs((1.0_KREAL+beta)*hn(i,j,:)-hs(i,j,:)- &
                        beta*(eqh+rh-h(i,j,:)))),maxval(abs((1.0_KREAL+beta)*bn(i,j,:)-bs(i,j,:)- &
                        beta*(eqb+rb-b(i,j,:)))))/scale)
                    call moment_2d(eqh,eqb,u,v,weight,mom,lok);if(.not.lok)return
                    target_err=max(target_err,maxval(abs(mom-wn(:,i,j)))/max(1.0_KREAL,maxval(abs(wn(:,i,j)))))
                else
                    call heat(w(:,i,j),h(i,j,:)-mh(i,j,:),b(i,j,:)-mb(i,j,:),u,v,weight,rold,lok)
                    if(.not.lok)return
                    rh=hs(i,j,:)+beta_old*(mh(i,j,:)+kxh(i,j,:)*rold(1)+kyh(i,j,:)*rold(2)-h(i,j,:))
                    rb=bs(i,j,:)+beta_old*(mb(i,j,:)+kxb(i,j,:)*rold(1)+kyb(i,j,:)*rold(2)-b(i,j,:))
                    call heat(wn(:,i,j),rh,rb,u,v,weight,rhs,lok);if(.not.lok)return
                    call heat(wn(:,i,j),mhn(i,j,:),mbn(i,j,:),u,v,weight,qtmp,lok);if(.not.lok)return
                    rhs=rhs-qtmp
                    a22=-beta*kapn(:,:,i,j)
                    a22(1,1)=a22(1,1)+1.0_KREAL+beta
                    a22(2,2)=a22(2,2)+1.0_KREAL+beta
                    det=a22(1,1)*a22(2,2)-a22(1,2)*a22(2,1)
                    if(.not.ieee_is_finite(det).or.det<=0.0_KREAL)return
                    rnew(1)=(rhs(1)*a22(2,2)-a22(1,2)*rhs(2))/det
                    rnew(2)=(a22(1,1)*rhs(2)-rhs(1)*a22(2,1))/det
                    eqh=mhn(i,j,:)+kxhn(i,j,:)*rnew(1)+kyhn(i,j,:)*rnew(2)
                    eqb=mbn(i,j,:)+kxbn(i,j,:)*rnew(1)+kybn(i,j,:)*rnew(2)
                    hn(i,j,:)=(rh+beta*eqh)/(1.0_KREAL+beta)
                    bn(i,j,:)=(rb+beta*eqb)/(1.0_KREAL+beta)
                endif
                if(trim(adjustl(TimeIntegration_method))/='S2D1VL')then
                    scale=max(1.0_KREAL,maxval(abs(rh)),maxval(abs(rb)),maxval(abs(hn(i,j,:))), &
                        maxval(abs(bn(i,j,:))))
                    cnres=max(cnres,max(maxval(abs((1+beta)*hn(i,j,:)-rh-beta*eqh)), &
                        maxval(abs((1+beta)*bn(i,j,:)-rb-beta*eqb)))/scale)
                    call heat(wn(:,i,j),hn(i,j,:)-mhn(i,j,:),bn(i,j,:)-mbn(i,j,:),u,v,weight,rhat,lok)
                    if(.not.lok)return
                    cnres=max(cnres,maxval(abs(rhat-rnew))/max(1.0_KREAL,maxval(abs(rhat)),maxval(abs(rnew))))
                endif
                call moment_2d(hn(i,j,:),bn(i,j,:),u,v,weight,mom,lok);if(.not.lok)return
                mmerr=max(mmerr,maxval(abs(mom-wn(:,i,j)))/max(1.0_KREAL,maxval(abs(wn(:,i,j)))))
                diag(9)=min(diag(9),minval(hn(i,j,:)));diag(10)=max(diag(10),maxval(hn(i,j,:)))
                diag(11)=min(diag(11),minval(bn(i,j,:)));diag(12)=max(diag(12),maxval(bn(i,j,:)))
            enddo
        enddo
        call cpu_time(t1);collision_seconds=t1-t0
        diag(1)=face_err;diag(2)=target_err;diag(3)=cnres;diag(4)=mmerr
        diag(5)=real(limited,KREAL);diag(6)=real(signed,KREAL)
        diag(8)=dt*(maxval(abs(u))/dx+maxval(abs(v))/dy)
        if(.not.all(ieee_is_finite(wn)).or..not.all(ieee_is_finite(hn)).or..not.all(ieee_is_finite(bn)))return
        if(cnres>gate.or.mmerr>gate)return
        w=wn;h=hn;b=bn;mh=mhn;mb=mbn;kxh=kxhn;kxb=kxbn;kyh=kyhn;kyb=kybn
        eah=eahn;eab=eabn;kappa=kapn;ok=.true.
        deallocate(fxw,fyw,fxh,fxb,fyh,fyb,wn,hs,bs,hn,bn,mhn,mbn,kxhn,kxbn,kyhn,kybn, &
            eahn,eabn,kapn)
    end subroutine advance_step

    subroutine init_velocity(u,v,weight,raw_node,raw_weight)
        real(KREAL),intent(out)::u(:),v(:),weight(:),raw_node(:),raw_weight(:)
        real(KREAL),parameter::x(16)=[-4.68873893930581875_KREAL,-3.86944790486012291_KREAL, &
            -3.17699916197995602_KREAL,-2.54620215784748138_KREAL,-1.95178799091625388_KREAL, &
            -1.38025853919888086_KREAL,-0.82295144914465590_KREAL,-0.27348104613815250_KREAL, &
             0.27348104613815250_KREAL,0.82295144914465590_KREAL,1.38025853919888086_KREAL, &
             1.95178799091625388_KREAL,2.54620215784748138_KREAL,3.17699916197995602_KREAL, &
             3.86944790486012291_KREAL,4.68873893930581875_KREAL]
        real(KREAL),parameter::wr(16)=[2.65480747401116728e-10_KREAL,2.32098084486520325e-7_KREAL, &
            2.71186009253787598e-5_KREAL,9.32284008624180709e-4_KREAL,1.28803115355099890e-2_KREAL, &
            8.38100413989858323e-2_KREAL,2.80647458528533678e-1_KREAL,5.07929479016613561e-1_KREAL, &
            5.07929479016613561e-1_KREAL,2.80647458528533678e-1_KREAL,8.38100413989858323e-2_KREAL, &
            1.28803115355099890e-2_KREAL,9.32284008624180709e-4_KREAL,2.71186009253787598e-5_KREAL, &
            2.32098084486520325e-7_KREAL,2.65480747401116728e-10_KREAL]
        real(KREAL)::ue(16),we(16)
        integer::i,j,k
        if(VelocityMode()/='full16')then
            call FillVelocity(u,v,weight,raw_node,raw_weight);return
        endif
        raw_node=x;raw_weight=wr;ue=sqrt(2.0_KREAL)*x;we=sqrt(2.0_KREAL)*exp(x*x)*wr
        do j=1,16;do i=1,16;k=i+(j-1)*16
            u(k)=ue(i);v(k)=ue(j);weight(k)=we(i)*we(j)
        enddo;enddo
    end subroutine init_velocity

    subroutine state_totals(w,total)
        real(KREAL),intent(in)::w(:,:,:)
        real(KREAL),intent(out)::total(4)
        integer::m
        do m=1,4;total(m)=sum(w(m,:,:));enddo
    end subroutine state_totals

    subroutine write_state(path,w,h,b,u,v,weight,time)
        character(len=*),intent(in)::path
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),weight(:),time
        integer::unit
        open(newunit=unit,file=path,status='replace',access='stream',form='unformatted')
        write(unit)int(size(h,1),KINT),int(size(h,2),KINT),int(size(h,3),KINT),time
        write(unit)w;write(unit)h;write(unit)b;write(unit)u;write(unit)v;write(unit)weight
        close(unit)
    end subroutine write_state

    subroutine write_field(path,w,h,b,eah,eab,u,v,weight,dx,dy,time)
        character(len=*),intent(in)::path
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),eah(:,:,:),eab(:,:,:)
        real(KREAL),intent(in)::u(:),v(:),weight(:),dx,dy,time
        integer::unit,i,j
        real(KREAL)::rho,ux,vy,theta_gas,p,q(2),x,y
        logical::ok
        open(newunit=unit,file=path,status='replace',action='write')
        write(unit,'(A)')'x,y,rho,p,qx,qy,U,V,time'
        do j=1,size(h,2);y=(real(j,KREAL)-0.5_KREAL)*dy
            do i=1,size(h,1);x=(real(i,KREAL)-0.5_KREAL)*dx
                call primitive(w(:,i,j),rho,ux,vy,theta_gas,p,ok)
                call heat(w(:,i,j),h(i,j,:)-eah(i,j,:),b(i,j,:)-eab(i,j,:),u,v,weight,q,ok)
                write(unit,'(9(ES24.16,:,","))')x,y,rho,p,q(1),q(2),ux,vy,time
            enddo
        enddo
        close(unit)
    end subroutine write_field

    subroutine write_velocity(u,v,weight,raw_node,raw_weight)
        real(KREAL),intent(in)::u(:),v(:),weight(:),raw_node(16),raw_weight(16)
        integer::unit,k
        open(newunit=unit,file='velocity.csv',status='replace',action='write')
        write(unit,'(A)')'k,u,v,effective_weight'
        do k=1,size(u);write(unit,'(I0,3(",",ES24.16))')k,u(k),v(k),weight(k);enddo
        close(unit)
        open(newunit=unit,file='gh16_1d.csv',status='replace',action='write')
        write(unit,'(A)')'i,raw_node,raw_weight,u,effective_weight'
        do k=1,16
            write(unit,'(I0,4(",",ES24.16))')k,raw_node(k),raw_weight(k),sqrt(2.0_KREAL)*raw_node(k), &
                sqrt(2.0_KREAL)*exp(raw_node(k)**2)*raw_weight(k)
        enddo
        close(unit)
    end subroutine write_velocity

    subroutine write_quadrature_check(raw_node,raw_weight,u,v,weight)
        real(KREAL),intent(in)::raw_node(16),raw_weight(16),u(:),v(:),weight(:)
        real(KREAL)::exact(0:6),calc(0:6),pi0
        integer::n,unit
        pi0=sqrt(pi_d);exact=0.0_KREAL
        exact(0)=pi0;exact(2)=0.5_KREAL*pi0;exact(4)=0.75_KREAL*pi0
        exact(6)=1.875_KREAL*pi0
        do n=0,6;calc(n)=sum(raw_weight*raw_node**n);enddo
        open(newunit=unit,file='quadrature_checks.csv',status='replace',action='write')
        write(unit,'(A)')'moment,calculated,exact,abs_error'
        do n=0,6
            write(unit,'(I0,3(",",ES24.16))')n,calc(n),exact(n),abs(calc(n)-exact(n))
        enddo
        write(unit,'(A,ES24.16)')'min_raw_weight,',minval(raw_weight)
        write(unit,'(A,ES24.16)')'min_effective_weight,',minval(weight)
        write(unit,'(A,ES24.16)')'max_abs_u,',maxval(abs(u))
        write(unit,'(A,ES24.16)')'max_abs_v,',maxval(abs(v))
        close(unit)
    end subroutine write_quadrature_check

    subroutine write_reconstruction_check()
        integer,parameter::grids(3)=[16,32,64]
        real(KREAL),parameter::gn(3)=[-0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL),0.0_KREAL, &
            0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL)]
        real(KREAL),allocatable::f(:,:)
        real(KREAL)::vl(3),vr(3),nl(3),nr(3),tl(3),tr(3),dx,xface,yp,ex
        real(KREAL)::vl2(3),vr2(3),nl2(3),nr2(3),tl2(3),tr2(3)
        real(KREAL)::ev,edn,edt,constant_error,transpose_error
        logical::ok
        integer::g,n,i,j,a,unit
        open(newunit=unit,file='reconstruction_checks.csv',status='replace',action='write')
        write(unit,'(A)')'N,value_rms,normal_derivative_rms,tangent_derivative_rms,constant_max,transpose_max'
        do g=1,3;n=grids(g);dx=2.0_KREAL*pi_d/real(n,KREAL);allocate(f(n,n))
            do j=1,n;do i=1,n
                f(i,j)=(sin(0.5_KREAL*dx)/(0.5_KREAL*dx))**2* &
                    sin((real(i,KREAL)-0.5_KREAL)*dx+(real(j,KREAL)-0.5_KREAL)*dx)
            enddo;enddo
            ev=0.0_KREAL;edn=0.0_KREAL;edt=0.0_KREAL;constant_error=0.0_KREAL
            transpose_error=0.0_KREAL
            do j=1,n;do i=1,n
                call s2d1_face3(f,1,i,j,vl,vr,nl,nr,tl,tr,ok)
                if(.not.ok)call abort_s2('manufactured WENO face reconstruction failed',313)
                xface=real(i,KREAL)*dx
                do a=1,3
                    yp=(real(j,KREAL)-0.5_KREAL+gn(a))*dx;ex=sin(xface+yp)
                    ev=ev+(vl(a)-ex)**2+(vr(a)-ex)**2
                    edn=edn+(nl(a)/dx-cos(xface+yp))**2+(nr(a)/dx-cos(xface+yp))**2
                    edt=edt+(tl(a)/dx-cos(xface+yp))**2+(tr(a)/dx-cos(xface+yp))**2
                enddo
                call s2d1_face3(transpose(f),2,i,j,vl2,vr2,nl2,nr2,tl2,tr2,ok)
                if(.not.ok)call abort_s2('transpose reconstruction failed',314)
                transpose_error=max(transpose_error,maxval(abs(vl-vl2)),maxval(abs(vr-vr2)), &
                    maxval(abs(nl-nl2)),maxval(abs(nr-nr2)),maxval(abs(tl-tl2)),maxval(abs(tr-tr2)))
            enddo;enddo
            f=1.0_KREAL;call s2d1_face3(f,1,n,n,vl,vr,nl,nr,tl,tr,ok)
            constant_error=max(maxval(abs(vl-1.0_KREAL)),maxval(abs(vr-1.0_KREAL)), &
                maxval(abs(nl)),maxval(abs(nr)),maxval(abs(tl)),maxval(abs(tr)))
            ev=sqrt(ev/real(6*n*n,KREAL));edn=sqrt(edn/real(6*n*n,KREAL))
            edt=sqrt(edt/real(6*n*n,KREAL))
            write(unit,'(I0,5(",",ES24.16))')n,ev,edn,edt,constant_error,transpose_error
            deallocate(f)
        enddo
        close(unit)
    end subroutine write_reconstruction_check

    subroutine write_local_qualification(u,v,weight)
        real(KREAL),intent(in)::u(:),v(:),weight(:)
        real(KREAL),parameter::prs(2)=[1.0_KREAL,2.0_KREAL/3.0_KREAL]
        real(KREAL),parameter::zs(4)=[0.1_KREAL,1.0_KREAL,32.0_KREAL,128.0_KREAL]
        real(KREAL)::w(4),mh(size(u)),mb(size(u)),kxh(size(u)),kxb(size(u)),kyh(size(u)),kyb(size(u))
        real(KREAL)::eah(size(u)),eab(size(u)),kap(2,2),rh(size(u)),rb(size(u)),ph(size(u)),pb(size(u))
        real(KREAL)::rseed(2),rhs(2),rnew(2),rhat(2),a22(2,2),det,beta,res,eqh(size(u)),eqb(size(u))
        real(KREAL)::rawh(size(u)),rawb(size(u)),mom(4),saved_pr,err
        logical::ok
        integer::ip0,iz,k,unit
        saved_pr=Pr;w=[1.03_KREAL,0.41_KREAL,-0.17_KREAL,1.73_KREAL]
        open(newunit=unit,file='local_cn_qualification.csv',status='replace',action='write')
        write(unit,'(A)')'Pr,z,k11,k12,k21,k22,equation_residual,r_identity,moment_residual'
        do ip0=1,2;Pr=prs(ip0)
            call make_mke(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kap,ok,err)
            if(.not.ok)call abort_s2('local qualification target failed',315)
            do k=1,size(u)
                rawh(k)=1.0e-4_KREAL*mh(k)*(u(k)-w(2)/w(1))*(v(k)-w(3)/w(1))
                rawb(k)=1.0e-4_KREAL*mb(k)*(u(k)-w(2)/w(1))*(v(k)-w(3)/w(1))
            enddo
            call projector_2d(w,rawh,rawb,1,u,v,weight,ph,pb,ok)
            if(.not.ok)call abort_s2('local qualification projection failed',316)
            rseed=[0.013_KREAL,-0.021_KREAL]
            rh=mh+kxh*rseed(1)+kyh*rseed(2)+ph
            rb=mb+kxb*rseed(1)+kyb*rseed(2)+pb
            do iz=1,4;beta=0.5_KREAL*zs(iz)
                call heat(w,rh,rb,u,v,weight,rhs,ok);if(.not.ok)call abort_s2('local heat failed',317)
                call heat(w,mh,mb,u,v,weight,rhat,ok);rhs=rhs-rhat
                a22=-beta*kap;a22(1,1)=a22(1,1)+1+beta;a22(2,2)=a22(2,2)+1+beta
                det=a22(1,1)*a22(2,2)-a22(1,2)*a22(2,1)
                rnew(1)=(rhs(1)*a22(2,2)-a22(1,2)*rhs(2))/det
                rnew(2)=(a22(1,1)*rhs(2)-rhs(1)*a22(2,1))/det
                eqh=mh+kxh*rnew(1)+kyh*rnew(2);eqb=mb+kxb*rnew(1)+kyb*rnew(2)
                ph=(rh+beta*eqh)/(1+beta);pb=(rb+beta*eqb)/(1+beta)
                res=max(maxval(abs((1+beta)*ph-rh-beta*eqh)), &
                    maxval(abs((1+beta)*pb-rb-beta*eqb)))
                call heat(w,ph-mh,pb-mb,u,v,weight,rhat,ok)
                call moment_2d(ph,pb,u,v,weight,mom,ok)
                write(unit,'(2(ES24.16,","),4(ES24.16,","),3(ES24.16,:,","))') &
                    Pr,zs(iz),kap(1,1),kap(1,2),kap(2,1),kap(2,2),res,maxval(abs(rhat-rnew)), &
                    maxval(abs(mom-w))
                if(max(res,maxval(abs(rhat-rnew)),maxval(abs(mom-w)))>gate) &
                    call abort_s2('local CN qualification exceeded 2e-12',318)
            enddo
        enddo
        close(unit);Pr=saved_pr
    end subroutine write_local_qualification

    subroutine write_equilibrium_direction_check(u,v,weight)
        real(KREAL),intent(in)::u(:),v(:),weight(:)
        real(KREAL),parameter::epsv(5)=[1.0e-3_KREAL,1.0e-4_KREAL,1.0e-5_KREAL,1.0e-6_KREAL,1.0e-7_KREAL]
        real(KREAL)::w(4),dw(4),ah(size(u)),ab(size(u)),ph(size(u)),pb(size(u)),mh(size(u)),mb(size(u))
        real(KREAL)::theta(4),jac(4,4),mom(4),fd_error,moment_error,zero_error,res,ep
        logical::ok
        integer::n,unit
        w=[1.03_KREAL,0.41_KREAL,-0.17_KREAL,1.73_KREAL]
        dw=[0.019_KREAL,-0.031_KREAL,0.023_KREAL,0.047_KREAL]
        call m_direction(w,dw,u,v,weight,ah,ab,ok,res)
        if(.not.ok)call abort_s2('analytic discrete-M direction failed',320)
        call moment_2d(ah,ab,u,v,weight,mom,ok);moment_error=maxval(abs(mom-dw))
        call m_direction(w,0.0_KREAL*dw,u,v,weight,ph,pb,ok,res)
        zero_error=max(maxval(abs(ph)),maxval(abs(pb)))
        open(newunit=unit,file='equilibrium_direction_checks.csv',status='replace',action='write')
        write(unit,'(A)')'epsilon,distribution_fd_max,moment_residual,zero_direction_max'
        do n=1,size(epsv);ep=epsv(n)
            call discrete_equilibrium_q_2d(w+ep*dw,1,u,v,weight,ph,pb,theta,jac,ok)
            if(.not.ok)call abort_s2('positive FD equilibrium failed',321)
            call discrete_equilibrium_q_2d(w-ep*dw,1,u,v,weight,mh,mb,theta,jac,ok)
            if(.not.ok)call abort_s2('negative FD equilibrium failed',322)
            fd_error=max(maxval(abs((ph-mh)/(2*ep)-ah)),maxval(abs((pb-mb)/(2*ep)-ab)))
            write(unit,'(4(ES24.16,:,","))')ep,fd_error,moment_error,zero_error
        enddo
        close(unit)
        if(moment_error>gate.or.zero_error>gate)call abort_s2('discrete-M direction identity failed',323)
    end subroutine write_equilibrium_direction_check

    subroutine write_face_wiring_check(u,v,weight)
        real(KREAL),intent(in)::u(:),v(:),weight(:)
        real(KREAL)::w(4),mh(size(u)),mb(size(u)),kx(size(u)),kxb(size(u)),ky(size(u)),kyb(size(u))
        real(KREAL)::ea(size(u)),eab(size(u)),kap(2,2),z(size(u)),dn1(4),dn2(4),dtan(4)
        real(KREAL)::fh0(size(u)),fb0(size(u)),fw0(4),fh1(size(u)),fb1(size(u)),fw1(4)
        real(KREAL)::fh2(size(u)),fb2(size(u)),fw2(4),err,res,response,linearity
        logical::ok
        integer::axis,unit
        w=[1.03_KREAL,0.41_KREAL,-0.17_KREAL,1.73_KREAL]
        dn1=[0.019_KREAL,-0.031_KREAL,0.023_KREAL,0.047_KREAL]
        dn2=[-0.013_KREAL,0.017_KREAL,-0.029_KREAL,0.037_KREAL]
        dtan=[0.011_KREAL,0.021_KREAL,-0.014_KREAL,-0.033_KREAL]
        call make_mke(w,u,v,weight,mh,mb,kx,kxb,ky,kyb,ea,eab,kap,ok,res)
        if(.not.ok)call abort_s2('face wiring equilibrium failed',324)
        z=0.0_KREAL
        open(newunit=unit,file='face_wiring_checks.csv',status='replace',action='write')
        write(unit,'(A)')'axis,response_max,linearity_max,face_moment_max'
        do axis=1,2
            call ugks_node(axis,0.005_KREAL,mh,mh,mb,mb,z,z,z,z,z,z,z,z, &
                0.0_KREAL*dn1,0.0_KREAL*dn2,0.0_KREAL*dtan,u,v,weight,fh0,fb0,fw0,err,ok)
            if(.not.ok)call abort_s2('zero-direction face wiring failed',325)
            call ugks_node(axis,0.005_KREAL,mh,mh,mb,mb,z,z,z,z,z,z,z,z, &
                dn1,dn2,dtan,u,v,weight,fh1,fb1,fw1,err,ok)
            if(.not.ok)call abort_s2('macro-direction face wiring failed',326)
            call ugks_node(axis,0.005_KREAL,mh,mh,mb,mb,z,z,z,z,z,z,z,z, &
                2*dn1,2*dn2,2*dtan,u,v,weight,fh2,fb2,fw2,err,ok)
            if(.not.ok)call abort_s2('scaled macro-direction face wiring failed',327)
            response=max(maxval(abs(fh1-fh0)),maxval(abs(fb1-fb0)),maxval(abs(fw1-fw0)))
            linearity=max(maxval(abs((fh2-fh0)-2*(fh1-fh0))),maxval(abs((fb2-fb0)-2*(fb1-fb0))), &
                maxval(abs((fw2-fw0)-2*(fw1-fw0))))
            write(unit,'(I0,3(",",ES24.16))')axis,response,linearity,err
            if(response<=100*epsilon(1.0_KREAL).or.linearity>gate) &
                call abort_s2('independent macro direction is not wired linearly',328)
        enddo
        close(unit)
    end subroutine write_face_wiring_check

    subroutine write_limiter_check()
        integer,parameter::n=9
        real(KREAL)::hf(n,n),bf(n,n),hface(4,3),bface(4,3),hnd(4,3),htd(4,3),bnd(4,3),btd(4,3)
        real(KREAL)::hvol(3,3),bvol(3,3),theta,hmean,bmean,havg,bavg,derivative_difference
        logical::ok,was_signed
        integer::i,j,a,c,unit
        do j=1,n;do i=1,n
            hf(i,j)=0.08_KREAL+0.22_KREAL*sin(1.7_KREAL*real(i,KREAL)+0.9_KREAL*real(j,KREAL))
            bf(i,j)=0.11_KREAL+0.27_KREAL*cos(0.8_KREAL*real(i,KREAL)-1.4_KREAL*real(j,KREAL))
        enddo;enddo
        hf(5,5)=0.06_KREAL;bf(5,5)=0.07_KREAL;hmean=hf(5,5);bmean=bf(5,5)
        call reconstruct_cell_pair(hf,bf,5,5,1.0_KREAL,1.0_KREAL,hface,bface,hnd,htd,bnd,btd, &
            theta,was_signed,ok)
        if(.not.ok)call abort_s2('common limiter reconstruction failed',329)
        call s2d1_cell9(hf,5,5,hvol,ok);call s2d1_cell9(bf,5,5,bvol,ok)
        havg=0.0_KREAL;bavg=0.0_KREAL
        do c=1,3;do a=1,3
            havg=havg+gw(a)*gw(c)*(hmean+theta*(hvol(a,c)-hmean))
            bavg=bavg+gw(a)*gw(c)*(bmean+theta*(bvol(a,c)-bmean))
        enddo;enddo
        derivative_difference=max(maxval(abs(hnd(1,:)-hnd(2,:))),maxval(abs(bnd(1,:)-bnd(2,:))))
        open(newunit=unit,file='limiter_checks.csv',status='replace',action='write')
        write(unit,'(A)')'theta,signed,min_h,min_b,h_mean_error,b_mean_error,left_right_derivative_difference'
        write(unit,'(ES24.16,",",L1,5(",",ES24.16))')theta,was_signed,minval(hface),minval(bface), &
            abs(havg-hmean),abs(bavg-bmean),derivative_difference
        close(unit)
        if(theta>=1.0_KREAL.or.was_signed.or.minval(hface)<-gate.or.minval(bface)<-gate.or. &
            abs(havg-hmean)>gate.or.abs(bavg-bmean)>gate.or.derivative_difference<=100*epsilon(1.0_KREAL)) &
            call abort_s2('common limiter activation test failed',330)
    end subroutine write_limiter_check

    subroutine RunS2D1()
        integer::np,ierr,nx,ny,nvel,i,j,k,step,unit
        real(KREAL)::dx,dy,cfl,x,y,rho,sincx,sincy,time,tstart,tend,init_start,init_end
        real(KREAL)::diag(12),flux_s,collision_s,flux_acc,collision_acc,check_s,output_s,t0c,t1c
        real(KREAL)::total0(4),total(4),drift,scale,mm,raw_node(16),raw_weight(16),mom(4)
        real(KREAL),allocatable::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),weight(:)
        real(KREAL),allocatable::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:)
        real(KREAL),allocatable::eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
        logical::ok
        character(len=24)::method
        call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
        if(np/=1)call abort_s2('S2D1 supports exactly one MPI rank',301)
        if(Nx_all/=Ny_all.or.Nx_all<5)call abort_s2('uniform square N>=5 is required',302)
        if(Nu/=16.or.Nv/=16)call abort_s2('full GH16 x GH16 is required',303)
        if(ck/=1)call abort_s2('CK must equal 1',304)
        if(any(BC_Type/=BCconst_periodic))call abort_s2('all four boundaries must be periodic',305)
        if(abs(Pr-1.0_KREAL)>1e-14_KREAL.and.abs(Pr-2.0_KREAL/3.0_KREAL)>1e-14_KREAL) &
            call abort_s2('Pr must equal 1 or 2/3',306)
        if(dt_global<=0.0_KREAL.or.MAX_ITER<1)call abort_s2('positive fixed dt and steps required',307)
        nx=Nx_all;ny=Ny_all;nvel=256;dx=2.0_KREAL*pi_d/real(nx,KREAL);dy=2.0_KREAL*pi_d/real(ny,KREAL)
        allocate(w(4,nx,ny),h(nx,ny,nvel),b(nx,ny,nvel),u(nvel),v(nvel),weight(nvel), &
            mh(nx,ny,nvel),mb(nx,ny,nvel),kxh(nx,ny,nvel),kxb(nx,ny,nvel), &
            kyh(nx,ny,nvel),kyb(nx,ny,nvel),eah(nx,ny,nvel),eab(nx,ny,nvel),kappa(2,2,nx,ny))
        call cpu_time(init_start)
        call init_velocity(u,v,weight,raw_node,raw_weight)
        call write_quadrature_check(raw_node,raw_weight,u,v,weight)
        call write_reconstruction_check()
        call write_local_qualification(u,v,weight)
        call write_equilibrium_direction_check(u,v,weight)
        call write_face_wiring_check(u,v,weight)
        call write_limiter_check()
        cfl=dt_global*(maxval(abs(u))/dx+maxval(abs(v))/dy)
        if(cfl>0.4_KREAL)call abort_s2('frozen dt violates transport CFL <= 0.4',308)
        sincx=sin(0.5_KREAL*dx)/(0.5_KREAL*dx);sincy=sin(0.5_KREAL*dy)/(0.5_KREAL*dy)
        do j=1,ny;y=(real(j,KREAL)-0.5_KREAL)*dy
            do i=1,nx;x=(real(i,KREAL)-0.5_KREAL)*dx
                select case(S2D1_initial_mode)
                case(0);rho=rho0+amp*sincx*sincy*sin(x+y)
                case(1);rho=rho0+amp*sincx*sin(x)
                case(2);rho=rho0+amp*sincy*sin(y)
                case(3);rho=rho0
                case default;call abort_s2('unsupported S2D1_initial_mode',319)
                end select
                w(:,i,j)=[rho,rho*u0,rho*v0,rho*(0.5_KREAL*(u0*u0+v0*v0)+1.5_KREAL*t0)]
                call make_mke(w(:,i,j),u,v,weight,mh(i,j,:),mb(i,j,:),kxh(i,j,:),kxb(i,j,:), &
                    kyh(i,j,:),kyb(i,j,:),eah(i,j,:),eab(i,j,:),kappa(:,:,i,j),ok,mm)
                if(.not.ok)call abort_s2('discrete conservative initialization failed',309)
                h(i,j,:)=mh(i,j,:);b(i,j,:)=mb(i,j,:)
            enddo
        enddo
        call build_targets(w,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,mm,ok)
        if(.not.ok)call abort_s2('initial volume target construction failed',310)
        call state_totals(w,total0)
        call write_velocity(u,v,weight,raw_node,raw_weight)
        call write_state('state_initial.bin',w,h,b,u,v,weight,0.0_KREAL)
        call write_field('field_initial.csv',w,h,b,eah,eab,u,v,weight,dx,dy,0.0_KREAL)
        call cpu_time(init_end)
        open(newunit=unit,file='steps.csv',status='replace',action='write')
        write(unit,'(A)')'step,time,face_moment,target_moment,cn_residual,macro_micro,limited,signed,'// &
            'dt_tau,cfl,h_min,h_max,b_min,b_max,global_drift,flux_s,collision_s'
        method=trim(adjustl(TimeIntegration_method));time=0.0_KREAL;flux_acc=0.0_KREAL
        collision_acc=0.0_KREAL;check_s=0.0_KREAL;output_s=0.0_KREAL
        call cpu_time(tstart)
        do step=1,MAX_ITER
            call advance_step(w,h,b,u,v,weight,dx,dy,dt_global,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa, &
                diag,flux_s,collision_s,ok)
            if(.not.ok)then
                close(unit);call write_state('last_accepted.bin',w,h,b,u,v,weight,time)
                call abort_s2('native S2/VL step rejected; last accepted state saved',311)
            endif
            time=time+dt_global;flux_acc=flux_acc+flux_s;collision_acc=collision_acc+collision_s
            call cpu_time(t0c);call state_totals(w,total)
            drift=maxval(abs(total-total0)/max(1.0_KREAL,abs(total0)))
            call cpu_time(t1c);check_s=check_s+t1c-t0c
            write(unit,'(I0,16(",",ES16.8))')step,time,diag(1:12),drift,flux_s,collision_s
            if(step==1.or.step==MAX_ITER)then
                call cpu_time(t0c)
                if(step==1)call write_state('state_prefix.bin',w,h,b,u,v,weight,time)
                call cpu_time(t1c);output_s=output_s+t1c-t0c
            endif
            if(drift>gate)then
                close(unit);call write_state('last_accepted.bin',w,h,b,u,v,weight,time)
                call abort_s2('periodic global conservation drift exceeded 2e-12',312)
            endif
            if(Nstep_show>0.and.mod(step,Nstep_show)==0.and.myID==0) &
                write(*,'(A,I0,A,F12.8,A,ES12.4)')'S2D1 step ',step,' t=',time,' drift=',drift
        enddo
        close(unit);call cpu_time(tend)
        call cpu_time(t0c)
        call write_state('state_final.bin',w,h,b,u,v,weight,time)
        call write_field('field_final.csv',w,h,b,eah,eab,u,v,weight,dx,dy,time)
        call cpu_time(t1c);output_s=output_s+t1c-t0c
        open(newunit=unit,file='timing.csv',status='replace',action='write')
        write(unit,'(A)')'method,N,Nv,Pr,dt,steps,time,cfl,initialization_s,flux_s,collision_s,check_s,output_s,total_loop_s'
        write(unit,'(A,",",I0,",",I0,11(",",ES16.8))')trim(method),nx,nvel,Pr,dt_global, &
            real(MAX_ITER,KREAL),time,cfl,init_end-init_start,flux_acc,collision_acc,check_s,output_s,tend-tstart
        close(unit)
        if(myID==0)write(*,'(A,1X,A,1X,A,ES12.4,A,ES12.4)') &
            'PASS_S2D1_PERIODIC_DIRECT_CN2',trim(method),'t=',time,' cfl=',cfl
        deallocate(w,h,b,u,v,weight,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa)
    end subroutine RunS2D1
end module mod_S2D1

module mod_S2D1Reconstruction
    use mod_S2MPI1Layout
    use Precision, only: KREAL
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: s2d1_face3, s2d1_side3, s2d1_cell9, s2d1_weno1d, s2d1_cell_faces9
    public :: FirstLayerPlane,prepare_first_layer_plane,s2d1_cell_faces9_cached,first_layer_plane_bytes
    public :: s2d1_cell_faces9_wall
    public :: VelocityTileWorkspace,VT_VELOCITIES,prepare_velocity_tile,reconstruct_velocity_tile_cell

    real(KREAL), parameter :: eps_weno=1.0e-8_KREAL
    real(KREAL), parameter :: gamma(4)=[0.01125_KREAL,0.1275_KREAL,0.01125_KREAL,0.85_KREAL]
    real(KREAL), parameter :: gnode(3)=[-0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL),0.0_KREAL, &
        0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL)]
    real(KREAL), save :: cmap(0:4,5,4)=0.0_KREAL
    real(KREAL), save :: wallmap(0:4,5,3,2)=0.0_KREAL
    logical, save :: maps_ready=.false.
    real(KREAL), save :: fixed_pow(8,0:4),fixed_dpow(8,0:4)

    ! Owned by one field and one velocity node during a single compute_flux.
    ! The caller prepares it again for each new state; invalid preparation
    ! cannot leave a previously valid plane consumable.
    type :: FirstLayerPlane
        private
        logical :: valid=.false.
        integer :: nx=0,ny=0
        logical :: scoped=.false., wall_scoped=.false.
        integer :: targets(4)=0
        logical,allocatable :: need_x(:,:),need_y(:,:)
        real(KREAL),allocatable :: xv(:,:,:),xd(:,:,:),rows(:,:,:),yv(:,:,:),yd(:,:,:)
    end type FirstLayerPlane

    integer,parameter :: VT_VELOCITIES=4,VT_FIELDS=2*VT_VELOCITIES,VT_LANES=8*VT_FIELDS
    type :: VelocityTileWorkspace
        type(FirstLayerPlane) :: cache(VT_FIELDS)
        real(KREAL) :: q(VT_LANES,5),comp(VT_LANES,5),fw(VT_LANES,4)
        real(KREAL) :: coef(VT_LANES,0:4,4),ccoef(VT_LANES,0:4,4)
        real(KREAL) :: cand(VT_LANES,4),dcand(VT_LANES,4),ccand(VT_LANES,4)
        real(KREAL) :: val(VT_LANES,3),der(VT_LANES,3),cval(VT_LANES,3)
        real(KREAL) :: face(4,3,VT_FIELDS),dn(4,3,VT_FIELDS),dt(4,3,VT_FIELDS),vol(3,3,VT_FIELDS)
        real(KREAL) :: pack_s=0,kernel_s=0,unpack_s=0,fallback_s=0
        integer(kind=8) :: packed_values=0,stored_values=0,cells=0,fallback_cells=0
        integer :: nf=0
        logical :: valid=.false.,profile=.false.,wall_batch=.false.
        real(KREAL) :: wallq(VT_LANES,5),wallc(VT_LANES,0:4)
        real(KREAL) :: wallv(VT_LANES,3),walld(VT_LANES,3)
    end type

contains

    integer pure function ip(i,n)
        integer, intent(in) :: i,n
        ip=modulo(i-1,n)+1
    end function ip

    subroutine init_maps()
        integer :: offsets5(5),offsets3(3),c,j,d,a
        real(KREAL) :: xi(8)
        if (maps_ready) return
        offsets3=[-2,-1,0]; call build_map(offsets3,cmap(:,:,1))
        offsets3=[-1,0,1];  call build_map(offsets3,cmap(:,:,2))
        offsets3=[0,1,2];   call build_map(offsets3,cmap(:,:,3))
        offsets5=[-2,-1,0,1,2]; call build_map(offsets5,cmap(:,:,4))
        do j=1,3
            offsets5=[1,2,3,4,5]-j
            call build_map(offsets5,wallmap(:,:,j,1),.true.)
            offsets5=[-4,-3,-2,-1,0]+j-1
            call build_map(offsets5,wallmap(:,:,j,2),.true.)
        enddo
        ! Same powers and operation order as prepared_eval, once per process.
        xi=[-0.5_KREAL,0.5_KREAL,gnode,0.0_KREAL,0.0_KREAL,0.0_KREAL]
        fixed_pow(:,0)=1.0_KREAL;fixed_dpow(:,0)=0.0_KREAL
        do d=1,4;do a=1,8
            fixed_pow(a,d)=fixed_pow(a,d-1)*xi(a)
            fixed_dpow(a,d)=real(d,KREAL)*xi(a)**(d-1)
        enddo;enddo
        maps_ready=.true.
    end subroutine init_maps

    subroutine build_map(offsets,map,ordered)
        integer, intent(in) :: offsets(:)
        real(KREAL), intent(out) :: map(0:4,5)
        logical,intent(in),optional::ordered
        real(KREAL) :: a(size(offsets),size(offsets)),rhs(size(offsets)),x(size(offsets))
        integer :: r,d,c
        map=0.0_KREAL
        do r=1,size(offsets)
            do d=0,size(offsets)-1
                a(r,d+1)=((real(offsets(r),KREAL)+0.5_KREAL)**(d+1)- &
                    (real(offsets(r),KREAL)-0.5_KREAL)**(d+1))/real(d+1,KREAL)
            end do
        end do
        do c=1,size(offsets)
            rhs=0.0_KREAL; rhs(c)=1.0_KREAL
            call solve_small(a,rhs,x)
            do d=0,size(offsets)-1
                if(present(ordered))then
                    if(ordered)then
                        map(d,c)=x(d+1)
                    else
                        map(d,offsets(c)+3)=x(d+1)
                    endif
                else
                    map(d,offsets(c)+3)=x(d+1)
                endif
            end do
        end do
    end subroutine build_map

    ! Five true cell averages on one side of the wall, expressed about cell j.
    ! The upper map index is distance from the upper wall (1,2,3).
    subroutine wall_coeff(q,j,ny,coef)
        real(KREAL),intent(in)::q(5)
        integer,intent(in)::j,ny
        real(KREAL),intent(out)::coef(0:4)
        call init_maps()
        if(global_j(j)<=3)then
            coef=matmul(wallmap(:,:,global_j(j),1),q)
        else
            coef=matmul(wallmap(:,:,global_ny(ny)-global_j(j)+1,2),q)
        endif
    end subroutine wall_coeff

    subroutine wall_eval(q,j,ny,xi,value,derivative)
        real(KREAL),intent(in)::q(5),xi
        integer,intent(in)::j,ny
        real(KREAL),intent(out)::value,derivative
        real(KREAL)::c(0:4)
        integer::d
        call wall_coeff(q,j,ny,c)
        value=c(4)
        do d=3,0,-1;value=value*xi+c(d);enddo
        derivative=4.0_KREAL*c(4)
        do d=3,1,-1;derivative=derivative*xi+real(d,KREAL)*c(d);enddo
    end subroutine wall_eval

    subroutine wall_samples(field,i,j,q)
        real(KREAL),intent(in)::field(:,:)
        integer,intent(in)::i,j
        real(KREAL),intent(out)::q(5)
        integer::s,first
        first=wall_first(j,size(field,2))
        do s=1,5;q(s)=field(i,first+s-1);enddo
    end subroutine wall_samples

    subroutine solve_small(a,b,x)
        real(KREAL), intent(in) :: a(:,:),b(:)
        real(KREAL), intent(out) :: x(:)
        real(KREAL) :: w(size(b),size(b)),r(size(b)),fac,tmp
        integer :: i,j,p,n,m
        n=size(b); w=a; r=b
        do p=1,n-1
            m=p-1+maxloc(abs(w(p:n,p)),dim=1)
            if (m/=p) then
                do j=1,n
                    tmp=w(p,j);w(p,j)=w(m,j);w(m,j)=tmp
                end do
                tmp=r(p);r(p)=r(m);r(m)=tmp
            end if
            do i=p+1,n
                fac=w(i,p)/w(p,p);w(i,p:n)=w(i,p:n)-fac*w(p,p:n);r(i)=r(i)-fac*r(p)
            end do
        end do
        x(n)=r(n)/w(n,n)
        do i=n-1,1,-1
            x(i)=(r(i)-sum(w(i,i+1:n)*x(i+1:n)))/w(i,i)
        end do
    end subroutine solve_small

    subroutine weno_coefficients(q,fw,coef,ok)
        real(KREAL),intent(in)::q(5)
        real(KREAL),intent(out)::fw(4),coef(0:4,4)
        logical,intent(out)::ok
        real(KREAL)::beta(4),tau,alpha(4),omega(4)
        integer::c,k,d
        call init_maps()
        beta(1)=13.0_KREAL/12.0_KREAL*(q(1)-2*q(2)+q(3))**2+0.25_KREAL*(q(1)-4*q(2)+3*q(3))**2
        beta(2)=13.0_KREAL/12.0_KREAL*(q(2)-2*q(3)+q(4))**2+0.25_KREAL*(q(2)-q(4))**2
        beta(3)=13.0_KREAL/12.0_KREAL*(q(3)-2*q(4)+q(5))**2+0.25_KREAL*(3*q(3)-4*q(4)+q(5))**2
        beta(4)=(231153*q(3)**2+104963*q(2)**2+6908*q(1)**2-38947*q(1)*q(4)+ &
            104963*q(4)**2+q(2)*(-51001*q(1)+179098*q(4)-38947*q(5))- &
            3*q(3)*(99692*q(2)-22641*q(1)+99692*q(4)-22641*q(5))+ &
            8209*q(1)*q(5)-51001*q(4)*q(5)+6908*q(5)**2)/5040.0_KREAL
        tau=sum(abs(beta(4)-beta(1:3)))/3.0_KREAL
        alpha=gamma*(1.0_KREAL+(tau/(beta+eps_weno))**2);omega=alpha/sum(alpha)
        fw(4)=omega(4)/gamma(4);fw(1:3)=omega(1:3)-fw(4)*gamma(1:3)
        ! The three quadratic candidates have exactly three active cells
        ! and zero cubic/quartic coefficients by build_map construction.
        coef=0.0_KREAL
        do c=1,3
            do k=c,c+2
                do d=0,2
                    coef(d,c)=coef(d,c)+cmap(d,k,c)*q(k)
                enddo
            enddo
        enddo
        do k=1,5
            do d=0,4
                coef(d,4)=coef(d,4)+cmap(d,k,4)*q(k)
            enddo
        enddo
        ok=all(ieee_is_finite(fw)).and.all(ieee_is_finite(coef))
    end subroutine weno_coefficients

    subroutine s2d1_weno1d(q,xi,value,derivative,ok)
        real(KREAL), intent(in) :: q(5),xi
        real(KREAL), intent(out) :: value,derivative
        logical, intent(out) :: ok
        real(KREAL) :: fw(4),coef(0:4,4),pow(0:4),dpow(0:4),cand(4),dcand(4)
        integer :: c,d
        call weno_coefficients(q,fw,coef,ok);if(.not.ok)return
        pow(0)=1.0_KREAL;dpow(0)=0.0_KREAL
        do d=1,4
            pow(d)=pow(d-1)*xi;dpow(d)=real(d,KREAL)*xi**(d-1)
        end do
        do c=1,4
            cand(c)=sum(coef(:,c)*pow);dcand(c)=sum(coef(:,c)*dpow)
        end do
        value=sum(fw*cand);derivative=sum(fw*dcand)
        ok=ieee_is_finite(value).and.ieee_is_finite(derivative)
    end subroutine s2d1_weno1d

    ! Evaluate a companion field with the nonlinear weights frozen from q.
    ! This is the tensor derivative of the same reconstructed expression;
    ! the companion does not compute a second set of WENO weights.

    ! The prepared object is local to one five-cell template and one field.
    ! In particular, no nonlinear weights are shared between distinct
    ! transverse templates, velocity nodes, fields, cells, or time levels.
    subroutine prepared_eval(fw,coef,xi,value,derivative,ok)
        real(KREAL),intent(in)::fw(4),coef(0:4,4),xi
        real(KREAL),intent(out)::value,derivative
        logical,intent(out)::ok
        real(KREAL)::pow(0:4),dpow(0:4),cand(4),dcand(4)
        integer::c,d
        pow(0)=1.0_KREAL;dpow(0)=0.0_KREAL
        do d=1,4
            pow(d)=pow(d-1)*xi;dpow(d)=real(d,KREAL)*xi**(d-1)
        enddo
        do c=1,4
            cand(c)=sum(coef(:,c)*pow);dcand(c)=sum(coef(:,c)*dpow)
        enddo
        value=sum(fw*cand);derivative=sum(fw*dcand)
        ok=ieee_is_finite(value).and.ieee_is_finite(derivative)
    end subroutine prepared_eval

    ! Independent fixed-node evaluations: SIMD across nodes, never across
    ! the degree/candidate reductions. Keep derivative checks at GQ3 too.
    subroutine prepared_fixed5(fw,coef,value,derivative,ok)
        real(KREAL),intent(in)::fw(4),coef(0:4,4)
        real(KREAL),intent(out)::value(5),derivative(5)
        logical,intent(out)::ok
        real(KREAL)::v(8),dv(8),cand(8),dcand(8)
        integer::a,c,d
        v=0.0_KREAL;dv=0.0_KREAL
        do c=1,4
            cand=0.0_KREAL;dcand=0.0_KREAL
            do d=0,4
                do a=1,8
                    cand(a)=cand(a)+coef(d,c)*fixed_pow(a,d)
                    dcand(a)=dcand(a)+coef(d,c)*fixed_dpow(a,d)
                enddo
            enddo
            do a=1,8
                v(a)=v(a)+fw(c)*cand(a)
                dv(a)=dv(a)+fw(c)*dcand(a)
            enddo
        enddo
        value=v(1:5);derivative=dv(1:5)
        ok=all(ieee_is_finite(value)).and.all(ieee_is_finite(derivative))
    end subroutine prepared_fixed5

    subroutine companion_eval3(q,companion,value,derivative,companion_value,ok)
        real(KREAL),intent(in)::q(5),companion(5)
        real(KREAL),intent(out)::value(3),derivative(3),companion_value(3)
        logical,intent(out)::ok
        real(KREAL)::fw(4),coef(0:4,4),ccoef(0:4,4),pow(0:4),dpow(0:4)
        real(KREAL)::cand(4),dcand(4),ccand(4)
        integer::s,c,d,k
        call weno_coefficients(q,fw,coef,ok);if(.not.ok)return
        ccoef=0.0_KREAL
        do c=1,3;do k=c,c+2;do d=0,2
            ccoef(d,c)=ccoef(d,c)+cmap(d,k,c)*companion(k)
        enddo;enddo;enddo
        do k=1,5;do d=0,4
            ccoef(d,4)=ccoef(d,4)+cmap(d,k,4)*companion(k)
        enddo;enddo
        do s=1,3
            pow(0)=1.0_KREAL;dpow(0)=0.0_KREAL
            do d=1,4;pow(d)=pow(d-1)*gnode(s);dpow(d)=real(d,KREAL)*gnode(s)**(d-1);enddo
            do c=1,4
                cand(c)=sum(coef(:,c)*pow);dcand(c)=sum(coef(:,c)*dpow)
                ccand(c)=sum(ccoef(:,c)*pow)
            enddo
            value(s)=sum(fw*cand);derivative(s)=sum(fw*dcand);companion_value(s)=sum(fw*ccand)
            if(.not.all(ieee_is_finite([value(s),derivative(s),companion_value(s)])))then
                ok=.false.;return
            endif
        enddo
        ok=.true.
    end subroutine companion_eval3

    subroutine s2d1_side3(field,axis,cell,tangent,side,value,dn,dt,ok,wall_y)
        real(KREAL), intent(in) :: field(:,:),side
        integer, intent(in) :: axis,cell,tangent
        real(KREAL), intent(out) :: value(3),dn(3),dt(3)
        logical, intent(out) :: ok
        logical, intent(in),optional :: wall_y
        real(KREAL) :: normal_value(5),normal_deriv(5),q(5)
        real(KREAL) :: f(4,3),n(4,3),t(4,3),vol(3,3)
        integer :: r,s,nx,ny,j,slot
        logical :: lok
        nx=size(field,1);ny=size(field,2);ok=.false.
        j=merge(tangent,cell,axis==1)
        if(present(wall_y))then
            if(wall_y.and.near_wall(j,ny))then
                if(axis==1)then
                    call s2d1_cell_faces9_wall(field,cell,j,f,n,t,vol,ok)
                    slot=merge(1,2,side<0.0_KREAL)
                else
                    call s2d1_cell_faces9_wall(field,tangent,j,f,n,t,vol,ok)
                    slot=merge(3,4,side<0.0_KREAL)
                endif
                if(ok)then;value=f(slot,:);dn=n(slot,:);dt=t(slot,:);endif
                return
            endif
        endif
        do r=-2,2
            do s=-2,2
                if(axis==1) q(s+3)=field(ip(cell+s,nx),ip(tangent+r,ny))
                if(axis==2) q(s+3)=field(ip(tangent+r,nx),ip(cell+s,ny))
            end do
            call s2d1_weno1d(q,side,normal_value(r+3),normal_deriv(r+3),lok)
            if(.not.lok)return
        end do
        call companion_eval3(normal_value,normal_deriv,value,dt,dn,lok)
        if(.not.lok)return
        ok=all(ieee_is_finite(value)).and.all(ieee_is_finite(dn)).and.all(ieee_is_finite(dt))
    end subroutine s2d1_side3

    subroutine s2d1_face3(field,axis,face,tangent,left,right,dnl,dnr,dtl,dtr,ok,wall_y)
        real(KREAL), intent(in) :: field(:,:)
        integer, intent(in) :: axis,face,tangent
        real(KREAL), intent(out) :: left(3),right(3),dnl(3),dnr(3),dtl(3),dtr(3)
        logical, intent(out) :: ok
        logical, intent(in),optional :: wall_y
        integer :: n,right_cell
        logical :: lok
        n=merge(size(field,1),size(field,2),axis==1);right_cell=ip(face+1,n)
        call s2d1_side3(field,axis,face,tangent,0.5_KREAL,left,dnl,dtl,lok,wall_y);if(.not.lok)then;ok=.false.;return;endif
        call s2d1_side3(field,axis,right_cell,tangent,-0.5_KREAL,right,dnr,dtr,lok,wall_y);ok=lok
    end subroutine s2d1_face3

    subroutine s2d1_cell9(field,i,j,value,ok,wall_y)
        real(KREAL), intent(in) :: field(:,:)
        integer, intent(in) :: i,j
        real(KREAL), intent(out) :: value(3,3)
        logical, intent(out) :: ok
        logical, intent(in),optional :: wall_y
        real(KREAL) :: rows(5,3),q(5),dummy,fw(4),coef(0:4,4)
        real(KREAL) :: face(4,3),dn(4,3),dt(4,3)
        integer :: r,s,a,b,nx,ny
        logical :: lok
        nx=size(field,1);ny=size(field,2);ok=.false.
        if(present(wall_y))then
            if(wall_y.and.near_wall(j,ny))then
                call s2d1_cell_faces9_wall(field,i,j,face,dn,dt,value,ok)
                return
            endif
        endif
        do r=-2,2
            do s=-2,2;q(s+3)=field(ip(i+s,nx),ip(j+r,ny));end do
            call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
            do a=1,3
                call prepared_eval(fw,coef,gnode(a),rows(r+3,a),dummy,lok);if(.not.lok)return
            end do
        end do
        do a=1,3
            call weno_coefficients(rows(:,a),fw,coef,lok);if(.not.lok)return
            do b=1,3
                call prepared_eval(fw,coef,gnode(b),value(a,b),dummy,lok);if(.not.lok)return
            end do
        end do
        ok=all(ieee_is_finite(value))
    end subroutine s2d1_cell9

    ! Fused cell-local x/y faces and x-then-y volume.  The normal x template
    ! used by the two x faces is exactly the first-stage volume template.
    ! The y-face templates are different and therefore prepared separately.
    subroutine s2d1_cell_faces9(field,i,j,face,dn,dt,volume,ok)
        real(KREAL),intent(in)::field(:,:)
        integer,intent(in)::i,j
        real(KREAL),intent(out)::face(4,3),dn(4,3),dt(4,3),volume(3,3)
        logical,intent(out)::ok
        real(KREAL)::xv(5,2),xd(5,2),yv(5,2),yd(5,2),rows(5,3)
        real(KREAL)::q(5),fw(4),coef(0:4,4),dummy
        integer::r,s,a,b,nx,ny
        logical::lok
        nx=size(field,1);ny=size(field,2);ok=.false.
        do r=-2,2
            do s=-2,2;q(s+3)=field(ip(i+s,nx),ip(j+r,ny));enddo
            call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
            call prepared_eval(fw,coef,-0.5_KREAL,xv(r+3,1),xd(r+3,1),lok);if(.not.lok)return
            call prepared_eval(fw,coef, 0.5_KREAL,xv(r+3,2),xd(r+3,2),lok);if(.not.lok)return
            do a=1,3
                call prepared_eval(fw,coef,gnode(a),rows(r+3,a),dummy,lok);if(.not.lok)return
            enddo
        enddo
        do r=-2,2
            do s=-2,2;q(s+3)=field(ip(i+r,nx),ip(j+s,ny));enddo
            call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
            call prepared_eval(fw,coef,-0.5_KREAL,yv(r+3,1),yd(r+3,1),lok);if(.not.lok)return
            call prepared_eval(fw,coef, 0.5_KREAL,yv(r+3,2),yd(r+3,2),lok);if(.not.lok)return
        enddo
        do s=1,2
            call companion_eval3(xv(:,s),xd(:,s),face(s,:),dt(s,:),dn(s,:),lok);if(.not.lok)return
            call companion_eval3(yv(:,s),yd(:,s),face(s+2,:),dt(s+2,:),dn(s+2,:),lok);if(.not.lok)return
        enddo
        do a=1,3
            call weno_coefficients(rows(:,a),fw,coef,lok);if(.not.lok)return
            do b=1,3
                call prepared_eval(fw,coef,gnode(b),volume(a,b),dummy,lok);if(.not.lok)return
            enddo
        enddo
        ok=all(ieee_is_finite(face)).and.all(ieee_is_finite(dn)).and. &
            all(ieee_is_finite(dt)).and.all(ieee_is_finite(volume))
    end subroutine s2d1_cell_faces9

    ! Tensor reconstruction with a one-sided quartic in y and frozen O3 in x.
    ! No ghost values or opposite-wall averages enter this routine.
    subroutine s2d1_cell_faces9_wall(field,i,j,face,dn,dt,volume,ok)
        real(KREAL),intent(in)::field(:,:)
        integer,intent(in)::i,j
        real(KREAL),intent(out)::face(4,3),dn(4,3),dt(4,3),volume(3,3)
        logical,intent(out)::ok
        real(KREAL)::xv(5,2),xd(5,2),rows(5,3),yv(5,2),yd(5,2)
        real(KREAL)::q(5),fw(4),coef(0:4,4),dummy
        integer::r,s,a,b,nx,ny,first
        logical::lok
        nx=size(field,1);ny=size(field,2);ok=.false.
        if(ny<5.or.j<1.or.j>ny)return
        if(.not.near_wall(j,ny))then
            call s2d1_cell_faces9(field,i,j,face,dn,dt,volume,ok)
            return
        endif
        first=wall_first(j,ny)
        do r=1,5
            do s=-2,2;q(s+3)=field(ip(i+s,nx),first+r-1);enddo
            call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
            call prepared_eval(fw,coef,-0.5_KREAL,xv(r,1),xd(r,1),lok);if(.not.lok)return
            call prepared_eval(fw,coef, 0.5_KREAL,xv(r,2),xd(r,2),lok);if(.not.lok)return
            do a=1,3
                call prepared_eval(fw,coef,gnode(a),rows(r,a),dummy,lok);if(.not.lok)return
            enddo
        enddo
        do r=-2,2
            call wall_samples(field,ip(i+r,nx),j,q)
            call wall_eval(q,j,ny,-0.5_KREAL,yv(r+3,1),yd(r+3,1))
            call wall_eval(q,j,ny, 0.5_KREAL,yv(r+3,2),yd(r+3,2))
        enddo
        do s=1,2
            do a=1,3
                call wall_eval(xv(:,s),j,ny,gnode(a),face(s,a),dt(s,a))
                call wall_eval(xd(:,s),j,ny,gnode(a),dn(s,a),dummy)
            enddo
            call companion_eval3(yv(:,s),yd(:,s),face(s+2,:),dt(s+2,:),dn(s+2,:),lok)
            if(.not.lok)return
        enddo
        do a=1,3;do b=1,3
            call wall_eval(rows(:,a),j,ny,gnode(b),volume(a,b),dummy)
        enddo;enddo
        ok=all(ieee_is_finite(face)).and.all(ieee_is_finite(dn)).and. &
            all(ieee_is_finite(dt)).and.all(ieee_is_finite(volume))
    end subroutine s2d1_cell_faces9_wall

    subroutine prepare_first_layer_plane(cache,field,ok,wall_y,target_box,wall_scope)
        type(FirstLayerPlane),intent(inout)::cache
        real(KREAL),intent(in)::field(:,:)
        logical,intent(out)::ok
        logical,intent(in),optional::wall_y
        integer,intent(in),optional::target_box(4)
        logical,intent(in),optional::wall_scope
        real(KREAL)::q(5),fw(4),coef(0:4,4),dummy,vfixed(5),dfixed(5)
        integer::nx,ny,i,j,s,a,ii,jj,first
        logical::lok,expand_wall
        expand_wall=.false.;if(present(wall_scope))expand_wall=wall_scope
        nx=size(field,1);ny=size(field,2);ok=.false.;cache%valid=.false.
        if(nx<1.or.ny<1)return
        if(cache%nx/=nx.or.cache%ny/=ny.or..not.allocated(cache%xv))then
            if(allocated(cache%xv))deallocate(cache%xv)
            if(allocated(cache%xd))deallocate(cache%xd)
            if(allocated(cache%rows))deallocate(cache%rows)
            if(allocated(cache%yv))deallocate(cache%yv)
            if(allocated(cache%yd))deallocate(cache%yd)
            allocate(cache%xv(nx,ny,2),cache%xd(nx,ny,2),cache%rows(nx,ny,3), &
                cache%yv(nx,ny,2),cache%yd(nx,ny,2))
            cache%nx=nx;cache%ny=ny
            if(allocated(cache%need_x))deallocate(cache%need_x,cache%need_y)
            allocate(cache%need_x(nx,ny),cache%need_y(nx,ny))
            cache%targets=0
        endif
        cache%scoped=present(target_box)
        if(present(target_box))then
            if(target_box(1)<1.or.target_box(2)>nx.or.target_box(1)>target_box(2))return
            if(target_box(3)<1.or.target_box(4)>ny.or.target_box(3)>target_box(4))return
            ! Validate raw input, including unconsumed padded entries.
            if(.not.all(ieee_is_finite(field)))return
            if(any(cache%targets/=target_box).or.(cache%wall_scoped.neqv.expand_wall))then
                cache%need_x=.false.;cache%need_y=.false.
                do jj=target_box(3),target_box(4);do ii=target_box(1),target_box(2)
                    do s=-2,2
                        cache%need_x(ii,ip(jj+s,ny))=.true.
                        cache%need_y(ip(ii+s,nx),jj)=.true.
                    enddo
                    if(expand_wall.and.near_wall(jj,ny))then
                        first=wall_first(jj,ny)
                        cache%need_x(ii,first:first+4)=.true.
                    endif
                enddo;enddo
                cache%targets=target_box;cache%wall_scoped=expand_wall
            endif
        endif
        do j=1,ny
            do i=1,nx
                if(cache%scoped)then
                    if(.not.cache%need_x(i,j))goto 110
                endif
                do s=-2,2;q(s+3)=field(ip(i+s,nx),j);enddo
                call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
                call prepared_fixed5(fw,coef,vfixed,dfixed,lok);if(.not.lok)return
                cache%xv(i,j,:)=vfixed(1:2);cache%xd(i,j,:)=dfixed(1:2)
                cache%rows(i,j,:)=vfixed(3:5)
110             continue
                if(cache%scoped)then
                    if(.not.cache%need_y(i,j))cycle
                endif
                if((slab_enabled.or.cart_context).and.(j<yrecon_lo().or.j>recon_hi(ny)))cycle
                if(present(wall_y))then
                    if(wall_y.and.near_wall(j,ny))then
                        call wall_samples(field,i,j,q)
                        call wall_eval(q,j,ny,-0.5_KREAL,cache%yv(i,j,1),cache%yd(i,j,1))
                        call wall_eval(q,j,ny, 0.5_KREAL,cache%yv(i,j,2),cache%yd(i,j,2))
                        cycle
                    endif
                endif
                do s=-2,2;q(s+3)=field(i,ip(j+s,ny));enddo
                call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
                call prepared_eval(fw,coef,-0.5_KREAL,cache%yv(i,j,1),cache%yd(i,j,1),lok);if(.not.lok)return
                call prepared_eval(fw,coef, 0.5_KREAL,cache%yv(i,j,2),cache%yd(i,j,2),lok);if(.not.lok)return
            enddo
        enddo
        cache%valid=.true.;ok=.true.
    end subroutine prepare_first_layer_plane

    ! Seven independent tensor templates: four face/companion pairs and
    ! three volume columns. Lane 8 is padding. No shared nonlinear weights.
    ! Lane is contiguous; all stencil/coefficient reductions retain order.
    ! Scalar nonlinear arithmetic is copied verbatim from weno_coefficients.
    ! Only independent polynomial coefficient construction is batched below.
    subroutine weno_weights_scalar(q,fw,ok)
        real(KREAL),intent(in)::q(5)
        real(KREAL),intent(out)::fw(4)
        logical,intent(out)::ok
        real(KREAL)::beta(4),tau,alpha(4),omega(4)
        beta(1)=13.0_KREAL/12.0_KREAL*(q(1)-2*q(2)+q(3))**2+0.25_KREAL*(q(1)-4*q(2)+3*q(3))**2
        beta(2)=13.0_KREAL/12.0_KREAL*(q(2)-2*q(3)+q(4))**2+0.25_KREAL*(q(2)-q(4))**2
        beta(3)=13.0_KREAL/12.0_KREAL*(q(3)-2*q(4)+q(5))**2+0.25_KREAL*(3*q(3)-4*q(4)+q(5))**2
        beta(4)=(231153*q(3)**2+104963*q(2)**2+6908*q(1)**2-38947*q(1)*q(4)+ &
            104963*q(4)**2+q(2)*(-51001*q(1)+179098*q(4)-38947*q(5))- &
            3*q(3)*(99692*q(2)-22641*q(1)+99692*q(4)-22641*q(5))+ &
            8209*q(1)*q(5)-51001*q(4)*q(5)+6908*q(5)**2)/5040.0_KREAL
        tau=sum(abs(beta(4)-beta(1:3)))/3.0_KREAL
        alpha=gamma*(1.0_KREAL+(tau/(beta+eps_weno))**2);omega=alpha/sum(alpha)
        fw(4)=omega(4)/gamma(4);fw(1:3)=omega(1:3)-fw(4)*gamma(1:3)
        ok=all(ieee_is_finite(fw))
    end subroutine weno_weights_scalar

    subroutine tensor_batch7(q,comp,value,derivative,comp_value,ok)
        real(KREAL),intent(in)::q(8,5),comp(8,5)
        real(KREAL),intent(out)::value(8,3),derivative(8,3),comp_value(8,3)
        logical,intent(out)::ok
        real(KREAL)::fw(8,4),lane_fw(4),lane_coef(0:4,4)
        real(KREAL)::coef(8,0:4,4),ccoef(8,0:4,4)
        real(KREAL)::cand(8,4),dcand(8,4),ccand(8,4),pow(0:4),dpow(0:4)
        integer::l,c,d,k,a
        ! Each lane keeps its own original scalar weights. Coefficient
        ! accumulation retains k order; independent lanes are contiguous.
        call init_maps()
        fw=0.0_KREAL;coef=0.0_KREAL;ccoef=0.0_KREAL
        do l=1,7
            call weno_weights_scalar(q(l,:),lane_fw,ok)
            if(.not.ok)return
            fw(l,:)=lane_fw
        enddo
        do c=1,3;do k=c,c+2;do d=0,2
            do l=1,8
                coef(l,d,c)=coef(l,d,c)+cmap(d,k,c)*q(l,k)
            enddo
        enddo;enddo;enddo
        do k=1,5;do d=0,4
            do l=1,8
                coef(l,d,4)=coef(l,d,4)+cmap(d,k,4)*q(l,k)
            enddo
        enddo;enddo
        ok=all(ieee_is_finite(coef));if(.not.ok)return
        do c=1,3;do k=c,c+2;do d=0,2
            do l=1,4
                ccoef(l,d,c)=ccoef(l,d,c)+cmap(d,k,c)*comp(l,k)
            enddo
        enddo;enddo;enddo
        do k=1,5;do d=0,4
            do l=1,4
                ccoef(l,d,4)=ccoef(l,d,4)+cmap(d,k,4)*comp(l,k)
            enddo
        enddo;enddo
        do a=1,3
            pow(0)=1.0_KREAL;dpow(0)=0.0_KREAL
            do d=1,4
                pow(d)=pow(d-1)*gnode(a);dpow(d)=real(d,KREAL)*gnode(a)**(d-1)
            enddo
            cand=0.0_KREAL;dcand=0.0_KREAL;ccand=0.0_KREAL
            do c=1,4;do d=0,merge(2,4,c<=3)
                do l=1,8
                    cand(l,c)=cand(l,c)+coef(l,d,c)*pow(d)
                    dcand(l,c)=dcand(l,c)+coef(l,d,c)*dpow(d)
                enddo
                do l=1,4
                    ccand(l,c)=ccand(l,c)+ccoef(l,d,c)*pow(d)
                enddo
            enddo;enddo
            value(:,a)=0.0_KREAL;derivative(:,a)=0.0_KREAL;comp_value(:,a)=0.0_KREAL
            do c=1,4
                do l=1,8
                    value(l,a)=value(l,a)+fw(l,c)*cand(l,c)
                    derivative(l,a)=derivative(l,a)+fw(l,c)*dcand(l,c)
                enddo
                do l=1,4
                    comp_value(l,a)=comp_value(l,a)+fw(l,c)*ccand(l,c)
                enddo
            enddo
        enddo
        ok=all(ieee_is_finite(value)).and.all(ieee_is_finite(derivative)).and.all(ieee_is_finite(comp_value))
    end subroutine tensor_batch7

    subroutine s2d1_cell_faces9_cached(cache,i,j,face,dn,dt,volume,ok,field,wall_y)
        type(FirstLayerPlane),intent(in)::cache
        integer,intent(in)::i,j
        real(KREAL),intent(out)::face(4,3),dn(4,3),dt(4,3),volume(3,3)
        logical,intent(out)::ok
        real(KREAL),intent(in),optional::field(:,:)
        logical,intent(in),optional::wall_y
        real(KREAL)::xv(5,2),xd(5,2),yv(5,2),yd(5,2),rows(5,3)
        real(KREAL)::qbatch(8,5),cbatch(8,5),values(8,3),derivs(8,3),cvalues(8,3)
        integer::r,s,a,b
        logical::lok
        ok=.false.
        if(.not.cache%valid.or.i<1.or.i>cache%nx.or.j<1.or.j>cache%ny)return
        if(cache%scoped)then
            if(i<cache%targets(1).or.i>cache%targets(2).or. &
               j<cache%targets(3).or.j>cache%targets(4))return
        endif
        if(present(wall_y).and.present(field))then
            if(wall_y.and.near_wall(j,cache%ny))then
                call s2d1_cell_faces9_wall(field,i,j,face,dn,dt,volume,ok)
                return
            endif
        endif
        do r=-2,2
            xv(r+3,:)=cache%xv(i,ip(j+r,cache%ny),:)
            xd(r+3,:)=cache%xd(i,ip(j+r,cache%ny),:)
            rows(r+3,:)=cache%rows(i,ip(j+r,cache%ny),:)
        enddo
        do r=-2,2
            yv(r+3,:)=cache%yv(ip(i+r,cache%nx),j,:)
            yd(r+3,:)=cache%yd(ip(i+r,cache%nx),j,:)
        enddo
        qbatch=0.0_KREAL;cbatch=0.0_KREAL
        do s=1,2
            qbatch(2*s-1,:)=xv(:,s);cbatch(2*s-1,:)=xd(:,s)
            qbatch(2*s,:)=yv(:,s);cbatch(2*s,:)=yd(:,s)
        enddo
        do a=1,3
            qbatch(4+a,:)=rows(:,a)
        enddo
        call tensor_batch7(qbatch,cbatch,values,derivs,cvalues,lok)
        if(.not.lok)return
        do s=1,2
            face(s,:)=values(2*s-1,:);dt(s,:)=derivs(2*s-1,:);dn(s,:)=cvalues(2*s-1,:)
            face(s+2,:)=values(2*s,:);dt(s+2,:)=derivs(2*s,:);dn(s+2,:)=cvalues(2*s,:)
        enddo
        do a=1,3
            volume(a,:)=values(4+a,:)
        enddo
        ok=all(ieee_is_finite(face)).and.all(ieee_is_finite(dn)).and. &
            all(ieee_is_finite(dt)).and.all(ieee_is_finite(volume))
    end subroutine s2d1_cell_faces9_cached

    subroutine tensor_velocity_tile(nf,q,comp,value,derivative,comp_value,fw,coef,ccoef,cand,dcand,ccand,ok,groups)
        integer,intent(in)::nf
        real(KREAL),intent(in)::q(VT_LANES,5),comp(VT_LANES,5)
        real(KREAL),intent(out)::value(VT_LANES,3),derivative(VT_LANES,3),comp_value(VT_LANES,3)
        real(KREAL),intent(out)::fw(VT_LANES,4),coef(VT_LANES,0:4,4),ccoef(VT_LANES,0:4,4)
        real(KREAL),intent(out)::cand(VT_LANES,4),dcand(VT_LANES,4),ccand(VT_LANES,4)
        logical,intent(out)::ok
        real(KREAL)::lane_fw(4),pow(0:4),dpow(0:4)
        integer::l,c,d,k,a,group,field,ng
        integer,intent(in),optional::groups
        ! Each lane keeps its own original scalar weights. Coefficient
        ! accumulation retains k order; independent lanes are contiguous.
        call init_maps()
        fw=0.0_KREAL;coef=0.0_KREAL;ccoef=0.0_KREAL
        ng=7;if(present(groups))ng=groups
        do group=0,ng-1;do field=1,nf
            l=group*VT_FIELDS+field
            call weno_weights_scalar(q(l,:),lane_fw,ok)
            if(.not.ok)return
            fw(l,:)=lane_fw
        enddo;enddo
        do c=1,3;do k=c,c+2;do d=0,2
            do l=1,VT_LANES
                coef(l,d,c)=coef(l,d,c)+cmap(d,k,c)*q(l,k)
            enddo
        enddo;enddo;enddo
        do k=1,5;do d=0,4
            do l=1,VT_LANES
                coef(l,d,4)=coef(l,d,4)+cmap(d,k,4)*q(l,k)
            enddo
        enddo;enddo
        ok=all(ieee_is_finite(coef));if(.not.ok)return
        do c=1,3;do k=c,c+2;do d=0,2
            do l=1,4*VT_FIELDS
                ccoef(l,d,c)=ccoef(l,d,c)+cmap(d,k,c)*comp(l,k)
            enddo
        enddo;enddo;enddo
        do k=1,5;do d=0,4
            do l=1,4*VT_FIELDS
                ccoef(l,d,4)=ccoef(l,d,4)+cmap(d,k,4)*comp(l,k)
            enddo
        enddo;enddo
        do a=1,3
            pow(0)=1.0_KREAL;dpow(0)=0.0_KREAL
            do d=1,4
                pow(d)=pow(d-1)*gnode(a);dpow(d)=real(d,KREAL)*gnode(a)**(d-1)
            enddo
            cand=0.0_KREAL;dcand=0.0_KREAL;ccand=0.0_KREAL
            do c=1,4;do d=0,merge(2,4,c<=3)
                do l=1,VT_LANES
                    cand(l,c)=cand(l,c)+coef(l,d,c)*pow(d)
                    dcand(l,c)=dcand(l,c)+coef(l,d,c)*dpow(d)
                enddo
                do l=1,4*VT_FIELDS
                    ccand(l,c)=ccand(l,c)+ccoef(l,d,c)*pow(d)
                enddo
            enddo;enddo
            value(:,a)=0.0_KREAL;derivative(:,a)=0.0_KREAL;comp_value(:,a)=0.0_KREAL
            do c=1,4
                do l=1,VT_LANES
                    value(l,a)=value(l,a)+fw(l,c)*cand(l,c)
                    derivative(l,a)=derivative(l,a)+fw(l,c)*dcand(l,c)
                enddo
                do l=1,4*VT_FIELDS
                    comp_value(l,a)=comp_value(l,a)+fw(l,c)*ccand(l,c)
                enddo
            enddo
        enddo
        ok=all(ieee_is_finite(value)).and.all(ieee_is_finite(derivative)).and.all(ieee_is_finite(comp_value))
    end subroutine tensor_velocity_tile

    ! First layer remains state-local and scalar-exact. Each field is independent.
    subroutine prepare_velocity_tile(ws,h,b,k0,nt,box,wall,ok)
        type(VelocityTileWorkspace),intent(inout)::ws
        real(KREAL),intent(in)::h(:,:,:),b(:,:,:)
        integer,intent(in)::k0,nt,box(4)
        logical,intent(in)::wall
        logical,intent(out)::ok
        integer::v
        ws%valid=.false.;ok=.false.;ws%nf=0
        if(nt<1.or.nt>VT_VELOCITIES.or.k0<1.or.k0+nt-1>size(h,3))return
        if(any(shape(h)/=shape(b)))return
        do v=1,nt
            call prepare_first_layer_plane(ws%cache(2*v-1),h(:,:,k0+v-1),ok,wall,box,ws%wall_batch.and.wall)
            if(.not.ok)return
            call prepare_first_layer_plane(ws%cache(2*v),b(:,:,k0+v-1),ok,wall,box,ws%wall_batch.and.wall)
            if(.not.ok)return
        enddo
        ws%nf=2*nt;ws%valid=.true.;ok=.true.
    end subroutine

    subroutine reconstruct_velocity_tile_cell(ws,h,b,k0,i,j,wall,ok)
        type(VelocityTileWorkspace),intent(inout)::ws
        real(KREAL),intent(in)::h(:,:,:),b(:,:,:)
        integer,intent(in)::k0,i,j
        logical,intent(in)::wall
        logical,intent(out)::ok
        integer::f,r,s,a,l,k,nx,ny
        real(KREAL)::t0,t1
        ok=.false.
        if(.not.ws%valid)return
        nx=ws%cache(1)%nx;ny=ws%cache(1)%ny
        if(i<ws%cache(1)%targets(1).or.i>ws%cache(1)%targets(2).or. &
           j<ws%cache(1)%targets(3).or.j>ws%cache(1)%targets(4))return
        if(wall.and.near_wall(j,ny))then
            if(ws%wall_batch)then
                call reconstruct_wall_tile(ws,i,j,ok)
                return
            endif
            if(ws%profile)call cpu_time(t0)
            do f=1,ws%nf
                k=k0+(f-1)/2
                if(mod(f,2)==1)then
                    call s2d1_cell_faces9_cached(ws%cache(f),i,j,ws%face(:,:,f),ws%dn(:,:,f), &
                         ws%dt(:,:,f),ws%vol(:,:,f),ok,h(:,:,k),wall)
                else
                    call s2d1_cell_faces9_cached(ws%cache(f),i,j,ws%face(:,:,f),ws%dn(:,:,f), &
                         ws%dt(:,:,f),ws%vol(:,:,f),ok,b(:,:,k),wall)
                endif
                if(.not.ok)return
            enddo
            if(ws%profile)then
                call cpu_time(t1);ws%fallback_s=ws%fallback_s+t1-t0
                ws%fallback_cells=ws%fallback_cells+1
            endif
            return
        endif
        if(ws%profile)call cpu_time(t0)
        ws%q=0;ws%comp=0
        ! Lane index: template*VT_FIELDS+field. Velocity/h/b is contiguous.
        do r=-2,2
            do s=1,2
                do f=1,ws%nf
                    l=(2*s-2)*VT_FIELDS+f
                    ws%q(l,r+3)=ws%cache(f)%xv(i,ip(j+r,ny),s)
                    ws%comp(l,r+3)=ws%cache(f)%xd(i,ip(j+r,ny),s)
                    l=(2*s-1)*VT_FIELDS+f
                    ws%q(l,r+3)=ws%cache(f)%yv(ip(i+r,nx),j,s)
                    ws%comp(l,r+3)=ws%cache(f)%yd(ip(i+r,nx),j,s)
                enddo
            enddo
            do a=1,3;do f=1,ws%nf
                ws%q((3+a)*VT_FIELDS+f,r+3)=ws%cache(f)%rows(i,ip(j+r,ny),a)
            enddo;enddo
        enddo
        if(ws%profile)then
            call cpu_time(t1);ws%pack_s=ws%pack_s+t1-t0;t0=t1
            ws%packed_values=ws%packed_values+55*ws%nf
        endif
        call tensor_velocity_tile(ws%nf,ws%q,ws%comp,ws%val,ws%der,ws%cval, &
            ws%fw,ws%coef,ws%ccoef,ws%cand,ws%dcand,ws%ccand,ok);if(.not.ok)return
        if(ws%profile)then
            call cpu_time(t1);ws%kernel_s=ws%kernel_s+t1-t0;t0=t1
        endif
        do s=1,2;do a=1,3;do f=1,ws%nf
            l=(2*s-2)*VT_FIELDS+f
            ws%face(s,a,f)=ws%val(l,a);ws%dt(s,a,f)=ws%der(l,a);ws%dn(s,a,f)=ws%cval(l,a)
            l=(2*s-1)*VT_FIELDS+f
            ws%face(s+2,a,f)=ws%val(l,a);ws%dt(s+2,a,f)=ws%der(l,a);ws%dn(s+2,a,f)=ws%cval(l,a)
        enddo;enddo;enddo
        do s=1,3;do a=1,3;do f=1,ws%nf
            ws%vol(s,a,f)=ws%val((3+s)*VT_FIELDS+f,a)
        enddo;enddo;enddo
        if(ws%profile)then
            call cpu_time(t1);ws%unpack_s=ws%unpack_s+t1-t0
            ws%stored_values=ws%stored_values+45*ws%nf;ws%cells=ws%cells+1
        endif
    end subroutine


    ! The wall map and its reduction are unchanged. Compute each unique map
    ! product once; Horner evaluation vectorizes across independent fields.
    subroutine wall_tile_kernel(nf,q,j,ny,c,value,derivative)
        integer,intent(in)::nf,j,ny
        real(KREAL),intent(in)::q(VT_LANES,5)
        real(KREAL),intent(out)::c(VT_LANES,0:4),value(VT_LANES,3),derivative(VT_LANES,3)
        real(KREAL)::cq(0:4)
        integer::g,f,l,a,d
        c=0.0_KREAL
        do g=0,6;do f=1,nf
            l=g*VT_FIELDS+f
            call wall_coeff(q(l,:),j,ny,cq)
            c(l,:)=cq
        enddo;enddo
        do a=1,3
            value(:,a)=c(:,4);derivative(:,a)=4.0_KREAL*c(:,4)
            do d=3,0,-1;do l=1,VT_LANES
                value(l,a)=value(l,a)*gnode(a)+c(l,d)
            enddo;enddo
            do d=3,1,-1;do l=1,VT_LANES
                derivative(l,a)=derivative(l,a)*gnode(a)+real(d,KREAL)*c(l,d)
            enddo;enddo
        enddo
    end subroutine

    subroutine reconstruct_wall_tile(ws,i,j,ok)
        type(VelocityTileWorkspace),intent(inout)::ws
        integer,intent(in)::i,j
        logical,intent(out)::ok
        integer::first,r,f,s,a,l,nx,ny
        real(KREAL)::t0,t1
        ok=.false.;nx=ws%cache(1)%nx;ny=ws%cache(1)%ny
        first=wall_first(j,ny)
        if(ws%profile)call cpu_time(t0)
        ws%q=0.0_KREAL;ws%comp=0.0_KREAL;ws%wallq=0.0_KREAL
        ! Use exactly the same five true wall-side rows; no opposite wall,
        ! no x-uniform assumption and no cross-field nonlinear weight sharing.
        do r=1,5
            do s=1,2;do f=1,ws%nf
                l=(s-1)*VT_FIELDS+f
                ws%wallq(l,r)=ws%cache(f)%xv(i,first+r-1,s)
                ws%wallq((s+1)*VT_FIELDS+f,r)=ws%cache(f)%xd(i,first+r-1,s)
                ws%q(l,r)=ws%cache(f)%yv(ip(i+r-3,nx),j,s)
                ws%comp(l,r)=ws%cache(f)%yd(ip(i+r-3,nx),j,s)
            enddo;enddo
            do a=1,3;do f=1,ws%nf
                ws%wallq((3+a)*VT_FIELDS+f,r)=ws%cache(f)%rows(i,first+r-1,a)
            enddo;enddo
        enddo
        if(ws%profile)then
            call cpu_time(t1);ws%pack_s=ws%pack_s+t1-t0;t0=t1
            ws%packed_values=ws%packed_values+55*ws%nf
        endif
        call wall_tile_kernel(ws%nf,ws%wallq,j,ny,ws%wallc,ws%wallv,ws%walld)
        call tensor_velocity_tile(ws%nf,ws%q,ws%comp,ws%val,ws%der,ws%cval, &
            ws%fw,ws%coef,ws%ccoef,ws%cand,ws%dcand,ws%ccand,ok,2)
        if(.not.ok)return
        if(.not.all(ieee_is_finite(ws%wallv)).or..not.all(ieee_is_finite(ws%walld)))then
            ok=.false.;return
        endif
        if(ws%profile)then
            call cpu_time(t1);ws%kernel_s=ws%kernel_s+t1-t0;t0=t1
        endif
        do s=1,2;do a=1,3;do f=1,ws%nf
            l=(s-1)*VT_FIELDS+f
            ws%face(s,a,f)=ws%wallv(l,a);ws%dt(s,a,f)=ws%walld(l,a)
            ws%dn(s,a,f)=ws%wallv((s+1)*VT_FIELDS+f,a)
            ws%face(s+2,a,f)=ws%val(l,a);ws%dt(s+2,a,f)=ws%der(l,a)
            ws%dn(s+2,a,f)=ws%cval(l,a)
        enddo;enddo;enddo
        do s=1,3;do a=1,3;do f=1,ws%nf
            ws%vol(s,a,f)=ws%wallv((3+s)*VT_FIELDS+f,a)
        enddo;enddo;enddo
        if(ws%profile)then
            call cpu_time(t1);ws%unpack_s=ws%unpack_s+t1-t0
            ws%stored_values=ws%stored_values+45*ws%nf;ws%cells=ws%cells+1
        endif
        ok=.true.
    end subroutine

    integer(kind=8) function first_layer_plane_bytes(cache) result(bytes)
        type(FirstLayerPlane),intent(in)::cache
        bytes=0_8
        if(allocated(cache%xv))bytes=bytes+int(size(cache%xv),8)*storage_size(0.0_KREAL)/8
        if(allocated(cache%xd))bytes=bytes+int(size(cache%xd),8)*storage_size(0.0_KREAL)/8
        if(allocated(cache%rows))bytes=bytes+int(size(cache%rows),8)*storage_size(0.0_KREAL)/8
        if(allocated(cache%yv))bytes=bytes+int(size(cache%yv),8)*storage_size(0.0_KREAL)/8
        if(allocated(cache%yd))bytes=bytes+int(size(cache%yd),8)*storage_size(0.0_KREAL)/8
    end function first_layer_plane_bytes
end module mod_S2D1Reconstruction

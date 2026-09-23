module mod_S2D1Reconstruction
    use Precision, only: KREAL
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: s2d1_face3, s2d1_side3, s2d1_cell9, s2d1_weno1d, s2d1_cell_faces9
    public :: FirstLayerPlane,prepare_first_layer_plane,s2d1_cell_faces9_cached,first_layer_plane_bytes

    real(KREAL), parameter :: eps_weno=1.0e-8_KREAL
    real(KREAL), parameter :: gamma(4)=[0.01125_KREAL,0.1275_KREAL,0.01125_KREAL,0.85_KREAL]
    real(KREAL), parameter :: gnode(3)=[-0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL),0.0_KREAL, &
        0.5_KREAL*sqrt(3.0_KREAL/5.0_KREAL)]
    real(KREAL), save :: cmap(0:4,5,4)=0.0_KREAL
    logical, save :: maps_ready=.false.

    ! Owned by one field and one velocity node during a single compute_flux.
    ! The caller prepares it again for each new state; invalid preparation
    ! cannot leave a previously valid plane consumable.
    type :: FirstLayerPlane
        private
        logical :: valid=.false.
        integer :: nx=0,ny=0
        real(KREAL),allocatable :: xv(:,:,:),xd(:,:,:),rows(:,:,:),yv(:,:,:),yd(:,:,:)
    end type FirstLayerPlane

contains

    integer pure function ip(i,n)
        integer, intent(in) :: i,n
        ip=modulo(i-1,n)+1
    end function ip

    subroutine init_maps()
        integer :: offsets5(5),offsets3(3),c
        if (maps_ready) return
        offsets3=[-2,-1,0]; call build_map(offsets3,cmap(:,:,1))
        offsets3=[-1,0,1];  call build_map(offsets3,cmap(:,:,2))
        offsets3=[0,1,2];   call build_map(offsets3,cmap(:,:,3))
        offsets5=[-2,-1,0,1,2]; call build_map(offsets5,cmap(:,:,4))
        maps_ready=.true.
    end subroutine init_maps

    subroutine build_map(offsets,map)
        integer, intent(in) :: offsets(:)
        real(KREAL), intent(out) :: map(0:4,5)
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
                map(d,offsets(c)+3)=x(d+1)
            end do
        end do
    end subroutine build_map

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
        integer::c
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
        do c=1,4
            coef(:,c)=matmul(cmap(:,:,c),q)
        end do
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
    subroutine s2d1_weno1d_companion(q,companion,xi,value,derivative,companion_value,ok)
        real(KREAL),intent(in)::q(5),companion(5),xi
        real(KREAL),intent(out)::value,derivative,companion_value
        logical,intent(out)::ok
        real(KREAL)::fw(4),coef(0:4,4),ccoef(0:4,4),pow(0:4),dpow(0:4)
        real(KREAL)::cand(4),dcand(4),ccand(4)
        integer::c,d
        call weno_coefficients(q,fw,coef,ok);if(.not.ok)return
        do c=1,4;ccoef(:,c)=matmul(cmap(:,:,c),companion);enddo
        pow(0)=1.0_KREAL;dpow(0)=0.0_KREAL
        do d=1,4;pow(d)=pow(d-1)*xi;dpow(d)=real(d,KREAL)*xi**(d-1);enddo
        do c=1,4
            cand(c)=sum(coef(:,c)*pow);dcand(c)=sum(coef(:,c)*dpow)
            ccand(c)=sum(ccoef(:,c)*pow)
        enddo
        value=sum(fw*cand);derivative=sum(fw*dcand);companion_value=sum(fw*ccand)
        ok=all(ieee_is_finite([value,derivative,companion_value]))
    end subroutine s2d1_weno1d_companion

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

    subroutine companion_eval3(q,companion,value,derivative,companion_value,ok)
        real(KREAL),intent(in)::q(5),companion(5)
        real(KREAL),intent(out)::value(3),derivative(3),companion_value(3)
        logical,intent(out)::ok
        real(KREAL)::fw(4),coef(0:4,4),ccoef(0:4,4),pow(0:4),dpow(0:4)
        real(KREAL)::cand(4),dcand(4),ccand(4)
        integer::s,c,d
        call weno_coefficients(q,fw,coef,ok);if(.not.ok)return
        do c=1,4;ccoef(:,c)=matmul(cmap(:,:,c),companion);enddo
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

    subroutine s2d1_side3(field,axis,cell,tangent,side,value,dn,dt,ok)
        real(KREAL), intent(in) :: field(:,:),side
        integer, intent(in) :: axis,cell,tangent
        real(KREAL), intent(out) :: value(3),dn(3),dt(3)
        logical, intent(out) :: ok
        real(KREAL) :: normal_value(5),normal_deriv(5),q(5)
        integer :: r,s,nx,ny
        logical :: lok
        nx=size(field,1);ny=size(field,2);ok=.false.
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

    subroutine s2d1_face3(field,axis,face,tangent,left,right,dnl,dnr,dtl,dtr,ok)
        real(KREAL), intent(in) :: field(:,:)
        integer, intent(in) :: axis,face,tangent
        real(KREAL), intent(out) :: left(3),right(3),dnl(3),dnr(3),dtl(3),dtr(3)
        logical, intent(out) :: ok
        integer :: n,right_cell
        logical :: lok
        n=merge(size(field,1),size(field,2),axis==1);right_cell=ip(face+1,n)
        call s2d1_side3(field,axis,face,tangent,0.5_KREAL,left,dnl,dtl,lok);if(.not.lok)then;ok=.false.;return;endif
        call s2d1_side3(field,axis,right_cell,tangent,-0.5_KREAL,right,dnr,dtr,lok);ok=lok
    end subroutine s2d1_face3

    subroutine s2d1_cell9(field,i,j,value,ok)
        real(KREAL), intent(in) :: field(:,:)
        integer, intent(in) :: i,j
        real(KREAL), intent(out) :: value(3,3)
        logical, intent(out) :: ok
        real(KREAL) :: rows(5,3),q(5),dummy,fw(4),coef(0:4,4)
        integer :: r,s,a,b,nx,ny
        logical :: lok
        nx=size(field,1);ny=size(field,2);ok=.false.
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

    subroutine prepare_first_layer_plane(cache,field,ok)
        type(FirstLayerPlane),intent(inout)::cache
        real(KREAL),intent(in)::field(:,:)
        logical,intent(out)::ok
        real(KREAL)::q(5),fw(4),coef(0:4,4),dummy
        integer::nx,ny,i,j,s,a
        logical::lok
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
        endif
        do j=1,ny
            do i=1,nx
                do s=-2,2;q(s+3)=field(ip(i+s,nx),j);enddo
                call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
                call prepared_eval(fw,coef,-0.5_KREAL,cache%xv(i,j,1),cache%xd(i,j,1),lok);if(.not.lok)return
                call prepared_eval(fw,coef, 0.5_KREAL,cache%xv(i,j,2),cache%xd(i,j,2),lok);if(.not.lok)return
                do a=1,3
                    call prepared_eval(fw,coef,gnode(a),cache%rows(i,j,a),dummy,lok);if(.not.lok)return
                enddo
                do s=-2,2;q(s+3)=field(i,ip(j+s,ny));enddo
                call weno_coefficients(q,fw,coef,lok);if(.not.lok)return
                call prepared_eval(fw,coef,-0.5_KREAL,cache%yv(i,j,1),cache%yd(i,j,1),lok);if(.not.lok)return
                call prepared_eval(fw,coef, 0.5_KREAL,cache%yv(i,j,2),cache%yd(i,j,2),lok);if(.not.lok)return
            enddo
        enddo
        cache%valid=.true.;ok=.true.
    end subroutine prepare_first_layer_plane

    subroutine s2d1_cell_faces9_cached(cache,i,j,face,dn,dt,volume,ok)
        type(FirstLayerPlane),intent(in)::cache
        integer,intent(in)::i,j
        real(KREAL),intent(out)::face(4,3),dn(4,3),dt(4,3),volume(3,3)
        logical,intent(out)::ok
        real(KREAL)::xv(5,2),xd(5,2),yv(5,2),yd(5,2),rows(5,3)
        real(KREAL)::fw(4),coef(0:4,4),dummy
        integer::r,s,a,b
        logical::lok
        ok=.false.
        if(.not.cache%valid.or.i<1.or.i>cache%nx.or.j<1.or.j>cache%ny)return
        do r=-2,2
            xv(r+3,:)=cache%xv(i,ip(j+r,cache%ny),:)
            xd(r+3,:)=cache%xd(i,ip(j+r,cache%ny),:)
            rows(r+3,:)=cache%rows(i,ip(j+r,cache%ny),:)
        enddo
        do r=-2,2
            yv(r+3,:)=cache%yv(ip(i+r,cache%nx),j,:)
            yd(r+3,:)=cache%yd(ip(i+r,cache%nx),j,:)
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
    end subroutine s2d1_cell_faces9_cached

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

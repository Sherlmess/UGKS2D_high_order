! High-order dedicated slab metadata. LAP and legacy paths remain unchanged.
module mod_S2MPI1Layout
    use Precision, only: KREAL
    use flow_para, only: S2W1_wall_mode
    use mpi
    implicit none
    real(KREAL)::communication_seconds=0.0_KREAL,communication_clock=0.0_KREAL
    logical :: slab_enabled=.false.
    integer,parameter :: HO_HALO=3
    integer :: cart_halo=3
    logical :: cart_forward_halo=.false.
    integer :: slab_lo=1,slab_hi=0,slab_origin=1,slab_global_ny=0
    integer :: slab_minus=MPI_PROC_NULL,slab_plus=MPI_PROC_NULL
    logical :: cart_context=.false.
    integer :: failure_i=0,failure_j=0
    integer :: cart_xlo=4,cart_xhi=0,cart_ylo=4,cart_yhi=0
    integer :: cart_oi=1,cart_oj=1,cart_nx=0,cart_ny=0
    integer :: cart_xminus=MPI_PROC_NULL,cart_xplus=MPI_PROC_NULL
    integer :: cart_yminus=MPI_PROC_NULL,cart_yplus=MPI_PROC_NULL
    real(KREAL)::owner_reconstruction_s=0,owner_pack_s=0,owner_mpi_s=0,owner_unpack_s=0,owner_sync_s=0,owner_exchange_s=0
    integer(kind=8)::owner_bytes=0,owner_cells=0,owner_calls=0
    real(KREAL)::owner_buffer_setup_s=0
    integer(kind=8)::owner_buffer_allocations=0,owner_buffer_capacity_bytes=0
    ! Contents are overwritten on every use; capacity alone survives calls.
    real(KREAL),allocatable,private,save::owner_sm(:),owner_sp(:),owner_rm(:),owner_rp(:)
contains

    ! Research-only owner-computes exchange. Raw state halos are unchanged.
    ! Whole reconstructed boundary cells are sent to preserve every existing
    ! face consumer, including the redundant corner faces and limiter counts.
    logical function owner_enabled()
        character(len=16)::value
        integer::status
        call get_environment_variable('S2_OWNER_RECON',value,status=status)
        owner_enabled=cart_context.and.status==0.and.trim(value)=='1'
    end function
    subroutine owner_ready(ok)
        logical,intent(inout)::ok
        integer::bad,allbad,e
        real(KREAL)::t
        bad=merge(0,1,ok);t=MPI_Wtime()
        call MPI_Allreduce(bad,allbad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,e)
        owner_sync_s=owner_sync_s+MPI_Wtime()-t
        communication_seconds=communication_seconds+MPI_Wtime()-t
        if(e/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,781,e)
        ok=allbad==0
    end subroutine
    subroutine owner_exchange_full(a,b,c,d,f,g,counts)
        real(KREAL),intent(inout)::a(:,:,:,:,:),b(:,:,:,:,:),c(:,:,:,:,:),d(:,:,:,:,:),f(:,:,:,:,:),g(:,:,:,:,:)
        integer,intent(inout)::counts(:,:,:)
        real(KREAL),allocatable::sm(:,:),sp(:,:),rm(:,:),rp(:,:)
        integer::nv,nc,nx,ny,i,j,e,st(MPI_STATUS_SIZE)
        real(KREAL)::t,tt
        tt=MPI_Wtime();nv=size(a,3);nc=72*nv+2
        nx=cart_xhi-cart_xlo+1;ny=cart_yhi-cart_ylo+1
        allocate(sm(nc,ny),sp(nc,ny),rm(nc,ny),rp(nc,ny))
        t=MPI_Wtime()
        do j=cart_ylo,cart_yhi
            call pack(cart_xlo,j,sm(:,j-cart_ylo+1));call pack(cart_xhi,j,sp(:,j-cart_ylo+1))
        enddo
        owner_pack_s=owner_pack_s+MPI_Wtime()-t;t=MPI_Wtime()
        call MPI_Sendrecv(sm,size(sm),MPI_DOUBLE_PRECISION,cart_xminus,22400,rp,size(rp),MPI_DOUBLE_PRECISION,cart_xplus,22400,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        call MPI_Sendrecv(sp,size(sp),MPI_DOUBLE_PRECISION,cart_xplus,22401,rm,size(rm),MPI_DOUBLE_PRECISION,cart_xminus,22401,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        owner_mpi_s=owner_mpi_s+MPI_Wtime()-t;t=MPI_Wtime()
        do j=cart_ylo,cart_yhi
            if(cart_xminus/=MPI_PROC_NULL)call unpack(cart_xlo-1,j,rm(:,j-cart_ylo+1))
            if(cart_xplus/=MPI_PROC_NULL)call unpack(cart_xhi+1,j,rp(:,j-cart_ylo+1))
        enddo
        owner_unpack_s=owner_unpack_s+MPI_Wtime()-t
        owner_bytes=owner_bytes+int(nc,8)*ny*8*(merge(1,0,cart_xminus/=MPI_PROC_NULL)+merge(1,0,cart_xplus/=MPI_PROC_NULL))
        deallocate(sm,sp,rm,rp)
        ! X exchange completes first; Y includes both X neighbours for corners.
        allocate(sm(nc,nx+2),sp(nc,nx+2),rm(nc,nx+2),rp(nc,nx+2))
        t=MPI_Wtime()
        do i=cart_xlo-1,cart_xhi+1
            call pack(i,cart_ylo,sm(:,i-cart_xlo+2));call pack(i,cart_yhi,sp(:,i-cart_xlo+2))
        enddo
        owner_pack_s=owner_pack_s+MPI_Wtime()-t;t=MPI_Wtime()
        call MPI_Sendrecv(sm,size(sm),MPI_DOUBLE_PRECISION,cart_yminus,22402,rp,size(rp),MPI_DOUBLE_PRECISION,cart_yplus,22402,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        call MPI_Sendrecv(sp,size(sp),MPI_DOUBLE_PRECISION,cart_yplus,22403,rm,size(rm),MPI_DOUBLE_PRECISION,cart_yminus,22403,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        owner_mpi_s=owner_mpi_s+MPI_Wtime()-t;t=MPI_Wtime()
        do i=cart_xlo-1,cart_xhi+1
            if(cart_yminus/=MPI_PROC_NULL)call unpack(i,cart_ylo-1,rm(:,i-cart_xlo+2))
            if(cart_yplus/=MPI_PROC_NULL)call unpack(i,cart_yhi+1,rp(:,i-cart_xlo+2))
        enddo
        owner_unpack_s=owner_unpack_s+MPI_Wtime()-t
        owner_bytes=owner_bytes+int(nc,8)*(nx+2)*8*(merge(1,0,cart_yminus/=MPI_PROC_NULL)+merge(1,0,cart_yplus/=MPI_PROC_NULL))
        owner_exchange_s=owner_exchange_s+MPI_Wtime()-tt
        communication_seconds=communication_seconds+MPI_Wtime()-tt
        owner_calls=owner_calls+1
    contains
        subroutine check_mpi(code)
            integer,intent(in)::code
            integer::err
            if(code/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,782,err)
        end subroutine
        subroutine pack(ix,jy,z)
            integer,intent(in)::ix,jy
            real(KREAL),intent(out)::z(:)
            integer::n
            n=12*nv
            z(1:n)=reshape(a(:,:,:,ix,jy),[n]);z(n+1:2*n)=reshape(b(:,:,:,ix,jy),[n])
            z(2*n+1:3*n)=reshape(c(:,:,:,ix,jy),[n]);z(3*n+1:4*n)=reshape(d(:,:,:,ix,jy),[n])
            z(4*n+1:5*n)=reshape(f(:,:,:,ix,jy),[n]);z(5*n+1:6*n)=reshape(g(:,:,:,ix,jy),[n])
            z(6*n+1:6*n+2)=real(counts(:,ix,jy),KREAL)
        end subroutine
        subroutine unpack(ix,jy,z)
            integer,intent(in)::ix,jy
            real(KREAL),intent(in)::z(:)
            integer::n
            n=12*nv
            a(:,:,:,ix,jy)=reshape(z(1:n),[4,3,nv]);b(:,:,:,ix,jy)=reshape(z(n+1:2*n),[4,3,nv])
            c(:,:,:,ix,jy)=reshape(z(2*n+1:3*n),[4,3,nv]);d(:,:,:,ix,jy)=reshape(z(3*n+1:4*n),[4,3,nv])
            f(:,:,:,ix,jy)=reshape(z(4*n+1:5*n),[4,3,nv]);g(:,:,:,ix,jy)=reshape(z(5*n+1:6*n),[4,3,nv])
            counts(:,ix,jy)=nint(z(6*n+1:6*n+2))
        end subroutine
    end subroutine

    subroutine owner_exchange(a,b,c,d,f,g,counts)
        real(KREAL),intent(inout)::a(:,:,:,:,:),b(:,:,:,:,:),c(:,:,:,:,:),d(:,:,:,:,:),f(:,:,:,:,:),g(:,:,:,:,:)
        integer,intent(inout)::counts(:,:,:)
        character(len=16)::value
        integer::status
        call get_environment_variable('S2_OWNER_COMPACT',value,status=status)
        if(status==0.and.trim(value)=='1')then
            call owner_exchange_compact(a,b,c,d,f,g,counts)
        else
            call owner_exchange_full(a,b,c,d,f,g,counts)
        endif
    end subroutine

    subroutine owner_exchange_compact(a,b,c,d,f,g,counts)
        real(KREAL),intent(inout)::a(:,:,:,:,:),b(:,:,:,:,:),c(:,:,:,:,:),d(:,:,:,:,:),f(:,:,:,:,:),g(:,:,:,:,:)
        integer,intent(inout)::counts(:,:,:)
        integer::nv,nx,ny,i,j,e,st(MPI_STATUS_SIZE),xm,xp,ym,yp,nm,np,p,q,mask
        real(KREAL)::t,tt
        tt=MPI_Wtime();nv=size(a,3)
        nx=cart_xhi-cart_xlo+1;ny=cart_yhi-cart_ylo+1
        ! X-minus -> right halo: side 1. X-plus -> left halo: sides 2/3/4.
        xm=(18*nv+2)*ny;xp=2*ny
        do j=cart_ylo,cart_yhi
            xp=xp+18*nv*popcnt(left_mask(j))
        enddo
        ! Y-minus -> top halo: side 3 except right corner (counts only).
        ! Y-plus -> bottom halo: 1/2/4, left corner 2/4, right corner 1.
        ym=18*nv*(nx+1)+2*(nx+2);yp=18*nv*(3*nx+3)+2*(nx+2)
        t=MPI_Wtime();nm=max(xm,ym);np=max(xp,yp)
        if(allocated(owner_sm))then
            if(size(owner_sm)<nm.or.size(owner_sp)<np)deallocate(owner_sm,owner_sp,owner_rm,owner_rp)
        endif
        if(.not.allocated(owner_sm))then
            allocate(owner_sm(nm),owner_rp(nm),owner_sp(np),owner_rm(np))
            owner_buffer_allocations=owner_buffer_allocations+1
            owner_buffer_capacity_bytes=16_8*(int(nm,8)+int(np,8))
        endif
        owner_buffer_setup_s=owner_buffer_setup_s+MPI_Wtime()-t
        t=MPI_Wtime();p=0;q=0
        do j=cart_ylo,cart_yhi
            if(cart_xminus/=MPI_PROC_NULL)call pack_cell(cart_xlo,j,1,owner_sm,p)
            if(cart_xplus/=MPI_PROC_NULL)call pack_cell(cart_xhi,j,left_mask(j),owner_sp,q)
        enddo
        owner_pack_s=owner_pack_s+MPI_Wtime()-t;t=MPI_Wtime()
        call MPI_Sendrecv(owner_sm,merge(xm,0,cart_xminus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_xminus,22400, &
            owner_rp,merge(xm,0,cart_xplus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_xplus,22400,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        call MPI_Sendrecv(owner_sp,merge(xp,0,cart_xplus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_xplus,22401, &
            owner_rm,merge(xp,0,cart_xminus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_xminus,22401,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        owner_mpi_s=owner_mpi_s+MPI_Wtime()-t;t=MPI_Wtime();p=0;q=0
        do j=cart_ylo,cart_yhi
            if(cart_xminus/=MPI_PROC_NULL)call unpack_cell(cart_xlo-1,j,left_mask(j),owner_rm,p)
            if(cart_xplus/=MPI_PROC_NULL)call unpack_cell(cart_xhi+1,j,1,owner_rp,q)
        enddo
        owner_unpack_s=owner_unpack_s+MPI_Wtime()-t
        owner_bytes=owner_bytes+8_8*(int(merge(xm,0,cart_xminus/=MPI_PROC_NULL),8)+int(merge(xp,0,cart_xplus/=MPI_PROC_NULL),8))
        ! Completed X strips supply the two required Y corner subsets.
        t=MPI_Wtime();p=0;q=0
        do i=cart_xlo-1,cart_xhi+1
            if(cart_yminus/=MPI_PROC_NULL)call pack_cell(i,cart_ylo,top_mask(i),owner_sm,p)
            if(cart_yplus/=MPI_PROC_NULL)call pack_cell(i,cart_yhi,bottom_mask(i),owner_sp,q)
        enddo
        owner_pack_s=owner_pack_s+MPI_Wtime()-t;t=MPI_Wtime()
        call MPI_Sendrecv(owner_sm,merge(ym,0,cart_yminus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_yminus,22402, &
            owner_rp,merge(ym,0,cart_yplus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_yplus,22402,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        call MPI_Sendrecv(owner_sp,merge(yp,0,cart_yplus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_yplus,22403, &
            owner_rm,merge(yp,0,cart_yminus/=MPI_PROC_NULL),MPI_DOUBLE_PRECISION,cart_yminus,22403,MPI_COMM_WORLD,st,e)
        call check_mpi(e)
        owner_mpi_s=owner_mpi_s+MPI_Wtime()-t;t=MPI_Wtime();p=0;q=0
        do i=cart_xlo-1,cart_xhi+1
            if(cart_yminus/=MPI_PROC_NULL)call unpack_cell(i,cart_ylo-1,bottom_mask(i),owner_rm,p)
            if(cart_yplus/=MPI_PROC_NULL)call unpack_cell(i,cart_yhi+1,top_mask(i),owner_rp,q)
        enddo
        owner_unpack_s=owner_unpack_s+MPI_Wtime()-t
        owner_bytes=owner_bytes+8_8*(int(merge(ym,0,cart_yminus/=MPI_PROC_NULL),8)+int(merge(yp,0,cart_yplus/=MPI_PROC_NULL),8))
        owner_exchange_s=owner_exchange_s+MPI_Wtime()-tt
        communication_seconds=communication_seconds+MPI_Wtime()-tt
        owner_calls=owner_calls+1
    contains
        integer function left_mask(jy) result(m)
            integer,intent(in)::jy
            m=14
            if(lower_wall().and.jy==cart_ylo)m=ibclr(m,2)
            if(upper_wall().and.jy==cart_yhi)m=ibclr(m,3)
        end function
        integer function top_mask(ix) result(m)
            integer,intent(in)::ix
            m=4
            if(ix==cart_xhi+1)m=0
        end function
        integer function bottom_mask(ix) result(m)
            integer,intent(in)::ix
            m=11
            if(ix==cart_xlo-1)m=10
            if(ix==cart_xhi+1)m=1
        end function
        subroutine check_mpi(code)
            integer,intent(in)::code
            integer::err
            if(code/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,782,err)
        end subroutine
        subroutine pack_cell(ix,jy,m,z,pos)
            integer,intent(in)::ix,jy,m
            real(KREAL),intent(inout)::z(:)
            integer,intent(inout)::pos
            integer::s,k,node
            do s=1,4
                if(.not.btest(m,s-1))cycle
                do k=1,nv;do node=1,3
                    z(pos+1)=a(s,node,k,ix,jy);z(pos+2)=b(s,node,k,ix,jy)
                    z(pos+3)=c(s,node,k,ix,jy);z(pos+4)=d(s,node,k,ix,jy)
                    z(pos+5)=f(s,node,k,ix,jy);z(pos+6)=g(s,node,k,ix,jy)
                    pos=pos+6
                enddo;enddo
            enddo
            z(pos+1:pos+2)=real(counts(:,ix,jy),KREAL);pos=pos+2
        end subroutine
        subroutine unpack_cell(ix,jy,m,z,pos)
            integer,intent(in)::ix,jy,m
            real(KREAL),intent(in)::z(:)
            integer,intent(inout)::pos
            integer::s,k,node
            do s=1,4
                if(.not.btest(m,s-1))cycle
                do k=1,nv;do node=1,3
                    a(s,node,k,ix,jy)=z(pos+1);b(s,node,k,ix,jy)=z(pos+2)
                    c(s,node,k,ix,jy)=z(pos+3);d(s,node,k,ix,jy)=z(pos+4)
                    f(s,node,k,ix,jy)=z(pos+5);g(s,node,k,ix,jy)=z(pos+6)
                    pos=pos+6
                enddo;enddo
            enddo
            counts(:,ix,jy)=nint(z(pos+1:pos+2));pos=pos+2
        end subroutine
    end subroutine
    integer function xlo()
        xlo=1
        if(cart_context)xlo=cart_xlo
    end function
    integer function xhi(n)
        integer,intent(in)::n
        xhi=n
        if(cart_context)xhi=cart_xhi
    end function
    integer function xrecon_lo()
        xrecon_lo=xlo()
        if(cart_context.and..not.left_wall())xrecon_lo=xlo()-1
    end function
    integer function xrecon_hi(n)
        integer,intent(in)::n
        xrecon_hi=xhi(n)
        if(cart_context.and..not.right_wall())xrecon_hi=xhi(n)+1
    end function
    integer function yrecon_lo()
        yrecon_lo=active_lo()
        if(cart_context.and..not.lower_wall())yrecon_lo=active_lo()-1
    end function
    logical function cavity_mode()
        cavity_mode=S2W1_wall_mode==2
    end function
    integer function global_i(i)
        integer,intent(in)::i
        global_i=i
        if(cart_context)global_i=cart_oi+i-cart_xlo
    end function
    integer function global_nx(n)
        integer,intent(in)::n
        global_nx=n
        if(cart_context)global_nx=cart_nx
    end function
    logical function near_x_wall(i,n)
        integer,intent(in)::i,n
        near_x_wall=cavity_mode().and.(global_i(i)<=3.or.global_i(i)>=global_nx(n)-2)
    end function
    integer function wall_first_x(i,n)
        integer,intent(in)::i,n
        wall_first_x=merge(1,n-4,global_i(i)<=3)
        if(cart_context)wall_first_x=merge(cart_xlo-cart_oi+1,cart_xlo+cart_nx-cart_oi-4,global_i(i)<=3)
    end function
    logical function left_wall()
        left_wall=.false.
        if(cavity_mode())then
            left_wall=.true.
            if(cart_context)left_wall=cart_xminus==MPI_PROC_NULL
        endif
    end function
    logical function right_wall()
        right_wall=.false.
        if(cavity_mode())then
            right_wall=.true.
            if(cart_context)right_wall=cart_xplus==MPI_PROC_NULL
        endif
    end function
    integer function wall_first(j,n)
        integer,intent(in)::j,n
        wall_first=merge(1,n-4,global_j(j)<=3)
        if(cart_context)wall_first=merge(cart_ylo-cart_oj+1,cart_ylo+cart_ny-cart_oj-4,global_j(j)<=3)
    end function
    subroutine exchange_cart(field)
        real(KREAL),intent(inout)::field(:,:,:)
        real(KREAL),allocatable::sm(:,:,:),sp(:,:,:),rm(:,:,:),rp(:,:,:)
        integer::nc,nx,ny,e,st(MPI_STATUS_SIZE)
        if(cart_forward_halo)then
            call exchange_cart_forward(field);return
        endif
        nc=size(field,1);nx=cart_xhi-cart_xlo+1;ny=cart_yhi-cart_ylo+1
        allocate(sm(nc,HO_HALO,ny),sp(nc,HO_HALO,ny),rm(nc,HO_HALO,ny),rp(nc,HO_HALO,ny))
        sm=field(:,cart_xlo:cart_xlo+2,cart_ylo:cart_yhi)
        sp=field(:,cart_xhi-2:cart_xhi,cart_ylo:cart_yhi)
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(sm,size(sm),MPI_DOUBLE_PRECISION,cart_xminus,22000,rp,size(rp),MPI_DOUBLE_PRECISION,cart_xplus,22000,MPI_COMM_WORLD,st,e)
        call MPI_Sendrecv(sp,size(sp),MPI_DOUBLE_PRECISION,cart_xplus,22001,rm,size(rm),MPI_DOUBLE_PRECISION,cart_xminus,22001,MPI_COMM_WORLD,st,e)
        if(cart_xminus/=MPI_PROC_NULL)field(:,1:3,cart_ylo:cart_yhi)=rm
        if(cart_xplus/=MPI_PROC_NULL)field(:,cart_xhi+1:cart_xhi+3,cart_ylo:cart_yhi)=rp
        if(cart_xminus==MPI_PROC_NULL)field(:,1:3,cart_ylo:cart_yhi)=spread(field(:,cart_xlo,cart_ylo:cart_yhi),2,3)
        if(cart_xplus==MPI_PROC_NULL)field(:,cart_xhi+1:cart_xhi+3,cart_ylo:cart_yhi)=spread(field(:,cart_xhi,cart_ylo:cart_yhi),2,3)
        deallocate(sm,sp,rm,rp)
        ! Include the freshly exchanged x strips: this supplies all tensor corners.
        allocate(sm(nc,nx+6,3),sp(nc,nx+6,3),rm(nc,nx+6,3),rp(nc,nx+6,3))
        sm=field(:,:,cart_ylo:cart_ylo+2);sp=field(:,:,cart_yhi-2:cart_yhi)
        call MPI_Sendrecv(sm,size(sm),MPI_DOUBLE_PRECISION,cart_yminus,22002,rp,size(rp),MPI_DOUBLE_PRECISION,cart_yplus,22002,MPI_COMM_WORLD,st,e)
        call MPI_Sendrecv(sp,size(sp),MPI_DOUBLE_PRECISION,cart_yplus,22003,rm,size(rm),MPI_DOUBLE_PRECISION,cart_yminus,22003,MPI_COMM_WORLD,st,e)
        if(cart_yminus/=MPI_PROC_NULL)field(:,:,1:3)=rm
        if(cart_yplus/=MPI_PROC_NULL)field(:,:,cart_yhi+1:cart_yhi+3)=rp
        ! Physical wall ghost strips are never consumed by the one-sided stencil.
        if(cart_yminus==MPI_PROC_NULL)field(:,:,1:3)=spread(field(:,:,cart_ylo),3,3)
        if(cart_yplus==MPI_PROC_NULL)field(:,:,cart_yhi+1:cart_yhi+3)=spread(field(:,:,cart_yhi),3,3)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end subroutine

    ! Forward one layer per round. A thin neighbour may forward a layer it
    ! received in an earlier round; no global gather or extrapolation is used.
    subroutine exchange_cart_forward(field)
        real(KREAL),intent(inout)::field(:,:,:)
        real(KREAL),allocatable::sm(:,:),sp(:,:),rm(:,:),rp(:,:)
        integer::k,nc,ny,e,st(MPI_STATUS_SIZE)
        nc=size(field,1);ny=cart_yhi-cart_ylo+1
        allocate(sm(nc,ny),sp(nc,ny),rm(nc,ny),rp(nc,ny))
        communication_clock=MPI_Wtime()
        do k=1,cart_halo
            sm=field(:,cart_xlo+k-1,cart_ylo:cart_yhi)
            sp=field(:,cart_xhi-k+1,cart_ylo:cart_yhi)
            call MPI_Sendrecv(sm,size(sm),MPI_DOUBLE_PRECISION,cart_xminus,22200, &
                rp,size(rp),MPI_DOUBLE_PRECISION,cart_xplus,22200,MPI_COMM_WORLD,st,e)
            call MPI_Sendrecv(sp,size(sp),MPI_DOUBLE_PRECISION,cart_xplus,22201, &
                rm,size(rm),MPI_DOUBLE_PRECISION,cart_xminus,22201,MPI_COMM_WORLD,st,e)
            if(cart_xminus/=MPI_PROC_NULL)field(:,cart_xlo-k,cart_ylo:cart_yhi)=rm
            if(cart_xplus/=MPI_PROC_NULL)field(:,cart_xhi+k,cart_ylo:cart_yhi)=rp
            if(cart_xminus==MPI_PROC_NULL)field(:,cart_xlo-k,cart_ylo:cart_yhi)=field(:,cart_xlo,cart_ylo:cart_yhi)
            if(cart_xplus==MPI_PROC_NULL)field(:,cart_xhi+k,cart_ylo:cart_yhi)=field(:,cart_xhi,cart_ylo:cart_yhi)
        enddo
        deallocate(sm,sp,rm,rp)
        allocate(sm(nc,size(field,2)),sp(nc,size(field,2)),rm(nc,size(field,2)),rp(nc,size(field,2)))
        ! Include completed x halo, so diagonal stencils receive identical corners.
        do k=1,cart_halo
            sm=field(:,:,cart_ylo+k-1);sp=field(:,:,cart_yhi-k+1)
            call MPI_Sendrecv(sm,size(sm),MPI_DOUBLE_PRECISION,cart_yminus,22202, &
                rp,size(rp),MPI_DOUBLE_PRECISION,cart_yplus,22202,MPI_COMM_WORLD,st,e)
            call MPI_Sendrecv(sp,size(sp),MPI_DOUBLE_PRECISION,cart_yplus,22203, &
                rm,size(rm),MPI_DOUBLE_PRECISION,cart_yminus,22203,MPI_COMM_WORLD,st,e)
            if(cart_yminus/=MPI_PROC_NULL)then
                field(:,:,cart_ylo-k)=rm
            else
                field(:,:,cart_ylo-k)=field(:,:,cart_ylo)
            endif
            if(cart_yplus/=MPI_PROC_NULL)then
                field(:,:,cart_yhi+k)=rp
            else
                field(:,:,cart_yhi+k)=field(:,:,cart_yhi)
            endif
        enddo
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end subroutine exchange_cart_forward

    integer function active_lo()
        active_lo=1
        if(slab_enabled)active_lo=slab_lo
        if(cart_context)active_lo=cart_ylo
    end function
    integer function active_hi(n)
        integer,intent(in)::n
        active_hi=n
        if(slab_enabled)active_hi=slab_hi
        if(cart_context)active_hi=cart_yhi
    end function
    integer function recon_hi(n)
        integer,intent(in)::n
        recon_hi=n
        if(slab_enabled)recon_hi=min(slab_hi+1,n)
        if(cart_context)recon_hi=cart_yhi+merge(0,1,upper_wall())
    end function
    integer function global_j(j)
        integer,intent(in)::j
        global_j=j
        if(slab_enabled)global_j=slab_origin+j-1
        if(cart_context)global_j=cart_oj+j-cart_ylo
    end function
    integer function global_ny(n)
        integer,intent(in)::n
        global_ny=n
        if(slab_enabled)global_ny=slab_global_ny
        if(cart_context)global_ny=cart_ny
    end function
    logical function near_wall(j,n)
        integer,intent(in)::j,n
        near_wall=global_j(j)<=3.or.global_j(j)>=global_ny(n)-2
    end function
    logical function lower_wall()
        lower_wall=.true.
        if(slab_enabled)lower_wall=slab_minus==MPI_PROC_NULL
        if(cart_context)lower_wall=cart_yminus==MPI_PROC_NULL
    end function
    logical function upper_wall()
        upper_wall=.true.
        if(slab_enabled)upper_wall=slab_plus==MPI_PROC_NULL
        if(cart_context)upper_wall=cart_yplus==MPI_PROC_NULL
    end function
    ! Same two Sendrecv directions as original mod_MPI; x periodic corners
    ! are implicit because NP_X=1 and every exchanged row contains all x.
    subroutine exchange_w(w)
        real(KREAL),intent(inout)::w(:,:,:)
        real(KREAL),allocatable::sm(:,:,:),sp(:,:,:),rm(:,:,:),rp(:,:,:)
        integer::e,st(MPI_STATUS_SIZE),n
        if(cart_context)then
            call exchange_cart(w);return
        endif
        if(.not.slab_enabled)return
        sm=w(:,:,slab_lo:slab_lo+HO_HALO-1);sp=w(:,:,slab_hi-HO_HALO+1:slab_hi)
        allocate(rm, mold=sm);allocate(rp,mold=sp);n=size(sm)
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(sm,n,MPI_DOUBLE_PRECISION,slab_minus,21000,rp,n,MPI_DOUBLE_PRECISION, &
            slab_plus,21000,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(sp,n,MPI_DOUBLE_PRECISION,slab_plus,21001,rm,n,MPI_DOUBLE_PRECISION, &
            slab_minus,21001,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(slab_minus/=MPI_PROC_NULL)w(:,:,1:HO_HALO)=rm
        if(slab_plus/=MPI_PROC_NULL)w(:,:,slab_hi+1:slab_hi+HO_HALO)=rp
    end subroutine
    subroutine exchange_micro(h)
        real(KREAL),intent(inout)::h(:,:,:)
        real(KREAL),allocatable::sm(:,:,:),sp(:,:,:),rm(:,:,:),rp(:,:,:)
        integer::e,st(MPI_STATUS_SIZE),n
        if(.not.slab_enabled)return
        sm=h(:,slab_lo:slab_lo+HO_HALO-1,:);sp=h(:,slab_hi-HO_HALO+1:slab_hi,:)
        allocate(rm,mold=sm);allocate(rp,mold=sp);n=size(sm)
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(sm,n,MPI_DOUBLE_PRECISION,slab_minus,21002,rp,n,MPI_DOUBLE_PRECISION, &
            slab_plus,21002,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(sp,n,MPI_DOUBLE_PRECISION,slab_plus,21003,rm,n,MPI_DOUBLE_PRECISION, &
            slab_minus,21003,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(slab_minus/=MPI_PROC_NULL)h(:,1:HO_HALO,:)=rm
        if(slab_plus/=MPI_PROC_NULL)h(:,slab_hi+1:slab_hi+HO_HALO,:)=rp
    end subroutine
end module

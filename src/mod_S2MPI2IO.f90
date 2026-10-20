! MPI-2 I/O and measurement only. All physical kernels remain MPI-1.
module mod_S2MPI2IO
    use iso_fortran_env, only: int32,int64
    use iso_c_binding, only: c_char,c_null_char,c_int,c_int64_t
    use Precision, only: KREAL,KINT
    use mpi
    use flow_para
    use flow_data, only: gridx_all,gridy_all
    use mpi_para
    use control_para, only: dt_global,TimeIntegration_method,TimeStep_method
    use mod_S2MPI1Layout
    use mod_S2D1, only: write_state,heat
    use mod_CorrectedWF2DEquilibrium, only: analytic_maxwellian_2d
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer :: checkpoint_stride=500,field_stride=0,diagnostic_stride=1,warmup_steps=0,audit_stride=0,stop_after=-1
    logical :: dump_initial=.false.,audit_faces=.false.,write_fields=.true.
    character(len=1024)::restart_path='',stop_file='STOP_REQUEST'
    real(KREAL)::io_seconds=0,io_comm_seconds=0
    integer :: checkpoint_gathers=0
    character(len=8),parameter::magic='S2MPI2V1',tail='COMPLETE'
    interface
        integer(c_int) function c_rename(old,new) bind(C,name='rename')
            import c_int,c_char
            character(kind=c_char),intent(in)::old(*),new(*)
        end function
        integer(c_int) function c_getpid() bind(C,name='getpid')
            import c_int
        end function
        integer(c_int64_t) function peak_rss() bind(C,name='mpi2_peak_rss')
            import c_int64_t
        end function
        integer(c_int64_t) function current_rss() bind(C,name='mpi2_current_rss')
            import c_int64_t
        end function
    end interface
contains
    subroutine io_require(ok,message)
        logical,intent(in)::ok
        character(len=*),intent(in)::message
        logical::allok
        integer::e,fu
        character(len=64)::filename
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(ok,allok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(allok)return
        write(filename,'("input_rejected_rank",I2.2,".txt")')myID
        open(newunit=fu,file=trim(filename),status='replace');write(fu,'(A)')message;close(fu)
        if(myID==0)write(*,'(A)')'MPI2 INPUT/IO REJECTED BEFORE PUBLISH: '//message
        communication_clock=MPI_Wtime()
        call MPI_Abort(MPI_COMM_WORLD,92,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end subroutine
    subroutine configure_io()
        integer::fu,ios,e,options(7),flags(3)
        logical::exists,ok
        namelist /mpi2_io/ checkpoint_stride,field_stride,diagnostic_stride,warmup_steps,audit_stride,stop_after, &
            dump_initial,audit_faces,write_fields,restart_path,stop_file
        ok=.true.
        if(myID==0)then
            inquire(file='MPI2_IO.nml',exist=exists)
            if(exists)then
                open(newunit=fu,file='MPI2_IO.nml',status='old',iostat=ios)
                ok=ios==0
                if(ok)then;read(fu,nml=mpi2_io,iostat=ios);close(fu);ok=ios==0;endif
            endif
            ok=ok.and.checkpoint_stride>=0.and.field_stride>=0.and.diagnostic_stride>=1.and.warmup_steps>=0.and.audit_stride>=0
        endif
        call io_require(ok,'invalid MPI2_IO.nml')
        options=[checkpoint_stride,field_stride,diagnostic_stride,warmup_steps,audit_stride,stop_after,0]
        flags=[merge(1,0,dump_initial),merge(1,0,audit_faces),merge(1,0,write_fields)]
        communication_clock=MPI_Wtime()
        call MPI_Bcast(options,7,MPI_INTEGER,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(flags,3,MPI_INTEGER,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(restart_path,len(restart_path),MPI_CHARACTER,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(stop_file,len(stop_file),MPI_CHARACTER,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        checkpoint_stride=options(1);field_stride=options(2);diagnostic_stride=options(3);warmup_steps=options(4)
        audit_stride=options(5);stop_after=options(6);dump_initial=flags(1)==1;audit_faces=flags(2)==1;write_fields=flags(3)==1
    end subroutine
    function physical_parameters() result(p)
        real(KREAL)::p(16)
        p=[real(ck,KREAL),Pr,S2_mu_ref,omega,Tcm,Tinf,Tw,Ma,Re,AOA, &
            S2W1_wall_ux(1),S2W1_wall_ux(2),S2W1_wall_theta(1),S2W1_wall_theta(2),ACC,dt_global]
    end function
    function option_identity() result(p)
        integer(int32)::p(12)
        p=int([BC_Type,FluxWall_type,TimeStep_method,S2PhysicalGas,S2W1_wall_mode,S2W2R_shared_predictor, &
            S2W1_time_centered,S2W1D_face_time_consistent,HO_HALO],int32)
    end function
    subroutine distribute(gw,gh,gb,w,h,b)
        real(KREAL),intent(in)::gw(:,:,:),gh(:,:,:),gb(:,:,:)
        real(KREAL),intent(out)::w(:,:,:),h(:,:,:),b(:,:,:)
        real(KREAL),allocatable::sw(:,:,:),sh(:,:,:),sb(:,:,:)
        integer::p,o,e,st(MPI_STATUS_SIZE),nl
        real(KREAL)::start
        start=MPI_Wtime();nl=slab_hi-slab_lo+1
        allocate(sw(4,Nx_all,nl),sh(Nx_all,nl,size(h,3)),sb(Nx_all,nl,size(h,3)))
        w=0;h=0;b=0
        if(myID==0)then
            do p=0,NP_XY-1
                o=Block_XY(p)%Origin_J;sw=gw(:,:,o:o+nl-1);sh=gh(:,o:o+nl-1,:);sb=gb(:,o:o+nl-1,:)
                if(p==0)then
                    w(:,:,slab_lo:slab_hi)=sw;h(:,slab_lo:slab_hi,:)=sh;b(:,slab_lo:slab_hi,:)=sb
                else
        communication_clock=MPI_Wtime()
                    call MPI_Send(sw,size(sw),MPI_DOUBLE_PRECISION,p,24000,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
                    call MPI_Send(sh,size(sh),MPI_DOUBLE_PRECISION,p,24001,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
                    call MPI_Send(sb,size(sb),MPI_DOUBLE_PRECISION,p,24002,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                endif
            enddo
        else
        communication_clock=MPI_Wtime()
            call MPI_Recv(sw,size(sw),MPI_DOUBLE_PRECISION,0,24000,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Recv(sh,size(sh),MPI_DOUBLE_PRECISION,0,24001,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Recv(sb,size(sb),MPI_DOUBLE_PRECISION,0,24002,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            w(:,:,slab_lo:slab_hi)=sw;h(:,slab_lo:slab_hi,:)=sh;b(:,slab_lo:slab_hi,:)=sb
        endif
        io_comm_seconds=io_comm_seconds+MPI_Wtime()-start
    end subroutine
    subroutine gather(w,h,b,gw,gh,gb)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:)
        real(KREAL),allocatable,intent(out)::gw(:,:,:),gh(:,:,:),gb(:,:,:)
        real(KREAL),allocatable::sw(:,:,:),sh(:,:,:),sb(:,:,:),packed(:)
        integer::nl,p,o,k,e
        real(KREAL)::start
        start=MPI_Wtime();nl=slab_hi-slab_lo+1
        sw=w(:,:,slab_lo:slab_hi);sh=h(:,slab_lo:slab_hi,:);sb=b(:,slab_lo:slab_hi,:)
        if(myID==0)then
            allocate(gw(4,Nx_all,Ny_all),gh(Nx_all,Ny_all,size(h,3)),gb(Nx_all,Ny_all,size(h,3)),packed(size(sh)*NP_XY))
        else
            allocate(gw(1,1,1),gh(1,1,1),gb(1,1,1),packed(1))
        endif
        communication_clock=MPI_Wtime()
        call MPI_Gather(sw,size(sw),MPI_DOUBLE_PRECISION,gw,size(sw),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Gather(sh,size(sh),MPI_DOUBLE_PRECISION,packed,size(sh),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(myID==0)then
            do p=0,NP_XY-1;o=Block_XY(p)%Origin_J;k=p*size(sh);gh(:,o:o+nl-1,:)=reshape(packed(k+1:k+size(sh)),shape(sh));enddo
        endif
        communication_clock=MPI_Wtime()
        call MPI_Gather(sb,size(sb),MPI_DOUBLE_PRECISION,packed,size(sb),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(myID==0)then
            do p=0,NP_XY-1;o=Block_XY(p)%Origin_J;k=p*size(sb);gb(:,o:o+nl-1,:)=reshape(packed(k+1:k+size(sb)),shape(sb));enddo
        endif
        io_comm_seconds=io_comm_seconds+MPI_Wtime()-start;checkpoint_gathers=checkpoint_gathers+1
    end subroutine
    subroutine load_input(w,h,b,u,v,wt,t,accepted,total0,balance,cumbottom,cumtop,resumed)
        real(KREAL),intent(out)::w(:,:,:),h(:,:,:),b(:,:,:),t,total0(4),balance(4),cumbottom(4),cumtop(4)
        real(KREAL),intent(in)::u(:),v(:),wt(:)
        integer(int64),intent(out)::accepted
        logical,intent(out)::resumed
        real(KREAL),allocatable::gw(:,:,:),gh(:,:,:),gb(:,:,:),uv(:),vv(:),wv(:),gx(:,:),gy(:,:)
        real(KREAL)::params(16),rdt
        integer(int32)::header(8),opts(12),rx,ry,rv
        character(len=8)::head,endmark
        character(len=32)::method,table
        character(len=1024)::path
        character(len=1)::extra
        integer::fu,ios,e,i
        logical::ok,is_checkpoint,file_open
        t=0;accepted=0;total0=0;balance=0;cumbottom=0;cumtop=0;resumed=.false.;ok=.true.;file_open=.false.
        allocate(gw(1,1,1),gh(1,1,1),gb(1,1,1))
        if(myID==0)then
            path=restart_path;is_checkpoint=len_trim(path)>0
            if(.not.is_checkpoint)call get_environment_variable('S2DO1_RESTART_PATH',path)
            open(newunit=fu,file=trim(path),status='old',form='unformatted',access='stream',iostat=ios);ok=ios==0;file_open=ok
            header=0;opts=0;params=0;rx=0;ry=0;rv=0;head='';endmark='';method='';table=''
            if(ok.and.is_checkpoint)then
                read(fu,iostat=ios)head,header,accepted,t,rdt,params,opts,method,table,total0,balance,cumbottom,cumtop
                ok=ios==0
                if(ok)ok=head==magic.and.all(header==int([1,16909060,storage_size(0.0_KREAL),Nx_all,Ny_all,size(u),ck,flag_vislaw],int32)) &
                    .and.rdt==dt_global.and.all(params==physical_parameters()).and.all(opts==option_identity()) &
                    .and.trim(method)==trim(TimeIntegration_method).and.trim(table)=='half32/scale2/1024' &
                    .and.accepted>=0.and.ieee_is_finite(t).and.t>=0.and.all(ieee_is_finite(total0)) &
                    .and.all(ieee_is_finite(balance)).and.all(ieee_is_finite(cumbottom)).and.all(ieee_is_finite(cumtop))
                if(ok)then
                    allocate(gx(Nx_all+1,Ny_all+1),gy(Nx_all+1,Ny_all+1));read(fu,iostat=ios)gx,gy
                    ok=ios==0.and.all(gx==gridx_all).and.all(gy==gridy_all)
                endif
                resumed=ok
            elseif(ok)then
                read(fu,iostat=ios)rx,ry,rv,t
                ok=ios==0.and.rx==Nx_all.and.ry==Ny_all.and.rv==size(u).and.t==0.0_KREAL
            endif
            if(ok)then
                deallocate(gw,gh,gb);allocate(gw(4,Nx_all,Ny_all),gh(Nx_all,Ny_all,size(u)),gb(Nx_all,Ny_all,size(u)), &
                    uv(size(u)),vv(size(u)),wv(size(u)))
                read(fu,iostat=ios)gw,gh,gb,uv,vv,wv
                ok=ios==0.and.all(uv==u).and.all(vv==v).and.all(wv==wt).and.all(ieee_is_finite(gw)) &
                    .and.all(ieee_is_finite(gh)).and.all(ieee_is_finite(gb))
                if(ok.and.is_checkpoint)then
                    read(fu,iostat=ios)endmark;ok=ios==0.and.endmark==tail
                endif
                if(ok)then;read(fu,iostat=ios)extra;ok=ios<0;endif
            endif
            if(file_open)close(fu)
        endif
        call io_require(ok,'checkpoint/initial file incomplete or physical, grid, table, dt, method identity mismatch')
        communication_clock=MPI_Wtime()
        call MPI_Bcast(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(accepted,1,MPI_INTEGER8,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(total0,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(balance,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(cumbottom,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(cumtop,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Bcast(resumed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        call distribute(gw,gh,gb,w,h,b)
        call exchange_w(w);call exchange_micro(h);call exchange_micro(b)
    end subroutine
    subroutine write_field(path,w,h,b,u,v,wt,t)
        character(len=*),intent(in)::path
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),wt(:),t
        real(KREAL)::rho,ux,vy,pres,theta,q(2),eqh(size(u)),eqb(size(u))
        integer::i,j,fu
        logical::ok
        open(newunit=fu,file=path,status='replace');write(fu,'(A)')'x,y,rho,U,V,p,T,qx,qy,time'
        do j=1,size(h,2);do i=1,size(h,1)
            rho=w(1,i,j);ux=w(2,i,j)/rho;vy=w(3,i,j)/rho
            pres=2.0_KREAL/real(ck+2,KREAL)*(w(4,i,j)-0.5_KREAL*rho*(ux*ux+vy*vy));theta=pres/rho
            call analytic_maxwellian_2d(w(:,i,j),ck,u,v,eqh,eqb,ok)
            call heat(w(:,i,j),h(i,j,:)-eqh,b(i,j,:)-eqb,u,v,wt,q,ok)
            write(fu,'(10(ES24.16,:,","))')(i-0.5_KREAL)*2*pi/Nx_all,(j-0.5_KREAL)/Ny_all,rho,ux,vy,pres,theta,q,t
        enddo;enddo
        close(fu)
    end subroutine
    subroutine parallel_field(prefix,w,h,b,u,v,wt,t)
        character(len=*),intent(in)::prefix
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),wt(:),t
        real(KREAL)::rho,ux,vy,pres,theta,q(2),eqh(size(u)),eqb(size(u))
        integer::fu,i,j
        logical::ok
        character(len=128)::path
        write(path,'(A,"_rank",I2.2,".csv")')trim(prefix),myID
        open(newunit=fu,file=trim(path),status='replace')
        write(fu,'(A)')'x,y,rho,U,V,p,T,qx,qy,global_j,owned'
        do j=1,size(h,2);do i=1,size(h,1)
            rho=w(1,i,j);ux=w(2,i,j)/rho;vy=w(3,i,j)/rho
            pres=2.0_KREAL/real(ck+2,KREAL)*(w(4,i,j)-0.5_KREAL*rho*(ux*ux+vy*vy));theta=pres/rho
            call analytic_maxwellian_2d(w(:,i,j),ck,u,v,eqh,eqb,ok)
            call heat(w(:,i,j),h(i,j,:)-eqh,b(i,j,:)-eqb,u,v,wt,q,ok)
            write(fu,'(9(ES24.16,","),I0,",",I0)')(i-0.5_KREAL)*2*pi/Nx_all, &
                (global_j(j)-0.5_KREAL)/Ny_all,rho,ux,vy,pres,theta,q,global_j(j),merge(1,0,j>=slab_lo.and.j<=slab_hi)
        enddo;enddo
        close(fu)
    end subroutine
    subroutine publish_checkpoint(w,h,b,u,v,wt,t,accepted,total0,balance,cumbottom,cumtop,final)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),wt(:),t,total0(4),balance(4),cumbottom(4),cumtop(4)
        integer(int64),intent(in)::accepted
        logical,intent(in)::final
        real(KREAL),allocatable::gw(:,:,:),gh(:,:,:),gb(:,:,:)
        integer::fu,ios,e,ren
        real(KREAL)::start
        logical::ok
        start=MPI_Wtime();call gather(w,h,b,gw,gh,gb);ok=.true.
        if(myID==0)then
            open(newunit=fu,file='checkpoint.bin.tmp',status='replace',form='unformatted',access='stream',iostat=ios);ok=ios==0
            if(ok)then
                write(fu,iostat=ios)magic,int([1,16909060,storage_size(0.0_KREAL),Nx_all,Ny_all,size(u),ck,flag_vislaw],int32), &
                    accepted,t,dt_global,physical_parameters(),option_identity()
                ok=ios==0
                if(ok)then
                    block
                        character(len=32)::method,table
                        method=TimeIntegration_method;table='half32/scale2/1024'
                        write(fu,iostat=ios)method,table,total0,balance,cumbottom,cumtop,gridx_all,gridy_all,gw,gh,gb,u,v,wt,tail
                        ok=ios==0
                    end block
                endif
                flush(fu,iostat=ios);ok=ok.and.ios==0;close(fu,iostat=ios);ok=ok.and.ios==0
            endif
            if(ok)then
                ren=c_rename('checkpoint.bin.tmp'//c_null_char,'checkpoint.bin'//c_null_char);ok=ren==0
            endif
            if(ok.and.final)then
                call write_state('native_state_final.bin.tmp',gw,gh,gb,u,v,wt,t)
                ren=c_rename('native_state_final.bin.tmp'//c_null_char,'native_state_final.bin'//c_null_char);ok=ren==0
            endif
            if(ok.and.write_fields)call write_field('native_field.csv',gw,gh,gb,u,v,wt,t)
        endif
        call io_require(ok,'atomic checkpoint publication failed; prior complete checkpoint remains recoverable')
        io_seconds=io_seconds+MPI_Wtime()-start
    end subroutine
    subroutine initial_output(w,h,b,u,v,wt,t)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),wt(:),t
        real(KREAL),allocatable::gw(:,:,:),gh(:,:,:),gb(:,:,:)
        real(KREAL)::start
        if(.not.dump_initial)return
        start=MPI_Wtime();call gather(w,h,b,gw,gh,gb)
        if(myID==0)call write_state('native_state_initial.bin',gw,gh,gb,u,v,wt,t)
        io_seconds=io_seconds+MPI_Wtime()-start
    end subroutine
    logical function stop_requested(segment_steps)
        integer,intent(in)::segment_steps
        integer::e
        stop_requested=.false.
        if(myID==0)then
            inquire(file=trim(stop_file),exist=stop_requested)
            if(stop_after>=0.and.segment_steps>=stop_after)stop_requested=.true.
        endif
        communication_clock=MPI_Wtime()
        call MPI_Bcast(stop_requested,1,MPI_LOGICAL,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end function
    subroutine audit_cache(step,w,h,b,mh,mb,xh,xb,yh,yb,eh,eb,kap)
        integer(int64),intent(in)::step
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),mh(:,:,:),mb(:,:,:),xh(:,:,:),xb(:,:,:),yh(:,:,:),yb(:,:,:),eh(:,:,:),eb(:,:,:),kap(:,:,:,:)
        character(len=128)::path
        integer::fu
        if(audit_stride==0)return
        if(mod(step,int(audit_stride,int64))/=0)return
        write(path,'("cache_step",I10.10,"_rank",I2.2,".bin")')step,myID
        open(newunit=fu,file=trim(path),status='replace',form='unformatted',access='stream')
        write(fu)int([Nx,Ny,size(h,3),Block_XY(myID)%Origin_J],int32)
        write(fu)w(:,:,slab_lo:slab_hi),h(:,slab_lo:slab_hi,:),b(:,slab_lo:slab_hi,:), &
            mh(:,slab_lo:slab_hi,:),mb(:,slab_lo:slab_hi,:),xh(:,slab_lo:slab_hi,:),xb(:,slab_lo:slab_hi,:), &
            yh(:,slab_lo:slab_hi,:),yb(:,slab_lo:slab_hi,:),eh(:,slab_lo:slab_hi,:),eb(:,slab_lo:slab_hi,:),kap(:,:, :,slab_lo:slab_hi)
        close(fu)
        write(path,'("halo_step",I10.10,"_rank",I2.2,".bin")')step,myID
        open(newunit=fu,file=trim(path),status='replace',form='unformatted',access='stream')
        write(fu)int([Nx,size(h,2),size(h,3),slab_origin],int32);write(fu)w,h,b;close(fu)
    end subroutine
end module

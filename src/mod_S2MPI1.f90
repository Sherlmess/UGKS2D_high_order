! Bounded native physical-space MPI driver for the frozen S2CN2 wall kernels.
module mod_S2MPI1
    use iso_fortran_env, only: int64
    use mod_S2MPI2IO
    use Precision, only: KREAL,KINT
    use mpi
    use mpi_para
    use flow_para
    use control_para, only: MAX_ITER,dt_global
    use InitializeFlow, only: InitAllocation,InitMesh
    use flow_data, only: ctr,vface,hface
    use mpi_wrapper, only: partion_mpi
    use mod_S2MPI1Layout
    use mod_S2D1, only: init_velocity,compute_flux,build_targets,advance_step,write_state,heat,w2_wall_calls
    use mod_CorrectedWF2DEquilibrium, only: analytic_maxwellian_2d
    use mod_S2PhysicalGas, only: gas_tau_calls,gas_tau_min,gas_tau_max
    use mod_VelocityTables, only: VelocityCount
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: RunS2MPI1
    real(KREAL),parameter::gate=2e-12_KREAL
contains
    logical function collective(ok) result(allok)
        logical,intent(in)::ok
        integer::e
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(ok,allok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end function
    subroutine require(ok,message)
        logical,intent(in)::ok
        character(len=*),intent(in)::message
        integer::e
        if(.not.collective(ok))then
            if(myID==0)write(*,'(A)')'MPI1 FAILURE: '//message
        communication_clock=MPI_Wtime()
            call MPI_Abort(MPI_COMM_WORLD,91,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        endif
    end subroutine
    subroutine sum4(a)
        real(KREAL),intent(inout)::a(4)
        real(KREAL)::tmp(4)
        integer::e
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(a,tmp,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,e);a=tmp
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end subroutine
    ! Lower-y owner sends last physical face. Upper-y rank uses precisely that payload.
    subroutine face_exchange(fw,fh,fb)
        real(KREAL),intent(inout)::fw(:,:,:),fh(:,:,:),fb(:,:,:)
        real(KREAL),allocatable::send(:),recv(:),echo(:)
        integer::n,e,st(MPI_STATUS_SIZE),a,c,nx,nv
        nx=size(fh,1);nv=size(fh,3);a=4*nx;c=nx*nv;n=a+2*c
        send=[reshape(fw(:,:,slab_hi),[a]),reshape(fh(:,slab_hi,:),[c]),reshape(fb(:,slab_hi,:),[c])]
        allocate(recv(n),echo(n));recv=0;echo=0
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(send,n,MPI_DOUBLE_PRECISION,slab_plus,22000,recv,n,MPI_DOUBLE_PRECISION, &
            slab_minus,22000,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(slab_minus/=MPI_PROC_NULL)then
            fw(:,:,slab_lo-1)=reshape(recv(:a),[4,nx])
            fh(:,slab_lo-1,:)=reshape(recv(a+1:a+c),[nx,nv])
            fb(:,slab_lo-1,:)=reshape(recv(a+c+1:),[nx,nv])
            ! Echo from the arrays consumed by the divergence, not merely the receive buffer.
            recv=[reshape(fw(:,:,slab_lo-1),[a]),reshape(fh(:,slab_lo-1,:),[c]), &
                reshape(fb(:,slab_lo-1,:),[c])]
        endif
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(recv,n,MPI_DOUBLE_PRECISION,slab_minus,22001,echo,n,MPI_DOUBLE_PRECISION, &
            slab_plus,22001,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        call require(slab_plus==MPI_PROC_NULL.or.all(echo==send),'shared W/h/b face payload mismatch')
    end subroutine
    subroutine save_state(path,w,h,b,u,v,wt,t)
        character(len=*),intent(in)::path
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),wt(:),t
        real(KREAL),allocatable::gw(:,:,:),gh(:,:,:),gb(:,:,:),sw(:,:,:),sh(:,:,:),sb(:,:,:),rh(:,:,:)
        integer::e,p,o,nx,nv,nloc,k,fu,i,j
        real(KREAL)::rho,ux,vy,theta,pressure,q(2),x,y,eqh(size(u)),eqb(size(u))
        logical::ok
        nx=size(h,1);nv=size(h,3);nloc=slab_hi-slab_lo+1
        sw=w(:,:,slab_lo:slab_hi);sh=h(:,slab_lo:slab_hi,:);sb=b(:,slab_lo:slab_hi,:)
        if(myID==0)then
            allocate(gw(4,nx,Ny_all),gh(nx,Ny_all,nv),gb(nx,Ny_all,nv),rh(nx,nloc,nv))
        else
            allocate(gw(1,1,1),gh(1,1,1),gb(1,1,1))
        endif
        communication_clock=MPI_Wtime()
        call MPI_Gather(sw,size(sw),MPI_DOUBLE_PRECISION,gw,size(sw),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        ! h/b gather rank-major slabs then root unpack into Fortran x/y/velocity ordering.
        block
            real(KREAL),allocatable::packed(:)
            allocate(packed(max(1,merge(size(sh)*NP_XY,1,myID==0))))
        communication_clock=MPI_Wtime()
            call MPI_Gather(sh,size(sh),MPI_DOUBLE_PRECISION,packed,size(sh),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            if(myID==0)then
                do p=0,NP_XY-1
                    o=Block_XY(p)%Origin_J;k=p*size(sh)
                    gh(:,o:o+nloc-1,:)=reshape(packed(k+1:k+size(sh)),shape(sh))
                enddo
            endif
        communication_clock=MPI_Wtime()
            call MPI_Gather(sb,size(sb),MPI_DOUBLE_PRECISION,packed,size(sb),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            if(myID==0)then
                do p=0,NP_XY-1
                    o=Block_XY(p)%Origin_J;k=p*size(sb)
                    gb(:,o:o+nloc-1,:)=reshape(packed(k+1:k+size(sb)),shape(sb))
                enddo
            endif
        end block
        if(myID/=0)return
        call write_state(path,gw,gh,gb,u,v,wt,t)
        if(path/='native_state_final.bin')return
        open(newunit=fu,file='native_field.csv',status='replace')
        write(fu,'(A)')'x,y,rho,U,V,p,T,qx,qy,time'
        do j=1,Ny_all;do i=1,nx
            rho=gw(1,i,j);ux=gw(2,i,j)/rho;vy=gw(3,i,j)/rho
            pressure=(2.0_KREAL/real(ck+2,KREAL))*(gw(4,i,j)-0.5_KREAL*rho*(ux*ux+vy*vy));theta=pressure/rho
            call analytic_maxwellian_2d(gw(:,i,j),ck,u,v,eqh,eqb,ok)
            call heat(gw(:,i,j),gh(i,j,:)-eqh,gb(i,j,:)-eqb,u,v,wt,q,ok)
            x=(real(i,KREAL)-0.5_KREAL)*2*pi/nx;y=(real(j,KREAL)-0.5_KREAL)/Ny_all
            write(fu,'(10(ES24.16,:,","))')x,y,rho,ux,vy,pressure,theta,q,t
        enddo;enddo
        close(fu)
    end subroutine
    subroutine load_state(w,h,b,u,v,wt,t)
        real(KREAL),intent(out)::w(:,:,:),h(:,:,:),b(:,:,:),t
        real(KREAL),intent(in)::u(:),v(:),wt(:)
        real(KREAL),allocatable::gw(:,:,:),gh(:,:,:),gb(:,:,:),uv(:),vv(:),wv(:)
        real(KREAL),allocatable::sw(:,:,:),sh(:,:,:),sb(:,:,:)
        character(len=1024)::path
        integer::fu,e,ios,rx,ry,rv,p,o,nloc,st(MPI_STATUS_SIZE)
        logical::ok
        w=0;h=0;b=0;t=0;ok=.true.;nloc=slab_hi-slab_lo+1
        if(myID==0)then
            call get_environment_variable('S2DO1_RESTART_PATH',path)
            open(newunit=fu,file=trim(path),status='old',form='unformatted',access='stream',iostat=ios)
            ok=ios==0
            if(ok)then
                read(fu,iostat=ios)rx,ry,rv,t
                ok=ios==0.and.rx==Nx_all.and.ry==Ny_all.and.rv==size(u)
            endif
        endif
        call require(ok,'complete initial state missing or wrong dimensions')
        if(myID==0)then
            allocate(gw(4,Nx_all,Ny_all),gh(Nx_all,Ny_all,size(u)),gb(Nx_all,Ny_all,size(u)), &
                uv(size(u)),vv(size(u)),wv(size(u)))
            read(fu,iostat=ios)gw,gh,gb,uv,vv,wv;close(fu)
            ok=ios==0.and.all(uv==u).and.all(vv==v).and.all(wv==wt).and.all(ieee_is_finite(gw)) &
                .and.all(ieee_is_finite(gh)).and.all(ieee_is_finite(gb))
        endif
        call require(ok,'initial payload/table mismatch')
        allocate(sw(4,Nx_all,nloc),sh(Nx_all,nloc,size(u)),sb(Nx_all,nloc,size(u)))
        if(myID==0)then
            do p=0,NP_XY-1
                o=Block_XY(p)%Origin_J;sw=gw(:,:,o:o+nloc-1);sh=gh(:,o:o+nloc-1,:);sb=gb(:,o:o+nloc-1,:)
                if(p==0)then
                    w(:,:,slab_lo:slab_hi)=sw;h(:,slab_lo:slab_hi,:)=sh;b(:,slab_lo:slab_hi,:)=sb
                else
        communication_clock=MPI_Wtime()
                    call MPI_Send(sw,size(sw),MPI_DOUBLE_PRECISION,p,23000,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
                    call MPI_Send(sh,size(sh),MPI_DOUBLE_PRECISION,p,23001,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
                    call MPI_Send(sb,size(sb),MPI_DOUBLE_PRECISION,p,23002,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                endif
            enddo
        else
        communication_clock=MPI_Wtime()
            call MPI_Recv(sw,size(sw),MPI_DOUBLE_PRECISION,0,23000,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Recv(sh,size(sh),MPI_DOUBLE_PRECISION,0,23001,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Recv(sb,size(sb),MPI_DOUBLE_PRECISION,0,23002,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            w(:,:,slab_lo:slab_hi)=sw;h(:,slab_lo:slab_hi,:)=sh;b(:,slab_lo:slab_hi,:)=sb
        endif
        communication_clock=MPI_Wtime()
        call MPI_Bcast(t,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        call exchange_w(w);call exchange_micro(h);call exchange_micro(b)
    end subroutine
    subroutine contract(w,h,b,fw,fh,fb)
        real(KREAL),intent(inout)::w(:,:,:),h(:,:,:),b(:,:,:),fw(:,:,:),fh(:,:,:),fb(:,:,:)
        real(KREAL)::value,total(4),err
        integer::i,j,k,ipx,jg,fu,countlocal,countall,e
        logical::ok
        w=-999;h=-999;b=-999
        do j=slab_lo,slab_hi;do i=1,size(h,1)
            value=real(i+100*global_j(j),KREAL)
            do k=1,4;w(k,i,j)=value+10000*k;enddo
            do k=1,size(h,3);h(i,j,k)=value+10000*k;b(i,j,k)=-value-20000*k;enddo
        enddo;enddo
        call exchange_w(w);call exchange_micro(h);call exchange_micro(b);err=0
        ! Includes both x periodic corners and all 3 neighbour rows, not x-uniform input.
        do j=1,size(h,2);do i=-2,size(h,1)+3
            ipx=modulo(i-1,size(h,1))+1;jg=global_j(j);value=real(ipx+100*jg,KREAL)
            do k=1,4;err=max(err,abs(w(k,ipx,j)-value-10000*k));enddo
            do k=1,size(h,3)
                err=max(err,abs(h(ipx,j,k)-value-10000*k),abs(b(ipx,j,k)+value+20000*k))
            enddo
        enddo;enddo
        call require(err==0,'encoded cell halo/corner index contract')
        countlocal=merge(1,0,lower_wall())+merge(1,0,upper_wall())
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(countlocal,countall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        call require(countall==2.and.(lower_wall().eqv.(Block_XY(myID)%Origin_J==1)).and. &
            (upper_wall().eqv.(Block_XY(myID)%Origin_J+Ny-1==Ny_all)),'physical wall ownership')
        fw=w;fh=h;fb=b;call face_exchange(fw,fh,fb)
        countlocal=size(h,1)*(slab_hi-slab_lo+1)-merge(size(h,1),0,upper_wall())
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(countlocal,countall,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        call require(countall==Nx_all*(Ny_all-1),'interior face unique count')
        ok=collective(myID/=NP_XY-1);call require(.not.ok,'directed collective reject contract')
        if(myID==0)then
            open(newunit=fu,file='contract.txt',status='replace')
            write(fu,*)'PASS encoded W/h/b all 3 halo layers, x/y corners, two walls, faces=',countall
            write(fu,*)'PASS single-rank veto reaches all ranks; real-step rollback is separate'
            close(fu)
        endif
    end subroutine
    subroutine RunS2MPI1()
        real(KREAL),allocatable::w(:,:,:),h(:,:,:),b(:,:,:),u(:),v(:),wt(:),rn(:),rw(:)
        real(KREAL),allocatable::pw(:,:,:),ph(:,:,:),pb(:,:,:),sw(:,:,:),sh(:,:,:),sb(:,:,:)
        real(KREAL),allocatable::rsw(:,:,:),rsh(:,:,:),rsb(:,:,:)
        real(KREAL),allocatable::vxw(:,:,:),vyw(:,:,:),vxh(:,:,:),vxb(:,:,:),vyh(:,:,:),vyb(:,:,:),blw(:,:),blh(:,:),blb(:,:)
        real(KREAL),allocatable::mh(:,:,:),mb(:,:,:),xh(:,:,:),xb(:,:,:),yh(:,:,:),yb(:,:,:),eh(:,:,:),eb(:,:,:),kap(:,:,:,:)
        real(KREAL),allocatable::tmh(:,:,:),tmb(:,:,:),txh(:,:,:),txb(:,:,:),tyh(:,:,:),tyb(:,:,:),teh(:,:,:),teb(:,:,:),tkap(:,:,:,:)
        real(KREAL),allocatable::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:),lw(:,:),lh(:,:),lb(:,:)
        real(KREAL)::dx,dy,t,diag(12),faceerr,targeterr,fs,cs,total0(4),total(4),bottom(4),top(4),balance(4),nextbalance(4)
        real(KREAL)::drift,start,localcost,cost,qerr,gd(12),tau_stats(3),gtau(3)
        real(KREAL)::cumbottom(4),cumtop(4),initial_start,init_s,loop_start,loop_s,step_start,step_s
        real(KREAL)::warm_s,measured_s,comm0,io0,step_comm0,step_io0,warm_comm,measured_comm,warm_io,measured_io
        real(KREAL)::times(12),global_times(12),out_start,output_s,local_field_start,loop_comm,warm_end,loop_end
        character(len=64)::field_prefix
        integer(int64)::accepted,first_accepted,rss_peak,rss_max,rss_sum
        integer::segment_accepted,time_unit,pidunit
        logical::resumed,stopped
        integer::n,nv,i,j,jlocal,im,jm,step,lim,signed,fu,ledger,counts,e,wallcalls,totalwalls,failrank,ios
        integer(kind=8)::gcalls(3),tcalls(3)
        character(len=32)::flag,rankstr
        logical::ok,rejected
        initial_start=MPI_Wtime();call configure_io()
        call partion_mpi()
        write(rankstr,'(I2.2)')myID
        open(newunit=pidunit,file='pid_rank'//trim(rankstr)//'.txt',status='replace')
        write(pidunit,*)c_getpid();close(pidunit)
        call InitAllocation()
        call InitMesh()
        call require(Ny>=max(5,HO_HALO+2),'local slab too small for audited templates')
        slab_minus=ID_Y_Minus;slab_plus=ID_Y_Plus
        slab_lo=1+merge(HO_HALO,0,slab_minus/=MPI_PROC_NULL);slab_hi=slab_lo+Ny-1
        n=slab_hi+merge(HO_HALO,0,slab_plus/=MPI_PROC_NULL)
        slab_origin=Block_XY(myID)%Origin_J-slab_lo+1;slab_global_ny=Ny_all
        dx=ctr(1,1)%length(1);dy=ctr(1,1)%length(2);nv=VelocityCount()**2
        ! Uniform Cartesian geometry indexed by global origins; root verifies the actual mesh input.
        ok=.true.
        if(myID==0)then
            block
                integer::mu,ni,nj,ii,jj
                real(KREAL),allocatable::gx(:,:),gy(:,:)
                open(newunit=mu,file='Mesh2d.dat',form='unformatted',status='old',iostat=ios)
                ok=ios==0
                if(ok)then
                    read(mu,iostat=ios)ni,nj;ok=ios==0.and.ni==Nx_all+1.and.nj==Ny_all+1
                    if(ok)then
                        allocate(gx(ni,nj),gy(ni,nj));read(mu,iostat=ios)gx,gy;ok=ios==0
                        do jj=1,nj;do ii=1,ni
                            ok=ok.and.abs(gx(ii,jj)-(ii-1)*dx)<1e-14_KREAL.and.abs(gy(ii,jj)-(jj-1)*dy)<1e-14_KREAL
                        enddo;enddo
                    endif
                    close(mu)
                endif
            end block
        endif
        call require(ok,'uniform Cartesian mesh/global geometry check')
        allocate(w(4,Nx,n),h(Nx,n,nv),b(Nx,n,nv),u(nv),v(nv),wt(nv),rn(VelocityCount()),rw(VelocityCount()), &
            mh(Nx,n,nv),mb(Nx,n,nv),xh(Nx,n,nv),xb(Nx,n,nv),yh(Nx,n,nv),yb(Nx,n,nv),eh(Nx,n,nv),eb(Nx,n,nv), &
            kap(2,2,Nx,n),fxw(4,Nx,n),fyw(4,Nx,n),fxh(Nx,n,nv),fxb(Nx,n,nv),fyh(Nx,n,nv),fyb(Nx,n,nv), &
            lw(4,Nx),lh(Nx,nv),lb(Nx,nv))
        call init_velocity(u,v,wt,rn,rw)
        call require(abs(maxval(abs(u))-16.051374658065264_KREAL)<1e-14_KREAL.and. &
            abs(maxval(abs(v))-16.051374658065264_KREAL)<1e-14_KREAL,'MPI1 requires frozen half32 scale=2 table')
        call require(dt_global*(maxval(abs(u))/dx+maxval(abs(v))/dy)<=0.4_KREAL,'transport CFL')
        flag='';call get_environment_variable('S2_MPI1_CONTRACT',flag)
        if(trim(flag)=='1')then;call contract(w,h,b,fyw,fyh,fyb);return;endif
        start=initial_start
        call load_input(w,h,b,u,v,wt,t,accepted,total0,balance,cumbottom,cumtop,resumed)
        call initial_output(w,h,b,u,v,wt,t)
        call build_targets(w,u,v,wt,mh,mb,xh,xb,yh,yb,eh,eb,kap,targeterr,ok)
        call require(ok.and.targeterr<=gate,'initial targets')
        if(.not.resumed)then
            do i=1,4;total0(i)=sum(w(i,:,slab_lo:slab_hi));enddo;call sum4(total0)
        endif
        first_accepted=accepted;segment_accepted=0;wallcalls=0;rejected=.false.;stopped=.false.
        out_start=MPI_Wtime();call audit_cache(accepted,w,h,b,mh,mb,xh,xb,yh,yb,eh,eb,kap)
        io_seconds=io_seconds+MPI_Wtime()-out_start
        open(newunit=time_unit,file='step_time_rank'//trim(rankstr)//'.csv',status='replace')
        write(time_unit,'(A)')'segment_step,cumulative_step,warmup,whole_step_s,MPI_s,IO_s,peak_RSS_bytes'
        if(myID==0)then
            open(newunit=counts,file='initial_metadata.txt',status='replace')
            write(counts,*)'step',accepted;write(counts,*)'time',t
            write(counts,*)'total0',total0;write(counts,*)'balance',balance
            write(counts,*)'cumbottom',cumbottom;write(counts,*)'cumtop',cumtop;close(counts)
        endif
        failrank=-1;flag='';call get_environment_variable('S2_MPI1_FAIL_RANK',flag)
        if(len_trim(flag)>0)read(flag,*,iostat=ios)failrank
        if(failrank>=0)then
            write(rankstr,'(I2.2)')myID
            open(newunit=counts,file='accepted_targets_before_rank'//trim(rankstr)//'.bin', &
                status='replace',access='stream',form='unformatted')
            write(counts)mh,mb,xh,xb,yh,yb,eh,eb,kap
            close(counts)
        endif
        if(myID==0)then
            open(newunit=ledger,file='wall_balance.csv',status='replace')
            write(ledger,'(A)')'step,time,delta1,delta2,delta3,delta4,bottom1,bottom2,bottom3,bottom4,top1,top2,top3,top4,balance1,balance2,balance3,balance4,drift'
            open(newunit=fu,file='native_steps.csv',status='replace')
            write(fu,'(A)')'step,time,face_moment,target_moment,cn_residual,macro_micro,limited,signed,beta,cfl,hmin,hmax,bmin,bmax,global_drift'
        endif
        init_s=MPI_Wtime()-initial_start-io_seconds
        warm_s=0;measured_s=0;warm_comm=0;measured_comm=0;warm_io=0;measured_io=0
        comm0=communication_seconds;io0=io_seconds;loop_start=MPI_Wtime();warm_end=loop_start
        stopped=stop_requested(0)
        do step=1,MAX_ITER
            if(stopped)exit
            step_start=MPI_Wtime();step_comm0=communication_seconds;step_io0=io_seconds
            call exchange_w(w);call exchange_micro(h);call exchange_micro(b)
            call compute_flux(w,h,b,u,v,wt,dx,dy,dt_global,fxw,fyw,fxh,fxb,fyh,fyb,faceerr,lim,signed,ok,lw,lh,lb)
            call require(ok.and.faceerr<=gate,'face kernel')
            wallcalls=wallcalls+w2_wall_calls
            call face_exchange(fyw,fyh,fyb)
            ! Capture every y face once for comparison with serial, including both walls.
            if(audit_faces)then
            out_start=MPI_Wtime()
            block
                integer::yu
                character(len=80)::ypath
                write(ypath,'("yfaces_step",I10.10,"_rank",I2.2,".bin")')accepted+1,myID
                open(newunit=yu,file=trim(ypath),form='unformatted',access='stream',status='replace')
                write(yu)int(Nx,KINT),int(Ny,KINT),int(nv,KINT),int(Block_XY(myID)%Origin_J,KINT)
                write(yu)fyw(:,:,slab_lo:slab_hi),fyh(:,slab_lo:slab_hi,:),fyb(:,slab_lo:slab_hi,:),lw,lh,lb
                close(yu)
            end block
            io_seconds=io_seconds+MPI_Wtime()-out_start
            endif
            ! Preserve native face mapping: each face integral is multiplied by
            ! its metric and stored BEFORE residual summation (as original L_flux).
            vxw=fxw*0;vyw=fyw*0;vxh=fxh*0;vxb=fxb*0;vyh=fyh*0;vyb=fyb*0
            blw=lw;blh=lh;blb=lb
            do j=slab_lo,slab_hi;do i=1,Nx
                jlocal=j-slab_lo+1;im=modulo(i,Nx)+1
                vxw(:,i,j)=fxw(:,i,j)*vface(im,jlocal)%length
                vxh(i,j,:)=fxh(i,j,:)*vface(im,jlocal)%length
                vxb(i,j,:)=fxb(i,j,:)*vface(im,jlocal)%length
                vyw(:,i,j)=fyw(:,i,j)*hface(i,jlocal+1)%length
                vyh(i,j,:)=fyh(i,j,:)*hface(i,jlocal+1)%length
                vyb(i,j,:)=fyb(i,j,:)*hface(i,jlocal+1)%length
            enddo;enddo
            do i=1,Nx
                if(lower_wall())then
                    blw(:,i)=lw(:,i)*hface(i,1)%length
                    blh(i,:)=lh(i,:)*hface(i,1)%length;blb(i,:)=lb(i,:)*hface(i,1)%length
                else
                    vyw(:,i,slab_lo-1)=fyw(:,i,slab_lo-1)*hface(i,1)%length
                    vyh(i,slab_lo-1,:)=fyh(i,slab_lo-1,:)*hface(i,1)%length
                    vyb(i,slab_lo-1,:)=fyb(i,slab_lo-1,:)*hface(i,1)%length
                endif
            enddo
            rsw=w*0;rsh=h*0;rsb=b*0
            do j=slab_lo,slab_hi;do i=1,Nx
                im=modulo(i-2,Nx)+1;jm=j-1;jlocal=j-slab_lo+1
                if(j==1)then
                    rsw(:,i,j)=(vxw(:,im,j)-vxw(:,i,j)+blw(:,i)-vyw(:,i,j))/ctr(i,jlocal)%area
                    rsh(i,j,:)=(vxh(im,j,:)-vxh(i,j,:)+blh(i,:)-vyh(i,j,:))*(1.0_KREAL/ctr(i,jlocal)%area)
                    rsb(i,j,:)=(vxb(im,j,:)-vxb(i,j,:)+blb(i,:)-vyb(i,j,:))*(1.0_KREAL/ctr(i,jlocal)%area)
                else
                    rsw(:,i,j)=(vxw(:,im,j)-vxw(:,i,j)+vyw(:,i,jm)-vyw(:,i,j))/ctr(i,jlocal)%area
                    rsh(i,j,:)=(vxh(im,j,:)-vxh(i,j,:)+vyh(i,jm,:)-vyh(i,j,:))*(1.0_KREAL/ctr(i,jlocal)%area)
                    rsb(i,j,:)=(vxb(im,j,:)-vxb(i,j,:)+vyb(i,jm,:)-vyb(i,j,:))*(1.0_KREAL/ctr(i,jlocal)%area)
                endif
            enddo;enddo
            sw=w;sh=h;sb=b
            do j=slab_lo,slab_hi;do i=1,Nx
                sw(:,i,j)=w(:,i,j)+rsw(:,i,j)
                sh(i,j,:)=h(i,j,:)+rsh(i,j,:)
                sb(i,j,:)=b(i,j,:)+rsb(i,j,:)
            enddo;enddo
            pw=w;ph=h;pb=b;tmh=mh;tmb=mb;txh=xh;txb=xb;tyh=yh;tyb=yb;teh=eh;teb=eb;tkap=kap
            call advance_step(pw,ph,pb,u,v,wt,dx,dy,dt_global,tmh,tmb,txh,txb,tyh,tyb,teh,teb,tkap, &
                diag,fs,cs,ok,sw,sh,sb,faceerr,lim,signed)
            ! All ranks reach this vote even on a local collision-kernel failure.
            ok=collective(ok.and.all(diag(1:4)<=gate))
            do i=1,4;total(i)=sum(pw(i,:,slab_lo:slab_hi));enddo;call sum4(total)
            bottom=0;top=0
            do i=1,Nx
                if(lower_wall())bottom=bottom+lw(:,i)/dy
                if(upper_wall())top=top+fyw(:,i,slab_hi)/dy
            enddo
            call sum4(bottom);call sum4(top);nextbalance=balance+top-bottom
            drift=maxval(abs(total+nextbalance-total0)/max(1.0_KREAL,abs(total0)))
            ok=collective(ok.and.drift<=gate.and.(myID/=failrank.or.step/=1))
            if(.not.ok)then
                out_start=MPI_Wtime()
                call save_state('native_failed_proposal.bin',pw,ph,pb,u,v,wt,t+dt_global)
                ! w/h/b and accepted targets, time and ledger have never been mutated.
                call save_state('native_rejected_accepted.bin',w,h,b,u,v,wt,t)
                write(rankstr,'(I2.2)')myID
                open(newunit=counts,file='accepted_targets_after_rank'//trim(rankstr)//'.bin', &
                    status='replace',access='stream',form='unformatted')
                write(counts)mh,mb,xh,xb,yh,yb,eh,eb,kap
                close(counts)
                open(newunit=counts,file='rejection_rank' //trim(rankstr)//'.txt',status='replace')
                write(counts,*)'rejected, accepted=',accepted,' time=',t,' ledger=',balance
                close(counts);rejected=.true.;io_seconds=io_seconds+MPI_Wtime()-out_start;exit
            endif
            w=pw;h=ph;b=pb;mh=tmh;mb=tmb;xh=txh;xb=txb;yh=tyh;yb=tyb;eh=teh;eb=teb;kap=tkap
            t=t+dt_global;accepted=accepted+1;segment_accepted=segment_accepted+1;balance=nextbalance
            cumbottom=cumbottom+bottom;cumtop=cumtop+top
            call exchange_w(w);call exchange_micro(h);call exchange_micro(b)
        communication_clock=MPI_Wtime()
            call MPI_Allreduce(diag,gd,12,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Allreduce(diag(5:6),gd(5:6),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Allreduce(diag(9),gd(9),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
            call MPI_Allreduce(diag(11),gd(11),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            out_start=MPI_Wtime()
            if(myID==0.and.mod(segment_accepted,diagnostic_stride)==0)then
                write(fu,'(I0,14(",",ES24.16))')accepted,t,gd,drift
                write(ledger,'(I0,18(",",ES24.16))')accepted,t,total-total0,bottom,top,balance,drift
                flush(fu);flush(ledger)
            endif
            call audit_cache(accepted,w,h,b,mh,mb,xh,xb,yh,yb,eh,eb,kap)
            if(field_stride>0)then
                if(mod(accepted,int(field_stride,int64))==0)then
                    write(field_prefix,'("field_step",I10.10)')accepted
                    call parallel_field(trim(field_prefix),w,h,b,u,v,wt,t)
                endif
            endif
            io_seconds=io_seconds+MPI_Wtime()-out_start
            if(checkpoint_stride>0)then
                if(mod(accepted,int(checkpoint_stride,int64))==0) &
                    call publish_checkpoint(w,h,b,u,v,wt,t,accepted,total0,balance,cumbottom,cumtop,.false.)
            endif
            stopped=stop_requested(segment_accepted)
            step_s=MPI_Wtime()-step_start
            if(segment_accepted<=warmup_steps)then
                warm_s=warm_s+step_s;warm_comm=warm_comm+communication_seconds-step_comm0
                warm_io=warm_io+io_seconds-step_io0
            else
                measured_s=measured_s+step_s;measured_comm=measured_comm+communication_seconds-step_comm0
                measured_io=measured_io+io_seconds-step_io0
            endif
            write(time_unit,'(I0,",",I0,",",I0,3(",",ES24.16),",",I0)')segment_accepted,accepted, &
                merge(1,0,segment_accepted<=warmup_steps),step_s,communication_seconds-step_comm0, &
                io_seconds-step_io0,peak_rss()
            if(segment_accepted==warmup_steps)warm_end=MPI_Wtime()
        enddo
        loop_end=MPI_Wtime();loop_s=loop_end-loop_start;loop_comm=communication_seconds-comm0
        if(segment_accepted<=warmup_steps.and.warmup_steps>0)warm_end=loop_end
        warm_s=warm_end-loop_start;measured_s=loop_end-warm_end
        out_start=MPI_Wtime()
        call publish_checkpoint(w,h,b,u,v,wt,t,accepted,total0,balance,cumbottom,cumtop,.true.)
        local_field_start=MPI_Wtime()
        if(write_fields)call parallel_field('local_field',w,h,b,u,v,wt,t)
        io_seconds=io_seconds+MPI_Wtime()-local_field_start
        output_s=MPI_Wtime()-out_start
        close(time_unit)
        times=[init_s,warm_s,measured_s,loop_s,warm_comm,measured_comm,loop_comm, &
            warm_io,measured_io,output_s,io_seconds,MPI_Wtime()-initial_start]
        communication_clock=MPI_Wtime()
        call MPI_Reduce(times,global_times,12,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        rss_peak=peak_rss()
        communication_clock=MPI_Wtime()
        call MPI_Reduce(rss_peak,rss_max,1,MPI_INTEGER8,MPI_MAX,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Reduce(rss_peak,rss_sum,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(myID==0)then
            open(newunit=counts,file='timing.csv',status='replace')
            write(counts,'(A)')'ranks,segment_steps,warm_steps,measured_steps,init_read_s,warmup_loop_s,measured_loop_s,all_loop_s,warm_MPI_s,measured_MPI_s,all_MPI_s,warm_IO_s,measured_IO_s,final_output_s,total_IO_s,native_total_s,max_rank_peak_RSS_bytes,sum_rank_peak_RSS_upper_bytes,full_state_gathers'
            write(counts,'(4(I0,","),12(ES24.16,","),I0,",",I0,",",I0)')NP_XY,segment_accepted, &
                min(segment_accepted,warmup_steps),max(0,segment_accepted-warmup_steps),global_times,rss_max,rss_sum,checkpoint_gathers
            close(counts)
            open(newunit=counts,file='final_metadata.txt',status='replace')
            write(counts,*)'step',accepted;write(counts,*)'time',t
            write(counts,*)'total0',total0;write(counts,*)'balance',balance
            write(counts,*)'cumbottom',cumbottom;write(counts,*)'cumtop',cumtop
            write(counts,*)'controlled_stop',stopped;write(counts,*)'resumed',resumed;close(counts)
        endif
        localcost=MPI_Wtime()-start
        communication_clock=MPI_Wtime()
        call MPI_Reduce(localcost,cost,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        communication_clock=MPI_Wtime()
        call MPI_Reduce(wallcalls,totalwalls,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        write(rankstr,'(I2.2)')myID
        open(newunit=counts,file='work_rank'//trim(rankstr)//'.csv',status='replace')
        write(counts,'(A)')'rank,origin_j,local_ny,halo_ny,accepted,face_stages,cn_stages,wall_nodes'
        write(counts,'(8(I0,:,","))')myID,Block_XY(myID)%Origin_J,Ny,n-Ny,accepted, &
            step-merge(0,1,rejected),step-merge(0,1,rejected),wallcalls
        close(counts)
        gcalls=gas_tau_calls
        communication_clock=MPI_Wtime()
        call MPI_Reduce(gcalls,tcalls,3,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        tau_stats=gas_tau_min
        communication_clock=MPI_Wtime()
        call MPI_Reduce(tau_stats,gtau,3,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(myID==0)then
            close(fu);close(ledger)
            open(newunit=counts,file='native_call_counts.csv',status='replace')
            write(counts,'(A)')'ranks,accepted,face_stages,cn_stages,wall_node_calls,collective_rejected'
            write(counts,'(5(I0,","),L1)')NP_XY,accepted,step-merge(0,1,rejected),step-merge(0,1,rejected),totalwalls,rejected
            close(counts)
            open(newunit=counts,file='physical_gas_calls.txt',status='replace');write(counts,*)tcalls;write(counts,*)'tau_min',gtau;close(counts)
            open(newunit=counts,file='short_cost.txt',status='replace');write(counts,*)cost;close(counts)
            write(*,*)'MPI1 completed ranks=',NP_XY,' accepted=',accepted,' rejected=',rejected,' seconds=',cost
        endif
    end subroutine
end module

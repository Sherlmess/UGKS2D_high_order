! Local high-order workspace over the native Block_XY decomposition.
! ctr remains the accepted state. Three strips plus corners; no global state gather here.
module mod_S2CartesianMPI
    use iso_c_binding,only:c_int64_t
    use Flow_DIS,only:time_dis
    use flow_data
    use mpi
    use mod_S2MPI1Layout
    use mod_S2D1, only: compute_flux,build_targets,advance_step
    implicit none
    private
    public :: CartesianSetup,CartesianTargets,CartesianEvolution,CartesianAdvance,CollectiveRequire
    public :: GlobalSum4,DumpLocalAccepted,cart_native,CartesianHaloProbe,CartesianStats,mesh_dx,mesh_dy,CartesianLedger
    logical,save::cart_native=.false.
    real(KREAL),save::mesh_dx=0.0_KREAL,mesh_dy=0.0_KREAL
    real(KREAL),save::accepted_ledger(4,4)=0.0_KREAL
    real(KREAL),save::seam_maxabs=0,seam_scaled=0
    real(KREAL),save::step_diag(12)=0.0_KREAL
    interface
        integer(c_int64_t) function peak_rss() bind(C,name='mpi2_peak_rss')
            import c_int64_t
        end function
    end interface
contains
    subroutine CollectiveRequire(ok,reason,code)
        logical,intent(in)::ok
        character(len=*),intent(in)::reason
        integer,intent(in)::code
        integer::bad,first,e
        bad=huge(1);if(.not.ok)bad=myID
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(bad,first,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(first==huge(1))return
        if(allocated(ctr))call SaveCtrBeforeAbort()
        communication_clock=MPI_Wtime()
        call MPI_Barrier(MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        if(.not.ok)write(*,'(A,I0,2A)')'S2 MPI failure rank=',myID,' : ',trim(reason)
        call MPI_Abort(MPI_COMM_WORLD,code,e)
    end subroutine
    subroutine CartesianLedger(total0,balance,bottom,top)
        real(KREAL),intent(in)::total0(4),balance(4),bottom(4),top(4)
        accepted_ledger(:,1)=total0;accepted_ledger(:,2)=balance;accepted_ledger(:,3)=bottom;accepted_ledger(:,4)=top
    end subroutine
    subroutine SaveCtrBeforeAbort()
        real(KREAL),allocatable::w(:,:,:),h(:,:,:),b(:,:,:)
        integer::i,j,fu
        allocate(w(4,nx,ny),h(nx,ny,Nvel),b(nx,ny,Nvel))
        do j=1,ny;do i=1,nx
            w(:,i,j)=ctr(i,j)%conVars;h(i,j,:)=ctr(i,j)%h_vec;b(i,j,:)=ctr(i,j)%b_vec
        enddo;enddo
        call DumpLocalAccepted(w,h,b,time_dis,'last_accepted_emergency')
        if(myID==0)then
            open(newunit=fu,file='last_accepted_emergency_ledger.txt',status='replace')
            write(fu,*)iter,time_dis,accepted_ledger
            close(fu)
        endif
    end subroutine
    subroutine CartesianStats(loop_s,output_s,comm_s)
        real(KREAL),intent(in)::loop_s,output_s,comm_s
        integer::fu
        character(len=80)::path
        if(.not.cart_native)return
        write(path,'("runtime_rank",I4.4,".csv")')myID
        open(newunit=fu,file=trim(path),status='replace')
        write(fu,'(A)')'rank,loop_s,communication_s,output_s,peak_rss_bytes,seam_maxabs,seam_scaled'
        write(fu,'(I0,3(",",ES24.16),",",I0,2(",",ES24.16))')myID,loop_s,comm_s,output_s,peak_rss(),seam_maxabs,seam_scaled
        close(fu)
        write(path,'("owner_cost_rank",I4.4,".csv")')myID
        open(newunit=fu,file=trim(path),status='replace')
        write(fu,'(A)')'rank,reconstruction_s,pack_s,mpi_wait_transfer_s,unpack_s,sync_s,exchange_total_s,sent_payload_bytes,reconstructed_cells,exchange_calls,buffer_setup_s,buffer_allocations,buffer_capacity_bytes'
        write(fu,'(I0,6(",",ES24.16),3(",",I0),",",ES24.16,2(",",I0))')myID,owner_reconstruction_s,owner_pack_s,owner_mpi_s,owner_unpack_s,owner_sync_s,owner_exchange_s,owner_bytes,owner_cells,owner_calls,owner_buffer_setup_s,owner_buffer_allocations,owner_buffer_capacity_bytes
        close(fu)
        write(path,'("gates_rank",I4.4,".txt")')myID
        open(newunit=fu,file=trim(path),status='replace');write(fu,'(12ES24.16)')step_diag;close(fu)
    end subroutine
    subroutine CartesianSetup()
        integer::r,e
        real(KREAL)::spacing(2)
        character(len=160)::msg
        cart_native=(trim(TimeIntegration_method)=='S2CN2'.or.IsHighOrderImplicit()).and.S2W2R_shared_predictor==0.and.S2W1_time_centered==0.and.S2W1D_face_time_consistent==0
        if(.not.cart_native.and.NP_XY/=1)call CollectiveRequire(.false.,'research midpoint/shared predictors remain single rank',475)
        if(.not.cart_native)return
        ! No rank-count whitelist. Every rank must own at least one real cell.
        call CollectiveRequire(Nx_all>=5.and.Ny_all>=5,'global reconstruction grid requires at least five cells per axis',470)
        cart_halo=3;cart_forward_halo=.false.
        do r=0,NP_XY-1
            call CollectiveRequire(Block_XY(r)%NI>=1.and.Block_XY(r)%NJ>=1,'empty Cartesian block',470)
            if(Block_XY(r)%NI<3.or.Block_XY(r)%NJ<3)cart_forward_halo=.true.
            ! A wall cell may need the fifth physical row across several ranks.
            if(S2W1_wall_mode==1.and.Block_XY(r)%NJ==1)cart_halo=4
        enddo
        if(cart_halo==4)cart_forward_halo=.true.
        if(myID==0)spacing=ctr(1,1)%length
        communication_clock=MPI_Wtime()
        call MPI_Bcast(spacing,2,MPI_KREAL,0,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        mesh_dx=spacing(1);mesh_dy=spacing(2)
        cart_xlo=cart_halo+1;cart_xhi=nx+cart_halo;cart_ylo=cart_halo+1;cart_yhi=ny+cart_halo
        cart_oi=Block_XY(myID)%Origin_I;cart_oj=Block_XY(myID)%Origin_J
        cart_nx=Nx_all;cart_ny=Ny_all
        cart_xminus=ID_X_minus;cart_xplus=ID_X_plus;cart_yminus=ID_Y_minus;cart_yplus=ID_Y_plus
        write(*,'(A,I0,A,2I6,A,2I6)')'S2 native block rank=',myID,' origin=',cart_oi,cart_oj,' cells=',nx,ny
    end subroutine
    function GlobalSum4(local) result(global)
        real(KREAL),intent(in)::local(4)
        real(KREAL)::global(4)
        integer::e
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(local,global,4,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
    end function
    subroutine PadState(w,h,b,pw,ph,pb)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:)
        real(KREAL),allocatable,intent(out)::pw(:,:,:),ph(:,:,:),pb(:,:,:)
        real(KREAL),allocatable::tmp(:,:,:)
        integer::k
        allocate(pw(4,nx+2*cart_halo,ny+2*cart_halo),ph(nx+2*cart_halo,ny+2*cart_halo,Nvel),pb(nx+2*cart_halo,ny+2*cart_halo,Nvel),tmp(Nvel,nx+2*cart_halo,ny+2*cart_halo))
        pw=0;pw(:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)=w
        cart_context=.true.;call exchange_cart(pw)
        tmp=0
        do k=1,Nvel;tmp(k,cart_xlo:cart_xhi,cart_ylo:cart_yhi)=h(:,:,k);enddo
        call exchange_cart(tmp)
        do k=1,Nvel;ph(:,:,k)=tmp(k,:,:);enddo
        do k=1,Nvel;tmp(k,cart_xlo:cart_xhi,cart_ylo:cart_yhi)=b(:,:,k);enddo
        call exchange_cart(tmp)
        do k=1,Nvel;pb(:,:,k)=tmp(k,:,:);enddo
    end subroutine
    subroutine CartesianTargets(w,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,error,ok)
        real(KREAL),intent(in)::w(:,:,:)
        real(KREAL),intent(out)::mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:),eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
        real(KREAL),intent(out)::error
        logical,intent(out)::ok
        real(KREAL),allocatable::pw(:,:,:),z(:,:,:,:),kap(:,:,:,:)
        integer::px,py
        px=nx+2*cart_halo;py=ny+2*cart_halo
        allocate(pw(4,px,py),z(px,py,Nvel,8),kap(2,2,px,py));pw=0
        pw(:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)=w
        cart_context=.true.;call exchange_cart(pw)
        call build_targets(pw,uVel,vVel,wVel,z(:,:,:,1),z(:,:,:,2),z(:,:,:,3),z(:,:,:,4),z(:,:,:,5),z(:,:,:,6),z(:,:,:,7),z(:,:,:,8),kap,error,ok)
        cart_context=.false.
        call CollectiveRequire(ok,'volume target rejected',471)
        mh=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,1);mb=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,2)
        kxh=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,3);kxb=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,4)
        kyh=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,5);kyb=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,6)
        eah=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,7);eab=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,8);kappa=kap(:,:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)
    end subroutine
    subroutine SharedFace(axis,macro,microh,microb)
        integer,intent(in)::axis
        real(KREAL),intent(inout)::macro(:,:,:),microh(:,:,:),microb(:,:,:)
        real(KREAL),allocatable::send(:,:),recv(:,:),left(:,:)
        real(KREAL)::err,scaled,global(2),local(2)
        integer::n,m,i,j,lo,hi,minus,plus,e,st(MPI_STATUS_SIZE)
        n=merge(ny,nx,axis==1);lo=cart_halo;hi=merge(cart_xhi,cart_yhi,axis==1)
        minus=merge(ID_X_minus,ID_Y_minus,axis==1);plus=merge(ID_X_plus,ID_Y_plus,axis==1)
        allocate(send(4+2*Nvel,n),recv(4+2*Nvel,n),left(4+2*Nvel,n))
        do m=1,n
            i=hi;j=m+cart_halo;if(axis==2)then;i=m+cart_halo;j=hi;endif
            send(:,m)=[macro(:,i,j),microh(i,j,:),microb(i,j,:)]
            i=lo;j=m+cart_halo;if(axis==2)then;i=m+cart_halo;j=lo;endif
            left(:,m)=[macro(:,i,j),microh(i,j,:),microb(i,j,:)]
        enddo
        communication_clock=MPI_Wtime()
        call MPI_Sendrecv(send,size(send),MPI_KREAL,plus,22100+axis,recv,size(recv),MPI_KREAL,minus,22100+axis,MPI_COMM_WORLD,st,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        err=0;scaled=0
        if(minus/=MPI_PROC_NULL)then
            err=maxval(abs(left-recv))
            scaled=maxval(abs(left-recv))/max(1.0_KREAL,maxval(abs(left)),maxval(abs(recv)))
        endif
        local=[err,scaled];call MPI_Allreduce(local,global,2,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,e)
        seam_maxabs=max(seam_maxabs,global(1));seam_scaled=max(seam_scaled,global(2))
        call CollectiveRequire(scaled<=2.0e-12_KREAL,'complete shared macro/h/b face discrepancy',472)
        if(minus/=MPI_PROC_NULL)then
            do m=1,n
                i=lo;j=m+cart_halo;if(axis==2)then;i=m+cart_halo;j=lo;endif
                macro(:,i,j)=recv(1:4,m);microh(i,j,:)=recv(5:4+Nvel,m);microb(i,j,:)=recv(5+Nvel:,m)
            enddo
        endif
    end subroutine
    subroutine CartesianEvolution(w,h,b,fxw,fyw,fxh,fxb,fyh,fyb,lower_w,lower_h,lower_b,error,limited,signed,ok)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:)
        real(KREAL),intent(out)::fxw(:,:,:),fyw(:,:,:),fxh(:,:,:),fxb(:,:,:),fyh(:,:,:),fyb(:,:,:)
        real(KREAL),allocatable,intent(inout)::lower_w(:,:),lower_h(:,:),lower_b(:,:)
        real(KREAL),intent(out)::error
        integer,intent(out)::limited,signed
        logical,intent(out)::ok
        real(KREAL),allocatable::pw(:,:,:),ph(:,:,:),pb(:,:,:),xw(:,:,:),yw(:,:,:),xh(:,:,:),xb(:,:,:),yh(:,:,:),yb(:,:,:),lw(:,:),lh(:,:),lb(:,:)
        integer::i,j,px,py,k,fu
        character(len=80)::path
        px=nx+2*cart_halo;py=ny+2*cart_halo
        call PadState(w,h,b,pw,ph,pb)
        allocate(xw(4,px,py),yw(4,px,py),xh(px,py,Nvel),xb(px,py,Nvel),yh(px,py,Nvel),yb(px,py,Nvel),lw(4,px),lh(px,Nvel),lb(px,Nvel))
        call compute_flux(pw,ph,pb,uVel,vVel,wVel,mesh_dx,mesh_dy,dt_global,xw,yw,xh,xb,yh,yb,error,limited,signed,ok,lw,lh,lb)
        cart_context=.false.
        call CollectiveRequire(ok,'face rejected before publication',473)
        call SharedFace(1,xw,xh,xb);call SharedFace(2,yw,yh,yb)
        do j=1,ny;do i=1,nx+1
            vface(i,j)%L_flux=xw(:,i+cart_halo-1,j+cart_halo)*vface(i,j)%length
            vface(i,j)%L_flux_h_vec=xh(i+cart_halo-1,j+cart_halo,:)*vface(i,j)%length
            vface(i,j)%L_flux_b_vec=xb(i+cart_halo-1,j+cart_halo,:)*vface(i,j)%length
            vface(i,j)%R_flux=vface(i,j)%L_flux
            vface(i,j)%R_flux_h_vec=vface(i,j)%L_flux_h_vec;vface(i,j)%R_flux_b_vec=vface(i,j)%L_flux_b_vec
        enddo;enddo
        do j=1,ny+1;do i=1,nx
            hface(i,j)%L_flux=yw(:,i+cart_halo,j+cart_halo-1)*hface(i,j)%length
            hface(i,j)%L_flux_h_vec=yh(i+cart_halo,j+cart_halo-1,:)*hface(i,j)%length
            hface(i,j)%L_flux_b_vec=yb(i+cart_halo,j+cart_halo-1,:)*hface(i,j)%length
            if(j==1.and.ID_Y_minus==MPI_PROC_NULL)then
                hface(i,j)%L_flux=lw(:,i+cart_halo)*hface(i,j)%length
                hface(i,j)%L_flux_h_vec=lh(i+cart_halo,:)*hface(i,j)%length
                hface(i,j)%L_flux_b_vec=lb(i+cart_halo,:)*hface(i,j)%length
            endif
            hface(i,j)%R_flux=hface(i,j)%L_flux
            hface(i,j)%R_flux_h_vec=hface(i,j)%L_flux_h_vec;hface(i,j)%R_flux_b_vec=hface(i,j)%L_flux_b_vec
        enddo;enddo
        fxw=xw(:,cart_xlo:cart_xhi,cart_ylo:cart_yhi);fyw=yw(:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)
        fxh=xh(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:);fxb=xb(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:)
        fyh=yh(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:);fyb=yb(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:)
        if(allocated(lower_w))then;lower_w=lw(:,cart_xlo:cart_xhi);lower_h=lh(cart_xlo:cart_xhi,:);lower_b=lb(cart_xlo:cart_xhi,:);endif
        if(S2Diagnostics)then
            write(path,'("faces_rank",I4.4,"_step",I4.4,".bin")')myID,iter
            open(newunit=fu,file=trim(path),access='stream',form='unformatted',status='replace')
            write(fu)nx,ny,Nvel,cart_oi,cart_oj
            do j=1,ny;do i=1,nx+1
                write(fu)vface(i,j)%L_flux,vface(i,j)%L_flux_h_vec,vface(i,j)%L_flux_b_vec
            enddo;enddo
            do j=1,ny+1;do i=1,nx
                write(fu)hface(i,j)%L_flux,hface(i,j)%L_flux_h_vec,hface(i,j)%L_flux_b_vec
            enddo;enddo
            close(fu)
        endif
    end subroutine
    subroutine CartesianAdvance(w,h,b,mh,mb,kxh,kxb,kyh,kyb,eah,eab,kappa,diag,flux_s,collision_s,ok,wn,hs,bs,error,limited,signed)
        real(KREAL),intent(inout)::w(:,:,:),h(:,:,:),b(:,:,:),mh(:,:,:),mb(:,:,:),kxh(:,:,:),kxb(:,:,:),kyh(:,:,:),kyb(:,:,:),eah(:,:,:),eab(:,:,:),kappa(:,:,:,:)
        real(KREAL),intent(in)::wn(:,:,:),hs(:,:,:),bs(:,:,:),error
        integer,intent(in)::limited,signed
        real(KREAL),intent(out)::diag(12),flux_s,collision_s
        logical,intent(out)::ok
        real(KREAL),allocatable::pw(:,:,:),ph(:,:,:),pb(:,:,:),z(:,:,:,:),kap(:,:,:,:),sw(:,:,:),sh(:,:,:),sb(:,:,:)
        integer::px,py,e,fu,st,inject,inject_step,bad,first
        character(len=32)::flag
        px=nx+2*cart_halo;py=ny+2*cart_halo
        call PadState(w,h,b,pw,ph,pb)
        allocate(z(px,py,Nvel,8),kap(2,2,px,py),sw(4,px,py),sh(px,py,Nvel),sb(px,py,Nvel));z=0;kap=0
        z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,1)=mh;z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,2)=mb
        z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,3)=kxh;z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,4)=kxb
        z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,5)=kyh;z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,6)=kyb
        z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,7)=eah;z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,8)=eab;kap(:,:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)=kappa
        sw=pw;sh=ph;sb=pb;sw(:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)=wn;sh(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:)=hs;sb(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:)=bs
        call advance_step(pw,ph,pb,uVel,vVel,wVel,mesh_dx,mesh_dy,dt_global,z(:,:,:,1),z(:,:,:,2),z(:,:,:,3),z(:,:,:,4),z(:,:,:,5),z(:,:,:,6),z(:,:,:,7),z(:,:,:,8),kap,diag,flux_s,collision_s,ok,sw,sh,sb,error,limited,signed)
        cart_context=.false.
        flag='';call get_environment_variable('S2E3_TEST_REJECT_RANK',flag,status=st)
        if(st==0)then
            read(flag,*,iostat=st)inject
            if(st==0)then
                inject_step=1
                flag='';call get_environment_variable('S2E3_TEST_REJECT_STEP',flag,status=st)
                if(st==0)read(flag,*,iostat=st)inject_step
                if(inject==myID.and.iter+1==inject_step)ok=.false.
            endif
        endif
        bad=huge(1);if(.not.ok)bad=myID
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(bad,first,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        ok=first==huge(1)
        if(.not.ok)then
            if(myID==first)write(*,'(A,I0,A,2I8)')'S2 collision/target rejected rank=',first,' global cell=',cart_oi+failure_i-cart_xlo,cart_oj+failure_j-cart_ylo
            return
        endif
        step_diag=diag
        w=pw(:,cart_xlo:cart_xhi,cart_ylo:cart_yhi);h=ph(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:);b=pb(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:)
        mh=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,1);mb=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,2);kxh=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,3);kxb=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,4)
        kyh=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,5);kyb=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,6);eah=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,7);eab=z(cart_xlo:cart_xhi,cart_ylo:cart_yhi,:,8);kappa=kap(:,:,cart_xlo:cart_xhi,cart_ylo:cart_yhi)
    end subroutine
    subroutine DumpLocalAccepted(w,h,b,time,stem)
        real(KREAL),intent(in)::w(:,:,:),h(:,:,:),b(:,:,:),time
        character(len=*),intent(in)::stem
        character(len=160)::path
        integer::fu,state_step
        state_step=iter
        if(trim(stem)=='accepted')state_step=iter+1
        write(path,'(A,"_rank",I4.4,".bin")')trim(stem),myID
        open(newunit=fu,file=trim(path),status='replace',access='stream',form='unformatted')
        write(fu)nx,ny,Nvel,cart_oi,cart_oj,state_step,time,w,h,b,uVel,vVel,wVel
        close(fu)
    end subroutine
    subroutine CartesianHaloProbe()
        real(KREAL),allocatable::f(:,:,:)
        real(KREAL)::expected,error,globalerror
        integer::i,j,k,gi,gj,e,fu
        character(len=80)::path
        allocate(f(4+2*Nvel,nx+2*cart_halo,ny+2*cart_halo));f=-999
        do j=cart_ylo,cart_yhi;do i=cart_xlo,cart_xhi;do k=1,4+2*Nvel
            f(k,i,j)=real(k+10000*(cart_oi+i-cart_xlo)+1000000*(cart_oj+j-cart_ylo),KREAL)
        enddo;enddo;enddo
        cart_context=.true.;call exchange_cart(f);cart_context=.false.;error=0
        do j=1,ny+2*cart_halo;do i=1,nx+2*cart_halo
            gi=modulo(cart_oi+i-cart_xlo-1,Nx_all)+1;gj=modulo(cart_oj+j-cart_ylo-1,Ny_all)+1
            if(S2W1_wall_mode==1.and.(cart_oj+j-cart_ylo<1.or.cart_oj+j-cart_ylo>Ny_all))cycle
            do k=1,4+2*Nvel
                expected=real(k+10000*gi+1000000*gj,KREAL);error=max(error,abs(f(k,i,j)-expected))
            enddo
        enddo;enddo
        communication_clock=MPI_Wtime()
        call MPI_Allreduce(error,globalerror,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,e)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        write(path,'("halo_rank",I4.4,".txt")')myID
        open(newunit=fu,file=trim(path),status='replace');write(fu,*)myID,nx,ny,cart_oi,cart_oj,error,globalerror;close(fu)
        call CollectiveRequire(error==0.0_KREAL,'complete W/h/b halo and corner encoding',474)
    end subroutine
end module

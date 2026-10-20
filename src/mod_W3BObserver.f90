module mod_W3BObserver
    use flow_para, only: ck
    use control_para, only: MAX_ITER
    use flow_data, only: KREAL,KINT,nx,ny,Nvel,ctr,uVel,vVel,wVel,dt_global
    use mpi, only: MPI_Wtime, MPI_Finalize
    use,intrinsic::ieee_arithmetic
    implicit none
    private
    public::W3B_Init,W3B_Observe
    logical,save::enabled=.false.
    integer,save::ru,su,checkpoint_stride=0,latest_stride=1
    real(KREAL),save::origin=0.1_KREAL,obs_seconds=0.0_KREAL,output_seconds=0.0_KREAL
    real(KREAL),save::scales(14)=1.0_KREAL
    real(KREAL),allocatable,save::w(:,:,:),h(:,:,:),b(:,:,:),oh(:,:,:),ob(:,:,:),q(:,:,:),oq(:,:,:)
    character(len=18),parameter::names(14)=[character(len=18)::'rho','Ux','Uy','theta','theta_perturbation', &
        'Pxy','qx','qy','rhoU','rhoV','rhoE','W_density','h','b']
contains
    subroutine capture()
        integer::i,j,k
        real(KREAL)::rho,ux,uy,theta,cx,cy,e
        do j=1,ny;do i=1,nx
            w(:,i,j)=ctr(i,j)%conVars;h(i,j,:)=ctr(i,j)%h_vec;b(i,j,:)=ctr(i,j)%b_vec
            rho=w(1,i,j);ux=w(2,i,j)/rho;uy=w(3,i,j)/rho
            theta=(2.0_KREAL/real(ck+2,KREAL))*(w(4,i,j)-0.5_KREAL*rho*(ux*ux+uy*uy))/rho
            q(:,i,j)=[rho,ux,uy,theta,theta-1.0_KREAL,0.0_KREAL,0.0_KREAL,0.0_KREAL,w(2:4,i,j),rho]
            do k=1,Nvel
                cx=uVel(k)-ux;cy=vVel(k)-uy;e=(cx*cx+cy*cy)*h(i,j,k)+b(i,j,k)
                q(6,i,j)=q(6,i,j)+wVel(k)*cx*cy*h(i,j,k)
                q(7,i,j)=q(7,i,j)+0.5_KREAL*wVel(k)*cx*e
                q(8,i,j)=q(8,i,j)+0.5_KREAL*wVel(k)*cy*e
            enddo
        enddo;enddo
    end subroutine
    subroutine snapshot(path,step)
        character(len=*),intent(in)::path
        integer,intent(in)::step
        integer::fu,err
        real(KREAL)::t0
        t0=MPI_Wtime()
        open(newunit=fu,file=path//'.tmp',status='replace',access='stream',form='unformatted')
        write(fu)int(nx,KINT),int(ny,KINT),int(Nvel,KINT),origin+real(step,KREAL)*dt_global,w,h,b,uVel,vVel,wVel
        close(fu)
        call rename(path//'.tmp',path,err)
        if(err/=0)error stop 'W3B output rename failed'
        output_seconds=output_seconds+MPI_Wtime()-t0
    end subroutine
    subroutine W3B_Init()
        character(len=128)::env
        integer::status
        real(KREAL)::t0
        call get_environment_variable('W3B_OBSERVE',env,status=status)
        enabled=status==0.and.trim(env)=='1'
        if(.not.enabled)return
        t0=MPI_Wtime()
        call get_environment_variable('W3_ORIGIN',env,status=status)
        if(status==0)read(env,*)origin
        call get_environment_variable('W3_CHECKPOINT_STRIDE',env,status=status)
        if(status==0)read(env,*)checkpoint_stride
        call get_environment_variable('W3_LATEST_STRIDE',env,status=status)
        if(status==0)read(env,*)latest_stride
        if(latest_stride<1)error stop 'W3_LATEST_STRIDE must be positive'
        call get_environment_variable('W3B_H_RMS',env);read(env,*)scales(13)
        call get_environment_variable('W3B_B_RMS',env);read(env,*)scales(14)
        scales(2:3)=.05_KREAL;scales(9:10)=.05_KREAL
        allocate(w(4,nx,ny),h(nx,ny,Nvel),b(nx,ny,Nvel),oh(nx,ny,Nvel),ob(nx,ny,Nvel),q(12,nx,ny),oq(12,nx,ny))
        call capture();oh=h;ob=b;oq=q
        call snapshot('checkpoint_0000.bin',0);call snapshot('latest_accepted.bin',0)
        open(newunit=ru,file='update_rates.csv',status='replace')
        write(ru,'(A)')'step,time,field,absolute_rms_rate,absolute_max_rate,fixed_scale,normalized_rms_rate,normalized_max_rate,state_rms,state_maxabs,x_nonuniform'
        open(newunit=su,file='observer_steps.csv',status='replace')
        write(su,'(A)')'step,time,h_min,b_min,finite,native_loop_s,observer_s,checkpoint_output_s'
        obs_seconds=MPI_Wtime()-t0
    end subroutine
    subroutine scalar_row(step,m,a,old)
        integer,intent(in)::step,m
        real(KREAL),intent(in)::a(:,:),old(:,:)
        real(KREAL)::rms,mx,nu,avg(size(a,2))
        integer::i
        rms=sqrt(sum(((a-old)/dt_global)**2)/real(size(a),KREAL));mx=maxval(abs((a-old)/dt_global))
        avg=sum(a,dim=1)/real(size(a,1),KREAL);nu=0.0_KREAL
        do i=1,size(a,1);nu=max(nu,maxval(abs(a(i,:)-avg)));enddo
        write(ru,'(I0,",",ES24.16,",",A,8(",",ES24.16))')step,origin+step*dt_global,trim(names(m)), &
            rms,mx,scales(m),rms/scales(m),mx/scales(m),sqrt(sum(a*a)/real(size(a),KREAL)),maxval(abs(a)),nu
    end subroutine
    subroutine micro_row(step,m,a,old)
        integer,intent(in)::step,m
        real(KREAL),intent(in)::a(:,:,:),old(:,:,:)
        real(KREAL)::rms,mx,nu
        integer::j,k
        rms=sqrt(sum(((a-old)/dt_global)**2)/real(size(a),KREAL));mx=maxval(abs((a-old)/dt_global));nu=0
        do k=1,Nvel;do j=1,ny
            nu=max(nu,maxval(abs(a(:,j,k)-sum(a(:,j,k))/real(nx,KREAL))))
        enddo;enddo
        write(ru,'(I0,",",ES24.16,",",A,8(",",ES24.16))')step,origin+step*dt_global,trim(names(m)), &
            rms,mx,scales(m),rms/scales(m),mx/scales(m),sqrt(sum(a*a)/real(size(a),KREAL)),maxval(abs(a)),nu
    end subroutine
    subroutine W3B_Observe(step,loop_seconds)
        integer,intent(in)::step
        real(KREAL),intent(in)::loop_seconds
        integer::m,stop_unit,stop_error
        real(KREAL)::t0
        character(len=64)::path
        logical::finite,save_now,stop_requested
        if(.not.enabled)return
        t0=MPI_Wtime();call capture()
        do m=1,12;call scalar_row(step,m,q(m,:,:),oq(m,:,:));enddo
        call micro_row(step,13,h,oh);call micro_row(step,14,b,ob)
        finite=all(ieee_is_finite(w)).and.all(ieee_is_finite(h)).and.all(ieee_is_finite(b))
        if(mod(step,latest_stride)==0.or.step==MAX_ITER)call snapshot('latest_accepted.bin',step)
        save_now=any(step==[80,240,400,560,720])
        if(checkpoint_stride>0)save_now=mod(step,checkpoint_stride)==0
        if(save_now)then
            write(path,'("checkpoint_",I8.8,".bin")')step
            call snapshot(trim(path),step)
        endif
        oh=h;ob=b;oq=q
        obs_seconds=obs_seconds+MPI_Wtime()-t0
        write(su,'(I0,3(",",ES24.16),",",L1,3(",",ES24.16))') &
            step,origin+step*dt_global,minval(h),minval(b),finite,loop_seconds,obs_seconds,output_seconds
        flush(ru);flush(su)
        if(.not.finite)error stop 'W3B nonfinite accepted state'
        inquire(file='RESOURCE_STOP_REQUEST',exist=stop_requested)
        if(stop_requested)then
            call snapshot('latest_accepted.bin',step)
            open(newunit=stop_unit,file='RESOURCE_STOP_ACCEPTED.txt',status='replace')
            write(stop_unit,*)step,origin+step*dt_global
            close(stop_unit)
            call MPI_Finalize(stop_error)
            stop 0
        endif
    end subroutine
end module

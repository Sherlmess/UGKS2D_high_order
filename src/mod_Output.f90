
module Output
    use mod_S2MPI1Layout,only:communication_seconds,communication_clock
    use flow_data
    use Flow_DIS
    use mod_ClosedCouetteMass,only:WriteMassRestart,restart_step
    use mod_ImplicitDiagnostics,only:GetImplicitReferences
    use Tools
    use mpi
    use mod_S2D1, only: IsS2CN2,IsS2Steady
    use mod_S2NativeAdapter, only: S2NativeWriteOutput,S2NativeOutputHeat,S2NativeCommitRestart,S2NativeRunSummary,S2NativePublishWallDiagnostics
    implicit none
    
contains
    !--------------------------------------------------
    !>Write result
    !--------------------------------------------------
    subroutine Output_result(write_restart_now)
        logical,optional,intent(in)::write_restart_now
        logical::save_restart
        save_restart=.true.
        if(present(write_restart_now))save_restart=write_restart_now
        if(IsS2Steady().or.trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau') then
            call S2NativeWriteOutput()
            return
        endif
        if(IsS2CN2().and.UnifiedPhysicalOperators.and.S2W1_wall_mode==1)call S2NativePublishWallDiagnostics()
        if(Iflag_WriteSurf==1)      call write_wallboundary
        call write_flow
        if(IsS2CN2().and..not.IsHighOrderImplicit())then
            if(Iflag_WriteSurf==1.and.S2W1_wall_mode==1)call write_wallboundary(2)
            if(S2Diagnostics)call S2NativeWriteOutput()
            call S2NativeRunSummary()
        endif
        if(save_restart)call Output_restart_only()
    end subroutine Output_result

    subroutine Output_restart_only()
        integer::unit
        real(KREAL)::ref_h,ref_b
        if(Iflag_WriteRestart/=1)return
        if(IsS2CN2().and..not.IsHighOrderImplicit())then
            call write_restart('flow2d.dat.next','vdf2d.dat.next','grid2d.grd.next')
            call S2NativeCommitRestart()
        else
            call write_restart()
        endif
        if(ClosedCouetteMassRescale)then
            call WriteMassRestart()
            call GetImplicitReferences(ref_h,ref_b)
            if(myID==0)then
                ! Preserve the ORIGINAL reference epoch across all complete restarts.
                call WriteRestartReference(ref_h,ref_b)
            endif
        endif
    end subroutine Output_restart_only

    subroutine WriteRestartReference(ref_h,ref_b)
        real(KREAL),intent(in)::ref_h,ref_b
        integer::unit,ios,epoch
        real(KREAL)::old_h,old_b,rtol,mtol
        logical::exists
        epoch=0
        inquire(file='residual_reference.in',exist=exists)
        if(exists)then
            open(newunit=unit,file='residual_reference.in',status='old',iostat=ios)
            if(ios/=0)error stop 'restart reference read failed'
            read(unit,*,iostat=ios)epoch,old_h,old_b,rtol,mtol
            close(unit)
            if(ios/=0)error stop 'restart reference parse failed'
        endif
        open(newunit=unit,file='residual_reference.in',status='replace',iostat=ios)
        if(ios/=0)error stop 'restart reference write failed'
        write(unit,'(i0,1x,4(es25.16e3,1x))')epoch,ref_h,ref_b,Implicit_ResidualTol,Implicit_MomentTol
        close(unit)
    end subroutine WriteRestartReference

    
    
    !---------------------------------------
    !>Output flow field: conVars based
    !---------------------------------------
    subroutine write_flow
        real(KREAL)                                     :: prim(4),prim0(4),solu0(10)
        integer(KINT)                                   :: i,j,k,pi,pj,np_recv,i1,j1,Ndata,tag
        character(len=20)                               :: str
        real(KREAL), allocatable, dimension(:,:,:)      :: solution,solution_all
        real(KREAL), allocatable, dimension(:,:,:)      :: buff1,buff2
        real(KREAL), allocatable, dimension(:)        :: surfx,surfy
        real(KREAL), allocatable, dimension(:,:)        :: FlowVars_all


        !Flow variables to be output:
        !1-10: density, u, v, T, p, txx, txy, tyy, qx, qy
        !11-14: per-cell residuals of conservative variables
        Ndata = 14
        
        !--------------------------------------------------
        !Prepare solutions
        !--------------------------------------------------
        allocate(solution(1:nx,1:ny,Ndata))
        if(Iflag_DIS==0) then       !Steady flow: output base flow: r0, u0, v0, T0, p0, sxx0,sxy0,syy0,qx0,qy0
            do j=1,ny
                do i=1,nx
                    prim = GetPrimary(ctr(i,j)%conVars)
                    solution(i,j,1:3)   = prim(1:3)                                                 !r,u,v
                    solution(i,j,4)     = 1.0/prim(4)                                               !Temperature
                    solution(i,j,5)     = 0.5*solution(i,j,1)*solution(i,j,4)                       !Pressure=0.5*rho*T
                    !Evaluate macro output moments directly from flattened cell-center VDF views.
                    solution(i,j,6:8)   = Get_Stress_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim)      !Shear stress
                    if(IsS2CN2().and..not.IsHighOrderImplicit())then
                        solution(i,j,9:10)=S2NativeOutputHeat(i,j)
                    else
                        solution(i,j,9:10)=GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim)
                    endif
                    solution(i,j,11:14) = ctr(i,j)%res_cell(1:4)                                      !Cell residual
                end do
            end do
        elseif(Iflag_DIS==1) then   !Unsteady flow: output disturbance: r', u', v', T', p', sxx',sxy',syy',qx',qy'
            do j=1,ny
                do i=1,nx
                    !Initial base flow
                    prim0 = GetPrimary(ctr(i,j)%conVars0)
                    solu0(1:3)   = prim0(1:3)                                               !r0,u0,v0
                    solu0(4)     = 1.0/prim0(4)                                             !T0
                    solu0(5)     = 0.5*solu0(1)*solu0(4)                                    !p0=0.5*r0*T0
                    !Evaluate initial macro output moments from flattened initial VDF views.
                    solu0(6:8)   = Get_Stress_Vec(ctr(i,j)%h0_vec,ctr(i,j)%b0_vec,prim0)  !sxx0,sxy0,syy0
                    solu0(9:10)  = GetHeatFlux_Vec(ctr(i,j)%h0_vec,ctr(i,j)%b0_vec,prim0) !qx0,qy0
                    
                    prim = GetPrimary(ctr(i,j)%conVars)
                    solution(i,j,1:3)   = prim(1:3)-solu0(1:3)                                              !r', u', v'
                    solution(i,j,4)     = 1.0/prim(4)-solu0(4)                                              !T'
                    !solution(i,j,5)     = 0.5*solution(i,j,1)*solution(i,j,4)-solu0(5)                      !p'
                    solution(i,j,5)     = 0.5*prim(1)/prim(4)-solu0(5)                      
                    !Subtract initial moments after evaluating the current flattened VDF views.
                    solution(i,j,6:8)   = Get_Stress_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim)-solu0(6:8)   !sxx',sxy',syy'
                    solution(i,j,9:10)  = GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim)-solu0(9:10) !qx',qy'
                    solution(i,j,11:14) = ctr(i,j)%res_cell(1:4)                                              !Cell residual
                end do
            end do
        endif

        !Processor 0 gathers data from other processors
        if(myid==0) then 
            allocate(solution_all(Nx_all,Ny_all,Ndata))
            do pj=0,NP_Y-1 
                do pi=0,NP_X-1
                    np_recv=pj*NP_X+pi
                    if(np_recv.eq.0) then
                        do j1=1,ny
                            do i1=1,nx
                                solution_all(i1,j1,1:Ndata) = solution(i1,j1,1:Ndata)
                            enddo
                        enddo
                    else
                        allocate(buff1(Block_XY(np_recv)%NI,Block_XY(np_recv)%NJ,Ndata))
                        tag=np_recv     !tag equals to block index
                        communication_clock=MPI_Wtime()
                        call MPI_Recv(buff1,Block_XY(np_recv)%NI*Block_XY(np_recv)%NJ*Ndata,MPI_KREAL,np_recv,tag,MPI_COMM_WORLD,Status,ierr)
                        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                        do j1=1,Block_XY(np_recv)%NJ
                            do i1=1,Block_XY(np_recv)%NI
                                i=Block_XY(np_recv)%Origin_I+i1-1
                                j=Block_XY(np_recv)%Origin_J+j1-1
                                solution_all(i,j,1:Ndata)=buff1(i1,j1,1:Ndata)
                            enddo
                        enddo
                        deallocate(buff1)
                    endif
                enddo
            enddo
        else
            allocate(buff2(nx,ny,Ndata))
            tag=IP_Y*NP_X+IP_X      !tag equals to block index
            do j1=1,ny
                do i1=1,nx
                    buff2(i1,j1,1:Ndata)=solution(i1,j1,1:Ndata)
                enddo
            enddo
            communication_clock=MPI_Wtime()
            call MPI_Send(buff2,nx*ny*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,ierr)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            deallocate(buff2)
        endif
        
        deallocate(solution)
        communication_clock=MPI_Wtime()
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        
        !--------------------------------------------------
        !Write to file
        !--------------------------------------------------
        if(myID==0) then
            !-------------------------------------------
            !>Output gener flow field
            !-------------------------------------------
            write(str , *) iter
            !====unformatted file
            !call write_wallboundary
            ! Legacy external-flow binary output needs surf.dat and nonzero Mach scaling.
            ! Native S2 uses the same Tecplot writer below, including stationary walls.
            if(.not.IsS2CN2().and..not.UnifiedPhysicalOperators)then
            allocate(surfx(nx_all),surfy(nx_all),FlowVars_all(nx_all,5))
            open(1,file='surf.dat',form='unformatted')
            read(1) 
            read(1) (surfx(i),i=1,nx_all),(surfy(i),i=1,nx_all),(FlowVars_all(i,1:5),i=1,nx_all)
            close(1)

        open(1,file='grid.grd',form='unformatted')
            write(1) 1
            write(1) nx_all,ny_all+1
            write(1) (0.5*(gridx_all(i,1)+gridx_all(i+1,1)),i=1,nx_all),((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),(0.5*(gridy_all(i,1)+gridy_all(i+1,1)),i=1,nx_all),((gridyc_all(i,j),i=1,nx_all),j=1,ny_all)
			!((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),((gridyc_all(i,j),i=1,nx_all),j=1,ny_all)
            close(1)

        open(1,file='grid_new.grd',form='unformatted')
            write(1) 1
            write(1) nx_all,ny_all+1
            write(1) (0.5*(gridx_all(i,1)+gridx_all(i+1,1)),i=1,nx_all),((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),((0.5*(gridy_all(i,j)+gridy_all(i+1,j)),i=1,nx_all),j=1,ny_all+1)
			!((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),((gridyc_all(i,j),i=1,nx_all),j=1,ny_all)
            close(1)

            open(1,file='uvw2D.dat',form='unformatted')
            write(1) 1
            write(1) nx_all,ny_all+1,1,6
            write(1) (FlowVars_all(i,1),i=1,nx_all),(((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,1),&
            (FlowVars_all(i,2)/Ma*sqrt(2.d0/gamma),i=1,nx_all),(((solution_all(i,j,k)/Ma*sqrt(2.d0/gamma),i=1,nx_all),j=1,ny_all),k=2,2),&
            (FlowVars_all(i,3)/Ma*sqrt(2.d0/gamma),i=1,nx_all),(((solution_all(i,j,k)/Ma*sqrt(2.d0/gamma),i=1,nx_all),j=1,ny_all),k=3,3),&
            (0.d0,i=1,nx_all),((0.d0,i=1,nx_all),j=1,ny_all),&
            (FlowVars_all(i,4),i=1,nx_all),(((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=4,4),&
            (FlowVars_all(i,5)/Ma**2*2.d0/gamma,i=1,nx_all),(((solution_all(i,j,k)/Ma**2*2.d0/gamma,i=1,nx_all),j=1,ny_all),k=5,5)
            !(((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,1),(((solution_all(i,j,k)/Ma*sqrt(2.d0/gamma),i=1,nx_all),j=1,ny_all),k=2,2),&
            !(((solution_all(i,j,k)/Ma*sqrt(2.d0/gamma),i=1,nx_all),j=1,ny_all),k=3,3),((0.d0,i=1,nx_all),j=1,ny_all),&
            !(((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=4,4),(((solution_all(i,j,k)/Ma**2*2.d0/gamma,i=1,nx_all),j=1,ny_all),k=5,5)
            close(1)
            deallocate(surfx,surfy,FlowVars_all)
            endif
            !open(1,file='Flow2d'//trim(fileName)//'_'//trim(adjustl(str))//'.dat',form='unformatted')
            !write(1) 1
            !write(1) nx_all,ny_all
            !write(1) (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,5)
            !close(1)
            
            !open(1,file='Flow2d.dat',form='unformatted')
            !write(1) 1
            !write(1) nx_all,ny_all
            !write(1) (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,5)
            !close(1)
            
            !====formatted file
            !history file: cell-centered form
            !open(1,file='Flow2dC'//trim(fileName)//'_'//trim(adjustl(str))//'.plt')
            open(1,file='Flow2dC.plt')
            write(1,*) "VARIABLES = x, y, r, u, v, T, p, sxx, sxy, syy, qx, qy"
            write(1,*) "ZONE  I = ",nx_all-1+2,", J = ",ny_all-1+2,"DATAPACKING=BLOCK, VARLOCATION=([3-12]=CELLCENTERED)"
            write(1,*) ((gridx_all(i,j),i=1,nx_all+1),j=1,ny_all+1),((gridy_all(i,j),i=1,nx_all+1),j=1,ny_all+1), &
                       (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,10)
            close(1)

            !final result (restart file): cell center value and point form
            !open(1,file='Flow2d.plt')
            open(1,file='Flow2d_'//trim(adjustl(str))//'.plt')
            write(1,*) "VARIABLES = x, y, r, u, v, T, p, sxx, sxy, syy, qx, qy"
            write(1,'(a20,I8,a20,I8,a20)') "ZONE  I = ",nx_all-1+1,", J = ",ny_all-1+1,", DATAPACKING=BLOCK"
            write(1,*) ((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),((gridyc_all(i,j),i=1,nx_all),j=1,ny_all), &
                       (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,10)
            close(1)

            !Per-cell residual cloud (4 conservative components)
            open(1,file='Residual2d_'//trim(adjustl(str))//'.plt')
            write(1,*) "VARIABLES = x, y, resRho, resRhoU, resRhoV, resE, logResRho, logResRhoU, logResRhoV, logResE"
            write(1,'(a20,I8,a20,I8,a20)') "ZONE  I = ",nx_all,", J = ",ny_all,", DATAPACKING=BLOCK"
            write(1,*) ((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),((gridyc_all(i,j),i=1,nx_all),j=1,ny_all), &
                       (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=11,14), &
                       (((log10(max(solution_all(i,j,k),SMV)),i=1,nx_all),j=1,ny_all),k=11,14)
            close(1)
            
            !-------------------------------------------
            !>Output specified line for some geometries
            !-------------------------------------------
            !=====Write wall-normal profile for Couette flow
            !mid = (1+nx_all)/2
            !open(1,file='profile.plt')
            !write(1,*) 'variables=y,density,u,v,T,p,txx,txy,tyy,qx,qy'
            !do j=1,ny_all
            !    write(1,'(11e24.16)') gridyc_all(mid,j),(solution_all(mid,j,i),i=1,10)
            !enddo
            !close(1)
            
            !=====Write central line profile for Cavity
            !mid = (1+nx_all)/2
            !open(unit=1,file='Y-U.plt',status="replace",action="write")
            !do j=1,ny_all
            !    write(1,*)  gridyc_all(mid,j),solution_all(mid,j,2)/0.15
            !end do
            !close(1)
            !
            !mid = (1+ny_all)/2
            !open(unit=1,file='X-V.plt',status="replace",action="write")
            !do i=1,nx_all
            !    write(1,*)  gridxc_all(i,mid),solution_all(i,mid,3)/0.15
            !end do
            !close(1)
            
            deallocate(solution_all)
        endif
        
    end subroutine write_flow
    
    !-------------------------------------------------------------------
    !>Output velocity distribution function for a cell
    !-------------------------------------------------------------------
    
    
    !-------------------------------------------------------------------
    !>Output surface surface flowVars and aerodynamics
    !-------------------------------------------------------------------
    !The surface quantities has been computed in the CalFluxBoundary, here only output BoundaryFace(3)
    subroutine write_wallboundary(side_arg)
        integer,intent(in),optional::side_arg
        integer::side,face_j,global_j
        character(len=32)::surf_name,surf_binary,coefficient_name
        integer(KINT)                                   :: i,pi,i1,np_recv,Ndata
        real(KREAL)                                     :: pinf
        real(KREAL), allocatable, dimension(:,:)        :: FlowVars,FlowVars_all         !rho,vn,vt,T,p,shxx,shxy,shyy,qx,qy (by interface distribution) 
        real(KREAL), allocatable, dimension(:,:)        :: Aerodynamics,Aerodynamics_all !p_wall, tau_wall, qn_wall (by particle-surface collision)
        real(KREAL), allocatable, dimension(:,:)        :: buff1,buff2
        
        side=1
        if(present(side_arg))side=side_arg
        face_j=1;global_j=1;surf_name='Surf.plt';surf_binary='surf.dat';coefficient_name='Surf_CpCfCh.plt'
        if(side==2)then
            face_j=ny+1;global_j=ny_all+1;surf_name='Surf_top.plt';surf_binary='surf_top.dat'
            coefficient_name='Surf_top_CpCfCh.plt'
        endif
        !Only the bottom boundary is post-processed; other boundary can be treated in the same way
        if(bc_type(side+2)==BCconst_wall_isothermal) then    !bottom boundary
            !assign surface quantities
            allocate(FlowVars(nx,5),Aerodynamics(nx,3))
            if(Iflag_DIS==0) then
                do i=1,nx
                    FlowVars(i,1:5)        = hface(i,face_j)%FlowVars(1:5)
                    Aerodynamics(i,1:3)     = hface(i,face_j)%Aerodynamics(1:3)
                end do
            elseif(Iflag_DIS==1) then
                do i=1,nx
                    FlowVars(i,1:5)        = hface(i,face_j)%FlowVars(1:5)-hface(i,face_j)%FlowVars0(1:5)
                    Aerodynamics(i,1:3)     = hface(i,face_j)%Aerodynamics(1:3)-hface(i,face_j)%Aerodynamics0(1:3)
                end do
            endif
            
            !Processor 0 gathers all data
            Ndata = 8           !MPI_Sendrecv: FlowVars(1:10), Aerodynamics(1:3)
            if(myid==0)then
                allocate(FlowVars_all(nx_all,5),Aerodynamics_all(nx_all,3))
                do pi=0,NP_X-1
                    np_recv=pi
                    if(side==2)np_recv=(NP_Y-1)*NP_X+pi
                    allocate(buff1(Block_XY(np_recv)%NI,Ndata))
                    if(np_recv==0)then
                        buff1(:,1:5)=FlowVars;buff1(:,6:8)=Aerodynamics
                    else
                        communication_clock=MPI_Wtime()
                        call MPI_Recv(buff1,size(buff1),MPI_KREAL,np_recv,8000+side,MPI_COMM_WORLD,Status,ierr)
                        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                    endif
                    do i1=1,Block_XY(np_recv)%NI
                        i=Block_XY(np_recv)%Origin_I+i1-1
                        FlowVars_all(i,:)=buff1(i1,1:5);Aerodynamics_all(i,:)=buff1(i1,6:8)
                    enddo
                    deallocate(buff1)
                enddo
            elseif((side==1.and.IP_Y==0).or.(side==2.and.IP_Y==NP_Y-1))then
                allocate(buff2(nx,Ndata));buff2(:,1:5)=FlowVars;buff2(:,6:8)=Aerodynamics
                communication_clock=MPI_Wtime()
                call MPI_Send(buff2,size(buff2),MPI_KREAL,0,8000+side,MPI_COMM_WORLD,ierr)
                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                deallocate(buff2)
            endif
            deallocate(FlowVars,Aerodynamics)
            
            communication_clock=MPI_Wtime()
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                
            !write to file
            if(myID==0) then
                !!Output surface flow variables
                !open(1,file='Surf_FlowVars.plt')
                !write(1,'(a)') 'variables=x,y,density,u,v,T,p'
                !do i=1,nx_all
                !    write(1,'(7e16.6)') 0.5*(gridx_all(i,global_j)+gridx_all(i+1,global_j)),0.5*(gridy_all(i,global_j)+gridy_all(i+1,global_j)),FlowVars_all(i,1:5)
                !end do   
                !close(1)
                !    
                !!Output surface aerodynamics
                !open(1,file='Surf_Aerodynamics.plt')
                !write(1,'(a)') 'variables=x,y,pnn,pnt,qn'
                !do i=1,nx_all
                !    write(1,'(5e16.6)') 0.5*(gridx_all(i,global_j)+gridx_all(i+1,global_j)),0.5*(gridy_all(i,global_j)+gridy_all(i+1,global_j)),Aerodynamics_all(i,1:3)
                !enddo
                !close(1)
                
                !Output surface quantities
                open(1,file=trim(surf_name))
                write(1,'(a)') 'variables=x,y,density,u,v,T,p,pnn_w,pnt_w,qn_w'
                do i=1,nx_all
                    write(1,'(10e16.6)') 0.5*(gridx_all(i,global_j)+gridx_all(i+1,global_j)),0.5*(gridy_all(i,global_j)+gridy_all(i+1,global_j)), &
                                        FlowVars_all(i,1:5), Aerodynamics_all(i,1:3)
                end do   
                close(1)
                
        open(1,file=trim(surf_binary),form='unformatted')
                    write(1) nx_all,7
                    write(1) (0.5*(gridx_all(i,global_j)+gridx_all(i+1,global_j)),i=1,nx_all),(0.5*(gridy_all(i,global_j)+gridy_all(i+1,global_j)),i=1,nx_all),(FlowVars_all(i,1:5),i=1,nx_all)
                    close(1)

    
                !Output surface coefficients for external flow
                pinf = 0.5*Prim_inf(1)/Prim_inf(4)          !p=0.5*rho*T
                if(Ma/=0.0_KREAL)then
                open(1,file=trim(coefficient_name))
                write(1,'(a)') 'variables=x,y,Cp,Cf,Ch'
                do i=1,nx_all
                    write(1,'(5es24.16)') 0.5*(gridx_all(i,global_j)+gridx_all(i+1,global_j)),  0.5*(gridy_all(i,global_j)+gridy_all(i+1,global_j)),&
                                        (Aerodynamics_all(i,1)-pinf)/(0.25*gamma*Ma**2), &  !surface pressure coefficient
                                        Aerodynamics_all(i,2)/(0.25*gamma*Ma**2), &         !skin friction coefficient
                                        -Aerodynamics_all(i,3)*2.0*(2.0/gamma)**1.5/Ma**3   !surface heat flux coefficient
                enddo
                close(1)
                endif
                deallocate(FlowVars_all,Aerodynamics_all)
            endif
        endif
    
    end subroutine write_wallboundary
    

    !-------------------------------------------------------------------
    !>Write restart file, including: grid2d.grd, flow2d.dat, and vdf2d.dat (single-block)
    !>Logic Updated: 
    !>  - Main VDF: Unstruct -> vdf2d_unstruct_full.dat (h only)
    !>  - Pert VDF: Unstruct -> vdf2d_unstruct_full_pert.dat (h' only)
    !-------------------------------------------------------------------
    subroutine write_restart(flow_path_arg,vdf_path_arg,grid_path_arg)
        character(len=*),intent(in),optional::flow_path_arg,vdf_path_arg,grid_path_arg
        character(len=1024)::flow_path,vdf_path,grid_path
        integer(KINT)                                   :: i,j,k,l,m,pi,pj,np_recv,i1,j1,Ndata,tag
        integer(KINT)                                   :: flat_k !Flattened velocity-space index for VDF output copies
        real(KREAL), allocatable, dimension(:,:,:)      :: solution,solution_all
        real(KREAL), allocatable, dimension(:,:,:)      :: buff1,buff2
        real(KREAL), allocatable, dimension(:,:,:,:)    :: solution_vdf,solution_vdf_all !Structured restart VDF cache with flat Nvel velocity dimension
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff1_vdf,buff2_vdf !Structured restart MPI buffers with flat Nvel velocity dimension
        ! perturbation VDF containers (h' and b')
        real(KREAL), allocatable, dimension(:,:,:,:)    :: solution_vdf_pert,solution_vdf_pert_all !Structured perturbation VDF cache with flat Nvel velocity dimension
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff1_vdf_pert,buff2_vdf_pert !Structured perturbation MPI buffers with flat Nvel velocity dimension
        ! Unstructured Grid VDF Variables
        ! For Main VDF (h)
        real(KREAL), allocatable, dimension(:,:,:,:) :: buff_vdf_local  !Unstructured VDF cache stored as (nx, ny, Nvel, 2)
        real(KREAL), allocatable, dimension(:,:,:,:) :: buff_vdf_all    !Global unstructured VDF cache stored as (Nx_all, Ny_all, Nvel, 2)
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff_vdf_pert_local  !Local perturbation VDF cache stored as (NI, NJ, Nvel, 2)
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff_vdf_pert_all    !Global perturbation VDF cache stored as (Nx_all, Ny_all, Nvel, 2)
        real(KREAL), allocatable, dimension(:,:,:,:)      :: recv_buffer   ! Shared temp buffer
            
        flow_path='flow2d.dat';vdf_path='vdf2d.dat';grid_path='grid2d.grd'
        if(present(flow_path_arg))flow_path=flow_path_arg
        if(present(vdf_path_arg))vdf_path=vdf_path_arg
        if(present(grid_path_arg))grid_path=grid_path_arg
        !---------------------------------------------------------------------
        !>Write grid: grid2d.grd
        !---This file will not be read when restart, but for visualization and interpolation
        !---------------------------------------------------------------------
        if(myid==0) then
            open(1,file=trim(grid_path),form='unformatted')
            write(1) 1
            write(1) nx_all,ny_all
            write(1) ((gridxc_all(i,j),i=1,nx_all),j=1,ny_all),((gridyc_all(i,j),i=1,nx_all),j=1,ny_all)
            close(1)
        endif

            
        !---------------------------------------------------------------------
        !>Write macro-conserved variables: flow2d.dat-->ctr(:,:)%conVars(1:4)
        !---------------------------------------------------------------------
        Ndata = 4
        
        !--prepare solution
        allocate(solution(1:nx,1:ny,Ndata))
        do j=1,ny
            do i=1,nx
                solution(i,j,1:4)   = ctr(i,j)%conVars
            end do
        end do

        !--processor 0 gathers data from other processors
        if(myid==0) then 
            allocate(solution_all(Nx_all,Ny_all,Ndata))
            do pj=0,NP_Y-1 
                do pi=0,NP_X-1
                    np_recv=pj*NP_X+pi
                    if(np_recv.eq.0) then
                        do j1=1,ny
                            do i1=1,nx
                                solution_all(i1,j1,1:Ndata) = solution(i1,j1,1:Ndata)
                            enddo
                        enddo
                    else
                        allocate(buff1(Block_XY(np_recv)%NI,Block_XY(np_recv)%NJ,Ndata))
                        tag=np_recv     !tag equals to block index
                        communication_clock=MPI_Wtime()
                        call MPI_Recv(buff1,Block_XY(np_recv)%NI*Block_XY(np_recv)%NJ*Ndata,MPI_KREAL,np_recv,tag,MPI_COMM_WORLD,Status,ierr)
                        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                        do j1=1,Block_XY(np_recv)%NJ
                            do i1=1,Block_XY(np_recv)%NI
                                i=Block_XY(np_recv)%Origin_I+i1-1
                                j=Block_XY(np_recv)%Origin_J+j1-1
                                solution_all(i,j,1:Ndata)=buff1(i1,j1,1:Ndata)
                            enddo
                        enddo
                        deallocate(buff1)
                    endif
                enddo
            enddo
        else
            allocate(buff2(nx,ny,Ndata))
            tag=IP_Y*NP_X+IP_X      !tag equals to block index
            do j1=1,ny
                do i1=1,nx
                    buff2(i1,j1,1:Ndata)=solution(i1,j1,1:Ndata)
                enddo
            enddo
            communication_clock=MPI_Wtime()
            call MPI_Send(buff2,nx*ny*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,ierr)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            deallocate(buff2)
        endif
        deallocate(solution)
        communication_clock=MPI_Wtime()
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
        
        !--write to file
        if(myid==0) then
            open(1,file=trim(flow_path),form='unformatted')
            write(1) 1
            write(1) nx_all,ny_all,Ndata
            write(1) (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,Ndata)
            close(1)
            deallocate(solution_all)
        endif
            
            
        !=====================================================================
        !> SECTION 1: Main VDF Output (h)
        !=====================================================================
        
        !>>> BRANCH A: Unstructured (Optimized) <<<
        if (Quadrature_method == Quad_Unstruct) then
            
            allocate(buff_vdf_local(nx, ny, Nvel, 2))
            do j = 1, ny
                do i = 1, nx
                    !Copy unstructured VDFs through flattened cell-center views.
                    do k=1,Nvel
                        buff_vdf_local(i, j, k, 1) = ctr(i, j)%h_vec(k)
                        buff_vdf_local(i, j, k, 2) = ctr(i, j)%b_vec(k)
                    enddo
                end do
            end do
            if (myid == 0) then
                allocate(buff_vdf_all(Nx_all, Ny_all, Nvel, 2))
                do pj = 0, NP_Y - 1
                    do pi = 0, NP_X - 1
                        np_recv = pj * NP_X + pi
                        if (np_recv == 0) then
                            do j1 = 1, ny
                                do i1 = 1, nx
                                    buff_vdf_all(i1, j1, :, :) = buff_vdf_local(i1, j1, :, :)
                                end do
                            end do
                        else
                            allocate(recv_buffer(Block_XY(np_recv)%NI, Block_XY(np_recv)%NJ, Nvel, 2))
                            tag = np_recv * 20 
                            communication_clock=MPI_Wtime()
                            call MPI_Recv(recv_buffer, Block_XY(np_recv)%NI * Block_XY(np_recv)%NJ * Nvel * 2, &
                                          MPI_KREAL, np_recv, tag, MPI_COMM_WORLD, Status, ierr)
                            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                            do j1 = 1, Block_XY(np_recv)%NJ
                                do i1 = 1, Block_XY(np_recv)%NI
                                    i = Block_XY(np_recv)%Origin_I + i1 - 1
                                    j = Block_XY(np_recv)%Origin_J + j1 - 1
                                    buff_vdf_all(i, j, :, :) = recv_buffer(i1, j1, :, :)
                                end do
                            end do
                            deallocate(recv_buffer)
                        endif
                    end do
                end do
            else
                tag = (IP_Y * NP_X + IP_X) * 20
                communication_clock=MPI_Wtime()
                call MPI_Send(buff_vdf_local, nx * ny * Nvel * 2, MPI_KREAL, 0, tag, MPI_COMM_WORLD, ierr)
                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            endif
            deallocate(buff_vdf_local)
            communication_clock=MPI_Wtime()
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock

            if (myid == 0) then
                open(1, file='vdf2d_unstruct_full.dat', form='unformatted')
                write(1) Nx_all, Ny_all, Nvel, 1
                write(1) (uSpace(k), k=1,Nvel), (vSpace(k), k=1,Nvel), (weight(k), k=1,Nvel)
                write(1) buff_vdf_all
                close(1)
                deallocate(buff_vdf_all)
            endif

        !>>> BRANCH B: Structured (Standard) <<<
        else 
            Ndata = 2
            allocate(solution_vdf(nx,ny,Nvel,Ndata))
            do j=1,ny
                do i=1,nx
                    !Copy structured VDFs through flattened cell-center views.
                    do flat_k=1,Nvel
                        solution_vdf(i,j,flat_k,1) = ctr(i,j)%h_vec(flat_k)
                        solution_vdf(i,j,flat_k,2) = ctr(i,j)%b_vec(flat_k)
                    enddo
                end do
            end do
            
            if(myid==0) then 
                allocate(solution_vdf_all(Nx_all,Ny_all,Nvel,Ndata))
                do pj=0,NP_Y-1 
                    do pi=0,NP_X-1
                        np_recv=pj*NP_X+pi
                        if(np_recv.eq.0) then
                            do j1=1,ny
                                do i1=1,nx
                                    solution_vdf_all(i1,j1,1:Nvel,1:Ndata) = solution_vdf(i1,j1,1:Nvel,1:Ndata)
                                enddo
                            enddo
                        else
                            allocate(buff1_vdf(Block_XY(np_recv)%NI,Block_XY(np_recv)%NJ,Nvel,Ndata))
                            tag=np_recv*10 
                            communication_clock=MPI_Wtime()
                            call MPI_Recv(buff1_vdf,Block_XY(np_recv)%NI*Block_XY(np_recv)%NJ*Nvel*Ndata, &
                                          MPI_KREAL,np_recv,tag,MPI_COMM_WORLD,Status,ierr)
                            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                            do j1=1,Block_XY(np_recv)%NJ
                                do i1=1,Block_XY(np_recv)%NI
                                    i=Block_XY(np_recv)%Origin_I+i1-1
                                    j=Block_XY(np_recv)%Origin_J+j1-1
                                    solution_vdf_all(i,j,1:Nvel,1:Ndata)=buff1_vdf(i1,j1,1:Nvel,1:Ndata)
                                enddo
                            enddo
                            deallocate(buff1_vdf)
                        endif
                    enddo
                enddo
            else
                allocate(buff2_vdf(nx,ny,Nvel,Ndata))
                tag=(IP_Y*NP_X+IP_X)*10 
                do j1=1,ny
                    do i1=1,nx
                        buff2_vdf(i1,j1,1:Nvel,1:Ndata)=solution_vdf(i1,j1,1:Nvel,1:Ndata)
                    enddo
                enddo
                communication_clock=MPI_Wtime()
                call MPI_Send(buff2_vdf,nx*ny*Nvel*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,ierr)
                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                deallocate(buff2_vdf)
            endif
            deallocate(solution_vdf)
            communication_clock=MPI_Wtime()
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
            
            if(myid==0) then
                open(1,file=trim(vdf_path),form='unformatted')
                write(1) 1
                write(1) nx_all,ny_all,Nvel*Ndata
                write(1) ((((solution_vdf_all(i,j,flat_k,k),i=1,nx_all),j=1,ny_all),flat_k=1,Nvel),k=1,Ndata)
                close(1)
                deallocate(solution_vdf_all)
            endif
        endif    
            
        !=====================================================================
        !> SECTION 2: Perturbation Output (Only if Iflag_DIS == 1)
        !=====================================================================
        if (Iflag_DIS == 1) then
            
            !>>> BRANCH A: Unstructured Perturbation (h' only) <<<
            if (Quadrature_method == Quad_Unstruct) then
                
                ! 1. Calculate Local Perturbation h' = h - h0
                
                Ndata = 2

                allocate(buff_vdf_pert_local(nx, ny, Nvel, Ndata))
                do j = 1, ny
                    do i = 1, nx
                        !Calculate unstructured perturbation VDFs through flattened cell-center views.
                        do k=1,Nvel
                            buff_vdf_pert_local(i, j, k, 1) = ctr(i, j)%h_vec(k)-ctr(i, j)%h0_vec(k)
                            buff_vdf_pert_local(i, j, k, 2) = ctr(i, j)%b_vec(k)-ctr(i, j)%b0_vec(k)
                        enddo
                    end do
                end do

                ! 2. MPI Gather (Use tag + 1 logic to avoid collision with Main VDF)
                if (myid == 0) then
                    allocate(buff_vdf_pert_all(Nx_all, Ny_all, Nvel, Ndata))
                    
                    do pj = 0, NP_Y - 1
                        do pi = 0, NP_X - 1
                            np_recv = pj * NP_X + pi
                            
                            if (np_recv == 0) then
                                do j1 = 1, ny
                                    do i1 = 1, nx
                                        buff_vdf_pert_all(i1, j1, :, :) = buff_vdf_pert_local(i1, j1, :, :)
                                    end do
                                end do
                            else
                                allocate(recv_buffer(Block_XY(np_recv)%NI, Block_XY(np_recv)%NJ, Nvel, Ndata))
                                tag = np_recv * 20 + 1   !Base tag + 1 for perturbation
                                
                                communication_clock=MPI_Wtime()
                                call MPI_Recv(recv_buffer, Block_XY(np_recv)%NI * Block_XY(np_recv)%NJ * Nvel * Ndata, &
                                              MPI_KREAL, np_recv, tag, MPI_COMM_WORLD, Status, ierr)
                                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                                
                                do j1 = 1, Block_XY(np_recv)%NJ
                                    do i1 = 1, Block_XY(np_recv)%NI
                                        i = Block_XY(np_recv)%Origin_I + i1 - 1
                                        j = Block_XY(np_recv)%Origin_J + j1 - 1
                                        buff_vdf_pert_all(i, j, :, :) = recv_buffer(i1, j1, :, :)
                                    end do
                                end do
                                deallocate(recv_buffer)
                            endif
                        end do
                    end do
                else
                    ! Send
                    tag = (IP_Y * NP_X + IP_X) * 20 + 1
                    communication_clock=MPI_Wtime()
                    call MPI_Send(buff_vdf_pert_local, nx * ny * Nvel * Ndata, MPI_KREAL, 0, tag, MPI_COMM_WORLD, ierr)
                    communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                endif
                
                deallocate(buff_vdf_pert_local)
                communication_clock=MPI_Wtime()
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock

                ! 3. Write New Format Perturbation File
                if (myid == 0) then
                    open(1, file='vdf2d_unstruct_full_pert.dat', form='unformatted')
                    ! Header
                    write(1) Nx_all, Ny_all, Nvel, 1
                    ! Coordinates (Useful for standalone plotting)
                    write(1) (uSpace(k), k=1,Nvel), (vSpace(k), k=1,Nvel), (weight(k), k=1,Nvel)
                    ! Data h'
                    write(1) buff_vdf_pert_all
                    close(1)
                    deallocate(buff_vdf_pert_all)
                endif

            !>>> BRANCH B: Structured Perturbation (Standard) <<<
            else 
                Ndata = 2
                allocate(solution_vdf_pert(nx,ny,Nvel,Ndata))
                do j=1,ny
                    do i=1,nx
                        !Calculate structured perturbation VDFs through flattened cell-center views.
                        do flat_k=1,Nvel
                            solution_vdf_pert(i,j,flat_k,1) = ctr(i,j)%h_vec(flat_k)-ctr(i,j)%h0_vec(flat_k)
                            solution_vdf_pert(i,j,flat_k,2) = ctr(i,j)%b_vec(flat_k)-ctr(i,j)%b0_vec(flat_k)
                        enddo
                    end do
                end do

                if(myid==0) then 
                    allocate(solution_vdf_pert_all(Nx_all,Ny_all,Nvel,Ndata))
                    do pj=0,NP_Y-1 
                        do pi=0,NP_X-1
                            np_recv=pj*NP_X+pi
                            if(np_recv.eq.0) then
                                do j1=1,ny
                                    do i1=1,nx
                                        solution_vdf_pert_all(i1,j1,1:Nvel,1:Ndata) = solution_vdf_pert(i1,j1,1:Nvel,1:Ndata)
                                    enddo
                                enddo
                            else
                                allocate(buff1_vdf_pert(Block_XY(np_recv)%NI,Block_XY(np_recv)%NJ,Nvel,Ndata))
                                tag=np_recv*10 + 1 
                                communication_clock=MPI_Wtime()
                                call MPI_Recv(buff1_vdf_pert,Block_XY(np_recv)%NI*Block_XY(np_recv)%NJ*Nvel*Ndata, &
                                              MPI_KREAL,np_recv,tag,MPI_COMM_WORLD,Status,ierr)
                                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                                do j1=1,Block_XY(np_recv)%NJ
                                    do i1=1,Block_XY(np_recv)%NI
                                        i=Block_XY(np_recv)%Origin_I+i1-1
                                        j=Block_XY(np_recv)%Origin_J+j1-1
                                        solution_vdf_pert_all(i,j,1:Nvel,1:Ndata)=buff1_vdf_pert(i1,j1,1:Nvel,1:Ndata)
                                    enddo
                                enddo
                                deallocate(buff1_vdf_pert)
                            endif
                        enddo
                    enddo
                else
                    allocate(buff2_vdf_pert(nx,ny,Nvel,Ndata))
                    tag=(IP_Y*NP_X+IP_X)*10 + 1 
                    do j1=1,ny
                        do i1=1,nx
                            buff2_vdf_pert(i1,j1,1:Nvel,1:Ndata)=solution_vdf_pert(i1,j1,1:Nvel,1:Ndata)
                        enddo
                    enddo
                    communication_clock=MPI_Wtime()
                    call MPI_Send(buff2_vdf_pert,nx*ny*Nvel*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,ierr)
                    communication_seconds=communication_seconds+MPI_Wtime()-communication_clock
                    deallocate(buff2_vdf_pert)
                endif
                deallocate(solution_vdf_pert)
                communication_clock=MPI_Wtime()
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                communication_seconds=communication_seconds+MPI_Wtime()-communication_clock

                if(myid==0) then
                    open(1,file='vdf2d_pert.dat',form='unformatted')
                    write(1) 1
                    write(1) nx_all,ny_all,Nvel*Ndata
                    write(1) ((((solution_vdf_pert_all(i,j,flat_k,k),i=1,nx_all),j=1,ny_all),flat_k=1,Nvel),k=1,Ndata)
                    close(1)
                    deallocate(solution_vdf_pert_all)
                endif
            endif
        endif
        
    end subroutine write_restart
    

end module Output

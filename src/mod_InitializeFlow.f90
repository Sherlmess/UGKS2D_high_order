!--------------------------------------------------
!>Initialization of mesh and intial flow field
!--------------------------------------------------
module InitializeFlow
    use flow_para
    use control_para, only: IsImplicitSteady, IsImplicitUnsteady, PeriodicComparisonProfile
    use mod_S2D1, only: IsS2CN2,IsS2Steady,S2InitVelocity => init_velocity
    use mod_S2NativeAdapter, only: S2NativeInitializeState,S2NativeReadMetadata
    use mod_S2SteadyIteration, only: InitializeS2Steady
    use Flow_DIS
    use tools
    use Boundary
    use mod_ImplicitBoundary, only: BuildStructuredMirrorMaps, ValidateImplicitMaxwellFaces
    use mod_ClosedCouetteMass, only: ReadMassRestart,InitializeClosedMass
    use mod_ImplicitMicro, only: EnforceAcceptedMomentCompatibility
    use mod_SteadyImplicit, only: InitializeSteadyImplicit
    use mod_UnsteadyImplicit, only: AllocateUnsteadyImplicitState, FillFringeAtStep
    use mpi
    use mpi_wrapper
    implicit none

contains
    !--------------------------------------------------
    !>Main initialization subroutine
    !--------------------------------------------------
    subroutine InitFlow()
        call InitAllocation()                   !Allocate arrays
        call InitMesh()                         !Initialize mesh
        call InitVelocity()                     !Initialize uSpace, vSpace and weights
        if(IsImplicitSteady() .or. IsImplicitUnsteady()) call ValidateImplicitMaxwellFaces()
        call InitFlowField()                    !Set the initial value
        if(ClosedCouetteMassRescale.and.Iflag_InitFlow==Init_Restart)call ReadMassRestart()
        call Init_fringe_x()                    !Initialize fringe function
        if(bc_type(1)==BCconst_inflow .or. &
           (Iflag_DIS==1 .and. DISBC_type==DISBC_inflow_TS)) then !Initialize inflow/TS boundary data
            call Init_inflow()
        endif
        !Apply physical boundaries only after their inflow data have been allocated.
        call Set_Ghost_PhysBoundary()
        call MPIUpdate_CellConVars()
        call MPIUpdate_CellVDF()
        call MPIUpdate_Conner_conVars()
        call MPIUpdate_Conner_VDF()
        if(IsImplicitSteady() .or. IsImplicitUnsteady().or.UnifiedPhysicalOperators) call InitializeSteadyImplicit()
        if((IsImplicitSteady().or.UnifiedPhysicalOperators).and.Implicit_InitHaloResync)then
            call MPIUpdate_CellConVars();call MPIUpdate_CellVDF()
            call MPIUpdate_Conner_conVars();call MPIUpdate_Conner_VDF()
        endif
        if(ClosedCouetteMassRescale)call InitializeClosedMass()
        if(IsS2Steady()) call InitializeS2Steady()
        if(IsImplicitUnsteady()) call AllocateUnsteadyImplicitState()
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine InitFlow

    
    !--------------------------------------------------
    !>Allocate arrays
    !--------------------------------------------------
    subroutine InitAllocation()
        integer(KINT)                                   :: i,j
        
        !------------------------------------------------------
        !>global mesh and flow field
        !------------------------------------------------------
        if(myID==0) then
            allocate(gridx_all(nx_all+1,ny_all+1),gridy_all(nx_all+1,ny_all+1))
            allocate(gridxc_all(nx_all,ny_all),gridyc_all(nx_all,ny_all))
        endif
        
        !------------------------------------------------------
        !>flow field
        !------------------------------------------------------
        allocate(gridx(1:nx+1,1:ny+1),gridy(1:nx+1,1:ny+1))         !geometry mesh
        allocate(ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))            !Cell center
        allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))             !Vertical and horizontal interfaces
        !Store macro residual pieces separately so future implicit/high-order paths can reuse R(U).
        allocate(res_macro_flux_store(4,1:nx,1:ny))
        allocate(res_macro_store(4,1:nx,1:ny))
        res_macro_flux_store = 0.0_KREAL
        res_macro_store = 0.0_KREAL
        do j=1-LAP,ny+LAP
            do i=1-LAP,nx+LAP
                ctr(i,j)%res_cell = 0.0_KREAL
            end do
        end do
        
        !------------------------------------------------------
        !>velocity space
        !------------------------------------------------------
        !Velocity-space arrays are stored as pure 1D buffers with k=i+(j-1)*Nu.
        Nvel = Nu*Nv
        allocate(uSpace(Nvel))
        allocate(vSpace(Nvel))
        allocate(weight(Nvel))
        allocate(uVel(Nvel))
        allocate(vVel(Nvel))
        allocate(wVel(Nvel))

        !Store micro face-flux residual pieces after Nvel is known.
        allocate(res_h_flux_store(Nvel,1:nx,1:ny))
        allocate(res_b_flux_store(Nvel,1:nx,1:ny))
        res_h_flux_store = 0.0_KREAL
        res_b_flux_store = 0.0_KREAL
        
        !------------------------------------------------------
        !>velocity distribution function
        !------------------------------------------------------
        !Cell center micro variables use pure Nvel storage; file I/O keeps Nu,Nv loops only at the edges.
        allocate(cell_h_store(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_b_store(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_sh_store(Nvel,2,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_sb_store(Nvel,2,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_h0_store(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_b0_store(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_res_h0_store(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(cell_res_b0_store(Nvel,1-LAP:nx+LAP,1-LAP:ny+LAP))
        ! Wall-exterior slots are not physical states, but complete R6 restart
        ! snapshots include their storage. Initialize it deterministically.
        if(UnifiedPhysicalOperators)then
            cell_h_store=0;cell_b_store=0;cell_h0_store=0;cell_b0_store=0
            cell_res_h0_store=0;cell_res_b0_store=0
            do j=1-LAP,ny+LAP;do i=1-LAP,nx+LAP
                ctr(i,j)%conVars=0;ctr(i,j)%conVars0=0;ctr(i,j)%sw=0
                ctr(i,j)%res_conVars0=0;ctr(i,j)%dt_local=dt_global
            enddo;enddo
        endif
        cell_sh_store = 0.0_KREAL
        cell_sb_store = 0.0_KREAL

        !Bind each cell directly to pure 1D velocity-space backing arrays.
        do j=1-LAP,ny+LAP
            do i=1-LAP,nx+LAP
                ctr(i,j)%h_vec(1:Nvel)      => cell_h_store(:,i,j)
                ctr(i,j)%b_vec(1:Nvel)      => cell_b_store(:,i,j)
                ctr(i,j)%sh_vec(1:Nvel,1:2) => cell_sh_store(:,:,i,j)
                ctr(i,j)%sb_vec(1:Nvel,1:2) => cell_sb_store(:,:,i,j)
                ctr(i,j)%h0_vec(1:Nvel)     => cell_h0_store(:,i,j)
                ctr(i,j)%b0_vec(1:Nvel)     => cell_b0_store(:,i,j)
                ctr(i,j)%res_h0_vec(1:Nvel) => cell_res_h0_store(:,i,j)
                ctr(i,j)%res_b0_vec(1:Nvel) => cell_res_b0_store(:,i,j)
            end do
        end do
        call CheckCellMicroPointers()

        !Cell interface micro fluxes use pure Nvel backing arrays.
        allocate(v_l_flux_h_store(Nvel,1:nx+1,1:ny))
        allocate(v_l_flux_b_store(Nvel,1:nx+1,1:ny))
        allocate(v_r_flux_h_store(Nvel,1:nx+1,1:ny))
        allocate(v_r_flux_b_store(Nvel,1:nx+1,1:ny))
        allocate(h_l_flux_h_store(Nvel,1:nx,1:ny+1))
        allocate(h_l_flux_b_store(Nvel,1:nx,1:ny+1))
        allocate(h_r_flux_h_store(Nvel,1:nx,1:ny+1))
        allocate(h_r_flux_b_store(Nvel,1:nx,1:ny+1))

        !Bind vertical-face flux arrays directly to 1D Nvel views.
        do j=1,ny
            do i=1,nx+1
                vface(i,j)%L_flux_h_vec(1:Nvel) => v_l_flux_h_store(:,i,j)
                vface(i,j)%L_flux_b_vec(1:Nvel) => v_l_flux_b_store(:,i,j)
                vface(i,j)%R_flux_h_vec(1:Nvel) => v_r_flux_h_store(:,i,j)
                vface(i,j)%R_flux_b_vec(1:Nvel) => v_r_flux_b_store(:,i,j)
            end do
        end do

        !Bind horizontal-face flux arrays directly to 1D Nvel views.
        do j=1,ny+1
            do i=1,nx
                hface(i,j)%L_flux_h_vec(1:Nvel) => h_l_flux_h_store(:,i,j)
                hface(i,j)%L_flux_b_vec(1:Nvel) => h_l_flux_b_store(:,i,j)
                hface(i,j)%R_flux_h_vec(1:Nvel) => h_r_flux_h_store(:,i,j)
                hface(i,j)%R_flux_b_vec(1:Nvel) => h_r_flux_b_store(:,i,j)
            end do
        end do
        call CheckFaceFluxPointers()

    end subroutine InitAllocation


    !--------------------------------------------------
    !>Verify cell micro pointers are associated with contiguous backing arrays.
    !--------------------------------------------------
    subroutine CheckCellMicroPointers()
        integer(KINT)                                   :: i,j

        do j=1-LAP,ny+LAP
            do i=1-LAP,nx+LAP
                !Verify each cell's pure 1D VDF views point at the expected backing arrays.
                if(.not.associated(ctr(i,j)%h_vec, cell_h_store(:,i,j)) .or. size(ctr(i,j)%h_vec)/=Nvel) &
                    call AbortCellPointerCheck(i,j,'h_vec')
                if(.not.associated(ctr(i,j)%b_vec, cell_b_store(:,i,j)) .or. size(ctr(i,j)%b_vec)/=Nvel) &
                    call AbortCellPointerCheck(i,j,'b_vec')
                if(.not.associated(ctr(i,j)%sh_vec, cell_sh_store(:,:,i,j)) .or. any(shape(ctr(i,j)%sh_vec)/=[Nvel,2])) &
                    call AbortCellPointerCheck(i,j,'sh_vec')
                if(.not.associated(ctr(i,j)%sb_vec, cell_sb_store(:,:,i,j)) .or. any(shape(ctr(i,j)%sb_vec)/=[Nvel,2])) &
                    call AbortCellPointerCheck(i,j,'sb_vec')
                if(.not.associated(ctr(i,j)%h0_vec, cell_h0_store(:,i,j)) .or. size(ctr(i,j)%h0_vec)/=Nvel) &
                    call AbortCellPointerCheck(i,j,'h0_vec')
                if(.not.associated(ctr(i,j)%b0_vec, cell_b0_store(:,i,j)) .or. size(ctr(i,j)%b0_vec)/=Nvel) &
                    call AbortCellPointerCheck(i,j,'b0_vec')
                if(.not.associated(ctr(i,j)%res_h0_vec, cell_res_h0_store(:,i,j)) .or. size(ctr(i,j)%res_h0_vec)/=Nvel) &
                    call AbortCellPointerCheck(i,j,'res_h0_vec')
                if(.not.associated(ctr(i,j)%res_b0_vec, cell_res_b0_store(:,i,j)) .or. size(ctr(i,j)%res_b0_vec)/=Nvel) &
                    call AbortCellPointerCheck(i,j,'res_b0_vec')
            end do
        end do
    end subroutine CheckCellMicroPointers


    subroutine AbortCellPointerCheck(i,j,name)
        integer(KINT), intent(in)                       :: i,j
        character(len=*), intent(in)                    :: name

        print '(A,I6,2I8,1X,A)', 'Cell micro pointer mismatch (rank,i,j,name):', myID,i,j,trim(name)
        call MPI_Abort(MPI_COMM_WORLD, 31, ierr)
    end subroutine AbortCellPointerCheck


    !--------------------------------------------------
    !>Verify face micro flux pointers are associated with backing arrays.
    !--------------------------------------------------
    subroutine CheckFaceFluxPointers()
        integer(KINT)                                   :: i,j

        do j=1,ny
            do i=1,nx+1
                !Verify vertical-face pure 1D flux views point at the expected backing arrays.
                if(.not.associated(vface(i,j)%L_flux_h_vec, v_l_flux_h_store(:,i,j)) .or. size(vface(i,j)%L_flux_h_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'vLh_vec')
                if(.not.associated(vface(i,j)%L_flux_b_vec, v_l_flux_b_store(:,i,j)) .or. size(vface(i,j)%L_flux_b_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'vLb_vec')
                if(.not.associated(vface(i,j)%R_flux_h_vec, v_r_flux_h_store(:,i,j)) .or. size(vface(i,j)%R_flux_h_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'vRh_vec')
                if(.not.associated(vface(i,j)%R_flux_b_vec, v_r_flux_b_store(:,i,j)) .or. size(vface(i,j)%R_flux_b_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'vRb_vec')
            end do
        end do

        do j=1,ny+1
            do i=1,nx
                !Verify horizontal-face pure 1D flux views point at the expected backing arrays.
                if(.not.associated(hface(i,j)%L_flux_h_vec, h_l_flux_h_store(:,i,j)) .or. size(hface(i,j)%L_flux_h_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'hLh_vec')
                if(.not.associated(hface(i,j)%L_flux_b_vec, h_l_flux_b_store(:,i,j)) .or. size(hface(i,j)%L_flux_b_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'hLb_vec')
                if(.not.associated(hface(i,j)%R_flux_h_vec, h_r_flux_h_store(:,i,j)) .or. size(hface(i,j)%R_flux_h_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'hRh_vec')
                if(.not.associated(hface(i,j)%R_flux_b_vec, h_r_flux_b_store(:,i,j)) .or. size(hface(i,j)%R_flux_b_vec)/=Nvel) &
                    call AbortFacePointerCheck(i,j,'hRb_vec')
            end do
        end do
    end subroutine CheckFaceFluxPointers


    subroutine AbortFacePointerCheck(i,j,name)
        integer(KINT), intent(in)                       :: i,j
        character(len=*), intent(in)                    :: name

        print '(A,I6,2I8,1X,A)', 'Face flux pointer mismatch (rank,i,j,name):', myID,i,j,trim(name)
        call MPI_Abort(MPI_COMM_WORLD, 32, ierr)
    end subroutine AbortFacePointerCheck


    !--------------------------------------------------
    !>Initialize mesh
    !--------------------------------------------------
    subroutine InitMesh()
        logical                 :: ext1
        integer(KINT)           :: nx0,ny0,i,j,k,pi,pj,np_send,i1,j1
        real(KREAL)             :: dx,dy,len1,len2
        real(KREAL),allocatable :: buff1(:,:,:),buff2(:,:,:)
        
        !---------------------------------------------------------
        !>Read main mesh and send to other processors (myID=0)
        !---------------------------------------------------------
        if(myID==0) then
            !check whether the grid file exists
            inquire(file="Mesh2d.dat",exist=ext1)
            if(.not.ext1.and.S2CartesianMesh.and.IsS2CN2())then
                do j=1,Ny_all+1;do i=1,Nx_all+1
                    gridx_all(i,j)=real(i-1,KREAL)*(2.0_KREAL*acos(-1.0_KREAL)/real(Nx_all,KREAL))
                    if(S2W1_wall_mode==1)then
                        gridy_all(i,j)=real(j-1,KREAL)/real(Ny_all,KREAL)
                    else
                        gridy_all(i,j)=real(j-1,KREAL)*(2.0_KREAL*acos(-1.0_KREAL)/real(Ny_all,KREAL))
                    endif
                enddo;enddo
                open(1,file='Mesh2d.dat',form='unformatted',status='new')
                write(1)Nx_all+1,Ny_all+1
                write(1)gridx_all,gridy_all
                close(1);ext1=.true.
            endif
            if(.not. ext1) then           ! the file does not exist
                call MPI_Abort(MPI_COMM_WORLD,473,ierr)
            endif
            
            !read global mesh
            open(1,file='Mesh2d.dat',form='unformatted')
            !read(1)
            read(1) nx0,ny0
            if(nx0/=nx_all+1 .or. ny0/=ny_all+1)call MPI_Abort(MPI_COMM_WORLD,474,ierr)
            read(1) ((gridx_all(i,j),i=1,nx0),j=1,ny0),((gridy_all(i,j),i=1,nx0),j=1,ny0) 
            ! x  y  gridx_all  gridy_all 
            close(1)
            
            !calculate global cell-center coordiantes
            do j=1,ny_all
                do i=1,nx_all
                    gridxc_all(i,j) = 0.25*(gridx_all(i,j)+gridx_all(i+1,j)+gridx_all(i,j+1)+gridx_all(i+1,j+1))
                    gridyc_all(i,j) = 0.25*(gridy_all(i,j)+gridy_all(i+1,j)+gridy_all(i,j+1)+gridy_all(i+1,j+1))
                enddo
            enddo
            
            !send local mesh to each processor
            do pj=0,NP_Y-1 
                do pi=0,NP_X-1
                    np_send=pj*NP_X+pi
                    if(np_send.eq.0) then
                        do j1=1,ny+1                        !total grid point is larger than total cell number +1 
                            do i1=1,nx+1
                                gridx(i1,j1)=gridx_all(i1,j1)
                                gridy(i1,j1)=gridy_all(i1,j1)
                            enddo
                        enddo
                    else
                        allocate(buff1(Block_XY(np_send)%NI+1,Block_XY(np_send)%NJ+1,2))
                        do j1=1,Block_XY(np_send)%NJ+1
                            do i1=1,Block_XY(np_send)%NI+1
                                i=Block_XY(np_send)%Origin_I+i1-1
                                j=Block_XY(np_send)%Origin_J+j1-1
                                buff1(i1,j1,1)=gridx_all(i,j)
                                buff1(i1,j1,2)=gridy_all(i,j)
                            enddo
                        enddo
                        call MPI_Send(buff1,(Block_XY(np_send)%NI+1)*(Block_XY(np_send)%NJ+1)*2,MPI_KREAL,np_send,9000,MPI_COMM_WORLD,ierr) 
                        !cut mesh data into blocks and send to corresponding processor
                        deallocate(buff1)
                    endif
                enddo
            enddo
        else
            allocate(buff2(nx+1,ny+1,2))
            call MPI_Recv(buff2,(nx+1)*(ny+1)*2,MPI_KREAL,0,9000,MPI_COMM_WORLD,Status,ierr)
            !receive local mesh data
            do j1=1,ny+1
                do i1=1,nx+1
                    gridx(i1,j1)=buff2(i1,j1,1)
                    gridy(i1,j1)=buff2(i1,j1,2)
                enddo
            enddo
            deallocate(buff2)
        endif
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        !-------------------------------------------------------------------------------------------
        !>For each block: calculate cell coordinate, length, area, and interface length and direction
        !-------------------------------------------------------------------------------------------
        !Cell: x and y coordinates, length_dx, length_dy, cell area      
        do j=1,ny
            do i=1,nx
                !cell center coordinate
                ctr(i,j)%x  =  0.25*(gridx(i,j)+gridx(i+1,j)+gridx(i,j+1)+gridx(i+1,j+1))
                ctr(i,j)%y  =  0.25*(gridy(i,j)+gridy(i+1,j)+gridy(i,j+1)+gridy(i+1,j+1))
                !length of horizontal direction
                len1        =  sqrt( (gridx(i+1,j)-gridx(i,j))**2+(gridy(i+1,j)-gridy(i,j))**2 )
                len2        =  sqrt( (gridx(i+1,j+1)-gridx(i,j+1))**2+(gridy(i+1,j+1)-gridy(i,j+1))**2 )
                ctr(i,j)%length(1)  =  0.5*(len1+len2)
                !length of vertical direction
                len1        =  sqrt( (gridx(i,j+1)-gridx(i,j))**2+(gridy(i,j+1)-gridy(i,j))**2 )
                len2        =  sqrt( (gridx(i+1,j+1)-gridx(i+1,j))**2+(gridy(i+1,j+1)-gridy(i+1,j))**2 )
                ctr(i,j)%length(2)  =  0.5*(len1+len2)
                !area
                ctr(i,j)%area   =  abs((gridx(i,j)-gridx(i+1,j+1))*(gridy(i+1,j)-gridy(i,j+1)) -  &
                                   (gridx(i+1,j)-gridx(i,j+1))*(gridy(i,j)-gridy(i+1,j+1)) )*0.5
            enddo
        enddo
        
        !Vertical interface: length and normal vector     
        do j=1,ny
            do i=1,nx+1
                dx = gridx(i,j+1)-gridx(i,j)
                dy = gridy(i,j+1)-gridy(i,j)
                vface(i,j)%length = sqrt(dx*dx+dy*dy)       ! length
                vface(i,j)%cosx   = dy/vface(i,j)%length    ! normal vector at (i,j) or (I-1/2,J)
                vface(i,j)%cosy   = -dx/vface(i,j)%length
            enddo
        enddo

        !Horizontal interface: length and normal vector
        do j=1,ny+1
            do i=1,nx
                dx=gridx(i+1,j)-gridx(i,j)
                dy=gridy(i+1,j)-gridy(i,j)
                hface(i,j)%length = sqrt(dx*dx+dy*dy)      ! length 
                hface(i,j)%cosx   = -dy/hface(i,j)%length  ! normal vector at i, j+1/2
                hface(i,j)%cosy   = dx/hface(i,j)%length    
            enddo
        enddo
        
        !-------------------------------------------------------
        !>Set ghost cell value: cell coordinates, length and area
        !-------------------------------------------------------
        !Update overlaped region
        call MPIUpdate_CellCoordinates()
        call MPIUpdate_Conner_Coordinates()
        
        !Physical boundary: interpolation
        if(IP_X==0) then        !left boundary
            do k=1,ghost
                ctr(1-k,1:ny)%x          = 2.0*ctr(1,1:ny)%x - ctr(1+k,1:ny)%x
                ctr(1-k,1:ny)%y          = 2.0*ctr(1,1:ny)%y - ctr(1+k,1:ny)%y
                ctr(1-k,1:ny)%length(1)  = ctr(1+k-1,1:ny)%length(1)
                ctr(1-k,1:ny)%length(2)  = ctr(1+k-1,1:ny)%length(2)
                ctr(1-k,1:ny)%area       = ctr(1+k-1,1:ny)%area
            enddo
        endif
        
        if(IP_X==NP_X-1) then   !right boundary
            do k=1,ghost
                ctr(nx+k,1:ny)%x          = 2.0*ctr(nx,1:ny)%x - ctr(nx-k,1:ny)%x
                ctr(nx+k,1:ny)%y          = 2.0*ctr(nx,1:ny)%y - ctr(nx-k,1:ny)%y
                ctr(nx+k,1:ny)%length(1)  = ctr(nx-k+1,1:ny)%length(1)
                ctr(nx+k,1:ny)%length(2)  = ctr(nx-k+1,1:ny)%length(2)
                ctr(nx+k,1:ny)%area       = ctr(nx-k+1,1:ny)%area
            enddo
        endif
        
        if(IP_Y==0) then        !bottom boundary
            do k=1,ghost
                ctr(1:nx,1-k)%x          = 2.0*ctr(1:nx,1)%x - ctr(1:nx,1+k)%x
                ctr(1:nx,1-k)%y          = 2.0*ctr(1:nx,1)%y - ctr(1:nx,1+k)%y
                ctr(1:nx,1-k)%length(1)  = ctr(1:nx,1+k-1)%length(1)
                ctr(1:nx,1-k)%length(2)  = ctr(1:nx,1+k-1)%length(2)
                ctr(1:nx,1-k)%area       = ctr(1:nx,1+k-1)%area 
            enddo
        endif
        
        if(IP_Y==NP_Y-1) then   !upper boundary
            do k=1,ghost
                ctr(1:nx,ny+k)%x          = 2.0*ctr(1:nx,ny)%x - ctr(1:nx,ny-k)%x
                ctr(1:nx,ny+k)%y          = 2.0*ctr(1:nx,ny)%y - ctr(1:nx,ny-k)%y
                ctr(1:nx,ny+k)%length(1)  = ctr(1:nx,ny-k+1)%length(1)
                ctr(1:nx,ny+k)%length(2)  = ctr(1:nx,ny-k+1)%length(2)
                ctr(1:nx,ny+k)%area       = ctr(1:nx,ny-k+1)%area
            enddo
        endif
              
    end subroutine InitMesh

    !--------------------------------------------------
    !>Initialize velocity
    !--------------------------------------------------
    subroutine InitVelocity()
        real(KREAL) :: raw_node(Nu),raw_weight(Nu)

        if(len_trim(StructuredVelocityFile)>0)then
            call InitVelocityFromStructuredFile()
        elseif (IsS2CN2().or.IsS2Steady().or.trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau') then
            ! S2CN2 uses the exact effective GH16 rule frozen by S2D-1-R1.
            ! uSpace/vSpace/weight remain the native authoritative velocity arrays.
            call S2InitVelocity(uSpace,vSpace,weight,raw_node,raw_weight)
        elseif (Quadrature_method==Quad_NewtonCotes4) then          !Initialize discrete velocity space using Newton-Cotes formulas
            call InitVelocityNewtonCotes_4th(Nu,Nv)         
        elseif (Quadrature_method==Quad_GaussHermite) then      !Initialize discrete velocity space using Gaussian-Hermite type quadrature
            call InitVelocityGaussHermite(Nu,Nv)
        elseif (Quadrature_method==Quad_Unstruct) then          !Initialize discrete velocity space using unstructured grid
            call InitVelocity_Unstruct(Nu, Nv)  
        else
            stop "Error: the Quadrature_method is illegal."
        end if
        call FlattenVelocitySpace()
        if (IsImplicitSteady() .or. IsImplicitUnsteady()) call BuildStructuredMirrorMaps()
        call CheckFlattenedVelocitySpace()
    end subroutine InitVelocity


    !--------------------------------------------------
    !>Synchronize primary Nvel velocity arrays from the velocity-space cache.
    !--------------------------------------------------
    subroutine InitVelocityFromStructuredFile()
        integer :: iu,ios,n1,n2,k,i,j,l
        logical :: valid
        if (myID==0) then
            open(newunit=iu,file=trim(StructuredVelocityFile),status='old',action='read',iostat=ios)
            if (ios/=0) call MPI_Abort(MPI_COMM_WORLD,201,ierr)
            read(iu,*,iostat=ios) n1,n2
            if (ios/=0) call MPI_Abort(MPI_COMM_WORLD,202,ierr)
            if(n1/=Nu.or.n2/=Nv) call MPI_Abort(MPI_COMM_WORLD,203,ierr)
            do k=1,Nvel
                read(iu,*,iostat=ios) uSpace(k),vSpace(k),weight(k)
                if (ios/=0) call MPI_Abort(MPI_COMM_WORLD,204,ierr)
            enddo
            close(iu)
        endif
        call MPI_Bcast(uSpace,Nvel,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(vSpace,Nvel,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(weight,Nvel,MPI_KREAL,0,MPI_COMM_WORLD,ierr)
        valid=all(ieee_is_finite(uSpace)).and.all(ieee_is_finite(vSpace)) &
            .and.all(ieee_is_finite(weight)).and.all(weight>0.0_KREAL)
        do j=1,Nv; do i=1,Nu
            k=i+(j-1)*Nu
            valid=valid.and.uSpace(k)==uSpace(i).and.vSpace(k)==vSpace(1+(j-1)*Nu)
            l=(Nu+1-i)+(Nv-j)*Nu
            valid=valid.and.uSpace(k)==-uSpace(l).and.vSpace(k)==-vSpace(l).and.weight(k)==weight(l)
        enddo; enddo
        if (.not.valid) call MPI_Abort(MPI_COMM_WORLD,205,ierr)
        U_MIN=minval(uSpace); U_MAX=maxval(uSpace)
        V_MIN=minval(vSpace); V_MAX=maxval(vSpace)
    end subroutine InitVelocityFromStructuredFile


    subroutine FlattenVelocitySpace()
        integer(KINT)                                   :: k

        Nvel = Nu*Nv
        if(.not.allocated(uVel) .or. size(uVel)/=Nvel) then
            if(allocated(uVel)) deallocate(uVel)
            if(allocated(vVel)) deallocate(vVel)
            if(allocated(wVel)) deallocate(wVel)
            allocate(uVel(Nvel))
            allocate(vVel(Nvel))
            allocate(wVel(Nvel))
        end if

        !Keep uVel/vVel/wVel as the primary names used by compute kernels.
        do k=1,Nvel
            uVel(k) = uSpace(k)
            vVel(k) = vSpace(k)
            wVel(k) = weight(k)
        end do
    end subroutine FlattenVelocitySpace


    !--------------------------------------------------
    !>Verify compute velocity arrays mirror the pure 1D velocity-space cache.
    !--------------------------------------------------
    subroutine CheckFlattenedVelocitySpace()
        real(KREAL)                                     :: max_error
        integer(KINT)                                   :: k

        max_error = 0.0_KREAL
        do k=1,Nvel
            max_error = max(max_error, abs(uVel(k)-uSpace(k)))
            max_error = max(max_error, abs(vVel(k)-vSpace(k)))
            max_error = max(max_error, abs(wVel(k)-weight(k)))
        end do

        if(max_error>0.0_KREAL) then
            print '(A,I6,1P,E15.7)', 'Flattened velocity-space mismatch on rank:', myID, max_error
            call MPI_Abort(MPI_COMM_WORLD, 30, ierr)
        end if
    end subroutine CheckFlattenedVelocitySpace

    !--------------------------------------------------
    !>Initialize discrete velocity space using Newton-Cotes formulas
    !--------------------------------------------------
    subroutine InitVelocityNewtonCotes_4th(num_u,num_v)
        integer(KINT), intent(inout)                    :: num_u,num_v
        real(KREAL)                                     :: du,dv !Spacing in u and v velocity space
        integer(KINT)                                   :: i,j,k

        !spacing in u and v velocity space
        du = (U_MAX-U_MIN)/(num_u-1)
        dv = (V_MAX-V_MIN)/(num_v-1)

        !Fill pure 1D velocity-space arrays while preserving the old Nu,Nv traversal.
        do j=1,num_v
            do i=1,num_u
                k = i+(j-1)*num_u
                uSpace(k) = U_MIN+(i-1)*du
                vSpace(k) = V_MIN+(j-1)*dv
                weight(k) = (newtonCoeff(i,num_u)*du)*(newtonCoeff(j,num_v)*dv)
            enddo
        enddo

        contains
            !--------------------------------------------------
            !>Calculate the coefficient for newton-cotes formula
            !>@param[in] idx          :index in velocity space
            !>@param[in] num          :total number in velocity space
            !>@return    newtonCoeff  :coefficient for newton-cotes formula
            !--------------------------------------------------
            pure function newtonCoeff(idx,num)
                integer(KINT), intent(in)               :: idx,num
                real(KREAL)                             :: newtonCoeff

                if (idx==1 .or. idx==num) then 
                    newtonCoeff = 14.0/45.0
                else if (mod(idx-5,4)==0) then
                    newtonCoeff = 28.0/45.0
                else if (mod(idx-3,4)==0) then
                    newtonCoeff = 24.0/45.0
                else
                    newtonCoeff = 64.0/45.0
                end if
            end function newtonCoeff
    end subroutine InitVelocityNewtonCotes_4th
        
    !--------------------------------------------------
    !>Set discrete velocity space using Gaussian-Hermite type quadrature
    !>@param[in] umid,vmid :middle value of the velocity space, zero or macroscopic velocity
    !--------------------------------------------------
    subroutine InitVelocityGaussHermite(num_u,num_v)
        integer(KINT), intent(inout)                    :: num_u,num_v
        real(KREAL)                                     :: umid,vmid
        real(KREAL),allocatable                         :: vcoords(:), weights(:) !Velocity points and weight for 28 points (symmetry)
        integer(KINT)                                   :: i,j,k

        allocate(vcoords(num_u), weights(num_u))
        if(num_u==28) then   
            !Set 28*28 velocity points and weight
            vcoords = [ -0.5392407922630E+01, -0.4628038787602E+01, -0.3997895360339E+01, -0.3438309154336E+01,&
                        -0.2926155234545E+01, -0.2450765117455E+01, -0.2007226518418E+01, -0.1594180474269E+01,&
                        -0.1213086106429E+01, -0.8681075880846E+00, -0.5662379126244E+00, -0.3172834649517E+00,&
                        -0.1331473976273E+00, -0.2574593750171E-01, +0.2574593750171E-01, +0.1331473976273E+00,&
                        +0.3172834649517E+00, +0.5662379126244E+00, +0.8681075880846E+00, +0.1213086106429E+01,&
                        +0.1594180474269E+01, +0.2007226518418E+01, +0.2450765117455E+01, +0.2926155234545E+01,&
                        +0.3438309154336E+01, +0.3997895360339E+01, +0.4628038787602E+01, +0.5392407922630E+01 ]

            weights = [ +0.2070921821819E-12, +0.3391774320172E-09, +0.6744233894962E-07, +0.3916031412192E-05,&
                        +0.9416408715712E-04, +0.1130613659204E-02, +0.7620883072174E-02, +0.3130804321888E-01,&
                        +0.8355201801999E-01, +0.1528864568113E+00, +0.2012086859914E+00, +0.1976903952423E+00,&
                        +0.1450007948865E+00, +0.6573088665062E-01, +0.6573088665062E-01, +0.1450007948865E+00,&
                        +0.1976903952423E+00, +0.2012086859914E+00, +0.1528864568113E+00, +0.8355201801999E-01,&
                        +0.3130804321888E-01, +0.7620883072174E-02, +0.1130613659204E-02, +0.9416408715712E-04,&
                        +0.3916031412192E-05, +0.6744233894962E-07, +0.3391774320172E-09, +0.2070921821819E-12 ]
        elseif(num_u==16) then
            !Set 16*16 velocity points and weight
            vcoords = [ -0.3686007162724397E+1, -0.2863133883708075E+1, -0.2183921153095858E+1, -0.1588855862270055E+1,&
                        -0.1064246312116224E+1, -0.6163028841823999, -0.2673983721677653, -0.5297864393185113E-1,&
                        0.5297864393185113E-1, 0.2673983721677653, 0.6163028841823999, 0.1064246312116224E+1,&
                        0.1588855862270055E+1, 0.2183921153095858E+1, 0.2863133883708075E+1, 0.3686007162724397E+1]

            weights = [ 0.1192596926595344E-5, 0.2020636491324107E-3, 0.5367935756025333E-2, 0.4481410991746290E-1,&
                        0.1574482826187903, 0.2759533979884218, 0.2683307544726388, 0.1341091884533595,&
                        0.1341091884533595, 0.2683307544726388, 0.2759533979884218, 0.1574482826187903,&
                        0.4481410991746290E-1, 0.5367935756025333E-2, 0.2020636491324107E-3, 0.1192596926595344E-5]
        elseif(num_u==8) then
            !Set 8*8 velocity points and weight
            vcoords = [ -0.2262664477010362E+1, -0.1342537825644992E+1, -0.6243246901871900E0, -0.1337764469960676E0,&
                        0.1337764469960676E0, 0.6243246901871900E0, 0.1342537825644992E+1, 0.2262664477010362E+1]
            
            weights = [ 0.6374323486257276E-2, 0.1334425003575195E0, 0.4211071018520622E0, 0.3253029997569190E0,&
                        0.3253029997569190E0, 0.4211071018520622E0, 0.1334425003575195E0, 0.6374323486257276E-2]
        endif

        umid = 0.0
        vmid = 0.0
        !Fill pure 1D velocity-space arrays while preserving the old Nu,Nv traversal.
        do j=1,num_v
            do i=1,num_u
                k = i+(j-1)*num_u
                uSpace(k) = umid+vcoords(i)
                vSpace(k) = vmid+vcoords(j)
                weight(k) = (weights(i)*exp(vcoords(i)**2))*(weights(j)*exp(vcoords(j)**2))
            enddo
        enddo
        
        deallocate(vcoords, weights)
    end subroutine InitVelocityGaussHermite


    !--------------------------------------------------
    !>Initialize velocity space using unstructured data
    !>Optimized for Symmetric Mesh: Adjacent Pairing (k, k+1)
    !--------------------------------------------------
    subroutine InitVelocity_Unstruct(num_u, num_v)
        integer(KINT), intent(in) :: num_u, num_v
        integer(KINT)             :: k, N_check
        character(len=100)        :: buffer
        logical                   :: file_exists

        if (myID == 0) then
            inquire(file="microProperties2d.dat", exist=file_exists)
            if (.not. file_exists) stop "Error: microProperties2d.dat not found!"

            open(20, file='microProperties2d.dat', status='old', action='read')
            read(20, '(A)') buffer
            read(20, *) N_check
            
            if (N_check /= num_u) then
                print*, "Error: Nv mismatch! File:", N_check, " Code:", num_u
                stop
            end if
            
            ! check Nv is even
            if (mod(num_u, 2) /= 0) stop "Error: Nv is odd, expected even pairs."

            ! read coordinates of velocity space
            read(20, '(A)') buffer
            do k = 1, num_u
                read(20, *) uSpace(k), vSpace(k)
            end do

            ! read weights
            read(20, '(A)') buffer
            do k = 1, num_u
                read(20, *) weight(k)
            end do
            close(20)
            
            !===================================================================
            ! calculate new U_MAX, U_MIN, V_MAX, V_MIN based on loaded grid
            !===================================================================

            
            U_MAX = maxval(uSpace(1:num_u))
            U_MIN = minval(uSpace(1:num_u))
            V_MAX = maxval(vSpace(1:num_u))
            V_MIN = minval(vSpace(1:num_u))
            
            ! Print confirmation on screen for verification
            write(*,*) "--------------------------------------------------------"
            write(*,*) " [INFO] Unstructured Velocity Grid Loaded."
            write(*,*) " [INFO] Overwriting Config.ini velocity bounds:"
            write(*,*) "        New U_Range: [", U_MIN, ", ", U_MAX, "]"
            write(*,*) "        New V_Range: [", V_MIN, ", ", V_MAX, "]"
            write(*,*) "--------------------------------------------------------"

            ! --- 4. Establish symmetric mapping (core optimization) ---
            ! Since the data in the file is already paired as (v>0, v<0)
            ! We directly bind k with k+1
            allocate(SymMap_V(num_u))
            do k = 1, num_u, 2
                SymMap_V(k)     = k + 1  ! positive v points to negative v
                SymMap_V(k+1)   = k      ! negative v points to positive v
            end do
            
            print*, "Unstruct Grid: Symmetric pairs mapped successfully."
        end if

        ! --- 5. Broadcast data ---
        call MPI_BCAST(uSpace, num_u*num_v, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(vSpace, num_u*num_v, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(weight, num_u*num_v, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)
        
        if (myID /= 0) allocate(SymMap_V(num_u))
        call MPI_BCAST(SymMap_V, num_u, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
        !===================================================================
        ! bcast updated scalar U_MAX to all other processes
        !===================================================================
        ! only myID=0 process reads unstructured grid from file and updates U_MAX, U_MIN, V_MAX, V_MIN,need to broadcast to others.

        
        call MPI_BCAST(U_MAX, 1, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(U_MIN, 1, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(V_MAX, 1, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(V_MIN, 1, MPI_KREAL, 0, MPI_COMM_WORLD, ierr)

    end subroutine InitVelocity_Unstruct

    
    !--------------------------------------------------
    !>Set the initial flow
    !--------------------------------------------------
    subroutine InitFlowField()
        !integer(KINT)   :: i,j
        
        if(IsS2CN2())then
            select case(Iflag_InitFlow)
            case(Init_Freestream)
                call InitFlow_Freestream
            case(Init_Macroflow)
                call InitFlow_Macroflow
            case(Init_Restart)
                if(.not.UnifiedPhysicalOperators)call S2NativeReadMetadata(.true.)
                call InitFlow_Restart
            case default
                error stop 'unsupported native S2 Iflag_InitFlow'
            end select
            call S2NativeInitializeState(.true.,Iflag_InitFlow/=Init_Restart)
            if(Iflag_InitFlow==Init_Restart.and..not.UnifiedPhysicalOperators)call S2NativeReadMetadata(.false.)
        elseif(IsS2Steady().or.trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau') then
            call S2NativeInitializeState()
        elseif(Iflag_DIS==0) then
            select case(Iflag_InitFlow)
            case(Init_Freestream)
                call InitFlow_Freestream
            case(Init_Macroflow)
                call InitFlow_Macroflow
            case(Init_Restart)
                call InitFlow_Restart
            end select
        elseif(Iflag_DIS==1) then
            call InitFlow_Restart
        endif
        
        !Update the overlaped region using MPI_Sendrecv
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPIUpdate_CellConVars()
        call MPIUpdate_CellVDF()
        call MPIUpdate_Conner_conVars()
        call MPIUpdate_Conner_VDF()
        call SnapshotInitialState()
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
    end subroutine InitFlowField


    !--------------------------------------------------
    !>Store the initial state for both steady and DIS runs
    !--------------------------------------------------
    subroutine SnapshotInitialState()
        integer(KINT) :: i,j,k !Cell loop indices and flattened velocity-space index

        do j=1-LAP,ny+LAP
            do i=1-LAP,nx+LAP
                ctr(i,j)%conVars0(1:4)   = ctr(i,j)%conVars(1:4)
                !Snapshot micro state and clear micro residuals through flattened velocity views.
                do k=1,Nvel
                    ctr(i,j)%h0_vec(k) = ctr(i,j)%h_vec(k)
                    ctr(i,j)%b0_vec(k) = ctr(i,j)%b_vec(k)
                    ctr(i,j)%res_h0_vec(k) = 0._KREAL
                    ctr(i,j)%res_b0_vec(k) = 0._KREAL
                enddo
                ctr(i,j)%res_conVars0    = 0._KREAL
            end do
        end do
    end subroutine SnapshotInitialState
    
    
    !--------------------------------------------------
    !>Initialize flow using freestream or initial parameter
    !--------------------------------------------------
    subroutine InitFlow_Freestream
        real(KREAL), allocatable, dimension(:)          :: H_vec,B_vec !Pure 1D freestream distributions
        real(KREAL)                                     :: conVars(4) !Conservative variables
        integer(KINT)                                   :: i,j,k !Cell loop indices and flattened velocity-space index

        !Allocate pure 1D Maxwellian arrays for freestream initialization.
        allocate(H_vec(Nvel))
        allocate(B_vec(Nvel))

        !Get conservative variables and Maxwellian distribution function
        conVars = GetConserved(Prim_inf)
        call DiscreteMaxwell_Vec(H_vec,B_vec,Prim_inf)

        !Initialize field
        do j=1,ny
            do i=1,nx
                ctr(i,j)%dt_local=dt_global
                ctr(i,j)%conVars = conVars
                !Copy freestream VDFs through flattened cell-center views.
                do k=1,Nvel
                    ctr(i,j)%h_vec(k) = H_vec(k)
                    ctr(i,j)%b_vec(k) = B_vec(k)
                enddo
            enddo
        enddo

        !Deallocation
        deallocate(H_vec)
        deallocate(B_vec)
    
    end subroutine InitFlow_Freestream
    
    
    !--------------------------------------------------
    !>Initialize flow using macroflow field
    !--------------------------------------------------
    subroutine InitFlow_Macroflow
        real(KREAL), allocatable, dimension(:,:,:)      :: flow_all
        real(KREAL), allocatable, dimension(:,:,:)      :: buff1,buff2
        real(KREAL), allocatable, dimension(:)          :: H_vec,B_vec !Pure 1D local Maxwellian arrays
        logical                                         :: ext1
        character(len=20)                               :: char
        real(KREAL)                                     :: conVars(4),prim(4),tmp
        integer(KINT)                                   :: nx0,ny0,i,j,k,pi,pj,np_send,i1,j1,Ndata
        
        !---------------------------------------------------------
        !>Init macroflow: read flow field and send to other processors
        !---------------------------------------------------------
        Ndata = 4       !ctr(:,:)%conVars(1:4)
        
        if(myID==0) then
            !check whether the macroflow file exists
            inquire(file="Flow2d.plt",exist=ext1)
            if(.not. ext1) then           ! the file does not exist
                call MPI_Abort(MPI_COMM_WORLD,475,ierr)
            endif
        
            !read main flow
            allocate(flow_all(nx_all,ny_all,Ndata))
            !---------------unformatted data-------------
            !open(1,file='Flow2d.dat',form='unformatted')
            !read(1)
            !read(1) nx0,ny0
            !if(nx0/=nx_all .or. ny0/=ny_all) stop 'Error: wrong cell numbers'
            !read(1) (((flow_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,Ndata)             !flow_all: density, u, v, T
            !close(1)
            !---------------formatted data-------------
            open(1,file='Flow2d.plt')
            read(1,*)
            read(1,'(a20,I8,a20,I8)') char,nx0,char,ny0
            if(nx0/=nx_all .or. ny0/=ny_all) stop 'Error: wrong cell numbers'
            read(1,*) ((tmp,i=1,nx_all),j=1,ny_all),((tmp,i=1,nx_all),j=1,ny_all),(((flow_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,Ndata)  !flow_all: density, u, v, T
            close(1)
            
            !transform to conserved variables
            do j=1,ny_all
                do i=1,nx_all
                    flow_all(i,j,4) = 1.0/flow_all(i,j,4)                       !flow_all: density, u, v, 1/T                             
                    conVars(1:Ndata)= GetConserved(flow_all(i,j,1:Ndata))  
                    flow_all(i,j,:) = conVars(:)                                !flow_all: density, x-momentum, y-momentum, total energy                    
                enddo
            enddo
            
            !send to each processor
            do pj=0,NP_Y-1 
                do pi=0,NP_X-1
                    np_send=pj*NP_X+pi
                    if(np_send.eq.0) then
                        do j1=1,ny
                            do i1=1,nx
                                ctr(i1,j1)%conVars(1:Ndata) = flow_all(i1,j1,1:Ndata)
                            enddo
                        enddo
                    else
                        allocate(buff1(Block_XY(np_send)%NI,Block_XY(np_send)%NJ,Ndata))
                        do j1=1,Block_XY(np_send)%NJ
                            do i1=1,Block_XY(np_send)%NI
                                i=Block_XY(np_send)%Origin_I+i1-1
                                j=Block_XY(np_send)%Origin_J+j1-1
                                buff1(i1,j1,1:Ndata)=flow_all(i,j,1:Ndata)
                            enddo
                        enddo
                        call MPI_Send(buff1,Block_XY(np_send)%NI*Block_XY(np_send)%NJ*Ndata,MPI_KREAL,np_send,9000,MPI_COMM_WORLD,ierr)
                        deallocate(buff1)
                    endif
                enddo
            enddo
            deallocate(flow_all)
        else
            allocate(buff2(nx,ny,Ndata))
            call MPI_Recv(buff2,nx*ny*Ndata,MPI_KREAL,0,9000,MPI_COMM_WORLD,Status,ierr)
            do j1=1,ny
                do i1=1,nx
                    ctr(i1,j1)%conVars(1:Ndata) = buff2(i1,j1,1:Ndata)
                enddo
            enddo
            deallocate(buff2)
        endif
        
        !---------------------------------------------------------
        !>Init micro-VDF: assume local Maxwell distribution at each cell
        !---------------------------------------------------------
        !Allocate pure 1D Maxwellian arrays for macroflow initialization.
        allocate(H_vec(Nvel))
        allocate(B_vec(Nvel))

        !Initialize field
        do j=1,ny
            do i=1,nx
                prim = GetPrimary(ctr(i,j)%conVars)
                call DiscreteMaxwell_Vec(H_vec,B_vec,Prim)
                !Copy local Maxwellian VDFs through flattened cell-center views.
                do k=1,Nvel
                    ctr(i,j)%h_vec(k) = H_vec(k)
                    ctr(i,j)%b_vec(k) = B_vec(k)
                enddo
            enddo
        enddo

        !Deallocation
        deallocate(H_vec)
        deallocate(B_vec)
    
    end subroutine
    

    !---------------------------------------------------------
    !>Initialize flow from restart file: flow2d.dat, vdf2d.dat
    !---------------------------------------------------------
    subroutine InitFlow_Restart
        integer(KINT)                                   :: nx0,ny0,nvar0,i,j,k,l,m,pi,pj,np_send,i1,j1,Ndata,tag
        integer(KINT)                                   :: flat_k !Flattened velocity-space index for restart VDF copies
        real(KREAL), allocatable, dimension(:,:,:)      :: solution_all
        real(KREAL), allocatable, dimension(:,:,:)      :: buff1,buff2
        
        ! --- structured mesh VDF variables: keep file order, store velocity dimension as flat Nvel ---
        real(KREAL), allocatable, dimension(:,:,:,:)    :: solution_vdf_all !Restart VDF cache with flat Nvel velocity dimension
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff1_vdf,buff2_vdf !MPI buffers with flat Nvel velocity dimension
        
        ! --- unstructured mesh VDF variables ---
        integer(KINT)                                   :: dummy_nv
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff_vdf_all_unstruct     !Global unstructured VDF cache stored as (Nx_all, Ny_all, Nvel, 2)
        real(KREAL), allocatable, dimension(:,:,:,:)    :: buff_vdf_local_unstruct   !Local unstructured VDF cache stored as (NI, NJ, Nvel, 2)
        real(KREAL), allocatable, dimension(:)          :: temp_coords               ! used to skip coordinate reading
        
        logical                                         :: ext1
        !character(len=20)                               :: char
        
        !---------------------------------------------------------
        !>Init macroflow: read flow2d.dat and send to other processors 
        ! ==> ctr(:,:)%conVars(1:4)
        !---------------------------------------------------------
        Ndata = 4
        
        if(myID==0) then
            !--read flow2d.dat
            allocate(solution_all(nx_all,ny_all,Ndata))
            inquire(file="flow2d.dat",exist=ext1)       !check whether the file exists
            if(.not. ext1) then
                stop 'Restart error: file flow2d.dat does not exist.'
            endif
            open(1,file='flow2d.dat',form='unformatted')
            read(1)
            read(1) nx0,ny0,nvar0
            
            if(nx0/=nx_all .or. ny0/=ny_all .or. nvar0/=Ndata) then
                print *, 'nx0 =', nx0, ', ny0 =', ny0, ', nvar0 =', nvar0
                stop 'Restart error: data in flow2d.dat does not conform.'
            endif
            read(1) (((solution_all(i,j,k),i=1,nx_all),j=1,ny_all),k=1,Ndata)
            close(1)
            
            !--send to each processor
            do pj=0,NP_Y-1 
                do pi=0,NP_X-1
                    np_send=pj*NP_X+pi
                    if(np_send.eq.0) then
                        do j1=1,ny
                            do i1=1,nx
                                ctr(i1,j1)%conVars(1:Ndata) = solution_all(i1,j1,1:Ndata)
                            enddo
                        enddo
                    else
                        allocate(buff1(Block_XY(np_send)%NI,Block_XY(np_send)%NJ,Ndata))
                        do j1=1,Block_XY(np_send)%NJ
                            do i1=1,Block_XY(np_send)%NI
                                i=Block_XY(np_send)%Origin_I+i1-1
                                j=Block_XY(np_send)%Origin_J+j1-1
                                buff1(i1,j1,1:Ndata)=solution_all(i,j,1:Ndata)
                            enddo
                        enddo
                        tag=np_send     !tag equals to block index
                        call MPI_Send(buff1,Block_XY(np_send)%NI*Block_XY(np_send)%NJ*Ndata,MPI_KREAL,np_send,tag,MPI_COMM_WORLD,ierr)
                        deallocate(buff1)
                    endif
                enddo
            enddo
            deallocate(solution_all)
        else
            allocate(buff2(nx,ny,Ndata))
            tag=IP_Y*NP_X+IP_X      !tag equals to block index
            call MPI_Recv(buff2,nx*ny*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,Status,ierr)
            do j1=1,ny
                do i1=1,nx
                    ctr(i1,j1)%conVars(1:Ndata) = buff2(i1,j1,1:Ndata)
                enddo
            enddo
            deallocate(buff2)
        endif
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        !---------------------------------------------------------
        !>Init micro-vdf: read vdf2d.dat (or unstruct) and send to other processors.
        !  VDF values are stored through h_vec/b_vec; legacy h/b 2D views share the same backing storage.
        !---------------------------------------------------------
        Ndata = 2   ! reset to 2 (h, b)
        
        if(myID==0) then
            !============================================================
            !  Master Process: Read and Distribute
            !============================================================
            
            !>>> Branch A: Unstructured Grid <<<
            if (Quadrature_method == Quad_Unstruct) then
                
                ! 1. Check if file exists
                inquire(file="vdf2d_unstruct_full.dat", exist=ext1)
                if(.not. ext1) stop 'Restart error: file vdf2d_unstruct_full.dat does not exist.'

                open(1, file='vdf2d_unstruct_full.dat', form='unformatted')

                ! 2. Read file header
                ! Output write logic: write(1) Nx_all, Ny_all,Nvel,1
                read(1) nx0, ny0, nvar0, dummy_nv

                ! 3. Dimension check (nvar0 stores flattened velocity count for unstructured files)
                if(nx0/=nx_all .or. ny0/=ny_all .or. nvar0/=Nvel) then
                    print*, "Restart Mismatch: Physical Grid Mismatch! File(", nx0, ny0, ") vs Code(", nx_all, ny_all, ")"
                    stop 'Restart error: Dimension mismatch in vdf2d_unstruct_full.dat'
                endif

                ! 4. Skip coordinate information (uSpace, vSpace, weight)
                allocate(temp_coords(Nvel))
                read(1) temp_coords, temp_coords, temp_coords  ! Skip uSpace, vSpace, weight
                deallocate(temp_coords)

                ! 5. Read full field data (h and b)
                ! Memory order: (i, j, flat velocity, h/b) -> Fortran array (Nx, Ny, Nvel, 2)
                allocate(buff_vdf_all_unstruct(nx_all, ny_all, Nvel, Ndata))
                read(1) buff_vdf_all_unstruct
                close(1)

                ! 6. Scatter data (Scatter)
                do pj=0,NP_Y-1 
                    do pi=0,NP_X-1
                        np_send = pj*NP_X + pi
                        
                        if(np_send == 0) then
                            ! Copy to self (Root)
                            do j1=1,ny
                                do i1=1,nx
                                    !Copy unstructured restart VDFs through flattened cell-center views.
                                    do k=1,Nvel
                                        ctr(i1,j1)%h_vec(k) = buff_vdf_all_unstruct(i1, j1, k, 1)
                                        ctr(i1,j1)%b_vec(k) = buff_vdf_all_unstruct(i1, j1, k, 2)
                                    enddo
                                enddo
                            enddo
                        else
                            ! Send to child process
                            allocate(buff_vdf_local_unstruct(Block_XY(np_send)%NI, Block_XY(np_send)%NJ, Nvel, Ndata))
                            
                            do j1=1,Block_XY(np_send)%NJ
                                do i1=1,Block_XY(np_send)%NI
                                    i = Block_XY(np_send)%Origin_I + i1 - 1
                                    j = Block_XY(np_send)%Origin_J + j1 - 1
                                    buff_vdf_local_unstruct(i1,j1,:,:) = buff_vdf_all_unstruct(i,j,:,:)
                                enddo
                            enddo

                            ! Use Tag = 20 to distinguish unstructured data
                            tag = np_send * 20
                            call MPI_Send(buff_vdf_local_unstruct, &
                                          Block_XY(np_send)%NI * Block_XY(np_send)%NJ * Nvel * Ndata, &
                                          MPI_KREAL, np_send, tag, MPI_COMM_WORLD, ierr)
                            
                            deallocate(buff_vdf_local_unstruct)
                        endif
                    enddo
                enddo
                deallocate(buff_vdf_all_unstruct)

            !>>> Branch B: Structured Grid - Keep original logic <<<
            else 
                !--read vdf2d.dat
                allocate(solution_vdf_all(nx_all,ny_all,Nvel,Ndata))
                inquire(file="vdf2d.dat",exist=ext1)       !check whether the file exists
                if(.not. ext1) then
                    stop 'Restart error: file vdf2d.dat does not exist.'
                endif
                open(1,file='vdf2d.dat',form='unformatted')
                read(1)
                read(1) nx0,ny0,nvar0
                
                if(nx0/=nx_all .or. ny0/=ny_all .or. nvar0/=Nvel*Ndata) then
                    print *, 'nx0 =', nx0, ', ny0 =', ny0, ', nvar0 =', nvar0
                    stop 'Restart error: data in vdf2d.dat does not conform.'
                endif
                read(1) ((((solution_vdf_all(i,j,flat_k,k),i=1,nx_all),j=1,ny_all),flat_k=1,Nvel),k=1,Ndata)
                close(1)
                
                !--send to each processor
                do pj=0,NP_Y-1 
                    do pi=0,NP_X-1
                        np_send=pj*NP_X+pi
                        if(np_send.eq.0) then
                            do j1=1,ny
                                do i1=1,nx
                                    !Copy structured restart VDFs through flattened cell-center views.
                                    do flat_k=1,Nvel
                                        ctr(i1,j1)%h_vec(flat_k) = solution_vdf_all(i1,j1,flat_k,1)
                                        ctr(i1,j1)%b_vec(flat_k) = solution_vdf_all(i1,j1,flat_k,2)
                                    enddo
                                enddo
                            enddo
                        else
                            allocate(buff1_vdf(Block_XY(np_send)%NI,Block_XY(np_send)%NJ,Nvel,Ndata))
                            do j1=1,Block_XY(np_send)%NJ
                                do i1=1,Block_XY(np_send)%NI
                                    i=Block_XY(np_send)%Origin_I+i1-1
                                    j=Block_XY(np_send)%Origin_J+j1-1
                                    buff1_vdf(i1,j1,1:Nvel,1:Ndata)=solution_vdf_all(i,j,1:Nvel,1:Ndata)
                                enddo
                            enddo
                            tag=np_send*10      !tag equals to block index multiply 10
                            call MPI_Send(buff1_vdf,Block_XY(np_send)%NI*Block_XY(np_send)%NJ*Nvel*Ndata,MPI_KREAL,np_send,tag,MPI_COMM_WORLD,ierr)
                            deallocate(buff1_vdf)
                        endif
                    enddo
                enddo
                deallocate(solution_vdf_all)
            endif

        else 
            !============================================================
            !  Worker Process: Receive Data
            !============================================================
            
            if (Quadrature_method == Quad_Unstruct) then
                ! --- Receive unstructured data (Tag * 20) ---
                allocate(buff_vdf_local_unstruct(nx, ny, Nvel, Ndata))

                tag = (IP_Y*NP_X + IP_X) * 20
                call MPI_Recv(buff_vdf_local_unstruct, nx*ny*Nvel*Ndata, MPI_KREAL, 0, tag, MPI_COMM_WORLD, Status, ierr)

                do j1=1,ny
                    do i1=1,nx
                        !Copy unstructured restart VDFs through flattened cell-center views.
                        do k=1,Nvel
                            ctr(i1,j1)%h_vec(k) = buff_vdf_local_unstruct(i1, j1, k, 1)
                            ctr(i1,j1)%b_vec(k) = buff_vdf_local_unstruct(i1, j1, k, 2)
                        enddo
                    enddo
                enddo
                deallocate(buff_vdf_local_unstruct)
                
            else
                ! --- Receive structured grid data (Tag * 10) ---
                allocate(buff2_vdf(nx,ny,Nvel,Ndata))
                tag=(IP_Y*NP_X+IP_X)*10      !tag equals to block index multiply 10
                call MPI_Recv(buff2_vdf,nx*ny*Nvel*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,Status,ierr)
                do j1=1,ny
                    do i1=1,nx
                        !Copy structured restart VDFs through flattened cell-center views.
                        do flat_k=1,Nvel
                            ctr(i1,j1)%h_vec(flat_k) = buff2_vdf(i1,j1,flat_k,1)
                            ctr(i1,j1)%b_vec(flat_k) = buff2_vdf(i1,j1,flat_k,2)
                        enddo
                    enddo
                enddo
                deallocate(buff2_vdf)
            endif
        endif
        
        !Assign the initial flow field snapshot for both steady and DIS cases
        do j=1,ny
            do i=1,nx
                ctr(i,j)%conVars0(1:4) = ctr(i,j)%conVars(1:4)
                !Snapshot restart VDFs through flattened cell-center views.
                do k=1,Nvel
                    ctr(i,j)%h0_vec(k) = ctr(i,j)%h_vec(k)
                    ctr(i,j)%b0_vec(k) = ctr(i,j)%b_vec(k)
                enddo
            end do
        end do
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
    end subroutine InitFlow_Restart
    
    
    !---------------------------------------------------------
    !>Initialize inflow boundary
    !   ==> steady simulation: inflow.dat & invdf.dat
    !   ==> unsteady simulation: BaseFlow: flow2d.dat & vdf2d.dat; DIS: eev1.dat ... eevn.dat 
    !---------------------------------------------------------
    subroutine Init_inflow()
        integer(KINT)                                   :: ny0,nvar0,j,k,l,m,pj,np_send,j1,Ndata,tag
        integer(KINT)                                   :: flat_k !Flattened velocity-space index for inflow VDF copies
        real(KREAL), allocatable, dimension(:,:)        :: solution_all
        real(KREAL), allocatable, dimension(:,:)        :: buff1,buff2
        real(KREAL), allocatable, dimension(:,:,:)      :: solution_vdf_all !Inflow VDF cache with flat Nvel velocity dimension
        real(KREAL), allocatable, dimension(:,:,:)      :: buff1_vdf,buff2_vdf !MPI buffers with flat Nvel velocity dimension
        real(KREAL), allocatable, dimension(:)          :: H_vec,B_vec              !Pure 1D local Maxwellian arrays
        real(KREAL), allocatable, dimension(:,:,:)      :: ts_all,buff1_ts,buff2_ts
        real(KREAL)                                     :: prim(4)
        logical                                         :: ext1
        !character(len=20)                               :: char

        !Allocate inflow array for base flow
        allocate(inflow_conVars(ny,4))
        !Store inflow VDFs as pure 1D velocity-space caches; file/MPI buffers keep their legacy Nu,Nv layout.
        allocate(inflow_h(ny,Nvel),inflow_b(ny,Nvel))
        
        !Steady flow simulation: inflow array is read from file, only for base flow
        if(Iflag_DIS==0) then           
            !----------------read inflow.dat-------------------
            Ndata = 4
            if(myid==0) then
                inquire(file="inflow.dat",exist=ext1)
                if(.not. ext1) then
                    stop 'Error in read inflow: file inflow.dat does not exist.'
                endif
                open(1,file='inflow.dat',form='unformatted')
                read(1)
                read(1) ny0,nvar0
                if(ny0/=ny_all .or. nvar0/=Ndata) then
                    stop 'Error in read inflow: data in inflow.dat does not conform.'
                endif
                allocate(solution_all(ny_all,Ndata))
                read(1) ((solution_all(j,k),j=1,ny_all),k=1,Ndata)
                close(1)
                
                !--send to each processor
                do pj=0,NP_Y-1 
                    np_send=pj*NP_X
                    if(np_send.eq.0) then
                        do j1=1,ny
                            inflow_conVars(j1,1:Ndata) = solution_all(j1,1:Ndata)
                        enddo
                    else
                        allocate(buff1(Block_XY(np_send)%NJ,Ndata))
                        do j1=1,Block_XY(np_send)%NJ
                            j=Block_XY(np_send)%Origin_J+j1-1
                            buff1(j1,1:Ndata)=solution_all(j,1:Ndata)
                        enddo
                        tag=np_send     !tag equals to block index
                        call MPI_Send(buff1,Block_XY(np_send)%NJ*Ndata,MPI_KREAL,np_send,tag,MPI_COMM_WORLD,ierr)
                        deallocate(buff1)
                    endif
                enddo
            elseif(myid/=0 .and. IP_X==0) then
                allocate(buff2(ny,Ndata))
                tag=IP_Y*NP_X           !tag equals to block index
                call MPI_Recv(buff2,ny*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,Status,ierr)
                do j1=1,ny
                    inflow_conVars(j1,1:Ndata) = buff2(j1,1:Ndata)
                enddo
                deallocate(buff2)
            endif
        
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
            !--------------------read invdf.dat-----------------------
            Ndata = 2
        
            if(myID==0) then
                allocate(solution_vdf_all(ny_all,Nvel,Ndata))
                inquire(file="invdf.dat",exist=ext1)        !check whether the file exists
                if(.not. ext1) then                         !file not exist: vdf is initialized by Maxwell distribution
                    !Allocate pure 1D Maxwellian arrays for generated inflow VDFs.
                    allocate(H_vec(Nvel),B_vec(Nvel))
                    do j=1,ny_all
                        prim=GetPrimary(solution_all(j,1:4))
                        call DiscreteMaxwell_Vec(H_vec,B_vec,prim)
                        !Store generated inflow VDFs directly in flattened velocity order.
                        solution_vdf_all(j,1:Nvel,1)=H_vec(1:Nvel)
                        solution_vdf_all(j,1:Nvel,2)=B_vec(1:Nvel)
                    enddo
                    deallocate(H_vec,B_vec)
                else                                        !file exist: read the file    
                    open(1,file='invdf.dat',form='unformatted')
                    read(1)
                    read(1) ny0,nvar0
                    if(ny0/=ny_all .or. nvar0/=Nvel*Ndata) then
                        stop 'Error in read vdf: data in invdf.dat does not conform.'
                    endif
                    read(1) (((solution_vdf_all(j,flat_k,k),j=1,ny_all),flat_k=1,Nvel),k=1,Ndata)
                    close(1)
                endif
            
                !--send to each processor
                do pj=0,NP_Y-1 
                        np_send=pj*NP_X
                        if(np_send.eq.0) then
                            do j1=1,ny
                                !Copy root inflow VDFs as contiguous flattened velocity arrays.
                                inflow_h(j1,1:Nvel) = solution_vdf_all(j1,1:Nvel,1)
                                inflow_b(j1,1:Nvel) = solution_vdf_all(j1,1:Nvel,2)
                            enddo
                    else
                        allocate(buff1_vdf(Block_XY(np_send)%NJ,Nvel,Ndata))
                        do j1=1,Block_XY(np_send)%NJ
                            j=Block_XY(np_send)%Origin_J+j1-1
                            buff1_vdf(j1,1:Nvel,1:Ndata)=solution_vdf_all(j,1:Nvel,1:Ndata)
                        enddo
                        tag=np_send*10      !tag equals to block index multiply 10
                        call MPI_Send(buff1_vdf,Block_XY(np_send)%NJ*Nvel*Ndata,MPI_KREAL,np_send,tag,MPI_COMM_WORLD,ierr)
                        deallocate(buff1_vdf)
                    endif
                enddo
                deallocate(solution_all)
                deallocate(solution_vdf_all)
            elseif(myid/=0 .and. IP_X==0) then
                allocate(buff2_vdf(ny,Nvel,Ndata))
                tag=IP_Y*NP_X*10      !tag equals to block index multiply 10
                call MPI_Recv(buff2_vdf,ny*Nvel*Ndata,MPI_KREAL,0,tag,MPI_COMM_WORLD,Status,ierr)
                do j1=1,ny
                    !Copy received inflow VDFs through explicit velocity loops.
                    do flat_k=1,Nvel
                        inflow_h(j1,flat_k) = buff2_vdf(j1,flat_k,1)
                        inflow_b(j1,flat_k) = buff2_vdf(j1,flat_k,2)
                    enddo
                enddo
                deallocate(buff2_vdf)
            endif
        
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
        !Unsteady flow simulation: inflow array is initialized by restart flow, disturbance eigfun. is read from file   
        elseif(Iflag_DIS==1) then
            !Inflow array for base flow
            do j=1,ny
                inflow_conVars(j,1:4) = ctr(1,j)%conVars0(1:4)
                !Copy base inflow VDFs from flattened initial-state views.
                inflow_h(j,1:Nvel) = ctr(1,j)%h0_vec(1:Nvel)
                inflow_b(j,1:Nvel) = ctr(1,j)%b0_vec(1:Nvel)
            enddo
            
            !Inflow array for disturbance
            if(DISBC_type==DISBC_inflow_TS) then
                Ndata = 5
                allocate(inflow_tsr(ny,Ndata,Num_dis_ts),inflow_tsi(ny,Ndata,Num_dis_ts))       !rho, u, v, T, p
                if(myid==0) then
                    allocate(ts_all(ny_all,Ndata*2,Num_dis_ts))
                    do m=1,Num_dis_ts
                        inquire(file='eev'//char(ichar('0')+m)//'.dat',exist=ext1)
                        if(.not. ext1) then
                            stop 'Error in read eev: file eev.dat does not exist.'
                        endif
                        open(1,file='eev'//char(ichar('0')+m)//'.dat',form='unformatted')
                        read(1)
                        read(1) ny0,nvar0
                        if(ny0/=ny_all .or. nvar0/=Ndata) then
                            stop 'Error in read eev: data in eev.dat does not conform.'
                        endif
                        read(1) ((ts_all(j,l,m),j=1,ny_all),l=1,Ndata),((ts_all(j,l,m),j=1,ny_all),l=Ndata+1,Ndata*2)
                        close(1)
                    enddo
                    
                    !--send to each processor
                    do pj=0,NP_Y-1 
                        np_send=pj*NP_X
                        if(np_send.eq.0) then
                            do j1=1,ny
                                inflow_tsr(j1,1:Ndata,1:Num_dis_ts) = ts_all(j1,1:Ndata,1:Num_dis_ts) 
                                inflow_tsi(j1,1:Ndata,1:Num_dis_ts) = ts_all(j1,(Ndata+1):Ndata*2,1:Num_dis_ts)
                            enddo
                        else
                            allocate(buff1_ts(Block_XY(np_send)%NJ,Ndata*2,Num_dis_ts))
                            do j1=1,Block_XY(np_send)%NJ
                                j=Block_XY(np_send)%Origin_J+j1-1
                                buff1_ts(j1,1:Ndata*2,1:Num_dis_ts)=ts_all(j,1:Ndata*2,1:Num_dis_ts)
                            enddo
                            tag=np_send     !tag equals to block index
                            call MPI_Send(buff1_ts,Block_XY(np_send)%NJ*Ndata*2*Num_dis_ts,MPI_KREAL,np_send,tag,MPI_COMM_WORLD,ierr)
                            deallocate(buff1_ts)
                        endif
                    enddo
                    deallocate(ts_all)
                elseif(myid/=0 .and. IP_X==0) then
                    allocate(buff2_ts(ny,Ndata*2,Num_dis_ts))
                    tag=IP_Y*NP_X           !tag equals to block index
                    call MPI_Recv(buff2_ts,ny*Ndata*2*Num_dis_ts,MPI_KREAL,0,tag,MPI_COMM_WORLD,Status,ierr)
                    do j1=1,ny
                        inflow_tsr(j1,1:Ndata,1:Num_dis_ts) = buff2_ts(j1,1:Ndata,1:Num_dis_ts)
                        inflow_tsi(j1,1:Ndata,1:Num_dis_ts) = buff2_ts(j1,(Ndata+1):Ndata*2,1:Num_dis_ts)
                    enddo
                    deallocate(buff2_ts)
                endif

                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                
            endif
            
        endif
            
    end subroutine Init_inflow
    

    !---------------------------------------------------------
    !>fringe function
    !---------------------------------------------------------
    subroutine Init_fringe_x()
        implicit none

        !Allocate once and refresh on every accepted physical step.
        if(.not.allocated(xlamd)) allocate(xlamd(nx,ny))
        call FillFringeAtStep(iter,xlamd)
        
    end subroutine
    
end module InitializeFlow

!-----------------------------------------------------------------
!>Boundary conditions: set ghost cell values for physical boundary
!-----------------------------------------------------------------
!Notice: 1) isothermal wall is realized through Cal_FluxBoundary instead of ghost cell;
!        2) the present symmetry boundary only applies to x- or y-axis

module Boundary
    use flow_data
    use Flow_DIS
    implicit none
    
contains

    !--------------------------------------------------
    !>Return the free-wave number for an angular frequency expressed on the
    !>UGKS thermal-speed time scale, sqrt(2*R*T_inf)/L_ref.
    !--------------------------------------------------
    pure real(KREAL) function FreeWaveNumberUGKS(wave_name,omer,theta,aoa,ma_value,gamma_value) result(kinf)
        character(len=*), intent(in) :: wave_name
        real(KREAL), intent(in) :: omer,theta,aoa,ma_value,gamma_value
        real(KREAL) :: projected_mean_speed,thermal_sound_speed

        projected_mean_speed=ma_value*cosd(theta-aoa)
        thermal_sound_speed=sqrt(gamma_value/2.0_KREAL)
        select case(trim(adjustl(wave_name)))
        case('F')
            kinf=omer/(thermal_sound_speed*(projected_mean_speed+1.0_KREAL))
        case('S')
            kinf=omer/(thermal_sound_speed*(projected_mean_speed-1.0_KREAL))
        case('E','V')
            kinf=omer/(thermal_sound_speed*projected_mean_speed)
        case default
            kinf=0.0_KREAL
        end select
    end function FreeWaveNumberUGKS

    !--------------------------------------------------
    !>Mirror a ghost-cell VDF by reversing both velocity components in flat Nvel storage.
    !--------------------------------------------------
    subroutine SetMirrorGhostBoth(dstCell,srcCell)
        type(CellCenter), intent(inout) :: dstCell !Ghost cell receiving mirrored VDFs
        type(CellCenter), intent(in)    :: srcCell !Interior cell providing source VDFs
        integer(KINT)                   :: k       !Flat velocity-space index
        integer(KINT)                   :: vel_i   !Structured u-index recovered from k
        integer(KINT)                   :: vel_j   !Structured v-index recovered from k
        integer(KINT)                   :: src_k   !Flat source index after u/v reflection

        !Reset VDF slopes once per ghost cell and mirror h/b with the existing flat ordering.
        dstCell%sh_vec = 0.0_KREAL
        dstCell%sb_vec = 0.0_KREAL
        do k=1,Nvel
            vel_i = mod(k-1,Nu)+1
            vel_j = (k-1)/Nu+1
            src_k = (Nu-vel_i+1)+(Nv-vel_j)*Nu
            dstCell%h_vec(k) = srcCell%h_vec(src_k)
            dstCell%b_vec(k) = srcCell%b_vec(src_k)
        enddo
    end subroutine SetMirrorGhostBoth

    !--------------------------------------------------
    !>Mirror a ghost-cell VDF by reversing the u velocity component in flat Nvel storage.
    !--------------------------------------------------
    subroutine SetMirrorGhostU(dstCell,srcCell)
        type(CellCenter), intent(inout) :: dstCell !Ghost cell receiving mirrored VDFs
        type(CellCenter), intent(in)    :: srcCell !Interior cell providing source VDFs
        integer(KINT)                   :: k       !Flat velocity-space index
        integer(KINT)                   :: vel_i   !Structured u-index recovered from k
        integer(KINT)                   :: vel_j   !Structured v-index recovered from k
        integer(KINT)                   :: src_k   !Flat source index after u reflection

        !Reset VDF slopes once per ghost cell and mirror h/b with the existing flat ordering.
        dstCell%sh_vec = 0.0_KREAL
        dstCell%sb_vec = 0.0_KREAL
        do k=1,Nvel
            vel_i = mod(k-1,Nu)+1
            vel_j = (k-1)/Nu+1
            src_k = (Nu-vel_i+1)+(vel_j-1)*Nu
            dstCell%h_vec(k) = srcCell%h_vec(src_k)
            dstCell%b_vec(k) = srcCell%b_vec(src_k)
        enddo
    end subroutine SetMirrorGhostU

    !--------------------------------------------------
    !>Mirror a ghost-cell VDF by reversing the v velocity component in flat Nvel storage.
    !--------------------------------------------------
    subroutine SetMirrorGhostV(dstCell,srcCell)
        type(CellCenter), intent(inout) :: dstCell !Ghost cell receiving mirrored VDFs
        type(CellCenter), intent(in)    :: srcCell !Interior cell providing source VDFs
        integer(KINT)                   :: k       !Flat velocity-space index
        integer(KINT)                   :: vel_i   !Structured u-index recovered from k
        integer(KINT)                   :: vel_j   !Structured v-index recovered from k
        integer(KINT)                   :: src_k   !Flat source index after v reflection

        !Unstructured velocity grids already provide a direct flat v-symmetry lookup table.
        dstCell%sh_vec = 0.0_KREAL
        dstCell%sb_vec = 0.0_KREAL
        if (Quadrature_method == Quad_Unstruct) then
            do k=1,Nvel
                src_k = SymMap_V(k)
                dstCell%h_vec(k) = srcCell%h_vec(src_k)
                dstCell%b_vec(k) = srcCell%b_vec(src_k)
            enddo
        else
            do k=1,Nvel
                vel_i = mod(k-1,Nu)+1
                vel_j = (k-1)/Nu+1
                src_k = vel_i+(Nv-vel_j)*Nu
                dstCell%h_vec(k) = srcCell%h_vec(src_k)
                dstCell%b_vec(k) = srcCell%b_vec(src_k)
            enddo
        endif
    end subroutine SetMirrorGhostV

    !--------------------------------------------------
    !>Set ghost cell of physical boundary
    !--------------------------------------------------
    subroutine Set_Ghost_PhysBoundary()   
        integer(KINT)               :: dirc_index,face_index, bc_index
        
        !---------------general geometry---------------------
            !              North
            !          ------------
            !          |           |
            !   West   |           |   East
            !          |           |
            !          ------------
            !              South
        !----------------------------------------------------
        if(IP_X==0) then    !left (west) boundary
            bc_index =1; dirc_index=Idirc; face_index=Lside
            select case(bc_type(bc_index))
            case(BCconst_wall_insulate)
                call bc_insulate(dirc_index,face_index)
            case(BCconst_symmetry)
                call bc_symmetry(dirc_index,face_index)
            case(BCconst_freestream)
                call bc_freestream(dirc_index,face_index)
            case(BCconst_extrapolation)
                call bc_extrapolation(dirc_index,face_index)
            case(BCconst_inflow)
                call bc_inflow()
            case(BCconst_periodic)
                !The ghost valued is updated through MPI_sendrecv, and nothing is needed here.
            case(BCconst_wall_isothermal)
                !isothermal wall: don't need ghost cell value
            end select
        endif
        
        if(IP_X==NP_X-1) then    !right (east) boundary
            bc_index =2; dirc_index=Idirc; face_index=Rside
            select case(bc_type(bc_index))
            case(BCconst_wall_insulate)
                call bc_insulate(dirc_index,face_index)
            case(BCconst_symmetry)
                call bc_symmetry(dirc_index,face_index)
            case(BCconst_freestream)
                call bc_freestream(dirc_index,face_index)
            case(BCconst_extrapolation)
                call bc_extrapolation(dirc_index,face_index)
            case(BCconst_inflow)
                write(*,'(a)') "Error in boundary setup: Inflow boundary only applied to left boundary, stop!!!"
                stop
            case(BCconst_periodic)
                !The ghost valued is updated through MPI_sendrecv, and nothing is needed here.
            case(BCconst_wall_isothermal)
                !isothermal wall: don't need ghost cell value
            end select
        endif
        
        if(IP_Y==0) then    !bottom (south) boundary
            bc_index =3; dirc_index=Jdirc; face_index=Lside
            select case(bc_type(bc_index))
            case(BCconst_wall_insulate)
                call bc_insulate(dirc_index,face_index)
            case(BCconst_symmetry)
                call bc_symmetry(dirc_index,face_index)
            case(BCconst_freestream)
                call bc_freestream(dirc_index,face_index)
            case(BCconst_extrapolation)
                call bc_extrapolation(dirc_index,face_index)
            case(BCconst_inflow)
                write(*,'(a)') "Error in boundary setup: Inflow boundary only applied to left boundary, stop!!!"
                stop
            case(BCconst_periodic)
                !The ghost valued is updated through MPI_sendrecv, and nothing is needed here.
            case(BCconst_wall_isothermal)
                !isothermal wall: don't need ghost cell value
            end select
        endif
        
        if(IP_Y==NP_Y-1) then    !upper (north) boundary
            bc_index =4; dirc_index=Jdirc; face_index=Rside
            select case(bc_type(bc_index))
            case(BCconst_wall_insulate)
                call bc_insulate(dirc_index,face_index)
            case(BCconst_symmetry)
                call bc_symmetry(dirc_index,face_index)
            case(BCconst_freestream)
                call bc_freestream(dirc_index,face_index)
            case(BCconst_extrapolation)
                call bc_extrapolation(dirc_index,face_index)
            case(BCconst_inflow)
                write(*,'(a)') "Error in boundary setup: Inflow boundary only applied to left boundary, stop!!!"
                stop
            case(BCconst_periodic)
                !The ghost valued is updated through MPI_sendrecv, and nothing is needed here.
            case(BCconst_wall_isothermal)
                !isothermal wall: don't need ghost cell value
            end select
        endif
        
        !DIS simulation: re-write ghost boundary value
        if(Iflag_DIS==1) then
            select case(DISBC_type)
            case(DISBC_inflow_TS)
                if(IP_X==0) call bc_inflow_DIS()        !only apply to left (inflow) boundary
            case(DISBC_free_wave)
                if((IP_X==0 .and. bc_type(1)==BCconst_freestream) .or. (IP_Y==NP_Y-1 .and. bc_type(4)==BCconst_freestream)) then
                    call bc_freestream_DIS()            !only apply to left or upper (farfield) boundary
                endif
            end select
        endif

    end subroutine Set_Ghost_PhysBoundary
    

    !-------------------------------------------------------------------
    
    
    !-------------------------------------------------------------------
    subroutine bc_extrapolation(dirc,side)
        integer,intent(in) :: dirc,side
        integer(KINT) :: i,j,k !Cell loop indices and flattened velocity-space index
        
        select case(dirc)
        case(Idirc)
            if(side==Lside) then
                do j=1,ny
                    do i=1,ghost
                        ctr(1-i,j)%conVars = ctr(1,j)%conVars
                        !Copy extrapolated VDFs through flattened ghost-cell views.
                        do k=1,Nvel
                            ctr(1-i,j)%h_vec(k) = ctr(1,j)%h_vec(k)
                            ctr(1-i,j)%b_vec(k) = ctr(1,j)%b_vec(k)
                        enddo
                        ctr(1-i,j)%sh_vec = ctr(1,j)%sh_vec
                        ctr(1-i,j)%sb_vec = ctr(1,j)%sb_vec
                    enddo
                end do
            elseif(side==Rside) then
                do j=1,ny
                    do i=1,ghost
                        ctr(nx+i,j)%conVars = ctr(nx,j)%conVars
                        !Copy extrapolated VDFs through flattened ghost-cell views.
                        do k=1,Nvel
                            ctr(nx+i,j)%h_vec(k) = ctr(nx,j)%h_vec(k)
                            ctr(nx+i,j)%b_vec(k) = ctr(nx,j)%b_vec(k)
                        enddo
                        ctr(nx+i,j)%sh_vec = ctr(nx,j)%sh_vec
                        ctr(nx+i,j)%sb_vec = ctr(nx,j)%sb_vec
                    enddo
                end do
            endif
        case(Jdirc)
            if(side==Lside) then
                do i=1,nx
                    do j=1,ghost
                        ctr(i,1-j)%conVars = ctr(i,1)%conVars
                        !Copy extrapolated VDFs through flattened ghost-cell views.
                        do k=1,Nvel
                            ctr(i,1-j)%h_vec(k) = ctr(i,1)%h_vec(k)
                            ctr(i,1-j)%b_vec(k) = ctr(i,1)%b_vec(k)
                        enddo
                        ctr(i,1-j)%sh_vec = ctr(i,1)%sh_vec
                        ctr(i,1-j)%sb_vec = ctr(i,1)%sb_vec
                    enddo
                end do
            elseif(side==Rside) then
                do i=1,nx
                    do j=1,ghost
                        ctr(i,ny+j)%conVars = ctr(i,ny)%conVars
                        !Copy extrapolated VDFs through flattened ghost-cell views.
                        do k=1,Nvel
                            ctr(i,ny+j)%h_vec(k) = ctr(i,ny)%h_vec(k)
                            ctr(i,ny+j)%b_vec(k) = ctr(i,ny)%b_vec(k)
                        enddo
                        ctr(i,ny+j)%sh_vec = ctr(i,ny)%sh_vec
                        ctr(i,ny+j)%sb_vec = ctr(i,ny)%sb_vec
                    enddo
                end do
            endif
        end select
        
    end subroutine bc_extrapolation

    
    !-------------------------------------------------------------------
    subroutine bc_insulate(dirc,side)           !no-slip and insulate boundary

        integer,intent(in) :: dirc,side
        integer(KINT) :: i,j !Cell and ghost-layer loop indices
    
        if(dirc==Idirc) then
            if(side==Lside) then
                do j=1,ny
                    do i=1,ghost
                        ctr(1-i,j)%conVars(1) = ctr(1+i-1,j)%conVars(1)
                        ctr(1-i,j)%conVars(2) = -ctr(1+i-1,j)%conVars(2)
                        ctr(1-i,j)%conVars(3) = -ctr(1+i-1,j)%conVars(3)
                        ctr(1-i,j)%conVars(4) = ctr(1+i-1,j)%conVars(4)
    
                        !Mirror ghost VDFs by reversing both velocity components in flat Nvel storage.
                        call SetMirrorGhostBoth(ctr(1-i,j),ctr(1+i-1,j))
                    enddo
                enddo
            else if(side==Rside) then
                do j=1,ny
                    do i=1,ghost
                        ctr(nx+i,j)%conVars(1) = ctr(nx-i+1,j)%conVars(1)
                        ctr(nx+i,j)%conVars(2) = -ctr(nx-i+1,j)%conVars(2)
                        ctr(nx+i,j)%conVars(3) = -ctr(nx-i+1,j)%conVars(3)
                        ctr(nx+i,j)%conVars(4) = ctr(nx-i+1,j)%conVars(4)
    
                        !Mirror ghost VDFs by reversing both velocity components in flat Nvel storage.
                        call SetMirrorGhostBoth(ctr(nx+i,j),ctr(nx-i+1,j))
                    enddo
                enddo
            end if
        else if(dirc==Jdirc) then
            if(side==Lside) then
                do i=1,nx
                    do j=1,ghost
                        ctr(i,1-j)%conVars(1) = ctr(i,1+j-1)%conVars(1)
                        ctr(i,1-j)%conVars(2) = -ctr(i,1+j-1)%conVars(2)
                        ctr(i,1-j)%conVars(3) = -ctr(i,1+j-1)%conVars(3)
                        ctr(i,1-j)%conVars(4) = ctr(i,1+j-1)%conVars(4)
    
                        !Mirror ghost VDFs by reversing both velocity components in flat Nvel storage.
                        call SetMirrorGhostBoth(ctr(i,1-j),ctr(i,1+j-1))
                    enddo
                enddo
            else if(side==Rside) then
                do i=1,nx
                    do j=1,ghost
                        ctr(i,ny+j)%conVars(1) = ctr(i,ny-j+1)%conVars(1)
                        ctr(i,ny+j)%conVars(2) = -ctr(i,ny-j+1)%conVars(1)
                        ctr(i,ny+j)%conVars(3) = -ctr(i,ny-j+1)%conVars(1)
                        ctr(i,ny+j)%conVars(4) = ctr(i,ny-j+1)%conVars(4)
    
                        !Mirror ghost VDFs by reversing both velocity components in flat Nvel storage.
                        call SetMirrorGhostBoth(ctr(i,ny+j),ctr(i,ny-j+1))
                    enddo
                enddo
            end if
        end if
    
    end subroutine bc_insulate
        
        
        !-------------------------------------------------------------------
    subroutine bc_symmetry(dirc,side)           !symmetry along x- or y-axis
    
    integer,intent(in) :: dirc,side
    integer(KINT) :: i,j !Cell and ghost-layer loop indices
    
    ! Check for unstructured mesh
    !if (dirc == Idirc .and. Quadrature_method == Quad_Unstruct) then
    !    print *, "Error: Left/Right Symmetry BC is not supported for this Unstructured Mesh!"
    !    stop
    !endif
    
    if(dirc==Idirc .and. side==Lside) then
        !BoundaryFace(1): firstly, judge symmetry along which axis
        if(abs(vface(1,1)%cosy)<1e-12 .and. abs(vface(1,ny)%cosy)<1e-12) then   !normal vector (1, 0) or (-1, 0): symmeter along y axis
            do j=1,ny
                do i=1,GHOST
                    ctr(1-i,j)%conVars(1) = ctr(1+i-1,j)%conVars(1)
                    ctr(1-i,j)%conVars(2) = -ctr(1+i-1,j)%conVars(2)    !Only reverse u velocity
                    ctr(1-i,j)%conVars(3) = ctr(1+i-1,j)%conVars(3) 
                    ctr(1-i,j)%conVars(4) = ctr(1+i-1,j)%conVars(4)
                    !Mirror the ghost VDF in u through the flat Nvel helper.
                    call SetMirrorGhostU(ctr(1-i,j),ctr(1+i-1,j))
                end do
            end do
        elseif(abs(vface(1,1)%cosx)<1e-12 .and. abs(vface(1,ny)%cosx)<1e-12) then   !normal vector (0, 1) or (0, -1): symmeter along x axis
            do j=1,ny
                do i=1,GHOST
                    ctr(1-i,j)%conVars(1) = ctr(1+i-1,j)%conVars(1)
                    ctr(1-i,j)%conVars(2) = ctr(1+i-1,j)%conVars(2) 
                    ctr(1-i,j)%conVars(3) = -ctr(1+i-1,j)%conVars(3)        !Only reverse v velocity 
                    ctr(1-i,j)%conVars(4) = ctr(1+i-1,j)%conVars(4)
                    !Mirror the ghost VDF in v through the flat Nvel helper.
                    call SetMirrorGhostV(ctr(1-i,j),ctr(1+i-1,j))
                end do
            end do
        else
            write(*,'(a)') 'Warning: the boundary face is not symmetry along x- or y-axis.' 
        endif
        
    elseif(dirc==Idirc .and. side==Rside) then
        !BoundaryFace(2): firstly, judge symmetry along which axis
        if(abs(vface(nx+1,1)%cosy)<1e-12 .and. abs(vface(nx+1,ny)%cosy)<1e-12) then   !normal vector (1, 0) or (-1, 0): symmeter along y axis
            do j=1,ny
                do i=1,GHOST
                    ctr(nx+i,j)%conVars(1) = ctr(nx-i+1,j)%conVars(1)
                    ctr(nx+i,j)%conVars(2) = -ctr(nx-i+1,j)%conVars(2)    !Only reverse u velocity
                    ctr(nx+i,j)%conVars(3) = ctr(nx-i+1,j)%conVars(3) 
                    ctr(nx+i,j)%conVars(4) = ctr(nx-i+1,j)%conVars(4)
                    !Mirror the ghost VDF in u through the flat Nvel helper.
                    call SetMirrorGhostU(ctr(nx+i,j),ctr(nx-i+1,j))
                end do
            end do
        elseif(abs(vface(nx+1,1)%cosx)<1e-12 .and. abs(vface(nx+1,ny)%cosx)<1e-12) then   !normal vector (0, 1) or (0, -1): symmeter along x axis
            do j=1,ny
                do i=1,GHOST
                    ctr(nx+i,j)%conVars(1) = ctr(nx-i+1,j)%conVars(1)
                    ctr(nx+i,j)%conVars(2) = ctr(nx-i+1,j)%conVars(2)
                    ctr(nx+i,j)%conVars(3) = -ctr(nx-i+1,j)%conVars(3)  !Only reverse v velocity
                    ctr(nx+i,j)%conVars(4) = ctr(nx-i+1,j)%conVars(4)
                    !Mirror the ghost VDF in v through the flat Nvel helper.
                    call SetMirrorGhostV(ctr(nx+i,j),ctr(nx-i+1,j))
                end do
            end do
        else
            write(*,'(a)') 'Warning: the boundary face is not symmetry along x- or y-axis.' 
        endif
        ! Jdirc means up and bottom boundary, Lside means bottom boundary
    elseif(dirc==Jdirc .and. side==Lside) then
        !BoundaryFace(3): firstly, judge symmetry along which axis
        if(abs(hface(1,1)%cosy)<1e-12 .and. abs(hface(nx,1)%cosy)<1e-12) then   !normal vector (1, 0) or (-1, 0): symmeter along y axis
            do i=1,nx
                do j=1,GHOST
                    ctr(i,1-j)%conVars(1) = ctr(i,1+j-1)%conVars(1)
                    ctr(i,1-j)%conVars(2) = -ctr(i,1+j-1)%conVars(2)     !Only reverse u velocity
                    ctr(i,1-j)%conVars(3) = ctr(i,1+j-1)%conVars(3) 
                    ctr(i,1-j)%conVars(4) = ctr(i,1+j-1)%conVars(4)
                    !Mirror the ghost VDF in u through the flat Nvel helper.
                    call SetMirrorGhostU(ctr(i,1-j),ctr(i,1+j-1))
                end do
            end do
        elseif(abs(hface(1,1)%cosx)<1e-12 .and. abs(hface(nx,1)%cosx)<1e-12) then   !normal vector (0, 1) or (0, -1): symmeter along x axis
            do i=1,nx
                do j=1,GHOST
                    ctr(i,1-j)%conVars(1) = ctr(i,1+j-1)%conVars(1)
                    ctr(i,1-j)%conVars(2) = ctr(i,1+j-1)%conVars(2)
                    ctr(i,1-j)%conVars(3) = -ctr(i,1+j-1)%conVars(3) !Only reverse v velocity
                    ctr(i,1-j)%conVars(4) = ctr(i,1+j-1)%conVars(4)
                    !Mirror the ghost VDF in v through the flat Nvel helper.
                    call SetMirrorGhostV(ctr(i,1-j),ctr(i,1+j-1))
                end do
            end do
        else
            write(*,'(a)') 'Warning: the boundary face is not symmetry along x- or y-axis.' 
        endif
        ! Jdirc means up and bottom boundary, Rside means upper boundary
    elseif(dirc==Jdirc .and. side==Rside) then
        !BoundaryFace(3): firstly, judge symmetry along which axis
        if(abs(hface(1,ny+1)%cosy)<1e-12 .and. abs(hface(nx,ny+1)%cosy)<1e-12) then   !normal vector (1, 0) or (-1, 0): symmeter along y axis
            do i=1,nx
                do j=1,GHOST
                    ctr(i,ny+j)%conVars(1) = ctr(i,ny-j+1)%conVars(1)
                    ctr(i,ny+j)%conVars(2) = -ctr(i,ny-j+1)%conVars(2)      !Only reverse u velocity
                    ctr(i,ny+j)%conVars(3) = ctr(i,ny-j+1)%conVars(3) 
                    ctr(i,ny+j)%conVars(4) = ctr(i,ny-j+1)%conVars(4)
                    !Mirror the ghost VDF in u through the flat Nvel helper.
                    call SetMirrorGhostU(ctr(i,ny+j),ctr(i,ny-j+1))
                end do
            end do
        elseif(abs(hface(1,ny+1)%cosx)<1e-12 .and. abs(hface(nx,ny+1)%cosx)<1e-12) then   !normal vector (0, 1) or (0, -1): symmeter along x axis
            do i=1,nx
                do j=1,GHOST
                    ctr(i,ny+j)%conVars(1) = ctr(i,ny-j+1)%conVars(1)
                    ctr(i,ny+j)%conVars(2) = ctr(i,ny-j+1)%conVars(2)
                    ctr(i,ny+j)%conVars(3) = -ctr(i,ny-j+1)%conVars(3) !Only reverse v velocity
                    ctr(i,ny+j)%conVars(4) = ctr(i,ny-j+1)%conVars(4)
                    !Mirror the ghost VDF in v through the flat Nvel helper.
                    call SetMirrorGhostV(ctr(i,ny+j),ctr(i,ny-j+1))
                end do
            end do
        else
            write(*,'(a)') 'Warning: the boundary face is not symmetry along x- or y-axis.' 
        endif
        
    endif
    end subroutine bc_symmetry
    
    
    !-------------------------------------------------------------------
    subroutine bc_freestream(dirc,side)
        use Tools,only              :  GetConserved,DiscreteMaxwell_Vec
        integer,intent(in)          :: dirc,side
        integer(KINT)               :: i,j,k
        real(KREAL)                 :: conVars(4) !Conservative variables
        real(KREAL), allocatable, dimension(:) :: H_vec,B_vec !Pure 1D freestream distributions

        !Allocate pure 1D Maxwellian arrays for freestream ghost cells.
        allocate(H_vec(Nvel))
        allocate(B_vec(Nvel))

        !Get conservative variables and Maxwellian distribution function
        conVars = GetConserved(Prim_inf)
        call DiscreteMaxwell_Vec(H_vec,B_vec,Prim_inf)

        select case(dirc)
        case(Idirc)
        if(side==Lside) then
            do j=1,ny
                do i=1,GHOST
                    ctr(1-i,j)%conVars = conVars
                    !Copy freestream VDFs through flattened ghost-cell views.
                    do k=1,Nvel
                        ctr(1-i,j)%h_vec(k) = H_vec(k)
                        ctr(1-i,j)%b_vec(k) = B_vec(k)
                    enddo
                    ctr(1-i,j)%sh_vec = 0.0
                    ctr(1-i,j)%sb_vec = 0.0
                enddo
            enddo
        elseif(side==Rside) then
            do j=1,ny
                do i=1,GHOST
                    ctr(nx+i,j)%conVars = conVars
                    !Copy freestream VDFs through flattened ghost-cell views.
                    do k=1,Nvel
                        ctr(nx+i,j)%h_vec(k) = H_vec(k)
                        ctr(nx+i,j)%b_vec(k) = B_vec(k)
                    enddo
                    ctr(nx+i,j)%sh_vec = 0.0
                    ctr(nx+i,j)%sb_vec = 0.0
                enddo
            enddo
        endif
        case(Jdirc)
            if(side==Lside) then
                do i=1,nx
                    do j=1,GHOST
                        ctr(i,1-j)%conVars = conVars
                        !Copy freestream VDFs through flattened ghost-cell views.
                        do k=1,Nvel
                            ctr(i,1-j)%h_vec(k) = H_vec(k)
                            ctr(i,1-j)%b_vec(k) = B_vec(k)
                        enddo
                        ctr(i,1-j)%sh_vec = 0.0
                        ctr(i,1-j)%sb_vec = 0.0
                    enddo
                enddo
            elseif(side==Rside) then
                do i=1,nx
                    do j=1,GHOST
                        ctr(i,ny+j)%conVars = conVars
                        !Copy freestream VDFs through flattened ghost-cell views.
                        do k=1,Nvel
                            ctr(i,ny+j)%h_vec(k) = H_vec(k)
                            ctr(i,ny+j)%b_vec(k) = B_vec(k)
                        enddo
                        ctr(i,ny+j)%sh_vec = 0.0
                        ctr(i,ny+j)%sb_vec = 0.0
                    enddo
                enddo
            endif
        end select
        
        !Deallocation
        deallocate(H_vec)
        deallocate(B_vec)
    end subroutine bc_freestream

    
    !-------------------------------------------------------------------
    subroutine bc_inflow()
        integer(KINT)               :: i,j,k !Cell and flattened velocity-space loop indices

        do j=1,ny
            do i=1,GHOST
                ctr(1-i,j)%conVars = inflow_conVars(j,1:4)
                !Copy stored inflow VDFs through flattened ghost-cell views.
                do k=1,Nvel
                    ctr(1-i,j)%h_vec(k) = inflow_h(j,k)
                    ctr(1-i,j)%b_vec(k) = inflow_b(j,k)
                enddo
                ctr(1-i,j)%sh_vec = 0.0
                ctr(1-i,j)%sb_vec = 0.0
            enddo
        enddo
        
    end subroutine bc_inflow
    
    !--------------------------------------------------
    !>DIS boundary for inflow TS wave
    !       : only apply to left (inflow) boundary
    !--------------------------------------------------
    subroutine bc_inflow_DIS()
        use Tools,only              :  GetPrimary,GetConserved,DiscreteMaxwell_Vec
        real(KREAL)                 :: conVars0(4),prim0(4),prim_dis(4),prim(4),conVars(4)
        real(KREAL)                 :: phasetmp
        integer(KINT)               :: i,j,k,l,m !Cell, velocity, and disturbance loop indices
        integer(KINT)               :: flat_k !Flattened velocity-space index
        real(KREAL), allocatable, dimension(:) :: H0_vec,B0_vec,H_vec,B_vec !Pure 1D Maxwellian arrays

        !Allocate pure 1D Maxwellian arrays for disturbed inflow ghost cells.
        allocate(H0_vec(Nvel),B0_vec(Nvel),H_vec(Nvel),B_vec(Nvel))
        
        do j=1,ny
            do i=1,GHOST
                conVars0    = inflow_conVars(j,1:4)
                prim0       = GetPrimary(conVars0)                  !r0, u0, v0, 1/T0
                call DiscreteMaxwell_Vec(H0_vec,B0_vec,prim0)     !H0,B0
                prim0(4)    = 1.d0/prim0(4)                         !r0, u0, v0, T0
                prim_dis = 0.d0
                do m=1,Num_dis_ts
                    phasetmp = -wave_ts(m)%omer*time_dis
                    do k=1,4
                        prim_dis(k)=prim_dis(k)+wave_ts(m)%amp*(inflow_tsr(j,k,m)*dcos(phasetmp)-inflow_tsi(j,k,m)*dsin(phasetmp))     !r', u', v', T'
                    enddo
                enddo
                prim        = prim0+prim_dis                        !r, u, v, T
                prim(4)     = 1.d0/prim(4)                          !r, u, v, 1/T
                conVars     = GetConserved(prim)
                call DiscreteMaxwell_Vec(H_vec,B_vec,prim)        !H, B

                ctr(1-i,j)%conVars = conVars
                !Copy disturbed inflow VDFs through flattened ghost-cell views.
                do flat_k=1,Nvel
                    ctr(1-i,j)%h_vec(flat_k) = inflow_h(j,flat_k)+H_vec(flat_k)-H0_vec(flat_k)
                    ctr(1-i,j)%b_vec(flat_k) = inflow_b(j,flat_k)+B_vec(flat_k)-B0_vec(flat_k)
                enddo
                ctr(1-i,j)%sh_vec = 0.0
                ctr(1-i,j)%sb_vec = 0.0
            enddo
        enddo
        
        deallocate(H0_vec,B0_vec,H_vec,B_vec)

    end subroutine bc_inflow_DIS
    
    
    !--------------------------------------------------
    !>DIS boundary for free-stream sound/entropy/vortex wave 
    !       : only apply to left or upper (farfield) boundary
    !--------------------------------------------------
    subroutine bc_freestream_DIS()
        use Tools,only              :  GetPrimary,GetConserved,DiscreteMaxwell_Vec
        integer(KINT)               :: i,j,k,m
        real(KREAL)                 :: phasetmp,omer,theta,kinf,kx,ky,xx,yy
        real(KREAL)                 :: prim0(4),prim_dis(4),prim(4),conVars(4)
        real(KREAL), allocatable    :: wave_coef(:,:)      !Eigenfunction: r, u, v, T, p
        real(KREAL), allocatable, dimension(:) :: H_vec,B_vec  !Pure 1D freestream disturbance distributions

        !Allocation
        allocate(H_vec(Nvel),B_vec(Nvel))
        allocate(wave_coef(7,Num_dis_free))     !r, u, v, T, p
        
        !Free-stream state
        Prim0 = Prim_inf
        Prim0(4) = 1.d0/Prim0(4)
        
        !Free-stream disturbance
        do m=1,Num_dis_free
            omer = wave_free(m)%omer
            theta = wave_free(m)%theta
            select case(trim(adjustl(wave_free(m)%name)))
            case('F')   !fast sound wave
                kinf           = FreeWaveNumberUGKS('F',omer,theta,AOA,Ma,gamma)
                wave_coef(1,m) = 2.d0/gamma                         !r
                wave_coef(2,m) = sqrt(2.d0/gamma)*cosd(theta)       !u
                wave_coef(3,m) = sqrt(2.d0/gamma)*sind(theta)       !v
                wave_coef(4,m) = 2.d0*(gamma-1.d0)/gamma            !T
                wave_coef(5,m) = 1.d0                               !p
                wave_coef(6,m) = kinf*cosd(theta)                   !kx
                wave_coef(7,m) = kinf*sind(theta)                   !ky
            case('S')   !Slow sound wave
                kinf           = FreeWaveNumberUGKS('S',omer,theta,AOA,Ma,gamma)
                wave_coef(1,m) = 2.d0/gamma                         !r
                wave_coef(2,m) = -sqrt(2.d0/gamma)*cosd(theta)      !u
                wave_coef(3,m) = -sqrt(2.d0/gamma)*sind(theta)      !v
                wave_coef(4,m) = 2.d0*(gamma-1.d0)/gamma            !T
                wave_coef(5,m) = 1.d0                               !p
                wave_coef(6,m) = kinf*cosd(theta)                   !kx
                wave_coef(7,m) = kinf*sind(theta)                   !ky
            case('E')
                kinf           = FreeWaveNumberUGKS('E',omer,theta,AOA,Ma,gamma)
                wave_coef(1,m) = 1.d0                               !r
                wave_coef(2,m) = 0.d0                               !u
                wave_coef(3,m) = 0.d0                               !v
                wave_coef(4,m) = -1.d0                              !T
                wave_coef(5,m) = 0.d0                               !p
                wave_coef(6,m) = kinf*cosd(theta)                   !kx
                wave_coef(7,m) = kinf*sind(theta)                   !ky
            case('V')
                kinf           = FreeWaveNumberUGKS('V',omer,theta,AOA,Ma,gamma)
                wave_coef(1,m) =  0.d0                              !r
                wave_coef(2,m) =  -tand(theta)                      !u
                wave_coef(3,m) =  1.d0                              !v
                wave_coef(4,m) =  0.d0                              !T
                wave_coef(5,m) =  0.d0                              !p
                wave_coef(6,m) = kinf*cosd(theta)                   !kx
                wave_coef(7,m) = kinf*sind(theta)                   !ky
            end select
        enddo
        
        if(IP_X==0) then                !Left boundary is freestream (flat-plate case)  
            do j=1,ny
                do i=1,GHOST
                    xx=ctr(1,j)%x;    yy=ctr(1,j)%y
                    prim_dis = 0.d0
                    do m=1,Num_dis_free
                        omer=wave_free(m)%omer;  kx=wave_coef(6,m);    ky=wave_coef(7,m)
                        phasetmp = kx*xx+ky*yy-omer*time_dis
                        do k=1,4
                           prim_dis(k)=prim_dis(k)+wave_free(m)%amp * wave_coef(k,m) * dcos(phasetmp)   !r', u', v', T'
                        enddo
                    enddo
                    prim        = prim0+prim_dis                        !r, u, v, T
                    prim(4)     = 1.d0/prim(4)                          !r, u, v, 1/T
                    conVars     = GetConserved(prim)
                    call DiscreteMaxwell_Vec(H_vec,B_vec,prim)        !H, B

                    ctr(1-i,j)%conVars = conVars
                    !Copy disturbed freestream VDFs through flattened ghost-cell views.
                    do k=1,Nvel
                        ctr(1-i,j)%h_vec(k) = H_vec(k)
                        ctr(1-i,j)%b_vec(k) = B_vec(k)
                    enddo
                    ctr(1-i,j)%sh_vec = 0.0
                    ctr(1-i,j)%sb_vec = 0.0
                enddo
            enddo
            
        endif

        if(IP_Y==NP_Y-1) then            !Upper boundary is freestream        
            do i=1,nx
                do j=1,GHOST
                    xx=ctr(i,ny)%x;    yy=ctr(i,ny)%y
                    prim_dis = 0.d0
                    do m=1,Num_dis_free
                        omer=wave_free(m)%omer;  kx=wave_coef(6,m);    ky=wave_coef(7,m)
                        phasetmp = kx*xx+ky*yy-omer*time_dis
                        do k=1,4
                           prim_dis(k)=prim_dis(k)+wave_free(m)%amp * wave_coef(k,m) * dcos(phasetmp)  !r', u', v', T'
                        enddo
                    enddo
                    prim        = prim0+prim_dis                        !r, u, v, T
                    prim(4)     = 1.d0/prim(4)                          !r, u, v, 1/T
                    conVars     = GetConserved(prim)
                    call DiscreteMaxwell_Vec(H_vec,B_vec,prim)        !H, B

                    ctr(i,ny+j)%conVars = conVars
                    !Copy disturbed freestream VDFs through flattened ghost-cell views.
                    do k=1,Nvel
                        ctr(i,ny+j)%h_vec(k) = H_vec(k)
                        ctr(i,ny+j)%b_vec(k) = B_vec(k)
                    enddo
                    ctr(i,ny+j)%sh_vec = 0.0
                    ctr(i,ny+j)%sb_vec = 0.0
                enddo
            enddo
        endif
        
        !Deallocation
        deallocate(H_vec,B_vec)
        
    end subroutine bc_freestream_DIS
    
    
    !-------------------------------------------------------------------
    subroutine bc_user
    
    
    end subroutine bc_user
    
end module Boundary

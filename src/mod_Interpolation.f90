module mod_Interpolation
    use flow_data
    use Tools
    implicit none
    
    !--------------------------------------------------
    !>Reconstruct the slope of initial distribution function
    !>Index method
    !---------------------------------
    !           (i,j+1)              |
    !      ----------------          |
    !      |              |          |
    !      |              |          |
    !      |              |          |
    ! (i,j)|     (i,j)    |(i+1,j)   |
    !      |      Cell    |          |
    !      |              |          |
    !      |              |          |
    !      ----------------          |
    !            (i,j)               |
    !---------------------------------
    !--------------------------------------------------
    
    contains
    
    subroutine Interpolation()
        integer(KINT)                                   :: i,j
        !real(KREAL)                                     :: conVars(4),sw1(4),sw2(4),len1,len2

        !-------------------------------------------------------------------------------------------
        !>Calculate micro- and macro slopes: for distribution function and conservative variables
        !-------------------------------------------------------------------------------------------

        !-----------------------------------------
        !1) Inner part
        !-----------------------------------------
        !i direction
        do j=1-GHOST+1,ny+GHOST-1
            do i=1-GHOST+1,nx+GHOST-1
                call InterpInner(ctr(i-1,j),ctr(i,j),ctr(i+1,j),IDIRC)
            end do
        end do

        !j direction
        do j=1-GHOST+1,ny+GHOST-1
            do i=1-GHOST+1,nx+GHOST-1
                call InterpInner(ctr(i,j-1),ctr(i,j),ctr(i,j+1),JDIRC)
            end do
        end do
        
        !--------------------------------------------------------------------------------------------------
        !2) Boundary part: the ghost cell has not been set for solid wall and one-side interpolation is adopted
        !--------------------------------------------------------------------------------------------------    
        if(IP_X==0 .and. bc_type(1)==BCconst_wall_isothermal) then          !left boundary
            do j=1-GHOST,ny+GHOST
                call InterpBoundary(ctr(1,j),ctr(1,j),ctr(1+1,j),IDIRC)
            end do  
        endif
                
        if(IP_X==NP_X-1 .and. bc_type(2)==BCconst_wall_isothermal) then     !right boundary
            do j=1-GHOST,ny+GHOST
                call InterpBoundary(ctr(nx-1,j),ctr(nx,j),ctr(nx,j),IDIRC)
            end do
        endif
        
        if(IP_Y==0 .and. bc_type(3)==BCconst_wall_isothermal) then          !bottom boundary
            do i=1-GHOST,nx+GHOST
                call InterpBoundary(ctr(i,1),ctr(i,1),ctr(i,1+1),JDIRC)
            end do
        endif
        
        if(IP_Y==NP_Y-1 .and. bc_type(4)==BCconst_wall_isothermal) then     !upper boundary
            do i=1-GHOST,nx+GHOST
                call InterpBoundary(ctr(i,ny-1),ctr(i,ny),ctr(i,ny),JDIRC)
            end do
        endif 

    end subroutine Interpolation
    

    !--------------------------------------------------
    !>Interpolation of the boundary cell
    !>@param[inout] targetCell    :the target boundary cell
    !>@param[inout] leftCell      :the left cell
    !>@param[inout] rightCell     :the right cell
    !>@param[in]    idx           :the index indicating i or j direction
    !--------------------------------------------------
    subroutine InterpBoundary(leftCell,targetCell,rightCell,idx)
        type(CellCenter), intent(inout)                 :: targetCell
        type(CellCenter), intent(inout)                 :: leftCell,rightCell
        integer(KINT), intent(in)                       :: idx
        real(KREAL)                                     :: denom
        integer(KINT)                                   :: k        !Flattened velocity-space index

        !Micro velocity distribution function.
        !Use 1D velocity views to avoid creating Nu-by-Nv temporaries.
        denom = 0.5*(rightCell%length(idx)+leftCell%length(idx))
        do k=1,Nvel
            targetCell%sh_vec(k,idx) = (rightCell%h_vec(k)-leftCell%h_vec(k))/denom
            targetCell%sb_vec(k,idx) = (rightCell%b_vec(k)-leftCell%b_vec(k))/denom
        end do

        !Macro conservative variables
        targetCell%sw(:,idx)    = (rightCell%conVars-leftCell%conVars)/denom
    end subroutine InterpBoundary
    

    !--------------------------------------------------
    !>Interpolation of the inner cells
    !>@param[in]    leftCell      :the left cell
    !>@param[inout] targetCell    :the target cell
    !>@param[in]    rightCell     :the right cell
    !>@param[in]    idx           :the index indicating i or j direction
    !--------------------------------------------------
    subroutine InterpInner(leftCell,targetCell,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellCenter), intent(inout)                 :: targetCell
        integer(KINT), intent(in)                       :: idx

        call Interp_VanLeer(leftCell,targetCell,rightCell,idx)
    end subroutine InterpInner

    !--------------------------------------------------
    !>VanLeerLimiter for reconstruction of distrubution function
    !>@param[in]    leftCell  :the left cell
    !>@param[inout] midCell   :the middle cell
    !>@param[in]    rightCell :the right cell
    !>@param[in]    idx       :the index indicating i or j direction
    !--------------------------------------------------
    subroutine Interp_VanLeer(leftCell,targetCell,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellCenter), intent(inout)                 :: targetCell
        integer(KINT), intent(in)                       :: idx
        real(KREAL)                                     :: denomL,denomR
        real(KREAL)                                     :: sL,sR
        real(KREAL)                                     :: sLw(4),sRw(4)
        integer(KINT)                                   :: k        !Flattened velocity-space index

        !>Micro velocity distribution function.
        !>Use 1D velocity views for Van Leer slopes; the limiter formula is unchanged.
        denomL = 0.5*(targetCell%length(idx)+leftCell%length(idx))
        denomR = 0.5*(rightCell%length(idx)+targetCell%length(idx))

        do k=1,Nvel
            sL = (targetCell%h_vec(k)-leftCell%h_vec(k))/denomL
            sR = (rightCell%h_vec(k)-targetCell%h_vec(k))/denomR
            targetCell%sh_vec(k,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

            sL = (targetCell%b_vec(k)-leftCell%b_vec(k))/denomL
            sR = (rightCell%b_vec(k)-targetCell%b_vec(k))/denomR
            targetCell%sb_vec(k,idx) = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)
        end do
        
        !>Macro conservative variables
        sLw = (targetCell%conVars-leftCell%conVars)/denomL
        sRw = (rightCell%conVars-targetCell%conVars)/denomR
        targetCell%sw(:,idx) = (sign(UP,sRw)+sign(UP,sLw))*abs(sRw)*abs(sLw)/(abs(sRw)+abs(sLw)+SMV)
    end subroutine Interp_VanLeer
    
end module mod_Interpolation
    

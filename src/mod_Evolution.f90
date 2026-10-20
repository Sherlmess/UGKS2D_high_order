module mod_Evolution
    use flow_data
    use tools
    use mod_ImplicitCommon, only: ProjectToConservedMoments, MomentConserved
    use mod_Interpolation,only: InterpInner
    use mod_S2D1, only: IsS2CN2,IsS2Steady
    use mod_S2NativeAdapter, only: S2NativeEvolution
    implicit none
    integer(KREAL), parameter                           :: MNUM = 6 !Number of normal velocity moments
    integer(KREAL), parameter                           :: MTUM = 5 !Number of tangential velocity moments
    
contains

    ! Reuse the low-order body equilibrium preparation, not the S2 Q solver.
    ! This prepares a local target; it never changes an accepted state or correction.
    subroutine PrepareNativeFaceEquilibrium(h,b,vn,vt,prim)
        real(KREAL),intent(out)::h(:),b(:)
        real(KREAL),intent(in)::vn(:),vt(:),prim(4)
        real(KREAL)::reference_h(size(h)),reference_b(size(b)),target(4)
        logical::ok
        integer::abort_error
        call DiscreteMaxwell_LocalVec(h,b,vn,vt,prim)
        if(.not.UnifiedPhysicalOperators)return
        reference_h=h;reference_b=b;target=GetConserved(prim)
        call ProjectToConservedMoments(h,b,target,vn,vt,wVel,reference_h,reference_b,ok)
        if(.not.ok)then
            call MPI_Abort(MPI_COMM_WORLD,446,abort_error)
            error stop 'native face equilibrium preparation failed'
        endif
    end subroutine PrepareNativeFaceEquilibrium


    ! Keep the frozen native Evolution body in a branch-free procedure so
    ! selecting S2CN2 cannot perturb code generation for the default path.
    subroutine Evolution()
        if(IsS2CN2().or.IsS2Steady()) then
            call S2NativeEvolution()
        else
            call EvolutionNative()
        endif
    end subroutine Evolution

    !--------------------------------------------------
    !>Calculate the flux across the interfaces
    !--------------------------------------------------
    subroutine EvolutionNative()
        integer(KINT)                                   :: i,j,k,vel_i,vel_j,sym_idx
        integer(KINT)                                   :: flat_k     !Flattened velocity-space index for symmetry ghost VDFs
        integer(KINT)                                   :: sym_flat_k !Flattened mirrored velocity-space index for symmetry ghost VDFs

        !--------------------------------------------------
        !>Inner part
        !--------------------------------------------------
        !i direction
        if(TimeStep_Method==Time_local_varying .or. IsImplicitUnsteady()) then
            do j=1,ny
                do i=1,nx+1
                    call CalcFlux_MultiDimension_LocalTime(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
                end do
            end do
        else
            do j=1,ny
                do i=1,nx+1
                    call CalcFlux_MultiDimension(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
                end do
            end do
        endif
        
        !j direction
        if(TimeStep_Method==Time_local_varying .or. IsImplicitUnsteady()) then
            do j=1,ny+1
                do i=1,nx
                    call CalcFlux_MultiDimension_LocalTime(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
                end do
            end do
        else
            do j=1,ny+1
                do i=1,nx
                    call CalcFlux_MultiDimension(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
                end do
            end do
        endif
        
        !---------------------------------------------------------------------------
        !>Boundary part
        !>Wall boundary:  the interface flux is recomputed.
        !>Other boundary: the ghost cell has been adopted before and nothing is needed.
        !---------------------------------------------------------------------------
        if(IP_X==0 .and. bc_type(1)==BCconst_wall_isothermal) then    !I direction and lower boundaryface
            do j=1,ny
                call CalcFluxBoundary(WallPrimary(1),vface(1,j),ctr(1,j),IDIRC,Lside) !RN means no frame rotation
            enddo
        endif
            
        if(IP_X==NP_X-1 .and. bc_type(2)==BCconst_wall_isothermal) then   !I direction and upper boundaryface
            do j=1,ny
                call CalcFluxBoundary(WallPrimary(2),vface(nx+1,j),ctr(nx,j),IDIRC,Rside) !RY means frame rotation
            enddo
        endif
            
        if(IP_Y==0 .and. bc_type(3)==BCconst_wall_isothermal) then   !J direction and lower boundaryface
            do i=1,nx
                call CalcFluxBoundary(WallPrimary(3),hface(i,1),ctr(i,1),JDIRC,Lside) !RN means no frame rotation
            enddo   
        endif
                
        if(IP_Y==NP_Y-1 .and. bc_type(4)==BCconst_wall_isothermal) then   !J direction and upper boundaryface
            do i=1,nx
                call CalcFluxBoundary(WallPrimary(4),hface(i,ny+1),ctr(i,ny),JDIRC,Rside) !RY means frame rotation
            enddo    
        endif
        
        !For flate-plate geometry: the symmetry boundary is implemented in front of the plate tip. 
        if(BG_geometry==BGconst_flatplate .and. IP_Y==0) then
            !---------------flat-plate geometry---------------------
            !-----------------boundary setup-----------
            !            farfield (extrapolation)
            !          ---------------------------
            !         |        |                 |
            !   inflow|        |                 |extrapolation
            !         |        |                 |
            !          -------- ------------------
            !          symmetry         wall
            !-----------------------------------------
            !Here: the symmetry boundary at the bottom is implemented (x<0).
            
            if(ctr(1,1)%x>0) return     !the region locates at the wall surface
            
            do i=1,nx
                if(ctr(i,1)%x<0) then   !the cell is located in front of the tip: symmetry
                    !set ghost value for symmetry boundary
                    do j=1,GHOST
                        ctr(i,1-j)%conVars(1) = ctr(i,1+j-1)%conVars(1)
                        ctr(i,1-j)%conVars(2) = ctr(i,1+j-1)%conVars(2)
                        ctr(i,1-j)%conVars(3) = -ctr(i,1+j-1)%conVars(3)    !Only reverse v velocity
                        ctr(i,1-j)%conVars(4) = ctr(i,1+j-1)%conVars(4)
                        if (Quadrature_method == Quad_Unstruct) then
                            !Copy unstructured symmetry ghost VDFs through flattened cell-center views.
                            do k=1,Nvel
                                sym_idx = SymMap_V(k)
                                ctr(i,1-j)%h_vec(k) = ctr(i,1+j-1)%h_vec(sym_idx)
                                ctr(i,1-j)%b_vec(k) = ctr(i,1+j-1)%b_vec(sym_idx)
                            end do
                        else
                            !Mirror structured symmetry ghost VDFs while keeping the storage loop flat.
                            do flat_k=1,Nvel
                                vel_i = mod(flat_k-1,Nu)+1
                                vel_j = (flat_k-1)/Nu+1
                                sym_flat_k = vel_i+(Nv-vel_j)*Nu
                                ctr(i,1-j)%h_vec(flat_k) = ctr(i,1+j-1)%h_vec(sym_flat_k)
                                ctr(i,1-j)%b_vec(flat_k) = ctr(i,1+j-1)%b_vec(sym_flat_k)
                            end do
                        end if
                        !Do the same thing for slope
                        ctr(i,1-j)%sh_vec = 0.0
                        ctr(i,1-j)%sb_vec = 0.0
                    end do
                    
                    !interpolation
                    do j=1-GHOST+1,1
                        call InterpInner(ctr(i,j-1),ctr(i,j),ctr(i,j+1),JDIRC)
                    end do
                    
                    !evolution: calculate flux
                    if(TimeStep_Method==Time_local_varying .or. IsImplicitUnsteady()) then
                        call CalcFlux_MultiDimension_LocalTime(ctr(i,1-1),hface(i,1),ctr(i,1),JDIRC)
                    else
                        call CalcFlux_MultiDimension(ctr(i,1-1),hface(i,1),ctr(i,1),JDIRC)
                    endif
                    
                elseif(ctr(i,1)%x>0) then
                    exit
                endif
            enddo
        endif

    end subroutine EvolutionNative

    !-------------------------------------------------------------------------
    !>Calculate flux of inner interface with multi-dimension (for global time step)
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !>@param[in]    idx       :index indicating i or j direction
    !-------------------------------------------------------------------------
    subroutine CalcFlux_MultiDimension(leftCell,face,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellInterface), intent(inout)              :: face
        integer(KINT), intent(in)                       :: idx !indicator for direction and boundary
        integer(KINT)                                   :: k
        real(KREAL), allocatable, dimension(:)          :: vn_vec,vt_vec !Pure 1D local-frame velocities
        real(KREAL), allocatable, dimension(:)          :: h_vec,b_vec !Pure 1D distribution function at the interface
        real(KREAL), allocatable, dimension(:)          :: H0_vec,B0_vec !Pure 1D Maxwellian distribution function
        real(KREAL), allocatable, dimension(:)          :: H_plus_vec,B_plus_vec !Pure 1D Shakhov correction
        real(KREAL), allocatable, dimension(:)          :: shn_vec,sbn_vec !Pure 1D normal distribution slopes
        real(KREAL), allocatable, dimension(:)          :: sht_vec,sbt_vec !Pure 1D tangential distribution slopes
        real(KREAL), allocatable, dimension(:)          :: delta_vec !Pure 1D Heaviside step function
        real(KREAL)                                     :: deltaL,deltaR
        real(KREAL)                                     :: wk,vnK,vtK,v2
        real(KREAL)                                     :: slope_h,slope_b !Velocity-space directional slopes used in flux integrals
        real(KREAL)                                     :: flux_plus(4) !Accumulated conservative flux from Shakhov correction
        real(KREAL)                                     :: flux_free(4) !Accumulated conservative flux from reconstructed f0
        real(KREAL)                                     :: flux_slope(4) !Accumulated conservative flux from spatial slopes
        real(KREAL)                                     :: h_aL,h_aR,h_b,h_aT !H-flux Maxwellian slope terms: left, right, tangential, time
        real(KREAL)                                     :: b_aL,b_aR,b_b,b_aT !B-flux Maxwellian slope terms: left, right, tangential, time
        real(KREAL)                                     :: prim(4),conVars(4) !Primary and conservative variables at the interface
        real(KREAL)                                     :: qf(2) !Heat flux in normal and tangential direction
        real(KREAL)                                     :: sw(4),swt(4) !Slope of conVars at normal and tangential direction
        real(KREAL)                                     :: aL(4),aR(4),b_slope(4),aT(4) !Micro slope of Maxwellian distribution at normal tangential and time.
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
        real(KREAL)                                     :: Mau0(4),MauL(4),MauR(4),Mbv(4),MauT(4) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<b*u*v*\psi>,<A*u*\psi>
        real(KREAL)                                     :: tau !Collision time
        real(KREAL)                                     :: Mt(5) !Some time integration terms
        real(KREAL)                                     :: dt!global or local timestep
    
        !--------------------------------------------------
        !Prepare
        !--------------------------------------------------
        !Allocate pure 1D velocity-space temporaries for the interface flux.
        allocate(vn_vec(Nvel))
        allocate(vt_vec(Nvel))
        allocate(delta_vec(Nvel))
        allocate(h_vec(Nvel))
        allocate(b_vec(Nvel))
        allocate(shn_vec(Nvel))
        allocate(sbn_vec(Nvel))
        allocate(sht_vec(Nvel))
        allocate(sbt_vec(Nvel))
        allocate(H0_vec(Nvel))
        allocate(B0_vec(Nvel))
        allocate(H_plus_vec(Nvel))
        allocate(B_plus_vec(Nvel))
        
        !Time step
        dt = dt_global
    
        !Convert the velocity space to local frame and set the Heaviside function.
        do k=1,Nvel
            vn_vec(k) = uVel(k)*face%cosx+vVel(k)*face%cosy
            vt_vec(k) =-uVel(k)*face%cosy+vVel(k)*face%cosx
            delta_vec(k) = (sign(UP,vn_vec(k))+1.0_KREAL)/2.0_KREAL
        enddo
    
        !--------------------------------------------------
        !Reconstruct initial distribution at interface
        !--------------------------------------------------
        do k=1,Nvel
            deltaL = delta_vec(k)
            deltaR = 1.0_KREAL-deltaL

            h_vec(k) = (leftCell%h_vec(k)+0.5_KREAL*leftCell%length(idx)*leftCell%sh_vec(k,idx))*deltaL+&
                       (rightCell%h_vec(k)-0.5_KREAL*rightCell%length(idx)*rightCell%sh_vec(k,idx))*deltaR
            b_vec(k) = (leftCell%b_vec(k)+0.5_KREAL*leftCell%length(idx)*leftCell%sb_vec(k,idx))*deltaL+&
                       (rightCell%b_vec(k)-0.5_KREAL*rightCell%length(idx)*rightCell%sb_vec(k,idx))*deltaR
            shn_vec(k) = leftCell%sh_vec(k,idx)*deltaL+rightCell%sh_vec(k,idx)*deltaR
            sbn_vec(k) = leftCell%sb_vec(k,idx)*deltaL+rightCell%sb_vec(k,idx)*deltaR

            ! Tangential part (notice the negative sign for J-direction rotation).
            if (idx==1) then
                sht_vec(k) = leftCell%sh_vec(k,2)*deltaL+rightCell%sh_vec(k,2)*deltaR
                sbt_vec(k) = leftCell%sb_vec(k,2)*deltaL+rightCell%sb_vec(k,2)*deltaR
            else
                sht_vec(k) = -(leftCell%sh_vec(k,1)*deltaL+rightCell%sh_vec(k,1)*deltaR)
                sbt_vec(k) = -(leftCell%sb_vec(k,1)*deltaL+rightCell%sb_vec(k,1)*deltaR)
            endif
        enddo
        
        !--------------------------------------------------
        !Obtain macroscopic variables
        !--------------------------------------------------
        !Conservative variables conVars at interface
        conVars = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            conVars(1) = conVars(1)+wk*h_vec(k)
            conVars(2) = conVars(2)+wk*vnK*h_vec(k)
            conVars(3) = conVars(3)+wk*vtK*h_vec(k)
            conVars(4) = conVars(4)+0.5_KREAL*wk*(v2*h_vec(k)+b_vec(k))
        enddo
        
        !Convert to primary variables
        prim = GetPrimary(conVars)
    
        !--------------------------------------------------
        !Calculate slope aL, aR
        !--------------------------------------------------
        sw  = (conVars-LocalFrame(leftCell%conVars,face%cosx,face%cosy))/(0.5*leftCell%length(idx)) !left slope of conVars
        aL  = MicroSlope(prim,sw) !Calculate a^L

        sw  = (LocalFrame(rightCell%conVars,face%cosx,face%cosy)-conVars)/(0.5*rightCell%length(idx)) !right slope of conVars
        aR  = MicroSlope(prim,sw) !Calculate a^R
        
        !--------------------------------------------------
        !Calculate slope b
        !--------------------------------------------------
        !b_slope is calculated through integrating sb
        swt = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            swt(1) = swt(1)+wk*sht_vec(k)
            swt(2) = swt(2)+wk*vnK*sht_vec(k)
            swt(3) = swt(3)+wk*vtK*sht_vec(k)
            swt(4) = swt(4)+0.5_KREAL*wk*(v2*sht_vec(k)+sbt_vec(k))
        enddo
        b_slope = MicroSlope(prim,swt)      !The negative sign has been considered in sht and sbt 
    
        !b_slope is calculated through interpolation of macroflow (notice the negative sign)
        !if(idx==Idirc) then     ! Idirc: vface
        !    swt = LocalFrame(0.5*(leftCell%sw(:,2-idx+1)+rightCell%sw(:,2-idx+1)),face%cosx,face%cosy)
        !elseif(idx==Jdirc) then ! Jdirc: hface, rotation for derivative
        !    swt = LocalFrame(-0.5*(leftCell%sw(:,2-idx+1)+rightCell%sw(:,2-idx+1)),face%cosx,face%cosy)
        !endif
        !b_slope = MicroSlope(prim,swt)
        
        !--------------------------------------------------
        !Calculate time slope of conVars and A
        !--------------------------------------------------
        !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR) 
        
        MauL = Moment_auvxi(aL,MuL,Mv,Mxi,1,0) !<aL*u*\psi>_{>0}
        MauR = Moment_auvxi(aR,MuR,Mv,Mxi,1,0) !<aR*u*\psi>_{<0}
        Mbv  = Moment_auvxi(b_slope,Mu,Mv,Mxi,0,1) !<b*v*\psi>
    
        sw = -prim(1)*(MauL+MauR+Mbv) !Time slope of conVars
        aT = MicroSlope(prim,sw) !Calculate A
    
        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim)
        if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau')then
            ComparisonFaceTauMin=min(ComparisonFaceTauMin,tau)
            ComparisonFaceTauMax=max(ComparisonFaceTauMax,tau)
            ComparisonFaceTauCalls=ComparisonFaceTauCalls+1
        endif
    
        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)
    
        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g0
        !--------------------------------------------------
        Mau0 = Moment_uvxi(Mu,Mv,Mxi,1,0,0) !<u*\psi>
        MauL = Moment_auvxi(aL,MuL,Mv,Mxi,2,0) !<aL*u^2*\psi>_{>0}
        MauR = Moment_auvxi(aR,MuR,Mv,Mxi,2,0) !<aR*u^2*\psi>_{<0}
        Mbv  = Moment_auvxi(b_slope,Mu,Mv,Mxi,1,1) !<b*u*v*\psi>
        MauT = Moment_auvxi(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>
    
        face%L_flux = Mt(1)*prim(1)*Mau0+Mt(2)*prim(1)*(MauL+MauR+Mbv)+Mt(3)*prim(1)*MauT
    
        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g+ and f0
        !--------------------------------------------------
        !Maxwellian distribution H0 and B0
        call PrepareNativeFaceEquilibrium(H0_vec,B0_vec,vn_vec,vt_vec,prim)
    
        !Calculate heat flux
        qf = GetHeatFlux_LocalVec(h_vec,b_vec,vn_vec,vt_vec,prim)
    
        !Shakhov part H+ and B+
        call ShakhovPart_LocalVec(H0_vec,B0_vec,vn_vec,vt_vec,qf,prim,H_plus_vec,B_plus_vec)
    
        !Conservative flux related to g+ and f0.
        !Use one Nvel loop for the three original sum groups:
        !flux_plus -> Mt(1)*g+, flux_free -> Mt(4)*f0, flux_slope -> -Mt(5)*grad(f0).
        flux_plus  = 0.0_KREAL
        flux_free  = 0.0_KREAL
        flux_slope = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            slope_h = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_b = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            flux_plus(1) = flux_plus(1)+wk*vnK*H_plus_vec(k)
            flux_plus(2) = flux_plus(2)+wk*vnK*vnK*H_plus_vec(k)
            flux_plus(3) = flux_plus(3)+wk*vtK*vnK*H_plus_vec(k)
            flux_plus(4) = flux_plus(4)+0.5_KREAL*wk*vnK*(v2*H_plus_vec(k)+B_plus_vec(k))

            flux_free(1) = flux_free(1)+wk*vnK*h_vec(k)
            flux_free(2) = flux_free(2)+wk*vnK*vnK*h_vec(k)
            flux_free(3) = flux_free(3)+wk*vtK*vnK*h_vec(k)
            flux_free(4) = flux_free(4)+0.5_KREAL*wk*vnK*(v2*h_vec(k)+b_vec(k))

            flux_slope(1) = flux_slope(1)+wk*vnK*slope_h
            flux_slope(2) = flux_slope(2)+wk*vnK*vnK*slope_h
            flux_slope(3) = flux_slope(3)+wk*vtK*vnK*slope_h
            flux_slope(4) = flux_slope(4)+0.5_KREAL*wk*vnK*(v2*slope_h+slope_b)
        enddo
        face%L_flux = face%L_flux+Mt(1)*flux_plus+Mt(4)*flux_free-Mt(5)*flux_slope
    
        !--------------------------------------------------
        !Calculate flux of distribution function.
        !Evaluate each velocity point through 1D views to avoid large Nu-by-Nv temporaries.
        do k=1,Nvel
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaL  = delta_vec(k)
            deltaR  = 1.0_KREAL-deltaL
            slope_h = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_b = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            h_aL = aL(1)*H0_vec(k)+aL(2)*vnK*H0_vec(k)+aL(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aL(4)*(v2*H0_vec(k)+B0_vec(k))
            h_aR = aR(1)*H0_vec(k)+aR(2)*vnK*H0_vec(k)+aR(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aR(4)*(v2*H0_vec(k)+B0_vec(k))
            h_b  = b_slope(1)*H0_vec(k)+b_slope(2)*vnK*H0_vec(k)+b_slope(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*b_slope(4)*(v2*H0_vec(k)+B0_vec(k))
            h_aT = aT(1)*H0_vec(k)+aT(2)*vnK*H0_vec(k)+aT(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aT(4)*(v2*H0_vec(k)+B0_vec(k))

            b_aL = aL(1)*B0_vec(k)+aL(2)*vnK*B0_vec(k)+aL(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aL(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_aR = aR(1)*B0_vec(k)+aR(2)*vnK*B0_vec(k)+aR(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aR(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_b  = b_slope(1)*B0_vec(k)+b_slope(2)*vnK*B0_vec(k)+b_slope(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*b_slope(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_aT = aT(1)*B0_vec(k)+aT(2)*vnK*B0_vec(k)+aT(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aT(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))

            face%L_flux_h_vec(k) = Mt(1)*vnK*(H0_vec(k)+H_plus_vec(k))+&
                                   Mt(2)*vnK*vnK*(h_aL*deltaL+h_aR*deltaR)+&
                                   Mt(2)*vnK*vtK*h_b+&
                                   Mt(3)*vnK*h_aT+&
                                   Mt(4)*vnK*h_vec(k)-Mt(5)*vnK*slope_h

            face%L_flux_b_vec(k) = Mt(1)*vnK*(B0_vec(k)+B_plus_vec(k))+&
                                   Mt(2)*vnK*vnK*(b_aL*deltaL+b_aR*deltaR)+&
                                   Mt(2)*vnK*vtK*b_b+&
                                   Mt(3)*vnK*b_aT+&
                                   Mt(4)*vnK*b_vec(k)-Mt(5)*vnK*slope_b
        enddo
        
        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        face%L_flux = GlobalFrame(face%L_flux,face%cosx,face%cosy) 
        
        !Total flux
        face%L_flux = face%length*face%L_flux

        !Scale micro fluxes and copy L->R through 1D views.
        do k=1,Nvel
            face%L_flux_h_vec(k) = face%length*face%L_flux_h_vec(k)
            face%L_flux_b_vec(k) = face%length*face%L_flux_b_vec(k)
            face%R_flux_h_vec(k) = face%L_flux_h_vec(k)
            face%R_flux_b_vec(k) = face%L_flux_b_vec(k)
        enddo

        !R_flux == L_flux
        face%R_flux     = face%L_flux
    
        ! Match the native body discrete moment convention. No flux/state projection.
        if(UnifiedPhysicalOperators)then
            call MomentConserved(face%L_flux_h_vec,face%L_flux_b_vec,uVel,vVel,wVel,face%L_flux)
            call MomentConserved(face%R_flux_h_vec,face%R_flux_b_vec,uVel,vVel,wVel,face%R_flux)
        endif

        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Deallocate array
        deallocate(vn_vec)
        deallocate(vt_vec)
        deallocate(delta_vec)
        deallocate(h_vec)
        deallocate(b_vec)
        deallocate(shn_vec)
        deallocate(sbn_vec)
        deallocate(sht_vec)
        deallocate(sbt_vec)
        deallocate(H0_vec)
        deallocate(B0_vec)
        deallocate(H_plus_vec)
        deallocate(B_plus_vec)
    end subroutine CalcFlux_MultiDimension
    
    
    !-----------------------------------------------------------------------
    !>Calculate flux of inner interface with multi-dimension
    !>@param[in]    leftCell  :cell left to the target interface
    !>@param[inout] face      :the target interface
    !>@param[in]    rightCell :cell right to the target interface
    !>@param[in]    idx       :index indicating i or j direction
    !-----------------------------------------------------------------------
    subroutine CalcFlux_MultiDimension_LocalTime(leftCell,face,rightCell,idx)
        type(CellCenter), intent(in)                    :: leftCell,rightCell
        type(CellInterface), intent(inout)              :: face
        integer(KINT), intent(in)                       :: idx !indicator for direction and boundary
        integer(KINT)                                   :: k !Flattened velocity-space index
        real(KREAL), allocatable, dimension(:)          :: vn_vec,vt_vec !Pure 1D local-frame velocities
        real(KREAL), allocatable, dimension(:)          :: h_vec,b_vec !Pure 1D distribution function at the interface
        real(KREAL), allocatable, dimension(:)          :: H0_vec,B0_vec !Pure 1D Maxwellian distribution function
        real(KREAL), allocatable, dimension(:)          :: H_plus_vec,B_plus_vec !Pure 1D Shakhov correction
        real(KREAL), allocatable, dimension(:)          :: shn_vec,sbn_vec !Pure 1D normal distribution slopes
        real(KREAL), allocatable, dimension(:)          :: sht_vec,sbt_vec !Pure 1D tangential distribution slopes
        real(KREAL), allocatable, dimension(:)          :: delta_vec !Pure 1D Heaviside step function
        real(KREAL)                                     :: flux(4,5) !row: flux(1:4,:)->flux for conVars; column: flux(:,1:5)->Mt(1:5)
        real(KREAL)                                     :: deltaL,deltaR !Left/right upwind weights at one velocity point
        real(KREAL)                                     :: wk,vnK,vtK,v2 !Reusable velocity-space scalars for Nvel loops
        real(KREAL)                                     :: slope_h,slope_b !Velocity-space directional slopes used in flux integrals
        real(KREAL)                                     :: flux_plus(4) !Accumulated conservative flux from Shakhov correction
        real(KREAL)                                     :: flux_free(4) !Accumulated conservative flux from reconstructed f0
        real(KREAL)                                     :: flux_slope(4) !Accumulated conservative flux from spatial slopes
        real(KREAL)                                     :: h_aL,h_aR,h_b,h_aT !H-flux Maxwellian slope terms: left, right, tangential, time
        real(KREAL)                                     :: b_aL,b_aR,b_b,b_aT !B-flux Maxwellian slope terms: left, right, tangential, time
        real(KREAL)                                     :: prim(4),conVars(4) !Primary and conservative variables at the interface
        real(KREAL)                                     :: qf(2) !Heat flux in normal and tangential direction
        real(KREAL)                                     :: sw(4),swt(4) !Slope of conVars at normal and tangential direction
        real(KREAL)                                     :: aL(4),aR(4),b_slope(4),aT(4) !Micro slope of Maxwellian distribution at normal tangential and time.
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM),Mv(0:MTUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<v^m>,<\xi^l>
        real(KREAL)                                     :: Mau0(4),MauL(4),MauR(4),Mbv(4),MauT(4) !<u\psi>,<aL*u^n*\psi>,<aR*u^n*\psi>,<b*u*v*\psi>,<A*u*\psi>
        real(KREAL)                                     :: tau !Collision time
        real(KREAL)                                     :: Mt(5) !Some time integration terms
        real(KREAL)                                     :: dt!global or local timestep
    
        !--------------------------------------------------
        !Prepare
        !--------------------------------------------------
        !Allocate pure 1D velocity-space temporaries for the local-time interface flux.
        allocate(vn_vec(Nvel))
        allocate(vt_vec(Nvel))
        allocate(delta_vec(Nvel))
        allocate(h_vec(Nvel))
        allocate(b_vec(Nvel))
        allocate(shn_vec(Nvel))
        allocate(sbn_vec(Nvel))
        allocate(sht_vec(Nvel))
        allocate(sbt_vec(Nvel))
        allocate(H0_vec(Nvel))
        allocate(B0_vec(Nvel))
        allocate(H_plus_vec(Nvel))
        allocate(B_plus_vec(Nvel))

        !Convert the velocity space to local frame and set the Heaviside function.
        do k=1,Nvel
            vn_vec(k) = uVel(k)*face%cosx+vVel(k)*face%cosy
            vt_vec(k) =-uVel(k)*face%cosy+vVel(k)*face%cosx
            delta_vec(k) = (sign(UP,vn_vec(k))+1.0_KREAL)/2.0_KREAL
        enddo
    
        !--------------------------------------------------
        !Reconstruct initial distribution at interface
        !--------------------------------------------------
        !Reconstruct interface distributions through 1D velocity views.
        do k=1,Nvel
            deltaL = delta_vec(k)
            deltaR = 1.0_KREAL-deltaL

            h_vec(k) = (leftCell%h_vec(k)+0.5_KREAL*leftCell%length(idx)*leftCell%sh_vec(k,idx))*deltaL+&
                       (rightCell%h_vec(k)-0.5_KREAL*rightCell%length(idx)*rightCell%sh_vec(k,idx))*deltaR
            b_vec(k) = (leftCell%b_vec(k)+0.5_KREAL*leftCell%length(idx)*leftCell%sb_vec(k,idx))*deltaL+&
                       (rightCell%b_vec(k)-0.5_KREAL*rightCell%length(idx)*rightCell%sb_vec(k,idx))*deltaR
            shn_vec(k) = leftCell%sh_vec(k,idx)*deltaL+rightCell%sh_vec(k,idx)*deltaR
            sbn_vec(k) = leftCell%sb_vec(k,idx)*deltaL+rightCell%sb_vec(k,idx)*deltaR

            !Tangential part (notice the negative sign for J-direction rotation).
            if (idx==1) then
                sht_vec(k) = leftCell%sh_vec(k,2)*deltaL+rightCell%sh_vec(k,2)*deltaR
                sbt_vec(k) = leftCell%sb_vec(k,2)*deltaL+rightCell%sb_vec(k,2)*deltaR
            else
                sht_vec(k) = -(leftCell%sh_vec(k,1)*deltaL+rightCell%sh_vec(k,1)*deltaR)
                sbt_vec(k) = -(leftCell%sb_vec(k,1)*deltaL+rightCell%sb_vec(k,1)*deltaR)
            endif
        enddo
        
        !--------------------------------------------------
        !Obtain macroscopic variables
        !--------------------------------------------------
        !Conservative variables conVars at interface.
        conVars = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            conVars(1) = conVars(1)+wk*h_vec(k)
            conVars(2) = conVars(2)+wk*vnK*h_vec(k)
            conVars(3) = conVars(3)+wk*vtK*h_vec(k)
            conVars(4) = conVars(4)+0.5_KREAL*wk*(v2*h_vec(k)+b_vec(k))
        enddo
        
        !Convert to primary variables
        prim = GetPrimary(conVars)
    
        !--------------------------------------------------
        !Calculate slope aL, aR
        !--------------------------------------------------
        sw  = (conVars-LocalFrame(leftCell%conVars,face%cosx,face%cosy))/(0.5*leftCell%length(idx)) !left slope of conVars
        aL  = MicroSlope(prim,sw) !Calculate a^L

        sw  = (LocalFrame(rightCell%conVars,face%cosx,face%cosy)-conVars)/(0.5*rightCell%length(idx)) !right slope of conVars
        aR  = MicroSlope(prim,sw) !Calculate a^R
        
        !--------------------------------------------------
        !Calculate slope b
        !--------------------------------------------------
        !b_slope is calculated through integrating sb
        !Tangential macro slope obtained by integrating the reconstructed VDF slopes.
        swt = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            swt(1) = swt(1)+wk*sht_vec(k)
            swt(2) = swt(2)+wk*vnK*sht_vec(k)
            swt(3) = swt(3)+wk*vtK*sht_vec(k)
            swt(4) = swt(4)+0.5_KREAL*wk*(v2*sht_vec(k)+sbt_vec(k))
        enddo
        b_slope = MicroSlope(prim,swt)
    
        !b_slope is calculated through interpolation of macroflow (notice the negative sign)
        !if(idx==Idirc) then     ! Idirc: vface
        !    swt = LocalFrame(0.5*(leftCell%sw(:,2-idx+1)+rightCell%sw(:,2-idx+1)),face%cosx,face%cosy)
        !elseif(idx==Jdirc) then ! Jdirc: hface, rotation for derivative
        !    swt = LocalFrame(-0.5*(leftCell%sw(:,2-idx+1)+rightCell%sw(:,2-idx+1)),face%cosx,face%cosy)
        !endif
        !b_slope = MicroSlope(prim,swt)
        
        !--------------------------------------------------
        !Calculate time slope of conVars and A
        !--------------------------------------------------
        !<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
        call CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR) 
        
        MauL = Moment_auvxi(aL,MuL,Mv,Mxi,1,0) !<aL*u*\psi>_{>0}
        MauR = Moment_auvxi(aR,MuR,Mv,Mxi,1,0) !<aR*u*\psi>_{<0}
        Mbv  = Moment_auvxi(b_slope,Mu,Mv,Mxi,0,1) !<b*u*v*\psi>
    
        sw = -prim(1)*(MauL+MauR+Mbv) !Time slope of conVars
        aT = MicroSlope(prim,sw) !Calculate A
    
        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g0
        !--------------------------------------------------
        Mau0 = Moment_uvxi(Mu,Mv,Mxi,1,0,0) !<u*\psi>
        MauL = Moment_auvxi(aL,MuL,Mv,Mxi,2,0) !<aL*u^2*\psi>_{>0}
        MauR = Moment_auvxi(aR,MuR,Mv,Mxi,2,0) !<aR*u^2*\psi>_{<0}
        Mbv  = Moment_auvxi(b_slope,Mu,Mv,Mxi,1,1) !<b*u*v*\psi>
        MauT = Moment_auvxi(aT,Mu,Mv,Mxi,1,0) !<A*u*\psi>
    
        !flux = Mt(1)*prim(1)*Mau0+Mt(2)*prim(1)*(MauL+MauR+Mbv)+Mt(3)*prim(1)*MauT
        flux(:,1) = prim(1)*Mau0
        flux(:,2) = prim(1)*(MauL+MauR+Mbv)
        flux(:,3) = prim(1)*MauT
    
        !--------------------------------------------------
        !Calculate the flux of conservative variables related to g+ and f0
        !--------------------------------------------------
        !Maxwellian distribution H0 and B0
        call DiscreteMaxwell_LocalVec(H0_vec,B0_vec,vn_vec,vt_vec,prim)
    
        !Calculate heat flux
        qf = GetHeatFlux_LocalVec(h_vec,b_vec,vn_vec,vt_vec,prim)
    
        !Shakhov part H+ and B+
        call ShakhovPart_LocalVec(H0_vec,B0_vec,vn_vec,vt_vec,qf,prim,H_plus_vec,B_plus_vec)
    
        !Conservative flux related to g+ and f0.
        !The local-time path keeps the original five-column flux decomposition:
        !column 1 receives g+, column 4 receives f0, and column 5 receives -grad(f0).
        flux_plus  = 0.0_KREAL
        flux_free  = 0.0_KREAL
        flux_slope = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            slope_h = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_b = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            flux_plus(1) = flux_plus(1)+wk*vnK*H_plus_vec(k)
            flux_plus(2) = flux_plus(2)+wk*vnK*vnK*H_plus_vec(k)
            flux_plus(3) = flux_plus(3)+wk*vtK*vnK*H_plus_vec(k)
            flux_plus(4) = flux_plus(4)+0.5_KREAL*wk*vnK*(v2*H_plus_vec(k)+B_plus_vec(k))

            flux_free(1) = flux_free(1)+wk*vnK*h_vec(k)
            flux_free(2) = flux_free(2)+wk*vnK*vnK*h_vec(k)
            flux_free(3) = flux_free(3)+wk*vtK*vnK*h_vec(k)
            flux_free(4) = flux_free(4)+0.5_KREAL*wk*vnK*(v2*h_vec(k)+b_vec(k))

            flux_slope(1) = flux_slope(1)+wk*vnK*slope_h
            flux_slope(2) = flux_slope(2)+wk*vnK*vnK*slope_h
            flux_slope(3) = flux_slope(3)+wk*vtK*vnK*slope_h
            flux_slope(4) = flux_slope(4)+0.5_KREAL*wk*vnK*(v2*slope_h+slope_b)
        enddo
        flux(:,1) = flux(:,1)+flux_plus
        flux(:,4) = flux_free
        flux(:,5) = -flux_slope
    
        !Micro distribution fluxes are evaluated directly in the L/R local-time loops below.
        !This removes the former flux_h(Nu,Nv,5)/flux_b(Nu,Nv,5) temporary arrays.
        
        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !collision time
        tau = GetTau(prim)

        !--------------------------------------------------
        !>L_flux of right cell
        !--------------------------------------------------
        !Time step
        dt = rightCell%dt_local
        
        !Time integration terms
        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)
        
        !Interface flux
        face%L_flux     = (Mt(1)*flux(:,1)+Mt(2)*flux(:,2)+Mt(3)*flux(:,3)+Mt(4)*flux(:,4)+Mt(5)*flux(:,5))*face%length

        !Right-cell local time step: combine micro flux components into L_flux through 1D views.
        do k=1,Nvel
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaL  = delta_vec(k)
            deltaR  = 1.0_KREAL-deltaL
            slope_h = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_b = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            h_aL = aL(1)*H0_vec(k)+aL(2)*vnK*H0_vec(k)+aL(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aL(4)*(v2*H0_vec(k)+B0_vec(k))
            h_aR = aR(1)*H0_vec(k)+aR(2)*vnK*H0_vec(k)+aR(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aR(4)*(v2*H0_vec(k)+B0_vec(k))
            h_b  = b_slope(1)*H0_vec(k)+b_slope(2)*vnK*H0_vec(k)+b_slope(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*b_slope(4)*(v2*H0_vec(k)+B0_vec(k))
            h_aT = aT(1)*H0_vec(k)+aT(2)*vnK*H0_vec(k)+aT(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aT(4)*(v2*H0_vec(k)+B0_vec(k))

            b_aL = aL(1)*B0_vec(k)+aL(2)*vnK*B0_vec(k)+aL(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aL(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_aR = aR(1)*B0_vec(k)+aR(2)*vnK*B0_vec(k)+aR(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aR(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_b  = b_slope(1)*B0_vec(k)+b_slope(2)*vnK*B0_vec(k)+b_slope(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*b_slope(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_aT = aT(1)*B0_vec(k)+aT(2)*vnK*B0_vec(k)+aT(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aT(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))

            face%L_flux_h_vec(k) = (Mt(1)*vnK*(H0_vec(k)+H_plus_vec(k))+&
                                   Mt(2)*(vnK*vnK*(h_aL*deltaL+h_aR*deltaR)+vnK*vtK*h_b)+&
                                   Mt(3)*vnK*h_aT+Mt(4)*vnK*h_vec(k)-Mt(5)*vnK*slope_h)*face%length

            face%L_flux_b_vec(k) = (Mt(1)*vnK*(B0_vec(k)+B_plus_vec(k))+&
                                   Mt(2)*(vnK*vnK*(b_aL*deltaL+b_aR*deltaR)+vnK*vtK*b_b)+&
                                   Mt(3)*vnK*b_aT+Mt(4)*vnK*b_vec(k)-Mt(5)*vnK*slope_b)*face%length
        enddo
        
        !--------------------------------------------------
        !>R_flux of left cell
        !--------------------------------------------------
        !Time step
        dt = leftCell%dt_local
        
        !Time integration terms
        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)
        
        !Interface flux
        face%R_flux     = (Mt(1)*flux(:,1)+Mt(2)*flux(:,2)+Mt(3)*flux(:,3)+Mt(4)*flux(:,4)+Mt(5)*flux(:,5))*face%length

        !Left-cell local time step: combine the same micro components into R_flux through 1D views.
        do k=1,Nvel
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaL  = delta_vec(k)
            deltaR  = 1.0_KREAL-deltaL
            slope_h = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_b = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            h_aL = aL(1)*H0_vec(k)+aL(2)*vnK*H0_vec(k)+aL(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aL(4)*(v2*H0_vec(k)+B0_vec(k))
            h_aR = aR(1)*H0_vec(k)+aR(2)*vnK*H0_vec(k)+aR(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aR(4)*(v2*H0_vec(k)+B0_vec(k))
            h_b  = b_slope(1)*H0_vec(k)+b_slope(2)*vnK*H0_vec(k)+b_slope(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*b_slope(4)*(v2*H0_vec(k)+B0_vec(k))
            h_aT = aT(1)*H0_vec(k)+aT(2)*vnK*H0_vec(k)+aT(3)*vtK*H0_vec(k)+&
                   0.5_KREAL*aT(4)*(v2*H0_vec(k)+B0_vec(k))

            b_aL = aL(1)*B0_vec(k)+aL(2)*vnK*B0_vec(k)+aL(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aL(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_aR = aR(1)*B0_vec(k)+aR(2)*vnK*B0_vec(k)+aR(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aR(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_b  = b_slope(1)*B0_vec(k)+b_slope(2)*vnK*B0_vec(k)+b_slope(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*b_slope(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))
            b_aT = aT(1)*B0_vec(k)+aT(2)*vnK*B0_vec(k)+aT(3)*vtK*B0_vec(k)+&
                   0.5_KREAL*aT(4)*(v2*B0_vec(k)+Mxi(2)*H0_vec(k))

            face%R_flux_h_vec(k) = (Mt(1)*vnK*(H0_vec(k)+H_plus_vec(k))+&
                                   Mt(2)*(vnK*vnK*(h_aL*deltaL+h_aR*deltaR)+vnK*vtK*h_b)+&
                                   Mt(3)*vnK*h_aT+Mt(4)*vnK*h_vec(k)-Mt(5)*vnK*slope_h)*face%length

            face%R_flux_b_vec(k) = (Mt(1)*vnK*(B0_vec(k)+B_plus_vec(k))+&
                                   Mt(2)*(vnK*vnK*(b_aL*deltaL+b_aR*deltaR)+vnK*vtK*b_b)+&
                                   Mt(3)*vnK*b_aT+Mt(4)*vnK*b_vec(k)-Mt(5)*vnK*slope_b)*face%length
        enddo
        
        !--------------------------------------------------
        !>Convert to global frame
        !--------------------------------------------------
        face%L_flux     = GlobalFrame(face%L_flux,face%cosx,face%cosy) 
        face%R_flux     = GlobalFrame(face%R_flux,face%cosx,face%cosy)
        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Release pure 1D velocity-space temporaries allocated above.
        deallocate(vn_vec)
        deallocate(vt_vec)
        deallocate(delta_vec)
        deallocate(h_vec)
        deallocate(b_vec)
        deallocate(shn_vec)
        deallocate(sbn_vec)
        deallocate(sht_vec)
        deallocate(sbt_vec)
        deallocate(H0_vec)
        deallocate(B0_vec)
        deallocate(H_plus_vec)
        deallocate(B_plus_vec)
    end subroutine CalcFlux_MultiDimension_LocalTime
    

    !--------------------------------------------------
    !>Calculate flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------
    subroutine CalcFluxBoundary(bc,face,cell,idx,side)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx,side

        if (FluxWall_type==FluxWall_Kinetic) then
            call KineticFluxBoundary(bc,face,cell,idx,side)
        elseif (FluxWall_type==FluxWall_Multiscale) then
            if(CouetteUseW2BWall)then
                call W2BMultiscaleFluxBoundary(bc,face,cell,idx,side)
            else
                call MultiscaleFluxBoundary(bc,face,cell,idx,side)
            endif
        elseif (FluxWall_type==FluxWall_Maxwell) then
            call KBC_Maxwell(bc,face,cell,idx,side)
        else
            stop "Error: illegal FluxWall_type!"
        end if
    end subroutine CalcFluxBoundary

    
    !-------------------------------------------------------------------
    !>Calculate kinetic flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------------------------
    subroutine KineticFluxBoundary(bc,face,cell,idx,side)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx,side
        integer(KINT)                                   :: k !Flattened velocity-space index
        real(KREAL)                                     :: rot
        real(KREAL), allocatable, dimension(:)          :: vn_vec,vt_vec !Pure 1D local-frame velocities
        real(KREAL), allocatable, dimension(:)          :: h_vec,b_vec !Pure 1D incident/boundary distributions
        real(KREAL), allocatable, dimension(:)          :: H0_vec,B0_vec !Pure 1D wall Maxwellian distributions
        real(KREAL), allocatable, dimension(:)          :: delta_vec !Pure 1D wall upwind selector
        real(KREAL)                                     :: flux(4)
        real(KREAL)                                     :: prim(4) !boundary condition in local frame
        real(KREAL)                                     :: prim_g(4),conVars_g(4)
        real(KREAL)                                     :: incidence,reflection
        real(KREAL)                                     :: wk,vnK,vtK,v2 !Reusable velocity-space scalars for Nvel loops
        real(KREAL)                                     :: dt   !local time step
        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !Allocate pure 1D velocity-space temporaries for kinetic wall flux loops.
        allocate(vn_vec(Nvel))
        allocate(vt_vec(Nvel))
        allocate(delta_vec(Nvel))
        allocate(h_vec(Nvel))
        allocate(b_vec(Nvel))
        allocate(H0_vec(Nvel))
        allocate(B0_vec(Nvel))
        
        !Time step
        dt = cell%dt_local
        
        !Frame ratation
        rot = -side

        !Convert velocity space to the wall frame and build the wall upwind selector.
        do k=1,Nvel
            vn_vec(k) = uVel(k)*face%cosx+vVel(k)*face%cosy
            vt_vec(k) =-uVel(k)*face%cosy+vVel(k)*face%cosx
            delta_vec(k) = (sign(UP,vn_vec(k))*rot+1.0_KREAL)/2.0_KREAL
        enddo

        !Boundary condition in local frame
        prim = LocalFrame(bc,face%cosx,face%cosy)

        !--------------------------------------------------
        !Obtain h^{in} and b^{in}, rotation accounts for the right wall
        !--------------------------------------------------
        !Reconstruct incident distributions at the wall through 1D velocity views.
        do k=1,Nvel
            h_vec(k) = cell%h_vec(k)-rot*0.5_KREAL*cell%length(idx)*cell%sh_vec(k,idx)
            b_vec(k) = cell%b_vec(k)-rot*0.5_KREAL*cell%length(idx)*cell%sb_vec(k,idx)
        enddo
        
        !--------------------------------------------------
        !Calculate wall density and Maxwellian distribution
        !--------------------------------------------------
        !Compute wall density coefficients using Nvel loops instead of array temporaries.
        incidence = 0.0_KREAL
        reflection = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            incidence = incidence+wk*vnK*h_vec(k)*(1.0_KREAL-delta_vec(k))
            reflection = reflection+wk*vnK*exp(-prim(4)*((vnK-prim(2))**2+(vtK-prim(3))**2))*delta_vec(k)
        enddo
        reflection = (prim(4)/PI)*reflection

        prim(1) = -incidence/reflection

        call DiscreteMaxwell_LocalVec(H0_vec,B0_vec,vn_vec,vt_vec,prim)
        
        !--------------------------------------------------
        !Distribution function at the boundary interface
        !--------------------------------------------------
        !Blend reflected Maxwellian and incident distribution at each velocity point.
        do k=1,Nvel
            h_vec(k) = H0_vec(k)*delta_vec(k)+h_vec(k)*(1.0_KREAL-delta_vec(k))
            b_vec(k) = B0_vec(k)*delta_vec(k)+b_vec(k)*(1.0_KREAL-delta_vec(k))
        enddo
        
        !--------------------------------------------------
        !Calculate flux
        !--------------------------------------------------
        !Accumulate conservative wall flux and boundary conservative state in one Nvel loop.
        flux = 0.0_KREAL
        conVars_g = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            flux(1) = flux(1)+wk*vnK*h_vec(k)
            flux(2) = flux(2)+wk*vnK*vnK*h_vec(k)
            flux(3) = flux(3)+wk*vnK*vtK*h_vec(k)
            flux(4) = flux(4)+0.5_KREAL*wk*vnK*(v2*h_vec(k)+b_vec(k))

            conVars_g(1) = conVars_g(1)+wk*h_vec(k)
            conVars_g(2) = conVars_g(2)+wk*vnK*h_vec(k)
            conVars_g(3) = conVars_g(3)+wk*vtK*h_vec(k)
            conVars_g(4) = conVars_g(4)+0.5_KREAL*wk*(v2*h_vec(k)+b_vec(k))
        enddo

        !--------------------------------------------------
        !Store surface FlowVars and Aerodynamics
        !--------------------------------------------------
        !Flow variables
        prim_g       = GlobalFrame(GetPrimary(conVars_g),face%cosx,face%cosy)
        face%FlowVars(1:3)   = prim_g(1:3)                      !density,u,v
        face%FlowVars(4)     = 1.0/prim_g(4)                    !temperature
        face%FlowVars(5)     = 0.5*prim_g(1)/prim_g(4)          !pressure
        !Aerodynamics (Notice: here assume Unw==0, i.e. prim(2)==0)
        face%Aerodynamics(1) = flux(2)                          !pnn: pressure
        face%Aerodynamics(2) = flux(3)                          !pnt: shear stress
        face%Aerodynamics(3) = flux(4)-prim(3)*flux(3)          !qn:  heat flux
        
        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        flux = GlobalFrame(flux,face%cosx,face%cosy)

        !Total flux
        flux = dt*face%length*flux
        if(side==Lside) then
            face%L_flux     = flux
            !Store micro flux through 1D face views to avoid flux_h/flux_b temporaries.
            do k=1,Nvel
                face%L_flux_h_vec(k) = dt*face%length*vn_vec(k)*h_vec(k)
                face%L_flux_b_vec(k) = dt*face%length*vn_vec(k)*b_vec(k)
            enddo
        elseif(side==Rside) then
            face%R_flux     = flux
            !Store micro flux through 1D face views to avoid flux_h/flux_b temporaries.
            do k=1,Nvel
                face%R_flux_h_vec(k) = dt*face%length*vn_vec(k)*h_vec(k)
                face%R_flux_b_vec(k) = dt*face%length*vn_vec(k)*b_vec(k)
            enddo
        endif
            
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Release pure 1D velocity-space temporaries allocated above.
        deallocate(vn_vec)
        deallocate(vt_vec)
        deallocate(delta_vec)
        deallocate(h_vec)
        deallocate(b_vec)
        deallocate(H0_vec)
        deallocate(B0_vec)
    end subroutine KineticFluxBoundary

    
    !--------------------------------------------------------------------
    !>Calculate multiscale flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------------------------
    subroutine MultiscaleFluxBoundary(bc,face,cell,idx,side)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx,side
        integer(KINT)                                   :: k !Flattened velocity-space index
        real(KREAL)                                     :: rot
        real(KREAL), allocatable, dimension(:)          :: vn_vec,vt_vec !Pure 1D local-frame velocities
        real(KREAL), allocatable, dimension(:)          :: h_vec,b_vec !Pure 1D reconstructed wall distributions
        real(KREAL), allocatable, dimension(:)          :: H_g_vec,B_g_vec !Pure 1D interface Maxwellian arrays
        real(KREAL), allocatable, dimension(:)          :: H_w_vec,B_w_vec !Pure 1D wall Maxwellian arrays
        real(KREAL), allocatable, dimension(:)          :: H_plus_vec,B_plus_vec !Pure 1D Shakhov correction arrays
        real(KREAL)                                     :: qf(2)
        real(KREAL), allocatable, dimension(:)          :: delta_vec !Pure 1D wall upwind selector
        real(KREAL)                                     :: flux(4)
        real(KREAL)                                     :: flux_wall(4) !Accumulated wall macro-flux correction from Nvel loop
        real(KREAL)                                     :: prim(4)
        real(KREAL)                                     :: prim_g(4) !Primary variables of g_{i+1/2}
        real(KREAL)                                     :: prim_w(4) !boundary condition in local frame
        real(KREAL)                                     :: prim_g0(4),conVars_g0(4),conVars_g(4),conVars_g1(4)
        real(KREAL)                                     :: incidence,reflection
        !real(KREAL)                                     :: tau,T1,T4,T5
        real(KREAL)                                     :: tau,Mt(5)
        real(KREAL)                                     :: wk,vnK,vtK,v2 !Reusable velocity-space scalars for Nvel loops
        real(KREAL)                                     :: deltaK,inflowK !Wall reflected/incident selectors for Nvel loops
        real(KREAL)                                     :: H_gK,B_gK,H_wK,B_wK !Per-velocity Maxwellian values for micro flux
        real(KREAL)                                     :: H_plusK,B_plusK !Per-velocity Shakhov correction values for micro flux
        real(KREAL)                                     :: hK,bK !Per-velocity reconstructed distribution values
        real(KREAL)                                     :: shnK,sbnK,shtK,sbtK !Per-velocity distribution slopes
        real(KREAL)                                     :: slope_hK,slope_bK !Per-velocity normal/tangential slope combinations
        real(KREAL)                                     :: wallExpK !Per-velocity Maxwellian kernel used by wall reflection
        real(KREAL)                                     :: Ha_gnK,Ha_gtK,Ha_TK !Per-velocity h-space expansion terms
        real(KREAL)                                     :: Ba_gnK,Ba_gtK,Ba_TK !Per-velocity b-space expansion terms
        real(KREAL)                                     :: dt   !local time step

        ! additional variables for slope calculation
        real(KREAL), allocatable, dimension(:)          :: shn_vec,sbn_vec !Pure 1D normal distribution slopes
        !real(KREAL), allocatable, dimension(:,:)        :: Ha_gn, Ba_gn
        !real(KREAL), allocatable, dimension(:,:)        :: Ha_gt, Ba_gt 
        !real(KREAL), allocatable, dimension(:,:)        :: Ha_T, Ba_T  
        real(KREAL)                                     :: a_gn(4), a_gt(4), a_T(4)    
        real(KREAL)                                     :: sw(4), sw_T(4)
        real(KREAL)                                     :: Mau(4), Mbv(4), Mau0(4), MauT(4)
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM)
        real(KREAL)                                     :: Mv(0:MTUM),Mxi(0:2)

        ! tangential slope variables
        real(KREAL), allocatable, dimension(:)          :: sht_vec,sbt_vec !Pure 1D tangential distribution slopes
        !limiter of tangential slope
        !real(KREAL)                                     :: tan_len, d_rho_t, limiter_phi
        !real(KREAL), parameter                          :: EPS = 1.0e-10
        !real(KREAL), parameter                          :: LIMIT_THRESHOLD = 0.5  ! allow max 50% change in tangential direction



        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !Allocate pure 1D velocity-space temporaries for multiscale wall flux loops.
        allocate(vn_vec(Nvel))
        allocate(vt_vec(Nvel))
        allocate(delta_vec(Nvel))
        allocate(h_vec(Nvel))
        allocate(b_vec(Nvel))
        allocate(H_plus_vec(Nvel))
        allocate(B_plus_vec(Nvel))
        allocate(H_g_vec(Nvel))
        allocate(B_g_vec(Nvel))
        allocate(H_w_vec(Nvel))
        allocate(B_w_vec(Nvel))
        allocate(shn_vec(Nvel))
        allocate(sbn_vec(Nvel))
        allocate(sht_vec(Nvel))
        allocate(sbt_vec(Nvel))

        !Time step
        dt = cell%dt_local
        
        !Frame rotation
        rot = -side

        !Convert velocity space to the wall frame and build the wall upwind selector.
        do k=1,Nvel
            vn_vec(k) = uVel(k)*face%cosx+vVel(k)*face%cosy
            vt_vec(k) =-uVel(k)*face%cosy+vVel(k)*face%cosx
            delta_vec(k) = (sign(UP,vn_vec(k))*rot+1.0_KREAL)/2.0_KREAL
        enddo

        !--------------------------------------------------
        !Reconstruct non-equlibrium distribution at interface
        !--------------------------------------------------
        !Load normal and tangential slopes through 1D velocity views.
        do k=1,Nvel
            shn_vec(k) = cell%sh_vec(k,idx)
            sbn_vec(k) = cell%sb_vec(k,idx)

            if (idx == 1) then ! normal is x, tangential is y
                sht_vec(k) = cell%sh_vec(k,2)
                sbt_vec(k) = cell%sb_vec(k,2)
                !tan_len = cell%length(2) ! tangential length
            else ! normal is y, tangential is x
                sht_vec(k) = -cell%sh_vec(k,1)
                sbt_vec(k) = -cell%sb_vec(k,1)
                !tan_len = cell%length(1) ! tangential length
            end if
        enddo

        !Limiter
        ! calculaate the density change due to tangential slope
        !d_rho_t = sum(weight*sht) * (0.5 * tan_len)
        
        ! If the density change > 50% of the center density, calculate a reduction factor phi
        ! prim(1) is the density of the current cell
        !prim = GetPrimary(cell%conVars) ! temporarily get for comparison
        !
        !if (abs(d_rho_t) > LIMIT_THRESHOLD * prim(1)) then
        !    limiter_phi = (LIMIT_THRESHOLD * prim(1)) / (abs(d_rho_t) + EPS)
        !else
        !    limiter_phi = 1.0
        !end if
        
        !sht = sht * limiter_phi
        !sbt = sbt * limiter_phi
        
        ! if use tangential slope sht,sbt, all Mt(5) need change vn*sh to (vn*sh + vt*sht),vn*sb to (vn*sb + vt*sbt)

        !Reconstruct incident non-equilibrium distribution at the wall.
        do k=1,Nvel
            h_vec(k) = cell%h_vec(k)-rot*0.5_KREAL*cell%length(idx)*shn_vec(k)
            b_vec(k) = cell%b_vec(k)-rot*0.5_KREAL*cell%length(idx)*sbn_vec(k)
        enddo

        !Calculate primary variable of g_{i+1/2} in local frame
        prim_w = LocalFrame(bc,face%cosx,face%cosy)
        prim_g = prim_w
        prim = LocalFrame(GetPrimary(cell%conVars),face%cosx,face%cosy)
        prim_g(1) = prim(1)/prim(4)*prim_w(4) !keep consistent pressure inside cell

        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim_g)

        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)

        !--------------------------------------------------
        !1. Calculate Maxwellian distribution and Shakhov part at g_{i+1/2}
        !--------------------------------------------------
        call DiscreteMaxwell_LocalVec(H_g_vec,B_g_vec,vn_vec,vt_vec,prim_g)

        qf = GetHeatFlux_LocalVec(h_vec,b_vec,vn_vec,vt_vec,prim_g)
        call ShakhovPart_LocalVec(H_g_vec,B_g_vec,vn_vec,vt_vec,qf,prim_g,H_plus_vec,B_plus_vec)

        !--------------------------------------------------
        !2. Spatial slope a_g: derived from local gradients
        !--------------------------------------------------
        !sw = LocalFrame(cell%sw(:,idx), face%cosx, face%cosy) 
        !a_g  = MicroSlope(prim_g, sw)
        
        !Alternative: derived from difference between cell center and interface
        conVars_g = GetConserved(prim_g)
        conVars_g1 = LocalFrame(cell%conVars, face%cosx, face%cosy)
        sw = (conVars_g - conVars_g1) / (0.5 * cell%length(idx))
        a_gn  = MicroSlope(prim_g, sw)
        
        !if (idx == 1) then
        !    sw_t = LocalFrame(cell%sw(:,2), face%cosx, face%cosy)
        !else
        !    sw_t = LocalFrame(-cell%sw(:,1), face%cosx, face%cosy) 
        !end if
        
        !Tangential macro slope obtained by integrating reconstructed VDF slopes.
        sw_t = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            sw_t(1) = sw_t(1)+wk*sht_vec(k)
            sw_t(2) = sw_t(2)+wk*vnK*sht_vec(k)
            sw_t(3) = sw_t(3)+wk*vtK*sht_vec(k)
            sw_t(4) = sw_t(4)+0.5_KREAL*wk*(v2*sht_vec(k)+sbt_vec(k))
        enddo

        a_gt = MicroSlope(prim_g, sw_t)
        !--------------------------------------------------
        !3. Temporal slope a_T
        !--------------------------------------------------
        call CalcMoment(prim_g, Mu, Mv, Mxi, MuL, MuR)
        Mau = Moment_auvxi(a_gn, Mu, Mv, Mxi, 1, 0) 
        Mbv = Moment_auvxi(a_gt, Mu, Mv, Mxi, 0, 1)
        sw = -prim_g(1)*(Mau + Mbv) 
        a_T = MicroSlope(prim_g, sw)
        
        !--------------------------------------------------
        !4. Prepare expansion terms for Flux calculation
        !--------------------------------------------------
        !Ha_gn = (a_gn(1) + a_gn(2)*vn + a_gn(3)*vt)*H_g + 0.5*a_gn(4)*((vn**2+vt**2)*H_g + B_g)
        !Ba_gn = (a_gn(1) + a_gn(2)*vn + a_gn(3)*vt)*B_g + 0.5*a_gn(4)*((vn**2+vt**2)*B_g + Mxi(2)*H_g)
        !Ha_gt = (a_gt(1) + a_gt(2)*vn + a_gt(3)*vt)*H_g + 0.5*a_gt(4)*((vn**2+vt**2)*H_g + B_g)
        !Ba_gt = (a_gt(1) + a_gt(2)*vn + a_gt(3)*vt)*B_g + 0.5*a_gt(4)*((vn**2+vt**2)*B_g + Mxi(2)*H_g)
        !Ha_T  = (a_T(1) + a_T(2)*vn + a_T(3)*vt)*H_g + 0.5*a_T(4)*((vn**2+vt**2)*H_g + B_g)
        !Ba_T  = (a_T(1) + a_T(2)*vn + a_T(3)*vt)*B_g + 0.5*a_T(4)*((vn**2+vt**2)*B_g + Mxi(2)*H_g)
        

        if (side == Lside) then
            Mu = MuR
        else
            Mu = MuL
        end if

        Mau0 = Moment_uvxi(Mu,Mv,Mxi,1,0,0) !<u*\psi>
        Mau = Moment_auvxi(a_gn,Mu,Mv,Mxi,2,0) !<a*u^2*\psi>
        Mbv  = Moment_auvxi(a_gt,Mu,Mv,Mxi,1,1) !<b*u*v*\psi>
        MauT = Moment_auvxi(a_T,Mu,Mv,Mxi,1,0) !<A*u*\psi>
        
        ! flux of Maxwellian part
        flux = Mt(1)*prim_g(1)*Mau0+Mt(2)*prim_g(1)*(Mau+Mbv)+Mt(3)*prim_g(1)!*MauT
        !--------------------------------------------------
        !Calculate wall density and Maxwellian distribution
        !--------------------------------------------------
        
        !incidence = Mt(1)*sum(weight*vn*(H_g+H_plus)*(1-delta)) + &
        !            Mt(2)*sum(weight*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + & 
        !            Mt(3)*sum(weight*vn*Ha_T*(1-delta)) + &      
        !            Mt(4)*sum(weight*vn*h*(1-delta)) - &
        !            Mt(5)*sum(weight*vn*(vn*shn + vt*sht)*(1-delta))
        
        !Compute wall incidence/reflection using 1D velocity views.
        incidence  = flux(1)
        reflection = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK
            slope_hK = vnK*shn_vec(k)+vtK*sht_vec(k)
            wallExpK = exp(-prim_w(4)*((vnK-prim_w(2))*(vnK-prim_w(2))+(vtK-prim_w(3))*(vtK-prim_w(3))))

            incidence  = incidence+wk*vnK*inflowK*(Mt(1)*H_plus_vec(k)+Mt(4)*h_vec(k)-Mt(5)*slope_hK)
            reflection = reflection+wk*vnK*wallExpK*deltaK
        enddo
        reflection = dt*(prim_w(4)/PI)*reflection
        prim_w(1) = -incidence/reflection

        call DiscreteMaxwell_LocalVec(H_w_vec,B_w_vec,vn_vec,vt_vec,prim_w)

        !--------------------------------------------------
        ! Calculate Fluxes 
        !--------------------------------------------------
        ! 1. Mass Flux
        !flux(1) = dt*sum(weight*vn*H_w*delta) + &
        !          Mt(1)*sum(weight*vn*(H_g+H_plus)*(1-delta)) + &
        !          Mt(2)*sum(weight*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta))+ &   
        !          Mt(3)*sum(weight*vn*Ha_T*(1-delta)) + &        
        !          Mt(4)*sum(weight*vn*h*(1-delta)) - &
        !          Mt(5)*sum(weight*vn*(vn*shn + vt*sht)*(1-delta)) 
        !Accumulate reflected and incident macro flux corrections through 1D velocity views.
        flux_wall = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK

            H_wK     = H_w_vec(k)
            B_wK     = B_w_vec(k)
            H_plusK  = H_plus_vec(k)
            B_plusK  = B_plus_vec(k)
            hK       = h_vec(k)
            bK       = b_vec(k)
            slope_hK = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_bK = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            flux_wall(1) = flux_wall(1)+wk*vnK*(dt*H_wK*deltaK+Mt(1)*H_plusK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            flux_wall(2) = flux_wall(2)+wk*vnK*vnK*(dt*H_wK*deltaK+Mt(1)*H_plusK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            flux_wall(3) = flux_wall(3)+wk*vnK*vtK*(dt*H_wK*deltaK+Mt(1)*H_plusK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            flux_wall(4) = flux_wall(4)+0.5_KREAL*wk*vnK*(dt*(v2*H_wK+B_wK)*deltaK+ &
                Mt(1)*(v2*H_plusK+B_plusK)*inflowK+Mt(4)*(v2*hK+bK)*inflowK- &
                Mt(5)*(v2*slope_hK+slope_bK)*inflowK)
        enddo
        flux = flux+flux_wall

        ! 2. Normal Momentum Flux
        !flux(2) = dt*sum(weight*vn*vn*H_w*delta) + &
        !          Mt(1)*sum(weight*vn*vn*(H_g+H_plus)*(1-delta)) + &
        !          Mt(2)*sum(weight*vn*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + & 
        !          Mt(3)*sum(weight*vn*vn*Ha_T*(1-delta)) + &     
        !          Mt(4)*sum(weight*vn*vn*h*(1-delta)) - &
        !          Mt(5)*sum(weight*vn*vn*(vn*shn + vt*sht)*(1-delta)) 
        ! The Nvel loop above replaces the four separate macro-flux sum expressions.
        ! 3. Tangential Momentum Flux
        !flux(3) = dt*sum(weight*vn*vt*H_w*delta) + &
        !          Mt(1)*sum(weight*vn*vt*(H_g+H_plus)*(1-delta)) + &
        !          Mt(2)*sum(weight*vn*vt*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + & 
        !          Mt(3)*sum(weight*vn*vt*Ha_T*(1-delta)) + &     
        !          Mt(4)*sum(weight*vt*vn*h*(1-delta)) - &
        !          Mt(5)*sum(weight*vt*vn*(vn*shn + vt*sht)*(1-delta)) 
        ! The Nvel loop above replaces the tangential momentum flux sum.
        ! 4. Energy Flux
        !flux(4) = dt*0.5*sum(weight*vn*((vn**2+vt**2)*H_w+B_w)*delta)+ &
        !          Mt(1)*0.5*sum(weight*vn*((vn**2+vt**2)*(H_g+H_plus)+(B_g+B_plus))*(1-delta))+ &
        !          Mt(2)*0.5*sum(weight*vn*((vn**2+vt**2)*(vn*Ha_gn+vt*Ha_gt)+(vn*Ba_gn+vt*Ba_gt))*(1-delta))+ & 
        !          Mt(3)*0.5*sum(weight*vn*((vn**2+vt**2)*Ha_T+Ba_T)*(1-delta))+ &      
        !          Mt(4)*0.5*sum(weight*vn*((vn**2+vt**2)*h+b)*(1-delta)) - &
        !          Mt(5)*0.5*sum(weight*vn*((vn**2+vt**2)*(vn*shn + vt*sht) + (vn*sbn + vt*sbt))*(1-delta))
        ! The Nvel loop above replaces the energy flux sum.



        !flux_h = dt*vn*H_w*delta + &
        !         Mt(1)*vn*(H_g+H_plus)*(1-delta) + &
        !         Mt(2)*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta) + &  
        !         Mt(3)*vn*Ha_T*(1-delta) + &       
        !         Mt(4)*vn*h*(1-delta) - &
        !         Mt(5)*vn*(vn*shn + vt*sht)*(1-delta) 
        
        !Store multiscale wall micro flux directly through 1D face views.
        do k=1,Nvel
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK

            H_gK     = H_g_vec(k)
            B_gK     = B_g_vec(k)
            H_wK     = H_w_vec(k)
            B_wK     = B_w_vec(k)
            H_plusK  = H_plus_vec(k)
            B_plusK  = B_plus_vec(k)
            hK       = h_vec(k)
            bK       = b_vec(k)
            shnK     = shn_vec(k)
            sbnK     = sbn_vec(k)
            shtK     = sht_vec(k)
            sbtK     = sbt_vec(k)

            Ha_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*H_gK+0.5_KREAL*a_gn(4)*(v2*H_gK+B_gK)
            Ha_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*H_gK+0.5_KREAL*a_gt(4)*(v2*H_gK+B_gK)
            Ha_TK  = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*H_gK+0.5_KREAL*a_T(4)*(v2*H_gK+B_gK)

            Ba_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*B_gK+0.5_KREAL*a_gn(4)*(v2*B_gK+Mxi(2)*H_gK)
            Ba_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*B_gK+0.5_KREAL*a_gt(4)*(v2*B_gK+Mxi(2)*H_gK)
            Ba_TK  = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*B_gK+0.5_KREAL*a_T(4)*(v2*B_gK+Mxi(2)*H_gK)

            if(side==Lside) then
                face%L_flux_h_vec(k) = face%length*(dt*vnK*H_wK*deltaK + &
                    Mt(1)*vnK*H_plusK*inflowK + Mt(2)*vnK*vnK*Ha_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ha_gtK*inflowK + Mt(3)*vnK*Ha_TK*inflowK + &
                    Mt(4)*vnK*hK*inflowK - Mt(5)*vnK*(vnK*shnK+vtK*shtK)*inflowK)

                face%L_flux_b_vec(k) = face%length*(dt*vnK*B_wK*deltaK + &
                    Mt(1)*vnK*B_plusK*inflowK + Mt(2)*vnK*vnK*Ba_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ba_gtK*inflowK + Mt(3)*vnK*Ba_TK*inflowK + &
                    Mt(4)*vnK*bK*inflowK - Mt(5)*vnK*(vnK*sbnK+vtK*sbtK)*inflowK)
            elseif(side==Rside) then
                face%R_flux_h_vec(k) = face%length*(dt*vnK*H_wK*deltaK + &
                    Mt(1)*vnK*H_plusK*inflowK + Mt(2)*vnK*vnK*Ha_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ha_gtK*inflowK + Mt(3)*vnK*Ha_TK*inflowK + &
                    Mt(4)*vnK*hK*inflowK - Mt(5)*vnK*(vnK*shnK+vtK*shtK)*inflowK)

                face%R_flux_b_vec(k) = face%length*(dt*vnK*B_wK*deltaK + &
                    Mt(1)*vnK*B_plusK*inflowK + Mt(2)*vnK*vnK*Ba_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ba_gtK*inflowK + Mt(3)*vnK*Ba_TK*inflowK + &
                    Mt(4)*vnK*bK*inflowK - Mt(5)*vnK*(vnK*sbnK+vtK*sbtK)*inflowK)
            endif
        enddo
        !flux_b = dt*vn*B_w*delta + &
        !         Mt(1)*vn*(B_g+B_plus)*(1-delta) + &
        !         Mt(2)*vn*(vn*Ba_gn + vt*Ba_gt)*(1-delta) + &   
        !         Mt(3)*vn*Ba_T*(1-delta) + &      
        !         Mt(4)*vn*b*(1-delta) - &
        !         Mt(5)*vn*(vn*sbn + vt*sbt)*(1-delta) 

        !--------------------------------------------------
        !Store surface FlowVars and Aerodynamics
        !--------------------------------------------------
        !Flow variables
        
        !conVars_g0(1) = dt*sum(weight*H_w*delta) + &
        !                Mt(1)*sum(weight*(H_g+H_plus)*(1-delta)) + &
        !                Mt(2)*sum(weight*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + &   
        !                Mt(3)*sum(weight*Ha_T*(1-delta)) + &       
        !                Mt(4)*sum(weight*h*(1-delta)) - &
        !                Mt(5)*sum(weight*(vn*shn + vt*sht)*(1-delta)) 

        !conVars_g0(2) = dt*sum(weight*vn*H_w*delta) + &
        !                Mt(1)*sum(weight*vn*(H_g+H_plus)*(1-delta)) + &
        !                Mt(2)*sum(weight*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + &
        !                Mt(3)*sum(weight*vn*Ha_T*(1-delta)) + &    
        !                Mt(4)*sum(weight*vn*h*(1-delta)) - &
        !                Mt(5)*sum(weight*vn*(vn*shn + vt*sht)*(1-delta)) 

        !conVars_g0(3) = dt*sum(weight*vt*H_w*delta) + &
        !                Mt(1)*sum(weight*vt*(H_g+H_plus)*(1-delta)) + &
        !                Mt(2)*sum(weight*vt*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + &
        !                Mt(3)*sum(weight*vt*Ha_T*(1-delta)) + &    
        !                Mt(4)*sum(weight*vt*h*(1-delta)) - &
        !                Mt(5)*sum(weight*vt*(vn*shn + vt*sht)*(1-delta)) 

        !conVars_g0(4) = dt*0.5*sum(weight*((vn**2+vt**2)*H_w+B_w)*delta)+ &
        !                Mt(1)*0.5*sum(weight*((vn**2+vt**2)*(H_g+H_plus)+(B_g+B_plus))*(1-delta))+ &
        !                Mt(2)*0.5*sum(weight*((vn**2+vt**2)*(vn*Ha_gn + vt*Ha_gt) + (vn*Ba_gn + vt*Ba_gt))*(1-delta))+ &
        !                Mt(3)*0.5*sum(weight*((vn**2+vt**2)*Ha_T+Ba_T)*(1-delta))+ &      
        !                Mt(4)*0.5*sum(weight*((vn**2+vt**2)*h+b)*(1-delta)) - &
        !                Mt(5)*0.5*sum(weight*((vn**2+vt**2)*(vn*shn + vt*sht) + (vn*sbn + vt*sbt))*(1-delta))

        !Build wall output conserved variables through a single 1D velocity loop.
        conVars_g0 = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK

            H_gK     = H_g_vec(k)
            B_gK     = B_g_vec(k)
            H_wK     = H_w_vec(k)
            B_wK     = B_w_vec(k)
            H_plusK  = H_plus_vec(k)
            B_plusK  = B_plus_vec(k)
            hK       = h_vec(k)
            bK       = b_vec(k)
            slope_hK = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_bK = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            Ha_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*H_gK+0.5_KREAL*a_gn(4)*(v2*H_gK+B_gK)
            Ha_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*H_gK+0.5_KREAL*a_gt(4)*(v2*H_gK+B_gK)
            Ha_TK  = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*H_gK+0.5_KREAL*a_T(4)*(v2*H_gK+B_gK)

            conVars_g0(1) = conVars_g0(1)+wk*(dt*H_wK*deltaK+Mt(1)*(H_gK+H_plusK)*inflowK+ &
                Mt(2)*vnK*Ha_gnK*inflowK+Mt(2)*vtK*Ha_gtK*inflowK+Mt(3)*Ha_TK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            conVars_g0(2) = conVars_g0(2)+wk*vnK*(dt*H_wK*deltaK+Mt(1)*(H_gK+H_plusK)*inflowK+ &
                Mt(2)*vnK*Ha_gnK*inflowK+Mt(2)*vtK*Ha_gtK*inflowK+Mt(3)*Ha_TK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            conVars_g0(3) = conVars_g0(3)+wk*vtK*(dt*H_wK*deltaK+Mt(1)*(H_gK+H_plusK)*inflowK+ &
                Mt(2)*vnK*Ha_gnK*inflowK+Mt(2)*vtK*Ha_gtK*inflowK+Mt(3)*Ha_TK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)

            !Keep the historical conVars_g0(4) algebra exactly; do not substitute the commented formula here.
            conVars_g0(4) = conVars_g0(4)+0.5_KREAL*wk*(dt*(v2*H_wK+B_wK)*deltaK+ &
                Mt(1)*(v2*(H_gK+H_plusK)+(B_gK+B_plusK))*inflowK+ &
                Mt(2)*vnK*(v2*Ha_gnK+Ha_gtK*vtK)*inflowK+ &
                Mt(2)*vtK*(v2*Ha_gnK+Ha_gtK*vtK)*inflowK+Mt(3)*v2*Ha_TK*inflowK+ &
                Mt(4)*(v2*hK+bK)*inflowK-Mt(5)*(v2*slope_hK+slope_bK)*inflowK)
        enddo
        

        prim_g0       = GlobalFrame(GetPrimary(conVars_g0/dt),face%cosx,face%cosy)
        face%FlowVars(1:3)   = prim_g0(1:3)                         !density,u,v
        face%FlowVars(4)     = 1.0/prim_g0(4)                       !temperature
        face%FlowVars(5)     = 0.5*prim_g0(1)/prim_g0(4)            !pressure
        !Aerodynamics (Notice: here assume Unw==0, i.e. prim_w(2)==0)
        face%Aerodynamics(1) = flux(2)/dt                           !pnn: pressure
        face%Aerodynamics(2) = flux(3)/dt                           !pnt: shear stress
        face%Aerodynamics(3) = flux(4)/dt-prim_w(3)*flux(3)/dt      !qn:  heat flux

        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        
        flux = GlobalFrame(flux,face%cosx,face%cosy)

        !Total flux
        flux = face%length*flux

        if(side==Lside) then
            face%L_flux     = flux
        elseif(side==Rside) then
            face%R_flux     = flux
        endif
        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Release pure 1D velocity-space temporaries allocated above.
        deallocate(vn_vec)
        deallocate(vt_vec)
        deallocate(delta_vec)
        deallocate(h_vec)
        deallocate(b_vec)
        deallocate(H_g_vec)
        deallocate(B_g_vec)
        deallocate(H_w_vec)
        deallocate(B_w_vec)
        deallocate(H_plus_vec)
        deallocate(B_plus_vec)
        deallocate(shn_vec)
        deallocate(sbn_vec)
        !deallocate(Ha_gn); deallocate(Ba_gn)
        !deallocate(Ha_gt); deallocate(Ba_gt)
        !deallocate(Ha_T); deallocate(Ba_T)
        deallocate(sht_vec)
        deallocate(sbt_vec)
    end subroutine MultiscaleFluxBoundary
    
    
    !-------------------------------------------------------------------
    !>Calculate kinetic flux of boundary interface, assuming left wall
    !>@param[in]    bc   :boundary condition
    !>@param[inout] face :the boundary interface
    !>@param[in]    cell :cell next to the boundary interface
    !>@param[in]    idx  :index indicating i or j direction
    !>@param[in]    rot  :indicating rotation
    !--------------------------------------------------------------------
    subroutine KBC_Maxwell(bc,face,cell,idx,side)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx,side
        integer(KINT)                                   :: k !Flattened velocity-space index
        real(KREAL)                                     :: rot
        real(KREAL), allocatable, dimension(:)          :: vn_vec,vt_vec !Pure 1D local-frame velocities
        real(KREAL), allocatable, dimension(:)          :: h_vec,b_vec !Pure 1D wall/interface distributions
        real(KREAL), allocatable, dimension(:)          :: Hw_vec,Bw_vec !Pure 1D wall Maxwellian arrays
        real(KREAL), allocatable, dimension(:)          :: Hr_vec,Br_vec !Pure 1D specular-reflection arrays
        real(KREAL), allocatable, dimension(:)          :: delta_vec !Pure 1D wall upwind selector
        real(KREAL)                                     :: flux(4)
        real(KREAL)                                     :: prim(4) !boundary condition in local frame
        real(KREAL)                                     :: prim_g(4),conVars_g(4)
        real(KREAL)                                     :: incidence,reflection
        real(KREAL)                                     :: wk,vnK,vtK,v2 !Reusable velocity-space scalars for Nvel loops
        real(KREAL)                                     :: deltaK,inflowK !Wall reflected/incident selectors for Nvel loops
        real(KREAL)                                     :: hK,bK,HwK,BwK,HrK,BrK !Per-velocity distribution values
        real(KREAL)                                     :: wallExpK !Per-velocity Maxwellian kernel used by wall reflection
        real(KREAL)                                     :: dt   !local time step
        integer(KINT)                                   :: vel_i,vel_j !Structured velocity indices recovered from flat_k when needed
        integer(KINT)                                   :: flat_k,sym_flat_k !Flattened index pair for specular reflection mapping
        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !Allocate pure 1D velocity-space temporaries for Maxwell wall flux loops.
        allocate(vn_vec(Nvel))
        allocate(vt_vec(Nvel))
        allocate(delta_vec(Nvel))
        allocate(h_vec(Nvel))
        allocate(b_vec(Nvel))
        allocate(Hw_vec(Nvel))
        allocate(Bw_vec(Nvel))
        allocate(Hr_vec(Nvel))
        allocate(Br_vec(Nvel))
        
        !Time step
        dt = cell%dt_local
        
        !Frame ratation
        rot = -side

        !Convert velocity space to the wall frame and build the wall upwind selector.
        do k=1,Nvel
            vn_vec(k) = uVel(k)*face%cosx+vVel(k)*face%cosy
            vt_vec(k) =-uVel(k)*face%cosy+vVel(k)*face%cosx
            delta_vec(k) = (sign(UP,vn_vec(k))*rot+1.0_KREAL)/2.0_KREAL
        enddo

        !Boundary condition in local frame
        prim = LocalFrame(bc,face%cosx,face%cosy)

        !--------------------------------------------------
        !Obtain h^{in} and b^{in}, rotation accounts for the right wall
        !--------------------------------------------------
        !Reconstruct incident distributions at the Maxwell wall through 1D views.
        do k=1,Nvel
            h_vec(k) = cell%h_vec(k)-rot*0.5_KREAL*cell%length(idx)*cell%sh_vec(k,idx)
            b_vec(k) = cell%b_vec(k)-rot*0.5_KREAL*cell%length(idx)*cell%sb_vec(k,idx)
        enddo
        
        !--------------------------------------------------
        !Obtain specular distribution: Hr, Br
        !--------------------------------------------------
        !Hr(vn,vt)=H_in(-vn,vt): build the same mirror relation with flattened indices.
        if(abs(face%cosy)<1e-12) then           !normal vector (1, 0) or (-1, 0): symmeter along y axis
            !Reverse u in flat Nvel storage; this preserves the old Nu+1-i mapping.
            do flat_k=1,Nvel
                vel_i = mod(flat_k-1,Nu)+1
                vel_j = (flat_k-1)/Nu+1
                sym_flat_k = (Nu-vel_i+1)+(vel_j-1)*Nu
                Hr_vec(flat_k) = h_vec(sym_flat_k)
                Br_vec(flat_k) = b_vec(sym_flat_k)
            enddo
        elseif(abs(face%cosx)<1e-12) then       !normal vector (0, 1) or (0, -1): symmeter along x axis
            if (Quadrature_method == Quad_Unstruct) then
                !Unstructured v reflection uses the precomputed flat SymMap_V table.
                do flat_k=1,Nvel
                    sym_flat_k = SymMap_V(flat_k)
                    Hr_vec(flat_k) = h_vec(sym_flat_k)
                    Br_vec(flat_k) = b_vec(sym_flat_k)
                end do
            else
                !Reverse v in flat Nvel storage; this preserves the old Nv+1-j mapping.
                do flat_k=1,Nvel
                    vel_i = mod(flat_k-1,Nu)+1
                    vel_j = (flat_k-1)/Nu+1
                    sym_flat_k = vel_i+(Nv-vel_j)*Nu
                    Hr_vec(flat_k) = h_vec(sym_flat_k)
                    Br_vec(flat_k) = b_vec(sym_flat_k)
                enddo
            end if
        else
            !General case: interpolation has to be used
            if (myID==0) print*, 'Warning: general interpolation code for Hr and Br in sub_KBC_Maxwell needs to be added.'
            ! temporarily set Hr = h to avoid crash, although not physically accurate
            Hr_vec = h_vec
            Br_vec = b_vec
        endif
        
        !--------------------------------------------------
        !Calculate wall density and Maxwellian distribution
        !--------------------------------------------------
        !Compute wall density terms through the flattened velocity-space view.
        incidence  = 0.0_KREAL
        reflection = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK
            wallExpK = exp(-prim(4)*((vnK-prim(2))*(vnK-prim(2))+(vtK-prim(3))*(vtK-prim(3))))

            incidence  = incidence+wk*vnK*h_vec(k)*inflowK
            reflection = reflection+wk*vnK*wallExpK*deltaK
        enddo
        reflection = (prim(4)/PI)*reflection

        prim(1) = -incidence/reflection

        call DiscreteMaxwell_LocalVec(Hw_vec,Bw_vec,vn_vec,vt_vec,prim)
        
        !--------------------------------------------------
        !Distribution function at the boundary interface
        !--------------------------------------------------
        !Mix diffuse and specular reflection directly in the flattened wall distributions.
        do k=1,Nvel
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK
            h_vec(k) = (ACC*Hw_vec(k)+(1.0_KREAL-ACC)*Hr_vec(k))*deltaK+h_vec(k)*inflowK
            b_vec(k) = (ACC*Bw_vec(k)+(1.0_KREAL-ACC)*Br_vec(k))*deltaK+b_vec(k)*inflowK
        enddo
        
        !--------------------------------------------------
        !Calculate flux
        !--------------------------------------------------
        !Accumulate Maxwell wall macro fluxes, output variables, and micro fluxes in one Nvel loop.
        flux = 0.0_KREAL
        conVars_g = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK
            hK  = h_vec(k)
            bK  = b_vec(k)

            flux(1) = flux(1)+wk*vnK*hK
            flux(2) = flux(2)+wk*vnK*vnK*hK
            flux(3) = flux(3)+wk*vnK*vtK*hK
            flux(4) = flux(4)+0.5_KREAL*wk*vnK*(v2*hK+bK)

            conVars_g(1) = conVars_g(1)+wk*hK
            conVars_g(2) = conVars_g(2)+wk*vnK*hK
            conVars_g(3) = conVars_g(3)+wk*vtK*hK
            conVars_g(4) = conVars_g(4)+0.5_KREAL*wk*(v2*hK+bK)

            if(side==Lside) then
                face%L_flux_h_vec(k) = dt*face%length*vnK*hK
                face%L_flux_b_vec(k) = dt*face%length*vnK*bK
            elseif(side==Rside) then
                face%R_flux_h_vec(k) = dt*face%length*vnK*hK
                face%R_flux_b_vec(k) = dt*face%length*vnK*bK
            endif
        enddo

        !--------------------------------------------------
        !Store surface FlowVars and Aerodynamics
        !--------------------------------------------------
        !Flow variables
        prim_g       = GlobalFrame(GetPrimary(conVars_g),face%cosx,face%cosy)
        face%FlowVars(1:3)   = prim_g(1:3)                      !density,u,v
        face%FlowVars(4)     = 1.0/prim_g(4)                    !temperature
        face%FlowVars(5)     = 0.5*prim_g(1)/prim_g(4)          !pressure
        !Aerodynamics (Notice: here assume Unw==0, i.e. prim(2)==0)
        face%Aerodynamics(1) = flux(2)                          !pnn: pressure
        face%Aerodynamics(2) = flux(3)                          !pnt: shear stress
        face%Aerodynamics(3) = flux(4)-prim(3)*flux(3)          !qn:  heat flux
        
        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        flux = GlobalFrame(flux,face%cosx,face%cosy)

        !Total flux
        flux = dt*face%length*flux
        if(side==Lside) then
            face%L_flux     = flux
        elseif(side==Rside) then
            face%R_flux     = flux
        endif
            
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Release pure 1D velocity-space temporaries allocated above.
        deallocate(vn_vec)
        deallocate(vt_vec)
        deallocate(delta_vec)
        deallocate(h_vec)
        deallocate(b_vec)
        deallocate(Hw_vec)
        deallocate(Bw_vec)
        deallocate(Hr_vec)
        deallocate(Br_vec)
    end subroutine KBC_Maxwell

    
    !--------------------------------------------------
    !>Calculate micro slope of Maxwellian distribution
    !>@param[in] prim        :primary variables
    !>@param[in] sw          :slope of conVars
    !>@return    MicroSlope  :slope of Maxwellian distribution
    !--------------------------------------------------
    function MicroSlope(prim,sw)
        real(KREAL), intent(in)                         :: prim(4),sw(4)
        real(KREAL)                                     :: MicroSlope(4)

        MicroSlope(4) = 4.0*prim(4)**2/((CK+2)*prim(1))*(2.0*sw(4)-2.0*prim(2)*sw(2)-2.0*prim(3)*sw(3)+sw(1)*(prim(2)**2+prim(3)**2-0.5*(CK+2)/prim(4)))
        MicroSlope(3) = 2.0*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*MicroSlope(4)
        MicroSlope(2) = 2.0*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*MicroSlope(4)
        MicroSlope(1) = sw(1)/prim(1)-prim(2)*MicroSlope(2)-prim(3)*MicroSlope(3)-0.5*(prim(2)**2+prim(3)**2+0.5*(CK+2)/prim(4))*MicroSlope(4)
    end function MicroSlope
    
    
    !--------------------------------------------------
    !>calculate moments of velocity and \xi
    !>@param[in] prim       :primary variables
    !>@param[out] Mu,Mv     :<u^n>,<v^m>
    !>@param[out] Mxi       :<\xi^2n>
    !>@param[out] MuL,MuR   :<u^n>_{>0},<u^n>_{<0}
    !--------------------------------------------------
    subroutine CalcMoment(prim,Mu,Mv,Mxi,MuL,MuR)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL), intent(out)                        :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM)
        real(KREAL), intent(out)                        :: Mv(0:MTUM)
        real(KREAL), intent(out)                        :: Mxi(0:2)
        integer :: i

        !Moments of normal velocity
        MuL(0) = 0.5*erfc(-sqrt(prim(4))*prim(2))
        MuL(1) = prim(2)*MuL(0)+0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))
        MuR(0) = 0.5*erfc(sqrt(prim(4))*prim(2))
        MuR(1) = prim(2)*MuR(0)-0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))

        do i=2,MNUM
            MuL(i) = prim(2)*MuL(i-1)+0.5*(i-1)*MuL(i-2)/prim(4)
            MuR(i) = prim(2)*MuR(i-1)+0.5*(i-1)*MuR(i-2)/prim(4)
        end do

        Mu = MuL+MuR

        !Moments of tangential velocity
        Mv(0) = 1.0
        Mv(1) = prim(3)

        do i=2,MTUM
            Mv(i) = prim(3)*Mv(i-1)+0.5*(i-1)*Mv(i-2)/prim(4)
        end do

        !Moments of \xi
        Mxi(0) = 1.0 !<\xi^0>
        Mxi(1) = 0.5*CK/prim(4) !<\xi^2>
        Mxi(2) = CK*(CK+2.0)/(4.0*prim(4)**2) !<\xi^4>
    end subroutine CalcMoment

    
    !--------------------------------------------------
    !>Calculate <u^\alpha*v^\beta*\xi^\delta*\psi>
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^delta>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@param[in] delta      :exponential index of \xi
    !>@return    Moment_uvxi :moment of <u^\alpha*v^\beta*\xi^\delta*\psi>
    !--------------------------------------------------
    function Moment_uvxi(Mu,Mv,Mxi,alpha,beta,delta)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,beta,delta
        real(KREAL)                                     :: Moment_uvxi(4)

        Moment_uvxi(1) = Mu(alpha)*Mv(beta)*Mxi(delta/2)
        Moment_uvxi(2) = Mu(alpha+1)*Mv(beta)*Mxi(delta/2)
        Moment_uvxi(3) = Mu(alpha)*Mv(beta+1)*Mxi(delta/2)
        Moment_uvxi(4) = 0.5*(Mu(alpha+2)*Mv(beta)*Mxi(delta/2)+Mu(alpha)*Mv(beta+2)*Mxi(delta/2)+Mu(alpha)*Mv(beta)*Mxi((delta+2)/2))
    end function Moment_uvxi

    
    !--------------------------------------------------
    !>Calculate <a*u^\alpha*v^\beta>
    !>@param[in] a          :micro slope of Maxwellian
    !>@param[in] Mu,Mv      :<u^\alpha>,<v^\beta>
    !>@param[in] Mxi        :<\xi^l>
    !>@param[in] alpha,beta :exponential index of u and v
    !>@return    Moment_auvxi  :moment of <a*u^\alpha*v^\beta>
    !--------------------------------------------------
    function Moment_auvxi(a,Mu,Mv,Mxi,alpha,beta)
        real(KREAL), intent(in)                         :: a(4)
        real(KREAL), intent(in)                         :: Mu(0:MNUM),Mv(0:MTUM),Mxi(0:2)
        integer(KINT), intent(in)                       :: alpha,beta
        real(KREAL)                                     :: Moment_auvxi(4)

        Moment_auvxi = a(1)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+0,0)+&
                    a(2)*Moment_uvxi(Mu,Mv,Mxi,alpha+1,beta+0,0)+&
                    a(3)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+1,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+2,beta+0,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+2,0)+&
                    0.5*a(4)*Moment_uvxi(Mu,Mv,Mxi,alpha+0,beta+0,2)
    end function Moment_auvxi
    
    
    subroutine W2BMultiscaleFluxBoundary(bc,face,cell,idx,side)
        real(KREAL), intent(in)                         :: bc(4) !Primary variables at boundary
        type(CellInterface), intent(inout)              :: face
        type(CellCenter), intent(in)                    :: cell
        integer(KINT), intent(in)                       :: idx,side
        integer(KINT)                                   :: k !Flattened velocity-space index
        real(KREAL)                                     :: rot
        real(KREAL), allocatable, dimension(:)          :: vn_vec,vt_vec !Pure 1D local-frame velocities
        real(KREAL), allocatable, dimension(:)          :: h_vec,b_vec !Pure 1D reconstructed wall distributions
        real(KREAL), allocatable, dimension(:)          :: H_g_vec,B_g_vec !Pure 1D interface Maxwellian arrays
        real(KREAL), allocatable, dimension(:)          :: H_w_vec,B_w_vec !Pure 1D wall Maxwellian arrays
        real(KREAL), allocatable, dimension(:)          :: H_plus_vec,B_plus_vec !Pure 1D Shakhov correction arrays
        real(KREAL)                                     :: qf(2)
        real(KREAL), allocatable, dimension(:)          :: delta_vec !Pure 1D wall upwind selector
        real(KREAL)                                     :: flux(4)
        real(KREAL)                                     :: flux_wall(4) !Accumulated wall macro-flux correction from Nvel loop
        real(KREAL)                                     :: prim(4)
        real(KREAL)                                     :: prim_g(4) !Primary variables of g_{i+1/2}
        real(KREAL)                                     :: prim_shape(4)
        real(KREAL)                                     :: prim_w(4) !boundary condition in local frame
        real(KREAL)                                     :: prim_g0(4),conVars_g0(4),conVars_g(4),conVars_g1(4)
        real(KREAL)                                     :: incidence,reflection
        !real(KREAL)                                     :: tau,T1,T4,T5
        real(KREAL)                                     :: tau,Mt(5)
        real(KREAL)                                     :: wk,vnK,vtK,v2 !Reusable velocity-space scalars for Nvel loops
        real(KREAL)                                     :: deltaK,inflowK !Wall reflected/incident selectors for Nvel loops
        real(KREAL)                                     :: H_gK,B_gK,H_wK,B_wK !Per-velocity Maxwellian values for micro flux
        real(KREAL)                                     :: H_plusK,B_plusK !Per-velocity Shakhov correction values for micro flux
        real(KREAL)                                     :: hK,bK !Per-velocity reconstructed distribution values
        real(KREAL)                                     :: shnK,sbnK,shtK,sbtK !Per-velocity distribution slopes
        real(KREAL)                                     :: slope_hK,slope_bK !Per-velocity normal/tangential slope combinations
        real(KREAL)                                     :: wallExpK !Per-velocity Maxwellian kernel used by wall reflection
        real(KREAL)                                     :: Ha_gnK,Ha_gtK,Ha_TK !Per-velocity h-space expansion terms
        real(KREAL)                                     :: Ba_gnK,Ba_gtK,Ba_TK !Per-velocity b-space expansion terms
        real(KREAL)                                     :: dt   !local time step

        real(KREAL),parameter :: PI=4.0_KREAL*atan(1.0_KREAL) ! Only this new wall kernel.

        ! additional variables for slope calculation
        real(KREAL), allocatable, dimension(:)          :: shn_vec,sbn_vec !Pure 1D normal distribution slopes
        !real(KREAL), allocatable, dimension(:,:)        :: Ha_gn, Ba_gn
        !real(KREAL), allocatable, dimension(:,:)        :: Ha_gt, Ba_gt 
        !real(KREAL), allocatable, dimension(:,:)        :: Ha_T, Ba_T  
        real(KREAL)                                     :: a_gn(4), a_gt(4), a_T(4)    
        real(KREAL)                                     :: sw(4), sw_T(4)
        real(KREAL)                                     :: Mau(4), Mbv(4), Mau0(4), MauT(4)
        real(KREAL)                                     :: Mu(0:MNUM),MuL(0:MNUM),MuR(0:MNUM)
        real(KREAL)                                     :: Mv(0:MTUM),Mxi(0:2)

        ! tangential slope variables
        real(KREAL), allocatable, dimension(:)          :: sht_vec,sbt_vec !Pure 1D tangential distribution slopes
        !limiter of tangential slope
        !real(KREAL)                                     :: tan_len, d_rho_t, limiter_phi
        !real(KREAL), parameter                          :: EPS = 1.0e-10
        !real(KREAL), parameter                          :: LIMIT_THRESHOLD = 0.5  ! allow max 50% change in tangential direction



        !--------------------------------------------------
        !prepare
        !--------------------------------------------------
        !Allocate pure 1D velocity-space temporaries for multiscale wall flux loops.
        allocate(vn_vec(Nvel))
        allocate(vt_vec(Nvel))
        allocate(delta_vec(Nvel))
        allocate(h_vec(Nvel))
        allocate(b_vec(Nvel))
        allocate(H_plus_vec(Nvel))
        allocate(B_plus_vec(Nvel))
        allocate(H_g_vec(Nvel))
        allocate(B_g_vec(Nvel))
        allocate(H_w_vec(Nvel))
        allocate(B_w_vec(Nvel))
        allocate(shn_vec(Nvel))
        allocate(sbn_vec(Nvel))
        allocate(sht_vec(Nvel))
        allocate(sbt_vec(Nvel))

        !Time step
        dt = cell%dt_local
        
        !Frame rotation
        rot = -side

        !Convert velocity space to the wall frame and build the wall upwind selector.
        do k=1,Nvel
            vn_vec(k) = uVel(k)*face%cosx+vVel(k)*face%cosy
            vt_vec(k) =-uVel(k)*face%cosy+vVel(k)*face%cosx
            delta_vec(k) = (sign(UP,vn_vec(k))*rot+1.0_KREAL)/2.0_KREAL
        enddo

        !--------------------------------------------------
        !Reconstruct non-equlibrium distribution at interface
        !--------------------------------------------------
        !Load normal and tangential slopes through 1D velocity views.
        do k=1,Nvel
            shn_vec(k) = cell%sh_vec(k,idx)
            sbn_vec(k) = cell%sb_vec(k,idx)

            if (idx == 1) then ! normal is x, tangential is y
                sht_vec(k) = cell%sh_vec(k,2)
                sbt_vec(k) = cell%sb_vec(k,2)
                !tan_len = cell%length(2) ! tangential length
            else ! normal is y, tangential is x
                sht_vec(k) = -cell%sh_vec(k,1)
                sbt_vec(k) = -cell%sb_vec(k,1)
                !tan_len = cell%length(1) ! tangential length
            end if
        enddo

        !Limiter
        ! calculaate the density change due to tangential slope
        !d_rho_t = sum(weight*sht) * (0.5 * tan_len)
        
        ! If the density change > 50% of the center density, calculate a reduction factor phi
        ! prim(1) is the density of the current cell
        !prim = GetPrimary(cell%conVars) ! temporarily get for comparison
        !
        !if (abs(d_rho_t) > LIMIT_THRESHOLD * prim(1)) then
        !    limiter_phi = (LIMIT_THRESHOLD * prim(1)) / (abs(d_rho_t) + EPS)
        !else
        !    limiter_phi = 1.0
        !end if
        
        !sht = sht * limiter_phi
        !sbt = sbt * limiter_phi
        
        ! if use tangential slope sht,sbt, all Mt(5) need change vn*sh to (vn*sh + vt*sht),vn*sb to (vn*sb + vt*sbt)

        !Reconstruct incident non-equilibrium distribution at the wall.
        do k=1,Nvel
            h_vec(k) = cell%h_vec(k)-rot*0.5_KREAL*cell%length(idx)*shn_vec(k)
            b_vec(k) = cell%b_vec(k)-rot*0.5_KREAL*cell%length(idx)*sbn_vec(k)
        enddo

        !Calculate primary variable of g_{i+1/2} in local frame
        prim_w = LocalFrame(bc,face%cosx,face%cosy)
        prim_g = prim_w
        prim = LocalFrame(GetPrimary(cell%conVars),face%cosx,face%cosy)
        prim_g(1) = prim(1)/prim(4)*prim_w(4) !keep consistent pressure inside cell

        !--------------------------------------------------
        !Calculate collision time and some time integration terms
        !--------------------------------------------------
        tau = GetTau(prim_g)

        Mt(4) = tau*(1.0-exp(-dt/tau))
        Mt(5) = -tau*dt*exp(-dt/tau)+tau*Mt(4)
        Mt(1) = dt-Mt(4)
        Mt(2) = -tau*Mt(1)+Mt(5) 
        Mt(3) = 0.5*dt**2-tau*Mt(1)

        !--------------------------------------------------
        !1. Calculate Maxwellian distribution and Shakhov part at g_{i+1/2}
        !--------------------------------------------------
        if(UnifiedPhysicalOperators)then
            call PrepareNativeFaceEquilibrium(H_g_vec,B_g_vec,vn_vec,vt_vec,prim_g)
        else
        do k=1,Nvel
            H_g_vec(k)=prim_g(1)*(prim_g(4)/PI)*exp(-prim_g(4)* &
                ((vn_vec(k)-prim_g(2))**2+(vt_vec(k)-prim_g(3))**2))
            B_g_vec(k)=H_g_vec(k)*real(CK,KREAL)/(2.0_KREAL*prim_g(4))
        enddo
        endif

        qf = GetHeatFlux_LocalVec(h_vec,b_vec,vn_vec,vt_vec,prim_g)
        call ShakhovPart_LocalVec(H_g_vec,B_g_vec,vn_vec,vt_vec,qf,prim_g,H_plus_vec,B_plus_vec)

        !--------------------------------------------------
        !2. Spatial slope a_g: derived from local gradients
        !--------------------------------------------------
        !sw = LocalFrame(cell%sw(:,idx), face%cosx, face%cosy) 
        !a_g  = MicroSlope(prim_g, sw)
        
        !Alternative: derived from difference between cell center and interface
        conVars_g = GetConserved(prim_g)
        conVars_g1 = LocalFrame(cell%conVars, face%cosx, face%cosy)
        sw = (conVars_g - conVars_g1) / (side * 0.5 * cell%length(idx))
        a_gn  = MicroSlope(prim_g, sw)
        
        !if (idx == 1) then
        !    sw_t = LocalFrame(cell%sw(:,2), face%cosx, face%cosy)
        !else
        !    sw_t = LocalFrame(-cell%sw(:,1), face%cosx, face%cosy) 
        !end if
        
        !Tangential macro slope obtained by integrating reconstructed VDF slopes.
        sw_t = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            vnK = vn_vec(k)
            vtK = vt_vec(k)
            v2  = vnK*vnK+vtK*vtK

            sw_t(1) = sw_t(1)+wk*sht_vec(k)
            sw_t(2) = sw_t(2)+wk*vnK*sht_vec(k)
            sw_t(3) = sw_t(3)+wk*vtK*sht_vec(k)
            sw_t(4) = sw_t(4)+0.5_KREAL*wk*(v2*sht_vec(k)+sbt_vec(k))
        enddo

        a_gt = MicroSlope(prim_g, sw_t)
        !--------------------------------------------------
        !3. Temporal slope a_T
        !--------------------------------------------------
        call CalcMoment(prim_g, Mu, Mv, Mxi, MuL, MuR)
        Mau = Moment_auvxi(a_gn, Mu, Mv, Mxi, 1, 0) 
        Mbv = Moment_auvxi(a_gt, Mu, Mv, Mxi, 0, 1)
        sw = -prim_g(1)*(Mau + Mbv) 
        a_T = MicroSlope(prim_g, sw)
        
        !--------------------------------------------------
        !4. Prepare expansion terms for Flux calculation
        !--------------------------------------------------
        !Ha_gn = (a_gn(1) + a_gn(2)*vn + a_gn(3)*vt)*H_g + 0.5*a_gn(4)*((vn**2+vt**2)*H_g + B_g)
        !Ba_gn = (a_gn(1) + a_gn(2)*vn + a_gn(3)*vt)*B_g + 0.5*a_gn(4)*((vn**2+vt**2)*B_g + Mxi(2)*H_g)
        !Ha_gt = (a_gt(1) + a_gt(2)*vn + a_gt(3)*vt)*H_g + 0.5*a_gt(4)*((vn**2+vt**2)*H_g + B_g)
        !Ba_gt = (a_gt(1) + a_gt(2)*vn + a_gt(3)*vt)*B_g + 0.5*a_gt(4)*((vn**2+vt**2)*B_g + Mxi(2)*H_g)
        !Ha_T  = (a_T(1) + a_T(2)*vn + a_T(3)*vt)*H_g + 0.5*a_T(4)*((vn**2+vt**2)*H_g + B_g)
        !Ba_T  = (a_T(1) + a_T(2)*vn + a_T(3)*vt)*B_g + 0.5*a_T(4)*((vn**2+vt**2)*B_g + Mxi(2)*H_g)
        

        if (side == Lside) then
            Mu = MuR
        else
            Mu = MuL
        end if

        Mau0 = Moment_uvxi(Mu,Mv,Mxi,1,0,0) !<u*\psi>
        Mau = Moment_auvxi(a_gn,Mu,Mv,Mxi,2,0) !<a*u^2*\psi>
        Mbv  = Moment_auvxi(a_gt,Mu,Mv,Mxi,1,1) !<b*u*v*\psi>
        MauT = Moment_auvxi(a_T,Mu,Mv,Mxi,1,0) !<A*u*\psi>
        
        ! flux of Maxwellian part
        flux = Mt(1)*prim_g(1)*Mau0+Mt(2)*prim_g(1)*(Mau+Mbv)+Mt(3)*prim_g(1)*MauT
        !--------------------------------------------------
        !Calculate wall density and Maxwellian distribution
        !--------------------------------------------------
        
        !incidence = Mt(1)*sum(weight*vn*(H_g+H_plus)*(1-delta)) + &
        !            Mt(2)*sum(weight*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + & 
        !            Mt(3)*sum(weight*vn*Ha_T*(1-delta)) + &      
        !            Mt(4)*sum(weight*vn*h*(1-delta)) - &
        !            Mt(5)*sum(weight*vn*(vn*shn + vt*sht)*(1-delta))
        
        ! Unit-density outgoing shapes use exactly the same low-order preparation.
        ! Wall density still comes only from the original zero-mass-flux equation.
        if(UnifiedPhysicalOperators)then
            prim_shape=prim_w;prim_shape(1)=1.0_KREAL
            call PrepareNativeFaceEquilibrium(H_w_vec,B_w_vec,vn_vec,vt_vec,prim_shape)
        endif
        !Compute wall incidence/reflection using 1D velocity views.
        incidence  = 0.0_KREAL
        reflection = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK
            slope_hK = vnK*shn_vec(k)+vtK*sht_vec(k)
            wallExpK = exp(-prim_w(4)*((vnK-prim_w(2))*(vnK-prim_w(2))+(vtK-prim_w(3))*(vtK-prim_w(3))))

            Ha_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*H_g_vec(k) &
                +0.5_KREAL*a_gn(4)*((vnK*vnK+vtK*vtK)*H_g_vec(k)+B_g_vec(k))
            Ha_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*H_g_vec(k) &
                +0.5_KREAL*a_gt(4)*((vnK*vnK+vtK*vtK)*H_g_vec(k)+B_g_vec(k))
            Ha_TK = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*H_g_vec(k) &
                +0.5_KREAL*a_T(4)*((vnK*vnK+vtK*vtK)*H_g_vec(k)+B_g_vec(k))
            incidence = incidence+wk*vnK*inflowK*(Mt(1)*(H_g_vec(k)+H_plus_vec(k)) &
                +Mt(2)*(vnK*Ha_gnK+vtK*Ha_gtK)+Mt(3)*Ha_TK+Mt(4)*h_vec(k)-Mt(5)*slope_hK)
            if(UnifiedPhysicalOperators)wallExpK=H_w_vec(k)
            reflection = reflection+wk*vnK*wallExpK*deltaK
        enddo
        if(UnifiedPhysicalOperators)then
            reflection=dt*reflection
        else
            reflection=dt*(prim_w(4)/PI)*reflection
        endif
        prim_w(1) = -incidence/reflection

        if(UnifiedPhysicalOperators)then
            H_w_vec=prim_w(1)*H_w_vec
            B_w_vec=prim_w(1)*B_w_vec
        else
        do k=1,Nvel
            H_w_vec(k)=prim_w(1)*(prim_w(4)/PI)*exp(-prim_w(4)* &
                ((vn_vec(k)-prim_w(2))**2+(vt_vec(k)-prim_w(3))**2))
            B_w_vec(k)=H_w_vec(k)*real(CK,KREAL)/(2.0_KREAL*prim_w(4))
        enddo
        endif

        !--------------------------------------------------
        ! Calculate Fluxes 
        !--------------------------------------------------
        ! 1. Mass Flux
        !flux(1) = dt*sum(weight*vn*H_w*delta) + &
        !          Mt(1)*sum(weight*vn*(H_g+H_plus)*(1-delta)) + &
        !          Mt(2)*sum(weight*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta))+ &   
        !          Mt(3)*sum(weight*vn*Ha_T*(1-delta)) + &        
        !          Mt(4)*sum(weight*vn*h*(1-delta)) - &
        !          Mt(5)*sum(weight*vn*(vn*shn + vt*sht)*(1-delta)) 
        !Accumulate reflected and incident macro flux corrections through 1D velocity views.
        flux_wall = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK

            H_wK     = H_w_vec(k)
            B_wK     = B_w_vec(k)
            H_plusK  = H_plus_vec(k)
            B_plusK  = B_plus_vec(k)
            hK       = h_vec(k)
            bK       = b_vec(k)
            slope_hK = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_bK = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            flux_wall(1) = flux_wall(1)+wk*vnK*(dt*H_wK*deltaK+Mt(1)*H_plusK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            flux_wall(2) = flux_wall(2)+wk*vnK*vnK*(dt*H_wK*deltaK+Mt(1)*H_plusK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            flux_wall(3) = flux_wall(3)+wk*vnK*vtK*(dt*H_wK*deltaK+Mt(1)*H_plusK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            flux_wall(4) = flux_wall(4)+0.5_KREAL*wk*vnK*(dt*(v2*H_wK+B_wK)*deltaK+ &
                Mt(1)*(v2*H_plusK+B_plusK)*inflowK+Mt(4)*(v2*hK+bK)*inflowK- &
                Mt(5)*(v2*slope_hK+slope_bK)*inflowK)
        enddo
        flux = flux+flux_wall

        ! 2. Normal Momentum Flux
        !flux(2) = dt*sum(weight*vn*vn*H_w*delta) + &
        !          Mt(1)*sum(weight*vn*vn*(H_g+H_plus)*(1-delta)) + &
        !          Mt(2)*sum(weight*vn*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + & 
        !          Mt(3)*sum(weight*vn*vn*Ha_T*(1-delta)) + &     
        !          Mt(4)*sum(weight*vn*vn*h*(1-delta)) - &
        !          Mt(5)*sum(weight*vn*vn*(vn*shn + vt*sht)*(1-delta)) 
        ! The Nvel loop above replaces the four separate macro-flux sum expressions.
        ! 3. Tangential Momentum Flux
        !flux(3) = dt*sum(weight*vn*vt*H_w*delta) + &
        !          Mt(1)*sum(weight*vn*vt*(H_g+H_plus)*(1-delta)) + &
        !          Mt(2)*sum(weight*vn*vt*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + & 
        !          Mt(3)*sum(weight*vn*vt*Ha_T*(1-delta)) + &     
        !          Mt(4)*sum(weight*vt*vn*h*(1-delta)) - &
        !          Mt(5)*sum(weight*vt*vn*(vn*shn + vt*sht)*(1-delta)) 
        ! The Nvel loop above replaces the tangential momentum flux sum.
        ! 4. Energy Flux
        !flux(4) = dt*0.5*sum(weight*vn*((vn**2+vt**2)*H_w+B_w)*delta)+ &
        !          Mt(1)*0.5*sum(weight*vn*((vn**2+vt**2)*(H_g+H_plus)+(B_g+B_plus))*(1-delta))+ &
        !          Mt(2)*0.5*sum(weight*vn*((vn**2+vt**2)*(vn*Ha_gn+vt*Ha_gt)+(vn*Ba_gn+vt*Ba_gt))*(1-delta))+ & 
        !          Mt(3)*0.5*sum(weight*vn*((vn**2+vt**2)*Ha_T+Ba_T)*(1-delta))+ &      
        !          Mt(4)*0.5*sum(weight*vn*((vn**2+vt**2)*h+b)*(1-delta)) - &
        !          Mt(5)*0.5*sum(weight*vn*((vn**2+vt**2)*(vn*shn + vt*sht) + (vn*sbn + vt*sbt))*(1-delta))
        ! The Nvel loop above replaces the energy flux sum.



        !flux_h = dt*vn*H_w*delta + &
        !         Mt(1)*vn*(H_g+H_plus)*(1-delta) + &
        !         Mt(2)*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta) + &  
        !         Mt(3)*vn*Ha_T*(1-delta) + &       
        !         Mt(4)*vn*h*(1-delta) - &
        !         Mt(5)*vn*(vn*shn + vt*sht)*(1-delta) 
        
        !Store multiscale wall micro flux directly through 1D face views.
        do k=1,Nvel
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK

            H_gK     = H_g_vec(k)
            B_gK     = B_g_vec(k)
            H_wK     = H_w_vec(k)
            B_wK     = B_w_vec(k)
            H_plusK  = H_plus_vec(k)
            B_plusK  = B_plus_vec(k)
            hK       = h_vec(k)
            bK       = b_vec(k)
            shnK     = shn_vec(k)
            sbnK     = sbn_vec(k)
            shtK     = sht_vec(k)
            sbtK     = sbt_vec(k)

            Ha_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*H_gK+0.5_KREAL*a_gn(4)*(v2*H_gK+B_gK)
            Ha_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*H_gK+0.5_KREAL*a_gt(4)*(v2*H_gK+B_gK)
            Ha_TK  = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*H_gK+0.5_KREAL*a_T(4)*(v2*H_gK+B_gK)

            Ba_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*B_gK+0.5_KREAL*a_gn(4)*(v2*B_gK+Mxi(2)*H_gK)
            Ba_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*B_gK+0.5_KREAL*a_gt(4)*(v2*B_gK+Mxi(2)*H_gK)
            Ba_TK  = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*B_gK+0.5_KREAL*a_T(4)*(v2*B_gK+Mxi(2)*H_gK)

            if(side==Lside) then
                face%L_flux_h_vec(k) = face%length*(dt*vnK*H_wK*deltaK + &
                    Mt(1)*vnK*(H_gK+H_plusK)*inflowK + Mt(2)*vnK*vnK*Ha_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ha_gtK*inflowK + Mt(3)*vnK*Ha_TK*inflowK + &
                    Mt(4)*vnK*hK*inflowK - Mt(5)*vnK*(vnK*shnK+vtK*shtK)*inflowK)

                face%L_flux_b_vec(k) = face%length*(dt*vnK*B_wK*deltaK + &
                    Mt(1)*vnK*(B_gK+B_plusK)*inflowK + Mt(2)*vnK*vnK*Ba_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ba_gtK*inflowK + Mt(3)*vnK*Ba_TK*inflowK + &
                    Mt(4)*vnK*bK*inflowK - Mt(5)*vnK*(vnK*sbnK+vtK*sbtK)*inflowK)
            elseif(side==Rside) then
                face%R_flux_h_vec(k) = face%length*(dt*vnK*H_wK*deltaK + &
                    Mt(1)*vnK*(H_gK+H_plusK)*inflowK + Mt(2)*vnK*vnK*Ha_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ha_gtK*inflowK + Mt(3)*vnK*Ha_TK*inflowK + &
                    Mt(4)*vnK*hK*inflowK - Mt(5)*vnK*(vnK*shnK+vtK*shtK)*inflowK)

                face%R_flux_b_vec(k) = face%length*(dt*vnK*B_wK*deltaK + &
                    Mt(1)*vnK*(B_gK+B_plusK)*inflowK + Mt(2)*vnK*vnK*Ba_gnK*inflowK + &
                    Mt(2)*vnK*vtK*Ba_gtK*inflowK + Mt(3)*vnK*Ba_TK*inflowK + &
                    Mt(4)*vnK*bK*inflowK - Mt(5)*vnK*(vnK*sbnK+vtK*sbtK)*inflowK)
            endif
        enddo
        !flux_b = dt*vn*B_w*delta + &
        !         Mt(1)*vn*(B_g+B_plus)*(1-delta) + &
        !         Mt(2)*vn*(vn*Ba_gn + vt*Ba_gt)*(1-delta) + &   
        !         Mt(3)*vn*Ba_T*(1-delta) + &      
        !         Mt(4)*vn*b*(1-delta) - &
        !         Mt(5)*vn*(vn*sbn + vt*sbt)*(1-delta) 

        !--------------------------------------------------
        !Store surface FlowVars and Aerodynamics
        !--------------------------------------------------
        !Flow variables
        
        !conVars_g0(1) = dt*sum(weight*H_w*delta) + &
        !                Mt(1)*sum(weight*(H_g+H_plus)*(1-delta)) + &
        !                Mt(2)*sum(weight*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + &   
        !                Mt(3)*sum(weight*Ha_T*(1-delta)) + &       
        !                Mt(4)*sum(weight*h*(1-delta)) - &
        !                Mt(5)*sum(weight*(vn*shn + vt*sht)*(1-delta)) 

        !conVars_g0(2) = dt*sum(weight*vn*H_w*delta) + &
        !                Mt(1)*sum(weight*vn*(H_g+H_plus)*(1-delta)) + &
        !                Mt(2)*sum(weight*vn*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + &
        !                Mt(3)*sum(weight*vn*Ha_T*(1-delta)) + &    
        !                Mt(4)*sum(weight*vn*h*(1-delta)) - &
        !                Mt(5)*sum(weight*vn*(vn*shn + vt*sht)*(1-delta)) 

        !conVars_g0(3) = dt*sum(weight*vt*H_w*delta) + &
        !                Mt(1)*sum(weight*vt*(H_g+H_plus)*(1-delta)) + &
        !                Mt(2)*sum(weight*vt*(vn*Ha_gn + vt*Ha_gt)*(1-delta)) + &
        !                Mt(3)*sum(weight*vt*Ha_T*(1-delta)) + &    
        !                Mt(4)*sum(weight*vt*h*(1-delta)) - &
        !                Mt(5)*sum(weight*vt*(vn*shn + vt*sht)*(1-delta)) 

        !conVars_g0(4) = dt*0.5*sum(weight*((vn**2+vt**2)*H_w+B_w)*delta)+ &
        !                Mt(1)*0.5*sum(weight*((vn**2+vt**2)*(H_g+H_plus)+(B_g+B_plus))*(1-delta))+ &
        !                Mt(2)*0.5*sum(weight*((vn**2+vt**2)*(vn*Ha_gn + vt*Ha_gt) + (vn*Ba_gn + vt*Ba_gt))*(1-delta))+ &
        !                Mt(3)*0.5*sum(weight*((vn**2+vt**2)*Ha_T+Ba_T)*(1-delta))+ &      
        !                Mt(4)*0.5*sum(weight*((vn**2+vt**2)*h+b)*(1-delta)) - &
        !                Mt(5)*0.5*sum(weight*((vn**2+vt**2)*(vn*shn + vt*sht) + (vn*sbn + vt*sbt))*(1-delta))

        ! Macro from the exact same time-and-length integrated micro flux.
        flux=0.0_KREAL
        do k=1,Nvel
            if(side==Lside)then
                hK=face%L_flux_h_vec(k)/face%length
                bK=face%L_flux_b_vec(k)/face%length
            else
                hK=face%R_flux_h_vec(k)/face%length
                bK=face%R_flux_b_vec(k)/face%length
            endif
            flux=flux+wVel(k)*[hK,vn_vec(k)*hK,vt_vec(k)*hK, &
                0.5_KREAL*((vn_vec(k)**2+vt_vec(k)**2)*hK+bK)]
        enddo
        !Build wall output conserved variables through a single 1D velocity loop.
        conVars_g0 = 0.0_KREAL
        do k=1,Nvel
            wk      = wVel(k)
            vnK     = vn_vec(k)
            vtK     = vt_vec(k)
            v2      = vnK*vnK+vtK*vtK
            deltaK  = delta_vec(k)
            inflowK = 1.0_KREAL-deltaK

            H_gK     = H_g_vec(k)
            B_gK     = B_g_vec(k)
            H_wK     = H_w_vec(k)
            B_wK     = B_w_vec(k)
            H_plusK  = H_plus_vec(k)
            B_plusK  = B_plus_vec(k)
            hK       = h_vec(k)
            bK       = b_vec(k)
            slope_hK = vnK*shn_vec(k)+vtK*sht_vec(k)
            slope_bK = vnK*sbn_vec(k)+vtK*sbt_vec(k)

            Ha_gnK = (a_gn(1)+a_gn(2)*vnK+a_gn(3)*vtK)*H_gK+0.5_KREAL*a_gn(4)*(v2*H_gK+B_gK)
            Ha_gtK = (a_gt(1)+a_gt(2)*vnK+a_gt(3)*vtK)*H_gK+0.5_KREAL*a_gt(4)*(v2*H_gK+B_gK)
            Ha_TK  = (a_T(1)+a_T(2)*vnK+a_T(3)*vtK)*H_gK+0.5_KREAL*a_T(4)*(v2*H_gK+B_gK)

            conVars_g0(1) = conVars_g0(1)+wk*(dt*H_wK*deltaK+Mt(1)*(H_gK+H_plusK)*inflowK+ &
                Mt(2)*vnK*Ha_gnK*inflowK+Mt(2)*vtK*Ha_gtK*inflowK+Mt(3)*Ha_TK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            conVars_g0(2) = conVars_g0(2)+wk*vnK*(dt*H_wK*deltaK+Mt(1)*(H_gK+H_plusK)*inflowK+ &
                Mt(2)*vnK*Ha_gnK*inflowK+Mt(2)*vtK*Ha_gtK*inflowK+Mt(3)*Ha_TK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)
            conVars_g0(3) = conVars_g0(3)+wk*vtK*(dt*H_wK*deltaK+Mt(1)*(H_gK+H_plusK)*inflowK+ &
                Mt(2)*vnK*Ha_gnK*inflowK+Mt(2)*vtK*Ha_gtK*inflowK+Mt(3)*Ha_TK*inflowK+ &
                Mt(4)*hK*inflowK-Mt(5)*slope_hK*inflowK)

            !Keep the historical conVars_g0(4) algebra exactly; do not substitute the commented formula here.
            conVars_g0(4) = conVars_g0(4)+0.5_KREAL*wk*(dt*(v2*H_wK+B_wK)*deltaK+ &
                Mt(1)*(v2*(H_gK+H_plusK)+(B_gK+B_plusK))*inflowK+ &
                Mt(2)*vnK*(v2*Ha_gnK+Ha_gtK*vtK)*inflowK+ &
                Mt(2)*vtK*(v2*Ha_gnK+Ha_gtK*vtK)*inflowK+Mt(3)*v2*Ha_TK*inflowK+ &
                Mt(4)*(v2*hK+bK)*inflowK-Mt(5)*(v2*slope_hK+slope_bK)*inflowK)
        enddo
        

        prim_g0       = GlobalFrame(GetPrimary(conVars_g0/dt),face%cosx,face%cosy)
        face%FlowVars(1:3)   = prim_g0(1:3)                         !density,u,v
        face%FlowVars(4)     = 1.0/prim_g0(4)                       !temperature
        face%FlowVars(5)     = 0.5*prim_g0(1)/prim_g0(4)            !pressure
        !Aerodynamics (Notice: here assume Unw==0, i.e. prim_w(2)==0)
        face%Aerodynamics(1) = flux(2)/dt                           !pnn: pressure
        face%Aerodynamics(2) = flux(3)/dt                           !pnt: shear stress
        face%Aerodynamics(3) = flux(4)/dt-prim_w(3)*flux(3)/dt      !qn:  heat flux

        !--------------------------------------------------
        !Final flux
        !--------------------------------------------------
        !Convert to global frame
        
        flux = GlobalFrame(flux,face%cosx,face%cosy)

        !Total flux
        flux = face%length*flux

        if(side==Lside) then
            face%L_flux     = flux
        elseif(side==Rside) then
            face%R_flux     = flux
        endif
        
        !--------------------------------------------------
        !Aftermath
        !--------------------------------------------------
        !Release pure 1D velocity-space temporaries allocated above.
        deallocate(vn_vec)
        deallocate(vt_vec)
        deallocate(delta_vec)
        deallocate(h_vec)
        deallocate(b_vec)
        deallocate(H_g_vec)
        deallocate(B_g_vec)
        deallocate(H_w_vec)
        deallocate(B_w_vec)
        deallocate(H_plus_vec)
        deallocate(B_plus_vec)
        deallocate(shn_vec)
        deallocate(sbn_vec)
        !deallocate(Ha_gn); deallocate(Ba_gn)
        !deallocate(Ha_gt); deallocate(Ba_gt)
        !deallocate(Ha_T); deallocate(Ba_T)
        deallocate(sht_vec)
        deallocate(sbt_vec)
    end subroutine W2BMultiscaleFluxBoundary

end module mod_Evolution

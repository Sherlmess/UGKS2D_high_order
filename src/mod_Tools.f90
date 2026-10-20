!--------------------------------------------------
!>Define some commonly used functions/subroutines
!--------------------------------------------------
module Tools
    use flow_data
    use mpi_para
    use mpi, only: MPI_COMM_WORLD, MPI_Abort
    use, intrinsic :: ieee_arithmetic
    implicit none

    contains
    
    !--------------------------------------------------
    !>Obtain discretized Maxwellian distribution on the global flattened velocity grid.
    !>@param[out] h_vec,b_vec :1D distribution functions
    !>@param[in]  prim        :primary variables
    !--------------------------------------------------
    subroutine DiscreteMaxwell_Vec(h_vec,b_vec,prim)
        real(KREAL), intent(out), dimension(:)             :: h_vec,b_vec !Flattened Maxwellian distributions
        real(KREAL), intent(in)                            :: prim(4)
        integer(KINT)                                      :: k !Flattened velocity-space index
        real(KREAL)                                        :: prefactor,b_factor,dvn,dvt,v2,hK !Reusable Nvel-loop scalars

        !Evaluate Maxwellian values directly on uVel/vVel without Nu-by-Nv temporary arrays.
        prefactor = prim(1)*(prim(4)/PI)
        b_factor = CK/(2.0_KREAL*prim(4))
        do k=1,Nvel
            dvn = uVel(k)-prim(2)
            dvt = vVel(k)-prim(3)
            v2  = dvn*dvn+dvt*dvt
            hK  = prefactor*exp(-prim(4)*v2)

            h_vec(k) = hK
            b_vec(k) = hK*b_factor
        enddo
    end subroutine DiscreteMaxwell_Vec

    !--------------------------------------------------
    !>Obtain discretized Maxwellian distribution on a supplied local velocity frame.
    !>@param[out] h_vec,b_vec   :1D distribution functions
    !>@param[in]  vn_vec,vt_vec :1D normal/tangential velocity coordinates
    !>@param[in]  prim          :primary variables in the local frame
    !--------------------------------------------------
    subroutine DiscreteMaxwell_LocalVec(h_vec,b_vec,vn_vec,vt_vec,prim)
        real(KREAL), intent(out), dimension(:)             :: h_vec,b_vec !Flattened Maxwellian distributions
        real(KREAL), intent(in), dimension(:)              :: vn_vec,vt_vec !Flattened local velocity coordinates
        real(KREAL), intent(in)                            :: prim(4)
        integer(KINT)                                      :: k !Flattened velocity-space index
        real(KREAL)                                        :: prefactor,b_factor,dvn,dvt,v2,hK !Reusable Nvel-loop scalars

        !Evaluate Maxwellian values directly on local-frame velocity arrays.
        prefactor = prim(1)*(prim(4)/PI)
        b_factor = CK/(2.0_KREAL*prim(4))
        do k=1,Nvel
            dvn = vn_vec(k)-prim(2)
            dvt = vt_vec(k)-prim(3)
            v2  = dvn*dvn+dvt*dvt
            hK  = prefactor*exp(-prim(4)*v2)

            h_vec(k) = hK
            b_vec(k) = hK*b_factor
        enddo
    end subroutine DiscreteMaxwell_LocalVec

    !--------------------------------------------------
    !>Calculate the Shakhov part H^+, B^+ on the global flattened velocity grid.
    !>@param[in]  H_vec,B_vec           :1D Maxwellian distribution functions
    !>@param[in]  qf                    :heat flux
    !>@param[in]  prim                  :primary variables
    !>@param[out] H_plus_vec,B_plus_vec :1D Shakhov correction
    !--------------------------------------------------
    subroutine ShakhovPart_Vec(H_vec,B_vec,qf,prim,H_plus_vec,B_plus_vec)
        real(KREAL), intent(in), dimension(:)              :: H_vec,B_vec !Flattened Maxwellian inputs
        real(KREAL), intent(in)                            :: qf(2)
        real(KREAL), intent(in)                            :: prim(4)
        real(KREAL), intent(out), dimension(:)             :: H_plus_vec,B_plus_vec !Flattened Shakhov outputs
        integer(KINT)                                      :: k !Flattened velocity-space index
        real(KREAL)                                        :: coeff,dvn,dvt,v2,qdot,shapeH,shapeB !Reusable Nvel-loop scalars

        !Evaluate Shakhov corrections directly on uVel/vVel without Nu-by-Nv temporary arrays.
        coeff = 0.8_KREAL*(1.0_KREAL-PR)*prim(4)*prim(4)/prim(1)
        if(CK>0) then
            do k=1,Nvel
                dvn = uVel(k)-prim(2)
                dvt = vVel(k)-prim(3)
                v2 = dvn*dvn+dvt*dvt
                qdot = dvn*qf(1)+dvt*qf(2)
                shapeH = 2.0_KREAL*prim(4)*v2-4.0_KREAL
                shapeB = shapeH+2.0_KREAL/CK

                H_plus_vec(k) = coeff*qdot*shapeH*H_vec(k)
                B_plus_vec(k) = coeff*qdot*shapeB*B_vec(k)
            enddo
        elseif(CK==0) then
            do k=1,Nvel
                dvn = uVel(k)-prim(2)
                dvt = vVel(k)-prim(3)
                v2 = dvn*dvn+dvt*dvt
                qdot = dvn*qf(1)+dvt*qf(2)
                shapeH = 2.0_KREAL*prim(4)*v2+CK-5.0_KREAL
                shapeB = 2.0_KREAL*prim(4)*v2+CK-3.0_KREAL

                H_plus_vec(k) = coeff*qdot*shapeH*H_vec(k)
                B_plus_vec(k) = coeff*qdot*shapeB*B_vec(k)
            enddo
        endif
    end subroutine ShakhovPart_Vec

    !--------------------------------------------------
    !>Calculate the Shakhov part H^+, B^+ on a supplied local velocity frame.
    !>@param[in]  H_vec,B_vec           :1D Maxwellian distribution functions
    !>@param[in]  vn_vec,vt_vec         :1D normal/tangential velocity coordinates
    !>@param[in]  qf                    :heat flux
    !>@param[in]  prim                  :primary variables in the local frame
    !>@param[out] H_plus_vec,B_plus_vec :1D Shakhov correction
    !--------------------------------------------------
    subroutine ShakhovPart_LocalVec(H_vec,B_vec,vn_vec,vt_vec,qf,prim,H_plus_vec,B_plus_vec)
        real(KREAL), intent(in), dimension(:)              :: H_vec,B_vec !Flattened Maxwellian inputs
        real(KREAL), intent(in), dimension(:)              :: vn_vec,vt_vec !Flattened local velocity coordinates
        real(KREAL), intent(in)                            :: qf(2)
        real(KREAL), intent(in)                            :: prim(4)
        real(KREAL), intent(out), dimension(:)             :: H_plus_vec,B_plus_vec !Flattened Shakhov outputs
        integer(KINT)                                      :: k !Flattened velocity-space index
        real(KREAL)                                        :: coeff,dvn,dvt,v2,qdot,shapeH,shapeB !Reusable Nvel-loop scalars

        !Evaluate Shakhov corrections directly on local-frame velocity arrays.
        coeff = 0.8_KREAL*(1.0_KREAL-PR)*prim(4)*prim(4)/prim(1)
        if(CK>0) then
            do k=1,Nvel
                dvn = vn_vec(k)-prim(2)
                dvt = vt_vec(k)-prim(3)
                v2 = dvn*dvn+dvt*dvt
                qdot = dvn*qf(1)+dvt*qf(2)
                shapeH = 2.0_KREAL*prim(4)*v2-4.0_KREAL
                shapeB = shapeH+2.0_KREAL/CK

                H_plus_vec(k) = coeff*qdot*shapeH*H_vec(k)
                B_plus_vec(k) = coeff*qdot*shapeB*B_vec(k)
            enddo
        elseif(CK==0) then
            do k=1,Nvel
                dvn = vn_vec(k)-prim(2)
                dvt = vt_vec(k)-prim(3)
                v2 = dvn*dvn+dvt*dvt
                qdot = dvn*qf(1)+dvt*qf(2)
                shapeH = 2.0_KREAL*prim(4)*v2+CK-5.0_KREAL
                shapeB = 2.0_KREAL*prim(4)*v2+CK-3.0_KREAL

                H_plus_vec(k) = coeff*qdot*shapeH*H_vec(k)
                B_plus_vec(k) = coeff*qdot*shapeB*B_vec(k)
            enddo
        endif
    end subroutine ShakhovPart_LocalVec
    
    !--------------------------------------------------
    !>Convert macro variables from global frame to local
    !>@param[in] w            :macro variables in global frame
    !>@param[in] cosx,cosy    :directional cosine
    !>@return    LocalFrame   :macro variables in local frame
    !--------------------------------------------------
    function LocalFrame(w,cosx,cosy)
        real(KREAL), intent(in)                         :: w(4)
        real(KREAL), intent(in)                         :: cosx,cosy
        real(KREAL)                                     :: LocalFrame(4)

        LocalFrame(1) = w(1)
        LocalFrame(2) = w(2)*cosx+w(3)*cosy
        LocalFrame(3) =-w(2)*cosy+w(3)*cosx
        LocalFrame(4) = w(4)
    end function LocalFrame

    !--------------------------------------------------
    !>Convert macro variables from local frame to global
    !>@param[in] w            :macro variables in local frame
    !>@param[in] cosx,cosy    :directional cosine
    !>@return    GlobalFrame  :macro variables in global frame
    !--------------------------------------------------
    function GlobalFrame(w,cosx,cosy)
        real(KREAL), intent(in)                         :: w(4)
        real(KREAL), intent(in)                         :: cosx,cosy
        real(KREAL)                                     :: GlobalFrame(4)

        GlobalFrame(1) = w(1)
        GlobalFrame(2) = w(2)*cosx-w(3)*cosy
        GlobalFrame(3) = w(2)*cosy+w(3)*cosx
        GlobalFrame(4) = w(4)
    end function GlobalFrame

    !--------------------------------------------------
    !>Convert primary variables to conservative variables in the equilibrium state
    !>@param[in] prim          :primary variables
    !>@return    GetConserved  :conservative variables
    !--------------------------------------------------
    function GetConserved(prim)
        real(KREAL), intent(in)                         :: prim(4) !Density, x-velocity, y-velocity, lambda=1/temperature
        real(KREAL)                                     :: GetConserved(4) !Density, x-momentum, y-momentum, total energy

        GetConserved(1) = prim(1)
        GetConserved(2) = prim(1)*prim(2)
        GetConserved(3) = prim(1)*prim(3)
        GetConserved(4) = 0.5*prim(1)/(prim(4)*(GAMMA-1.0))+0.5*prim(1)*(prim(2)**2+prim(3)**2)
    end function GetConserved

    !--------------------------------------------------
    !>Convert conservative variables to primary variables in the equilibrium state
    !>@param[in] w           :conservative variables
    !>@return    GetPrimary  :primary variables
    !--------------------------------------------------
    function GetPrimary(w)
        real(KREAL), intent(in)                         :: w(4) !Density, x-momentum, y-momentum, total energy
        real(KREAL)                                     :: GetPrimary(4) !Density, x-velocity, y-velocity, lambda=1/temperature
        real(KREAL)                                     :: denom
        integer                                         :: ierr_local

        if (.not.all(ieee_is_finite(w))) then
            print '(A,I6,1P,4E15.7)', 'Non-finite conservative vars entering GetPrimary (rank,w):', myID,w
            call MPI_Abort(MPI_COMM_WORLD, 101, ierr_local)
        end if

        if (w(1) <= 0._KREAL) then
            print '(A,I6,1P,4E15.7)', 'Non-positive density in GetPrimary (rank,w):', myID,w
            call MPI_Abort(MPI_COMM_WORLD, 102, ierr_local)
        end if

        denom = w(4)-0.5*(w(2)**2+w(3)**2)/w(1)
        if (.not.ieee_is_finite(denom) .or. denom <= 0._KREAL) then
            print '(A,I6,1P,5E15.7)', 'Invalid energy state in GetPrimary (rank,w,denom):', myID,w,denom
            call MPI_Abort(MPI_COMM_WORLD, 103, ierr_local)
        end if

        GetPrimary(1) = w(1)
        GetPrimary(2) = w(2)/w(1)
        GetPrimary(3) = w(3)/w(1)
        GetPrimary(4) = 0.5*w(1)/((GAMMA-1.0)*denom)
    end function GetPrimary

    !--------------------------------------------------
    !>obtain ratio of specific heat
    !>@param[in] ck        :internal degree of freedom
    !>@return    get_gamma :ratio of specific heat
    !--------------------------------------------------
    function get_gamma(ck)
        integer,intent(in) :: ck
        real(KREAL) :: get_gamma

        get_gamma = float(ck+4)/float(ck+2)
    end function get_gamma
    
    !--------------------------------------------------
    !>get the nondimensionalized viscosity coefficient at reference state
    !>@param[in] kn          :Knudsen number
    !>@param[in] alpha,omega :index related to HS/VHS/VSS model
    !>@return    get_mu_ref  :nondimensionalized viscosity coefficient
    !--------------------------------------------------
    function get_mu_ref(kn,alpha,omega)
        real(KREAL),intent(in) :: kn,alpha,omega
        real(KREAL) :: get_mu_ref

        get_mu_ref = 5*(alpha+1)*(alpha+2)*sqrt(PI)/(4*alpha*(5-2*omega)*(7-2*omega))*kn
    end function get_mu_ref
    
    !--------------------------------------------------
    !>Calculate collision time
    !>@param[in] prim    :primary variables
    !>@return    GetTau  :collision time
    !--------------------------------------------------
    function GetTau(prim)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL)                                     :: GetTau
        real(KREAL)                                     :: mu
        integer                                         :: ierr_local

        if (prim(1) <= 0._KREAL .or. prim(4) <= 0._KREAL .or. .not.all(ieee_is_finite(prim))) then
            print '(A,I6,1P,4E15.7)', 'Invalid primary vars entering GetTau (rank,prim):', myID,prim
            call MPI_Abort(MPI_COMM_WORLD, 104, ierr_local)
        end if

        if(trim(adjustl(PeriodicComparisonProfile))=='BGK_ConstTau')then
            GetTau=0.1_KREAL
            return
        endif

        if(flag_vislaw==PowerLaw) then
            mu = mu_ref*(1.0/prim(4))**omega
        elseif(flag_vislaw==Sutherland) then
            mu = mu_ref*sqrt(1.0/prim(4))*(1.0+Tcm)/(1.0+Tcm*prim(4))
        endif

        if (.not.ieee_is_finite(mu) .or. mu <= 0._KREAL) then
            print '(A,I6,1P,5E15.7)', 'Invalid viscosity in GetTau (rank,prim,mu):', myID,prim,mu
            call MPI_Abort(MPI_COMM_WORLD, 105, ierr_local)
        end if

        GetTau = mu/(0.5*prim(1)/prim(4))       !tau = mu/p

        if (.not.ieee_is_finite(GetTau) .or. GetTau <= 0._KREAL) then
            print '(A,I6,1P,5E15.7)', 'Invalid tau computed (rank,prim,tau):', myID,prim,GetTau
            call MPI_Abort(MPI_COMM_WORLD, 106, ierr_local)
        end if
    end function GetTau
    
    !--------------------------------------------------
    !>Obtain speed of sound
    !>@param[in] prim    :primary variables
    !>@return    GetSoundSpeed :speed of sound
    !--------------------------------------------------
    function GetSoundSpeed(prim)
        real(KREAL), intent(in)                         :: prim(4)
        real(KREAL)                                     :: GetSoundSpeed !Speed of sound
        integer                                         :: ierr_local

        if (.not.all(ieee_is_finite(prim)) .or. prim(4) <= 0._KREAL) then
            print '(A,I6,1P,4E15.7)', 'Invalid primary vars entering GetSoundSpeed (rank,prim):', myID,prim
            call MPI_Abort(MPI_COMM_WORLD, 107, ierr_local)
        end if

        GetSoundSpeed = sqrt(0.5*GAMMA/prim(4))
    end function GetSoundSpeed

    !--------------------------------------------------
    !>Get pressure from flattened velocity-space distributions
    !>@param[in] h_vec,b_vec    :1D distribution functions
    !>@param[in] prim           :primary variables
    !>@return    GetPressure_Vec:pressure
    !--------------------------------------------------
    function GetPressure_Vec(h_vec,b_vec,prim)
        real(KREAL), intent(in), dimension(:)             :: h_vec,b_vec !Flattened input distributions
        real(KREAL), intent(in)                           :: prim(4)
        real(KREAL)                                       :: GetPressure_Vec !pressure
        integer(KINT)                                     :: k !Flattened velocity-space index
        real(KREAL)                                       :: wk,dvn,dvt,v2,energy_sum,b_sum !Nvel-loop accumulators

        !Accumulate pressure moments directly from flattened velocity-space views.
        energy_sum = 0.0_KREAL
        b_sum = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            dvn = uVel(k)-prim(2)
            dvt = vVel(k)-prim(3)
            v2  = dvn*dvn+dvt*dvt

            energy_sum = energy_sum+wk*v2*h_vec(k)
            b_sum = b_sum+wk*b_vec(k)
        enddo

        GetPressure_Vec = (energy_sum+b_sum)/(CK+2.0_KREAL)
    end function GetPressure_Vec
    
    !--------------------------------------------------
    !>Calculate heat flux from flattened velocity-space distributions
    !>@param[in] h_vec,b_vec    :1D distribution functions
    !>@param[in] prim           :primary variables
    !>@return    GetHeatFlux_Vec:heat flux in x/y direction
    !--------------------------------------------------
    function GetHeatFlux_Vec(h_vec,b_vec,prim)
        real(KREAL), intent(in), dimension(:)             :: h_vec,b_vec !Flattened input distributions
        real(KREAL), intent(in)                           :: prim(4)
        real(KREAL)                                       :: GetHeatFlux_Vec(2) !Heat flux in x/y direction
        integer(KINT)                                     :: k !Flattened velocity-space index
        real(KREAL)                                       :: wk,dvn,dvt,v2 !Reusable Nvel-loop scalars

        !Accumulate heat flux directly from flattened velocity-space views.
        GetHeatFlux_Vec = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            dvn = uVel(k)-prim(2)
            dvt = vVel(k)-prim(3)
            v2  = dvn*dvn+dvt*dvt

            GetHeatFlux_Vec(1) = GetHeatFlux_Vec(1)+0.5_KREAL*wk*dvn*(v2*h_vec(k)+b_vec(k))
            GetHeatFlux_Vec(2) = GetHeatFlux_Vec(2)+0.5_KREAL*wk*dvt*(v2*h_vec(k)+b_vec(k))
        enddo
    end function GetHeatFlux_Vec

    !--------------------------------------------------
    !>Calculate heat flux from flattened distributions on a supplied local velocity frame.
    !>@param[in] h_vec,b_vec    :1D distribution functions
    !>@param[in] vn_vec,vt_vec  :1D normal/tangential velocity coordinates
    !>@param[in] prim           :primary variables in the local frame
    !>@return    GetHeatFlux_LocalVec:heat flux in normal/tangential direction
    !--------------------------------------------------
    function GetHeatFlux_LocalVec(h_vec,b_vec,vn_vec,vt_vec,prim)
        real(KREAL), intent(in), dimension(:)             :: h_vec,b_vec !Flattened input distributions
        real(KREAL), intent(in), dimension(:)             :: vn_vec,vt_vec !Flattened local velocity coordinates
        real(KREAL), intent(in)                           :: prim(4)
        real(KREAL)                                       :: GetHeatFlux_LocalVec(2) !Heat flux in normal/tangential direction
        integer(KINT)                                     :: k !Flattened velocity-space index
        real(KREAL)                                       :: wk,dvn,dvt,v2 !Reusable Nvel-loop scalars

        !Accumulate heat flux directly from flattened local-frame velocity views.
        GetHeatFlux_LocalVec = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            dvn = vn_vec(k)-prim(2)
            dvt = vt_vec(k)-prim(3)
            v2  = dvn*dvn+dvt*dvt

            GetHeatFlux_LocalVec(1) = GetHeatFlux_LocalVec(1)+0.5_KREAL*wk*dvn*(v2*h_vec(k)+b_vec(k))
            GetHeatFlux_LocalVec(2) = GetHeatFlux_LocalVec(2)+0.5_KREAL*wk*dvt*(v2*h_vec(k)+b_vec(k))
        enddo
    end function GetHeatFlux_LocalVec
    
    !--------------------------------------------------
    !>Calculate shear stress tensor from flattened velocity-space distributions
    !>@param[in] h_vec,b_vec   :1D distribution functions
    !>@param[in] prim          :primary variables
    !>@return    Get_Stress_Vec:shear stress of tau_xx, tau_xy, tau_yy
    !--------------------------------------------------
    function Get_Stress_Vec(h_vec,b_vec,prim)
        real(KREAL), intent(in), dimension(:)             :: h_vec,b_vec !Flattened input distributions
        real(KREAL), intent(in)                           :: prim(4)
        real(KREAL)                                       :: Get_Stress_Vec(3)
        real(KREAL)                                       :: Pressure !Scalar pressure from flattened moment
        integer(KINT)                                     :: k !Flattened velocity-space index
        real(KREAL)                                       :: wk,dvn,dvt,hK !Reusable Nvel-loop scalars

        !Accumulate stress directly from flattened velocity-space views.
        Pressure = GetPressure_Vec(h_vec,b_vec,prim)
        Get_Stress_Vec = 0.0_KREAL
        do k=1,Nvel
            wk  = wVel(k)
            dvn = uVel(k)-prim(2)
            dvt = vVel(k)-prim(3)
            hK  = h_vec(k)

            Get_Stress_Vec(1) = Get_Stress_Vec(1)+wk*dvn*dvn*hK
            Get_Stress_Vec(2) = Get_Stress_Vec(2)+wk*dvn*dvt*hK
            Get_Stress_Vec(3) = Get_Stress_Vec(3)+wk*dvt*dvt*hK
        enddo
        Get_Stress_Vec(1) = Get_Stress_Vec(1)-Pressure
        Get_Stress_Vec(3) = Get_Stress_Vec(3)-Pressure
    end function Get_Stress_Vec

    !--------------------------------------------------
    !>Use an absolute residual only when its reference is explicitly classified
    !>as near zero; otherwise preserve the relative residual definition.
    !--------------------------------------------------
    pure elemental real(KREAL) function RelativeOrAbsoluteResidual(value,reference,absolute_reference_threshold)
        real(KREAL), intent(in) :: value,reference,absolute_reference_threshold
        if (reference<=absolute_reference_threshold) then
            RelativeOrAbsoluteResidual=value
        else
            RelativeOrAbsoluteResidual=value/reference
        endif
    end function RelativeOrAbsoluteResidual

    pure elemental logical function IsNearZeroResidualScale(component_scale,ordinary_scale)
        real(KREAL), intent(in) :: component_scale,ordinary_scale
        IsNearZeroResidualScale=component_scale<=sqrt(epsilon(1.0_KREAL))* &
            max(ordinary_scale,tiny(1.0_KREAL))
    end function IsNearZeroResidualScale

end module Tools

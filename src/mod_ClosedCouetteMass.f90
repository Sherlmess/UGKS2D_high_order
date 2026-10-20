! Opt-in uniform rescaling; no mass collectives are called on the disabled path.
module mod_ClosedCouetteMass
    use flow_data
    use Flow_DIS
    use mod_ImplicitState
    use mod_ImplicitCommon, only: MomentConserved, IsPhysicalState
    use Tools, only: GetPrimary, GetTau, GetHeatFlux_Vec
    use mpi
    use Boundary, only: Set_Ghost_PhysBoundary
    use mpi_wrapper, only: MPIUpdate_CellConVars, MPIUpdate_CellVDF, MPIUpdate_Conner_conVars, MPIUpdate_Conner_VDF
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(KREAL), save :: closed_mass_target=0.0_KREAL
    real(KREAL), save :: mass_scale=1.0_KREAL, mass_scale_seconds=0.0_KREAL, mass_refresh_seconds=0.0_KREAL
    real(KREAL), save :: correction_ledger(4)=0.0_KREAL
    real(KREAL), save :: mass_before(5),mass_after(5),mass_before_metrics(5),mass_after_metrics(5),mass_shape_error(4)
    real(KREAL), save :: update_before(4),update_after(4),scale_increment(4)
    logical, save :: target_loaded=.false., snapshot_active=.false.
    integer, save :: restart_step=0
    real(KREAL), allocatable :: actual_tau(:,:)
    ! Scalar derived-type copies retain pointer associations. Backing storage is copied separately.
    type(CellCenter), allocatable :: saved_ctr(:,:)
    type(CellInterface), allocatable :: saved_vface(:,:),saved_hface(:,:)
    real(KREAL) :: saved_ledger(4),saved_diagnostics(23)
    real(KREAL),allocatable :: saved_actual_tau(:,:)
    real(KREAL), allocatable :: saved_cell_h_store(:,:,:)
    real(KREAL), allocatable :: saved_cell_b_store(:,:,:)
    real(KREAL), allocatable :: saved_cell_sh_store(:,:,:,:)
    real(KREAL), allocatable :: saved_cell_sb_store(:,:,:,:)
    real(KREAL), allocatable :: saved_v_l_flux_h_store(:,:,:)
    real(KREAL), allocatable :: saved_v_l_flux_b_store(:,:,:)
    real(KREAL), allocatable :: saved_v_r_flux_h_store(:,:,:)
    real(KREAL), allocatable :: saved_v_r_flux_b_store(:,:,:)
    real(KREAL), allocatable :: saved_h_l_flux_h_store(:,:,:)
    real(KREAL), allocatable :: saved_h_l_flux_b_store(:,:,:)
    real(KREAL), allocatable :: saved_h_r_flux_h_store(:,:,:)
    real(KREAL), allocatable :: saved_h_r_flux_b_store(:,:,:)
    real(KREAL), allocatable :: saved_res_macro_flux_store(:,:,:)
    real(KREAL), allocatable :: saved_res_macro_store(:,:,:)
    real(KREAL), allocatable :: saved_res_h_flux_store(:,:,:)
    real(KREAL), allocatable :: saved_res_b_flux_store(:,:,:)
    real(KREAL), allocatable :: saved_rate_w(:,:,:)
    real(KREAL), allocatable :: saved_stream_h(:,:,:)
    real(KREAL), allocatable :: saved_stream_b(:,:,:)
    real(KREAL), allocatable :: saved_delta_w(:,:,:)
    real(KREAL), allocatable :: saved_delta_h(:,:,:)
    real(KREAL), allocatable :: saved_delta_b(:,:,:)
    real(KREAL), allocatable :: saved_w_tilde(:,:,:)
    real(KREAL), allocatable :: saved_tau_tilde(:,:)
    real(KREAL), allocatable :: saved_q_frozen(:,:,:)
    real(KREAL), allocatable :: saved_equilibrium_h(:,:,:)
    real(KREAL), allocatable :: saved_equilibrium_b(:,:,:)
contains
    subroutine SaveMassTransaction()
        if (snapshot_active) call MassAbort("nested mass transaction")
        saved_ctr=ctr; saved_vface=vface; saved_hface=hface
        saved_diagnostics=[implicit_macro_residual,implicit_macro_relative,implicit_macro_absolute, &
            implicit_h_norm,implicit_b_norm,implicit_h_relative,implicit_b_relative,implicit_h_absolute,implicit_b_absolute, &
            implicit_moment_error,implicit_mass_defect,implicit_ncfl_active,implicit_alpha,mass_scale]
        saved_ledger=correction_ledger
        if (allocated(actual_tau)) saved_actual_tau=actual_tau
        saved_cell_h_store=cell_h_store
        saved_cell_b_store=cell_b_store
        saved_cell_sh_store=cell_sh_store
        saved_cell_sb_store=cell_sb_store
        saved_v_l_flux_h_store=v_l_flux_h_store
        saved_v_l_flux_b_store=v_l_flux_b_store
        saved_v_r_flux_h_store=v_r_flux_h_store
        saved_v_r_flux_b_store=v_r_flux_b_store
        saved_h_l_flux_h_store=h_l_flux_h_store
        saved_h_l_flux_b_store=h_l_flux_b_store
        saved_h_r_flux_h_store=h_r_flux_h_store
        saved_h_r_flux_b_store=h_r_flux_b_store
        saved_res_macro_flux_store=res_macro_flux_store
        saved_res_macro_store=res_macro_store
        saved_res_h_flux_store=res_h_flux_store
        saved_res_b_flux_store=res_b_flux_store
        saved_rate_w=rate_w
        saved_stream_h=stream_h
        saved_stream_b=stream_b
        saved_delta_w=delta_w
        saved_delta_h=delta_h
        saved_delta_b=delta_b
        saved_w_tilde=w_tilde
        saved_tau_tilde=tau_tilde
        saved_q_frozen=q_frozen
        saved_equilibrium_h=equilibrium_h
        saved_equilibrium_b=equilibrium_b
        snapshot_active=.true.
    end subroutine
    subroutine RestoreMassTransaction()
        ctr=saved_ctr; vface=saved_vface; hface=saved_hface
        implicit_macro_residual=saved_diagnostics(1:4);implicit_macro_relative=saved_diagnostics(5:8)
        implicit_macro_absolute=saved_diagnostics(9:12)
        implicit_h_norm=saved_diagnostics(13);implicit_b_norm=saved_diagnostics(14)
        implicit_h_relative=saved_diagnostics(15);implicit_b_relative=saved_diagnostics(16)
        implicit_h_absolute=saved_diagnostics(17);implicit_b_absolute=saved_diagnostics(18)
        implicit_moment_error=saved_diagnostics(19);implicit_mass_defect=saved_diagnostics(20)
        implicit_ncfl_active=saved_diagnostics(21);implicit_alpha=saved_diagnostics(22);mass_scale=saved_diagnostics(23)
        correction_ledger=saved_ledger
        if (allocated(saved_actual_tau)) actual_tau=saved_actual_tau
        cell_h_store=saved_cell_h_store
        cell_b_store=saved_cell_b_store
        cell_sh_store=saved_cell_sh_store
        cell_sb_store=saved_cell_sb_store
        v_l_flux_h_store=saved_v_l_flux_h_store
        v_l_flux_b_store=saved_v_l_flux_b_store
        v_r_flux_h_store=saved_v_r_flux_h_store
        v_r_flux_b_store=saved_v_r_flux_b_store
        h_l_flux_h_store=saved_h_l_flux_h_store
        h_l_flux_b_store=saved_h_l_flux_b_store
        h_r_flux_h_store=saved_h_r_flux_h_store
        h_r_flux_b_store=saved_h_r_flux_b_store
        res_macro_flux_store=saved_res_macro_flux_store
        res_macro_store=saved_res_macro_store
        res_h_flux_store=saved_res_h_flux_store
        res_b_flux_store=saved_res_b_flux_store
        rate_w=saved_rate_w
        stream_h=saved_stream_h
        stream_b=saved_stream_b
        delta_w=saved_delta_w
        delta_h=saved_delta_h
        delta_b=saved_delta_b
        w_tilde=saved_w_tilde
        tau_tilde=saved_tau_tilde
        q_frozen=saved_q_frozen
        equilibrium_h=saved_equilibrium_h
        equilibrium_b=saved_equilibrium_b
    end subroutine
    subroutine ReleaseMassTransaction()
        deallocate(saved_ctr,saved_vface,saved_hface)
        if (allocated(saved_actual_tau)) deallocate(saved_actual_tau)
        deallocate(saved_cell_h_store)
        deallocate(saved_cell_b_store)
        deallocate(saved_cell_sh_store)
        deallocate(saved_cell_sb_store)
        deallocate(saved_v_l_flux_h_store)
        deallocate(saved_v_l_flux_b_store)
        deallocate(saved_v_r_flux_h_store)
        deallocate(saved_v_r_flux_b_store)
        deallocate(saved_h_l_flux_h_store)
        deallocate(saved_h_l_flux_b_store)
        deallocate(saved_h_r_flux_h_store)
        deallocate(saved_h_r_flux_b_store)
        deallocate(saved_res_macro_flux_store)
        deallocate(saved_res_macro_store)
        deallocate(saved_res_h_flux_store)
        deallocate(saved_res_b_flux_store)
        deallocate(saved_rate_w)
        deallocate(saved_stream_h)
        deallocate(saved_stream_b)
        deallocate(saved_delta_w)
        deallocate(saved_delta_h)
        deallocate(saved_delta_b)
        deallocate(saved_w_tilde)
        deallocate(saved_tau_tilde)
        deallocate(saved_q_frozen)
        deallocate(saved_equilibrium_h)
        deallocate(saved_equilibrium_b)
        snapshot_active=.false.
    end subroutine
    subroutine MassAbort(message)
        character(len=*),intent(in) :: message
        if (myID==0) write(*,'(a)') 'ClosedCouetteMassRescale: '//message
        call MPI_Abort(MPI_COMM_WORLD,941,ierr)
    end subroutine

    pure logical function ValidMass(m)
        real(KREAL),intent(in) :: m
        ValidMass=ieee_is_finite(m) .and. m>0.0_KREAL
    end function

    subroutine MassStatistics(totals,metrics,ok)
        ! totals = integral W(1:4), integral discrete h density. Physical cells only.
        ! metrics = max moment error, max negative h/b fraction, minimum rho/internal energy.
        real(KREAL),intent(out) :: totals(5),metrics(5)
        logical,intent(out) :: ok
        real(KREAL) :: local(5),mx(3),mn(2),mm(4),e,nh,nb
        logical :: local_ok
        integer :: i,j
        local=0;mx=0;mn=huge(1.0_KREAL);local_ok=.true.
        do j=1,ny; do i=1,nx
            call MomentConserved(ctr(i,j)%h_vec,ctr(i,j)%b_vec,uVel,vVel,wVel,mm)
            local_ok=local_ok .and. all(ieee_is_finite(ctr(i,j)%h_vec)) .and. &
                all(ieee_is_finite(ctr(i,j)%b_vec)) .and. IsPhysicalState(ctr(i,j)%conVars) .and. &
                IsPhysicalState(mm) .and. ValidMass(ctr(i,j)%area)
            local(1:4)=local(1:4)+ctr(i,j)%area*ctr(i,j)%conVars
            local(5)=local(5)+ctr(i,j)%area*mm(1)
            mx(1)=max(mx(1),maxval(abs(ctr(i,j)%conVars-mm)/max(abs(ctr(i,j)%conVars),1.0_KREAL)))
            nh=sum(wVel*max(-ctr(i,j)%h_vec,0.0_KREAL))/max(sum(wVel*abs(ctr(i,j)%h_vec)),tiny(1.0_KREAL))
            nb=sum(wVel*max(-ctr(i,j)%b_vec,0.0_KREAL))/max(sum(wVel*abs(ctr(i,j)%b_vec)),tiny(1.0_KREAL))
            mx(2)=max(mx(2),nh);mx(3)=max(mx(3),nb)
            if (ValidMass(ctr(i,j)%conVars(1))) then
                e=ctr(i,j)%conVars(4)-0.5_KREAL*sum(ctr(i,j)%conVars(2:3)**2)/ctr(i,j)%conVars(1)
                mn=min(mn,[ctr(i,j)%conVars(1),e])
            endif
        enddo; enddo
        call MPI_Allreduce(local,totals,5,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(mx,metrics(1:3),3,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(mn,metrics(4:5),2,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(local_ok,ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=ok .and. all(ieee_is_finite(totals)) .and. all(ieee_is_finite(metrics)) .and. &
            metrics(1)<=Implicit_MomentTol
    end subroutine

    subroutine CheckMassState(ok)
        logical,intent(out) :: ok
        real(KREAL) :: t(5),m(5)
        call MassStatistics(t,m,ok)
    end subroutine

    subroutine InitializeClosedMass()
        logical :: ok
        real(KREAL) :: t(5),m(5),low,high
        integer :: i,j
        logical :: rectangular,global_rectangular
        rectangular=.true.
        do j=1,ny+1;do i=1,nx+1
            rectangular=rectangular .and. abs(gridx(i,j)-gridx(i,1))<=1.0e-14_KREAL .and. &
                abs(gridy(i,j)-gridy(1,j))<=1.0e-14_KREAL
        enddo;enddo
        call MPI_Allreduce(rectangular,global_rectangular,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global_rectangular) call MassAbort('requires verified Cartesian Couette geometry')
        call MassStatistics(t,m,ok)
        if (.not.ok) call MassAbort('initial macro/micro or physical-state gate failed')
        if (Iflag_InitFlow==Init_Freestream) then
            closed_mass_target=t(1)
            correction_ledger=0.0_KREAL
            target_loaded=.true.
        else if (.not.target_loaded) then
            call MassAbort('restart lacks persistent M0; current mass is not a target')
        endif
        call MPI_Allreduce(closed_mass_target,low,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(closed_mass_target,high,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        if (.not.ValidMass(low) .or. low/=high) call MassAbort('invalid or inconsistent M0 across ranks')
        allocate(actual_tau(nx,ny))
        do j=1,ny;do i=1,nx
            actual_tau(i,j)=GetTau(GetPrimary(ctr(i,j)%conVars))
        enddo;enddo
        if (myID==0) write(*,'(a,es25.16e3)') 'persistent_M0=',closed_mass_target
    end subroutine

    subroutine RescaleClosedMass(ok)
        logical,intent(out) :: ok
        logical :: local_ok,global_ok
        integer :: i,j
        real(KREAL) :: start,prim0(4),prim1(4),err(4),globerr(4),oldmass(4),change(4),global_change(4)
        real(KREAL) :: h0(Nvel),b0(Nvel),w0(4),tau
        start=MPI_Wtime()
        call MassStatistics(mass_before,mass_before_metrics,ok)
        if (.not.ok .or. .not.ValidMass(mass_before(1)) .or. .not.ValidMass(closed_mass_target)) then
            ok=.false.;return
        endif
        mass_scale=closed_mass_target/mass_before(1)
        if (.not.ValidMass(mass_scale)) then
            ok=.false.;return
        endif
        ! Preflight every scaled candidate before any GetPrimary call can abort on overflow.
        local_ok=.true.
        do j=1,ny;do i=1,nx
            local_ok=local_ok .and. IsPhysicalState(mass_scale*ctr(i,j)%conVars) .and. &
                all(ieee_is_finite(mass_scale*ctr(i,j)%h_vec)) .and. all(ieee_is_finite(mass_scale*ctr(i,j)%b_vec))
        enddo;enddo
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
        if (.not.ok) return
        err=0.0_KREAL;change=0.0_KREAL
        do j=1,ny; do i=1,nx
            w0=ctr(i,j)%conVars;h0=ctr(i,j)%h_vec;b0=ctr(i,j)%b_vec
            prim0=GetPrimary(w0)
            ctr(i,j)%conVars=mass_scale*w0
            ctr(i,j)%h_vec=mass_scale*h0
            ctr(i,j)%b_vec=mass_scale*b0
            prim1=GetPrimary(ctr(i,j)%conVars)
            err(1)=max(err(1),maxval(abs(prim1(2:3)-prim0(2:3))))
            err(2)=max(err(2),abs(1.0_KREAL/prim1(4)-1.0_KREAL/prim0(4)))
            err(3)=max(err(3),maxval(abs(ctr(i,j)%h_vec/ctr(i,j)%conVars(1)-h0/w0(1))))
            err(4)=max(err(4),maxval(abs(ctr(i,j)%b_vec/ctr(i,j)%conVars(1)-b0/w0(1))))
            change=change+ctr(i,j)%area*(w0-saved_ctr(i,j)%conVars)
        enddo; enddo
        call MPI_Allreduce(err,globerr,4,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        mass_shape_error=globerr
        call MPI_Allreduce(change,global_change,4,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        update_before=global_change
        call MassStatistics(mass_after,mass_after_metrics,ok)
        scale_increment=mass_after(1:4)-mass_before(1:4)
        update_after=update_before+scale_increment
        mass_scale_seconds=mass_scale_seconds+MPI_Wtime()-start
        if (.not.ok) return
        start=MPI_Wtime();local_ok=.true.
        do j=1,ny; do i=1,nx
            ! Recompute the actual pressure-dependent tau; do not retain predictor tau.
            prim1=GetPrimary(ctr(i,j)%conVars);tau=GetTau(prim1)
            local_ok=local_ok .and. ValidMass(tau)
            actual_tau(i,j)=tau;tau_tilde(i,j)=tau
            w_tilde(:,i,j)=ctr(i,j)%conVars
            q_frozen(:,i,j)=GetHeatFlux_Vec(ctr(i,j)%h_vec,ctr(i,j)%b_vec,prim1)
        enddo; enddo
        ! Predictor equilibrium belongs to the old direction and must be rebuilt before reuse.
        equilibrium_h=0.0_KREAL;equilibrium_b=0.0_KREAL
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        ok=global_ok
        if (ok) then
            call Set_Ghost_PhysBoundary()
            call MPIUpdate_CellConVars()
            call MPIUpdate_CellVDF()
            call MPIUpdate_Conner_conVars()
            call MPIUpdate_Conner_VDF()
            implicit_moment_error=mass_after_metrics(1)
            correction_ledger=correction_ledger+scale_increment
        endif
        mass_refresh_seconds=mass_refresh_seconds+MPI_Wtime()-start
    end subroutine

    function RestartIdentity() result(identity)
        real(KREAL) :: identity(31)
        identity=[real(Nx_all,KREAL),real(Ny_all,KREAL),real(NP_X,KREAL),real(NP_Y,KREAL), &
            real(Nu,KREAL),real(Nv,KREAL),CFL,Implicit_nCFL,real(ck,KREAL),Pr,omega,Tinf,Tw,CouetteMuRef, &
            WallUx,real(BC_Type,KREAL),real(flag_vislaw,KREAL),real(FluxWall_type,KREAL),ACC, &
            real(BG_geometry,KREAL),real(Quadrature_method,KREAL), &
            merge(1.0_KREAL,0.0_KREAL,Implicit_ZhuMovingWall),merge(1.0_KREAL,0.0_KREAL,Implicit_InitHaloResync), &
            sum(wVel),sum(wVel*uVel*uVel)]
    end function

    subroutine WriteMassRestart()
        integer :: unit,i,j,ios,version
        logical :: ok,global_ok
        character(len=64) :: name
        version=4
        write(name,'(a,i0,a)') 'closed_mass_rank',myID,'.bin'
        open(newunit=unit,file=trim(name),access='stream',form='unformatted',status='replace',iostat=ios)
        ok=ios==0
        if (ok) then
            write(unit,iostat=ios) version,nx,ny,Nvel,LAP,restart_step+iter,closed_mass_target,correction_ledger,RestartIdentity()
            ok=ios==0
            do j=1-LAP,ny+LAP; do i=1-LAP,nx+LAP
                if (ok) then
                    write(unit,iostat=ios) ctr(i,j)%conVars,ctr(i,j)%h_vec,ctr(i,j)%b_vec
                    ok=ios==0
                endif
            enddo; enddo
            close(unit)
        endif
        call MPI_Allreduce(ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global_ok) call MassAbort('cannot write complete M0-bearing restart')
    end subroutine

    subroutine ReadMassRestart()
        integer :: unit,i,j,ios,header(6)
        logical :: ok,global_ok
        real(KREAL) :: identity(31),w(4),h(Nvel),b(Nvel),lo,hi,lmin(4),lmax(4)
        integer :: smin,smax
        character(len=64) :: name
        write(name,'(a,i0,a)') 'closed_mass_rank',myID,'.bin'
        open(newunit=unit,file=trim(name),access='stream',form='unformatted',status='old',iostat=ios)
        ok=ios==0
        if (ok) then
            read(unit,iostat=ios) header,closed_mass_target,correction_ledger,identity
            ok=ios==0
            if (ok) ok=all(header(1:5)==[4,nx,ny,Nvel,LAP]) .and. &
                all(identity==RestartIdentity()) .and. ValidMass(closed_mass_target) .and. &
                all(ieee_is_finite(correction_ledger))
            if (ok) then
                restart_step=header(6)
                do j=1-LAP,ny+LAP; do i=1-LAP,nx+LAP
                    if (ok) then
                        read(unit,iostat=ios) w,h,b
                        ok=ios==0
                        if (ok .and. i>=1 .and. i<=nx .and. j>=1 .and. j<=ny) then
                            ! Bind the M0 metadata to the matching native physical restart, not a different flow pair.
                            ok=all(w==ctr(i,j)%conVars) .and. all(h==ctr(i,j)%h_vec) .and. all(b==ctr(i,j)%b_vec)
                        endif
                        if (ok) then
                            ctr(i,j)%conVars=w;ctr(i,j)%h_vec=h;ctr(i,j)%b_vec=b
                        endif
                    endif
                enddo; enddo
            endif
            close(unit)
        endif
        call MPI_Allreduce(ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global_ok) call MassAbort('missing, mismatched, or invalid M0-bearing restart; refusing to redefine M0')
        call MPI_Allreduce(closed_mass_target,lo,1,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(closed_mass_target,hi,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        if (lo/=hi) call MassAbort('restart M0 differs between ranks')
        call MPI_Allreduce(correction_ledger,lmin,4,MPI_KREAL,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(correction_ledger,lmax,4,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(restart_step,smin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(restart_step,smax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        if (any(lmin/=lmax) .or. smin/=smax) call MassAbort('restart ledger or origin step differs across ranks')
        target_loaded=.true.
        ! Original native restart iteration semantics are additional steps; preserve origin separately.
    end subroutine
end module mod_ClosedCouetteMass

! Optional bounded wiring observations. No numerical changes or extra sweeps.
module mod_NativeWiringChecks
    use flow_data
    use mod_ImplicitState
    use mod_ImplicitCommon,only:MomentConserved
    use mpi
    use,intrinsic::ieee_arithmetic
    implicit none
    logical,save::enabled=.false.,initialized=.false.,dumped=.false.
contains
    logical function WiringChecksEnabled()
        character(len=8)::value
        integer::status
        if(.not.initialized)then
            call get_environment_variable('NATIVE_WIRING_CHECKS',value,status=status)
            enabled=status==0.and.trim(value)=='1';initialized=.true.
        endif
        WiringChecksEnabled=enabled
    end function
    subroutine DumpPhysicalOperators()
        integer::fu,i,j
        character(len=64)::name
        if(.not.WiringChecksEnabled().or.dumped.or..not.UnifiedPhysicalOperators)return
        write(name,'(A,I0,A)')'physical_operator_rank',myID,'.bin'
        open(newunit=fu,file=trim(name),status='replace',access='stream',form='unformatted')
        write(fu)nx,ny,Nvel,Block_XY(myID)%Origin_I,Block_XY(myID)%Origin_J
        write(fu)ctr(1,1)%dt_local,uVel,vVel,wVel
        do j=1,ny;do i=1,nx
            write(fu)ctr(i,j)%area,ctr(i,j)%conVars,ctr(i,j)%h_vec,ctr(i,j)%b_vec
        enddo;enddo
        write(fu)rate_w,stream_h,stream_b,physical_equilibrium_h,physical_equilibrium_b, &
            physical_collision_h,physical_collision_b
        do j=1,ny;do i=1,nx+1
            write(fu)vface(i,j)%L_flux,vface(i,j)%R_flux,vface(i,j)%L_flux_h_vec, &
                vface(i,j)%R_flux_h_vec,vface(i,j)%L_flux_b_vec,vface(i,j)%R_flux_b_vec
        enddo;enddo
        do j=1,ny+1;do i=1,nx
            write(fu)hface(i,j)%L_flux,hface(i,j)%R_flux,hface(i,j)%L_flux_h_vec, &
                hface(i,j)%R_flux_h_vec,hface(i,j)%L_flux_b_vec,hface(i,j)%R_flux_b_vec
        enddo;enddo
        close(fu);dumped=.true.
    end subroutine
    subroutine RecordWiringState()
        real(KREAL)::local(4),global(4),moment(4),gap,global_gap,err,errg
        integer::i,j,fu,counts(3),lo(3),hi(3)
        logical::ok,allok
        if(.not.WiringChecksEnabled())return
        local=0;gap=0;ok=.true.;err=0
        do j=1,ny;do i=1,nx
            call MomentConserved(ctr(i,j)%h_vec,ctr(i,j)%b_vec,uVel,vVel,wVel,moment)
            local=local+ctr(i,j)%area*ctr(i,j)%conVars
            gap=max(gap,maxval(abs(moment-ctr(i,j)%conVars)/max(1.0_KREAL,abs(ctr(i,j)%conVars))))
            ok=ok.and.all(ieee_is_finite(ctr(i,j)%conVars)).and.all(ieee_is_finite(ctr(i,j)%h_vec)).and.all(ieee_is_finite(ctr(i,j)%b_vec))
        enddo;enddo
        if(allocated(delta_h))err=max(maxval(abs(delta_w)),maxval(abs(delta_h)),maxval(abs(delta_b)))
        call MPI_Allreduce(local,global,4,MPI_KREAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(gap,global_gap,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(err,errg,1,MPI_KREAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(ok,allok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        counts=[macro_pr_passes,r4_forward_scans,r4_backward_scans]
        call MPI_Allreduce(counts,lo,3,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(counts,hi,3,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        if(myID/=0)return
        if(iter==0)then
            open(newunit=fu,file='wiring_states.csv',status='replace')
            write(fu,'(A)')'iter,mass,Px,Py,E,moment_gap,max_correction,finite,PR_passes_min,forward_min,backward_min,PR_passes_max,forward_max,backward_max'
        else
            open(newunit=fu,file='wiring_states.csv',status='old',position='append')
        endif
        write(fu,'(I0,6(",",ES24.16),",",L1,6(",",I0))')iter,global,global_gap,errg,allok,lo,hi
        close(fu)
    end subroutine
end module

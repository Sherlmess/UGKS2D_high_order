module mod_S2SteadyIteration
    use flow_data
    use mpi
    use mod_ImplicitState
    use mod_ImplicitBoundary, only: BuildStructuredMirrorMaps
    use mod_ImplicitDiagnostics, only: InitializeImplicitDiagnostics
    use mod_ImplicitCommon, only: MomentConserved,ProjectToConservedMoments
    use mpi_wrapper, only: MPIExchangeDeltaMacro,MPIExchangeDeltaMicro
    use mod_ImplicitMacro, only: SolveMacroPR2,BuildPredictedMacroState,ComputeMacroDefectNorm
    use mod_ImplicitMicro, only: BuildPredictedEquilibrium,SolveMicroLUSGS,RecoverMacroFromKineticState, &
        ComputeMicroDefectNorm
    use mod_SteadyImplicit, only: FindAdmissibleMicroAlpha
    use mod_S2NativeAdapter, only: S2NativeInvalidateTargets
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: InitializeS2Steady,S2NativeApplySteadyUpdate
    real(KREAL),parameter::gate=2.0e-12_KREAL
    real(KREAL),save::total0(4)=0.0_KREAL
    real(KREAL),save::prediction_gap=0.0_KREAL,global_drift=0.0_KREAL
    integer,save::update_count=0
    logical,save::r1_diagnostic=.false.
contains
    subroutine InitializeS2Steady()
        real(KREAL)::moment(4),err
        integer::i,j,env_status
        character(len=8)::env_value
        ! The native low-order Jacobian uses GetTau(mu_ref).  The S2 target
        ! itself uses the frozen constant tau=0.1; Ma=0 in this periodic case
        ! would otherwise leave only the preconditioner's mu_ref at zero.
        mu_ref=0.1_KREAL
        call AllocateImplicitState()
        call BuildStructuredMirrorMaps()
        call InitializeImplicitDiagnostics()
        err=0.0_KREAL;total0=0.0_KREAL
        do j=1,ny;do i=1,nx
            call MomentConserved(ctr(i,j)%h_vec,ctr(i,j)%b_vec,uVel,vVel,wVel,moment)
            err=max(err,maxval(abs(moment-ctr(i,j)%conVars))/max(1.0_KREAL,maxval(abs(ctr(i,j)%conVars))))
            total0=total0+ctr(i,j)%conVars
        enddo;enddo
        if(err>gate)call AbortI1('initial macro/micro moment mismatch')
        env_value=''
        call get_environment_variable('S2I1_R1_DIAG',env_value,status=env_status)
        r1_diagnostic=env_status==0.and.trim(env_value)=='1'
        s2_steady_coupling_on=.true.
        if(myID==0)then
            open(unit=166,file='s2i1_updates.csv',status='replace',action='write')
            write(166,'(A)')'iter,ncfl,alpha,retries,prediction_gap,global_drift'
            if(r1_diagnostic)then
                open(unit=167,file='s2i1_stages.csv',status='replace',action='write')
                write(167,'(A)')'iter,retry,stage,dW_sum1,dW_sum2,dW_sum3,dW_sum4,'// &
                    'dMf_sum1,dMf_sum2,dMf_sum3,dMf_sum4,dG_sum1,dG_sum2,dG_sum3,dG_sum4,'// &
                    'local_moment_gap,macro_linear_defect,micro_h_linear_defect,'// &
                    'micro_b_linear_defect,equilibrium_increment_moment_gap'
            endif
            open(unit=168,file='s2i1_constraint.csv',status='replace',action='write')
            write(168,'(A)')'iter,retry,phase,norm_before,norm_after,correction_norm,'// &
                'global_before1,global_before2,global_before3,global_before4,'// &
                'global_after1,global_after2,global_after3,global_after4,'// &
                'local_moment_gap_before,local_moment_gap_after'
        endif
    end subroutine InitializeS2Steady

    subroutine S2NativeApplySteadyUpdate(accepted)
        logical,intent(out)::accepted
        real(KREAL),allocatable::oldw(:,:,:),oldh(:,:,:),oldb(:,:,:)
        real(KREAL)::macro_alpha,moment(4),total(4),gap,drift
        integer::i,j,retry,fu
        logical::ok
        accepted=.false.
        if(.not.s2_steady_rhs)call AbortI1('low-order solver did not receive S2 RHS')
        allocate(oldw(4,nx,ny),oldh(Nvel,nx,ny),oldb(Nvel,nx,ny))
        do j=1,ny;do i=1,nx
            oldw(:,i,j)=ctr(i,j)%conVars
            oldh(:,i,j)=ctr(i,j)%h_vec
            oldb(:,i,j)=ctr(i,j)%b_vec
        enddo;enddo
        implicit_ncfl_active=Implicit_nCFL
        do retry=0,Implicit_MaxRetries
            call ResetImplicitCorrections()
            equilibrium_increment_h=0.0_KREAL;equilibrium_increment_b=0.0_KREAL
            implicit_retry_count=retry
            call SolveMacroPR2(ok)
            if(ok)call RecordR1Stage('pr2',retry)
            if(ok)then
                call ConstrainS2MacroDirection(ok)
                if(ok)call RecordR1Stage('pr2_constrained',retry)
            endif
            if(.not.ok.and.myID==0)write(*,'(A,I0,A,ES16.8)')'S2I1 reject macro retry=',retry,' nCFL=',implicit_ncfl_active
            if(ok)then
                call BuildPredictedMacroState(macro_alpha,ok)
                if(ok)call RecordR1Stage('predicted_w',retry)
                if(.not.ok.and.myID==0)write(*,'(A,I0)')'S2I1 reject predicted macro retry=',retry
            endif
            if(ok)then
                call BuildPredictedEquilibrium(ok)
                if(ok)call RecordR1Stage('predicted_g',retry)
                if(.not.ok.and.myID==0)write(*,'(A,I0)')'S2I1 reject predicted equilibrium retry=',retry
            endif
            if(ok)then
                call SolveMicroLUSGS(ok,freeze_wall_closure=.true.)
                if(ok)call RecordR1Stage('micro_post',retry)
                if(.not.ok.and.myID==0)write(*,'(A,I0)')'S2I1 reject micro retry=',retry
            endif
            if(ok)then
                call ConstrainS2SteadyDirection(ok)
                if(ok)call RecordR1Stage('micro_constrained',retry)
            endif
            if(ok)then
                call FindAdmissibleMicroAlpha(implicit_alpha,ok,.false.)
                if(.not.ok.and.myID==0)write(*,'(A,I0)')'S2I1 reject line search retry=',retry
            endif
            if(ok)exit
            implicit_ncfl_active=0.5_KREAL*implicit_ncfl_active
        enddo
        if(.not.ok)return
        gap=0.0_KREAL
        do j=1,ny;do i=1,nx
            ctr(i,j)%h_vec=ctr(i,j)%h_vec+implicit_alpha*delta_h(:,i,j)
            ctr(i,j)%b_vec=ctr(i,j)%b_vec+implicit_alpha*delta_b(:,i,j)
            call MomentConserved(ctr(i,j)%h_vec,ctr(i,j)%b_vec,uVel,vVel,wVel,moment)
            gap=max(gap,maxval(abs(moment-(oldw(:,i,j)+implicit_alpha*delta_w(:,i,j)))/ &
                max(1.0_KREAL,abs(oldw(:,i,j)+implicit_alpha*delta_w(:,i,j)))))
        enddo;enddo
        call RecoverMacroFromKineticState(ok)
        if(ok)then
            total=0.0_KREAL
            do j=1,ny;do i=1,nx
                total=total+ctr(i,j)%conVars
            enddo;enddo
            drift=maxval(abs(total-total0)/max(1.0_KREAL,abs(total0)))
            ok=ieee_is_finite(drift).and.drift<=gate.and.gap<=gate
        else
            drift=huge(1.0_KREAL)
        endif
        if(.not.ok)then
            if(myID==0)then
                write(*,'(A,2ES24.16)')'S2I1 reject recovered/global drift: ',gap,drift
                open(newunit=fu,file='s2i1_failed_proposal.bin',status='replace',access='stream',form='unformatted')
                write(fu)nx,ny,Nvel,oldw,oldh,oldb
                do j=1,ny;do i=1,nx
                    write(fu)ctr(i,j)%conVars,ctr(i,j)%h_vec,ctr(i,j)%b_vec
                enddo;enddo
                close(fu)
            endif
            do j=1,ny;do i=1,nx
                ctr(i,j)%conVars=oldw(:,i,j)
                ctr(i,j)%h_vec=oldh(:,i,j)
                ctr(i,j)%b_vec=oldb(:,i,j)
            enddo;enddo
            call S2NativeInvalidateTargets()
            return
        endif
        prediction_gap=gap;global_drift=drift;update_count=update_count+1
        if(myID==0)then
            write(166,'(I0,",",ES24.16,",",ES24.16,",",I0,",",ES24.16,",",ES24.16)') &
                update_count,implicit_ncfl_active,implicit_alpha,retry,prediction_gap,global_drift
            flush(166)
        endif
        call S2NativeInvalidateTargets()
        accepted=.true.
    end subroutine S2NativeApplySteadyUpdate

    subroutine ConstrainS2MacroDirection(ok)
        logical,intent(out)::ok
        real(KREAL)::before(4),after(4),volume,mean(4),norm0,norm1,corr
        integer::i,j
        before=0.0_KREAL;after=0.0_KREAL;volume=0.0_KREAL
        norm0=0.0_KREAL;norm1=0.0_KREAL;corr=0.0_KREAL
        do j=1,ny;do i=1,nx
            before=before+ctr(i,j)%area*delta_w(:,i,j)
            volume=volume+ctr(i,j)%area
            norm0=norm0+ctr(i,j)%area*sum(delta_w(:,i,j)**2)
        enddo;enddo
        ok=ieee_is_finite(volume).and.volume>0.0_KREAL
        if(.not.ok)return
        mean=before/volume
        do j=1,ny;do i=1,nx
            delta_w(:,i,j)=delta_w(:,i,j)-mean
            after=after+ctr(i,j)%area*delta_w(:,i,j)
            norm1=norm1+ctr(i,j)%area*sum(delta_w(:,i,j)**2)
            corr=corr+ctr(i,j)%area*sum(mean**2)
        enddo;enddo
        call MPIExchangeDeltaMacro(delta_w)
        ok=all(ieee_is_finite(delta_w(:,1:nx,1:ny))).and. &
            maxval(abs(after))<=gate*max(1.0_KREAL,maxval(abs(total0*ctr(1,1)%area)))
        if(myID==0)then
            write(168,'(I0,",",I0,",",A,13(",",ES24.16))') &
                iter,implicit_retry_count,'macro',sqrt(norm0),sqrt(norm1),sqrt(corr), &
                before,after,0.0_KREAL,0.0_KREAL
            flush(168)
        endif
    end subroutine ConstrainS2MacroDirection

    subroutine ConstrainS2SteadyDirection(ok)
        logical,intent(out)::ok
        real(KREAL)::before(4),after(4),mb(4),ma(4),h0(Nvel),b0(Nvel)
        real(KREAL)::norm0,norm1,corr,gap0,gap1,vol
        integer::i,j
        logical::projection_ok
        before=0.0_KREAL;after=0.0_KREAL
        norm0=0.0_KREAL;norm1=0.0_KREAL;corr=0.0_KREAL;gap0=0.0_KREAL;gap1=0.0_KREAL
        ok=.true.
        do j=1,ny;do i=1,nx
            vol=ctr(i,j)%area
            h0=delta_h(:,i,j);b0=delta_b(:,i,j)
            call MomentConserved(h0,b0,uVel,vVel,wVel,mb)
            before=before+vol*mb
            gap0=max(gap0,maxval(abs(mb-delta_w(:,i,j))))
            norm0=norm0+vol*sum(wVel*(h0*h0+b0*b0))
            call ProjectToConservedMoments(delta_h(:,i,j),delta_b(:,i,j),delta_w(:,i,j), &
                uVel,vVel,wVel,ctr(i,j)%h_vec,ctr(i,j)%b_vec,projection_ok)
            if(.not.projection_ok)ok=.false.
            call MomentConserved(delta_h(:,i,j),delta_b(:,i,j),uVel,vVel,wVel,ma)
            after=after+vol*ma
            gap1=max(gap1,maxval(abs(ma-delta_w(:,i,j))))
            norm1=norm1+vol*sum(wVel*(delta_h(:,i,j)**2+delta_b(:,i,j)**2))
            corr=corr+vol*sum(wVel*((delta_h(:,i,j)-h0)**2+(delta_b(:,i,j)-b0)**2))
        enddo;enddo
        call MPIExchangeDeltaMicro(delta_h,delta_b)
        ok=ok.and.all(ieee_is_finite(delta_h(:,1:nx,1:ny))).and. &
            all(ieee_is_finite(delta_b(:,1:nx,1:ny))).and. &
            gap1<=gate.and.maxval(abs(after))<= &
            gate*max(1.0_KREAL,maxval(abs(total0*ctr(1,1)%area)))
        if(myID==0)then
            write(168,'(I0,",",I0,",",A,13(",",ES24.16))') &
                iter,implicit_retry_count,'micro',sqrt(norm0),sqrt(norm1),sqrt(corr), &
                before,after,gap0,gap1
            flush(168)
        endif
    end subroutine ConstrainS2SteadyDirection

    subroutine RecordR1Stage(stage,retry)
        character(len=*),intent(in)::stage
        integer,intent(in)::retry
        real(KREAL)::wd(4),md(4),gd(4),gm(4),em(4),gap,eqgap,mdef,hdef,bdef,vol
        integer::i,j
        if(.not.r1_diagnostic)return
        wd=0.0_KREAL;md=0.0_KREAL;gd=0.0_KREAL;gap=0.0_KREAL;eqgap=0.0_KREAL
        do j=1,ny;do i=1,nx
            vol=ctr(i,j)%area
            call MomentConserved(delta_h(:,i,j),delta_b(:,i,j),uVel,vVel,wVel,gm)
            call MomentConserved(equilibrium_increment_h(:,i,j),equilibrium_increment_b(:,i,j), &
                uVel,vVel,wVel,em)
            wd=wd+vol*delta_w(:,i,j);md=md+vol*gm;gd=gd+vol*em
            gap=max(gap,maxval(abs(gm-delta_w(:,i,j))))
            eqgap=max(eqgap,maxval(abs(em-delta_w(:,i,j))))
        enddo;enddo
        mdef=ComputeMacroDefectNorm();hdef=0.0_KREAL;bdef=0.0_KREAL
        if(stage=='micro_post'.or.stage=='micro_constrained')call ComputeMicroDefectNorm(hdef,bdef)
        if(myID==0)then
            write(167,'(I0,",",I0,",",A,17(",",ES24.16))') &
                iter,retry,trim(stage),wd,md,gd,gap,mdef,hdef,bdef,eqgap
            flush(167)
        endif
    end subroutine RecordR1Stage

    subroutine AbortI1(message)
        character(len=*),intent(in)::message
        if(myID==0)write(*,'(A)')'S2I1 FAILURE: '//trim(message)
        call MPI_Abort(MPI_COMM_WORLD,451,ierr)
    end subroutine AbortI1
end module mod_S2SteadyIteration

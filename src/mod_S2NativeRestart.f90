! Native W/h/b restart sidecar. Files use original writer/readers unchanged.
module mod_S2NativeRestart
    use iso_fortran_env, only:int64,int32
    use iso_c_binding, only:c_char,c_null_char,c_int64_t
    use flow_data
    use Flow_DIS, only:time_dis
    use mod_S2MPI2IO, only:physical_parameters,option_identity,c_rename
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    use mpi
    implicit none
    private
    public::ReadNativeMetadata,CommitNativeRestart
    character(len=16),parameter::schema='S2NATIVE-E2-V1'
    interface
        integer(c_int64_t) function file_hash(path) bind(C,name='s2_file_hash')
            import c_char,c_int64_t
            character(kind=c_char),intent(in)::path(*)
        end function
    end interface
contains
    subroutine require(ok,message)
        logical,intent(in)::ok
        character(len=*),intent(in)::message
        integer::e
        if(ok)return
        write(*,'(A,I0,A)')'NATIVE RESTART REJECTED BEFORE PUBLISH rank=',myID,': '//message
        call MPI_Abort(MPI_COMM_WORLD,458,e)
    end subroutine
    subroutine CommitNativeRestart(total0,balance,cumbottom,cumtop)
        real(KREAL),intent(in)::total0(4),balance(4),cumbottom(4),cumtop(4)
        integer::fu,ios,rc
        integer(int64)::hashes(2),meta_hash
        logical::exists
        if(myID/=0)return
        hashes=[file_hash('flow2d.dat.next'//c_null_char),file_hash('vdf2d.dat.next'//c_null_char)]
        call require(all(hashes/=0),'staged flow/vdf files missing')
        open(newunit=fu,file='s2_restart.meta.next',status='replace',access='stream',form='unformatted',iostat=ios)
        call require(ios==0,'metadata staging open failed')
        write(fu,iostat=ios)schema,int([storage_size(0.0_KREAL),Nx_all,Ny_all,Nvel,ck,flag_vislaw],int32), &
            iter,time_dis,dt_global,physical_parameters(),option_identity(),TimeIntegration_method, &
            PI,bz,SMV,mu_ref,uVel,vVel,wVel,gridx_all,gridy_all,total0,balance,cumbottom,cumtop, &
            res,res_relative,res_absolute,res_num,res_den,res_rms,res_max,hashes,schema
        call require(ios==0,'metadata staging write failed')
        close(fu,iostat=ios);call require(ios==0,'metadata close failed')
        meta_hash=file_hash('s2_restart.meta.next'//c_null_char)
        open(newunit=fu,file='s2_restart.valid.next',status='replace',iostat=ios)
        call require(ios==0,'marker staging open failed')
        write(fu,*,iostat=ios)schema,meta_hash
        call require(ios==0,'marker staging write failed')
        close(fu,iostat=ios);call require(ios==0,'marker close failed')
        ! Invalidate old marker only after every next-generation file has closed successfully.
        inquire(file='s2_restart.valid',exist=exists)
        if(exists)then
            open(newunit=fu,file='s2_restart.valid',status='old',iostat=ios)
            call require(ios==0,'old marker invalidation open failed')
            close(fu,status='delete',iostat=ios);call require(ios==0,'old marker invalidation failed')
        endif
        rc=c_rename('flow2d.dat.next'//c_null_char,'flow2d.dat'//c_null_char);call require(rc==0,'flow publish failed')
        rc=c_rename('vdf2d.dat.next'//c_null_char,'vdf2d.dat'//c_null_char);call require(rc==0,'vdf publish failed')
        rc=c_rename('grid2d.grd.next'//c_null_char,'grid2d.grd'//c_null_char);call require(rc==0,'grid publish failed')
        rc=c_rename('s2_restart.meta.next'//c_null_char,'s2_restart.meta'//c_null_char);call require(rc==0,'metadata publish failed')
        rc=c_rename('s2_restart.valid.next'//c_null_char,'s2_restart.valid'//c_null_char);call require(rc==0,'marker publish failed')
    end subroutine
    subroutine ReadNativeMetadata(validate_only,total0,balance,cumbottom,cumtop)
        logical,intent(in)::validate_only
        real(KREAL),intent(out)::total0(4),balance(4),cumbottom(4),cumtop(4)
        integer::fu,ios,step
        integer(int32)::header(6),opts(12)
        integer(int64)::hashes(2),meta_hash
        character(len=16)::label,endmark
        character(len=24)::method
        character(len=1)::extra
        real(KREAL)::t,rdt,params(16),constants(4),norms(4,7)
        real(KREAL),allocatable::uv(:),vv(:),wt(:),gx(:,:),gy(:,:)
        logical::exists
        call require(len_trim(S2CheckpointImport)==0,'native restart and legacy import are mutually exclusive')
        inquire(file='s2_restart.valid',exist=exists);call require(exists,'valid marker missing')
        open(newunit=fu,file='s2_restart.valid',status='old',iostat=ios);call require(ios==0,'marker open failed')
        read(fu,*,iostat=ios)label,meta_hash;close(fu)
        call require(ios==0.and.label==schema,'marker schema invalid')
        call require(meta_hash==file_hash('s2_restart.meta'//c_null_char),'metadata generation checksum mismatch')
        allocate(uv(Nvel),vv(Nvel),wt(Nvel),gx(Nx_all+1,Ny_all+1),gy(Nx_all+1,Ny_all+1))
        open(newunit=fu,file='s2_restart.meta',status='old',access='stream',form='unformatted',iostat=ios)
        call require(ios==0,'metadata open failed')
        read(fu,iostat=ios)label,header,step,t,rdt,params,opts,method,constants,uv,vv,wt,gx,gy, &
            total0,balance,cumbottom,cumtop,norms,hashes,endmark
        call require(ios==0,'metadata truncated')
        read(fu,iostat=ios)extra;close(fu);call require(ios<0,'metadata trailing bytes')
        call require(label==schema.and.endmark==schema.and.all(header== &
            int([storage_size(0.0_KREAL),Nx_all,Ny_all,Nvel,ck,flag_vislaw],int32)),'schema/dimension/model identity mismatch')
        call require(rdt==dt_global.and.all(params==physical_parameters()).and.all(opts==option_identity()) &
            .and.method==TimeIntegration_method.and.all(constants==[PI,bz,SMV,mu_ref]),'physical/config/constants identity mismatch')
        call require(all(uv==uVel).and.all(vv==vVel).and.all(wt==wVel), &
            'velocity table identity mismatch')
        if(myID==0)call require(all(gx==gridx_all).and.all(gy==gridy_all),'grid identity mismatch')
        call require(all(hashes==[file_hash('flow2d.dat'//c_null_char),file_hash('vdf2d.dat'//c_null_char)]), &
            'flow/vdf generation checksum mismatch')
        call require(step>=0.and.step<MAX_ITER.and.t>=0.0_KREAL,'restart step/time or target MAX_ITER invalid')
        call require(all(ieee_is_finite(norms)).and.all(ieee_is_finite(total0)).and.all(ieee_is_finite(balance)) &
            .and.all(ieee_is_finite(cumbottom)).and.all(ieee_is_finite(cumtop)),'nonfinite ledger/residual metadata')
        if(validate_only)return
        iter=step;time_dis=t
        res=norms(:,1);res_relative=norms(:,2);res_absolute=norms(:,3)
        res_num=norms(:,4);res_den=norms(:,5);res_rms=norms(:,6);res_max=norms(:,7)
        if(myID==0)write(*,'(A,I0,A,ES24.16)')'Native restart accepted step=',iter,' time=',time_dis
    end subroutine
end module

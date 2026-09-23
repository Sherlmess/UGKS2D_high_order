module mod_ImplicitBoundary
    use flow_data
    use Tools, only: DiscreteMaxwell_LocalVec, LocalFrame
    use mod_ImplicitState, only: mirror_u_index, mirror_v_index, &
        mirror_u_valid, mirror_v_valid
    use mpi
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: BuildStructuredMirrorMaps, GetMacroBoundaryDelta
    public :: GetSymmetryMicroDelta, GetSymmetryMicroValue, ComputeDiffuseWallFaceDelta
    public :: ComputeMaxwellWallFaceDelta
    public :: ResolveImplicitBoundaryType, IsAxisAlignedSymmetryNormal
    public :: ValidateImplicitSymmetryFaces, IsSupportedMaxwellFaceNormal
    public :: ValidateImplicitMaxwellFaces
contains
    subroutine BuildStructuredMirrorMaps()
        integer(KINT) :: i, j, k, ku, kv
        real(KREAL) :: tolerance

        if (Nu*Nv/=Nvel) call AbortMirrorMap('Nu*Nv does not equal Nvel')
        if (allocated(mirror_u_index)) then
            if (size(mirror_u_index)/=Nvel) deallocate(mirror_u_index, mirror_v_index)
        endif
        if (.not.allocated(mirror_u_index)) allocate(mirror_u_index(Nvel), mirror_v_index(Nvel))
        mirror_u_valid=.true.
        mirror_v_valid=.true.

        do j=1,Nv
            do i=1,Nu
                k=i+(j-1)*Nu
                ku=(Nu-i+1)+(j-1)*Nu
                kv=i+(Nv-j)*Nu
                mirror_u_index(k)=ku
                mirror_v_index(k)=kv

                tolerance=100.0_KREAL*epsilon(1.0_KREAL)* &
                    max(1.0_KREAL,abs(wVel(k)),abs(wVel(ku)))
                if (abs(wVel(ku)-wVel(k))>tolerance) mirror_u_valid=.false.
                tolerance=100.0_KREAL*epsilon(1.0_KREAL)* &
                    max(1.0_KREAL,abs(wVel(k)),abs(wVel(kv)))
                if (abs(wVel(kv)-wVel(k))>tolerance) mirror_v_valid=.false.

                tolerance=100.0_KREAL*epsilon(1.0_KREAL)*max(1.0_KREAL,abs(uVel(k)))
                if (abs(uVel(ku)+uVel(k))>tolerance) mirror_u_valid=.false.
                tolerance=100.0_KREAL*epsilon(1.0_KREAL)*max(1.0_KREAL,abs(vVel(k)))
                if (abs(vVel(ku)-vVel(k))>tolerance) mirror_u_valid=.false.
                if (abs(vVel(kv)+vVel(k))>tolerance) mirror_v_valid=.false.
                tolerance=100.0_KREAL*epsilon(1.0_KREAL)*max(1.0_KREAL,abs(uVel(k)))
                if (abs(uVel(kv)-uVel(k))>tolerance) mirror_v_valid=.false.
            enddo
        enddo
        do k=1,Nvel
            if (mirror_u_index(mirror_u_index(k))/=k) mirror_u_valid=.false.
            if (mirror_v_index(mirror_v_index(k))/=k) mirror_v_valid=.false.
        enddo
    end subroutine BuildStructuredMirrorMaps

    integer(KINT) function ResolveImplicitBoundaryType(i,j,dirc,side,configured_type) result(boundary_type)
        integer(KINT), intent(in) :: i,j,dirc,side,configured_type

        boundary_type=configured_type
        if (BG_geometry==BGconst_flatplate .and. dirc==JDIRC .and. side==Lside .and. &
            j==1 .and. ctr(i,1)%x<0.0_KREAL) boundary_type=BCconst_symmetry
    end function ResolveImplicitBoundaryType

    pure logical function IsAxisAlignedSymmetryNormal(nx_face,ny_face)
        real(KREAL), intent(in) :: nx_face,ny_face
        real(KREAL) :: tolerance

        tolerance=100.0_KREAL*epsilon(1.0_KREAL)
        IsAxisAlignedSymmetryNormal= &
            (abs(abs(nx_face)-1.0_KREAL)<=tolerance .and. abs(ny_face)<=tolerance) .or. &
            (abs(abs(ny_face)-1.0_KREAL)<=tolerance .and. abs(nx_face)<=tolerance)
    end function IsAxisAlignedSymmetryNormal

    logical function HasStructuredMirrorForNormal(nx_face,ny_face)
        real(KREAL), intent(in) :: nx_face,ny_face
        real(KREAL) :: tolerance

        tolerance=100.0_KREAL*epsilon(1.0_KREAL)
        HasStructuredMirrorForNormal=.false.
        if (abs(abs(nx_face)-1.0_KREAL)<=tolerance .and. abs(ny_face)<=tolerance) then
            HasStructuredMirrorForNormal=mirror_u_valid
        elseif (abs(abs(ny_face)-1.0_KREAL)<=tolerance .and. abs(nx_face)<=tolerance) then
            HasStructuredMirrorForNormal=mirror_v_valid
        endif
    end function HasStructuredMirrorForNormal

    logical function IsSupportedMaxwellFaceNormal(nx_face,ny_face,accommodation)
        real(KREAL), intent(in) :: nx_face,ny_face,accommodation
        real(KREAL) :: endpoint_tolerance

        IsSupportedMaxwellFaceNormal=ieee_is_finite(accommodation) .and. &
            accommodation>=0.0_KREAL .and. accommodation<=1.0_KREAL
        if (.not.IsSupportedMaxwellFaceNormal) return
        endpoint_tolerance=100.0_KREAL*epsilon(1.0_KREAL)
        if (abs(accommodation-1.0_KREAL)<=endpoint_tolerance) return
        IsSupportedMaxwellFaceNormal=HasStructuredMirrorForNormal(nx_face,ny_face)
    end function IsSupportedMaxwellFaceNormal

    subroutine ValidateImplicitSymmetryFaces()
        integer(KINT) :: i,j,boundary_type
        logical :: local_ok,global_ok

        local_ok=.true.
        if (IP_X==0) then
            do j=1,ny
                boundary_type=ResolveImplicitBoundaryType(1,j,IDIRC,Lside,BC_Type(1))
                if (boundary_type==BCconst_symmetry) local_ok=local_ok .and. &
                    HasStructuredMirrorForNormal(-vface(1,j)%cosx,-vface(1,j)%cosy)
            enddo
        endif
        if (IP_X==NP_X-1) then
            do j=1,ny
                boundary_type=ResolveImplicitBoundaryType(nx,j,IDIRC,Rside,BC_Type(2))
                if (boundary_type==BCconst_symmetry) local_ok=local_ok .and. &
                    HasStructuredMirrorForNormal(vface(nx+1,j)%cosx,vface(nx+1,j)%cosy)
            enddo
        endif
        if (IP_Y==0) then
            do i=1,nx
                boundary_type=ResolveImplicitBoundaryType(i,1,JDIRC,Lside,BC_Type(3))
                if (boundary_type==BCconst_symmetry) local_ok=local_ok .and. &
                    HasStructuredMirrorForNormal(-hface(i,1)%cosx,-hface(i,1)%cosy)
            enddo
        endif
        if (IP_Y==NP_Y-1) then
            do i=1,nx
                boundary_type=ResolveImplicitBoundaryType(i,ny,JDIRC,Rside,BC_Type(4))
                if (boundary_type==BCconst_symmetry) local_ok=local_ok .and. &
                    HasStructuredMirrorForNormal(hface(i,ny+1)%cosx,hface(i,ny+1)%cosy)
            enddo
        endif
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global_ok) call AbortMirrorMap('implicit symmetry face lacks its required velocity mirror')
    end subroutine ValidateImplicitSymmetryFaces

    subroutine ValidateImplicitMaxwellFaces()
        integer(KINT) :: i,j,boundary_type
        logical :: local_ok,global_ok

        if (.not.IsImplicitSteady() .and. .not.IsImplicitUnsteady()) return
        if (FluxWall_type/=FluxWall_Maxwell) return
        local_ok=.true.
        if (IP_X==0) then
            do j=1,ny
                boundary_type=ResolveImplicitBoundaryType(1,j,IDIRC,Lside,BC_Type(1))
                if (boundary_type==BCconst_wall_isothermal) local_ok=local_ok .and. &
                    IsSupportedMaxwellFaceNormal(-vface(1,j)%cosx,-vface(1,j)%cosy,ACC)
            enddo
        endif
        if (IP_X==NP_X-1) then
            do j=1,ny
                boundary_type=ResolveImplicitBoundaryType(nx,j,IDIRC,Rside,BC_Type(2))
                if (boundary_type==BCconst_wall_isothermal) local_ok=local_ok .and. &
                    IsSupportedMaxwellFaceNormal(vface(nx+1,j)%cosx,vface(nx+1,j)%cosy,ACC)
            enddo
        endif
        if (IP_Y==0) then
            do i=1,nx
                boundary_type=ResolveImplicitBoundaryType(i,1,JDIRC,Lside,BC_Type(3))
                if (boundary_type==BCconst_wall_isothermal) local_ok=local_ok .and. &
                    IsSupportedMaxwellFaceNormal(-hface(i,1)%cosx,-hface(i,1)%cosy,ACC)
            enddo
        endif
        if (IP_Y==NP_Y-1) then
            do i=1,nx
                boundary_type=ResolveImplicitBoundaryType(i,ny,JDIRC,Rside,BC_Type(4))
                if (boundary_type==BCconst_wall_isothermal) local_ok=local_ok .and. &
                    IsSupportedMaxwellFaceNormal(hface(i,ny+1)%cosx,hface(i,ny+1)%cosy,ACC)
            enddo
        endif
        call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
        if (.not.global_ok) call AbortMirrorMap('implicit Maxwell face lacks its required velocity mirror')
    end subroutine ValidateImplicitMaxwellFaces

    subroutine AbortMirrorMap(message)
        character(len=*), intent(in) :: message

        write(*,'(A,I0,2A)') 'Invalid structured velocity mirror map on rank ', myID, ': ', trim(message)
        call MPI_Abort(MPI_COMM_WORLD, 202, ierr)
    end subroutine AbortMirrorMap

    subroutine GetMacroBoundaryDelta(boundary_type,nx_face,ny_face,delta_inner,delta_boundary,ok)
        integer(KINT), intent(in) :: boundary_type
        real(KREAL), intent(in) :: nx_face, ny_face, delta_inner(4)
        real(KREAL), intent(out) :: delta_boundary(4)
        logical, intent(out) :: ok
        real(KREAL) :: normal_momentum

        ok=.true.
        select case(boundary_type)
        case(BCconst_freestream,BCconst_inflow)
            delta_boundary=0.0_KREAL
        case(BCconst_extrapolation)
            delta_boundary=delta_inner
        case(BCconst_symmetry,BCconst_wall_isothermal)
            delta_boundary=delta_inner
            normal_momentum=delta_inner(2)*nx_face+delta_inner(3)*ny_face
            delta_boundary(2)=delta_inner(2)-2.0_KREAL*normal_momentum*nx_face
            delta_boundary(3)=delta_inner(3)-2.0_KREAL*normal_momentum*ny_face
        case default
            delta_boundary=0.0_KREAL
            ok=.false.
        end select
    end subroutine GetMacroBoundaryDelta

    subroutine GetSymmetryMicroDelta(delta_inner,nx_face,ny_face,delta_boundary,ok)
        real(KREAL), intent(in) :: delta_inner(:), nx_face, ny_face
        real(KREAL), intent(out) :: delta_boundary(:)
        logical, intent(out) :: ok
        real(KREAL) :: tolerance

        ok=.false.
        delta_boundary=0.0_KREAL
        if (.not.allocated(mirror_u_index) .or. size(delta_inner)/=Nvel .or. &
            size(delta_boundary)/=Nvel) return
        tolerance=100.0_KREAL*epsilon(1.0_KREAL)
        if (abs(abs(nx_face)-1.0_KREAL)<=tolerance .and. abs(ny_face)<=tolerance) then
            if (.not.mirror_u_valid) return
            delta_boundary=delta_inner(mirror_u_index)
        elseif (abs(abs(ny_face)-1.0_KREAL)<=tolerance .and. abs(nx_face)<=tolerance) then
            if (.not.mirror_v_valid) return
            delta_boundary=delta_inner(mirror_v_index)
        else
            return
        endif
        ok=.true.
    end subroutine GetSymmetryMicroDelta

    subroutine GetSymmetryMicroValue(delta_inner,nx_face,ny_face,k,delta_boundary,ok)
        real(KREAL), intent(in) :: delta_inner(:),nx_face,ny_face
        integer(KINT), intent(in) :: k
        real(KREAL), intent(out) :: delta_boundary
        logical, intent(out) :: ok
        real(KREAL) :: tolerance

        ok=.false.
        delta_boundary=0.0_KREAL
        if (.not.allocated(mirror_u_index) .or. size(delta_inner)/=Nvel) return
        if (k<1 .or. k>Nvel) return
        tolerance=100.0_KREAL*epsilon(1.0_KREAL)
        if (abs(abs(nx_face)-1.0_KREAL)<=tolerance .and. abs(ny_face)<=tolerance) then
            if (.not.mirror_u_valid) return
            delta_boundary=delta_inner(mirror_u_index(k))
        elseif (abs(abs(ny_face)-1.0_KREAL)<=tolerance .and. abs(nx_face)<=tolerance) then
            if (.not.mirror_v_valid) return
            delta_boundary=delta_inner(mirror_v_index(k))
        else
            return
        endif
        ok=.true.
    end subroutine GetSymmetryMicroValue

    subroutine ComputeDiffuseWallFaceDelta(prim_wall_input,nx_face,ny_face,side, &
        delta_h_inner,delta_b_inner,delta_h_face,delta_b_face,delta_rho,ok)
        real(KREAL), intent(in) :: prim_wall_input(4), nx_face, ny_face
        integer(KINT), intent(in) :: side
        real(KREAL), intent(in) :: delta_h_inner(:), delta_b_inner(:)
        real(KREAL), intent(out) :: delta_h_face(:), delta_b_face(:)
        real(KREAL), intent(out) :: delta_rho
        logical, intent(out) :: ok
        real(KREAL) :: prim_local(4), incidence, reflection, outgoing_flux_scale
        real(KREAL) :: vn(Nvel), vt(Nvel), delta_selector(Nvel)
        real(KREAL) :: h_unit(Nvel), b_unit(Nvel)
        integer(KINT) :: k

        ok=.false.; delta_rho=0.0_KREAL
        delta_h_face=0.0_KREAL; delta_b_face=0.0_KREAL
        if (size(delta_h_inner)/=Nvel .or. size(delta_b_inner)/=Nvel .or. &
            size(delta_h_face)/=Nvel .or. size(delta_b_face)/=Nvel) return
        if (side/=Lside .and. side/=Rside) return

        do k=1,Nvel
            vn(k)=uVel(k)*nx_face+vVel(k)*ny_face
            vt(k)=-uVel(k)*ny_face+vVel(k)*nx_face
            if (vn(k)<0.0_KREAL) then
                delta_selector(k)=1.0_KREAL
            else
                delta_selector(k)=0.0_KREAL
            endif
        enddo

        prim_local=LocalFrame(prim_wall_input,nx_face,ny_face)
        prim_local(1)=1.0_KREAL
        call DiscreteMaxwell_LocalVec(h_unit,b_unit,vn,vt,prim_local)
        incidence=sum(wVel*vn*delta_h_inner*(1.0_KREAL-delta_selector))
        reflection=sum(wVel*vn*h_unit*delta_selector)
        outgoing_flux_scale=sum(abs(wVel*vn*h_unit*delta_selector))
        if (abs(reflection)<=1000.0_KREAL*epsilon(1.0_KREAL)* &
            max(outgoing_flux_scale,tiny(1.0_KREAL))) return

        delta_rho=-incidence/reflection
        delta_h_face=delta_rho*h_unit*delta_selector+ &
                     delta_h_inner*(1.0_KREAL-delta_selector)
        delta_b_face=delta_rho*b_unit*delta_selector+ &
                     delta_b_inner*(1.0_KREAL-delta_selector)
        ok=.true.
    end subroutine ComputeDiffuseWallFaceDelta

    subroutine ComputeMaxwellWallFaceDelta(prim_wall_input,nx_face,ny_face,side,accommodation, &
        delta_h_inner,delta_b_inner,delta_h_face,delta_b_face,delta_rho,ok)
        real(KREAL), intent(in) :: prim_wall_input(4),nx_face,ny_face,accommodation
        integer(KINT), intent(in) :: side
        real(KREAL), intent(in) :: delta_h_inner(:),delta_b_inner(:)
        real(KREAL), intent(out) :: delta_h_face(:),delta_b_face(:),delta_rho
        logical, intent(out) :: ok
        real(KREAL) :: diffuse_h(Nvel),diffuse_b(Nvel),endpoint_tolerance
        real(KREAL) :: specular_h(Nvel),specular_b(Nvel),vn(Nvel)
        logical :: mirror_ok
        integer(KINT) :: k

        ok=.false.
        delta_rho=0.0_KREAL
        delta_h_face=0.0_KREAL
        delta_b_face=0.0_KREAL
        if (.not.ieee_is_finite(accommodation)) return
        if (accommodation<0.0_KREAL .or. accommodation>1.0_KREAL) return

        call ComputeDiffuseWallFaceDelta(prim_wall_input,nx_face,ny_face,side, &
            delta_h_inner,delta_b_inner,diffuse_h,diffuse_b,delta_rho,ok)
        if (.not.ok) return
        endpoint_tolerance=100.0_KREAL*epsilon(1.0_KREAL)
        if (abs(accommodation-1.0_KREAL)<=endpoint_tolerance) then
            delta_h_face=diffuse_h
            delta_b_face=diffuse_b
            return
        endif

        call GetSymmetryMicroDelta(delta_h_inner,nx_face,ny_face,specular_h,mirror_ok)
        if (.not.mirror_ok) then
            ok=.false.
            return
        endif
        call GetSymmetryMicroDelta(delta_b_inner,nx_face,ny_face,specular_b,mirror_ok)
        if (.not.mirror_ok) then
            ok=.false.
            return
        endif

        do k=1,Nvel
            vn(k)=uVel(k)*nx_face+vVel(k)*ny_face
        enddo
        delta_h_face=delta_h_inner
        delta_b_face=delta_b_inner
        where(vn<0.0_KREAL)
            delta_h_face=accommodation*diffuse_h+(1.0_KREAL-accommodation)*specular_h
            delta_b_face=accommodation*diffuse_b+(1.0_KREAL-accommodation)*specular_b
        endwhere
        ok=.true.
    end subroutine ComputeMaxwellWallFaceDelta
end module mod_ImplicitBoundary

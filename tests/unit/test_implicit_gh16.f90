program test_implicit_gh16
    use mpi
    use flow_data
    use InitializeFlow, only: InitVelocityGaussHermite, FlattenVelocitySpace
    use Tools, only: DiscreteMaxwell_Vec
    use mod_ImplicitState, only: mirror_u_index, mirror_v_index
    use mod_ImplicitBoundary, only: BuildStructuredMirrorMaps, ComputeDiffuseWallFaceDelta, &
        ComputeMaxwellWallFaceDelta
    use mod_ImplicitCommon, only: MomentConserved, ProjectToConservedMoments
    use TestSupport
    use, intrinsic :: ieee_arithmetic
    implicit none

    integer :: ierr_local
    integer(KINT) :: i,j,k,ku,kv,expected_u,expected_v,normal_index,accommodation_index
    logical :: ok
    real(KREAL) :: target(4),moment(4),relative_error,mass_flux,delta_rho
    real(KREAL) :: accommodation,maxwell_mass_flux
    real(KREAL) :: prim_reference(4),prim_wall_test(4)
    real(KREAL), allocatable :: reference_h(:),reference_b(:),h(:),b(:)
    real(KREAL), allocatable :: delta_h_inner(:),delta_b_inner(:),delta_h_face(:),delta_b_face(:)
    real(KREAL), dimension(2,4) :: normals
    integer(KINT), dimension(4) :: sides
    real(KREAL), parameter :: accommodations(5)=[0.0_KREAL,0.25_KREAL,0.5_KREAL,0.75_KREAL,1.0_KREAL]

    call MPI_Init(ierr_local)
    call MPI_Comm_rank(MPI_COMM_WORLD,myID,ierr_local)

    Nu=16; Nv=16; Nvel=Nu*Nv; CK=1
    allocate(uSpace(Nvel),vSpace(Nvel),weight(Nvel))
    allocate(uVel(Nvel),vVel(Nvel),wVel(Nvel))
    call InitVelocityGaussHermite(Nu,Nv)
    call FlattenVelocitySpace()

    call AssertTrue(Nu==16 .and. Nv==16 .and. Nvel==256, 'GH16 velocity dimensions')
    call AssertTrue(all(ieee_is_finite(uVel)) .and. all(ieee_is_finite(vVel)) .and. &
                    all(ieee_is_finite(wVel)), 'GH16 velocity arrays finite')
    call AssertTrue(all(wVel>0.0_KREAL), 'GH16 quadrature weights positive')

    do j=1,Nv
        do i=1,Nu
            k=i+(j-1)*Nu
            ku=(Nu-i+1)+(j-1)*Nu
            kv=i+(Nv-j)*Nu
            call AssertTrue(uVel(k)*uVel(ku)<0.0_KREAL .and. &
                            abs(uVel(k)+uVel(ku))<1.0e-14_KREAL, 'GH16 u nodes are paired')
            call AssertTrue(vVel(k)*vVel(kv)<0.0_KREAL .and. &
                            abs(vVel(k)+vVel(kv))<1.0e-14_KREAL, 'GH16 v nodes are paired')
            call AssertTrue(abs(vVel(k)-vVel(ku))<1.0e-14_KREAL .and. &
                            abs(uVel(k)-uVel(kv))<1.0e-14_KREAL, 'GH16 tangential nodes are preserved')
        enddo
    enddo

    call BuildStructuredMirrorMaps()
    do k=1,Nvel
        i=mod(k-1,Nu)+1
        j=(k-1)/Nu+1
        expected_u=(Nu-i+1)+(j-1)*Nu
        expected_v=i+(Nv-j)*Nu
        call AssertTrue(mirror_u_index(k)==expected_u .and. mirror_v_index(k)==expected_v, &
                        'GH16 mirror indices match structured grid')
        call AssertTrue(mirror_u_index(mirror_u_index(k))==k .and. &
                        mirror_v_index(mirror_v_index(k))==k, 'GH16 mirror indices are involutions')
        call AssertTrue(abs(wVel(k)-wVel(mirror_u_index(k)))<1.0e-14_KREAL, &
            'GH16 u mirror preserves weight')
        call AssertTrue(abs(wVel(k)-wVel(mirror_v_index(k)))<1.0e-14_KREAL, &
            'GH16 v mirror preserves weight')
    enddo

    allocate(reference_h(Nvel),reference_b(Nvel),h(Nvel),b(Nvel))
    prim_reference=[1.0_KREAL,0.15_KREAL,-0.10_KREAL,1.0_KREAL]
    call DiscreteMaxwell_Vec(reference_h,reference_b,prim_reference)
    call MomentConserved(reference_h,reference_b,uVel,vVel,wVel,target)
    h=reference_h*(1.0_KREAL+1.0e-3_KREAL*(1.0_KREAL+0.1_KREAL*uVel))
    b=reference_b*(1.0_KREAL+2.0e-3_KREAL*(1.0_KREAL-0.05_KREAL*vVel))
    call ProjectToConservedMoments(h,b,target,uVel,vVel,wVel,reference_h,reference_b,ok)
    call MomentConserved(h,b,uVel,vVel,wVel,moment)
    relative_error=maxval(abs(moment-target))/max(1.0_KREAL,maxval(abs(target)))
    call AssertTrue(ok, 'GH16 moment projection accepted')
    call AssertTrue(relative_error<=1.0e-10_KREAL, 'GH16 moment projection relative error')

    allocate(delta_h_inner(Nvel),delta_b_inner(Nvel),delta_h_face(Nvel),delta_b_face(Nvel))
    delta_h_inner=0.2_KREAL*(1.0_KREAL+0.01_KREAL*uVel+0.01_KREAL*vVel)
    delta_b_inner=0.1_KREAL*(1.0_KREAL-0.01_KREAL*uVel+0.01_KREAL*vVel)
    prim_wall_test=[1.0_KREAL,0.0_KREAL,0.0_KREAL,1.0_KREAL]
    normals(:,1)=[-1.0_KREAL,0.0_KREAL]; sides(1)=Lside
    normals(:,2)=[ 1.0_KREAL,0.0_KREAL]; sides(2)=Rside
    normals(:,3)=[0.0_KREAL,-1.0_KREAL]; sides(3)=Lside
    normals(:,4)=[0.0_KREAL, 1.0_KREAL]; sides(4)=Rside
    do normal_index=1,4
        call ComputeDiffuseWallFaceDelta(prim_wall_test,normals(1,normal_index),normals(2,normal_index), &
            sides(normal_index),delta_h_inner,delta_b_inner,delta_h_face,delta_b_face,delta_rho,ok)
        mass_flux=sum(wVel*(uVel*normals(1,normal_index)+vVel*normals(2,normal_index))*delta_h_face)
        call AssertTrue(ok, 'GH16 diffuse wall correction accepted')
        call AssertTrue(abs(mass_flux)<1.0e-10_KREAL, 'GH16 diffuse wall zero mass flux')
    enddo

    do normal_index=1,4
        do accommodation_index=1,size(accommodations)
            accommodation=accommodations(accommodation_index)
            call ComputeMaxwellWallFaceDelta(prim_wall_test,normals(1,normal_index),normals(2,normal_index), &
                sides(normal_index),accommodation,delta_h_inner,delta_b_inner, &
                delta_h_face,delta_b_face,delta_rho,ok)
            maxwell_mass_flux=sum(wVel*(uVel*normals(1,normal_index)+vVel*normals(2,normal_index))*delta_h_face)
            call AssertTrue(ok,'GH16 Maxwell wall correction accepted')
            call AssertTrue(abs(maxwell_mass_flux)<1.0e-10_KREAL,'GH16 Maxwell wall zero mass flux')
        enddo
    enddo

    call FinishTests('implicit_gh16')
    call MPI_Finalize(ierr_local)
end program test_implicit_gh16

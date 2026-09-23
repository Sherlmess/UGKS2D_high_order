program test_implicit_boundary
    use mpi
    use flow_data
    use control_para, only: TimeIntegration_method, IsImplicitSteady, IsImplicitUnsteady, FluxWall_type
    use mpi_para, only: NP_X, NP_Y, IP_X, IP_Y
    use Boundary, only: FreeWaveNumberUGKS
    use mod_ImplicitState
    use mod_ImplicitBoundary
    use TestSupport
    implicit none
    integer :: ierr_local,normal_index
    integer(KINT) :: accommodation_index,k
    logical :: ok
    real(KREAL) :: delta_boundary(4), delta_rho, mass_delta,accommodation,tolerance
    real(KREAL) :: symmetry_value
    real(KREAL) :: expected_wave_number
    real(KREAL) :: delta_h_inner(4), delta_b_inner(4)
    real(KREAL) :: delta_h_face(4), delta_b_face(4), symmetry_delta(4)
    real(KREAL) :: diffuse_h(4),diffuse_b(4),specular_h(4),specular_b(4)
    real(KREAL) :: expected_h(4),expected_b(4),maxwell_h(4),maxwell_b(4)
    real(KREAL) :: normals(2,4),normal_velocity(4),upwind_mass_flux
    real(KREAL), parameter :: accommodations(5)=[0.0_KREAL,0.25_KREAL,0.5_KREAL,0.75_KREAL,1.0_KREAL]
    integer(KINT) :: sides(4)

    call MPI_Init(ierr_local)
    nx=1; ny=1; Nu=2; Nv=2; Nvel=4; ck=1
    allocate(uVel(Nvel), vVel(Nvel), wVel(Nvel))
    allocate(ctr(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vface(1:nx+1,1:ny),hface(1:nx,1:ny+1))
    uVel=[-1.0_KREAL,1.0_KREAL,-1.0_KREAL,1.0_KREAL]
    vVel=[-1.0_KREAL,-1.0_KREAL,1.0_KREAL,1.0_KREAL]
    wVel=0.25_KREAL
    call AllocateImplicitState()

    ! All prescribed disturbance frequencies use the UGKS thermal-speed time
    ! scale.  Check the matching dispersion relation for every free-wave type.
    expected_wave_number=0.8_KREAL/(sqrt(1.4_KREAL/2.0_KREAL)*(0.5_KREAL+1.0_KREAL))
    call AssertClose(FreeWaveNumberUGKS('F',0.8_KREAL,0.0_KREAL,0.0_KREAL,0.5_KREAL,1.4_KREAL), &
        expected_wave_number,1.0e-14_KREAL,'fast wave uses UGKS frequency scale')
    expected_wave_number=0.8_KREAL/(sqrt(1.4_KREAL/2.0_KREAL)*(0.5_KREAL-1.0_KREAL))
    call AssertClose(FreeWaveNumberUGKS('S',0.8_KREAL,0.0_KREAL,0.0_KREAL,0.5_KREAL,1.4_KREAL), &
        expected_wave_number,1.0e-14_KREAL,'slow wave uses UGKS frequency scale')
    expected_wave_number=0.8_KREAL/(sqrt(1.4_KREAL/2.0_KREAL)*0.5_KREAL)
    call AssertClose(FreeWaveNumberUGKS('E',0.8_KREAL,0.0_KREAL,0.0_KREAL,0.5_KREAL,1.4_KREAL), &
        expected_wave_number,1.0e-14_KREAL,'entropy wave uses UGKS frequency scale')
    call AssertClose(FreeWaveNumberUGKS('V',0.8_KREAL,0.0_KREAL,0.0_KREAL,0.5_KREAL,1.4_KREAL), &
        expected_wave_number,1.0e-14_KREAL,'vorticity wave uses UGKS frequency scale')

    call BuildStructuredMirrorMaps()
    BG_geometry=BGconst_flatplate
    ctr(1,1)%x=-0.1_KREAL
    call AssertTrue(ResolveImplicitBoundaryType(1,1,JDIRC,Lside,BCconst_wall_isothermal)== &
        BCconst_symmetry,'flat-plate pre-tip implicit symmetry')
    ctr(1,1)%x=0.1_KREAL
    call AssertTrue(ResolveImplicitBoundaryType(1,1,JDIRC,Lside,BCconst_wall_isothermal)== &
        BCconst_wall_isothermal,'flat-plate post-tip implicit wall')
    call AssertTrue(IsAxisAlignedSymmetryNormal(0.0_KREAL,-1.0_KREAL), &
        'axis-aligned symmetry normal accepted')
    call AssertTrue(.not.IsAxisAlignedSymmetryNormal(0.5_KREAL,sqrt(0.75_KREAL)), &
        'oblique symmetry normal rejected')
    call AssertTrue(IsSupportedMaxwellFaceNormal(0.0_KREAL,-1.0_KREAL,0.0_KREAL), &
        'axis-aligned specular face accepted')
    call AssertTrue(IsSupportedMaxwellFaceNormal(1.0_KREAL,0.0_KREAL,0.5_KREAL), &
        'axis-aligned partial face accepted')
    call AssertTrue(IsSupportedMaxwellFaceNormal(sqrt(0.5_KREAL),sqrt(0.5_KREAL),1.0_KREAL), &
        'oblique diffuse face accepted')
    call AssertTrue(.not.IsSupportedMaxwellFaceNormal(sqrt(0.5_KREAL),sqrt(0.5_KREAL),0.5_KREAL), &
        'oblique partial face rejected')
    call AssertTrue(all(mirror_u_index==[2,1,4,3]), 'u mirror map')
    call AssertTrue(all(mirror_v_index==[3,4,1,2]), 'v mirror map')
    call GetMacroBoundaryDelta(BCconst_symmetry,1.0_KREAL,0.0_KREAL, &
        [1.0_KREAL,2.0_KREAL,3.0_KREAL,4.0_KREAL],delta_boundary,ok)
    call AssertTrue(ok .and. maxval(abs(delta_boundary- &
        [1.0_KREAL,-2.0_KREAL,3.0_KREAL,4.0_KREAL]))<1.0e-13_KREAL, &
        'macro reflection')
    call GetSymmetryMicroDelta([1.0_KREAL,2.0_KREAL,3.0_KREAL,4.0_KREAL], &
                               1.0_KREAL,0.0_KREAL,symmetry_delta,ok)
    call AssertTrue(ok .and. all(symmetry_delta==[2.0_KREAL,1.0_KREAL,4.0_KREAL,3.0_KREAL]), &
                    'micro x reflection')

    Prim_wall=[1.0_KREAL,0.0_KREAL,0.0_KREAL,1.0_KREAL]
    delta_h_inner=[0.4_KREAL,-0.2_KREAL,0.1_KREAL,0.3_KREAL]
    delta_b_inner=[0.2_KREAL,-0.1_KREAL,0.05_KREAL,0.15_KREAL]
    normals(:,1)=[-1.0_KREAL,0.0_KREAL]; sides(1)=Lside
    normals(:,2)=[ 1.0_KREAL,0.0_KREAL]; sides(2)=Rside
    normals(:,3)=[0.0_KREAL,-1.0_KREAL]; sides(3)=Lside
    normals(:,4)=[0.0_KREAL, 1.0_KREAL]; sides(4)=Rside
    do normal_index=1,4
        call GetSymmetryMicroDelta(delta_h_inner,normals(1,normal_index),normals(2,normal_index), &
            symmetry_delta,ok)
        call AssertTrue(ok,'symmetry vector oracle accepted')
        call GetSymmetryMicroDelta(delta_b_inner,normals(1,normal_index),normals(2,normal_index), &
            specular_b,ok)
        call AssertTrue(ok,'b symmetry vector oracle accepted')
        do k=1,Nvel
            call GetSymmetryMicroValue(delta_h_inner,normals(1,normal_index), &
                normals(2,normal_index),k,symmetry_value,ok)
            call AssertTrue(ok .and. symmetry_value==symmetry_delta(k), &
                'scalar symmetry lookup matches vector oracle')
            call GetSymmetryMicroValue(delta_b_inner,normals(1,normal_index), &
                normals(2,normal_index),k,symmetry_value,ok)
            call AssertTrue(ok .and. symmetry_value==specular_b(k), &
                'scalar b symmetry lookup matches vector oracle')
        enddo
    enddo
    call GetSymmetryMicroValue(delta_h_inner,0.0_KREAL,-1.0_KREAL,0,symmetry_value,ok)
    call AssertTrue(.not.ok,'scalar symmetry lookup rejects invalid velocity index')
    do normal_index=1,4
        call ComputeDiffuseWallFaceDelta(Prim_wall,normals(1,normal_index),normals(2,normal_index), &
            sides(normal_index),delta_h_inner,delta_b_inner,delta_h_face,delta_b_face,delta_rho,ok)
        normal_velocity=uVel*normals(1,normal_index)+vVel*normals(2,normal_index)
        mass_delta=sum(wVel*normal_velocity*delta_h_face)
        upwind_mass_flux=sum(wVel*(max(normal_velocity,0.0_KREAL)*delta_h_inner+ &
            min(normal_velocity,0.0_KREAL)*delta_h_face))
        call AssertTrue(ok .and. abs(mass_delta)<1.0e-12_KREAL, &
                        'diffuse wall face zero mass increment')
        call AssertTrue(ok .and. abs(upwind_mass_flux)<1.0e-12_KREAL, &
                        'diffuse wall upwind zero mass increment')
    enddo

    tolerance=1.0e-12_KREAL
    do normal_index=1,4
        call ComputeDiffuseWallFaceDelta(Prim_wall,normals(1,normal_index),normals(2,normal_index), &
            sides(normal_index),delta_h_inner,delta_b_inner,diffuse_h,diffuse_b,delta_rho,ok)
        call AssertTrue(ok,'diffuse oracle accepted for Maxwell test')
        call GetSymmetryMicroDelta(delta_h_inner,normals(1,normal_index),normals(2,normal_index),specular_h,ok)
        call AssertTrue(ok,'specular h oracle accepted for Maxwell test')
        call GetSymmetryMicroDelta(delta_b_inner,normals(1,normal_index),normals(2,normal_index),specular_b,ok)
        call AssertTrue(ok,'specular b oracle accepted for Maxwell test')
        normal_velocity=uVel*normals(1,normal_index)+vVel*normals(2,normal_index)
        do accommodation_index=1,size(accommodations)
            accommodation=accommodations(accommodation_index)
            call ComputeMaxwellWallFaceDelta(Prim_wall,normals(1,normal_index),normals(2,normal_index), &
                sides(normal_index),accommodation,delta_h_inner,delta_b_inner,maxwell_h,maxwell_b,delta_rho,ok)
            expected_h=delta_h_inner
            expected_b=delta_b_inner
            do k=1,Nvel
                if (normal_velocity(k)<0.0_KREAL) then
                    expected_h(k)=accommodation*diffuse_h(k)+(1.0_KREAL-accommodation)*specular_h(k)
                    expected_b(k)=accommodation*diffuse_b(k)+(1.0_KREAL-accommodation)*specular_b(k)
                endif
            enddo
            mass_delta=sum(wVel*normal_velocity*maxwell_h)
            call AssertTrue(ok,'Maxwell wall correction accepted')
            call AssertTrue(maxval(abs(maxwell_h-expected_h))<tolerance,'Maxwell h convex combination')
            call AssertTrue(maxval(abs(maxwell_b-expected_b))<tolerance,'Maxwell b convex combination')
            call AssertTrue(abs(mass_delta)<tolerance,'Maxwell wall zero mass correction')
        enddo
    enddo

    call ComputeMaxwellWallFaceDelta(Prim_wall,sqrt(0.5_KREAL),sqrt(0.5_KREAL),Lside,0.5_KREAL, &
        delta_h_inner,delta_b_inner,maxwell_h,maxwell_b,delta_rho,ok)
    call AssertTrue(.not.ok,'oblique partial-accommodation wall rejected')

    ! A shifted streamwise range is valid for a horizontal flat plate: the
    ! wall correction needs only the still-symmetric v reflection.
    uVel=[-1.0_KREAL,2.0_KREAL,-1.0_KREAL,2.0_KREAL]
    vVel=[-1.0_KREAL,-1.0_KREAL,1.0_KREAL,1.0_KREAL]
    call BuildStructuredMirrorMaps()
    call GetSymmetryMicroDelta(delta_h_inner,0.0_KREAL,-1.0_KREAL,symmetry_delta,ok)
    call AssertTrue(ok .and. maxval(abs(symmetry_delta-delta_h_inner([3,4,1,2])))<tolerance, &
        'horizontal reflection accepts asymmetric unused u grid')
    call GetSymmetryMicroDelta(delta_h_inner,-1.0_KREAL,0.0_KREAL,symmetry_delta,ok)
    call AssertTrue(.not.ok,'vertical reflection rejects unavailable u mirror')
    call AssertTrue(IsSupportedMaxwellFaceNormal(0.0_KREAL,-1.0_KREAL,0.6_KREAL), &
        'horizontal partial Maxwell accepts valid v mirror')
    call AssertTrue(.not.IsSupportedMaxwellFaceNormal(-1.0_KREAL,0.0_KREAL,0.6_KREAL), &
        'vertical partial Maxwell rejects unavailable u mirror')
    call ComputeMaxwellWallFaceDelta(Prim_wall,0.0_KREAL,-1.0_KREAL,Lside,0.6_KREAL, &
        delta_h_inner,delta_b_inner,maxwell_h,maxwell_b,delta_rho,ok)
    normal_velocity=vVel
    mass_delta=sum(wVel*normal_velocity*maxwell_h)
    call AssertTrue(ok .and. abs(mass_delta)<tolerance, &
        'horizontal partial Maxwell conserves mass with asymmetric unused u grid')

    ! The same Maxwell-face preflight applies to both implicit modes.  This
    ! fixture keeps the horizontal ACC=0.6 face valid while the vertical
    ! direction remains intentionally unavailable in the assertions above.
    BG_geometry=BGconst_general
    NP_X=1; NP_Y=1; IP_X=0; IP_Y=0
    BC_Type=BCconst_freestream
    BC_Type(3)=BCconst_wall_isothermal
    FluxWall_type=FluxWall_Maxwell
    ACC=0.6_KREAL
    vface%cosx=1.0_KREAL; vface%cosy=0.0_KREAL
    hface%cosx=0.0_KREAL; hface%cosy=1.0_KREAL
    TimeIntegration_method='ImplicitSteady'
    call ValidateImplicitMaxwellFaces()
    call AssertTrue(IsImplicitSteady(),'Maxwell face preflight dispatches in steady implicit mode')
    TimeIntegration_method='ImplicitUnsteady'
    call ValidateImplicitMaxwellFaces()
    call AssertTrue(IsImplicitUnsteady(),'Maxwell face preflight dispatches in unsteady implicit mode')

    call FinishTests('implicit_boundary')
    call MPI_Finalize(ierr_local)
end program test_implicit_boundary

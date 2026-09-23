program test_implicit_common
    use mpi
    use flow_data
    use tools
    use mod_ImplicitCommon
    use TestSupport
    implicit none
    integer :: ierr_local
    real(KREAL) :: w(4), flux(4), moment(4)
    real(KREAL) :: h(2), b(2), u(2), v(2), weight_local(2)
    real(KREAL) :: h_proj(4),b_proj(4),u_proj(4),v_proj(4),w_proj(4)
    real(KREAL) :: ref_h_proj(4),ref_b_proj(4),target(4),projected(4)
    logical :: ok

    call MPI_Init(ierr_local)
    gamma=1.4_KREAL
    w=GetConserved([1.0_KREAL,2.0_KREAL,-1.0_KREAL,0.5_KREAL])
    call EulerFluxNormal(w,1.0_KREAL,0.0_KREAL,flux)
    call AssertTrue(maxval(abs(flux-[2.0_KREAL,5.0_KREAL,-2.0_KREAL,12.0_KREAL])) &
                    <1.0e-13_KREAL, 'x Euler flux')
    call EulerFluxNormal(w,0.0_KREAL,1.0_KREAL,flux)
    call AssertTrue(maxval(abs(flux-[-1.0_KREAL,-2.0_KREAL,2.0_KREAL,-6.0_KREAL])) &
                    <1.0e-13_KREAL, 'y Euler flux')
    call AssertTrue(IsPhysicalState(w), 'physical conservative state')
    call AssertTrue(.not.IsPhysicalState([-1.0_KREAL,0.0_KREAL,0.0_KREAL,1.0_KREAL]), &
                    'negative density rejected')
    u=[-1.0_KREAL,1.0_KREAL]; v=0.0_KREAL; weight_local=0.5_KREAL
    h=[1.0_KREAL,3.0_KREAL]; b=[2.0_KREAL,4.0_KREAL]
    call MomentConserved(h,b,u,v,weight_local,moment)
    call AssertTrue(maxval(abs(moment-[2.0_KREAL,1.0_KREAL,0.0_KREAL,2.5_KREAL])) &
                    <1.0e-13_KREAL, 'kinetic moments')
    call AssertClose(PseudoTimeStep(0.1_KREAL,100.0_KREAL,0.5_KREAL),20.0_KREAL, &
                     1.0e-13_KREAL,'pseudo dt')
    u_proj=[-1.0_KREAL,1.0_KREAL,-1.0_KREAL,1.0_KREAL]
    v_proj=[-1.0_KREAL,-1.0_KREAL,1.0_KREAL,1.0_KREAL]; w_proj=0.25_KREAL
    h_proj=[0.2_KREAL,0.3_KREAL,0.4_KREAL,0.5_KREAL]
    b_proj=[0.1_KREAL,0.2_KREAL,0.3_KREAL,0.4_KREAL]
    ref_h_proj=1.0_KREAL; ref_b_proj=1.0_KREAL
    target=[1.2_KREAL,0.15_KREAL,-0.08_KREAL,1.7_KREAL]
    call ProjectToConservedMoments(h_proj,b_proj,target,u_proj,v_proj,w_proj, &
                                   ref_h_proj,ref_b_proj,ok)
    call MomentConserved(h_proj,b_proj,u_proj,v_proj,w_proj,projected)
    call AssertTrue(ok .and. maxval(abs(projected-target))<1.0e-12_KREAL, &
                    'moment-compatible projection')
    call FinishTests('implicit_common')
    call MPI_Finalize(ierr_local)
end program test_implicit_common

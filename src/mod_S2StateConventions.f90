! E1: representation-only interfaces for E2. No projection or state repair.
! Caller sets gamma=(CK+4)/(CK+2), as InitPara.set_const_para already does.
module mod_S2StateConventions
    use Precision, only: KREAL
    use Tools, only: GetPrimary,GetConserved,GetSoundSpeed
    implicit none
    private
    public :: NativeTemperatureFromTheta,ThetaFromNativeTemperature,LambdaFromTheta
    public :: S2Thermodynamics,S2ConservedFromTheta,S2PressureMatchedWall
contains
    pure real(KREAL) function NativeTemperatureFromTheta(theta)
        real(KREAL),intent(in)::theta
        NativeTemperatureFromTheta=2.0_KREAL*theta
    end function
    pure real(KREAL) function ThetaFromNativeTemperature(T_native)
        real(KREAL),intent(in)::T_native
        ThetaFromNativeTemperature=0.5_KREAL*T_native
    end function
    pure real(KREAL) function LambdaFromTheta(theta)
        real(KREAL),intent(in)::theta
        LambdaFromTheta=0.5_KREAL/theta
    end function
    subroutine S2Thermodynamics(w,prim,T_native,theta,p,e_internal,sound)
        real(KREAL),intent(in)::w(4)
        real(KREAL),intent(out)::prim(4),T_native,theta,p,e_internal,sound
        prim=GetPrimary(w)
        T_native=1.0_KREAL/prim(4)
        theta=ThetaFromNativeTemperature(T_native)
        p=prim(1)*theta
        e_internal=w(4)/w(1)-0.5_KREAL*(prim(2)**2+prim(3)**2)
        sound=GetSoundSpeed(prim)
    end subroutine
    function S2ConservedFromTheta(rho,ux,vy,theta) result(w)
        real(KREAL),intent(in)::rho,ux,vy,theta
        real(KREAL)::w(4)
        w=GetConserved([rho,ux,vy,LambdaFromTheta(theta)])
    end function
    subroutine S2PressureMatchedWall(w,ux,theta_w,prim_wall)
        ! Cartesian wall primitive. Reflection density remains the wall kernel's job.
        real(KREAL),intent(in)::w(4),ux,theta_w
        real(KREAL),intent(out)::prim_wall(4)
        real(KREAL)::prim(4),T_native,theta,p,e_internal,sound
        call S2Thermodynamics(w,prim,T_native,theta,p,e_internal,sound)
        prim_wall=[p/theta_w,ux,0.0_KREAL,LambdaFromTheta(theta_w)]
    end subroutine
end module

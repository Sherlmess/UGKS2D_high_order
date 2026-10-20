module mod_S2PhysicalGas
    use Precision, only: KREAL
    use flow_para, only: ck,S2PhysicalGas
    use Tools, only: GetTau
    use,intrinsic::ieee_arithmetic
    implicit none
    real(KREAL),save::gas_tau_min(3)=huge(1.0_KREAL),gas_tau_max(3)=0.0_KREAL
    integer,save::gas_tau_calls(3)=0
contains
    function s2_tau(w,site) result(tau)
        real(KREAL),intent(in)::w(4)
        integer,intent(in),optional::site
        integer::n
        real(KREAL)::tau,p,prim(4)
        tau=0.1_KREAL
        if(S2PhysicalGas==0)return
        p=2.0_KREAL/real(ck+2,KREAL)*(w(4)-0.5_KREAL*(w(2)**2+w(3)**2)/w(1))
        prim=[w(1),w(2)/w(1),w(3)/w(1),0.5_KREAL*w(1)/p]
        tau=GetTau(prim)
        n=1;if(present(site))n=site
        gas_tau_min(n)=min(gas_tau_min(n),tau);gas_tau_max(n)=max(gas_tau_max(n),tau)
        gas_tau_calls(n)=gas_tau_calls(n)+1
    end function
end module

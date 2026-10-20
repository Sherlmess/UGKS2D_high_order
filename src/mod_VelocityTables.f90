module mod_VelocityTables
    use Precision, only: KREAL
    use flow_para, only: S2VelocityTable,S2VelocityScale
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public::VelocityMode,VelocityCount,FillVelocity,VelocityMax
contains
    function VelocityMode() result(mode)
        character(len=16)::mode
        integer::st
        mode=S2VelocityTable;st=0
        if(len_trim(mode)==0)call get_environment_variable('S2_VELOCITY_TABLE',mode,status=st)
        if(st==1.or.len_trim(mode)==0)mode='full16'
        if(st/=0.and.st/=1)error stop 'invalid/truncated S2_VELOCITY_TABLE'
        select case(trim(mode))
        case('full16','half16','half32')
        case default;error stop 'S2_VELOCITY_TABLE must be full16/half16/half32'
        end select
    end function
    integer function VelocityCount()
        VelocityCount=16
        if(VelocityMode()=='half32')VelocityCount=32
    end function
    real(KREAL) function VelocityMax()
        real(KREAL),allocatable::u(:),v(:),w(:),rn(:),rw(:)
        integer::n
        n=VelocityCount();allocate(u(n*n),v(n*n),w(n*n),rn(n),rw(n))
        call FillVelocity(u,v,w,rn,rw);VelocityMax=maxval(abs(u))
    end function
    subroutine FillVelocity(u,v,weight,raw_node,raw_weight)
        real(KREAL),intent(out)::u(:),v(:),weight(:),raw_node(:),raw_weight(:)
        integer::i,j,k,n,st,ios
        character(len=64)::scale_env
        real(KREAL)::velocity_scale
        real(KREAL)::ue(size(raw_node)),we(size(raw_node))
        n=VelocityCount()
        if(size(raw_node)/=n.or.size(raw_weight)/=n.or.size(u)/=n*n.or.size(v)/=n*n.or.size(weight)/=n*n) &
            error stop 'velocity table shape mismatch'
        select case(trim(VelocityMode()))
        case('full16');call FillFull16(u,v,weight,raw_node,raw_weight);return
        case('half16')
            raw_node=[ &
                -3.68600716272439710e+00_KREAL, -2.86313388370807509e+00_KREAL, &
                -2.18392115309585799e+00_KREAL, -1.58885586227005504e+00_KREAL, &
                -1.06424631211622400e+00_KREAL, -6.16302884182399868e-01_KREAL, &
                -2.67398372167765297e-01_KREAL, -5.29786439318511299e-02_KREAL, &
                5.29786439318511299e-02_KREAL, 2.67398372167765297e-01_KREAL, &
                6.16302884182399868e-01_KREAL, 1.06424631211622400e+00_KREAL, &
                1.58885586227005504e+00_KREAL, 2.18392115309585799e+00_KREAL, &
                2.86313388370807509e+00_KREAL, 3.68600716272439710e+00_KREAL]
            raw_weight=[ &
                1.19259692659534398e-06_KREAL, 2.02063649132410713e-04_KREAL, &
                5.36793575602533338e-03_KREAL, 4.48141099174629021e-02_KREAL, &
                1.57448282618790292e-01_KREAL, 2.75953397988421778e-01_KREAL, &
                2.68330754472638822e-01_KREAL, 1.34109188453359507e-01_KREAL, &
                1.34109188453359507e-01_KREAL, 2.68330754472638822e-01_KREAL, &
                2.75953397988421778e-01_KREAL, 1.57448282618790292e-01_KREAL, &
                4.48141099174629021e-02_KREAL, 5.36793575602533338e-03_KREAL, &
                2.02063649132410713e-04_KREAL, 1.19259692659534398e-06_KREAL]
        case('half32')
            raw_node=[ &
                -5.67501793404192334e+00_KREAL, -4.94637720404838621e+00_KREAL, &
                -4.34360634547017010e+00_KREAL, -3.80737711675589763e+00_KREAL, &
                -3.31576927503869756e+00_KREAL, -2.85826652854326690e+00_KREAL, &
                -2.42945049160214266e+00_KREAL, -2.02680815216886678e+00_KREAL, &
                -1.64985424039743411e+00_KREAL, -1.29978932127703506e+00_KREAL, &
                -9.79404170330729240e-01_KREAL, -6.93073720301999319e-01_KREAL, &
                -4.46696225961683790e-01_KREAL, -2.47397669452455921e-01_KREAL, &
                -1.02802245237917422e-01_KREAL, -1.97536584600771156e-02_KREAL, &
                1.97536584600771156e-02_KREAL, 1.02802245237917422e-01_KREAL, &
                2.47397669452455921e-01_KREAL, 4.46696225961683790e-01_KREAL, &
                6.93073720301999319e-01_KREAL, 9.79404170330729240e-01_KREAL, &
                1.29978932127703506e+00_KREAL, 1.64985424039743411e+00_KREAL, &
                2.02680815216886678e+00_KREAL, 2.42945049160214266e+00_KREAL, &
                2.85826652854326690e+00_KREAL, 3.31576927503869756e+00_KREAL, &
                3.80737711675589763e+00_KREAL, 4.34360634547017010e+00_KREAL, &
                4.94637720404838621e+00_KREAL, 5.67501793404192334e+00_KREAL]
            raw_weight=[ &
                8.67420445249457028e-15_KREAL, 1.53767791618982658e-11_KREAL, &
                3.61154913974275996e-09_KREAL, 2.59000761941504528e-07_KREAL, &
                7.95049571962246903e-06_KREAL, 1.25286221329562052e-04_KREAL, &
                1.13569531068878323e-03_KREAL, 6.40991442405013188e-03_KREAL, &
                2.39197096186835705e-02_KREAL, 6.19996960991567869e-02_KREAL, &
                1.16572490553503091e-01_KREAL, 1.65438637755609519e-01_KREAL, &
                1.83562801116245566e-01_KREAL, 1.62921292314545607e-01_KREAL, &
                1.13608556894151586e-01_KREAL, 5.05246320213775663e-02_KREAL, &
                5.05246320213775663e-02_KREAL, 1.13608556894151586e-01_KREAL, &
                1.62921292314545607e-01_KREAL, 1.83562801116245566e-01_KREAL, &
                1.65438637755609519e-01_KREAL, 1.16572490553503091e-01_KREAL, &
                6.19996960991567869e-02_KREAL, 2.39197096186835705e-02_KREAL, &
                6.40991442405013188e-03_KREAL, 1.13569531068878323e-03_KREAL, &
                1.25286221329562052e-04_KREAL, 7.95049571962246903e-06_KREAL, &
                2.59000761941504528e-07_KREAL, 3.61154913974275996e-09_KREAL, &
                1.53767791618982658e-11_KREAL, 8.67420445249457028e-15_KREAL]
        end select
        ue=sqrt(2.0_KREAL)*raw_node;we=sqrt(2.0_KREAL)*exp(raw_node*raw_node)*raw_weight
        ! Optional dimensional change of quadrature variable; defaults preserve old arithmetic.
        velocity_scale=1.0_KREAL
        st=1
        if(S2VelocityScale>0.0_KREAL)then
            velocity_scale=S2VelocityScale
        else
            call get_environment_variable('S2_VELOCITY_SCALE',scale_env,status=st)
        endif
        if(st==0)then
            read(scale_env,*,iostat=ios)velocity_scale
            if(ios/=0)error stop 'invalid velocity scale syntax'
            if(.not.ieee_is_finite(velocity_scale).or.velocity_scale<=0.0_KREAL)error stop 'invalid velocity scale'
        elseif(st/=1)then
            error stop 'truncated velocity scale'
        endif
        if(velocity_scale/=1.0_KREAL)then
            ue=velocity_scale*ue;we=velocity_scale*we
        endif
        do j=1,n;do i=1,n;k=i+(j-1)*n
            u(k)=ue(i);v(k)=ue(j);weight(k)=we(i)*we(j)
        enddo;enddo
    end subroutine
    subroutine FillFull16(u,v,weight,raw_node,raw_weight)
        real(KREAL),intent(out)::u(:),v(:),weight(:),raw_node(:),raw_weight(:)
        real(KREAL),parameter::x(16)=[-4.68873893930581875_KREAL,-3.86944790486012291_KREAL, &
            -3.17699916197995602_KREAL,-2.54620215784748138_KREAL,-1.95178799091625388_KREAL, &
            -1.38025853919888086_KREAL,-0.82295144914465590_KREAL,-0.27348104613815250_KREAL, &
             0.27348104613815250_KREAL,0.82295144914465590_KREAL,1.38025853919888086_KREAL, &
             1.95178799091625388_KREAL,2.54620215784748138_KREAL,3.17699916197995602_KREAL, &
             3.86944790486012291_KREAL,4.68873893930581875_KREAL]
        real(KREAL),parameter::wr(16)=[2.65480747401116728e-10_KREAL,2.32098084486520325e-7_KREAL, &
            2.71186009253787598e-5_KREAL,9.32284008624180709e-4_KREAL,1.28803115355099890e-2_KREAL, &
            8.38100413989858323e-2_KREAL,2.80647458528533678e-1_KREAL,5.07929479016613561e-1_KREAL, &
            5.07929479016613561e-1_KREAL,2.80647458528533678e-1_KREAL,8.38100413989858323e-2_KREAL, &
            1.28803115355099890e-2_KREAL,9.32284008624180709e-4_KREAL,2.71186009253787598e-5_KREAL, &
            2.32098084486520325e-7_KREAL,2.65480747401116728e-10_KREAL]
        real(KREAL)::ue(16),we(16)
        integer::i,j,k
        raw_node=x;raw_weight=wr;ue=sqrt(2.0_KREAL)*x;we=sqrt(2.0_KREAL)*exp(x*x)*wr
        do j=1,16;do i=1,16;k=i+(j-1)*16
            u(k)=ue(i);v(k)=ue(j);weight(k)=we(i)*we(j)
        enddo;enddo
    end subroutine FillFull16
end module mod_VelocityTables

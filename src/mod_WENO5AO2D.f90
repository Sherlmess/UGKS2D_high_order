module WENO5AO2D
    use Precision, only: KREAL
    use WENO5AOScalar, only: rk => rk, pidx, WENO5AOCellValueCache, &
        WENO5AOCellFirstDirectionCache, prepare_weno5ao_cell_value, &
        evaluate_weno5ao_cell_value, prepare_weno5ao_cell_first_direction, &
        evaluate_weno5ao_cell_first_direction
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(KREAL), parameter, public :: GQ3_NODE(3) = [-sqrt(3.0_KREAL/5.0_KREAL), &
        0.0_KREAL,sqrt(3.0_KREAL/5.0_KREAL)]
    real(KREAL), parameter, public :: GQ3_WEIGHT(3) = [5.0_KREAL/9.0_KREAL, &
        8.0_KREAL/9.0_KREAL,5.0_KREAL/9.0_KREAL]

contains

    subroutine reconstruct_weno5ao_tensor_point(field,i,j,xi,eta,value,ok)
        real(KREAL), intent(in) :: field(:,:),xi,eta
        integer, intent(in) :: i,j
        real(KREAL), intent(out) :: value
        logical, intent(out) :: ok
        integer :: nx,ny,s,col
        real(KREAL), allocatable :: line(:),outer(:)
        type(WENO5AOCellValueCache) :: cache

        value=0.0_KREAL
        ok=.false.
        nx=size(field,1); ny=size(field,2)
        if (nx<5 .or. ny<5 .or. i<1 .or. i>nx .or. j<1 .or. j>ny) return
        if (.not.all(ieee_is_finite(field))) return
        allocate(line(ny),outer(5))
        do s=1,5
            col=pidx(i+s-3,nx)
            line=field(col,:)
            call prepare_weno5ao_cell_value(line,j,cache)
            if (.not.cache%valid) then
                deallocate(line,outer)
                return
            end if
            call evaluate_weno5ao_cell_value(cache,eta,outer(s))
        end do
        call prepare_weno5ao_cell_value(outer,3,cache)
        if (.not.cache%valid) then
            deallocate(line,outer)
            return
        end if
        call evaluate_weno5ao_cell_value(cache,xi,value)
        ok=ieee_is_finite(value)
        deallocate(line,outer)
    end subroutine reconstruct_weno5ao_tensor_point

    subroutine reconstruct_weno5ao_tensor_point_direction(field,field_direction,i,j,xi,eta, &
        value,value_direction,ok)
        real(KREAL), intent(in) :: field(:,:),field_direction(:,:),xi,eta
        integer, intent(in) :: i,j
        real(KREAL), intent(out) :: value,value_direction
        logical, intent(out) :: ok
        integer :: nx,ny,s,col
        real(KREAL), allocatable :: line(:),line_direction(:),outer(:),outer_direction(:)
        type(WENO5AOCellFirstDirectionCache) :: cache
        logical :: differentiable

        value=0.0_KREAL; value_direction=0.0_KREAL; ok=.false.
        nx=size(field,1); ny=size(field,2)
        if (nx<5 .or. ny<5 .or. size(field_direction,1)/=nx .or. &
            size(field_direction,2)/=ny .or. i<1 .or. i>nx .or. j<1 .or. j>ny) return
        if (.not.all(ieee_is_finite(field)) .or. .not.all(ieee_is_finite(field_direction))) return
        allocate(line(ny),line_direction(ny),outer(5),outer_direction(5))
        do s=1,5
            col=pidx(i+s-3,nx)
            line=field(col,:)
            line_direction=field_direction(col,:)
            call prepare_weno5ao_cell_first_direction(line,line_direction,j,cache,differentiable)
            if (.not.differentiable .or. .not.cache%valid) then
                deallocate(line,line_direction,outer,outer_direction)
                return
            end if
            call evaluate_weno5ao_cell_first_direction(cache,eta,outer(s),outer_direction(s))
        end do
        call prepare_weno5ao_cell_first_direction(outer,outer_direction,3,cache,differentiable)
        if (.not.differentiable .or. .not.cache%valid) then
            deallocate(line,line_direction,outer,outer_direction)
            return
        end if
        call evaluate_weno5ao_cell_first_direction(cache,xi,value,value_direction)
        ok=ieee_is_finite(value) .and. ieee_is_finite(value_direction)
        deallocate(line,line_direction,outer,outer_direction)
    end subroutine reconstruct_weno5ao_tensor_point_direction

    subroutine reconstruct_weno5ao_tensor_face(field,axis,face,tangent,eta,left_value,right_value,ok)
        real(KREAL), intent(in) :: field(:,:),eta
        integer, intent(in) :: axis,face,tangent
        real(KREAL), intent(out) :: left_value,right_value
        logical, intent(out) :: ok
        integer :: nx,ny,left_cell,right_cell

        left_value=0.0_KREAL; right_value=0.0_KREAL; ok=.false.
        nx=size(field,1); ny=size(field,2)
        if (axis<1 .or. axis>2) return
        if (axis==1) then
            if (nx<5 .or. ny<5 .or. face<1 .or. face>nx .or. tangent<1 .or. tangent>ny) return
            left_cell=face
            right_cell=pidx(face+1,nx)
        else
            if (nx<5 .or. ny<5 .or. face<1 .or. face>ny .or. tangent<1 .or. tangent>nx) return
            left_cell=face
            right_cell=pidx(face+1,ny)
        end if
        call tensor_side_value(field,axis,left_cell,tangent,0.5_KREAL,eta,left_value,ok)
        if (.not.ok) return
        call tensor_side_value(field,axis,right_cell,tangent,-0.5_KREAL,eta,right_value,ok)
    end subroutine reconstruct_weno5ao_tensor_face

    subroutine reconstruct_weno5ao_tensor_face_direction(field,field_direction,axis,face,tangent,eta, &
        left_value,right_value,left_direction,right_direction,ok)
        real(KREAL), intent(in) :: field(:,:),field_direction(:,:),eta
        integer, intent(in) :: axis,face,tangent
        real(KREAL), intent(out) :: left_value,right_value,left_direction,right_direction
        logical, intent(out) :: ok
        integer :: nx,ny,left_cell,right_cell

        left_value=0.0_KREAL; right_value=0.0_KREAL
        left_direction=0.0_KREAL; right_direction=0.0_KREAL; ok=.false.
        nx=size(field,1); ny=size(field,2)
        if (axis<1 .or. axis>2) return
        if (size(field_direction,1)/=nx .or. size(field_direction,2)/=ny) return
        if (axis==1) then
            if (nx<5 .or. ny<5 .or. face<1 .or. face>nx .or. tangent<1 .or. tangent>ny) return
            left_cell=face
            right_cell=pidx(face+1,nx)
        else
            if (nx<5 .or. ny<5 .or. face<1 .or. face>ny .or. tangent<1 .or. tangent>nx) return
            left_cell=face
            right_cell=pidx(face+1,ny)
        end if
        call tensor_side_direction(field,field_direction,axis,left_cell,tangent,0.5_KREAL,eta, &
            left_value,left_direction,ok)
        if (.not.ok) return
        call tensor_side_direction(field,field_direction,axis,right_cell,tangent,-0.5_KREAL,eta, &
            right_value,right_direction,ok)
    end subroutine reconstruct_weno5ao_tensor_face_direction

    subroutine tensor_side_value(field,axis,outer_center,tangent,xi,eta,value,ok)
        real(KREAL), intent(in) :: field(:,:),xi,eta
        integer, intent(in) :: axis,outer_center,tangent
        real(KREAL), intent(out) :: value
        logical, intent(out) :: ok
        integer :: nx,ny,s,index,inner_size,inner_center
        real(KREAL), allocatable :: line(:),outer(:)
        type(WENO5AOCellValueCache) :: cache

        value=0.0_KREAL; ok=.false.
        nx=size(field,1); ny=size(field,2)
        if (axis==1) then
            inner_size=ny; inner_center=tangent
        else
            inner_size=nx; inner_center=tangent
        end if
        allocate(line(inner_size),outer(5))
        do s=1,5
            index=pidx(outer_center+s-3,merge(nx,ny,axis==1))
            if (axis==1) then
                line=field(index,:)
            else
                line=field(:,index)
            end if
            call prepare_weno5ao_cell_value(line,inner_center,cache)
            if (.not.cache%valid) then
                deallocate(line,outer)
                return
            end if
            call evaluate_weno5ao_cell_value(cache,eta,outer(s))
        end do
        call prepare_weno5ao_cell_value(outer,3,cache)
        if (.not.cache%valid) then
            deallocate(line,outer)
            return
        end if
        call evaluate_weno5ao_cell_value(cache,xi,value)
        ok=ieee_is_finite(value)
        deallocate(line,outer)
    end subroutine tensor_side_value

    subroutine tensor_side_direction(field,field_direction,axis,outer_center,tangent,xi,eta, &
        value,value_direction,ok)
        real(KREAL), intent(in) :: field(:,:),field_direction(:,:),xi,eta
        integer, intent(in) :: axis,outer_center,tangent
        real(KREAL), intent(out) :: value,value_direction
        logical, intent(out) :: ok
        integer :: nx,ny,s,index,inner_size,inner_center
        real(KREAL), allocatable :: line(:),line_direction(:),outer(:),outer_direction(:)
        type(WENO5AOCellFirstDirectionCache) :: cache
        logical :: differentiable

        value=0.0_KREAL; value_direction=0.0_KREAL; ok=.false.
        nx=size(field,1); ny=size(field,2)
        if (size(field_direction,1)/=nx .or. size(field_direction,2)/=ny) return
        if (axis==1) then
            inner_size=ny; inner_center=tangent
        else
            inner_size=nx; inner_center=tangent
        end if
        allocate(line(inner_size),line_direction(inner_size),outer(5),outer_direction(5))
        do s=1,5
            index=pidx(outer_center+s-3,merge(nx,ny,axis==1))
            if (axis==1) then
                line=field(index,:); line_direction=field_direction(index,:)
            else
                line=field(:,index); line_direction=field_direction(:,index)
            end if
            call prepare_weno5ao_cell_first_direction(line,line_direction,inner_center,cache,differentiable)
            if (.not.differentiable .or. .not.cache%valid) then
                deallocate(line,line_direction,outer,outer_direction)
                return
            end if
            call evaluate_weno5ao_cell_first_direction(cache,eta,outer(s),outer_direction(s))
        end do
        call prepare_weno5ao_cell_first_direction(outer,outer_direction,3,cache,differentiable)
        if (.not.differentiable .or. .not.cache%valid) then
            deallocate(line,line_direction,outer,outer_direction)
            return
        end if
        call evaluate_weno5ao_cell_first_direction(cache,xi,value,value_direction)
        ok=ieee_is_finite(value) .and. ieee_is_finite(value_direction)
        deallocate(line,line_direction,outer,outer_direction)
    end subroutine tensor_side_direction

end module WENO5AO2D



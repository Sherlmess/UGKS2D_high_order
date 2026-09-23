!#UGKS2D parallel simulation code: developed by JHOU 
!---------------------------------
!>To be added
!---------------------------------
!1) Axisymmetrical flow solver
!2) General Maxwell reflection boundary
!
!---------------------------------
!>Notice
!---------------------------------
!1)The symmetry conditions only apply to x- or y- axis
!2)The Maxwell reflection only apply to x- or y-parallel boundary
!
!------------------------------------------------------------------------------------------------
!added by SjZhang   
!2025/12/7 : add unstructured velocity space support(only for upper/lower symmetry boundary condition)
!2026/5/25 : simplified the code structure

!###==============================================================================================


program UGKS2D
    use ImportPara
    use mpi
    use mpi_wrapper
    use InitializeFlow
    use Solver
    use Output
    use mod_S2D1, only: IsS2D1,RunS2D1
    implicit none
    
    !Initial of MPI
    call MPI_Init(ierr) 
    
    !Get myid
    call MPI_Comm_Rank(MPI_COMM_WORLD,myID,ierr)
    
    if(myID==0) then
        write(*,'(a)')  '-------------UGKS2D-DIS simulation beginning: Parallel in XY, version 2.0 (by JHOU, oujihui@tju.edu.cn)-------------'
    endif
    
    !Import parameters
    call InitPara()

    if(IsS2D1()) then
        call RunS2D1()
        call MPI_Finalize(ierr)
        stop
    endif
    
    !Partion mpi
    call Partion_mpi()
    
    !Initialization
    call InitFlow()

    !Time integration
    call Solve_UGKS()

    !Output solution
    call Output_result()
    
    !Finalize MPI
    call MPI_Finalize(ierr)

end program UGKS2D

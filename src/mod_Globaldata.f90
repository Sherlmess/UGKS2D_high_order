    !--------------------------------------------------
    !>Specify the precision
    !--------------------------------------------------
    module Precision
        use mpi, only: MPI_DOUBLE_PRECISION
        implicit none
        ! include "mpif.h" 
        
        integer, parameter                  :: KINT     = kind(1)
        integer, parameter                  :: KREAL    = kind(0.d0)
        integer, parameter                  :: MPI_KREAL= MPI_DOUBLE_PRECISION
    end module Precision

    !--------------------------------------------------
    !>Define the global constant parameters
    !--------------------------------------------------
    module const_para
        use Precision
        implicit none
    
        integer(KINT),parameter             :: GHOST = 2                !Ghost cell number: applied to physical boundary
    
        real(KREAL), parameter              :: PI   = 4.0*atan(1.0)     !Pi
        real(KREAL), parameter              :: bz   = 1.3806488E-23     !Boltzmann constant
        real(KREAL), parameter              :: SMV  = tiny(0.0)         !Small value to avoid 0/0
        real(KREAL), parameter              :: UP   = 1.0               !Used in sign() function
        
        !Initialize flow
        integer(KINT), parameter            :: Init_Freestream = 0      !initialize flow using initial or freestream parameter
        integer(KINT), parameter            :: Init_Macroflow = 1       !initialize flow using macroflow file
        integer(KINT), parameter            :: Init_Restart = 2         !initialize flow using restart file
    
        !Time step
        integer(KINT), parameter            :: Time_global_fixed=0      !fixed global time step
        integer(KINT), parameter            :: Time_global_varying=1    !varying global time step
        integer(KINT), parameter            :: Time_local_varying=2     !varying local time step
    
        !Reconstruction
        integer(KINT), parameter            :: Reconst_Vanleer = 1      !VanleerLimiter reconstruction

        !Quadrature method: supported options
        integer(KINT), parameter            :: Quad_NewtonCotes4 = 1    !Newton-Cotes 4th order
        integer(KINT), parameter            :: Quad_GaussHermite = 2    !Gauss-Hermite
        integer(KINT), parameter            :: Quad_Unstruct = 3        !Unstructured velocity grid (Read from file)

        !Evolution: supported option
        integer(KINT), parameter            :: Flux_Multidimension = 1  !Flux with multi-dimension
        
        !Solid wall boundary type
        integer(KINT), parameter            :: FluxWall_Kinetic = 0     !Kinetic diffuse boundary condition
        integer(KINT), parameter            :: FluxWall_Multiscale = 1  !Multiscale boundary condition
        integer(KINT), parameter            :: FluxWall_Maxwell = 2     !Maxwell boundary condition
    
        !Direction
        integer(KINT), parameter            :: Idirc = 1                !I direction
        integer(KINT), parameter            :: Jdirc = 2                !J direction
        integer(KINT), parameter            :: Lside = -1               !Left (or lower) side (donot change the value)
        integer(KINT), parameter            :: Rside = 1                !Right (or upper) side (donot change the value)
        
        !Gas viscous law
        integer(KINT), parameter            :: Sutherland = 0           !Sutherland law
        integer(KINT), parameter            :: PowerLaw = 1             !HS, VHS or VSS model: viscosity-temperature power law
    
        !Boundary condition
        !Special geometries
        integer(KINT), parameter            :: BGconst_General=0        !general geometry
        integer(KINT), parameter            :: BGconst_FlatPlate=1      !flat-plate geometry
    
        !General boundary type
        integer(KINT), parameter            :: BCconst_wall_isothermal=2!isothermal wall (velocity slip and temperature jump)
        integer(KINT), parameter            :: BCconst_wall_insulate=21 !insulate wall (no-slip and adiabatic condition)
        integer(KINT), parameter            :: BCconst_symmetry=3       !symmetry boundary
        integer(KINT), parameter            :: BCconst_freestream=4     !free-stream-in boundary
        integer(KINT), parameter            :: BCconst_inflow=5         !inflow boundary
        integer(KINT), parameter            :: BCconst_extrapolation=6  !extrapolation boundary
        integer(KINT), parameter            :: BCconst_periodic=10      !periodic boundary
        
        !DIS boundary type
        integer(KINT), parameter            :: DISBC_None=0             !no B.C. disturbance
        integer(KINT), parameter            :: DISBC_Inflow_TS=1        !inflow TS wave
        integer(KINT), parameter            :: DISBC_Free_wave=2        !free-stream wave
    
    end module const_para
    
    
    !--------------------------------------------------
    !>Define global variables for MPI: shared by all nodes
    !--------------------------------------------------
    module mpi_para
        use Precision
        use mpi, only: MPI_STATUS_SIZE
        implicit none
        
        integer(KINT),parameter             :: LAP = 2                  !Overlap length for Block-Block (LAP >= GHOST)
        
        !num process
        integer(KINT)                       :: myID                     !Current process
        integer(KINT)                       :: IP_XY,NP_XY              !NP_ALL: total process, NP_XY: XY-process, NP_UV: UV-process
        
        !Commucation node
        integer(KINT)                       :: ID_X_minus,ID_X_plus,ID_Y_minus,ID_Y_plus
        
        !block type
        type block_type
            integer(KINT)                   :: NI,NJ                    !Mesh points in I- and J-dirc
            integer(KINT)                   :: Origin_I,Origin_J        !Global i j index of left-bottom corner
        end type block_type
        
        !physical (xy) space
        integer(KINT)                       :: NP_X,NP_Y                !Total block number in x- and y-dirc
        integer(KINT)                       :: IP_X,IP_Y                !Current IP index in XY-dirc
        type(block_type),allocatable,dimension(:) :: Block_XY
        
        !For MPI run
        integer(KINT)                       :: Status(MPI_status_Size),ierr 
    
    end module mpi_para
    
    !--------------------------------------------------
    !>Define control parameters for the simulation
    !--------------------------------------------------
    module control_para
        use const_para
        use mpi_para
        implicit none
        !scheme
        integer(KINT)                       :: Quadrature_method        !Quadrature_method
        integer(KINT)                       :: FluxWall_type            !Solid wall flux
        integer(KINT)                       :: TimeStep_method          !Time integration method
        integer(KINT)                       :: Iflag_InitFlow           !Initialize flow option
        integer(KINT)                       :: Iflag_WriteSurf          !Whether write surface quantities: 0 for no, 1 for yes
        integer(KINT)                       :: Iflag_WriteRestart       !Whether write restart file: 0 for no, 1 for yes
        character(len=24)                   :: TimeIntegration_method   !Explicit, ImplicitSteady, or ImplicitUnsteady
        character(len=24)                   :: PeriodicComparisonProfile !Optional S2D-3 BGK constant-tau comparison
        real(KREAL)                         :: ComparisonFaceTauMin,ComparisonFaceTauMax
        real(KREAL)                         :: ComparisonCellTauMin,ComparisonCellTauMax
        integer(KINT)                       :: ComparisonFaceTauCalls,ComparisonCellTauCalls
        integer(KINT)                       :: S2D1_initial_mode        !0 sin(x+y), 1 sin(x), 2 sin(y), 3 constant
        real(KREAL)                         :: Implicit_nCFL
        real(KREAL)                         :: Implicit_ResidualTol
        real(KREAL)                         :: Implicit_MomentTol
        real(KREAL)                         :: Implicit_MinAlpha
        integer(KINT)                       :: Implicit_MacroSweeps
        integer(KINT)                       :: Implicit_MicroSweeps
        integer(KINT)                       :: Implicit_MaxRetries
        integer(KINT)                       :: Implicit_MaxIter
        real(KREAL)                         :: Unsteady_Theta
        real(KREAL)                         :: Unsteady_InnerTol
        integer(KINT)                       :: Unsteady_MaxInner
        
        !Computation
        real(KREAL)                         :: dt_global                !Global time step
        real(KREAL)                         :: CFL                      !CFL number
        real(KREAL)                         :: EPS                      !Convergence criteria
        real(KREAL)                         :: simtime,dsimtime         !Simulation time (CPU time)
        real(KREAL)                         :: dcommtime                !Commucation time
        integer(KINT)                       :: iter                     !Current iteration
        integer(KINT)                       :: max_iter                 !Maximal iterations
        integer(KINT)                       :: Nstep_show               !Show residuals every Nstep_show step
        integer(KINT)                       :: Nstep_write              !Write files every Nstep_write step
        real(KREAL)                         :: res(4)                   !Global L2-normalized residual
        real(KREAL)                         :: res_relative(4)          !Relative residual before zero-scale selection
        real(KREAL)                         :: res_absolute(4)          !Area-weighted absolute residual rate
        real(KREAL)                         :: res_num(4)               !Global L2 numerator of residuals
        real(KREAL)                         :: res_den(4)               !Global normalization denominator
        real(KREAL)                         :: res_rms(4)               !Global unnormalized RMS residuals
        real(KREAL)                         :: res_max(4)               !Global max of per-cell residuals
        
        !For output
        character(len=8)                    :: date
        character(len=10)                   :: time
        character(len=100)                  :: fileName

    contains

        pure logical function IsImplicitSteady()
            IsImplicitSteady = trim(adjustl(TimeIntegration_method)) == 'ImplicitSteady'
        end function IsImplicitSteady

        pure logical function IsS2Steady()
            IsS2Steady = trim(adjustl(TimeIntegration_method)) == 'S2Steady'
        end function IsS2Steady

        pure logical function IsImplicitUnsteady()
            IsImplicitUnsteady = trim(adjustl(TimeIntegration_method)) == 'ImplicitUnsteady'
        end function IsImplicitUnsteady

    end module control_para
    
    
    !--------------------------------------------------
    !>Define flow parameters
    !--------------------------------------------------
    module flow_para
        use control_para
        implicit none
        
        !Nondimensional free-stream parameters
        real(KREAL)                         :: Ma                       !Mach number
        real(KREAL)                         :: Re                       !Reynolds number
        real(KREAL)                         :: AOA                      !Angle of attack
        real(KREAL)                         :: Tinf                     !Free-stream temperature
        real(KREAL)                         :: Tw                       !Wall temperature
        real(KREAL)                         :: Kn                       !Knudsen number (mean free path is defined by mu/rho*sqrt(pi/(2RT)))
    
        !Gas properties
        integer(KINT)                       :: ck                       !Internal degree of freedom, here 1 denotes 2D flow for monatomic gas
        real(KREAL)                         :: Pr                       !Prandtl number
        integer(KINT)                       :: flag_vislaw              !Viscous law
        real(KREAL)                         :: omega                    !Temperature dependence index in HS/VHS/VSS model
        real(KREAL)                         :: Tcm                      !Temperature constant in Sutherland law
        real(KREAL)                         :: gamma                    !Ratio of specific heat
        real(KREAL)                         :: mu_ref                   !Viscosity coefficient in reference state      
        
        !Physical cell center and velocity space
        integer(KINT)                       :: Nx_all,Ny_all            !Total cell number in x- and y-direction
        integer(KINT)                       :: Nx,Ny                    !Local-block cell number in x- and y-direction
        
        integer(KINT)                       :: Nu,Nv                    !Number of mesh points in velocity space for u and v
        real(KREAL)                         :: U_MIN,U_MAX,V_MIN,V_MAX  !Minimum and maximum micro velocity
    
        !Free-stream or wall state: density, u, v, lambda=1/T  
        real(KREAL)                         :: Prim_inf(4)              !Free-stream gas
        real(KREAL)                         :: Prim_wall(4)             !Wall condition

        !Boundary conditio
        integer(KINT)                       :: BC_Type(4)               !Boundary conditions         
        integer(KINT)                       :: BG_geometry              !0: general, 1:flat plate
        real(KREAL)                         :: ACC                      !Accommodation coefficient at the solid wall 
    
    end module flow_para
    
    
    !--------------------------------------------------
    !>Define flow data
    !--------------------------------------------------
    module flow_data
        use flow_para
        implicit none
        
        !--------------------------------------------------
        !Define mesh and flow variables
        !--------------------------------------------------
        !>Cell and interface index
        !---------------------------------
        !            (i,j+1)              |
        !       ----------------          |
        !       |              |          |
        !       |              |          |
        !       |              |          |
        !(i-1,j)|     (i,j)    |(i+1,j)   |
        !       |      Cell    |interface |
        !       |              |          |
        !       |              |          |
        !       ----------------          |
        !           (i,j-1)               |
        !---------------------------------
        !--------------------------------------------------
        !Define mesh data structure
        !--------------------------------------------------
        !Cell center
        type CellCenter
            !Time step
            real(KREAL)                                     :: dt_local                 !Local time step
            !Geometry
            real(KREAL)                                     :: x,y                      !Cell center coordinates
            real(KREAL)                                     :: area                     !Cell area
            real(KREAL)                                     :: length(2)                !Cell length in i and j direction
            !Flow field: conVars=[density, x-momentum, y-momentum, total energy]
            real(KREAL)                                     :: conVars(4)               !Conservative variables at cell center
            real(KREAL)                                     :: res_cell(4)              !Per-cell residual of conservative variables
            real(KREAL)                                     :: sw(4,2)                  !Slope of conVars in i and j direction
            real(KREAL), pointer, dimension(:)              :: h_vec=>null(),b_vec=>null()      !1D velocity-space view of h/b, k=i+(j-1)*Nu
            real(KREAL), pointer, dimension(:,:)            :: sh_vec=>null(),sb_vec=>null()    !1D velocity-space slope view, second index is i/j direction
            !##Initial flow data (for disturbance field computation)
            real(KREAL)                                     :: conVars0(4)
            real(KREAL)                                     :: res_conVars0(4)
            real(KREAL), pointer, dimension(:)              :: h0_vec=>null(),b0_vec=>null()    !1D view of initial h/b
            real(KREAL), pointer, dimension(:)              :: res_h0_vec=>null(),res_b0_vec=>null()    !1D view of initial residual h/b
        end type CellCenter

        !Cell interface
        type CellInterface
            !Geometry
            real(KREAL)                                     :: length                   !Length of cell interface
            real(KREAL)                                     :: cosx,cosy                !Directional cosine in global frame
            !Flow variables
            real(KREAL)                                     :: conVars(4),swt(4)        !Conservative variables and tangential derivative at cell interface
            real(KREAL)                                     :: FlowVars(5),Aerodynamics(3)  !Store surface flow variables and aerodynamics
            real(KREAL)                                     :: FlowVars0(5),Aerodynamics0(3)    !Initial state (for DIS simulation)
            !Flux
            real(KREAL)                                     :: L_flux(4)                !Conservative variables flux at cell interface 
            real(KREAL)                                     :: R_flux(4)                !Conservative variables flux at cell interface
            real(KREAL), pointer, dimension(:)              :: L_flux_h_vec=>null(),L_flux_b_vec=>null()    !1D left micro-flux view
            real(KREAL), pointer, dimension(:)              :: R_flux_h_vec=>null(),R_flux_b_vec=>null()    !1D right micro-flux view
        end type CellInterface
    
        !Mesh information
        real(KREAL),allocatable,dimension(:,:)              :: gridx_all,gridy_all      !cell vertice: total mesh
        real(KREAL),allocatable,dimension(:,:)              :: gridxc_all,gridyc_all    !cell center: total mesh
        real(KREAL),allocatable,dimension(:,:)              :: gridx,gridy              !grid vertice: current mesh
        
        !cell center and interface information 
        type(CellCenter),allocatable,dimension(:,:)         :: ctr                      !Cell center
        type(CellInterface),allocatable,dimension(:,:)      :: vface,hface              !Vertical and horizontal interfaces
        real(KREAL),allocatable,dimension(:,:,:)            :: res_macro_flux_store     !Macro face-flux contribution dU_flux, stored as (4,nx,ny)
        real(KREAL),allocatable,dimension(:,:,:)            :: res_macro_store          !Macro explicit update RHS after base-flow/fringe terms, stored as (4,nx,ny)
        real(KREAL),allocatable,dimension(:,:,:)            :: res_h_flux_store         !Micro h face-flux contribution dh_flux, stored as (Nvel,nx,ny)
        real(KREAL),allocatable,dimension(:,:,:)            :: res_b_flux_store         !Micro b face-flux contribution db_flux, stored as (Nvel,nx,ny)
        real(KREAL),allocatable,target,dimension(:,:,:)     :: cell_h_store,cell_b_store              !Pure 1D backing arrays for ctr%h_vec/b_vec
        real(KREAL),allocatable,target,dimension(:,:,:)     :: cell_h0_store,cell_b0_store            !Pure 1D backing arrays for initial h/b
        real(KREAL),allocatable,target,dimension(:,:,:)     :: cell_res_h0_store,cell_res_b0_store    !Pure 1D backing arrays for initial VDF residuals
        real(KREAL),allocatable,target,dimension(:,:,:,:)   :: cell_sh_store,cell_sb_store            !Pure 1D backing arrays for ctr%sh_vec/sb_vec slopes
        real(KREAL),allocatable,target,dimension(:,:,:)     :: v_l_flux_h_store,v_l_flux_b_store      !Pure 1D backing arrays for vertical-face left fluxes
        real(KREAL),allocatable,target,dimension(:,:,:)     :: v_r_flux_h_store,v_r_flux_b_store      !Pure 1D backing arrays for vertical-face right fluxes
        real(KREAL),allocatable,target,dimension(:,:,:)     :: h_l_flux_h_store,h_l_flux_b_store      !Pure 1D backing arrays for horizontal-face left fluxes
        real(KREAL),allocatable,target,dimension(:,:,:)     :: h_r_flux_h_store,h_r_flux_b_store      !Pure 1D backing arrays for horizontal-face right fluxes
        
        !Velocity space
        real(KREAL),allocatable,dimension(:)                :: weight                   !Pure 1D quadrature weights, k=i+(j-1)*Nu
        real(KREAL),allocatable,dimension(:)                :: uSpace,vSpace            !Pure 1D velocity coordinates, k=i+(j-1)*Nu
        integer(KINT)                                       :: Nvel                     !Flattened velocity-space point count, Nvel=Nu*Nv
        real(KREAL),allocatable,dimension(:)                :: wVel,uVel,vVel           !Flattened weight/u/v arrays using k=i+(j-1)*Nu
        ! symmetric boundary mapping table: SymMap_V(k) stores the index of the point with opposite v velocity to point k
        integer(KINT), allocatable :: SymMap_V(:)

        !Inflow boundary
        real(KREAL),allocatable,dimension(:,:)              :: inflow_conVars
        real(KREAL),allocatable,dimension(:,:)              :: inflow_h,inflow_b       !Pure 1D inflow VDF cache, indexed as (j,k)
    
    end module flow_data
    
    
    !--------------------------------------------------
    !>Define flow data and control paras for disturbance
    !--------------------------------------------------
    module Flow_DIS
        use Precision
        implicit none
        
        !Control para
        integer(KINT)                                       ::  Iflag_DIS               !0:steady base flow, 1:unsteady disturbance evolution (subtract initial Res.)
        integer(KINT)                                       ::  Iflag_Fringe            !0:not apply fringe, 1: apply fringe in the outflow
        integer(KINT)                                       ::  DISBC_Type              !0:no B.C. disturbance, 1:inflow TS wave, 2: free-stream wave
        
        !Global time for disturbance evolution
        real(KREAL)                                         ::  time_dis                !global time
        
        !Define disturbance Type
        type wave
            character(len=2)                                ::  name                    !TS (TS wave), F (fast wave), S (slow wave), E (entropy wave), V (vortex wave)
            real(KREAL)                                     ::  omer,omei               !frequency
            real(KREAL)                                     ::  amp                     !amplitude
            real(KREAL)                                     ::  alfr,alfi               !stream-wise wave number (used in TS wave)
            real(KREAL)                                     ::  theta                   !incident angle (used in free-stream wave)
        end type wave
        
        !Define inflow TS wave
        integer(KINT)                                       ::  Num_dis_TS              !number of TS wave
        type(wave),dimension(20)                            ::  wave_TS                 !TS wave (maximum 20, could be changed)
        real(KREAL),allocatable,dimension(:,:,:)            ::  inflow_tsr,inflow_tsi   !Eigenfunction: rho, u, v, T, p
        
        !Define free-stream wave
        integer(KINT)                                       ::  Num_dis_free            !number of free-stream wave
        type(wave),dimension(20)                            ::  wave_free               !free-stream wave (maximum 20, could be changed)
        
        !Fringe region
        integer(KINT)                                       ::  Iters_fringe
        real(KREAL)                                         ::  rfringe(4)
        real(KREAL),allocatable,dimension(:,:)              ::  Xlamd
        
    end module Flow_DIS

# Standalone native Fortran/MPI build. No Python required.
FC = mpifort
CC = cc
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
TARGET = $(BIN_DIR)/ugks_solver
MAIN_SRC = $(SRC_DIR)/main_UGKS2D.f90
MOD_SRCS := $(filter-out $(MAIN_SRC),$(wildcard $(SRC_DIR)/*.f90))
MODULE_OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MOD_SRCS)) $(OBJ_DIR)/mpi2_platform.o
MAIN_OBJ = $(OBJ_DIR)/main_UGKS2D.o
OBJS = $(MODULE_OBJS) $(MAIN_OBJ)
COMPILER ?= gnu
ifneq (,$(filter mpiifort mpiifx ifort ifx,$(notdir $(FC))))
COMPILER = intel
endif
ifeq ($(COMPILER),intel)
COMMON_FLAGS = -fpp -module $(OBJ_DIR)
RELEASE_FLAGS = -O2
DEBUG_FLAGS = -O0 -g -check all -check noarg_temp_created -traceback
else ifeq ($(COMPILER),gnu)
COMMON_FLAGS = -cpp -J$(OBJ_DIR) -ffree-line-length-none -fall-intrinsics
RELEASE_FLAGS = -O2
DEBUG_FLAGS = -O0 -g -fcheck=all -fbacktrace
else
$(error COMPILER must be gnu or intel)
endif
FLAGS ?= $(COMMON_FLAGS) $(RELEASE_FLAGS)
CFLAGS ?= -O2 -std=c99
LDFLAGS ?=
# SDK linker path required by some macOS Command Line Tools installations.
ifeq ($(shell uname -s),Darwin)
SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)
ifneq ($(strip $(SDKROOT)),)
export SDKROOT
LDFLAGS += -L$(SDKROOT)/usr/lib
endif
endif
.DEFAULT_GOAL := all
.PHONY: all release debug clean
all: $(TARGET)
release:
	$(MAKE) clean
	$(MAKE) FLAGS='$(COMMON_FLAGS) $(RELEASE_FLAGS)' $(TARGET)
debug:
	$(MAKE) clean
	$(MAKE) FLAGS='$(COMMON_FLAGS) $(DEBUG_FLAGS)' $(TARGET)
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FLAGS) $(LDFLAGS) -o $@ $^
$(MAIN_OBJ): $(MODULE_OBJS)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FLAGS) -c $< -o $@
$(OBJ_DIR)/mpi2_platform.o: $(SRC_DIR)/mpi2_platform.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(TARGET)

# Frozen E3 module dependency graph.
$(OBJ_DIR)/mod_MPI.o: $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_ImplicitCommon.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_ImplicitState.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_ImplicitBoundary.o: $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_ImplicitMacro.o: $(OBJ_DIR)/mod_ImplicitBoundary.o $(OBJ_DIR)/mod_ImplicitCommon.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_MPI.o
$(OBJ_DIR)/mod_ImplicitMicro.o: $(OBJ_DIR)/mod_ImplicitBoundary.o $(OBJ_DIR)/mod_ImplicitCommon.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_MPI.o $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_ImplicitDiagnostics.o: $(OBJ_DIR)/mod_ImplicitMicro.o $(OBJ_DIR)/mod_ImplicitCommon.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Tools.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_ImportPara.o: $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_Interpolation.o: $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_Evolution.o: $(OBJ_DIR)/mod_Interpolation.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_WENO5AOScalar.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_WENO5AO2D.o: $(OBJ_DIR)/mod_WENO5AOScalar.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_CorrectedWF2DTypes.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_CorrectedWF2DEquilibrium.o: $(OBJ_DIR)/mod_CorrectedWF2DTypes.o $(OBJ_DIR)/mod_WENO5AO2D.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_CorrectedWF2DClosure.o: $(OBJ_DIR)/mod_CorrectedWF2DEquilibrium.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_S2D1Reconstruction.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_S2D1.o: $(OBJ_DIR)/mod_S2D1Reconstruction.o $(OBJ_DIR)/mod_CorrectedWF2DClosure.o $(OBJ_DIR)/mod_CorrectedWF2DEquilibrium.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_S2NativeAdapter.o: $(OBJ_DIR)/mod_TimeIntegrator.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_S2SteadyTermination.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Residual.o: $(OBJ_DIR)/mod_S2NativeAdapter.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_SteadyImplicit.o: $(OBJ_DIR)/mod_ImplicitDiagnostics.o $(OBJ_DIR)/mod_ImplicitMacro.o $(OBJ_DIR)/mod_ImplicitMicro.o $(OBJ_DIR)/mod_Residual.o
$(OBJ_DIR)/mod_S2SteadyIteration.o: $(OBJ_DIR)/mod_SteadyImplicit.o $(OBJ_DIR)/mod_S2NativeAdapter.o $(OBJ_DIR)/mod_ImplicitDiagnostics.o $(OBJ_DIR)/mod_ImplicitMacro.o $(OBJ_DIR)/mod_ImplicitMicro.o $(OBJ_DIR)/mod_ImplicitBoundary.o $(OBJ_DIR)/mod_ImplicitState.o
$(OBJ_DIR)/mod_TimeIntegrator.o: $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_UnsteadyImplicit.o: $(OBJ_DIR)/mod_SteadyImplicit.o $(OBJ_DIR)/mod_ImplicitMacro.o $(OBJ_DIR)/mod_ImplicitMicro.o $(OBJ_DIR)/mod_Residual.o $(OBJ_DIR)/mod_Interpolation.o $(OBJ_DIR)/mod_Evolution.o $(OBJ_DIR)/mod_Boundary.o $(OBJ_DIR)/mod_MPI.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Boundary.o: $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_InitializeFlow.o: $(OBJ_DIR)/mod_S2SteadyIteration.o $(OBJ_DIR)/mod_S2NativeAdapter.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_UnsteadyImplicit.o $(OBJ_DIR)/mod_SteadyImplicit.o $(OBJ_DIR)/mod_ImplicitBoundary.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Boundary.o $(OBJ_DIR)/mod_MPI.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Output.o: $(OBJ_DIR)/mod_S2NativeAdapter.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Evolution.o: $(OBJ_DIR)/mod_S2NativeAdapter.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_Interpolation.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Solve.o: $(OBJ_DIR)/mod_S2SteadyIteration.o $(OBJ_DIR)/mod_S2NativeAdapter.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_UnsteadyImplicit.o $(OBJ_DIR)/mod_InitializeFlow.o $(OBJ_DIR)/mod_SteadyImplicit.o $(OBJ_DIR)/mod_Output.o $(OBJ_DIR)/mod_Evolution.o $(OBJ_DIR)/mod_Residual.o $(OBJ_DIR)/mod_TimeIntegrator.o $(OBJ_DIR)/mod_Interpolation.o $(OBJ_DIR)/mod_Boundary.o $(OBJ_DIR)/mod_Tools.o $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_MPI.o
$(OBJ_DIR)/main_UGKS2D.o: $(OBJ_DIR)/mod_S2D1.o
$(OBJ_DIR)/mod_S2MultiscaleWall.o: $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_CorrectedWF2DClosure.o
$(OBJ_DIR)/mod_S2D1.o: $(OBJ_DIR)/mod_S2MultiscaleWall.o
$(OBJ_DIR)/mod_W3BObserver.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_Solve.o: $(OBJ_DIR)/mod_W3BObserver.o
$(OBJ_DIR)/mod_VelocityTables.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_S2D1.o: $(OBJ_DIR)/mod_VelocityTables.o
$(OBJ_DIR)/mod_ImportPara.o: $(OBJ_DIR)/mod_VelocityTables.o
$(OBJ_DIR)/mod_InitializeFlow.o: $(OBJ_DIR)/mod_VelocityTables.o
$(OBJ_DIR)/mod_S2PhysicalGas.o: $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_S2MultiscaleWall.o $(OBJ_DIR)/mod_S2D1.o: $(OBJ_DIR)/mod_S2PhysicalGas.o
$(OBJ_DIR)/mod_S2MultiscaleWallLegacy.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_S2MultiscaleWall.o: $(OBJ_DIR)/mod_S2MultiscaleWallLegacy.o
$(OBJ_DIR)/mod_S2MPI1Layout.o: $(OBJ_DIR)/mod_Globaldata.o
$(OBJ_DIR)/mod_S2D1Reconstruction.o $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_ImportPara.o: $(OBJ_DIR)/mod_S2MPI1Layout.o
$(OBJ_DIR)/mod_S2MPI1.o: $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_MPI.o
$(OBJ_DIR)/main_UGKS2D.o: $(OBJ_DIR)/mod_S2MPI1.o
$(OBJ_DIR)/mod_S2MPI1.o: $(OBJ_DIR)/mod_InitializeFlow.o
$(OBJ_DIR)/mod_S2MPI2IO.o: $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_S2MPI1Layout.o
$(OBJ_DIR)/mod_S2MPI1.o: $(OBJ_DIR)/mod_S2MPI2IO.o
$(OBJ_DIR)/mod_MPI.o: $(OBJ_DIR)/mod_S2MPI1Layout.o
$(OBJ_DIR)/mod_S2StateConventions.o: $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_S2NativeAdapter.o: $(OBJ_DIR)/mod_S2MPI2IO.o
$(OBJ_DIR)/mod_S2NativeRestart.o: $(OBJ_DIR)/mod_S2MPI2IO.o
$(OBJ_DIR)/mod_S2NativeAdapter.o: $(OBJ_DIR)/mod_S2NativeRestart.o
$(OBJ_DIR)/mod_S2CartesianMPI.o: $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_Globaldata.o $(OBJ_DIR)/mod_S2MPI1Layout.o
$(OBJ_DIR)/mod_S2NativeAdapter.o: $(OBJ_DIR)/mod_S2CartesianMPI.o
$(OBJ_DIR)/mod_TimeIntegrator.o $(OBJ_DIR)/mod_Output.o $(OBJ_DIR)/mod_Solve.o: $(OBJ_DIR)/mod_S2MPI1Layout.o

$(OBJ_DIR)/mod_ClosedCouetteMass.o: $(OBJ_DIR)/mod_Boundary.o $(OBJ_DIR)/mod_MPI.o $(OBJ_DIR)/mod_ImplicitCommon.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Tools.o
$(OBJ_DIR)/mod_R5Observe.o: $(OBJ_DIR)/mod_ClosedCouetteMass.o $(OBJ_DIR)/mod_ImplicitDiagnostics.o $(OBJ_DIR)/mod_Evolution.o $(OBJ_DIR)/mod_Residual.o $(OBJ_DIR)/mod_ImplicitMicro.o
$(OBJ_DIR)/mod_SteadyImplicit.o: $(OBJ_DIR)/mod_R5Observe.o
$(OBJ_DIR)/mod_InitializeFlow.o $(OBJ_DIR)/mod_Output.o: $(OBJ_DIR)/mod_ClosedCouetteMass.o
$(OBJ_DIR)/mod_TimeIntegrator.o: $(OBJ_DIR)/mod_ImplicitCommon.o
$(OBJ_DIR)/mod_Residual.o: $(OBJ_DIR)/mod_ImplicitMicro.o

$(OBJ_DIR)/mod_NativeWiringChecks.o: $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_ImplicitCommon.o
$(OBJ_DIR)/mod_Solve.o: $(OBJ_DIR)/mod_NativeWiringChecks.o

$(OBJ_DIR)/mod_Evolution.o: $(OBJ_DIR)/mod_ImplicitCommon.o

$(OBJ_DIR)/mod_Output.o: $(OBJ_DIR)/mod_ImplicitDiagnostics.o

# IO-R2 residual coefficient serialization.
$(OBJ_DIR)/mod_Solve.o: $(OBJ_DIR)/mod_ClosedCouetteMass.o

$(OBJ_DIR)/mod_S2D1Reconstruction.o: $(OBJ_DIR)/mod_Interpolation.o

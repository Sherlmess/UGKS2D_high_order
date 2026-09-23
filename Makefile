# ==========================================
# 智能双编译器 Makefile (gfortran & ifort)
# 使用说明:
#   1. 本地开发 (Mac): 直接输入 make
#   2. 超算运行: 输入 make FC=ifort
# ==========================================

# 1. 默认编译器设置 (使用 MPI 封装编译器，可被外部 override)
FC = mpifort

# 2. 目录设置
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
TARGET = $(BIN_DIR)/ugks_solver
TEST_DIR = tests/unit
TEST_BIN_DIR = bin/tests

# 3. 源文件扫描（模块 + 主程序）
MAIN_SRC := $(SRC_DIR)/main_UGKS2D.f90
MOD_SRCS := $(filter-out $(MAIN_SRC), $(wildcard $(SRC_DIR)/*.f90))
SRCS := $(MOD_SRCS) $(MAIN_SRC)

# 将 .f90 替换为 .o 文件
MODULE_OBJS = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(MOD_SRCS))
MAIN_OBJ = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(MAIN_SRC))
OBJS = $(MODULE_OBJS) $(MAIN_OBJ)
TEST_SUPPORT_OBJ = $(OBJ_DIR)/mod_TestSupport.o
UNIT_TEST_NAMES = test_implicit_config test_implicit_state test_residual_rates test_implicit_common test_implicit_gh16 test_implicit_boundary test_implicit_macro test_implicit_micro test_steady_implicit_update test_unsteady_baseline test_unsteady_time_scales test_unsteady_macro test_unsteady_micro
UNIT_TEST_BINS = $(addprefix $(TEST_BIN_DIR)/,$(UNIT_TEST_NAMES))
MPI_TEST_NAMES = test_implicit_mpi
MPI_TEST_BINS = $(addprefix $(TEST_BIN_DIR)/,$(MPI_TEST_NAMES))

# ==========================================
# 4. 智能编译器配置
# ==========================================

# --- 配置 A: gfortran (你的 Mac) ---
ifeq ($(FC), gfortran)
    # 输出 .mod 到 obj
    MOD_FLAG = -J$(OBJ_DIR)
    # 严厉开发模式
    DEBUG_FLAGS = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -std=f2008 -pedantic -ffree-line-length-none -fall-intrinsics
    # 发布模式
    RELEASE_FLAGS = -O3 -march=native -funroll-loops
endif

# --- 配置 A-2: mpifort (实为 gfortran + MPI 封装) ---
ifeq ($(FC), mpifort)
    MOD_FLAG = -J$(OBJ_DIR)
    DEBUG_FLAGS = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -ffpe-summary=all \
                  -std=f2008 -pedantic -ffree-line-length-none -fall-intrinsics
    RELEASE_FLAGS = -O3 -march=native -funroll-loops
endif

# --- 配置 B: ifort (超算) ---
ifeq ($(FC), ifort)
    # 输出 .mod 到 obj (注意这里是 -module)
    MOD_FLAG = -module $(OBJ_DIR)
    DEBUG_FLAGS = -g -O0 -check all -traceback -fpe0
    RELEASE_FLAGS = -O3 -xHost -ipo -qopenmp
endif

# --- 配置 C: Intel MPI 封装编译器 (超算，纯 MPI) ---
ifeq ($(FC), mpiifort)
    MOD_FLAG = -module $(OBJ_DIR)
    DEBUG_FLAGS = -g -O0 -check all -traceback -fpe0
    RELEASE_FLAGS = -O3 -xHost -ipo
endif

# 通用预处理参数
COMMON_FLAGS = -cpp $(MOD_FLAG)

# 默认使用 Debug 模式
FLAGS = $(COMMON_FLAGS) $(DEBUG_FLAGS)

# ==========================================
# 5. 构建规则
# ==========================================

.PHONY: all clean debug release test-source test-unit test-regression test-performance

all: debug

debug:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory FC=$(FC) FLAGS='$(COMMON_FLAGS) $(DEBUG_FLAGS)' $(TARGET)

release:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory FC=$(FC) FLAGS='$(COMMON_FLAGS) $(RELEASE_FLAGS)' $(TARGET)

test-source:
	python3 tests/helpers/check_intel_mpi_module_import.py
	python3 tests/helpers/check_unsteady_baseline_source.py
	python3 tests/helpers/check_shakhov_admissibility_policy.py

test-unit: test-source debug
	@$(MAKE) --no-print-directory FC=$(FC) FLAGS='$(COMMON_FLAGS) $(DEBUG_FLAGS)' $(UNIT_TEST_BINS) $(MPI_TEST_BINS)
	@set -e; for test_bin in $(UNIT_TEST_BINS); do mpirun -np 1 $$test_bin; done
	@set -e; for test_bin in $(MPI_TEST_BINS); do mpirun -np 2 $$test_bin; done

test-regression: test-unit
	bash tests/regression/run_explicit_structured.sh
	bash tests/regression/run_implicit_uniform.sh
	bash tests/regression/run_implicit_unsteady_uniform.sh
	bash tests/regression/run_implicit_unsteady_fwave.sh
	bash tests/regression/run_implicit_unsteady_acc06.sh
	bash tests/regression/run_implicit_unsteady_temporal.sh
	bash tests/regression/run_implicit_unsteady_temporal_order.sh
	bash tests/regression/run_flatplate_steady_compare.sh
	bash tests/regression/run_flatplate_maxwell_accommodation.sh

test-performance: release
	bash tests/performance/run_ncfl_sweep.sh

# 链接
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	@echo "? [$(FC)] Linking..."
	$(FC) $(FLAGS) -o $@ $^
	@echo "? Build successful! Executable: $(TARGET)"

# main 程序在编译前需要已有模块文件
$(MAIN_OBJ): $(MODULE_OBJS)

# 模块依赖关系（确保 .mod 生成顺序正确）
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
$(OBJ_DIR)/mod_S2NativeAdapter.o: $(OBJ_DIR)/mod_S2D1.o $(OBJ_DIR)/mod_S2SteadyTermination.o $(OBJ_DIR)/mod_ImplicitState.o $(OBJ_DIR)/mod_Globaldata.o
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

$(TEST_SUPPORT_OBJ): $(TEST_DIR)/mod_TestSupport.f90 $(OBJ_DIR)/mod_Globaldata.o
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FLAGS) -I$(OBJ_DIR) -c $< -o $@

$(TEST_BIN_DIR)/%: $(TEST_DIR)/%.f90 $(MODULE_OBJS) $(TEST_SUPPORT_OBJ)
	@mkdir -p $(TEST_BIN_DIR)
	$(FC) $(FLAGS) -I$(OBJ_DIR) -o $@ $< $(TEST_SUPPORT_OBJ) $(MODULE_OBJS)
# 编译
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJ_DIR)
	@echo "? [$(FC)] Compiling $<..."
	$(FC) $(FLAGS) -c $< -o $@

# 清理
clean:
	@echo "? Cleaning up..."
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(TARGET) $(UNIT_TEST_BINS) $(MPI_TEST_BINS)

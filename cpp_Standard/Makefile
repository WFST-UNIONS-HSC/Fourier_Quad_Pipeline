CXX = mpicxx
MPIRUN ?= mpirun
STACK_PREFIX ?=
EIGEN_INCLUDE ?=

CPPFLAGS = -Iinclude -Iconfig -Iinclude/process_extcat -Iinclude/process_init -Iinclude/process_main -Iinclude/process_rearr -Iinclude/process_fd
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -MMD -MP
LDFLAGS =
LDLIBS = -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm

ifneq ($(strip $(STACK_PREFIX)),)
CPPFLAGS += -I$(STACK_PREFIX)/include -I$(STACK_PREFIX)/include/eigen3
LDFLAGS += -L$(STACK_PREFIX)/lib -Wl,-rpath,$(STACK_PREFIX)/lib
endif

ifneq ($(strip $(EIGEN_INCLUDE)),)
CPPFLAGS += -I$(EIGEN_INCLUDE)
endif

SRCS = main.cpp \
       src/process_extcat/process_extcat.cpp \
       src/process_rearr/process_rearr.cpp \
       src/process_rearr/CatalogRearranger.cpp \
       src/process_fd/process_fd.cpp \
       src/process_fd/ShearCatalogReader.cpp \
       src/process_fd/StarCutCalculator.cpp \
       src/process_fd/KMeansClusterer.cpp \
       src/process_fd/FDMeasurement.cpp \
       src/process_fd/QuadraticFitting.cpp \
       src/process_main/process_main.cpp \
       src/process_main/ExternalCatalogReader.cpp \
       src/process_init/process_init.cpp \
       src/process_init/Initializer.cpp \
       src/process_init/FitsExtractor.cpp \
       src/process_main/PreProcess.cpp \
       src/process_main/Astrometry.cpp \
       src/process_main/SourceExtractor.cpp \
       src/process_main/FourierTransformSt1.cpp \
       src/process_main/FourierTransformSt2.cpp \
       src/process_main/PSFModel.cpp \
       src/process_main/PSFRecons.cpp \
       src/process_main/ShearMeasurement.cpp \
       src/process_main/ExposureInfo.cpp \
       src/process_main/CatalogCombiner.cpp \
       src/process_main/FitsIO.cpp \
       src/process_main/LinearSolve.cpp \
       src/process_main/Universalblock.cpp \
       src/process_main/UniversalUtils.cpp \
       src/process_main/ImageProcessing.cpp \
       src/process_main/NumericalRecipes.cpp \
       src/process_main/MPIScheduler.cpp \
       src/process_main/ExStar.cpp

OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
TARGET = Fourier_Quad_Pipe
UNIVERSALBLOCK_TEST = tests/test_universalblock
CATALOG_ROW_COUNT_TEST = tests/test_catalog_row_count
PSF_MODEL_STATE_TEST = tests/test_psf_model_state
LEGACY_STAMP_CAPACITY_TEST = tests/test_legacy_stamp_capacity
MPI_FAILURE_TEST = tests/test_mpi_failure
ASTROMETRY_DYNAMIC_CAPACITY_TEST = tests/test_astrometry_dynamic_capacity
CATALOG_LIFECYCLE_TEST = tests/test_catalog_lifecycle
EXPOSURE_RUNTIME_SIZING_TEST = tests/test_exposure_runtime_sizing

all: $(TARGET)

test-universalblock: $(UNIVERSALBLOCK_TEST)
	./$(UNIVERSALBLOCK_TEST)

test-catalog-row-count: $(CATALOG_ROW_COUNT_TEST)
	./$(CATALOG_ROW_COUNT_TEST)

test-psf-model-state: $(PSF_MODEL_STATE_TEST)
	./$(PSF_MODEL_STATE_TEST)

test-legacy-stamp-capacity: $(LEGACY_STAMP_CAPACITY_TEST)
	./$(LEGACY_STAMP_CAPACITY_TEST)

test-astrometry-dynamic-capacity: $(ASTROMETRY_DYNAMIC_CAPACITY_TEST)
	./$(ASTROMETRY_DYNAMIC_CAPACITY_TEST)

test-mpi-failure: $(MPI_FAILURE_TEST)
	@status=0; \
	OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
	timeout 10 $(MPIRUN) -np 2 ./$(MPI_FAILURE_TEST) || status=$$?; \
	if [ $$status -eq 0 ] || [ $$status -eq 124 ]; then \
		echo "MPI failure propagation test failed with status $$status"; \
		exit 1; \
	fi; \
	echo "MPI failure propagation test passed with status $$status"

test-catalog-lifecycle: $(CATALOG_LIFECYCLE_TEST)
	./$(CATALOG_LIFECYCLE_TEST)

test-exposure-runtime-sizing: $(EXPOSURE_RUNTIME_SIZING_TEST)
	./$(EXPOSURE_RUNTIME_SIZING_TEST)

$(UNIVERSALBLOCK_TEST): tests/UniversalblockTest.cpp \
                       src/process_main/Universalblock.cpp \
                       src/process_main/UniversalUtils.cpp \
                       src/process_main/FitsIO.cpp \
                       src/process_main/LinearSolve.cpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/UniversalblockTest.cpp \
		src/process_main/Universalblock.cpp \
		src/process_main/UniversalUtils.cpp \
		src/process_main/FitsIO.cpp \
		src/process_main/LinearSolve.cpp \
		$(LDFLAGS) -o $@ $(LDLIBS)

$(CATALOG_ROW_COUNT_TEST): tests/CatalogRowCountTest.cpp \
                          include/process_main/CatalogRowCount.hpp \
                          include/process_main/MPIFailure.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/CatalogRowCountTest.cpp -o $@

$(PSF_MODEL_STATE_TEST): tests/PSFModelStateTest.cpp \
                        include/process_main/PSFModelState.hpp \
                        include/process_main/MPIFailure.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/PSFModelStateTest.cpp -o $@

$(LEGACY_STAMP_CAPACITY_TEST): tests/LegacyStampCapacityTest.cpp \
                              src/process_main/FitsIO.cpp \
                              include/process_main/FitsIO.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/LegacyStampCapacityTest.cpp \
		src/process_main/FitsIO.cpp $(LDFLAGS) -o $@ -lcfitsio

$(MPI_FAILURE_TEST): tests/MPIFailureTest.cpp \
                    include/process_main/MPIFailure.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
	tests/MPIFailureTest.cpp -o $@

$(ASTROMETRY_DYNAMIC_CAPACITY_TEST): tests/AstrometryDynamicCapacityTest.cpp \
                                    src/process_main/Astrometry.cpp \
                                    src/process_main/UniversalUtils.cpp \
                                    src/process_main/LinearSolve.cpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		-ffunction-sections -fdata-sections \
		tests/AstrometryDynamicCapacityTest.cpp \
		src/process_main/Astrometry.cpp \
		src/process_main/UniversalUtils.cpp \
		src/process_main/LinearSolve.cpp $(LDFLAGS) \
		-Wl,--gc-sections -o $@

$(CATALOG_LIFECYCLE_TEST): tests/CatalogCombinerLifecycleTest.cpp \
                          src/process_main/CatalogCombiner.cpp \
                          src/process_main/Universalblock.cpp \
                          src/process_main/UniversalUtils.cpp \
                          src/process_main/FitsIO.cpp \
                          src/process_main/LinearSolve.cpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/CatalogCombinerLifecycleTest.cpp \
		src/process_main/CatalogCombiner.cpp \
		src/process_main/Universalblock.cpp \
		src/process_main/UniversalUtils.cpp \
		src/process_main/FitsIO.cpp \
		src/process_main/LinearSolve.cpp \
		$(LDFLAGS) -o $@ $(LDLIBS)

$(EXPOSURE_RUNTIME_SIZING_TEST): tests/ExposureRuntimeSizingTest.cpp \
                                include/process_main/ExposureInfo.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/ExposureRuntimeSizingTest.cpp -o $@

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET) $(UNIVERSALBLOCK_TEST) \
	      $(CATALOG_ROW_COUNT_TEST) $(PSF_MODEL_STATE_TEST) \
	      $(LEGACY_STAMP_CAPACITY_TEST) $(MPI_FAILURE_TEST) \
	      $(ASTROMETRY_DYNAMIC_CAPACITY_TEST) $(CATALOG_LIFECYCLE_TEST) \
	      $(EXPOSURE_RUNTIME_SIZING_TEST)

-include $(DEPS)

.PHONY: all clean test-universalblock test-catalog-row-count \
	        test-psf-model-state test-legacy-stamp-capacity \
	        test-astrometry-dynamic-capacity test-mpi-failure \
	        test-catalog-lifecycle test-exposure-runtime-sizing

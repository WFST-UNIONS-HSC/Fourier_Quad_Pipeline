CXX = mpicxx
STACK_PREFIX ?=

CPPFLAGS =
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -MMD -MP
LDFLAGS =
LDLIBS = -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm

ifneq ($(strip $(STACK_PREFIX)),)
CPPFLAGS += -I$(STACK_PREFIX)/include -I$(STACK_PREFIX)/include/eigen3
LDFLAGS += -L$(STACK_PREFIX)/lib -Wl,-rpath,$(STACK_PREFIX)/lib
endif

SRCS = main.cpp \
       PreProcess.cpp \
       Astrometry.cpp \
       SourceExtractor.cpp \
       FourierTransformSt1.cpp \
       FourierTransformSt2.cpp \
       PSFModel.cpp \
       PSFRecons.cpp \
       ShearMeasurement.cpp \
       ExposureInfo.cpp \
       CatalogCombiner.cpp \
       FitsIO.cpp \
       LinearSolve.cpp \
       UniversalUtils.cpp \
       ImageProcessing.cpp \
       NumericalRecipes.cpp \
       MPIScheduler.cpp \
       ExStar.cpp

OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
TARGET = Fourier_Quad_Pipe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)

.PHONY: all clean

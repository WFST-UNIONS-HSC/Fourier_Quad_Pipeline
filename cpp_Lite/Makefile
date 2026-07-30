LIB_PATH_ROOT = /lustre/home/acct-phyzj/phyzj-m31/ysx/pack
INC_DIRS = -I$(LIB_PATH_ROOT)/include -I$(LIB_PATH_ROOT)/include/eigen3
LIB_DIRS = -L$(LIB_PATH_ROOT)/lib
RPATH = -Wl,-rpath,$(LIB_PATH_ROOT)/lib

CXX = mpicxx
CXXFLAGS = -O3 -std=c++17 $(INC_DIRS) -Wall
LIBS = -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm

SRCS = main.cpp \
       PreProcess.cpp \
       Astrometry.cpp \
       SourceExtractor.cpp \
       FourierTransformSt1.cpp \
       FourierTransformSt2.cpp \
       PSFModel.cpp \
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
TARGET = Fourier_Quad_Main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIB_DIRS) $(LIBS) $(RPATH)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean

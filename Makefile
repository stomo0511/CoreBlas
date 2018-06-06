#
UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
  TMATRIX_ROOT = /home/stomo/WorkSpace/TileMatrix
  CXX =	g++
endif
ifeq ($(UNAME),Darwin)
  TMATRIX_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/TileMatrix
  CXX =	/usr/local/bin/g++ 
endif
#
BLAS_ROOT = /opt/OpenBLAS
BLAS_INC_DIR = $(BLAS_ROOT)/include
BLAS_LIB_DIR = $(BLAS_ROOT)/lib
BLAS_LIBS = -lopenblas_seq
#
PLASMA_ROOT = /opt/plasma-17.1
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lcoreblas -lplasma
#
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTileMatrix
#
CXXFLAGS =	-fopenmp -m64 -O2 -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR)

LOBJS =		CoreBlasTile.o

LIBS = libCoreBlasTile.a

TARGET = test

$(LIBS):	$(LOBJS)
	$(AR) r $(LIBS) $(LOBJS)
	ranlib $(LIBS)

$(TARGET): CoreBlasTileTest.o $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ CoreBlasTileTest.o $(LIBS) -L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(BLAS_LIBS)
	
all:	$(LIBS) $(TARGET)

clean:
	rm -f *.o $(LOBJS) $(LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

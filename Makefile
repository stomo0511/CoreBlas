#
BLAS_ROOT = /opt/OpenBLAS
BLAS_INC_DIR = $(BLAS_ROOT)/include
BLAS_LIB_DIR = $(BLAS_ROOT)/lib
BLAS_LIBS = -lopenblas_seq
#
PLASMA_ROOT = /opt/PLASMA
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lplasma -lcoreblas -lquark -lpthread
#
TMATRIX_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/TileMatrix
#TMATRIX_ROOT = /home/stomo/WorkSpace/TileMatrix
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTileMatrix
#
CXX = /usr/local/bin/g++
#CXXFLAGS = -g -Wall -fopenmp -DDEBUG -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR)
CXXFLAGS = -O3 -fopenmp -DDEBUG -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR)

LOBJS =		CoreBlasTile.o

LIBS = libCoreBlasTile.a

TARGET = test

$(LIBS):	$(LOBJS)
	$(AR) r $(LIBS) $(LOBJS)
	ranlib $(LIBS)

$(TARGET): CoreBlasTileTest.o $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ CoreBlasTileTest.o $(LIBS) -L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(BLAS_LIBS)
	
all:	$(LIBS)

clean:
	rm -f *.o $(LOBJS) $(LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

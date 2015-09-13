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
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTMatrix
#
CXX = /usr/local/bin/g++
CXXFLAGS = -g -Wall -fopenmp -DDEBUG -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR)

LOBJS =		CoreBlas.o

LIBS = libCoreBlas.a

TARGET =	test

$(LIBS):	$(LOBJS)
	$(AR) r $(LIBS) $(LOBJS)
	ranlib $(LIBS)

$(TARGET):	CoreBlasTest.o $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ CoreBlasTest.o $(LIBS)

all:	$(LIBS) $(TARGET)

clean:
	rm -f $(LOBJS) $(LIBS) $(TARGET)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

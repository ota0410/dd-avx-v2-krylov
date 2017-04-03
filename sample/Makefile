CXX=g++

CXXFLAGS=-O3 -mavx2 -fopenmp -lddavx

SRC=./main.cpp
OBJ=main.o ../libddavx.a

INCLUDE=-I ../include
LIBDIR=-L ../

all: $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(INCLUDE) $(LIBDIR) -o $(OBJ)

clean:
	@rm -rf $(OBJ) 

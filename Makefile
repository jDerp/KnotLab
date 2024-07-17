CXX = g++
CXXFLAGS = -std=c++14 -O2 -Wall -Wl,-rpath,/usr/local/lib
LDFLAGS = -I${shell pwd}/iRRAMx/include -L${shell pwd}/iRRAMx/lib -liRRAMx -liRRAM -lmpfr -lgmp

all: eigen eigen_semi

%: %.cpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

CXX = g++
CXXFLAGS = -std=c++14 -O2 -Wall -Wl,-rpath,/usr/local/lib
LDFLAGS = -I${shell pwd}/iRRAMx/include -L${shell pwd}/iRRAMx/lib -liRRAMx -liRRAM -lmpfr -lgmp

all: burau burau_approx

%: %.cpp
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

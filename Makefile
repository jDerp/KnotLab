eigen: eigen.cpp
	g++ -std=c++14 -O2 -Wall -Wl,-rpath,/usr/local/lib eigen.cpp -I${shell pwd}/iRRAMx/include -L${shell pwd}/iRRAMx/lib -liRRAMx -liRRAM -lmpfr -lgmp -o eigen
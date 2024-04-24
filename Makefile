eigen: eigen.cpp
	g++ -std=c++14 -O2 -Wall -Wl,-rpath,/usr/local/lib eigen.cpp -I${HOME}/git_dir/iRRAMx/include -L${HOME}/git_dir/iRRAMx/lib -liRRAMx -liRRAM -lmpfr -lgmp -o eigen
# Braid using Exact Real Computation

## Purpose

- Precise checking for Eigenvalues of Braid's Burau representation 


## Library used

- fbrausse/iRRAM
- realcomputation/iRRAMx

## Files

- burau.cpp

Main code for Burau Spectral Test. The test is of arbitrary precision which is semi-decidable (halts when irreducible/unexchangeable, loops infinitely otherwise). 

- burau_approx.cpp

Code for soft Burau Spectral test. The test is of a given precision which will output whether there is any obstruction up to that precision.

- makefile

Run make to compile burau and burau_approx. It can also be used to compile other iRRAM/iRRAMx program by run ```make [filename]``` e.g. if your file name is code123.cpp, run ```make code123```.

- testcases

Testcases are formatted for burau.cpp in which "ir" is irreducible, "r" is reducible, "ne" is unexchangeable and "e" is exchangeable.




  

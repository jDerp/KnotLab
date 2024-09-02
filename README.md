# Braid using Exact Real Computation

## Purpose

- Precise checking for Eigenvalues of Braid's Burau representation 


## Library used

- fbrausse/iRRAM
- realcomputation/iRRAMx

## Installation

Recursively clone the Spectral Test Repository to the desired directory, in which the iRRAM and iRRAMx will be automatically be cloned by the submodule links. 

```
git clone --recursive https://github.com/jDerp/KnotLab
cd KnotLab
```

Since iRRAMx requires iRRAM to be installed system-wide, iRRAM is required to be installed system-wide. Note that you would require admin privileges to install. The following is a step-by-step guide on how to install iRRAM system-wide for ubuntu linux:

Install the required packages via the following command.

```
sudo apt-get install libgmp-dev libmpfr-dev libpng-dev autoconf gcc g++ make libtool
```

Move into the iRRAM directory and run quickinstall.

```
sudo ./path/to/KnotLab/iRRAM/QUICKINSTALL_run_me
```

Let the quickinstall download GMP and MPFR from Trier and install locally in the iRRAM folder.

To do this, select 2), use default directory, then n, y, y for both GMP and MPFR.

After this step, iRRAM should be installed in your system.

To install iRRAMx, move to iRRAMx directory and run make.

```
cd path/to/KnotLab/iRRAMx
make
```

After this step, iRRAMx should be installed locally.

Then, running make on the original directory should compile the Spectral program successfully.

```
cd path/to/KnotLab
make
```


## Files

- burau.cpp

Main code for Burau Spectral Test. The test is of arbitrary precision which is semi-decidable (halts when irreducible/unexchangeable, loops infinitely otherwise). 

- burau_approx.cpp

Code for soft Burau Spectral test. The test is of a given precision which will output whether there is any obstruction up to that precision.

- makefile

Run make to compile burau and burau_approx. It can also be used to compile other iRRAM/iRRAMx program by run ```make [filename]``` e.g. if your file name is code123.cpp, run ```make code123```.

- testcases

Testcases are formatted for burau.cpp in which "ir" is irreducible, "r" is reducible, "ne" is unexchangeable and "e" is exchangeable.




  

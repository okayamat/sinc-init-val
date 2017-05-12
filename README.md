# sinc-init-val
Numerical solvers of initial value problems by means of Sinc methods

## Overview
These programs solve the following four examples of initial value problems.

(1) Find y(t) and z(t) satisfying for 0<t<2 that
* y'(t) = z(t),   y(0) = Ai(0),
* z'(t) = t y(t), z(0) = Ai'(0),  
where Ai(t) is the Airy function.

(2) Find y(t) and z(t) satisfying for 0<t<1 that
* y'(t) = z(t),                 y(0) = 0,
* z'(t) = 2 y(t) / (1 + x^2)^2, z(0) = 1.

(3) Find y(t) and z(t) satisfying for 0<t<2 that
* y'(t) = - y(t) + z(t) / (2 sqrt(t)), y(0) = 0,
* z'(t) = - y(t) / sqrt(t),              z(0) = 1.

(4) Find y(t) and z(t) satisfying for -1<t<1 that
* y'(t) = - 2[t F^2(t) + sin(4 arctanh t)] z(t) / F(t), y(-1) = 0,
* z'(t) =   2[t F^2(t) + sin(4 arctanh t)] y(t) / F(t), z(-1) = 1,  
where F(t) = sqrt(cos(4 arctanh t) + cosh(pi)).

Those problems are solved by means of the following 4 methods:
* SE-Sinc-Nyström method [1]
* SE-Sinc-collocation method [2]
* DE-Sinc-Nyström method [1]
* DE-Sinc-collocation method [3]

Each program solves those problems increasing N as N = 1, 9, 17, 25, ...,
and outputs N, computation time, and maximum error over the target interval.

## How to compile
These programs need GSL (Gnu Scientific Library) and CPPLapack.
Install them first (and libraries needed to them). Then, modify
make files according to your installation.

## Results
Outputs by those programs are stored in data/ directory, with .dat extension.
Gnuplot programs and created eps graphs are also stored in the directory.

computation environment:
* OS: Mac OS X 10.6
* CPU: two 2.93 GHz 6-Core Intel Xeon
* Memory: 32 GB DDR3 ECC SDRAM
* Compiler: GCC 4.0.1
* Libraries: GSL 1.16, LAPACK 3.4.2, CPPLapack 20050325, installed by Fink

## References
1. A. Nurmuhammad, M. Muhammad, M. Mori:
 Numerical solution of initial value problems based on the double exponential
 transformation, Publ. Res. Inst. Math. Sci., Kyoto Univ., Vol. 41 (2005),
 pp. 937--948.
2. F. Stenger: Numerical Methods Based on Sinc and Analytic Functions,
 Springer-Verlag, New York, 1993.
3. T. Okayama: Theoretical analysis of Sinc-collocation methods and
 Sinc-Nystöm methods for initial value problems, BIT Numer. Math., to appear.

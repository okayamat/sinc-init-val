# sinc-init-val
Numerical solvers of initial value problems by means of Sinc methods

## Overview
These programs solve the following four initial value problems.

(1) Find y_1(t) and y_2(t) satisfying for 0<t<2 that
* y_1'(t) = y_2(t),   y_1(0) = Ai(0),
* y_2'(t) = t y_1(t), y_2(0) = Ai'(0),  
where Ai(t) is the Airy function.

(2) Find y_1(t) and y_2(t) satisfying for 0<t<1 that
* y_1'(t) = y_2(t),                 y_1(0) = 0,
* y_2'(t) = 2 y_1(t) / (1 + x^2)^2, y_2(0) = 1.

(3) Find y_1(t) and y_2(t) satisfying for 0<t<2 that
* y_1'(t) = - y_1(t) + y_2(t) / (2 sqrt(t)), y_1(0) = 0,
* y_2'(t) = - y_1(t) / sqrt(t),              y_2(0) = 1.

(4) Find y_1(t) and y_2(t) satisfying for -1<t<1 that
* y_1'(t) = - 2[t F^2(t) + sin(4 arctanh t)] y_2(t)/F(t), y_1(-1) = 0,
* y_2'(t) =   2[t F^2(t) + sin(4 arctanh t)] y_1(t)/F(t), y_2(-1) = 1,  
where F(t) = sqrt(cos(4 arctanh t) + cosh(pi)).

Those problems are solved by means of the following 4 methods:
* SE-Sinc-Nyström method [1]
* SE-Sinc-collocation method [2]
* DE-Sinc-Nyström method [1]
* DE-Sinc-collocation method [3]

Each program solves those problems increasing N as N = 1, 9, 17, 25, ...,
and outputs N, computation time, and maximum error over the target interval.

## How to compile
These programs needs GSL (Gnu Scientific Library) and CPPLapack.
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
* Libraries: GSL 1.16, CPPLapack 20050325, installed by Fink

## References
1. A. Nurmuhammad, M. Muhammad, M. Mori:
 Numerical solution of initial value problems based on the double exponential
 transformation, Publ. Res. Inst. Math. Sci., Kyoto Univ., Vol. 41 (2005),
 pp. 937--948.
2. F. Stenger: Numerical Methods Based on Sinc and Analytic Functions,
 Springer-Verlag, New York, 1993.
3. T. Okayama: Theoretical analysis of Sinc-collocation methods and
 Sinc-Nystöm methods for initial value problems, arXiv:1304.6508 [math.NA].

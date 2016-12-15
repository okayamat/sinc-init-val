# sinc-init-val
Numerical solvers of initial value problems by means of Sinc methods
====
## Overview
These programs solve the following four initial value problems.

[1] Find y_1(t) and y_2(t) satisfying for 0<t<2 that
* y_1'(t) = y_2(t),   y_1(0) = Ai(0),
* y_2'(t) = t y_1(t), y_2(0) = Ai'(0),
where Ai(t) is the Airy function.

[2] Find y_1(t) and y_2(t) satisfying for 0<t<1 that
* y_1'(t) = y_2(t),                 y_1(0) = 0,
* y_2'(t) = 2 y_1(t) / (1 + x^2)^2, y_2(0) = 1.

[3] Find y_1(t) and y_2(t) satisfying for 0<t<2 that
* y_1'(t) = - y_1(t) + y_2(t) / (2 sqrt(t)), y_1(0) = 0,
* y_2'(t) = - y_1(t) / sqrt(t),              y_2(0) = 1.

[4] Find y_1(t) and y_2(t) satisfying for -1<t<1 that
* y_1'(t) = - 2[t F^2(t) + sin(4 arctanh t)] y_2(t)/F(t), y_1(-1) = 0,
* y_2'(t) =   2[t F^2(t) + sin(4 arctanh t)] y_1(t)/F(t), y_2(-1) = 1,
where F(t) = sqrt(cos(4 arctanh t) + cosh(pi)).

Those problems are solved by means of the following 4 methods:
* SE-Sinc-Nystroem method
* SE-Sinc-collocation method
* DE-Sinc-Nystroem method
* DE-Sinc-collocation method

Each program solves those problems with increasing N, and outputs
N, computation time, and maximum error over the target interval.

## How to compile
These programs needs GSL (Gnu Scientific Library) and CPPLapack.
Install them first (and libraries needed to them). Then, modify
make files according to your installation.


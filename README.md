# Physics 4G03
Solutions to problems from Physics 4G03, taught at McMaster University during the Fall 2018 term. Most problems and solutions are implemented in C. Plots are generally generated using matplotlib and pandas in Python.

## Problems and solutions

1. Implementation of Lanczos algorithm in Matlab @ [[lanczos.m]](https://github.com/brindavenkataramani/4G03/blob/master/lanczos.m)
    * Written so I could make sure I wasn't losing my mind when trying to implement the same thing in C.
    * Super barebones; feel free to use it for your own implementations.
2. Utility functions to read/write matrices from/to .txt files in C @ [[matrix_io.c]](https://github.com/brindavenkataramani/4G03/blob/master/matrix_io.c)
    * Just for fun; it might be useful later on (hopefully).
3. Jacobi routine along with some utility functions and a driver program to test the routine @[[jacobi.cpp]](https://github.com/brindavenkataramani/4G03/blob/master/jacobi.cpp)
    * I still need to implement some validations after running the routine.
    * This is probably grossly inefficient, but it works well for Hermitian matrices (up to n = 25+).
    * I also need to implement a sorting method, to sort eigenvalues.

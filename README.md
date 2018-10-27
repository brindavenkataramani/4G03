# Physics 4G03
Solutions to problems from Physics 4G03, taught at McMaster University during the Fall 2018 term. Most problems and solutions are implemented in C. Plots are generally generated using matplotlib and pandas in Python.

## Problems and solutions

1. Implementation of Lanczos algorithm in Matlab @ [[lanczos.m]](https://github.com/brindavenkataramani/4G03/blob/master/lanczos.m)
    * Written so I could make sure I wasn't losing my mind when trying to implement the same thing in C.
    * Super barebones; feel free to use it for your own implementations.
2. Utility functions to read/write matrices from/to .txt files in C @ [[matrix_io.c]](https://github.com/brindavenkataramani/4G03/blob/master/matrix_io.c)
    * Just for fun; it might be useful later on (hopefully).
3. Jacobi routine along with some utility functions and a driver program to test the routine @ [[jacobi.cpp]](https://github.com/brindavenkataramani/4G03/blob/master/jacobi.cpp)
    * I still need to implement some validations after running the routine.
    * This is probably grossly inefficient, but it works well for Hermitian matrices (up to n = 25+).
    * I also need to implement a sorting method, to sort eigenvalues.
    * Also includes some functions to generate a special Hermitian matrix based on the potential $$V(x)=\frac{\hbar^2}{2\mu}\bigg(\frac{0.05}{\sin^2{x}}+\frac{5}{\cos^2{x}}\bigg)$$ (plug that into LaTeX).
4. Working implementation of Lanczos routine in C++ @ [[lanczos.cpp]](https://github.com/brindavenkataramani/4G03/blob/master/lanczos.cpp)
    * Includes driver program to test.
    * Can possibly improve it by moving some housekeeping from the main program to the Lanczos routine itself.
    * More possible improvements include fine-tuning some functions like swap_values() to be more useful in different situations.
    * Validated against Matlab results.
    * Works with Hermitian matrices.
    
    

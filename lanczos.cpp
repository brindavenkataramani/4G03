// ***********************
// LIBRARIES & DEFINITIONS
// ***********************

#include <iostream>
#include <cmath> // Note that we can just use abs instead of fabs in C++.
#include <vector> // For vector structure.
#include <cstdlib>
#include <cstdio>

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix; // Define a matrix structure as vector of rows. Indexing is basically same as 2-dimensional array.

// ******************************
// END OF LIBRARIES & DEFINITIONS
// ******************************

// ************************
// VECTOR UTILITY FUNCTIONS
// ************************

/*

Includes functions to do the following:
  * Compute dot product.
  * Compute matrix-vector product.
  * Copy a vector.
  * Perform scalar multiplication on a vector.

*/

// Working implementation of dot product for two vectors.

double dotprod(const vector<double>& y, const vector<double>& x) { 
  double sum = 0.0;

  for (int i = 0; i < y.size(); i++) {
    sum += x[i] * y[i];

  }

  return (sum);

}

// Working implementation of matrix vector product.

vector<double> matrixprod(Matrix &A, vector<double> &v, int N) {
  int i, j;
  vector<double> temp(N);
  vector<double> res(N);

  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {
      temp[j] = A[i][j];

    }

    res[i] = dotprod(temp, v);

  }

  return res;

}

// Copy a vector w to v.

void swap_vals(vector<double> &v, vector<double> &w, int N) {
  int i;

  for(i = 0; i < N; i++) {
    v[i] = w[i];

  }

}

// Replace v with cv.

void scalar_mult(vector<double> &v, double c, int N) {
  int i;

  for(i = 0; i < N; i++) {
    v[i] = c*v[i];

  }

}

// *******************************
// END OF VECTOR UTILITY FUNCTIONS
// *******************************

// **************************
// STANDARD UTILITY FUNCTIONS
// **************************

/*

Other utility functions:
  * Functions to make symmetric matrices/random vectors.
  * Printing functions.

*/

// Makes a symmetric matrix.

void init_matrix(Matrix &A, int N) {
  int i, j; 

  for(i = 0; i < N; i++) {
    for(j = i; j < N; j++) {
      A[i][j] = rand() % 10; // If we want integers.
      // A[i][j] = (double)rand()/(double)RAND_MAX; // If we want doubles.
      A[j][i] = A[i][j];

    }

  }

}

// Makes a vector.

void init_vector(vector<double> &v, int N) {
  int i;

  for(i = 0; i < N; i++) {
    v[i] = rand() % 10;
    // v[i] = (double)rand()/(double)RAND_MAX;

  }

}

// Prints a matrix.

void print_matrix(Matrix &A, int N) {
  int i, j; 

  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {
      printf("%7.3f ", A[i][j]);

    }

    printf("\n");

  }

}

// Prints a vector.

void print_vector(vector<double> &v, int N) {
  int i;

  for(i = 0; i < N; i++) {
    printf("%f\n", v[i]);

  }

}

// *********************************
// END OF STANDARD UTILITY FUNCTIONS
// *********************************

// ***************
// LANCZOS ROUTINE
// ***************

/* 

The Lanczos routine. Necessary to do some housekeeping in the main loop itself.

*/

void Lanczos(Matrix &A, Matrix &lan, vector<double> &v0, vector<double> &v1, vector<double> &w, vector<double> &f, vector<double> &temp, int N, int M) {
  int i, j;

   for(i = 0; i < M - 1; i++) {
    lan[i][i+1] = sqrt(dotprod(f,f));
    lan[i+1][i] = lan[i][i+1]; //
    scalar_mult(f,1/lan[i][i+1], N); // 
    swap_vals(v1, f, N); //
    swap_vals(temp, v1, N); //
    temp = matrixprod(A, temp, N); //
    
    for(j = 0; j < N; j++) {
      w[j] = temp[j] - lan[i][i+1]*v0[j]; //

    }

    lan[i+1][i+1] = dotprod(v1, w); //

    for(j = 0; j < N; j++) {
      f[j] = w[j] - lan[i+1][i+1]*v1[j];

    }

    swap_vals(v0, v1, N);

  }

}

// **********************
// END OF LANCZOS ROUTINE
// **********************

// **************
// DRIVER PROGRAM
// **************

int main() {
  srand(time(0));
  int N = 10;
  int M = 6;
  int i, j;
  Matrix mat(N, Row(N));
  Matrix lan(M, Row(M));
  vector<double> v0(N), v1(N), w(N), f(N), temp(N);
  
  init_matrix(mat, N);
  print_matrix(mat, N);

  // Housekeeping to start Lanczos.

  init_vector(v0, N);
  scalar_mult(v0, 1/sqrt(dotprod(v0,v0)), N); // Normalize.
  w = matrixprod(mat, v0, N);
  lan[0][0] = dotprod(v0, w);

  for(i = 0; i < N; i++) {
    f[i] = w[i] - lan[0][0]*v0[i];

  }

  // End of housekeeping to start Lanczos.

  Lanczos(mat, lan, v0, v1, w, f, temp, N, M); // Lanczos.

  printf("\n");
  print_matrix(lan, M); // Print tridiagonal matrix.

}

// *********************
// END OF DRIVER PROGRAM
// *********************

/*

This program contains all methods needed to run and test the Jacobi routine for the diagonalization of real symmetric matrices (i.e. Hermitian matrices). It was last updated on 2018-10-26 @ 16:53.

TO-DO:
  * Tune the Jacobi routine (i.e. do validations before returning eigenvector.)
  * Implement sorting method in utility methods section.

*/

// ***********************
// LIBRARIES & DEFINITIONS
// ***********************

#include <iostream>
#include <cmath> // Note that we can just use abs instead of fabs in C++.
#include <vector> // For vector structure.

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix; // Define a matrix structure as vector of rows. Indexing is basically same as 2-dimensional array.

// ******************************
// END OF LIBRARIES & DEFINITIONS
// ******************************

// ***************
// UTILITY METHODS
// ***************

/* 

These are all methods used to test the Jacobi routine. The methods include:
  * A method to initialize a random symmetric matrix [init_matrix()].
  * A method to print out a matrix [print_matrix()].
  * A method to print out a vector; this is for printing the vector of eigenvalues. [print_vector()].
  * A method to sort eigenvalues from least to greatest [eig_sort()].

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

void print_vector(vector<double> v, int N) {
  int i;

  printf("Eigenvalues are:\n");

  for(i = 0; i < N; i++) {
    printf("%f\n", v[i]);

  }

}

// **********************
// END OF UTILITY METHODS
// **********************

// **************
// JACOBI ROUTINE
// **************

// Jacobi diagonalization of a matrix.
// Everything outside this point is still working.
// Nevermind, this is also working.
// What I need to do is tune it up.

/*

TO-DO:
  * Figure out how to break out of a loop if it doesn't converge.
  * Flag that the matrix didn't converge basically.
  * Make sure it works in all cases.

*/

void Jacobi(Matrix &A, vector<double>& v, int N) {
  int maxit = 100; // Number of iterations to try the Jacobi routine..
  double eps = 1.0e-22; // Threshold. Any sum less than this is considered to be 
  double pi = M_PI; // For easier access to pi.
  double phi; // For later (basically rotation angle).
  
  int t, i, j, k; // Counting variables.

  double sum; // Used for validation at the end.

  // Do this until max iterations is reached.

  for(t = 0; t < maxit; t++) {
    double s = 0; // Sum of off-diagonal elements on upper-triangular part.
    
    // Compute the sum.

    for(i = 0; i < N - 1; i++) {
      for(j = i + 1; j < N; j++) {
          s = s + abs(A[i][j]); // Compute absolute sum because otherwise off-diagonal elements may destroy each other even if the matrix is not diagonalized.

      }

    }

    // If it's already diagonal, break out.
    // We can neglect below diagonal elements since the matrix is symmetric.

    if(s < eps) {
      printf("Converged after %i iterations.", t);
      break;

    }

    // Do the rotations otherwise.

    else {
      double limit = s/(N*(N-1)/2.0);

      for(i = 0; i < N - 1; i++) {
        for(j = i + 1; j < N; j++) {
          if(fabs(A[i][j]) > limit) {
            double denom = A[i][i] - A[j][j];
            double phi;

            if(fabs(denom) < eps) {
              phi = pi/4;

            }

            else {
              phi = 0.5*atan(2*A[i][j]/denom);

            }

            double si = sin(phi);
            double co = cos(phi);
            double store;

            for(k = i + 1; k < j; k++) {
              store = A[i][k];
              A[i][k] = A[i][k]*co + A[k][j]*si;
              A[k][j] = A[k][j]*co - store*si;

            }

            for(k = j + 1; k < N; k++) {
              store = A[i][k];
              A[i][k] = A[i][k]*co + A[j][k]*si;
              A[j][k] = A[j][k]*co - store*si;

            }

            for(k = 0; k < i; k++) {
              store = A[k][i];
              A[k][i] = A[k][i]*co + A[k][j]*si;
              A[k][j] = A[k][j]*co - store*si;

            }

            store = A[i][i];
            A[i][i] = A[i][i]*co*co + 2.0*A[i][j]*co*si + A[j][j]*si*si;
            A[j][j] = A[j][j]*co*co - 2.0*A[i][j]*co*si + store*si*si;
            A[i][j] = 0;

          }

        }
      
      }

    }

  }

  // This doesn't seem to work for some reason.
  // Need to figure that out.
  // NEVERMIND; THIS IS FIXED.
  // So what I need to do is return 0s as eigenvalues assuming we can't converge.

  for(i = 0; i < N; i++) {
    v[i] = A[i][i];

  }

}

// *********************
// END OF JACOBI ROUTINE
// *********************

// ******************************
// GENERATION OF SPECIAL MATRICES
// ******************************

/*

Methods to generate special matrices.

*/

// Quick utility function for Dirac delta-function implementation.

int delta_func(int i, int j) {
  if(i == j) {
    return 1;

  }

  else {
    return 0;

  }

}

// Quick utility function to return minimum between two integers. Returns either or if they are equal.

int min(int i, int j) {
  if(i > j) {
    return j;

  }

  else if(i < j) {
    return i;

  }

  else {
    return i;

  }

}

// Function to populate matrix for given potential.
// In units of 2\hbar^2/\mu, where \mu is just mass.

void populate_matrix(Matrix &A, int N) {
  int i, j;

  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {
      A[i][j] = (float)(i+1)*(float)(i+1)*(float)delta_func(i, j) + min(i+1,j+1)*(0.05 + pow(-1,fabs(i-j))*5);

    }

  }

}

// *************************************
// END OF GENERATION OF SPECIAL MATRICES
// *************************************

// **************
// DRIVER PROGRAM
// **************

/*

This contains the main method where I test the Jacobi routine I implemented. I test the validity of its results against the results found using Matlab for the same matrix.

*/

int main() {

  srand(time(0)); // Initialize random seed so we aren't always generating the same matrices over and over again.
  int N = 5; // Size of N x N matrix.
  Matrix mat(N, Row(N)); // Declare N x N matrix.
  vector<double> eig(N); // Declare vector to store eigenvalues.

  /*

  init_matrix(mat, N); // Initialize matrix with random elements.
  print_matrix(mat, N); // Print matrix (validation).
  Jacobi(mat, eig, N); // Diagonalize the matrix using the Jacobi routine.
  printf("\n");
  print_matrix(mat, N); // Print the matrix (validation).
  print_vector(eig, N); // Print the eigenvalues of the matrix.

  */

  populate_matrix(mat, N);
  print_matrix(mat, N);
  Jacobi(mat, eig, N);
  printf("\n");
  print_matrix(mat, N);
  print_vector(eig, N); // Note that the eigenvalues are 1/4 of what they should be. In units of 2\hbar^2/\mu, where \mu is just mass. So actually in units of 4*\hbar^2/2\mu (which is what we want). 

}

// *********************
// END OF DRIVER PROGRAM
// *********************

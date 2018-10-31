// ***********************
// LIBRARIES & DEFINITIONS
// ***********************

#define BIT_SET(a,b) ((a) |= (1U<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1U<<(b))) 
#define BIT_FLIP(a,b) ((a) ^= (1U<<(b)))
#define BIT_CHECK(a,b) ((bool)((a) & (1U<<(b)))) // Set on the condition f else clear
//bool f; // conditional flag
//unsigned int m; // the bit mask
//unsigned int w; // the word to modify: if (f) w |= m; else w &= ~m;
#define COND_BIT_SET(a,b,f) ((a) = ((a) & ~(1U<<(b))) | ((-(unsigned int)f) & (1U<<(b))))

#include <iostream>
#include <cmath> // Note that we can just use abs instead of fabs in C++.
#include <vector> // For vector structure.

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix; // Define a matrix structure as vector of rows. Indexing is basically same as 2-dimensional array.

// ************************
// SPIN HAMILTONIAN ACTION
// ************************

void hv(vector<double>& y, const vector<double>& x, int L) {
  bool b ;
  unsigned int k;
  for (unsigned int i=0;i<x.size();i++){ 
    if (abs(x[i])>2.2e-16) {
      int jm = L-1;
      double xov2=x[i]/2.0; 
      for (int j=0;j<L;j++){
        k=i; 
        COND_BIT_SET(k,jm,BIT_CHECK(i,j)); 
        COND_BIT_SET(k,j,BIT_CHECK(i,jm)); 
        y[k % L] += xov2;
        jm = j;
      } 
    }
  }

  for (unsigned int i=0;i<x.size();i++) {
    y[i]=y[i]-((double) L)/2.0*x[i]/2.0; 

  }

}

// ******************************
// END OF SPIN HAMILTONIAN ACTION
// ******************************

// ************************
// VECTOR UTILITY FUNCTIONS
// ************************

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

void swap_vals(vector<double> &v, vector<double> &w, int N) {
  int i;

  for(i = 0; i < N; i++) {
    v[i] = w[i];

  }

}

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

void Lanczos(Matrix &lan, vector<double> &v0, vector<double> &v1, vector<double> &w, vector<double> &f, vector<double> &temp, int N, int M) {
  int i, j;

   for(i = 0; i < M - 1; i++) {
    lan[i][i+1] = sqrt(dotprod(f,f));
    lan[i+1][i] = lan[i][i+1]; //
    scalar_mult(f,1/lan[i][i+1], N); // 
    swap_vals(v1, f, N); //
    hv(temp, v1, N); //
    
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

int main() {
  // srand(time(0));
  int N = 3;
  int M = 1;
  int i, j;
  Matrix mat(N, Row(N));
  Matrix lan(M, Row(M));
  vector<double> v0(N), v1(N), w(N), f(N), temp(N);
  
  init_matrix(mat, N);
  print_matrix(mat, N);
  init_vector(v0, N);
  print_vector(v0, N);
  scalar_mult(v0, 1/sqrt(dotprod(v0,v0)), N); // Normalize.
  printf("\n");
  print_vector(v0, N);
  hv(w, v0, N); // so the issue is here.
  printf("\n");
  print_vector(v0, N);
  printf("\n");
  print_vector(w, N); 
  lan[0][0] = dotprod(v0, w);
  printf("\n"); 
  printf("%f\n", lan[0][0]);

  for(i = 0; i < N; i++) {
    f[i] = w[i] - lan[0][0]*v0[i];

  }

  printf("\n");
  print_vector(f, N);

  // This is all fine.

  // Now for Lanczos inline.
  // Will implment as function later.

  /*

  for(i = 0; i < M - 1; i++) {
    lan[i][i+1] = sqrt(dotprod(f,f));
    lan[i+1][i] = lan[i][i+1]; //
    scalar_mult(f,1/lan[i][i+1], N); // 
    swap_vals(v1, f, N); //
    swap_vals(temp, v1, N); //
    temp = matrixprod(mat, temp, N); //
    
    for(j = 0; j < N; j++) {
      w[j] = temp[j] - lan[i][i+1]*v0[j]; //

    }

    lan[i+1][i+1] = dotprod(v1, w); //

    for(j = 0; j < N; j++) {
      f[j] = w[j] - lan[i+1][i+1]*v1[j];

    }

    swap_vals(v0, v1, N);

  }

  */

  Lanczos(lan, v0, v1, w, f, temp, N, M);

  printf("\n");
  print_matrix(lan, M);

}

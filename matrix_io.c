# Some functions to read/write matrices from/to .txt files.
# Utility functions for my Lanczos implementation in C.

#include <stdio.h>
#include <stdlib.h>

#define N 5 // So a 5 x 5 matrix.

double matrix[N][N]; // Global matrix variable.

// This is a utility function to initialize a matrix with random elements.

void init_matrix() {
  int i, j;

  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j ++) {
      matrix[i][j] = (double)rand()/(double)RAND_MAX;

    }

  }

}

// This will print a matrix to a file.

void print_matrix() {
  FILE *fp;
  fp = fopen("matrix.txt", "w");

  int i, j;

  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {

      if(j != N - 1) {
        fprintf(fp, "%f ", matrix[i][j]);

      }

      else {
        fprintf(fp, "%f", matrix[i][j]);
        
      }

    }
    
    if(i != N - 1) {
      fprintf(fp, "\n");

    }

  }

  fclose(fp);

}

// This will read a matrix from a file.

void read_matrix() {
  FILE *fp;
  fp = fopen("matrix.txt", "r");

  int i, j;

  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++) {
      
      if (!fscanf(fp, "%lf", &matrix[i][j])) break;

    }

  } 

  fclose(fp);

}

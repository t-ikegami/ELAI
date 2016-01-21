#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "blas.hpp"

using namespace std;

typedef elai::vector< double > Vector;
typedef elai::matrix< double > Matrix;

int main()
{
  using namespace elai;

  int n = 5;
  int nnz = 17;
  int ind[] = { 0, 5, 8, 11, 14, 17 };
  int col[] =
    { 0, 1, 2, 3, 4
    , 0, 1, 2
    , 0, 1, 2
    , 0, 3, 4
    , 0, 3, 4
    };
  Matrix A( n, n, nnz, ind, col );
  Vector x( n );

  for ( int i = 0; i < n; ++i )
  {
    x( i ) = 0.5;
    A( i, i ) = -2.;
    for ( int j = 0; j < n; ++j ) if ( abs( i - j ) == 1 ) A( i, j ) = 1.;
  }
  cout << "A = " << A;

  Vector b = A * x;
/*
  cout << "A x = " << b;

  ifstream ifs( "./matrix.mm" );
  Matrix B( ifs );
  cout << "B = " << B;

  Matrix C = A + B;
  cout << "A + B = " << C;
*/

  int ord[ 5 ] = { 3, 2, 1, 4, 0 };
  cout << "Before: " << A;
  A.reorder( ord );
  cout << "After: " << A;
}

#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"

using namespace std;

typedef elai::vector< float > Vector;
typedef elai::matrix< float > Matrix;

int main()
{
  using namespace elai;

  int n = 5;
  int nnz = 13;
  int ind[] = { 0, 2, 5, 8, 11, nnz };
  int col[] =
    { 0, 1
    , 0, 1, 2
    , 1, 2, 3
    , 2, 3, 4
    , 3, 4
    };
  float c[] =
    { -2., 1.
    , 1., -2., 1.
    , 1., -2., 1.
    , 1., -2., 1.
    , 1., -2
    };
  Matrix A( n, n, nnz, ind, col, c );
  Vector x( n ), b( n );

  { // a * b
    x( 0 ) = 1.;
    x( 1 ) = 2.;
    x( 2 ) = 3.;
    x( 3 ) = 4.;
    x( 4 ) = 5.;
    b( 0 ) = .5;
    b( 1 ) = 1.;
    b( 2 ) = 1.5;
    b( 3 ) = 2.;
    b( 4 ) = 2.5;
    x = x - 2. * b;
    cout << x;
  }
}

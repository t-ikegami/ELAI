#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"

using namespace std;

typedef elai::vector< float > Vector;
typedef elai::matrix< float > Matrix;

int main()
{
  using namespace elai;

  int n = 6;
  int nnz = 16;
  int ind[] = { 0, 2, 5, 8, 10, 13, nnz };
  int col[] =
    { 0, 3
    , 0, 1, 4
    , 0, 2, 5
    , 0, 3
    , 1, 3, 4
    , 2, 3, 5
    };
  float c[] =
    { 2., -1.
    , 1.5, .1, -1.
    , -1.5, .1, -1.
    , -1., 2.
    , 1., 1.5, .1
    , -1., -1.5, .1
    };
  Matrix A( n, n, nnz, ind, col, c );
  Vector x( n ), b( n );

  for ( int i = 0; i < n; ++i )
  {
    x( i ) = 0.;
    b( i ) = 1.;
  }

  A.normalize();
  b.scale( A.scaleRow() );
  cout << A;
  cout << b;
  A.unnormalize();
  b.unscale( A.scaleRow() );
  cout << A;
  cout << b;
}

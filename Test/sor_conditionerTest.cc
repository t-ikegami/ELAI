#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "sor_conditioner.hpp"

using namespace std;

typedef elai::vector< float > Vector;
typedef elai::matrix< float > Matrix;
typedef elai::sor_conditioner< float > SorP;

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
    { 2., -1.
    , -1., 2., -1.
    , -1., 2., -1.
    , -1., 2., -1.
    , -1., 2
    };
  Matrix A( n, n, nnz, ind, col, c );
  Vector x( n ), b( n );
  SorP prec( A, 1.07 );
  // elai::preconditioner< float > *prec = new SorP( A );

  for ( int i = 0; i < n; ++i )
  {
    x( i ) = 0.;
    b( i ) = 1.;
  }
  x = b;

  cout << A;
  prec.forward( x );
  prec.backward( x );
  cout << x;
}

#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "ic.hpp"
#include "cg.hpp"

using namespace std;

typedef elai::vector< float > Vector;
typedef elai::matrix< float > Matrix;
typedef elai::ic< float > IC;
typedef elai::cg< float > CG;

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
  IC prec( A );
  CG solver( A, b, &prec );

  prec.factor();
  for ( int i = 0; i < n; ++i )
  {
    x( i ) = 0.;
    b( i ) = 1.;
  }

  bool flg = solver.solve( x );
  if ( flg )
  {
    cout << "Solved.." << endl;
    cout << x;
  }
  else cout << "Diverged!!" << endl;
}

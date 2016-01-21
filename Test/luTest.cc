#include <iostream>
#include "mpi.h"
#include "vector.hpp"
#include "matrix.hpp"
#include "lu.hpp"
#include "mumps.hpp"

using namespace std;

typedef elai::vector< double > Vector;
typedef elai::matrix< double > Matrix;
typedef elai::mumps< double > LU;

int main( int argc, char **argv )
{
  using namespace elai;
  int rank;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

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
  double c[] =
    { 2., -1.
    , -1., 2., -1.
    , -1., 2., -1.
    , -1., 2., -1.
    , -1., 2
    };
  Matrix A( n, n, nnz, ind, col, c );
  Vector x( n ), b( n );
  LU *lu;
  //LU lu( A, MPI_COMM_WORLD );

  lu = new LU( A, MPI_COMM_WORLD );
  for ( int i = 0; i < n; ++i )
  {
    x( i ) = 0.;
    b( i ) = 1.;
  }

  lu->factor();
  lu->solve( b, x );
  if ( rank == 0 ) cout << x;

  // Exact solution x =
  //  2.5 4 4.5 4 2.5

  delete lu;

  /*
  // Replay
  lu = new LU( A, MPI_COMM_WORLD );
  lu->factor();
  lu->solve( b, x );
  if ( rank == 0 ) cout << x;
  delete lu;
  */

  MPI_Finalize();
}

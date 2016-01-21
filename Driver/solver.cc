#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "stdlib.h"
#include "elai.h"

using namespace std;

typedef double Scalar;
typedef elai::time_monitor< double > Timer;
typedef elai::generator< Scalar > Generator;
typedef Generator::Space Space;
typedef Generator::Family Family;
typedef elai::subjugator< Generator::Element, Generator::Neighbour > Subjugator;
typedef elai::sync Sync;
typedef elai::coherence Coherence;
typedef elai::vector< Scalar > Vector;
typedef elai::matrix< Scalar > Matrix;
typedef elai::ksp< Scalar > KSP;
typedef elai::bicgstab< Scalar > BCGSTAB;
typedef elai::bicgsafe< Scalar > BCGSAFE;
typedef elai::gmres< Scalar > GMRES;
typedef elai::ilu< Scalar > ILU;
#ifdef ELAI_USE_MUMPS
typedef elai::mumps< Scalar > LU;
#endif
typedef elai::linear_function< Generator::Element, Generator::Neighbour, Scalar > Function;
typedef elai::linear_operator< Generator::Element, Generator::Neighbour, Scalar > Operator;

//
// PARAMETERS
//
enum ksp_method
  { ELAI_BCGS
  , ELAI_BCGSA
  , ELAI_GMRES
#ifdef ELAI_USE_MUMPS
  , ELAI_LU
#endif
  };
ksp_method method;
Scalar cthres, fthres, sthres;
int flevel, imax;
int mysize, myrank;
bool scaled, preconditioned;

bool solve( Matrix& A, Vector& x, Vector& b, Coherence *coherent = NULL )
{
  bool flg = false;
  double fact_elapsed;
  Timer t0;
  KSP *solver = NULL;
  ILU *prec = NULL;
  Scalar ratio = static_cast< Scalar >( 1. );
  Scalar rnorm, cnorm = rnorm = static_cast< Scalar >( 1. );

  fact_elapsed = t0();

  if ( scaled )
  {
    A.normalize( sthres, coherent );
    b.scale( A.scaleRow() );
    x.unscale( A.scaleCol() );
    rnorm = A.scaleRowNorm();
    cnorm = A.scaleColNorm();
    ratio = A.scaleRatio();
  }

  if ( preconditioned ) prec = new ILU( A, flevel, fthres, false );

  if ( method == ELAI_BCGS ) solver = new BCGSTAB( A, b, prec, coherent );
  else if ( method == ELAI_BCGSA ) solver = new BCGSAFE( A, b, prec, coherent );
  else if ( method == ELAI_GMRES ) solver = new GMRES( A, b, prec, coherent );

  solver->iter_max( imax );
  solver->rel_thres( cthres );

  if ( preconditioned ) prec->factor( fthres );
  flg = solver->solve( x );

  if ( scaled )
  {
    x.scale( A.scaleCol() );
    b.unscale( A.scaleRow() );
    A.unnormalize();
  }

  if ( myrank == 0 ) {
    cerr << "FACT: " << fact_elapsed << "-sec." << endl;
    cerr << "PREC: " << solver->prec_elapsed() << "-sec." << endl;
    cerr << "ITER: " << solver->elapsed() - solver->prec_elapsed() << "-sec." << endl;
    cerr << "SCAL: " << ratio << " ( " << rnorm << ", " << cnorm << " ) " << endl;
  }

  return flg;
}

void solve_dist
  ( Operator& A
  , Function& U
  , Function& V
  , const vector< int >& ranks
  )
{
  Subjugator loc( A.dom(), A.topo(), ranks );
  Space base( loc( myrank ) );
  Family topo( loc( myrank, base ) );
  std::cout << base.size() << std::endl;
  Operator a( base, topo );
  Function u( base ), v( base );
  Coherence coherent( u, MPI_COMM_WORLD );
  Sync sync( U, MPI_COMM_WORLD );

  a.reflectIn( A, loc );
  u.reflectIn( U, loc );
  v.reflectIn( V, loc );
  solve( a.action(), u.ran(), v.ran(), &coherent );
  U.clear( 0e0 );
  U.reflect( u, loc );
  sync();
}

void run_dist
  ( Matrix& A
  , Vector& x
  , Vector& b
  )
{
  Generator gen( A );
  const Space& base = gen.space();
  const Family& topo = gen.family();
  Operator a( base, base, topo, A );
  Function u( base, x ), v( base, b );
  vector< int > ranks;

  for ( int i = 0; i < mysize; ++i ) ranks.push_back( i );
  solve_dist( a, u, v, ranks );

  if ( myrank == 0 )
  {
    Matrix& A = a.action();
    Vector& x = u.ran();
    Vector& b = v.ran();
    Vector res = A * x - b;
    Scalar r = res * res, r0;

    r = sqrt( r );
    cout << "||A x - b|| = " << setprecision( 15 ) << r << endl;
    r0 = x * x; r0 = sqrt( r0 );
    r /= r0;
    cout << "||A x - b||/||x|| = " << setprecision( 15 ) << r << endl;
  }
}

void run
  ( Matrix& A
  , Vector& x
  , Vector& b
  )
{
  solve( A, x, b );

  Vector v = A * x - b;
  Scalar r = v * v, r0;

  r = sqrt( r );
  cout << "||A x - b|| = " << setprecision( 15 ) << r << endl;
  r0 = x * x; r0 = sqrt( r0 );
  r /= r0;
  cout << "||A x - b||/||x|| = " << setprecision( 15 ) << r << endl;
}

#ifdef ELAI_USE_MUMPS
void direct
  ( Matrix& A
  , Vector& x
  , Vector& b
  )
{
  LU lu( A, MPI_COMM_WORLD );

  lu.factor();
  lu.solve( b, x );
  if ( myrank == 0 )
  {
    cout << "USE " << mysize << "-PROC." << endl;
    Vector res( x );
    double r, r0;
    //ofstream xfile( "x.mtx" );

    res = A * x - b;
    r = res * res; r = sqrt( r );
    cout << "||A x - b|| = " << setprecision( 15 ) << r << endl;
    r0 = x * x; r0 = sqrt( r0 );
    r /= r0;
    cout << "||A x - b||/||x|| = " << setprecision( 15 ) << r << endl;
    //xfile << setprecision( 15 ) << x;
  }
}
#endif

void demo()
{
  typedef Generator::Element Node;
  typedef Generator::Neighbour Link;

  int nnd = 16;
  int color = 0;

  Space mesh;
  Family topo;

  for ( int i = 0; i < nnd; ++i )
  {
    Node n( i, color );
    mesh.join( n );

    Link lnk( n );
    if ( 0 <= i - 4 ) lnk.join( Node( i - 4, color ) );
    if ( 0 <= i - 1 ) lnk.join( Node( i - 1, color ) );
    if ( i + 1 < nnd ) lnk.join( Node( i + 1, color ) );
    if ( i + 4 < nnd ) lnk.join( Node( i + 4, color ) );
    topo.join( lnk );
  }

  Operator problem( mesh, topo );
  Function solution( mesh ), inhomogeneous( mesh );

  for ( int i = 0; i < nnd; ++i )
  {
    problem( Node( i, color ) ) = -4;
    if ( 0 <= i - 4 ) problem( Node( i - 4, color ), Node( i, color ) ) = 1;
    if ( 0 <= i - 1 ) problem( Node( i - 1, color ), Node( i, color ) ) = 1;
    if ( i + 1 < nnd ) problem( Node( i, color ), Node( i + 1, color ) ) = 1;
    if ( i + 4 < nnd ) problem( Node( i, color ), Node( i + 4, color ) ) = 1;
    inhomogeneous( Node( i, color ) ) = 1;
  }

  {
    Matrix& A = problem.action();
    Vector& x = solution.ran();
    Vector& b = inhomogeneous.ran();

    if ( myrank == 0 )
    {
      cout << "DEMO SETTINGS.." << endl << endl;
      cout << "COEFFS=" << endl << A << endl;
      cout << "RHS   =" << endl << b << endl;
      cout << "USING MUMPS.." << endl;
    }
#ifdef ELAI_USE_MUMPS
    direct( A, x, b );
#endif
  }
}

int main( int argc, char *argv[] )
{
  using namespace elai;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mysize );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  const string CTHRES( getenv( "CTHRES" ) );
  const string FLEVEL( getenv( "FLEVEL" ) );
  const string FTHRES( getenv( "FTHRES" ) );
  const string STHRES( getenv( "STHRES" ) );
  const string KSP( getenv( "KSP" ) );

  istringstream( CTHRES ) >> cthres;
  istringstream( FTHRES ) >> fthres;
  istringstream( STHRES ) >> sthres;
  istringstream( FLEVEL ) >> flevel;
  if ( myrank == 0 )
  {
    cout << setprecision( 15 );
    cout << KSP << " " << cthres << " " << fthres << " " << sthres << " " << flevel << std::endl;
  }

  if ( argc < 2 )
  {
    demo();
    MPI_Finalize();

    return 0;
  }

  ifstream mfile( argv[ 1 ] );
  ifstream bfile( argv[ 2 ] );
  Matrix A( mfile );
  Vector b( bfile );
  Vector x( b.m() );

  imax = A.m() / 5;
  scaled = false;
  preconditioned = false;
  cerr << setprecision( 15 );
#ifdef ELAI_USE_MUMPS
  method = ELAI_LU;
  if ( !KSP.compare( "LU" ) ) method = ELAI_LU; else
#else
  method = ELAI_BCGS;
#endif
  if ( !KSP.compare( "BCGS" ) ) method = ELAI_BCGS;
  else if ( !KSP.compare( "SBCGS" ) ) { scaled = true; method = ELAI_BCGS; }
  else if ( !KSP.compare( "PBCGS" ) ) { preconditioned = true; method = ELAI_BCGS; }
  else if ( !KSP.compare( "SPBCGS" ) ) { scaled = true; preconditioned = true; method = ELAI_BCGS; }
  else if ( !KSP.compare( "BCGSA" ) ) method = ELAI_BCGSA;
  else if ( !KSP.compare( "SBCGSA" ) ) { scaled = true; method = ELAI_BCGSA; }
  else if ( !KSP.compare( "PBCGSA" ) ) { preconditioned = true; method = ELAI_BCGSA; }
  else if ( !KSP.compare( "SPBCGSA" ) ) { scaled = true; preconditioned = true; method = ELAI_BCGSA; }
  else if ( !KSP.compare( "GMRES" ) ) method = ELAI_GMRES;
  else if ( !KSP.compare( "SGMRES" ) ) { scaled = true; method = ELAI_GMRES; }
  else if ( !KSP.compare( "PGMRES" ) ) { preconditioned = true; method = ELAI_GMRES; }
  else if ( !KSP.compare( "SPGMRES" ) ) { scaled = true; preconditioned = true; method = ELAI_GMRES; }

#ifdef ELAI_USE_MUMPS
  if ( method == ELAI_LU ) direct( A, x, b ); else
#endif
  if ( 1 < mysize ) run_dist( A, x, b );
  else run( A, x, b );

  MPI_Finalize();
}

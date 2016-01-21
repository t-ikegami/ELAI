#include <iostream>
#include <fstream>
#include "blas.hpp"

using namespace std;

typedef elai::vector< double > Vector;

int main()
{
  using namespace elai;

  int n = 10;
  Vector u( n );
  for ( int i = 0; i < n; ++i ) u( i ) = 0.5 * i;

  cout << "u = " << u;

  ifstream ifs( "./vector.mm" );
  Vector v( ifs );
  cout << "v = " << v;

  double prod = u * v;
  cout << "u * v = " << prod << endl;

  Vector w = 0.5 * ( u + v );
  cout << "0.5 ( u + v ) = " << w;

  w = w - v;
  cout << "u - v = " << w;

  Vector dic( 4 );
  int ord[ 4 ] = { 3, 2, 1, 0 };
  for ( int i = 0; i < 4; ++i ) dic( i ) = i + 1;
  cout << "Before:" << dic << endl;
  dic.reorder( ord );
  cout << "After:" << dic << endl;
}

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "def.hpp"
#include "portal.hpp"
#include "space.hpp"
#include "family.hpp"
#include "linear_operator.hpp"
#include "entire_operator.hpp"
using namespace std;
int myrank;

typedef int NID;
const NID WHOLE = -1;

typedef int VAR;
const VAR ANY  =15;
const VAR PSI  = 1;
const VAR ELEC = 2;
const VAR HOLE = 4;
const VAR NUM  = 8;

class Element
{
  const NID i_;
  const VAR t_;
  const int r_;

public:
  typedef long long id_type;

  Element( id_type id )
    : i_( ( id & 0xfffffffcl ) >> 2 )
    , t_( id & 0x3l )
    , r_( ( id & 0xffffffff00000000l ) >> 32 )
  {
    cout << myrank << ":" << i_ << " " << t_ << " " << r_ << std::endl;
  }
  Element( NID i, VAR t, int r ) : i_( i ), t_( t ), r_( r ) {}
  Element( const Element& x ) : i_( x.i_ ), t_( x.t_ ), r_( x.r_ ) {}
  Element( const Element& x, int c ) : i_( x.i_ ), t_( x.t_ ), r_( c ) {}
  ~Element() {}

  NID id() const { return i_; }
  VAR var() const { return t_; }
  id_type operator()() const
  {
    id_type r = static_cast< id_type >( r_ ) & 0xffffffffl;
    id_type i = static_cast< id_type >( i_ ) & 0xcfffffffl;
    id_type t = static_cast< id_type >( t_ ) & 0x3l;
    //cout << myrank << ":" << r << "-" << i << "-" << t << " = ";
    //cout << ( ( r << 32 ) | ( i << 2 ) | t ) << endl;

    return ( r << 32 ) | ( i << 2 ) | t;
  }
  int color() const { return r_; }

  bool operator<( const Element& x ) const
  {
    if ( r_ != x.r_ ) return r_ < x.r_;
    if ( i_ < x.i_ ) return true;
    return i_ == x.i_ && t_ < x.t_;
  }
  bool operator==( const Element& x ) const
  {
    return r_ == x.r_ && i_ == x.i_ && t_ == x.t_;
  }
};

class Neighbour
{
  typedef set< Element > Adjacency;

  const Element x_;
  Adjacency adj_;

public:
  typedef Adjacency::iterator iterator;
  typedef Adjacency::const_iterator const_iterator;

  Neighbour( NID i, VAR t, int r ) : x_( Element( i, t, r ) ) {}
  Neighbour( const Element& x ) : x_( x )
  {
    if ( x_.id() != WHOLE && x_.var() != ANY ) adj_.insert( x_ );
  }
  Neighbour( const Neighbour& u ) : x_( u.x_ ), adj_( u.adj_ ) {}
  ~Neighbour() {}

  Neighbour& operator=( const Neighbour& src )
  {
    if ( this == &src ) return *this;
    if ( x_ == src.x_ ) adj_ = src.adj_;
    return *this;
  }

  const Element& element() const { return x_; }
  int size() const { return adj_.size(); }

  Neighbour& join( const Element& e )
  {
    assert( e.id() != WHOLE && e.var() != ANY );
    adj_.insert( e ); return *this;
  }
  Neighbour& join( NID i, VAR t, int r )
  {
    assert( i != WHOLE && t != ANY );
    adj_.insert( Element( i, t, r ) ); return *this;
  }

  Neighbour& erase( const Element& x )
  {
    adj_.erase( x ); return *this;
  }
  Neighbour& erase( const_iterator& it )
  {
    adj_.erase( it ); return *this;
  }

  bool operator()( const Element& x ) const
  {
    if ( x_.id() == WHOLE && x_.var() == ANY ) return true;
    else if ( x_.id() == WHOLE && x_.var() == x.var() ) return true;
    else if ( x_.id() == x.id() && x_.var() == ANY ) return true;
    for ( const_iterator it = adj_.begin(); it != adj_.end(); ++it ) if ( *it == x ) return true;
    return false;
  }

  iterator begin() { return adj_.begin(); }
  const_iterator begin() const { return adj_.begin(); }
  iterator end() { return adj_.end(); }
  const_iterator end() const { return adj_.end(); }
};

typedef elai::portal Portal;
typedef elai::space< Element > Space;
typedef struct elai::space< Element >::const_point Point;
typedef elai::family< Element, Neighbour > Family;
typedef elai::linear_operator< Element, Neighbour, float > Operator;
typedef elai::entire_operator< Element, Neighbour, float > OperatorDist;

int main( int argc, char **argv )
{
  int sz, rk;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &sz );
  MPI_Comm_rank( MPI_COMM_WORLD, &rk );

  if ( sz != 4 ) exit( 1 );

  if ( rk == 0 )
  {
    int rid = rk;

    Space s;
    s.join( Element( 0, PSI, rid ) );
    s.join( Element( 1, PSI, rid ) );
    s.join( Element( 2, PSI, rid ), Element( 1, PSI, 1 ) );
    s.join( Element( 3, PSI, rid ) );
    s.join( Element( 4, PSI, rid ) );
    s.join( Element( 5, PSI, rid ), Element( 4, PSI, 1 ) );
    s.join( Element( 6, PSI, rid ), Element( 3, PSI, 2 ) );
    s.join( Element( 7, PSI, rid ), Element( 4, PSI, 2 ) );
    s.join( Element( 8, PSI, rid ), Element( 4, PSI, 3 ) );

    Neighbour l0( Element( 0, PSI, rid ) );
    l0.join( Element( 1, PSI, rid ) );
    l0.join( Element( 3, PSI, rid ) );

    Neighbour l1( Element( 1, PSI, rid ) );
    l1.join( Element( 0, PSI, rid ) );
    l1.join( Element( 2, PSI, rid ) );
    l1.join( Element( 4, PSI, rid ) );

    Neighbour l2( Element( 2, PSI, rid ) );
    l2.join( Element( 1, PSI, rid ) );
    l2.join( Element( 5, PSI, rid ) );

    Neighbour l3( Element( 3, PSI, rid ) );
    l3.join( Element( 0, PSI, rid ) );
    l3.join( Element( 4, PSI, rid ) );
    l3.join( Element( 6, PSI, rid ) );

    Neighbour l4( Element( 4, PSI, rid ) );
    l4.join( Element( 1, PSI, rid ) );
    l4.join( Element( 3, PSI, rid ) );
    l4.join( Element( 5, PSI, rid ) );
    l4.join( Element( 7, PSI, rid ) );

    Neighbour l5( Element( 5, PSI, rid ) );
    l5.join( Element( 2, PSI, rid ) );
    l5.join( Element( 4, PSI, rid ) );
    l5.join( Element( 8, PSI, rid ) );

    Neighbour l6( Element( 6, PSI, rid ) );
    l6.join( Element( 3, PSI, rid ) );
    l6.join( Element( 7, PSI, rid ) );

    Neighbour l7( Element( 7, PSI, rid ) );
    l7.join( Element( 4, PSI, rid ) );
    l7.join( Element( 6, PSI, rid ) );
    l7.join( Element( 8, PSI, rid ) );

    Neighbour l8( Element( 8, PSI, rid ) );
    l8.join( Element( 5, PSI, rid ) );
    l8.join( Element( 7, PSI, rid ) );

    Family t( s );
    t.join( l0 ).join( l1 ).join( l2 ).join( l3 ).join( l4 ).join( l5 ).join( l6 ).join( l7 ).join( l8 );

    OperatorDist dist;
    Operator *a[ sz ];
    a[ rid ] = new Operator( s, t );
    dist.join( *a[ 0 ] );
    for ( int i = 1; i < sz; ++i )
    {
      Portal port( i, MPI_COMM_WORLD );

      cout << "Recv from " << i << endl;
      a[ i ] = new Operator( port );
      cout << "Merge with " << i << endl;
      dist.join( *a[ i ] );
    }

    Operator A( dist( true ) );
    A.action().clear( 1.01e0 );
    ofstream of( "entire.mtx" );
    of << setprecision( 15 ) << A.action();

    for ( int i = 0; i < sz; ++i ) delete a[ i ];
  }
  else if ( rk == 1 )
  {
    int rid = rk;

    Space s;
    s.join( Element( 0, PSI, rid ), Element( 1, PSI, 0 ) );
    s.join( Element( 1, PSI, rid ) );
    s.join( Element( 2, PSI, rid ) );
    s.join( Element( 3, PSI, rid ), Element( 4, PSI, 0 ) );
    s.join( Element( 4, PSI, rid ) );
    s.join( Element( 5, PSI, rid ) );
    s.join( Element( 6, PSI, rid ), Element( 4, PSI, 2 ) );
    s.join( Element( 7, PSI, rid ), Element( 4, PSI, 3 ) );
    s.join( Element( 8, PSI, rid ), Element( 5, PSI, 3 ) );

    Neighbour l0( Element( 0, PSI, rid ) );
    l0.join( Element( 1, PSI, rid ) );
    l0.join( Element( 3, PSI, rid ) );

    Neighbour l1( Element( 1, PSI, rid ) );
    l1.join( Element( 0, PSI, rid ) );
    l1.join( Element( 2, PSI, rid ) );
    l1.join( Element( 4, PSI, rid ) );

    Neighbour l2( Element( 2, PSI, rid ) );
    l2.join( Element( 1, PSI, rid ) );
    l2.join( Element( 5, PSI, rid ) );

    Neighbour l3( Element( 3, PSI, rid ) );
    l3.join( Element( 0, PSI, rid ) );
    l3.join( Element( 4, PSI, rid ) );
    l3.join( Element( 6, PSI, rid ) );

    Neighbour l4( Element( 4, PSI, rid ) );
    l4.join( Element( 1, PSI, rid ) );
    l4.join( Element( 3, PSI, rid ) );
    l4.join( Element( 5, PSI, rid ) );
    l4.join( Element( 7, PSI, rid ) );

    Neighbour l5( Element( 5, PSI, rid ) );
    l5.join( Element( 2, PSI, rid ) );
    l5.join( Element( 4, PSI, rid ) );
    l5.join( Element( 8, PSI, rid ) );

    Neighbour l6( Element( 6, PSI, rid ) );
    l6.join( Element( 3, PSI, rid ) );
    l6.join( Element( 7, PSI, rid ) );

    Neighbour l7( Element( 7, PSI, rid ) );
    l7.join( Element( 4, PSI, rid ) );
    l7.join( Element( 6, PSI, rid ) );
    l7.join( Element( 8, PSI, rid ) );

    Neighbour l8( Element( 8, PSI, rid ) );
    l8.join( Element( 5, PSI, rid ) );
    l8.join( Element( 7, PSI, rid ) );

    Family t( s );
    t.join( l0 ).join( l1 ).join( l2 ).join( l3 ).join( l4 ).join( l5 ).join( l6 ).join( l7 ).join( l8 );

    Operator A( s, t );
    Portal port( 0, MPI_COMM_WORLD );
    port( A );
    cout << rid << " has sent to root" << endl;
  }
  else if ( rk == 2 )
  {
    int rid = rk;

    Space s;
    s.join( Element( 0, PSI, rid ), Element( 3, PSI, 0 ) );
    s.join( Element( 1, PSI, rid ), Element( 4, PSI, 0 ) );
    s.join( Element( 2, PSI, rid ), Element( 4, PSI, 1 ) );
    s.join( Element( 3, PSI, rid ) );
    s.join( Element( 4, PSI, rid ) );
    s.join( Element( 5, PSI, rid ), Element( 4, PSI, 3 ) );
    s.join( Element( 6, PSI, rid ) );
    s.join( Element( 7, PSI, rid ) );
    s.join( Element( 8, PSI, rid ), Element( 7, PSI, 3 ) );

    Neighbour l0( Element( 0, PSI, rid ) );
    l0.join( Element( 1, PSI, rid ) );
    l0.join( Element( 3, PSI, rid ) );

    Neighbour l1( Element( 1, PSI, rid ) );
    l1.join( Element( 0, PSI, rid ) );
    l1.join( Element( 2, PSI, rid ) );
    l1.join( Element( 4, PSI, rid ) );

    Neighbour l2( Element( 2, PSI, rid ) );
    l2.join( Element( 1, PSI, rid ) );
    l2.join( Element( 5, PSI, rid ) );

    Neighbour l3( Element( 3, PSI, rid ) );
    l3.join( Element( 0, PSI, rid ) );
    l3.join( Element( 4, PSI, rid ) );
    l3.join( Element( 6, PSI, rid ) );

    Neighbour l4( Element( 4, PSI, rid ) );
    l4.join( Element( 1, PSI, rid ) );
    l4.join( Element( 3, PSI, rid ) );
    l4.join( Element( 5, PSI, rid ) );
    l4.join( Element( 7, PSI, rid ) );

    Neighbour l5( Element( 5, PSI, rid ) );
    l5.join( Element( 2, PSI, rid ) );
    l5.join( Element( 4, PSI, rid ) );
    l5.join( Element( 8, PSI, rid ) );

    Neighbour l6( Element( 6, PSI, rid ) );
    l6.join( Element( 3, PSI, rid ) );
    l6.join( Element( 7, PSI, rid ) );

    Neighbour l7( Element( 7, PSI, rid ) );
    l7.join( Element( 4, PSI, rid ) );
    l7.join( Element( 6, PSI, rid ) );
    l7.join( Element( 8, PSI, rid ) );

    Neighbour l8( Element( 8, PSI, rid ) );
    l8.join( Element( 5, PSI, rid ) );
    l8.join( Element( 7, PSI, rid ) );

    Family t( s );
    t.join( l0 ).join( l1 ).join( l2 ).join( l3 ).join( l4 ).join( l5 ).join( l6 ).join( l7 ).join( l8 );

    Operator A( s, t );
    Portal port( 0, MPI_COMM_WORLD );
    port( A );
    cout << rid << " has sent to root" << endl;
  }
  else if ( rk == 3 )
  {
    int rid = rk;

    Space s;
    s.join( Element( 0, PSI, rid ), Element( 4, PSI, 0 ) );
    s.join( Element( 1, PSI, rid ), Element( 4, PSI, 1 ) );
    s.join( Element( 2, PSI, rid ), Element( 5, PSI, 1 ) );
    s.join( Element( 3, PSI, rid ), Element( 4, PSI, 2 ) );
    s.join( Element( 4, PSI, rid ) );
    s.join( Element( 5, PSI, rid ) );
    s.join( Element( 6, PSI, rid ), Element( 7, PSI, 2 ) );
    s.join( Element( 7, PSI, rid ) );
    s.join( Element( 8, PSI, rid ) );

    Neighbour l0( Element( 0, PSI, rid ) );
    l0.join( Element( 1, PSI, rid ) );
    l0.join( Element( 3, PSI, rid ) );

    Neighbour l1( Element( 1, PSI, rid ) );
    l1.join( Element( 0, PSI, rid ) );
    l1.join( Element( 2, PSI, rid ) );
    l1.join( Element( 4, PSI, rid ) );

    Neighbour l2( Element( 2, PSI, rid ) );
    l2.join( Element( 1, PSI, rid ) );
    l2.join( Element( 5, PSI, rid ) );

    Neighbour l3( Element( 3, PSI, rid ) );
    l3.join( Element( 0, PSI, rid ) );
    l3.join( Element( 4, PSI, rid ) );
    l3.join( Element( 6, PSI, rid ) );

    Neighbour l4( Element( 4, PSI, rid ) );
    l4.join( Element( 1, PSI, rid ) );
    l4.join( Element( 3, PSI, rid ) );
    l4.join( Element( 5, PSI, rid ) );
    l4.join( Element( 7, PSI, rid ) );

    Neighbour l5( Element( 5, PSI, rid ) );
    l5.join( Element( 2, PSI, rid ) );
    l5.join( Element( 4, PSI, rid ) );
    l5.join( Element( 8, PSI, rid ) );

    Neighbour l6( Element( 6, PSI, rid ) );
    l6.join( Element( 3, PSI, rid ) );
    l6.join( Element( 7, PSI, rid ) );

    Neighbour l7( Element( 7, PSI, rid ) );
    l7.join( Element( 4, PSI, rid ) );
    l7.join( Element( 6, PSI, rid ) );
    l7.join( Element( 8, PSI, rid ) );

    Neighbour l8( Element( 8, PSI, rid ) );
    l8.join( Element( 5, PSI, rid ) );
    l8.join( Element( 7, PSI, rid ) );

    Family t( s );
    t.join( l0 ).join( l1 ).join( l2 ).join( l3 ).join( l4 ).join( l5 ).join( l6 ).join( l7 ).join( l8 );

    Operator A( s, t );
    Portal port( 0, MPI_COMM_WORLD );
    port( A );
    cout << rid << " has sent to root" << endl;
  }

  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
}

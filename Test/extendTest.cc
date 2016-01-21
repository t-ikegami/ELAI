#include <iostream>
#include "mpi.h"
#include "portal.hpp"
#include "space.hpp"
#include "family.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "linear_function.hpp"
#include "linear_operator.hpp"
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

typedef elai::space< Element > Space;
typedef struct elai::space< Element >::const_point Point;
typedef elai::family< Element, Neighbour > Family;
typedef elai::vector< float > Vector;
typedef elai::matrix< float > Matrix;
typedef elai::linear_function< Element, Neighbour, float > Function;
typedef elai::linear_operator< Element, Neighbour, float > Operator;
typedef elai::portal Portal;

int main( int argc, char **argv )
{
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank ); 

  if ( myrank == 0 )
  {
    Space s1;
    int pair = 1, n = 4;
    Portal port( pair, MPI_COMM_WORLD );

    s1.join( Element( 0, PSI, myrank ) );
    s1.join( Element( 1, PSI, myrank ) );
    s1.join( Element( 2, PSI, myrank ), Element( 0, PSI, pair ) );
    s1.join( Element( 3, PSI, myrank ), Element( 1, PSI, pair ) );

    Neighbour u0( Element( 0, PSI, myrank ) );
    u0.join( Element( 1, PSI, myrank ) );

    Neighbour u1( Element( 1, PSI, myrank ) );
    u1.join( Element( 0, PSI, myrank ) );
    u1.join( Element( 2, PSI, myrank ) );

    Neighbour u2( Element( 2, PSI, myrank ) );
    u2.join( Element( 1, PSI, myrank ) );
    u2.join( Element( 3, PSI, myrank ) );

    Neighbour u3( Element( 3, PSI, myrank ) );
    u3.join( Element( 2, PSI, myrank ) );

    Family T1( s1 );
    T1.join( u0 ).join( u1 ).join( u2 ).join( u3 );

    Operator A1( s1, T1 );
    {
      Matrix& a = A1.action();
      for ( int i = 0; i < n; ++i )
      {
        if ( i != 0 ) a( i, i-1 ) = 1.;
        a( i, i ) = -2.;
        if ( i != n-1 ) a( i, i+1 ) = 1.;
      }
    }
    Operator A2( port );
    Operator A( A1.extend( A2 ) );

    cout << A.action();
    MPI_Barrier( MPI_COMM_WORLD );
  }
  else if ( myrank == 1 )
  {
    int pair = 0, n = 4;
    Portal port( pair, MPI_COMM_WORLD );

    Space s1;
    s1.join( Element( 0, PSI, myrank ), Element( 2, PSI, pair ) );
    s1.join( Element( 1, PSI, myrank ), Element( 3, PSI, pair ) );
    s1.join( Element( 2, PSI, myrank ) );
    s1.join( Element( 3, PSI, myrank ) );

    Neighbour u0( Element( 0, PSI, myrank ) );
    u0.join( Element( 1, PSI, myrank ) );

    Neighbour u1( Element( 1, PSI, myrank ) );
    u1.join( Element( 0, PSI, myrank ) );
    u1.join( Element( 2, PSI, myrank ) );

    Neighbour u2( Element( 2, PSI, myrank ) );
    u2.join( Element( 1, PSI, myrank ) );
    u2.join( Element( 3, PSI, myrank ) );

    Neighbour u3( Element( 3, PSI, myrank ) );
    u3.join( Element( 2, PSI, myrank ) );

    Family T1( s1 );
    T1.join( u0 ).join( u1 ).join( u2 ).join( u3 );

    Operator A( s1, T1 );
    {
      Matrix& a = A.action();
      for ( int i = 0; i < n; ++i )
      {
        if ( i != 0 ) a( i, i-1 ) = 1.;
        a( i, i ) = -2.;
        if ( i != n-1 ) a( i, i+1 ) = 1.;
      }
    }
    port( A );

    MPI_Barrier( MPI_COMM_WORLD );
  }
  MPI_Finalize();
}

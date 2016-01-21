#include <iostream>
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"
#include "coherence.hpp"
#include "sync.hpp"
#include "portal.hpp"
#include "subjugator.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "linear_function.hpp"
#include "linear_operator.hpp"
using namespace std;
int myrank, mysize;

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
    //cout << myrank << ":" << i_ << " " << t_ << " " << r_ << std::endl;
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

typedef elai::coherence Coherence;
typedef elai::sync Sync;
typedef elai::space< Element > Space;
typedef elai::family< Element, Neighbour > Family;
typedef elai::subjugator< Element, Neighbour > Subjugator;
typedef elai::vector< float > Vector;
typedef elai::matrix< float > Matrix;
typedef elai::linear_function< Element, Neighbour, float > Function;
typedef elai::linear_operator< Element, Neighbour, float > Operator;

void ddc0( MPI_Comm intra, MPI_Comm inter )
{
  int rid = 0, pair = 1;
  std::vector< int > palette;
  palette.push_back( 0 );
  palette.push_back( 2 );

  Space all;
  all.join( Element( 0, PSI, rid ) );
  all.join( Element( 1, PSI, rid ) );
  all.join( Element( 2, PSI, rid ), Element( 1, PSI, pair ) );
  all.join( Element( 3, PSI, rid ) );
  all.join( Element( 4, PSI, rid ) );
  all.join( Element( 5, PSI, rid ), Element( 4, PSI, pair ) );

  Neighbour u0( Element( 0, PSI, rid ) );
  u0.join( Element( 1, PSI, rid ) );
  u0.join( Element( 3, PSI, rid ) );

  Neighbour u1( Element( 1, PSI, rid ) );
  u1.join( Element( 0, PSI, rid ) );
  u1.join( Element( 2, PSI, rid ) );
  u1.join( Element( 4, PSI, rid ) );

  Neighbour u2( Element( 2, PSI, rid ) );
  u2.join( Element( 1, PSI, rid ) );
  u2.join( Element( 5, PSI, rid ) );

  Neighbour u3( Element( 3, PSI, rid ) );
  u3.join( Element( 0, PSI, rid ) );
  u3.join( Element( 4, PSI, rid ) );

  Neighbour u4( Element( 4, PSI, rid ) );
  u4.join( Element( 1, PSI, rid ) );
  u4.join( Element( 3, PSI, rid ) );
  u4.join( Element( 5, PSI, rid ) );

  Neighbour u5( Element( 5, PSI, rid ) );
  u5.join( Element( 2, PSI, rid ) );
  u5.join( Element( 4, PSI, rid ) );

  Family T1( all );
  T1.join( u0 ).join( u1 ).join( u2 ).join( u3 ).join( u4 ).join( u5 );

  Subjugator subj( all, T1, palette, intra, inter );

  const Space sub = subj( myrank );
  const Family t = subj( myrank, sub );
  Operator A( sub, t );
  Function f( sub );
  Coherence coherent( f, MPI_COMM_WORLD );

  {
    Vector& x = f.ran();
    for ( int i = 0; i < x.m(); ++i ) x( i ) = myrank;
  }
  coherent( f.ran().val() );
  cout << "sub" << myrank << ":" << f.ran();

  Function F( all );
  Sync F_sync( F, intra );

  F.reflect( f, subj );
  F_sync();
}

void ddc1( MPI_Comm intra, MPI_Comm inter )
{
  int rid = 1, pair = 0;
  std::vector< int > palette;
  palette.push_back( 1 );
  palette.push_back( 3 );

  Space all;
  all.join( Element( 0, PSI, rid ), Element( 1, PSI, pair ) );
  all.join( Element( 1, PSI, rid ) );
  all.join( Element( 2, PSI, rid ) );
  all.join( Element( 3, PSI, rid ), Element( 4, PSI, pair ) );
  all.join( Element( 4, PSI, rid ) );
  all.join( Element( 5, PSI, rid ) );

  Neighbour u0( Element( 0, PSI, rid ) );
  u0.join( Element( 1, PSI, rid ) );
  u0.join( Element( 3, PSI, rid ) );

  Neighbour u1( Element( 1, PSI, rid ) );
  u1.join( Element( 0, PSI, rid ) );
  u1.join( Element( 2, PSI, rid ) );
  u1.join( Element( 4, PSI, rid ) );

  Neighbour u2( Element( 2, PSI, rid ) );
  u2.join( Element( 1, PSI, rid ) );
  u2.join( Element( 5, PSI, rid ) );

  Neighbour u3( Element( 3, PSI, rid ) );
  u3.join( Element( 0, PSI, rid ) );
  u3.join( Element( 4, PSI, rid ) );

  Neighbour u4( Element( 4, PSI, rid ) );
  u4.join( Element( 1, PSI, rid ) );
  u4.join( Element( 3, PSI, rid ) );
  u4.join( Element( 5, PSI, rid ) );

  Neighbour u5( Element( 5, PSI, rid ) );
  u5.join( Element( 2, PSI, rid ) );
  u5.join( Element( 4, PSI, rid ) );

  Family T1( all );
  T1.join( u0 ).join( u1 ).join( u2 ).join( u3 ).join( u4 ).join( u5 );

  Subjugator subj( all, T1, palette, intra, inter );

  const Space sub = subj( myrank );
  const Family t = subj( myrank, sub );
  Operator A( sub, t );
  Function f( sub );
  Coherence coherent( f, MPI_COMM_WORLD );

  {
    Vector& x = f.ran();
    for ( int i = 0; i < x.m(); ++i ) x( i ) = myrank;
  }
  coherent( f.ran().val() );
  cout << "sub" << myrank << ":" << f.ran();

  Function F( all );
  Sync F_sync( F, intra );

  F.reflect( f, subj );
  F_sync();
}

int main( int argc, char **argv )
{
  // Leader: 0, 1 => inter
  // RID=0 : 0, 2 => intra
  // RID=1 : 1, 3 => intra
  // WHOLE : 0, 1, 2, 3 => MPI_COMM_WORLD
  int rank_leader[] = { 0, 1 };
  int is_leader, rid;
  MPI_Comm intra, inter;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &mysize );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  rid = myrank % 2;
  is_leader = MPI_UNDEFINED;
  for ( int i = 0; i < 2; ++i )
    if ( rank_leader[ i ] == myrank ) is_leader = 1;
  MPI_Comm_split( MPI_COMM_WORLD, is_leader, rid, &inter );
  MPI_Comm_split( MPI_COMM_WORLD, rid, myrank, &intra );

  if ( mysize == 4 )
  {
    if ( rid == 0 ) ddc0( intra, inter );
    else ddc1( intra, inter );
  }

  MPI_Finalize();
}

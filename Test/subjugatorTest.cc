#include <iostream>
#include <vector>
#include "space.hpp"
#include "family.hpp"
#include "subjugator.hpp"
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
    cout << myrank << ":" << r << "-" << i << "-" << t << " = ";
    cout << ( ( r << 32 ) | ( i << 2 ) | t ) << endl;

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

typedef elai::space< Element > Space;
typedef struct elai::space< Element >::const_point Point;
typedef elai::family< Element, Neighbour > Family;
typedef elai::subjugator< Element, Neighbour > Subjugator;

int main( int argc, char **argv )
{
  Space all;
  all.join( Element( 0, PSI, 0 ) );
  all.join( Element( 1, PSI, 0 ) );
  all.join( Element( 2, PSI, 0 ) );
  all.join( Element( 3, PSI, 0 ) );
  all.join( Element( 4, PSI, 0 ) );
  all.join( Element( 5, PSI, 0 ) );
  all.join( Element( 6, PSI, 0 ) );

  Neighbour u0( Element( 0, PSI, 0 ) );
  u0.join( Element( 1, PSI, 0 ) );

  Neighbour u1( Element( 1, PSI, 0 ) );
  u1.join( Element( 0, PSI, 0 ) );
  u1.join( Element( 2, PSI, 0 ) );

  Neighbour u2( Element( 2, PSI, 0 ) );
  u2.join( Element( 1, PSI, 0 ) );
  u2.join( Element( 3, PSI, 0 ) );

  Neighbour u3( Element( 3, PSI, 0 ) );
  u3.join( Element( 2, PSI, 0 ) );
  u3.join( Element( 4, PSI, 0 ) );

  Neighbour u4( Element( 4, PSI, 0 ) );
  u4.join( Element( 3, PSI, 0 ) );
  u4.join( Element( 5, PSI, 0 ) );

  Neighbour u5( Element( 5, PSI, 0 ) );
  u5.join( Element( 4, PSI, 0 ) );
  u5.join( Element( 6, PSI, 0 ) );

  Neighbour u6( Element( 6, PSI, 0 ) );
  u6.join( Element( 5, PSI, 0 ) );

  Family T1( all );
  T1.join( u0 ).join( u1 ).join( u2 ).join( u3 ).join( u4 ).join( u5 ).join( u6 );

  std::vector< int > palette;
  palette.push_back( 0 );
  palette.push_back( 1 );
  Subjugator subj( all, T1, palette );

  Space sub0 = subj( 0 );
  Space sub1 = subj( 1 );

  std::cout << "sub0:" << std::endl;
  for ( Space::const_iterator it = sub0.begin(); it != sub0.end(); ++it )
  {
    Space::const_point pt( it );

    std::cout << "  ";
    if ( sub0.internal_contain( pt.element ) )
      std::cout << pt.index << " => external";
    else std::cout << pt.index;
    std::cout << std::endl;
  }

  std::cout << "sub1:" << std::endl;
  for ( Space::const_iterator it = sub1.begin(); it != sub1.end(); ++it )
  {
    Space::const_point pt( it );

    std::cout << "  ";
    if ( sub1.internal_contain( pt.element ) )
      std::cout << pt.index << " => external";
    else std::cout << pt.index;
    std::cout << std::endl;
  }

}

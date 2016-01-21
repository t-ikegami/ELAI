#include <iostream>
#include "space.hpp"
#include "family.hpp"
#include "metis.hpp"
using namespace std;

typedef int NID;
const NID WHOLE = -1;

typedef int VAR;
const VAR ANY  =15;
const VAR PSI  = 1;
//const VAR ELEC = 2;
const VAR HOLE = 4;
const VAR NUM  = 8;

class Element
{
  const NID i_;
  const VAR t_;

public:
  Element( NID i, VAR t ) : i_( i ), t_( t ) {}
  Element( const Element& x ) : i_( x.i_ ), t_( x.t_ ) {}
  ~Element() {}

  NID id() const { return i_; }
  VAR var() const { return t_; }

  bool operator<( const Element& x ) const
  {
    if ( i_ < x.i_ ) return true;
    return i_ == x.i_ && t_ < x.t_;
  }
  bool operator==( const Element& x ) const
  {
    return i_ == x.i_ && t_ == x.t_;
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

  Neighbour( NID i, VAR t ) : x_( Element( i, t ) ) {}
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
  Neighbour& join( NID i, VAR t )
  {
    assert( i != WHOLE && t != ANY );
    adj_.insert( Element( i, t ) ); return *this;
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
typedef elai::metis< Element, Neighbour > Metis;

int main()
{
  {
    Space s1;
    s1.join( Element( 0, PSI ) );
    s1.join( Element( 1, PSI ) );
    s1.join( Element( 2, PSI ) );
    s1.join( Element( 3, PSI ) );
    s1.join( Element( 4, PSI ) );
    s1.join( Element( 5, PSI ) );
    s1.join( Element( 6, PSI ) );
    s1.join( Element( 7, PSI ) );

    Neighbour u0psi( Element( 0, PSI ) );
    u0psi.join( Element( 1, PSI ) );
    u0psi.join( Element( 2, PSI ) );
    u0psi.join( Element( 4, PSI ) );

    Neighbour u1psi( Element( 1, PSI ) );
    u1psi.join( Element( 0, PSI ) );
    u1psi.join( Element( 3, PSI ) );
    u1psi.join( Element( 5, PSI ) );

    Neighbour u2psi( Element( 2, PSI ) );
    u2psi.join( Element( 0, PSI ) );
    u2psi.join( Element( 3, PSI ) );
    u2psi.join( Element( 6, PSI ) );

    Neighbour u3psi( Element( 3, PSI ) );
    u3psi.join( Element( 1, PSI ) );
    u3psi.join( Element( 2, PSI ) );
    u3psi.join( Element( 7, PSI ) );

    Neighbour u4psi( Element( 4, PSI ) );
    u4psi.join( Element( 0, PSI ) );
    u4psi.join( Element( 5, PSI ) );
    u4psi.join( Element( 6, PSI ) );

    Neighbour u5psi( Element( 5, PSI ) );
    u5psi.join( Element( 1, PSI ) );
    u5psi.join( Element( 4, PSI ) );
    u5psi.join( Element( 7, PSI ) );

    Neighbour u6psi( Element( 6, PSI ) );
    u6psi.join( Element( 2, PSI ) );
    u6psi.join( Element( 4, PSI ) );
    u6psi.join( Element( 7, PSI ) );

    Neighbour u7psi( Element( 7, PSI ) );
    u7psi.join( Element( 3, PSI ) );
    u7psi.join( Element( 5, PSI ) );
    u7psi.join( Element( 6, PSI ) );

    Family T1( s1 );
    T1.join( u0psi );
    T1.join( u1psi );
    T1.join( u2psi );
    T1.join( u3psi );
    T1.join( u4psi );
    T1.join( u5psi );
    T1.join( u6psi );
    T1.join( u7psi );

    Metis ord( s1, T1 );
    const Space& s = ord.reordered();

    cout << "Before(" << s1.size() << "):";
    for ( Space::const_iterator it = s1.begin(); it != s1.end(); ++it )
      cout << " " << Point( it ).index;
    cout << endl;
    cout << "After(" << s.size() << "):";
    for ( Space::const_iterator it = s.begin(); it != s.end(); ++it )
      cout << " " << Point( it ).index;
    cout << endl;
  }
}

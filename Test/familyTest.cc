#include <iostream>
#include "space.hpp"
#include "family.hpp"

using namespace std;

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

  const Element x_; // maybe opened point
  Adjacency adj_;   // only actual point

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

  Neighbour& join( const Element& e ) { adj_.insert( e ); return *this; }
  Neighbour& join( NID i, VAR t ) { adj_.insert( Element( i, t ) ); return *this; }

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
typedef elai::space< Element >::const_point Point;
typedef elai::family< Element, Neighbour > Family;

int main()
{
  {
    Space s1, s2, s3;

    cout << "compose s1.." << endl;
    s1.join( Element( 0, PSI ) );
    s1.join( Element( 0, NUM ) );
    s1.join( Element( 1, PSI ) );
    s1.join( Element( 1, ELEC) );
    s1.join( Element( 1, HOLE) );
    s1.join( Element( 2, PSI ) );
    s1.join( Element( 2, NUM ) );

    cout << "compose s2.." << endl;
    s2.join( Element( 0, PSI ) );
    s2.join( Element( 1, PSI ) );
    s2.join( Element( 2, PSI ) );

    Family psi( Neighbour( WHOLE, PSI ) );

    cout << "filter psi.." << endl;
    s3 = psi( s1 );
    cout << "psi:" << endl;
    for ( Space::const_iterator it = s3.begin(); it != s3.end(); ++it )
    {
      cout << Point( it ).element.id() << " " << Point( it ).element.var() << endl;
    }
    cout << endl;
    assert( s3 == s2 );
  }

  {
    Space s1;
    s1.join( Element( 0, PSI ) );
    s1.join( Element( 1, PSI ) );
    s1.join( Element( 1, ELEC) );
    s1.join( Element( 1, HOLE) );
    s1.join( Element( 2, PSI ) );
    s1.join( Element( 2, ELEC) );
    s1.join( Element( 2, HOLE) );
    s1.join( Element( 3, PSI ) );
    s1.join( Element( 3, ELEC) );
    s1.join( Element( 3, HOLE) );
    s1.join( Element( 4, PSI ) );
    s1.join( Element( 4, ELEC) );
    s1.join( Element( 4, HOLE) );

    Neighbour u0psi( Element( 0, PSI ) );
    u0psi.join( Element( 1, PSI ) );
    u0psi.join( Element( 1, ELEC) );
    u0psi.join( Element( 1, HOLE) );

    Neighbour u1psi( Element( 1, PSI ) );
    u1psi.join( Element( 0, PSI ) );
    u1psi.join( Element( 1, ELEC) );
    u1psi.join( Element( 1, HOLE) );
    u1psi.join( Element( 2, PSI ) );
    u1psi.join( Element( 2, ELEC) );
    u1psi.join( Element( 2, HOLE) );
    u1psi.join( Element( 3, PSI ) );
    u1psi.join( Element( 3, ELEC) );
    u1psi.join( Element( 3, HOLE) );
    u1psi.join( Element( 4, PSI ) );
    u1psi.join( Element( 4, ELEC) );
    u1psi.join( Element( 4, HOLE) );

    Neighbour u1elec( Element( 1, ELEC ) );
    u1elec.join( Element( 0, PSI ) );
    u1elec.join( Element( 1, PSI ) );
    u1elec.join( Element( 2, PSI ) );
    u1elec.join( Element( 2, ELEC) );
    u1elec.join( Element( 3, PSI ) );
    u1elec.join( Element( 3, ELEC) );
    u1elec.join( Element( 4, PSI ) );
    u1elec.join( Element( 4, ELEC) );

    Neighbour u1hole( Element( 1, HOLE ) );
    u1hole.join( Element( 0, PSI ) );
    u1hole.join( Element( 1, PSI ) );
    u1hole.join( Element( 2, PSI ) );
    u1hole.join( Element( 2, HOLE) );
    u1hole.join( Element( 3, PSI ) );
    u1hole.join( Element( 3, HOLE) );
    u1hole.join( Element( 4, PSI ) );
    u1hole.join( Element( 4, HOLE) );

    Neighbour u2psi( Element( 2, PSI ) );
    u2psi.join( Element( 1, PSI ) );
    u2psi.join( Element( 1, ELEC) );
    u2psi.join( Element( 1, HOLE) );
    u2psi.join( Element( 2, ELEC) );
    u2psi.join( Element( 2, HOLE) );
    u2psi.join( Element( 4, PSI ) );
    u2psi.join( Element( 4, ELEC) );
    u2psi.join( Element( 4, HOLE) );

    Neighbour u2elec( Element( 2, ELEC ) );
    u2elec.join( Element( 1, PSI ) );
    u2elec.join( Element( 1, ELEC) );
    u2elec.join( Element( 2, PSI ) );
    u2elec.join( Element( 4, PSI ) );
    u2elec.join( Element( 4, ELEC) );

    Neighbour u2hole( Element( 2, HOLE ) );
    u2hole.join( Element( 1, PSI ) );
    u2hole.join( Element( 1, HOLE) );
    u2hole.join( Element( 2, PSI ) );
    u2hole.join( Element( 4, PSI ) );
    u2hole.join( Element( 4, HOLE) );

    Neighbour u3psi( Element( 3, PSI ) );
    u3psi.join( Element( 1, PSI ) );
    u3psi.join( Element( 1, ELEC) );
    u3psi.join( Element( 1, HOLE) );
    u3psi.join( Element( 3, ELEC) );
    u3psi.join( Element( 3, HOLE) );
    u3psi.join( Element( 4, PSI ) );
    u3psi.join( Element( 4, ELEC) );
    u3psi.join( Element( 4, HOLE) );

    Neighbour u3elec( Element( 3, ELEC ) );
    u3elec.join( Element( 1, PSI ) );
    u3elec.join( Element( 1, ELEC) );
    u3elec.join( Element( 3, PSI ) );
    u3elec.join( Element( 4, PSI ) );
    u3elec.join( Element( 4, ELEC) );

    Neighbour u3hole( Element( 3, HOLE ) );
    u3hole.join( Element( 1, PSI ) );
    u3hole.join( Element( 1, HOLE) );
    u3hole.join( Element( 3, PSI ) );
    u3hole.join( Element( 4, PSI ) );
    u3hole.join( Element( 4, HOLE) );

    Neighbour u4psi( Element( 4, PSI ) );
    u4psi.join( Element( 1, PSI ) );
    u4psi.join( Element( 1, ELEC) );
    u4psi.join( Element( 1, HOLE) );
    u4psi.join( Element( 2, PSI ) );
    u4psi.join( Element( 2, ELEC) );
    u4psi.join( Element( 2, HOLE) );
    u4psi.join( Element( 3, PSI ) );
    u4psi.join( Element( 3, ELEC) );
    u4psi.join( Element( 3, HOLE) );
    u4psi.join( Element( 4, ELEC) );
    u4psi.join( Element( 4, HOLE) );

    Neighbour u4elec( Element( 4, ELEC ) );
    u4elec.join( Element( 1, PSI ) );
    u4elec.join( Element( 1, ELEC) );
    u4elec.join( Element( 2, PSI ) );
    u4elec.join( Element( 2, ELEC) );
    u4elec.join( Element( 3, PSI ) );
    u4elec.join( Element( 3, ELEC) );
    u4elec.join( Element( 4, PSI ) );

    Neighbour u4hole( Element( 4, HOLE ) );
    u4hole.join( Element( 1, PSI ) );
    u4hole.join( Element( 1, HOLE) );
    u4hole.join( Element( 2, PSI ) );
    u4hole.join( Element( 2, HOLE) );
    u4hole.join( Element( 3, PSI ) );
    u4hole.join( Element( 3, HOLE) );
    u4hole.join( Element( 4, PSI ) );

    Family T1( s1 );
    T1.join( u0psi );
    T1.join( u1psi ).join( u1elec ).join( u1hole );
    T1.join( u2psi ).join( u2elec ).join( u2hole );
    T1.join( u3psi ).join( u3elec ).join( u3hole );
    T1.join( u4psi ).join( u4elec ).join( u4hole );

    cout << "Check numbers of adjacencies in s1:" << endl;

    cout << T1( Element( 0, PSI ) ).size() << endl;
    assert( T1( Element( 0, PSI ) ).size() == 4 );

    cout << T1( Element( 1, PSI ) ).size() << endl;
    cout << T1( Element( 1, ELEC) ).size() << endl;
    cout << T1( Element( 1, HOLE) ).size() << endl;
    assert( T1( Element( 1, PSI ) ).size() == 13 );
    assert( T1( Element( 1, ELEC) ).size() == 9 );
    assert( T1( Element( 1, HOLE) ).size() == 9 );

    cout << T1( Element( 2, PSI ) ).size() << endl;
    cout << T1( Element( 2, ELEC) ).size() << endl;
    cout << T1( Element( 2, HOLE) ).size() << endl;
    assert( T1( Element( 2, PSI ) ).size() == 9 );
    assert( T1( Element( 2, ELEC) ).size() == 6 );
    assert( T1( Element( 2, HOLE) ).size() == 6 );

    cout << T1( Element( 3, PSI ) ).size() << endl;
    cout << T1( Element( 3, ELEC) ).size() << endl;
    cout << T1( Element( 3, HOLE) ).size() << endl;
    assert( T1( Element( 3, PSI ) ).size() == 9 );
    assert( T1( Element( 3, ELEC) ).size() == 6 );
    assert( T1( Element( 3, HOLE) ).size() == 6 );

    cout << T1( Element( 4, PSI ) ).size() << endl;
    cout << T1( Element( 4, ELEC) ).size() << endl;
    cout << T1( Element( 4, HOLE) ).size() << endl;
    assert( T1( Element( 4, PSI ) ).size() == 12 );
    assert( T1( Element( 4, ELEC) ).size() == 8 );
    assert( T1( Element( 4, HOLE) ).size() == 8 );

    cout << endl;
  }
}

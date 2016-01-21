#include <iostream>
#include <map>
#include <set>
#include "vector.hpp"
#include "linear_function.hpp"

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

  Neighbour& operator=( Neighbour& src )
  {
    if ( this == &src ) return *this;
    if ( x_ == src.x_ ) adj_ = src.adj_;
    return *this;
  }

  const Element& element() const { return x_; }

  Neighbour& join( NID i, VAR t )
  {
    if ( i != WHOLE && t != ANY ) adj_.insert( Element( i, t ) );
    return *this;
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
typedef elai::space< Element >::point Point;
typedef elai::family< Element, Neighbour > Family;
typedef elai::vector< float > Vector;
typedef elai::linear_function< Element, Neighbour, float > Function;

int main()
{
  {
    Space omega;

    omega.join( Element( 0, NUM  ) );
    omega.join( Element( 1, PSI  ) );
    omega.join( Element( 2, PSI  ) );
    omega.join( Element( 2, ELEC ) );
    omega.join( Element( 2, HOLE ) );
    omega.join( Element( 3, PSI  ) );
    omega.join( Element( 3, ELEC ) );
    omega.join( Element( 3, HOLE ) );
    omega.join( Element( 4, PSI  ) );
    omega.join( Element( 4, NUM  ) );

    Function f( omega );

    Function psi = f.localize( Neighbour( WHOLE, PSI  ) );
    Function elec= f.localize( Neighbour( WHOLE, ELEC ) );
    Function hole= f.localize( Neighbour( WHOLE, HOLE ) );
    Function num = f.localize( Neighbour( WHOLE, NUM  ) );

    { Vector& v = psi.ran();
      for ( int i = 0; i < v.m(); ++i ) v( i ) = 1.;
      cout << "psi =" << v;
    }
    { Vector& v = elec.ran();
      for ( int i = 0; i < v.m(); ++i ) v( i ) = 2.;
      cout << "elec =" << v;
    }
    { Vector& v = hole.ran();
      for ( int i = 0; i < v.m(); ++i ) v( i ) = 3.;
      cout << "hole =" << v;
    }
    { Vector& v = num.ran();
      for ( int i = 0; i < v.m(); ++i ) v( i ) = 4.;
      cout << "num =" << v;
    }

    f.reflectIn( psi ).reflectIn( elec ).reflectIn( hole ).reflectIn( num );
    cout << f.ran();
  }

  {
    Space s1;
    s1.join( Element( 0, PSI ) );
    s1.join( Element( 1, PSI ) );
    s1.join( Element( 2, PSI ) );
    s1.join( Element( 3, PSI ), Element( 2, NUM ) );
    s1.join( Element( 4, PSI ), Element( 3, NUM ) );
    Function f1( s1 );

    Space s2;
    s2.join( Element( 0, NUM ), Element( 1, PSI ) );
    s2.join( Element( 1, NUM ), Element( 2, PSI ) );
    s2.join( Element( 2, NUM ) );
    s2.join( Element( 3, NUM ) );
    s2.join( Element( 4, NUM ) );
    s2.join( Element( 5, NUM ) );
    Function f2( s2 );

    Function f = f1.extend( f2 );
    Vector &v = f.ran();
    for ( int i = 0; i < f.dim(); ++i ) v( i ) = i;
    cout << v;
    f1.reflect( f );
    cout << f1.ran();
    f2.reflect( f );
    cout << f2.ran();
  }
}

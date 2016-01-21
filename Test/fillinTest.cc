#include <iostream>
#include <map>
#include <set>
#include "vector.hpp"
#include "matrix.hpp"
#include "linear_function.hpp"
#include "linear_operator.hpp"
#include "fillin.hpp"

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

  Neighbour& erase( const Element& x )
  {
    adj_.erase( x );
    return *this;
  }
  Neighbour& erase( const_iterator& it )
  {
    adj_.erase( it );
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
typedef elai::matrix< float > Matrix;
typedef elai::linear_operator< Element, Neighbour, float > Operator;
typedef elai::fillin< float > Fillin;

int main()
{
  {
    Space s1;
    s1.join( Element( 1, PSI ) );
    s1.join( Element( 2, PSI ) );
    s1.join( Element( 3, PSI ) );
    s1.join( Element( 4, PSI ) );
    s1.join( Element( 5, PSI ) );
    s1.join( Element( 6, PSI ) );
    s1.join( Element( 7, PSI ) );

    Neighbour u1( Element( 1, PSI ) );
    u1.join( Element( 3, PSI ) ).join( Element( 4, PSI ) );

    Neighbour u2( Element( 2, PSI ) );
    u2.join( Element( 3, PSI ) ).join( Element( 5, PSI ) ).join( Element( 7, PSI ) );

    Neighbour u3( Element( 3, PSI ) );
    u3.join( Element( 1, PSI ) ).join( Element( 2, PSI ) );

    Neighbour u4( Element( 4, PSI ) );
    u4.join( Element( 1, PSI ) ).join( Element( 6, PSI ) ).join( Element( 7, PSI ) );

    Neighbour u5( Element( 5, PSI ) );
    u5.join( Element( 2, PSI ) );

    Neighbour u6( Element( 6, PSI ) );
    u6.join( Element( 4, PSI ) );

    Neighbour u7( Element( 7, PSI ) );
    u7.join( Element( 2, PSI ) ).join( Element( 4, PSI ) );

    Family T1( s1 );
    T1.join( u1 ).join( u2 ).join( u3 ).join( u4 ).join( u5 ).join( u6 ).join( u7 );

    Operator A( s1, T1 );
    Fillin f( A.action() );
    f( 2 );
    Matrix a( f.m(), f.n(), f.nnz(), f.xadj(), f.adjy(), f.coef() );
    cout << A.action();
    cout << a;
  }
}

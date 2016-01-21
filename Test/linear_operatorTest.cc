#include <iostream>
#include <map>
#include <set>
#include "vector.hpp"
#include "matrix.hpp"
#include "linear_function.hpp"
#include "linear_operator.hpp"

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

int main()
{
  {
    cout << "Test1" << endl;
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

    Operator A( s1, T1 );

    Operator poisson( A.localize( Neighbour( WHOLE, PSI ) ) );
    //Operator psi( A.localize( Neighbour( WHOLE, PSI ), Neighbour( WHOLE, PSI ) ) );
    {
      poisson( Element( 0, PSI ) ) = 1e0;
      poisson( Element( 1, PSI ) ) = 1e0;
      poisson( Element( 2, PSI ) ) = 1e0;
      poisson( Element( 3, PSI ) ) = 1e0;
      poisson( Element( 4, PSI ) ) = 1e0;
    }

    Operator connE( A.localize( Neighbour( WHOLE, ELEC ) ) );
    //Operator elec( A.localize( Neighbour( WHOLE, ELEC ), Neighbour( WHOLE, ELEC ) ) );
    {
      connE( Element( 1, ELEC) ) = 1e-2;
      connE( Element( 2, ELEC) ) = 1e-2;
      connE( Element( 3, ELEC) ) = 1e-2;
      connE( Element( 4, ELEC) ) = 1e-2;
    }
    //cout << connE.action();

    Operator connH( A.localize( Neighbour( WHOLE, HOLE ) ) );
    //Operator hole( A.localize( Neighbour( WHOLE, HOLE ), Neighbour( WHOLE, HOLE ) ) );
    {
      connH( Element( 1, HOLE) ) = 5e-2;
      connH( Element( 2, HOLE) ) = 5e-2;
      connH( Element( 3, HOLE) ) = 5e-2;
      connH( Element( 4, HOLE) ) = 5e-2;
    }

    A.reflect( poisson ).reflect( connE ).reflect( connH );
    cout << A.action();
    Operator thermal( A.localize( Neighbour( WHOLE, PSI ), Neighbour( WHOLE, PSI ) ) );
    cout << thermal.action();
    Operator ddE( A.localize( Neighbour( WHOLE, ELEC ), Neighbour( WHOLE, ELEC ) ) );
    cout << ddE.action();
    Operator ddH( A.localize( Neighbour( WHOLE, HOLE ), Neighbour( WHOLE, HOLE ) ) );
    cout << ddH.action();
  }

  {
    cout << "Test2" << endl;
    Space s1;
    s1.join( Element( 0, PSI ) );
    s1.join( Element( 1, PSI ) );
    s1.join( Element( 2, PSI ) );
    Neighbour u1_0( Element( 0, PSI ) );
    u1_0.join( Element( 1, PSI ) );
    Neighbour u1_1( Element( 1, PSI ) );
    u1_1.join( Element( 0, PSI ) ).join( Element( 2, PSI ) );
    Neighbour u1_2( Element( 2, PSI ) );
    u1_2.join( Element( 1, PSI ) );
    Family T1( s1 );
    T1.join( u1_0 ).join( u1_1 ).join( u1_2 );
    Function f1( s1 );
    Operator A1( s1, T1 );

    Space s2;
    s2.join( Element( 1, PSI ) );
    s2.join( Element( 2, PSI ) );
    s2.join( Element( 3, PSI ) );
    Neighbour u2_1( Element( 1, PSI ) );
    u2_1.join( Element( 2, PSI ) );
    Neighbour u2_2( Element( 2, PSI ) );
    u2_2.join( Element( 1, PSI ) ).join( Element( 3, PSI ) );
    Neighbour u2_3( Element( 3, PSI ) );
    u2_3.join( Element( 2, PSI ) );
    Family T2( s2 );
    T2.join( u2_1 ).join( u2_2 ).join( u2_3 );
    Function f2( s2 );
    Operator A2( s2, T2 );

    Space s1_in, s2_in;
    s1_in.join( Element( 0, PSI ) );
    s1_in.join( Element( 1, PSI ) );
    s2_in.join( Element( 2, PSI ) );
    s2_in.join( Element( 3, PSI ) );
    Function f1_ = f1.localize( s1_in );
    Function f2_ = f2.localize( s2_in );
    Function f = f1_.extend( f2_ );
    cout << f1_.ran();
    cout << f2_.ran();
    cout << f.ran();

    Operator A1_ = A1.localize( s1_in );
    Operator A2_ = A2.localize( s2_in );
    Operator A = A1_.extend( A2_ );
    cout << A1_.action();
    cout << A2_.action();
    cout << A.action();
  }
}

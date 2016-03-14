/*
 *
 * Enexss Linear Algebra Interface (ELAI)
 *
 * Copyright 2013-2015 H. KOSHIMOTO, AIST
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#ifndef __ELAI_LINEAR_OPERATOR__
#define __ELAI_LINEAR_OPERATOR__

#include <vector>
#include "space.hpp"
#include "family.hpp"
#include "portal.hpp"
#include "subjugator.hpp"
#include "linear_function.hpp"
#include "matrix.hpp"

namespace elai
{

/*
 * Infinite Vector Space:
 *  A( f ) = g
 *  A; Linear Operator: function -> function
 *  f; Solution: element -> scalar
 *  g; Inhomogeneous Term: function -> scalar
 *
 * Finite Vector Space:
 *  A f = g
 *  A; Matrix: i -> ( j -> scalar )
 *  f; Vector: j -> scalar
 *  g; Vector: i -> scalar
 *
 * Connection Space:
 *  for solution:
 *    element -> j
 *    j -> element
 *  for function:
 *    function -> i
 *    i -> function
 *
 * Special Case:
 *  same elements on solution and inhomogeneous terms.
 *  generally speaking, an elements decides it's value (=solution) and it's function.
 *
 */
template< class Element, class Neighbour, class Range >
class linear_operator
{
  friend class portal;

public:
  typedef Range range;
  typedef vector< range > Vector;
  typedef matrix< range > Matrix;
  typedef space< Element > Space;
  typedef family< Element, Neighbour > Family;

private:
  typedef struct Space::const_point Point;

  const Space f_;
  const Space x_;
  const Family tau_;
  Matrix A_;

#ifdef ELAI_USE_MPI
  void marshalize( const portal& port ) const
  {
    port( f_ );
    port( x_ );
    port( tau_ );
    port.send( A_.val(), mpi_< range >().type, A_.nnz() );
  }

  void reinforce( const portal& port ) const
  {
    port.send( A_.val(), mpi_< range >().type, A_.nnz() );
  }
#endif

  void setup()
  {
    int *ind, *col, m, n, nnz, off;
    std::vector< Space > adjs;
    range *c;

    m = f_.size();
    ind = new int[ m + 1 ];

    ind[ 0 ] = 0; off = 0;
    for ( typename Space::const_iterator it = f_.begin(); it != f_.end(); ++it )
    {
      adjs.push_back( tau_( Point( it ).element ) );
      ind[ Point( it ).index + 1 ] = adjs[ off++ ].size();
    }
    for ( int i = 0; i < m; ++i ) ind[ i + 1 ] += ind[ i ];

    nnz = ind[ m ];
    col = new int[ nnz ];
    c = new range[ nnz ];

    off = 0; n = 0;
    for ( typename Space::const_iterator it = f_.begin(); it != f_.end(); ++it )
    {
      const Space &adj = adjs[ off++ ];
      int i = Point( it ).index, k = 0;
      std::set< int > jset;

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        int j = x_.index( Point( jt ).element ); // 'j's are not ordered.

        if ( n < j ) n = j;
        jset.insert( j ); // 'j's are ordered in jset.
      }
      for ( std::set< int >::const_iterator j = jset.begin(); j != jset.end(); ++j )
      {
        col[ ind[ i ] + k ] = *j;
        c[ ind[ i ] + k ] = static_cast< range >( 0 );
        ++k;
      }
    }
    A_.setup( m, n + 1, nnz, ind, col, c );

    delete [] c;
    delete [] col;
    delete [] ind;
  }

public:
  linear_operator( const Space& f, const Family& tau )
    : f_( f ), x_( f ), tau_( tau )
  {
    setup();
  }

  linear_operator
    ( const Space& f
    , const Space& x
    , const Family& tau
    )
    : f_( f ), x_( x ), tau_( tau )
  {
    setup();
  }

  linear_operator
    ( const Space& f
    , const Space& x
    , const Family& tau
    , const Matrix& A
    )
    : f_( f ), x_( x ), tau_( tau ), A_( A )
  {}

  linear_operator( const linear_operator< Element, Neighbour, Range >& src )
    : f_( src.f_ ), x_( src.x_ ), tau_( src.tau_ ), A_( src.A_ ) {}

#ifdef ELAI_USE_MPI
  linear_operator( portal& port )
    : f_( port ), x_( port ), tau_( port ), A_()
  {
    setup();
    port.recv( A_.val(), mpi_< range >().type, A_.nnz() );
  }

  void operator()( portal& port )
  {
    port.recv( A_.val(), mpi_< range >().type, A_.nnz() );
  }
#endif

  ~linear_operator()
  {}

  int dim() const { return f_.size(); }
  int codim() const { return x_.size(); }
  const Space& dom() const { return x_; }
  const Space& ran() const { return f_; }
  const Family& topo() const { return tau_; }
  Matrix& action() { return A_; }
  const Matrix& action() const { return A_; }

  range& operator()( const Element& f )
  { return A_( f_.index( f ), x_.index( f ) ); }
  const range& operator()( const Element& f ) const
  { return A_( f_.index( f ), x_.index( f ) ); }
  range& operator()( const Element& f, const Element& x )
  { return A_( f_.index( f ), x_.index( x ) ); }
  const range& operator()( const Element& f, const Element& x ) const
  { return A_( f_.index( f ), x_.index( x ) ); }

#ifdef ELAI_USE_PYTHON
  void setMatrix(const Element& f, range val)
  {
     A_( f_.index( f ), x_.index( f ) ) = val;
  }

  void setMatrix(const Element& f, const Element& x, range val)
  {
    A_( f_.index( f ), x_.index( x ) ) = val;
  }
#endif

  linear_operator< Element, Neighbour, Range >& clear( const Range c )
  {
    A_.clear( c );

    return *this;
  }
  linear_operator< Element, Neighbour, Range >& clear( const Range c, const Space& s )
  {
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
      for ( typename Space::const_iterator jt = x_.begin(); jt != x_.end(); ++jt )
        A_( f_.index( Point( it ).element ), x_.index( Point( jt ).element ) ) = c;

    return *this;
  }
  linear_operator< Element, Neighbour, Range >& clear( const Range c, const Space& s, const Space& t )
  {
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
      for ( typename Space::const_iterator jt = t.begin(); jt != t.end(); ++jt )
        A_( f_.index( Point( it ).element ), x_.index( Point( jt ).element ) ) = c;

    return *this;
  }

  // Row-wise localizeation
  linear_operator< Element, Neighbour, Range > localize( const Space& s ) const
  {
    return localize( Family( s ) );
  }
  linear_operator< Element, Neighbour, Range > localize( const Neighbour& u ) const
  {
    return localize( Family( u ) );
  }
  linear_operator< Element, Neighbour, Range > localize( const Family& tau ) const
  {
#ifdef ELAI_USE_C11
    Space&& s = std::move( tau( f_ ) ); // reordered
    Family&& theta = std::move( tau_.localize( s ) );
#else
    Space s = tau( f_ ); // reordered
    Family theta = tau_.localize( s );
#endif
    std::vector< Space > adjs;
    int m = s.size(), n = 0;
    int *ind = new int[ m + 1 ], offset;

    offset = 0;
    ind[ 0 ] = 0;
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      adjs.push_back( theta( Point( it ).element ) );
      ind[ Point( it ).index + 1 ] = adjs[ offset++ ].size();
    }
    for ( int i = 0; i < m; ++i ) ind[ i + 1 ] += ind[ i ];

    int nnz = ind[ m ];
    int *col = new int[ nnz ];
    range *c = new range[ nnz ];

    offset = 0;
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      Space& adj = adjs[ offset++ ];
      int i = Point( it ).index; // local i
      int I = f_.index( Point( it ).element ); // global i
      std::set< int > jset;
      int k = 0;

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        int J = x_.index( Point( jt ).element ); // global j, are not ordered.

        if ( n < J ) n = J;
        jset.insert( J ); // global j, are ordered in jset.
      }
      for ( std::set < int >::const_iterator J = jset.begin(); J != jset.end(); ++J )
      {
        col[ ind[ i ] + k ] = *J;
        c[ ind[ i ] + k ] = A_( I, *J );
        ++k;
      }
    }
    Matrix a( m, n + 1, nnz, ind, col, c );

    delete [] c;
    delete [] col;
    delete [] ind;

    return linear_operator< Element, Neighbour, Range >( s, x_, theta, a );
  }

  // Row/Col localizeation
  linear_operator< Element, Neighbour, Range > localize
    ( const Space& s
    , const Space& t
    ) const
  {
    return localize( Family( s ), Family( t ) );
  }
  linear_operator< Element, Neighbour, Range > localize
    ( const Neighbour& u
    , const Neighbour& v
    ) const
  {
    return localize( Family( u ), Family( v ) );
  }
  linear_operator< Element, Neighbour, Range > localize
    ( const Family& sigma
    , const Family& tau
    ) const
  {
#ifdef ELAI_USE_C11
    Space&& s = std::move( sigma( f_ ) ); // reordered
    Space&& t = std::move( tau( x_ ) ); // reordered
    Family&& theta = std::move( tau_.localize( s ) );
#else
    Space s = sigma( f_ ); // reordered
    Space t = tau( x_ ); // reordered
    Family theta = tau_.localize( s );
#endif
    std::vector< Space > adjs;
    int m = s.size(), n = t.size();
    int *ind = new int[ m + 1 ], offset;

    offset = 0;
    ind[ 0 ] = 0;
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
#ifdef ELAI_USE_C11
      Space&& candidates = std::move( theta( Point( it ).element ) );
#else
      Space candidates = theta( Point( it ).element );
#endif
      adjs.push_back( candidates & t );
      ind[ Point( it ).index + 1 ] = adjs[ offset++ ].size();
    }
    for ( int i = 0; i < m; ++i ) ind[ i + 1 ] += ind[ i ];

    int nnz = ind[ m ];
    int *col = new int[ nnz ];
    range *c = new range[ nnz ];

    offset = 0;
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      Space& adj = adjs[ offset++ ];
      int i = Point( it ).index; // local i
      int I = f_.index( Point( it ).element ); // global i
      int k = 0;

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& e = Point( jt ).element;
        int j = t.index( e ); // local j, are ordered in adj.
        int J = x_.index( e ); // global j

        col[ ind[ i ] + k ] = j;
        c[ ind[ i ] + k ] = A_( I, J );
        ++k;
      }
    }
    Matrix a( m, n, nnz, ind, col, c );

    delete [] c;
    delete [] col;
    delete [] ind;

    return linear_operator< Element, Neighbour, Range >( s, t, theta, a );
  }

  linear_operator< Element, Neighbour, Range >& reflectIn
    ( const linear_operator< Element, Neighbour, Range >& B )
  {
    for ( typename Space::const_iterator it = B.f_.begin(); it != B.f_.end(); ++it )
    {
      const Element& d = Point( it ).element;
      int i = Point( it ).index;
      int I = -1;

      if ( f_.govern( d ) ) I = f_.index( d );
      else if ( f_.external_contain( d ) ) I = f_.external_index( d );

      if ( I < 0 ) continue;
      const Space& adj = B.tau_( d );

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& e = Point( jt ).element;
        int j = B.x_.index( e );
        int J = -1;

        if ( x_.govern( e ) ) J = x_.index( e );
        else if ( x_.external_contain( e ) ) J = x_.external_index( e );

        if ( 0 <= J ) A_( I, J ) = B.A_( i, j );
      }
    }

    return *this;
  }

  linear_operator< Element, Neighbour, Range >& reflectIn
    ( const linear_operator< Element, Neighbour, Range >& B
    , const subjugator< Element, Neighbour >& subj
    )
  {
    for ( typename Space::const_iterator it = B.f_.begin(); it != B.f_.end(); ++it )
    {
      const Element& d = Point( it ).element;
      int i = Point( it ).index;
      int I = -1;
      int c = subj.color( d );
      const Element x( d, c );

      if ( f_.govern( x ) ) I = f_.index( x );
      else if ( f_.external_contain( x ) ) I = f_.external_index( x );

      if ( I < 0 ) continue;
      const Space& adj = B.tau_( d );

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& e = Point( jt ).element;
        int j = B.x_.index( e );
        int J = -1;
        int c = subj.color( e );
        const Element y( e, c );

        if ( x_.govern( y ) ) J = x_.index( y );
        else if ( x_.external_contain( y ) ) J = x_.external_index( y );

        if ( 0 <= J ) A_( I, J ) = B.A_( i, j );
      }
    }

    return *this;
  }

  linear_operator< Element, Neighbour, Range >& reflect
    ( const linear_operator< Element, Neighbour, Range >& B )
  {
    for ( typename Space::const_iterator it = f_.begin(); it != f_.end(); ++it )
    {
      const Element& d = Point( it ).element;
      int I = Point( it ).index;
      int i = -1;

      if ( B.f_.govern( d ) ) i = B.f_.index( d );
      //else if ( B.f_.external_contain( d ) ) i = B.f_.external_index( d );

      if ( i < 0 ) continue;
      const Space& adj = tau_( d );

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& e = Point( jt ).element;
        int J = x_.index( e );
        int j = -1;

        if ( B.x_.govern( e ) ) j = B.x_.index( e );
        //else if ( B.x_.external_contain( e ) ) j = B.x_.external_index( e );

        if ( 0 <= j ) A_( I, J ) = B.A_( i, j );
      }
    }

    return *this;
  }

  linear_operator< Element, Neighbour, Range >& reflect
    ( const linear_operator< Element, Neighbour, Range >& B
    , const subjugator< Element, Neighbour >& subj
    )
  {
    for ( typename Space::const_iterator it = f_.begin(); it != f_.end(); ++it )
    {
      const Element& d = Point( it ).element;
      int I = Point( it ).index;
      int i = -1;
      int c = subj.color( d );
      const Element x( d, c );

      if ( B.f_.govern( x ) ) i = B.f_.index( x );
      //else if ( B.f_.external_contain( x ) ) i = B.f_.external_index( x );

      if ( i < 0 ) continue;
      const Space& adj = tau_( d );

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& e = Point( jt ).element;
        int J = x_.index( e );
        int j = -1;
        int c = subj.color( e );
        const Element y( e, c );

        if ( B.x_.govern( y ) ) j = B.x_.index( y );
        //else if ( B.x_.external_contain( y ) ) j = B.x_.external_index( y );

        if ( 0 <= j ) A_( I, J ) = B.A_( i, j );
      }
    }

    return *this;
  }

  // 'extend' method changes a non-zero structure.
  linear_operator< Element, Neighbour, Range > extend
    ( const linear_operator< Element, Neighbour, Range >& B ) const
  {
#ifdef ELAI_USE_C11
    Space&& s = std::move( f_ | B.f_ ); // reordered
    Space&& t = std::move( x_ | B.x_ ); // reordered
    Family&& tau = std::move( B.tau_ );
#else
    Space s = f_ | B.f_; // reordered
    Space t = x_ | B.x_; // reordered
    Family tau = B.tau_;
#endif
    std::vector< Space > adjs;
    int *ind = new int[ s.size() + 1 ], m = s.size(), offset;

    for ( typename Space::const_marginal_iterator it = B.f_.internal_begin()
        ; it != B.f_.internal_end(); ++it
        )
    { // exchange from internal elements to external elements
      if ( f_.internal_contain( typename Space::const_internal_point( it ).external ) )
        tau.flip
          ( typename Space::const_internal_point( it ).internal
          , typename Space::const_internal_point( it ).external
          );
    }
#ifdef ELAI_USE_C11
    tau = std::move( tau_ | tau );
#else
    tau = tau_ | tau;
#endif

    offset = 0;
    ind[ 0 ] = 0;
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      const Element& e = Point( it ).element;
      adjs.push_back( tau( e ) );
      ind[ Point( it ).index + 1] = adjs[ offset++ ].size();
    }
    for ( int i = 0; i < m; ++i ) ind[ i + 1 ] += ind[ i ];

    int nnz = ind[ m ];
    int *col = new int[ nnz ];
    range *c = new range[ nnz ];

    offset = 0;
    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      // ?? Space& adj = adjs[ offset++ ];
      const Space& adj = adjs[ offset++ ];
      const Element& e = Point( it ).element;
      int i = Point( it ).index; // global i
      int k = 0;

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        const Element& d = Point( jt ).element;
        int j = t.index( d ); // global j
        int I, J; // local i, j depend on spaces

        if ( f_.contain( e ) && x_.contain( d ) )
        {
          I = f_.index( e );
          J = x_.index( d );
          c[ ind[ i ] + k ] = A_( I, J );
        }
        else if ( f_.contain( e ) )
        {
          I = B.f_.external_index( e );
          J = B.x_.index( d );
          c[ ind[ i ] + k ] = B.A_( I, J );
        }
        else
        {
          if ( B.f_.external_contain( e ) ) I = B.f_.external_index( e );
          else I = B.f_.index( e );
          if ( B.x_.external_contain( d ) ) J = B.f_.external_index( d );
          else J = B.x_.index( d );
          c[ ind[ i ] + k ] = B.A_( I, J );
        }
        col[ ind[ i ] + k ] = j;
        ++k;
      }
    }

    Matrix A( s.size(), t.size(), nnz, ind, col, c );

    delete [] c;
    delete [] col;
    delete [] ind;

    return linear_operator< Element, Neighbour, Range >( s, t, tau, A );
  }
};

}

#endif//__ELAI_LINEAR_OPERATOR__

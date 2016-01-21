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
#ifndef __ELAI_GENERATOR__
#define __ELAI_GENERATOR__

#include <set>
#include "def.hpp"
#include "space.hpp"
#include "family.hpp"
#include "matrix.hpp"

namespace elai
{

template< class Coef >
class generator
{
public:
  class Element
  {
    const int i_, c_;

  public:
    typedef long long id_type;

    Element( const id_type id )
      : i_( id & 0x00000000ffffffff )
      , c_( ( id & 0xffffffff00000000 ) >> 32 )
    {}
    Element( int i, int c ) : i_( i ), c_( c ) {}
    Element( const Element& x ) : i_( x.i_ ), c_( x.c_ ) {}
    Element( const Element& x, int c ) : i_( x.i_ ), c_( c ) {}
    ~Element() {}

    id_type operator()() const
    {
      return ( static_cast< id_type >( c_ ) << 32 ) | ( static_cast< id_type >( i_ ) );
    }
    int color() const { return c_; }

    bool operator<( const Element& x ) const
    {
      if ( c_ == x.c_ ) return i_ < x.i_;
      return c_ < x.c_;
    }
    bool operator==( const Element& x ) const { return i_ == x.i_ && c_ == x.c_; }
  };

  class Neighbour
  {
    typedef std::set< Element > Adjacency;

    const Element x_;
    Adjacency adj_;

  public:
    typedef typename Adjacency::iterator iterator;
    typedef typename Adjacency::const_iterator const_iterator;

    Neighbour( int i, int c ) : x_( Element( i, c ) ) { adj_.insert( x_ ); }
    Neighbour( const Element& x ) : x_( x ) { adj_.insert( x_ ); }
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

    Neighbour& join( int i, int c )
    {
      adj_.insert( Element( i, c ) );
      return *this;
    }
    Neighbour& join( const Element& y )
    {
      adj_.insert( y );
      return *this;
    }

    Neighbour& erase( const Element& y )
    {
      adj_.erase( y );
      return *this;
    }
    Neighbour& erase( const_iterator& it )
    {
      adj_.erase( it );
      return *this;
    }

    bool operator()( const Element& y ) const
    {
      for ( const_iterator it = adj_.begin(); it != adj_.end(); ++it )
        if ( *it == y ) return true;
      return false;
    }

    iterator begin() { return adj_.begin(); }
    const_iterator begin() const { return adj_.begin(); }
    const_iterator end() const { return adj_.end(); }
  };

  typedef matrix< Coef > Matrix;
  typedef elai::space< Element > Space;
  typedef elai::family< Element, Neighbour > Family;

  Element element( const typename Element::id_type id ) const { return Element( id ); }
  Element element( int i, int c ) const { return Element( i, c ); }
  Element element( const Element& x ) const { return Element( x ); }
  Element element( const Element& x, int c ) const { return Element( x, c ); }
  Neighbour neighbour( int i, int c ) const { return Neighbour( i, c ); }
  Neighbour neighbour( const Element& x ) const { return Neighbour( x ); }
  Neighbour neighbour( const Neighbour& u ) const { return Neighbour( u ); }

private:
  const Matrix A_;
  Space base_;
  Family topo_;

public:
  generator( const Matrix& A ) : A_( A )
  {
    const int c = 0;
    const int *ind = A_.ind(), *col = A_.col();

    for ( int i = 0; i < A_.m(); ++i ) base_.join( Element( i, c ) );
    for ( int i = 0; i < A_.m(); ++i )
    {
      Neighbour u( i, c );

      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        u.join( j, c );
      }
      topo_.join( u );
    }
  }
  ~generator() {}

  const Space& space() const { return base_; }
  const Family& family() const { return topo_; }
};

}

#endif//__ELAI_GENERATOR__

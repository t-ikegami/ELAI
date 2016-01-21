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
#ifndef __ELAI_SPACE__
#define __ELAI_SPACE__

#include <map>
#include <iostream>
#include <utility>
#include "def.hpp"
#include "portal.hpp"

namespace elai
{

/*
 * Element Requirements:
 *  - No arguments constructor returns the minimum item.
 *  - Unique ID Type: typedef ??? id_type
 *  -                 id_type must be comparable.
 *  - Copy constructors: Element( const Element& )
 *  -                    Element( const Element&, int ) // with color
 *  -                    Element( const id_type& )
 *  - Unique ID: id_type operator()() const
 *  - Comparator: bool operator<( const Element& ) const
 *  -             bool operator==( const Element& ) const
 *  - Coloring: int color() const
 */

template< class Element >
class space
{
  friend class portal;

  typedef std::map< Element, int > Container; // Element x Index
  typedef std::map< Element, Element > Margin;

  Container xi_; // Element -> Index
  Margin i2e_;   // Internal Element -> External Element, Not Governed Element!
  Margin e2i_;   // External Element -> External Element, Not Governed Element!

  space( const Container& xi, const Margin& i2e, const Margin& e2i )
    : xi_( xi ), i2e_( i2e ), e2i_( e2i )
  {}

#ifdef ELAI_USE_MPI
  void marshalize( const portal& port ) const
  {
    int m = xi_.size(), n = i2e_.size(), off, *ind;
    typename Element::id_type *id;
    typename Element::id_type *iid;
    typename Element::id_type *eid;

    port.send( &m, mpi_< int >().type );
    port.send( &n, mpi_< int >().type );

    ind = new int[ m ];
    id = new typename Element::id_type[ m ];
    off = 0;
    for ( const_iterator it = xi_.begin(); it != xi_.end(); ++it )
    {
      id[ off ] = const_point( it ).element();
      ind[ off ] = const_point( it ).index;
      ++off;
    }
    port.send( id, mpi_< typename Element::id_type >().type, m );
    port.send( ind, mpi_< int >().type, m );

    iid = new typename Element::id_type[ n ];
    eid = new typename Element::id_type[ n ];
    off = 0;
    for ( const_marginal_iterator it = i2e_.begin(); it != i2e_.end(); ++it )
    {
      iid[ off ] = const_internal_point( it ).internal();
      eid[ off ] = const_internal_point( it ).external();
      ++off;
    }
    port.send( iid, mpi_< typename Element::id_type >().type, n );
    port.send( eid, mpi_< typename Element::id_type >().type, n );

    delete [] eid;
    delete [] iid;
    delete [] id;
    delete [] ind;
  }
#endif

public:
  typedef typename Container::value_type value_type;
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;

  struct point
  {
    point( value_type& v ) : element( v.first ), index( v.second ) {}
    point( iterator& it ) : element( it->first ), index( it->second ) {}
    const Element& element;
    int& index;
  };
  struct const_point
  {
    const_point( const value_type& v ) : element( v.first ), index( v.second ) {}
    const_point( const_iterator& it ) : element( it->first ), index( it->second ) {}
    const Element& element;
    const int& index;
  };

  typedef typename Margin::value_type marginal_value_type;
  typedef typename Margin::iterator marginal_iterator;
  typedef typename Margin::const_iterator const_marginal_iterator;

  struct internal_point
  {
    internal_point( marginal_value_type& v )
      : internal( v.first ), external( v.second ) {}
    internal_point( marginal_iterator& it )
      : internal( it->first ), external( it->second ) {}
    const Element& internal;
    Element& external;
  };
  struct const_internal_point
  {
    const_internal_point( const marginal_value_type& v )
      : internal( v.first ), external( v.second ) {}
    const_internal_point( const_marginal_iterator& it )
      : internal( it->first ), external( it->second ) {}
    const Element& internal;
    const Element& external;
  };

  struct external_point
  {
    external_point( marginal_value_type& v )
      : external( v.first ), internal( v.second ) {}
    external_point( marginal_iterator& it )
      : external( it->first ), internal( it->second ) {}
    const Element& external;
    Element& internal;
  };
  struct const_external_point
  {
    const_external_point( const marginal_value_type& v )
      : external( v.first ), internal( v.second ) {}
    const_external_point( const_marginal_iterator& it )
      : external( it->first ), internal( it->second ) {}
    const Element& external;
    const Element& internal;
  };

  space() {}
  space( const space< Element >& src )
    : xi_( src.xi_ ), i2e_( src.i2e_ ), e2i_( src.e2i_ )
  {}
#ifdef ELAI_USE_MPI
  space( portal& port )
  {
    int m, n, *ind;
    typename Element::id_type *id;
    typename Element::id_type *iid;
    typename Element::id_type *eid;

    port.recv( &m, mpi_< int >().type );
    port.recv( &n, mpi_< int >().type );

    ind = new int[ m ];
    id = new typename Element::id_type[ m ];
    port.recv( id, mpi_< typename Element::id_type >().type, m );
    port.recv( ind, mpi_< int >().type, m );
    for ( int i = 0; i < m; ++i )
      xi_.insert( value_type( Element( id[ i ] ), ind[ i ] ) );

    iid = new typename Element::id_type[ n ];
    eid = new typename Element::id_type[ n ];
    port.recv( iid, mpi_< typename Element::id_type >().type, n );
    port.recv( eid, mpi_< typename Element::id_type >().type, n );

    for ( int i = 0; i < n; ++i ) link( Element( iid[ i ] ), Element( eid[ i ] ) );

    delete [] eid;
    delete [] iid;
    delete [] id;
    delete [] ind;
  }
#endif
  ~space() {}

  space< Element >& operator=( const space< Element >& src )
  {
    xi_.clear(); xi_ = src.xi_;
    i2e_.clear(); i2e_ = src.i2e_;
    e2i_.clear(); e2i_ = src.e2i_;

    return *this;
  }

  int size() const { return xi_.size(); }
  int marginal_size() const { assert( i2e_.size() == e2i_.size() ); return i2e_.size(); }

  space< Element >& join( const Element& x )
  {
    int i = xi_.size();

    assert( xi_.find( x ) == xi_.end() );
    xi_.insert( value_type( x, i ) );

    return *this;
  }
  space< Element >& join( const Element& x, const Element& y )
  {
    int i = xi_.size();

    assert( xi_.find( x ) == xi_.end() );
    xi_.insert( value_type( x, i ) );
    i2e_.insert( marginal_value_type( x, y ) );
    e2i_.insert( marginal_value_type( y, x ) );
    assert( i2e_.size() == e2i_.size() );

    return *this;
  }

  space< Element >& link( const Element& x, const Element& y )
  {
    assert( xi_.find( x ) != xi_.end() );
    i2e_.insert( marginal_value_type( x, y ) );
    e2i_.insert( marginal_value_type( y, x ) );

    return *this;
  }

  // The find method requires operator== of Element.
  // intarnal_contain or external_contain mean internal-side or external-side of marginals.
  // If you want "internaly" contain, use govern.
  bool contain( const Element& x ) const { return xi_.find( x ) != xi_.end(); }
  bool govern( const Element& x ) const { return contain( x ) && !internal_contain( x ); }
  bool internal_contain( const Element& x ) const { return i2e_.find( x ) != i2e_.end(); }
  bool external_contain( const Element& y ) const { return e2i_.find( y ) != e2i_.end(); }

  int& index( const Element& x )
  {
    iterator it = xi_.find( x );
    assert( it != xi_.end() );

    return point( it ).index;
  }
  const int& index( const Element& x ) const
  {
    const_iterator it = xi_.find( x );
    assert( it != xi_.end() );

    return const_point( it ).index;
  }

  Element& external_element( const Element& x )
  {
    marginal_iterator it = i2e_.find( x );
    assert( it != i2e_.end() );

    return internal_point( it ).external;
  }
  const Element& external_element( const Element& x ) const
  {
    const_marginal_iterator it = i2e_.find( x );
    assert( it != i2e_.end() );

    return const_internal_point( it ).external;
  }

  Element& internal_element( const Element& y )
  {
    marginal_iterator it = e2i_.find( y );
    assert( it != e2i_.end() );

    return external_point( it ).internal;
  }
  const Element& internal_element( const Element& y ) const
  {
    const_marginal_iterator it = e2i_.find( y );
    assert( it != e2i_.end() );

    return const_external_point( it ).internal;
  }

  int& external_index( const Element& y )
  {
    marginal_iterator it = e2i_.find( y );
    assert( it != e2i_.end() );

    return index( external_point( it ).internal );
  }
  const int& external_index( const Element& y ) const
  {
    const_marginal_iterator it = e2i_.find( y );
    assert( it != e2i_.end() );

    return index( const_external_point( it ).internal );
  }

  iterator begin() { return xi_.begin(); }
  const_iterator begin() const { return xi_.begin(); }
  const_iterator end() const { return xi_.end(); }

  marginal_iterator internal_begin() { return i2e_.begin(); }
  const_marginal_iterator internal_begin() const { return i2e_.begin(); }
  const_marginal_iterator internal_end() const { return i2e_.end(); }

  marginal_iterator external_begin() { return e2i_.begin(); }
  const_marginal_iterator external_begin() const { return e2i_.begin(); }
  const_marginal_iterator external_end() const { return e2i_.end(); }

  space< Element > operator|( const space< Element >& rhs ) const
  {
    int i = 0;
    Container s;
    Margin i2e, e2i;

    {
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && rhs.govern( external_element( x ) ) ) { ++il; continue; }
        else if ( !rhs.govern( y ) && govern( rhs.external_element( y ) ) ) { ++ir; continue; }
        else if ( !govern( x ) && !rhs.govern( y ) && ( external_element( x ) == rhs.external_element( y ) ) )
        {
          const Element& xx( external_element( x ) );

          s.insert( value_type( x, i ) );
          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
          ++i; ++il; ++ir;
          continue;
        }
        if ( x < y )
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il;
        }
        else if ( y < x )
        {
          s.insert( value_type( y, i ) );
          if ( rhs.internal_contain( y ) )
          {
            const Element& yy( rhs.external_element( y ) );

            i2e.insert( marginal_value_type( y, yy ) );
            e2i.insert( marginal_value_type( yy, y ) );
          }
          ++i; ++ir;
        }
        else
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il; ++ir;
        }
      }
      while ( il != end() )
      {
        const Element& x = const_point( il ).element;

        if ( !govern( x ) )
        {
          const Element& xx( external_element( x ) );

          if ( rhs.contain( xx ) || rhs.external_contain( xx ) ) { ++il; continue; }
        }
        s.insert( value_type( x, i ) );
        if ( internal_contain( x ) )
        {
          const Element& xx( external_element( x ) );

          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
        }
        ++i; ++il;
      }
      while ( ir != rhs.end() )
      {
        const Element& y = const_point( ir ).element;

        if ( !rhs.govern( y ) )
        {
          const Element& yy( rhs.external_element( y ) );

          if ( contain( yy ) || external_contain( yy ) ) { ++ir; continue; }
        }
        s.insert( value_type( y, i ) );
        if ( rhs.internal_contain( y ) )
        {
          const Element& yy( rhs.external_element( y ) );

          i2e.insert( marginal_value_type( y, yy ) );
          e2i.insert( marginal_value_type( yy, y ) );
        }
        ++i; ++ir;
      }
    }

    assert( i2e.size() == e2i.size() );
    return space< Element >( s, i2e, e2i );
  }

  space< Element > operator&( const space< Element >& rhs ) const
  {
    int i = 0;
    Container s;
    Margin i2e, e2i;

    { // xi_ part
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && rhs.govern( external_element( x ) ) ) { ++il; continue; }
        else if ( !rhs.govern( y ) && govern( rhs.external_element( y ) ) ) { ++ir; continue; }
        else if ( !govern( x ) && !rhs.govern( y ) && ( external_element( x ) == rhs.external_element( y ) ) )
        {
          const Element& xx( external_element( x ) );

          s.insert( value_type( x, i ) );
          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
          ++i; ++il; ++ir;
          continue;
        }
        if ( x < y ) ++il;
        else if ( y < x ) ++ir;
        else
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il; ++ir;
        }
      }
    }

    assert( i2e.size() == e2i.size() );
    return space< Element >( s, i2e, e2i );
  }

  space< Element > operator/( const space< Element >& rhs ) const
  {
    int i = 0;
    Container s;
    Margin i2e, e2i;

    {
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && rhs.govern( external_element( x ) ) ) { ++il; continue; }
        else if ( !rhs.govern( y ) && govern( rhs.external_element( y ) ) ) { ++ir; continue; }
        else if ( !govern( x ) && !rhs.govern( y ) && ( external_element( x ) == rhs.external_element( y ) ) )
        { ++il; ++ir; continue; }
        if ( x < y )
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il;
        }
        else if ( y < x ) ++ir;
        else { ++il; ++ir; }
      }
      while ( il != end() )
      {
        const Element& x = const_point( il ).element;

        if ( !govern( x ) )
        {
          const Element& xx( external_element( x ) );

          if ( rhs.contain( xx ) || rhs.external_contain( xx ) ) { ++il; continue; }
        }
        s.insert( value_type( x, i ) );
        if ( internal_contain( x ) )
        {
          const Element& xx( external_element( x ) );

          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
        }
        ++i; ++il;
      }
    }

    assert( i2e.size() == e2i.size() );
    return space< Element >( s, i2e, e2i );
  }

  space< Element >& operator|=( const space< Element >& rhs )
  {
    int i = 0;
    Container s;
    Margin i2e, e2i;

    {
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && rhs.govern( external_element( x ) ) ) { ++il; continue; }
        else if ( !rhs.govern( y ) && govern( rhs.external_element( y ) ) ) { ++ir; continue; }
        else if ( !govern( x ) && !rhs.govern( y ) &&  ( external_element( x ) == rhs.external_element( y ) ) )
        {
            const Element& xx( external_element( x ) );

            s.insert( value_type( x, i ) );
            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
            ++i; ++il; ++ir;
            continue;
        }
        if ( x < y )
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il;
        }
        else if ( y < x )
        {
          s.insert( value_type( y, i ) );
          if ( rhs.internal_contain( y ) )
          {
            const Element& yy( rhs.external_element( y ) );

            i2e.insert( marginal_value_type( y, yy ) );
            e2i.insert( marginal_value_type( yy, y ) );
          }
          ++i; ++ir;
        }
        else
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          { // Is it impossible?
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il; ++ir;
        }
      }
      while ( il != end() )
      {
        const Element& x = const_point( il ).element;

        if ( !govern( x ) )
        {
          const Element& xx( external_element( x ) );

          if ( rhs.contain( xx ) || rhs.external_contain( xx ) ) { ++il; continue; }
        }
        s.insert( value_type( x, i ) );
        if ( internal_contain( x ) )
        {
          const Element& xx( external_element( x ) );

          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
        }
        ++i; ++il;
      }
      while ( ir != rhs.end() )
      {
        const Element& y = const_point( ir ).element;

        if ( !rhs.govern( y ) )
        {
          const Element& yy( rhs.external_element( y ) );

          if ( contain( yy ) || external_contain( yy ) ) { ++ir; continue; }
        }
        s.insert( value_type( y, i ) );
        if ( rhs.internal_contain( y ) )
        {
          const Element& yy( rhs.external_element( y ) );

          i2e.insert( marginal_value_type( y, yy ) );
          e2i.insert( marginal_value_type( yy, y ) );
        }
        ++i; ++ir;
      }
    }

    assert( i2e.size() == e2i.size() );
    xi_.swap( s );
    i2e_.swap( i2e );
    e2i_.swap( e2i );

    return *this;
  }

  space< Element >& operator&=( const space< Element >& rhs )
  {
    int i = 0;
    Container s;
    Margin i2e, e2i;

    { // xi_ part
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && rhs.govern( external_element( x ) ) ) { ++il; continue; }
        else if ( !rhs.govern( y ) && govern( rhs.external_element( y ) ) ) { ++ir; continue; }
        else if ( !govern( x ) && !rhs.govern( y ) && ( external_element( x ) == rhs.external_element( y ) ) )
        {
          const Element& xx( external_element( x ) );

          s.insert( value_type( x, i ) );
          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
          ++i; ++il; ++ir;
          continue;
        }
        if ( x < y ) ++il;
        else if ( y < x ) ++ir;
        else
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il; ++ir;
        }
      }
    }

    assert( i2e.size() == e2i.size() );
    xi_.swap( s );
    i2e_.swap( i2e );
    e2i_.swap( e2i );

    return *this;
  }

  space< Element > operator/=( const space< Element >& rhs )
  {
    int i = 0;
    Container s;
    Margin i2e, e2i;

    {
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && rhs.govern( external_element( x ) ) ) { ++il; continue; }
        else if ( !rhs.govern( y ) && govern( rhs.external_element( y ) ) ) { ++ir; continue; }
        else if ( !govern( x ) && !rhs.govern( y ) && ( external_element( x ) == rhs.external_element( y ) ) )
        { ++il; ++ir; continue; }
        if ( x < y )
        {
          s.insert( value_type( x, i ) );
          if ( internal_contain( x ) )
          {
            const Element& xx( external_element( x ) );

            i2e.insert( marginal_value_type( x, xx ) );
            e2i.insert( marginal_value_type( xx, x ) );
          }
          ++i; ++il;
        }
        else if ( y < x ) ++ir;
        else { ++il; ++ir; }
      }
      while ( il != end() )
      {
        const Element& x = const_point( il ).element;

        if ( !govern( x ) )
        {
          const Element& xx( external_element( x ) );

          if ( rhs.contain( xx ) || rhs.external_contain( xx ) ) { ++il; continue; }
        }
        s.insert( value_type( x, i ) );
        if ( internal_contain( x ) )
        {
          const Element& xx( external_element( x ) );

          i2e.insert( marginal_value_type( x, xx ) );
          e2i.insert( marginal_value_type( xx, x ) );
        }
        ++i; ++il;
      }
    }

    assert( i2e.size() == e2i.size() );
    xi_.swap( s );
    i2e_.swap( i2e );
    e2i_.swap( e2i );

    return *this;
  }

  bool operator==( const space< Element >& rhs ) const
  {
    { // xi_ part
      const_iterator il = begin();
      const_iterator ir = rhs.begin();

      while ( il != end() && ir != rhs.end() )
      {
        const Element& x = const_point( il ).element;
        const Element& y = const_point( ir ).element;

        if ( !govern( x ) && !rhs.govern( y ) && ( external_element( x ) == rhs.external_element( y ) ) )
        {
          ++il; ++ir;
          continue;
        }
        if ( !( x == y ) ) return false;
        else ++il; ++ir;
      }
      if ( il != end() ) return false;
      if ( ir != rhs.end() ) return false;
    }

    return true;
  }
};

}

#endif//__ELAI_SPACE__

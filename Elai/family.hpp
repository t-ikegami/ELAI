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
#ifndef __ELAI_FAMILY__
#define __ELAI_FAMILY__

#include <map>
#include <set>
#include "def.hpp"
#include "space.hpp"
#include "portal.hpp"

namespace elai
{

/*
 * Neighbour Requirements:
 *  - Constructors: Neighbour( const Element& ); This argument is the generating element.
 *  -               Neighbour( const Neighbour& ); This constructor needs the below helper.
 *  -               Neighbour& operator=( const Neighbour& ); // copy-constructor
 *  - Generating element: const Element& element(); Neighbour knows it's generating element.
 *  - Comparators: bool operator()( const Element& ) const; Is an Element contained?
 *  - Add an adjacent point: Neighbour& join( const Element& )
 *  - Remove an adjacent point: Neighbour& erase( const Element& )
 *  -                           Neighbour& erase( const_iterator& )
 *  - Iterators: iterator begin(); For adjacented elements
 *  -            const_iterator begin() const
 *  -            const_iterator end() const
 *  - Semantics: typedef Element *iterator
 *  -            typedef const Element *const_iterator
 *  -            iterators must contain the generating element with the except an opened point.
 */

template< class Element, class Neighbour >
class family
{
public:
  typedef space< Element > Space;
  typedef struct Space::const_point Point;

private:
  friend class portal;

  typedef std::map< Element, Neighbour > Family; // pair( generating-point, it's neighbour )

  // Be aware that a generating point might be not an actual point.
  // So you have to seek points in neighbourhoods.
  Family tau_;

  family( const Family& tau ) : tau_( tau ) {}

#ifdef ELAI_USE_MPI
  void marshalize( const portal& port ) const
  {
    typedef typename Element::id_type id_type;

    int m, nnz, *cnt;
    id_type *base_id;
    id_type *adj_id;

    m = tau_.size();
    port.send( &m, mpi_< int >().type );

    cnt = new int[ m ];
    base_id = new id_type[ m ];
    nnz = 0;
    for ( typename Family::const_iterator it = tau_.begin(); it != tau_.end(); ++it )
    {
      const Element& x = it->first;
      const Space& s = operator()( x );

      base_id[ nnz ] = x();
      cnt[ nnz ] = s.size() - 1;
      ++nnz;
    }
    nnz = 0;
    for ( int i = 0; i < m; ++i ) nnz += cnt[ i ];
    port.send( &nnz, mpi_< int >().type );
    port.send( base_id, mpi_< id_type >().type, m );
    port.send( cnt, mpi_< int >().type, m );

    adj_id = new id_type[ nnz ];
    nnz = 0;
    for ( typename Family::const_iterator it = tau_.begin(); it != tau_.end(); ++it )
    {
      const Element& x = it->first;
      const Space& s = operator()( x );

      for ( typename Space::const_iterator jt = s.begin(); jt != s.end(); ++jt )
      {
        const Element& y = typename Space::const_point( jt ).element;

        if ( x == y ) continue;
        adj_id[ nnz++ ] = y();
      }
    }
    port.send( adj_id, mpi_< id_type >().type, nnz );

    delete [] adj_id;
    delete [] base_id;
    delete [] cnt;
  }
#endif

public:
  family() {}
  family( const Space& s )
  {
    for ( typename Space::const_iterator it = s.begin()
        ; it != s.end(); ++it
        ) tau_.insert( typename Family::value_type( Point( it ).element, Neighbour( Point( it ).element ) ) );
    // This constructor needs Neighbour( const Element& ).
  }
  family( const Neighbour& u )
  {
    tau_.insert( typename Family::value_type( u.element(), u ) );
  }
  family( const family< Element, Neighbour >& src ) : tau_( src.tau_ ) {}
#ifdef ELAI_USE_MPI
  family( portal& port )
  {
    typedef typename Element::id_type id_type;

    int m, nnz, *cnt;
    id_type *base_id;
    id_type *adj_id;

    port.recv( &m, mpi_< int >().type );
    port.recv( &nnz, mpi_< int >().type );

    cnt = new int[ m ];
    base_id = new id_type[ m ];
    port.recv( base_id, mpi_< id_type >().type, m );
    port.recv( cnt, mpi_< int >().type, m );

    adj_id = new id_type[ nnz ];
    port.recv( adj_id, mpi_< id_type >().type, nnz );

    for ( int i = 0, off = 0; i < m; ++i )
    {
      const Element x( base_id[ i ] );
      Neighbour u( x );

      for ( int j = 0; j < cnt[ i ]; ++j )
      {
        const Element y( adj_id[ off++ ] );

        u.join( y );
      }
      join( u );
    }

    delete [] adj_id;
    delete [] base_id;
    delete [] cnt;
  }
#endif
  ~family() {}

  family< Element, Neighbour >& operator=( const family< Element, Neighbour >& src )
  {
    tau_ = src.tau_;
    return *this;
  }

  family< Element, Neighbour >& join( const Neighbour& u )
  {
    typename Family::iterator it = tau_.find( u.element() );

    if ( it == tau_.end() ) tau_.insert( typename Family::value_type( u.element(), u ) );
    //else it->second = u; // This code needs Neighbour::operator=( const Neighbour& )
    else
    {
      for ( typename Neighbour::const_iterator jt = u.begin(); jt != u.end(); ++jt )
        it->second.join( *jt );
    }

    return *this;
  }

  // Exchange the element and adjacent elements from x to y.
  family< Element, Neighbour >& flip( const Element& x, const Element& y, bool symmetric = false )
  {
    if ( symmetric )
    {
      typename Family::iterator it = tau_.find( x );
      Space cand;

      if ( it != tau_.end() )
      {
        const Neighbour& src = it->second;
        Neighbour tmp( y );

        for ( typename Neighbour::const_iterator p = src.begin(); p != src.end(); ++p )
        {
          if ( *p == x ) continue;
          tmp.join( *p );
          cand.join( *p );
        }
        tau_.erase( it );
        tau_.insert( typename Family::value_type( y, tmp ) );
      }

      for ( typename Space::const_iterator p = cand.begin(); p != cand.end(); ++p )
      {
        const Element& pt = typename Space::const_point( p ).element;
        typename Family::iterator it = tau_.find( pt );

        if ( it != tau_.end() )
        {
          Neighbour& u = it->second;

          for ( typename Neighbour::iterator jt = u.begin(); jt != u.end(); ++jt )
          {
            if ( *jt == x )
            {
              u.erase( jt );
              u.join( y );
              break;
            }
          }
        }
      }
    }
    else
    {
      typename Family::iterator it = tau_.find( x );

      if ( it != tau_.end() )
      {
        const Neighbour& src = it->second;
        Neighbour tmp( y );

        for ( typename Neighbour::const_iterator p = src.begin(); p != src.end(); ++p )
        {
          if ( *p == x ) continue;
          tmp.join( *p );
        }
        tau_.erase( it );
        tau_.insert( typename Family::value_type( y, tmp ) );
      }

      for ( typename Family::iterator it = tau_.begin(); it != tau_.end(); ++it )
      {
        Neighbour& u( it->second );

        for ( typename Neighbour::iterator p = u.begin(); p != u.end(); ++p )
        {
          if ( *p == x )
          {
            u.erase( p );
            u.join( y );
            break;
          }
        }
      }
    }

    return *this;
  }

  family< Element, Neighbour > operator|( const family< Element, Neighbour >& rhs ) const
  {
    Family tau = tau_;

    for ( typename Family::const_iterator m2 = rhs.tau_.begin(); m2 != rhs.tau_.end(); ++m2 )
    {
      const Element& y = m2->first;
      typename Family::iterator m1 = tau.find( y );
      const Neighbour& u( m2->second );

      if ( m1 == tau.end() ) // unkown element
      {
        tau.insert( typename Family::value_type( y, u ) );
      }
      else // kwown element
      {
        Neighbour& v = m1->second;

        // append elements in rhs's internal to this internal.
        for ( typename Neighbour::const_iterator p = u.begin(); p != u.end(); ++p )
          v.join( *p );
      }
    }

    return family< Element, Neighbour >( tau );
  }

  family< Element, Neighbour >& operator|=( const family< Element, Neighbour >& rhs )
  {
    Family tau = tau_;

    for ( typename Family::const_iterator m2 = rhs.tau_.begin(); m2 != rhs.tau_.end(); ++m2 )
    {
      const Element& y = m2->first;
      typename Family::iterator m1 = tau.find( y );
      const Neighbour& u( m2->second );

      if ( m1 == tau.end() ) // unkown element
      {
        tau.insert( typename Family::value_type( y, u ) );
      }
      else // kwown element
      {
        Neighbour& v = m1->second;

        // append elements in rhs's internal to this internal.
        for ( typename Neighbour::const_iterator p = u.begin(); p != u.end(); ++p )
          v.join( *p );
      }
    }
    std::swap( tau_, tau );

    return *this;
  }

  family< Element, Neighbour > localize( const Space& s ) const
  {
    Family tau;

    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      const Element& e = Point( it ).element;

      if ( tau_.find( e ) != tau_.end() )
      {
        const Neighbour& neigh = tau_.find( e )->second;

        tau.insert( typename Family::value_type( e, neigh ) );
      }
    }

    return family< Element, Neighbour >( tau );
  }

  family< Element, Neighbour > localize( const Space& s, const Space& t ) const
  {
    Family tau;

    for ( typename Space::const_iterator it = s.begin(); it != s.end(); ++it )
    {
      const Element& e = Point( it ).element;

      if ( tau_.find( e ) != tau_.end() )
      {
        const Neighbour& neigh0 = tau_.find( e )->second;
        Neighbour neigh( e );

        for ( typename Neighbour::const_iterator it = neigh0.begin()
            ; it != neigh0.end(); ++it
            ) if ( t.contain( *it ) ) neigh.join( *it );
        tau.insert( typename Family::value_type( e, neigh ) );
      }
    }

    return family< Element, Neighbour >( tau );
  }

  // A included-point-set in this family.
  // Iteration via points in Neighbour tau_.
  // This method expects that neighbourhoods < s.
  Space operator()( const Space& s ) const
  {
    Space sub_s;

    for ( typename Space::const_iterator p = s.begin(); p != s.end(); ++p )
    {
      const Element& e = Point( p ).element;

      //if ( tau_.find( e ) != tau_.end() ) sub_s.join( e );
      for ( typename Family::const_iterator it = tau_.begin(); it != tau_.end(); ++it )
      {
        const Neighbour& neigh = it->second;

        if ( neigh( e ) ) { sub_s.join( e ); break; }
      }
    }

    return sub_s;
  }

  // A k-level adjacent-point-set in this family from the point x.
  // Iteration via points in Neighbour tau_.
  // At first, searching the (maybe opened) point exactly matches the argument x.
  // Then, enumerating adjacent points of the above point.
  Space operator()( const Element& x, int k = 0 ) const
  {
    Space sub_s;
    std::set< Element > passed;
/*
    typename Family::const_iterator it = tau_.begin();

    // If x were a generating point, tau_.find( x ) returns tau_.end() always.
    // We have to seek the matching element over the tau_.
    while ( it != tau_.end() )
    {
      if ( Neighbour( it->first )( x ) ) break;
      ++it;
    }
*/
    typename Family::const_iterator it = tau_.find( x );

    if ( it != tau_.end() ) // On-hit.
    {
      // From the mathing element ( maybe a generating point ), we list adjacents up incrementaly.
      do
      {
        const Element& chk = it->first;
        const Neighbour& neigh = it->second;

        // Iterators of Neighbour are always pointing an actual element.
        // We can include them all.
        for ( typename Neighbour::const_iterator p = neigh.begin(); p != neigh.end(); ++p )
          if ( !sub_s.contain( *p ) ) sub_s.join( *p );

        passed.insert( chk );

        // Seeking the candidate which we start to list adjacents up from.
        for ( typename Space::const_iterator jt = sub_s.begin(); jt != sub_s.end(); ++jt )
        {
          const Element& e = Point( jt ).element;

          if ( passed.find( e ) == passed.end() )
          {
            it = tau_.find( e ); // Next actual element ( k-chain )
            break;
          }
        }
        --k;
      }
      while( 0 < k && it != tau_.end() );
    }

/*
    // Non sense: typename Family::const_iterator it = tau_.find( x );
    for ( typename Family::const_iterator it = tau_.begin(); it != tau_.end(); ++it )
    {
      // Is x in Neibour of an (opened or actual) point?
      if ( Neighbour( it->first )( x ) ) // needs bool Neighbour::operator()( const Element& )
      {
        int k = lv;
        do
        {
          const Neighbour& neigh = it->second;

          // Iterators of Neighbour are always pointing an actual point.
          for ( typename Neighbour::const_iterator p = neigh.begin(); p != neigh.end(); ++p ) sub_s.join( *p );

          // Only actual points are to be checked.
          // First time, it->first == x, if it->first is actual.
          if ( sub_s.contain( it->first ) ) passed.insert( it->first );

          for ( typename Space::const_iterator jt = sub_s.begin(); jt != sub_s.end(); ++jt )
          {
            if ( passed.find( Point( jt ).element ) == passed.end() )
            {
              it = tau_.find( Point( jt ).element ); // Next actual point (k-chain)
              break;
            }
          }
          --k;
        }
        while ( ( 0 < k ) && ( it != tau_.end() ) );
        break;
      }
    }
*/

    return sub_s;
  }
};

}

#endif//__ELAI_FAMILY__

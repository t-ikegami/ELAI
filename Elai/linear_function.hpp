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
#ifndef __ELAI_LINEAR_FUNCTION__
#define __ELAI_LINEAR_FUNCTION__

#include "def.hpp"
#include "coherence.hpp"
#include "sync.hpp"
#include "portal.hpp"
#include "space.hpp"
#include "family.hpp"
#include "subjugator.hpp"
#include "vector.hpp"

namespace elai
{

template< class Element, class Neighbour, class Range >
class linear_function
{
  friend class coherence;
  friend class sync;
  friend class portal;

public:
  typedef Range range;
  typedef vector< range > Vector;
  typedef space< Element > Space;
  typedef family< Element, Neighbour > Family;

private:
  typedef struct Space::const_point Point;

#ifdef ELAI_USE_MPI
  typedef struct Space::const_internal_point Marginal;
#endif

  int m_;
  const Space x_;
  Vector f_;

#ifdef ELAI_USE_MPI
  void coherence_setup( coherence& coherent )
  {
    MPI_Comm comm = coherent.comm();
    int myrank = coherent.myself();
    int comm_size = 0;

    MPI_Comm_size( comm, &comm_size );
    int MNN[ comm_size ][ comm_size ];

    for ( int i = 0; i < comm_size; ++i ) MNN[ myrank ][ i ] = 0;
    for ( typename Space::const_marginal_iterator it = x_.internal_begin()
        ; it != x_.internal_end(); ++it
        ) MNN[ myrank ][ Marginal( it ).external.color() ] += 1;
    MPI_Allgather
      ( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MNN, comm_size, MPI_INT, comm );

    // Ready for RECV
    int recv_num = 0, recv_rank_num = 0;
    for ( int i = 0; i < comm_size; ++i )
    {
      if ( MNN[ myrank ][ i ] <= 0 ) continue;
      ++recv_rank_num;
      recv_num += MNN[ myrank ][ i ];
    }
    int *recv_rank = new int[ recv_rank_num ];
    int *recv_id_num = new int[ recv_rank_num ];
    MPI_Request *recv_req = new MPI_Request[ recv_rank_num ];
    typename Element::id_type *recv_id = new typename Element::id_type[ recv_num ];
    for ( int i = 0, j = 0, k = 0, off = 0; i < comm_size; ++i )
    {
      if ( MNN[ myrank ][ i ] <= 0 ) continue;
      recv_rank[ j ] = i;
      recv_id_num[ j ] = MNN[ myrank ][ i ];
      for ( typename Space::const_marginal_iterator it = x_.internal_begin()
          ; it != x_.internal_end(); ++it
          )
      {
        if ( i != Marginal( it ).external.color() ) continue;
        recv_id[ off + k ] = Marginal( it ).external();
        ++k;
      }
      off += MNN[ myrank ][ i ];
      k = 0;
      ++j;
    }

    // Ready for SEND
    int send_num = 0, send_rank_num = 0;
    for ( int i = 0; i < comm_size; ++i )
    {
      if ( MNN[ i ][ myrank ] <= 0 ) continue;
      ++send_rank_num;
      send_num += MNN[ i ][ myrank ];
    }
    int *send_rank = new int[ send_rank_num ];
    int *send_id_num = new int[ send_rank_num ];
    MPI_Request *send_req = new MPI_Request[ send_rank_num ];
    typename Element::id_type *send_id = new typename Element::id_type[ send_num ];
    for ( int i = 0, j = 0; i < comm_size; ++i )
    {
      if ( MNN[ i ][ myrank ] <= 0 ) continue;
      send_rank[ j ] = i;
      send_id_num[ j ] = MNN[ i ][ myrank ];
      ++j;
    }

    // Send/Recv Element's id
    for ( int i = 0, off = 0; i < send_rank_num; ++i )
    {
      MPI_Irecv
        ( &send_id[ off ]
        , send_id_num[ i ]
        , mpi_< typename Element::id_type >().type
        , send_rank[ i ]
        , myrank
        , comm
        , &send_req[ i ]
        );
      off += send_id_num[ i ];
    }
    for ( int i = 0, off = 0; i < recv_rank_num; ++i )
    {
      MPI_Send
        ( &recv_id[ off ]
        , recv_id_num[ i ]
        , mpi_< typename Element::id_type >().type
        , recv_rank[ i ]
        , recv_rank[ i ]
        , comm
        );
      off += recv_id_num[ i ];
    }
    if ( 0 < send_rank_num ) MPI_Waitall( send_rank_num, send_req, MPI_STATUSES_IGNORE );

    MPI_Aint base;
    MPI_Get_address( f_.val(), &base );

    // Commit Recv Datatype
    MPI_Datatype *recv = new MPI_Datatype[ recv_rank_num ];
    for ( int i = 0, off = 0; i < recv_rank_num; ++i )
    {
      int n = recv_id_num[ i ];
      int *len = new int[ n ];
      MPI_Aint *addr = new MPI_Aint[ n ];
      MPI_Datatype *type = new MPI_Datatype[ n ];

      for ( int j = 0; j < n; ++j )
      {
        int idx = x_.external_index( Element( recv_id[ off++ ] ) );

        coherent.exclude( idx );
        len[ j ] = 1;
        MPI_Get_address( &f_( idx ), &addr[ j ] );
        //  ( &f_( x_.external_index( Element( recv_id[ off++ ] ) ) ), &addr[ j ] );
        addr[ j ] -= base;
        type[ j ] = mpi_< range >().type;
      }
      MPI_Type_create_struct( n, len, addr, type, &recv[ i ] );
      MPI_Type_commit( &recv[ i ] );
      coherent.read( recv[ i ], recv_rank[ i ] );

      delete [] type;
      delete [] addr;
      delete [] len;
    }

    // Commit Send Datatype
    MPI_Datatype *send = new MPI_Datatype[ send_rank_num ];
    for ( int i = 0, off = 0; i < send_rank_num; ++i )
    {
      int n = send_id_num[ i ];
      int *len = new int[ n ];
      MPI_Aint *addr = new MPI_Aint[ n ];
      MPI_Datatype *type = new MPI_Datatype[ n ];

      for ( int j = 0; j < n; ++j )
      {
        len[ j ] = 1;
        MPI_Get_address
          ( &f_( x_.index( Element( send_id[ off++ ] ) ) ), &addr[ j ] );
        addr[ j ] -= base;
        type[ j ] = mpi_< range >().type;
      }
      MPI_Type_create_struct( n, len, addr, type, &send[ i ] );
      MPI_Type_commit( &send[ i ] );
      coherent.write( send[ i ], send_rank[ i ] );

      delete [] type;
      delete [] addr;
      delete [] len;
    }

    delete [] send; delete [] recv;
    delete [] send_id; delete [] recv_id;
    delete [] send_req; delete [] recv_req;
    delete [] send_id_num; delete [] recv_id_num;
    delete [] send_rank; delete [] recv_rank;
  }

  void sync_setup( sync& s )
  {
    s.setup( f_.val(), f_.m(), mpi_< range >().type );
  }

  void marshalize( const portal& port ) const
  {
    port( x_ );
#ifdef ELAI_USE_MPI3
    // portal::send( const void *, MPI_Datatype, int ) const needs MPI-3
    port.send( f_.val(), mpi_< range >().type, f_.m() );
#else
    range f[ f_.m() ];

    for ( int i = 0; i < f_.m(); ++i ) f[ i ] = f_( i );
    port.send( f, mpi_< range >().type, f_.m() );
#endif
  }

  void reinforce( const portal& port ) const
  {
    // Asynchronous updates.
    port.send( f_.val(), mpi_< range >().type, f_.m() );
  }
#endif

public:
  linear_function
    ( const Space& x
    )
    : m_( x.size() ), x_( x ), f_( m_ )
  {}

  linear_function
    ( const Space& x, const range *v
    )
    : m_( x.size() ), x_( x ), f_()
  {
    f_.setup( m_ );
    for ( int i = 0; i < m_; ++i ) f_( i ) = v[ i ];
  }

  linear_function
    ( const Space& x, const Vector& f
    )
    : m_( x.size() ), x_( x ), f_( f )
  {}

  linear_function
    ( const linear_function< Element, Neighbour, Range >& f )
    : m_( f.m_ ), x_( f.x_ ), f_( f.f_ )
  {}

#ifdef ELAI_USE_MPI
  linear_function( portal& port )
    : m_( 0 ), x_( port ), f_()
  {
    m_ = x_.size();
    f_.setup( m_ );
    port.recv( f_.val(), mpi_< range >().type, f_.m() );
  }
#endif

  ~linear_function()
  {}

#ifdef ELAI_USE_MPI
  void operator()( portal& port )
  {
    port.recv( f_.val(), mpi_< range >().type, f_.m() );
  }
#endif

  int dim() const { return m_; }
  int codim() const { return 1; }
  const Space& dom() const { return x_; }
  Vector& ran() { return f_; }

  range& operator()( const Element& x ) { return f_( x_.index( x ) ); }
  const range& operator()( const Element& x ) const { return f_( x_.index( x ) ); }

#ifdef ELAI_USE_PYTHON
  void setVector(const Element& x, range val)
  {
    f_( x_.index( x ) ) = val;
  }
#endif

  linear_function< Element, Neighbour, Range >& clear( const Range v )
  {
    f_.clear( v );

    return *this;
  }
  linear_function< Element, Neighbour, Range >& clear( const Range v, const Space& s )
  {
    for ( typename Space::const_iterator it = s.begin()
        ; it != s.end(); ++it
        ) f_( x_.index( Point( it ).element ) ) = v;

    return *this;
  }

  // 'localize' methods actually mean restriction-maps.
  linear_function< Element, Neighbour, Range > localize( const Space& s ) const
  {
    return localize( Family( s ) );
  }
  linear_function< Element, Neighbour, Range > localize( const Neighbour& u ) const
  {
    return localize( Family( u ) );
  }
  linear_function< Element, Neighbour, Range > localize( const Family& tau ) const
  {
#ifdef ELAI_USE_C11
    Space&& s = std::move( tau( x_ ) );
#else
    Space s = tau( x_ );
#endif
    Vector f( s.size() );

    for ( typename Space::const_iterator it = s.begin()
        ; it != s.end(); ++it
        ) f( Point( it ).index ) = f_( x_.index( Point( it ).element ) );
        //----------------------   -------------------------------------
        //local  function-range  = global function-range

    return linear_function< Element, Neighbour, Range >( s, f );
  }

  linear_function< Element, Neighbour, Range >& reflectIn
    ( const linear_function< Element, Neighbour, Range >& f )
  {
    for ( typename Space::const_iterator it = f.x_.begin(); it != f.x_.end(); ++it )
    {
      const Element& x = Point( it ).element;
      int i = Point( it ).index;
      int I = -1;

      if ( x_.govern( x ) ) I = x_.index( x );
      else if ( x_.external_contain( x ) ) I = x_.external_index( x );

      if ( 0 <= I ) f_( I ) = f.f_( i );
    }

    return *this;
  }

  linear_function< Element, Neighbour, Range >& reflectIn
    ( const linear_function< Element, Neighbour, Range >& f
    , const subjugator< Element, Neighbour >& subj
    )
  {
    for ( typename Space::const_iterator it = f.x_.begin(); it != f.x_.end(); ++it )
    {
      const Element& x = Point( it ).element;
      int i = Point( it ).index;
      int I = -1;
      int c = subj.color( x );
      const Element y( x, c );

      if ( x_.govern( y ) ) I = x_.index( y );
      else if ( x_.external_contain( y ) ) I = x_.external_index( y );

      if ( 0 <= I ) f_( I ) = f.f_( i );
    }

    return *this;
  }

  linear_function< Element, Neighbour, Range >& reflect
    ( const linear_function< Element, Neighbour, Range >& f )
  {
    for ( typename Space::const_iterator it = x_.begin(); it != x_.end(); ++it )
    {
      const Element& x = Point( it ).element;
      int I = Point( it ).index;
      int i = -1;

      if ( f.x_.govern( x ) ) i = f.x_.index( x );
      //else if ( f.x_.external_contain( x ) ) i = f.x_.external_index( x );

      if ( 0 <= i ) f_( I ) = f.f_( i );
    }

    return *this;
  }

  linear_function< Element, Neighbour, Range >& reflect
    ( const linear_function< Element, Neighbour, Range >& f
    , const subjugator< Element, Neighbour >& subj
    )
  {
    for ( typename Space::const_iterator it = x_.begin(); it != x_.end(); ++it )
    {
      const Element& x = Point( it ).element;
      int I = Point( it ).index;
      int i = -1;
      int c = subj.color( x );
      const Element y( x, c );

      if ( f.x_.govern( y ) ) i = f.x_.index( y );
      //else if ( f.x_.external_contain( y ) ) i = f.x_.external_index( y );

      if ( 0 <= i ) f_( I ) = f.f_( i );
    }

    return *this;
  }

  linear_function< Element, Neighbour, Range > extend( const linear_function< Element, Neighbour, Range >& f ) const
  {
#ifdef ELAI_USE_C11
    Space&& s = std::move( x_ | f.x_ );
#else
    Space s = x_ | f.x_;
#endif
    Vector v( s.size() );

    for ( typename Space::const_iterator it = s.begin()
        ; it != s.end(); ++it
        )
    {
      const Element& e = Point( it ).element;
      const int i = Point( it ).index;

      if ( f.x_.contain( e ) ) v( i ) = f( e );
      else v( i ) = f_( x_.index( e ) );
    }

    return linear_function< Element, Neighbour, Range >( s, v );
  }

/*
  // this: o  oxx ooo   xo <- iterate
  //       |   || |||   | 
  //    f: ooo oooooo  oo  
  linear_function< Element, Neighbour, Range >& reflect( const linear_function< Element, Neighbour, Range >& f )
  {
    for ( typename Space::const_iterator it = x_.begin(); it != x_.end(); ++it )
    {
      const Element& x = Point( it ).element;
      const int i = Point( it ).index;

      if ( x_.govern( x ) ) f_( i ) = f( x );
      else f_( i ) = f( x_.external_element( x ) );
    }

    return *this;
  }
*/
};

}

#endif//__ELAI_LINEAR_FUNCTION__

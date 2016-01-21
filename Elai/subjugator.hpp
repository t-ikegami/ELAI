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
#ifndef __ELAI_SUBJUGATOR__
#define __ELAI_SUBJUGATOR__

#include "def.hpp"
#include "space.hpp"
#include "family.hpp"

#ifdef ELAI_USE_MPI
#include "mpi.h"
#endif

namespace elai
{

template< class Element, class Neighbour >
class subjugator
{
  typedef space< Element > Space;
  typedef family< Element, Neighbour > Family;
#ifdef ELAI_USE_MPI
  typedef struct Space::const_internal_point Marginal;
#endif

  const Space& base_;
  const Family& adjacent_;
  const std::vector< int >& palette_;
  int *color_;
#ifdef ELAI_USE_MPI
  MPI_Comm intra_, inter_;
#endif

  void setup()
  {
    if ( palette_.size() <= 0 ) return;
    unsigned int partition = ( base_.size() - base_.marginal_size() ) / palette_.size();
    unsigned int cnt = 0, i = 0;

    color_ = new int[ base_.size() ];
    for ( typename Space::const_iterator it = base_.begin(); it != base_.end(); ++it )
    {
      typename Space::const_point pt( it );

      if ( base_.govern( pt.element ) )
      {
        color_[ pt.index ] = palette_[ i ];
        if ( partition <= ++cnt && ( i + 1 ) < palette_.size() ) { ++i; cnt = 0; }
      }
      else color_[ pt.index ] = pt.element.color();
    }
#ifdef ELAI_USE_MPI
    if ( inter_ != NULL ) marginal_setup();
#endif
  }

#ifdef ELAI_USE_MPI
  void marginal_setup()
  {
    int myrank;
    MPI_Comm_rank( intra_, &myrank ); // LOCAL RANK

    if ( myrank == 0 )
    { // Inter Communicator is available only for leaders.
      int mysize, myrid;
      MPI_Comm_size( inter_, &mysize ); // LEADER SIZE == NUMBER of REGIONS
      MPI_Comm_rank( inter_, &myrid );

      int MNN[ mysize ][ mysize ];

      for ( int i = 0; i < mysize; ++i ) MNN[ myrid ][ i ] = 0;
      for ( typename Space::const_marginal_iterator it = base_.internal_begin()
          ; it != base_.internal_end(); ++it
          ) MNN[ myrid ][ Marginal( it ).external.color() ] += 1;
      MPI_Allgather( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MNN, mysize, MPI_INT, inter_ );

      // FOR SEND EXTERNAL ELEMENTS
      int cnt_send_regs = 0;
      int cnt_send_elms = 0;
      for ( int i = 0; i < mysize; ++i )
      {
        if ( MNN[ myrid ][ i ] <= 0 ) continue;
        cnt_send_regs += 1;
        cnt_send_elms += MNN[ myrid ][ i ];
      }
      int send_regs[ cnt_send_regs ];
      int send_nums[ cnt_send_regs ];
      typename Element::id_type send_elms[ cnt_send_elms ];
      for ( int i = 0, j = 0, k = 0, off = 0; i < mysize; ++i )
      {
        if ( MNN[ myrid ][ i ] <= 0 ) continue;
        send_regs[ j ] = i;
        send_nums[ j ] = MNN[ myrid ][ i ];
        for ( typename Space::const_marginal_iterator it = base_.internal_begin()
            ; it != base_.internal_end(); ++it
            )
        {
          if ( i != Marginal( it ).external.color() ) continue;
          send_elms[ off + k ] = Marginal( it ).external();
          ++k;
        }
        off += MNN[ myrid ][ i ];
        k = 0;
        ++j;
      }

      // FOR RECV INTERNAL ELEMENTS
      int cnt_recv_regs = 0;
      int cnt_recv_elms = 0;
      for ( int j = 0; j < mysize; ++j )
      {
        if ( MNN[ j ][ myrid ] <= 0 ) continue;
        cnt_recv_regs += 1;
        cnt_recv_elms += MNN[ j ][ myrid ];
      }
      int recv_regs[ cnt_recv_regs ];
      int recv_nums[ cnt_recv_regs ];
      typename Element::id_type recv_elms[ cnt_recv_elms ];
      MPI_Request recv_req[ cnt_recv_regs ];
      for ( int i = 0, j = 0; i < mysize; ++i )
      {
        if ( MNN[ i ][ myrid ] <= 0 ) continue;
        recv_regs[ j ] = i;
        recv_nums[ j ] = MNN[ i ][ myrid ];
        ++j;
      }

      // SEND EXTERNAL ELEMENTS & RECV INTERNAL ELEMENTS
      for ( int i = 0, off = 0; i < cnt_recv_regs; ++i )
      {
        MPI_Irecv
          ( &recv_elms[ off ]
          , recv_nums[ i ], mpi_< typename Element::id_type >().type
          , recv_regs[ i ], myrid, inter_, &recv_req[ i ]
          );
        off += recv_nums[ i ];
      }
      for ( int i = 0, off = 0; i < cnt_send_regs; ++i )
      {
        MPI_Send
          ( &send_elms[ off ]
          , send_nums[ i ], mpi_< typename Element::id_type >().type
          , send_regs[ i ], send_regs[ i ], inter_
          );
        off += send_nums[ i ];
      }
      if ( 0 < cnt_recv_regs ) MPI_Waitall( cnt_recv_regs, recv_req, MPI_STATUSES_IGNORE );

      // FOR SEND INTERNAL COLORS
      int send_cols[ cnt_recv_elms ];
      for ( int i = 0, off = 0; i < cnt_recv_regs; ++i )
      {
        int rid = recv_regs[ i ];

        for ( int j = 0; j < recv_nums[ i ]; ++j )
        {
          send_cols[ off ] = color_[ base_.index( Element( recv_elms[ off ] ) ) ];
          ++off;
        }
      }

      // FOR RECV EXTERNAL COLORS
      int recv_cols[ cnt_send_elms ];
      MPI_Request send_req[ cnt_send_elms ];

      // SEND INTERNAL COLORS & RECV EXTERNAL COLORS
      for ( int i = 0, off = 0; i < cnt_send_regs; ++i )
      {
        MPI_Irecv
          ( &recv_cols[ off ]
          , send_nums[ i ], MPI_INT
          , send_regs[ i ], myrid, inter_, &send_req[ i ]
          );
        off += send_nums[ i ];
      }
      for ( int i = 0, off = 0; i < cnt_recv_regs; ++i )
      {
        MPI_Send
          ( &send_cols[ off ]
          , recv_nums[ i ], MPI_INT
          , recv_regs[ i ], recv_regs[ i ], inter_
          );
        off += recv_nums[ i ];
      }
      if ( 0 < cnt_send_regs ) MPI_Waitall( cnt_send_regs, send_req, MPI_STATUSES_IGNORE );

      // UPDATE COLORS
      for ( int i = 0; i < cnt_send_elms; ++i )
        color_[ base_.external_index( Element( send_elms[ i ] ) ) ] = recv_cols[ i ];
    }

    // SYNC COLORS IN THE REGION
    MPI_Barrier( intra_ );
    MPI_Bcast( color_, base_.size(), MPI_INT, 0, intra_ );
  }
#endif

  void terminate()
  {
    if ( color_ != NULL ) { delete [] color_; color_ = NULL; }
  }

public:
  subjugator
    ( const Space& base, const Family& adjacent, const std::vector< int >& palette
#ifdef ELAI_USE_MPI
    , MPI_Comm intra = NULL, MPI_Comm inter = NULL
#endif
    )
    : base_( base ), adjacent_( adjacent ), palette_( palette ), color_( NULL )
#ifdef ELAI_USE_MPI
    , intra_( intra ), inter_( inter )
#endif
  {
    setup();
  }
  ~subjugator()
  {
    terminate();
  }

  int color( const Element& x ) const { return color_[ base_.index( x ) ]; }

  Space operator()( int color )
  {
    Space base, sub, margin;

    for ( typename Space::const_iterator it = base_.begin()
        ; it != base_.end(); ++it
        )
    {
      typename Space::const_point pt( it );

      if ( color_[ pt.index ] == color ) base.join( pt.element );
    }
    for ( typename Space::const_iterator it = base.begin()
        ; it != base.end(); ++it
        )
    {
      typename Space::const_point pt0( it );
#ifdef ELAI_USE_C11
      const Space&& adj = std::move( adjacent_( pt0.element ) );
#else
      const Space adj = adjacent_( pt0.element );
#endif

      sub.join( Element( pt0.element, color ) );
      for ( typename Space::const_iterator jt = adj.begin()
          ; jt != adj.end(); ++jt
          )
      {
        typename Space::const_point pt( jt );

        if ( base.contain( pt.element ) || margin.contain( Element( pt.element, color ) ) ) continue;

        int c = color_[ base_.index( pt.element ) ];

        margin.join( Element( pt.element, color ), Element( pt.element, c ) );
      }
    }

    return sub | margin;
  }

  Family operator()( int color, const Space& sub )
  {
    Space base, margin;
    Family adjs;

    for ( typename Space::const_iterator it = base_.begin()
        ; it != base_.end(); ++it
        )
    {
      typename Space::const_point pt( it );

      if ( color_[ pt.index ] == color ) base.join( pt.element );
    }
    for ( typename Space::const_iterator it = base.begin()
        ; it != base.end(); ++it
        )
    {
      typename Space::const_point pt0( it );
#ifdef ELAI_USE_C11
      const Space&& adj = std::move( adjacent_( pt0.element ) );
#else
      const Space adj = adjacent_( pt0.element );
#endif
      Neighbour neigh( Element( pt0.element, color ) );

      for ( typename Space::const_iterator jt = adj.begin()
          ; jt != adj.end(); ++jt
          )
      {
        typename Space::const_point pt( jt );
        const Element pt1( pt.element, color );

        if ( !base.contain( pt.element ) && !margin.contain( pt.element ) )
          margin.join( pt.element, Element( pt.element, color ) ); // ( BASE, RECOLORED )
        if ( !sub.contain( pt1 ) ) continue;
        neigh.join( pt1 );
      }
      adjs.join( neigh );
    }
    for ( typename Space::const_marginal_iterator it = margin.internal_begin()
        ; it != margin.internal_end(); ++it
        )
    {
      typename Space::const_internal_point pt0( it );
#ifdef ELAI_USE_C11
      const Space&& adj = std::move( adjacent_( pt0.element ) );
#else
      const Space adj = adjacent_( pt0.internal ); // BASE ELEMENT
#endif
      Neighbour neigh( pt0.external ); // RECOLORED ELEMENT

      for ( typename Space::const_iterator jt = adj.begin(); jt != adj.end(); ++jt )
      {
        typename Space::const_point pt( jt );
        const Element pt1( pt.element, color );

        if ( !sub.contain( pt1 ) ) continue;
        neigh.join( pt1 );
      }
      adjs.join( neigh );
    }

    return adjs;
  }
};

}

#endif//__ELAI_SUBJUGATOR__

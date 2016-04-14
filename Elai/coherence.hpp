/*
 *
 * Elastic Linear Algebra Interface (ELAI)
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
#ifndef __ELAI_COHERENCE__
#define __ELAI_COHERENCE__

#include <complex>
#include <vector>
#include <utility>
#include "def.hpp"

#ifdef ELAI_USE_MPI
#include "mpi.h"

namespace elai
{

class coherent_type
{
public:
  coherent_type( MPI_Datatype pair_type, int pair_rank )
    : type( pair_type ), pair( pair_rank )
  {}
  coherent_type( const coherent_type& v )
    : type( v.type ), pair( v.pair )
  {}
  coherent_type( std::pair< MPI_Datatype, int >& v )
    : type( v.first ), pair( v.second )
  {}
  coherent_type( std::pair< MPI_Datatype, int > *it )
    : type( it->first ), pair( it->second )
  {}
  MPI_Datatype type;
  int pair;
};

class coherence
{
  typedef std::vector< coherent_type > Slot;
  typedef std::vector< int > ESlot;

  int myself_;
  MPI_Comm comm_;
  Slot rslot_, wslot_;
  ESlot eslot_;
  MPI_Request *req_;

  template< class Coef >
  Coef fix_prod( const Coef& a, const Coef& b ) const { return a * b; }
  template< class Coef >
  Coef fix_prod( const std::complex< Coef >& a, const std::complex< Coef >& b ) const
  { return a * std::conj( b ); }

public:
  typedef Slot::iterator iterator;
  typedef Slot::const_iterator const_iterator;

  /*
  coherence()
    : myself_( -1 ), comm_( NULL ), rslot_(), wslot_(), req_( NULL )
  {}
  */
  template< class Target >
  coherence( Target& obj, MPI_Comm comm )
    : myself_( -1 ), comm_( comm ), rslot_(), wslot_(), eslot_(), req_( NULL )
  {
    MPI_Comm_rank( comm_, &myself_ );
    obj.coherence_setup( *this );
    if ( 0 < rslot_.size() ) req_ = new MPI_Request[ rslot_.size() ];
  }
  ~coherence()
  {
    if ( req_ != NULL ) { delete [] req_; req_ = NULL; }
  }

  int myself() const { return myself_; }
  MPI_Comm comm() { return comm_; }

  void operator()( void *base )
  {
    if ( comm_ == NULL ) return;
    for ( unsigned int i = 0; i < rslot_.size(); ++i )
    {
      coherent_type& slot( rslot_[ i ] );

      MPI_Irecv( base, 1, slot.type, slot.pair, myself_, comm_, &req_[ i ] );
    }
    for ( unsigned int i = 0; i < wslot_.size(); ++i )
    {
      coherent_type& slot( wslot_[ i ] );

      MPI_Send( base, 1, slot.type, slot.pair, slot.pair, comm_ );
    }
    if ( 0 < rslot_.size() ) MPI_Waitall( rslot_.size(), req_, MPI_STATUSES_IGNORE );
  }

  void read( MPI_Datatype pair_type, int pair_rank )
  { read( coherent_type( pair_type, pair_rank ) ); }
  void read( coherent_type type )
  { rslot_.push_back( type ); }

  void write( MPI_Datatype pair_type, int pair_rank )
  { write( coherent_type( pair_type, pair_rank ) ); }
  void write( coherent_type type )
  { wslot_.push_back( type ); }

  void exclude( int i )
  { eslot_.push_back( i ); }

  template< class Coef >
  void fix( Coef *ptr, const Coef *lhs, const Coef *rhs ) const
  {
    for ( typename ESlot::const_iterator it = eslot_.begin(); it != eslot_.end(); ++it )
      *ptr -= fix_prod( lhs[ *it ], rhs[ *it ] );
    MPI_Allreduce( MPI_IN_PLACE, static_cast< void * >( ptr ), 1, mpi_< Coef >().type, MPI_SUM, comm_ );
  }

  bool all_true( const bool flg ) const
  {
    int result = flg ? 0 : 1;

    MPI_Allreduce( MPI_IN_PLACE, &result, 1, mpi_< int >().type, MPI_SUM, comm_ );

    return result == 0;
  }
};

}
#endif

#endif//__ELAI_COHERENCE__

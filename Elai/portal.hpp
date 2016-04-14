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
#ifndef __ELAI_PORTAL__
#define __ELAI_PORTAL__

#include "def.hpp"

#ifdef ELAI_USE_MPI
#include "mpi.h"

namespace elai
{

class portal
{
  MPI_Comm comm_;
  MPI_Status status_;
  int pair_, self_;

public:
  portal( int pair, MPI_Comm comm )
    : comm_( comm ), status_(), pair_( pair ), self_( -1 )
  {
    MPI_Comm_rank( comm_, &self_ );
  }
  portal( const portal& src )
    : comm_( src.comm_ ), status_(), pair_( src.pair_ ), self_( src.self_ )
  {}

  void send( const void *ptr, MPI_Datatype type, int n = 1 ) const
  {
#ifdef ELAI_USE_MPI3
    MPI_Send( ptr, n, type, pair_, self_, comm_ );
#else
    MPI_Send( const_cast< void * >( ptr ), n, type, pair_, self_, comm_ );
#endif
  }

  void recv( void *ptr, MPI_Datatype type, int n = 1 )
  {
    MPI_Recv( ptr, n, type, pair_, pair_, comm_, &status_ );
  }

  template< class Target >
  const Target& operator()( const Target& obj ) const
  {
    obj.marshalize( *this );

    return obj;
  }

  template< class Target >
  void reinforce( const Target& obj ) const { obj.reinforce( *this ); }
};

}
#endif

#endif//__ELAI_PORTAL__

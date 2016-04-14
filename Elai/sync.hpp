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
#ifndef __ELAI_SYNC__
#define __ELAI_SYNC__

#include "def.hpp"

#ifdef ELAI_USE_MPI
#include "mpi.h"

namespace elai
{

class sync
{
  MPI_Comm comm_;
  void *ptr_;
  int cnt_;
  MPI_Datatype type_;

public:
  template< class Target >
  sync( Target& obj, MPI_Comm comm )
    : comm_( comm ), ptr_( NULL ), cnt_( 0 ), type_()
  {
    obj.sync_setup( *this );
  }
  ~sync()
  {
  }

  void setup( void *ptr, int cnt, MPI_Datatype type )
  {
    ptr_ = ptr;
    cnt_ = cnt;
    type_ = type;
  }

  void operator()()
  {
    MPI_Allreduce( MPI_IN_PLACE, ptr_, cnt_, type_, MPI_SUM, comm_ );
  }
};

}
#endif

#endif//__ELAI_SYNC__

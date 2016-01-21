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
#ifndef __ELAI_LU__
#define __ELAI_LU__

#include "def.hpp"
#include "vector.hpp"
#include "matrix.hpp"

#ifdef ELAI_USE_MPI
#include "mpi.h"
#define ELAI_USE_LU             \
public: /* FOR ODD COMPILERS */ \
  using lu< Coef >::A_;         \
  using lu< Coef >::mem_;       \
  using lu< Coef >::comm_;      \
  using lu< Coef >::size_;      \
  using lu< Coef >::rank_ // ; is missed advisedly.
#else
#define ELAI_USE_LU             \
public: /* FOR ODD COMPILERS */ \
  using lu< Coef >::A_;         \
  using lu< Coef >::mem_ // ; is missed advisedly.
#endif

namespace elai
{

template < class Coef >
class lu
{
protected:
  matrix< Coef >& A_;
  size_t mem_;

#ifdef ELAI_USE_MPI
  MPI_Comm comm_;
  int size_, rank_;
#endif

  virtual bool factor_() = 0;
  virtual bool solve_( const vector< Coef >& b, vector< Coef >& x ) = 0;

public:
  lu
    ( matrix< Coef >& A
#ifdef ELAI_USE_MPI
    , MPI_Comm comm
#endif
    )
    : A_( A ), mem_( 0 )
#ifdef ELAI_USE_MPI
    , comm_( comm )
#endif
  {
#ifdef ELAI_USE_MPI
    MPI_Comm_size( comm_, &size_ );
    MPI_Comm_rank( comm_, &rank_ );
#endif
  }
  ~lu()
  {}

  bool factor()
  {
    return factor_();
  }

  bool solve( const vector< Coef >& b, vector< Coef >& x )
  {
    return solve_( b, x );
  }

  size_t mem() const { return mem_; }
};

}

#endif//__ELAI_LU__

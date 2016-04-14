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
#ifndef __ELAI_PRECONDITIONER__
#define __ELAI_PRECONDITIONER__

#include "def.hpp"
#include "expression.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"

namespace elai
{

template< class Range >
class preconditioner
{
protected:
  const matrix< Range >& A_;

  virtual void forward_( vector< Range >& x ) const = 0;
  virtual void backward_( vector< Range >& x ) const = 0;
  virtual void forwardInv_( vector< Range >& x ) const = 0;
  virtual void backwardInv_( vector< Range >& x ) const = 0;

public:
  preconditioner( const matrix< Range >& A )
    : A_( A )
  {}
  virtual ~preconditioner() {}

  void forward( vector< Range >& x ) const     // x -> L^-1 x
  {
    forward_( x );
  }
  void backward( vector< Range >& x ) const    // x -> U^-1 x
  {
    backward_( x );
  }
  void forwardInv( vector< Range >& x ) const  // x -> U x
  {
    forwardInv_( x );
  }
  void backwardInv( vector< Range >& x ) const // x -> L x
  {
    backwardInv_( x );
  }
};

}

#endif//__ELAI_PRECONDITIONER__

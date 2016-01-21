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
#ifndef __ELAI_JACOBI_CONDITIONER__
#define __ELAI_JACOBI_CONDITIONER__

#include "def.hpp"
#include "expression.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"
#include "preconditioner.hpp"

namespace elai
{

// A = L + D + U
// Then M = D

template< class Range >
class jacobi_conditioner : public preconditioner< Range >
{
public:
  jacobi_conditioner( const matrix< Range >& A ) : preconditioner< Range >( A ) {}

protected:
  void forward_( vector< Range >& x ) const {}
  void backward_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    for ( int i = 0; i < x.m(); ++i ) x( i ) /= A( i, i );
  }
  void forwardInv_( vector< Range >& x ) const {}
  void backwardInv_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    for ( int i = 0; i < x.m(); ++i ) x( i ) *= A( i, i );
  }
};

}

#endif//__ELAI_JACOBI_CONDITIONER__

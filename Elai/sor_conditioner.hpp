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
#ifndef __ELAI_SOR_CONDITIONER__
#define __ELAI_SOR_CONDITIONER__

#include "def.hpp"
#include "expression.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"
#include "preconditioner.hpp"

namespace elai
{

// A = L + D + U
// Then M = ( D / w + L ) * w / ( 2 - w ) * D^-1 * ( D / w + U)
//          -------------   -----------------------------------
//        =      L'       *                  U'
// such that 0 < w < 2

template< class Range >
class sor_conditioner : public preconditioner< Range >
{
  int *diag_;
  Range omega_, iomega_;

protected:
  void forward_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = A.ind();
    const int *col = A.col();
    const Range *coef = A.val();

    for ( int i = 0; i < A.m(); ++i )
    {
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( i < j ) break;
        else if ( i == j ) x( i ) /= coef[ k ] * iomega_;
        else x( i ) -= x( j ) * coef[ k ];
      }
    }
  }
  void backward_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = A.ind();
    const int *col = A.col();
    const Range *coef = A.val();
    const Range omega = ( static_cast< Range >( 2. ) - omega_ );

    for ( int i = A.m() - 1; 0 <= i ; --i )
    {
      const Range diag = coef[ diag_[ i ] ];

      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( j <= i ) continue;
        else x( i ) -=  x( j ) * omega_ * coef[ k ] / ( omega * diag );
      }
      x( i ) *= omega;
    }
  }
  void forwardInv_( vector< Range >& x ) const
  {
    const vector< Range > tmp( x );
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = A.ind();
    const int *col = A.col();
    const Range *coef = A.val();

    for ( int i = 0; i < A.m(); ++i )
    {
      x( i ) = static_cast< Range >( 0. );
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( i < j ) break;
        else if ( i == j ) x( i ) += tmp( j ) * coef[ k ] * iomega_;
        else x( i ) += tmp( j ) * coef[ k ];
      }
    }
  }
  void backwardInv_( vector< Range >& x ) const
  {
    const vector< Range > tmp( x );
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = A.ind();
    const int *col = A.col();
    const Range *coef = A.val();
    const Range omega = ( static_cast< Range >( 2. ) - omega_ );

    for ( int i = 0; i < A.m(); ++i )
    {
      const Range diag = coef[ diag_[ i ] ];

      x( i ) = static_cast< Range >( 0. );
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( i < j ) continue;
        else if ( i == j ) x( i ) += tmp( j ) / omega;
        else x( i ) += tmp( j ) * omega_ * coef[ k ] / ( omega * diag );
      }
    }
  }

public:
  sor_conditioner( const matrix< Range >& A, Range omega )
    : preconditioner< Range >( A ), diag_( NULL ), omega_( omega )
  {
    diag_ = new int[ A.m() ];
    for ( int i = 0; i < A.m(); ++i )
      for ( int k = A.ind()[ i ]; k < A.ind()[ i + 1 ]; ++k )
        if ( A.col()[ k ] == i ) { diag_[ i ] = k; break; }
    if ( omega_ < static_cast< Range >( .0 ) ) omega_ = static_cast< Range >( 1.03 );
    else if ( static_cast< Range >( 2. ) < omega_ ) omega_ = static_cast< Range >( 1.03 );
    iomega_ = static_cast< Range >( 1. ) / omega_;
  }
  ~sor_conditioner()
  {
    delete [] diag_;
  }
};

}

#endif//__ELAI_SOR_CONDITIONER__

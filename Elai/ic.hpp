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
#ifndef __ELAI_IC__
#define __ELAI_IC__

#include "def.hpp"
#include "expression.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"
#include "preconditioner.hpp"
#include "fillin.hpp"

namespace elai
{

template< class Range >
class ic : public preconditioner< Range >
{
  fillin< Range > prec_;
  int *diag_;

protected:
  void forward_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = prec_.xadj();
    const int *col = prec_.adjy();
    const Range *coef = prec_.coef();

    for ( int i = 0; i < A.m(); ++i )
    {
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( i <= j ) break; // L( i, i ) = 1.;
        x( i ) -= coef[ k ] * x( j );
      }
    }
  }
  void backward_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = prec_.xadj();
    const int *col = prec_.adjy();
    const Range *coef = prec_.coef();

    for ( int i = A.m() - 1; 0 <= i; --i )
    {
      Range diag = static_cast< Range >( 1 );

      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( j < i ) continue;
        else if ( j == i ) diag = coef[ k ];
        else if ( i < j ) x( i ) -= diag * coef[ k ] * x( j );
      }
      x( i ) /= diag;
    }
  }
  void forwardInv_( vector< Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = prec_.xadj();
    const int *col = prec_.adjy();
    const Range *coef = prec_.coef();
    vector< Range > tmp( x );

    for ( int i = 0; i < A.m(); ++i )
    {
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( i <= j ) break;
        tmp( i ) += coef[ k ] * x( j );
      }
    }
    x = tmp;
  }
  void backwardInv_( vector<Range >& x ) const
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = prec_.xadj();
    const int *col = prec_.adjy();
    const Range *coef = prec_.coef();
    vector< Range > tmp( x );

    for ( int i = 0; i < A.m(); ++i )
    {
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        int j = col[ k ];

        if ( j < i ) continue;
        tmp( i ) += coef[ k ] * x( j );
      }
    }
    x = tmp;
  }

public:
  ic( const matrix< Range >& A, int lv = 0 )
    : preconditioner< Range >( A ), prec_( A ), diag_( NULL )
  {
    prec_( lv );
    diag_ = new int[ A.m() ];
  }
  ~ic()
  {
    // DO NOTHING! OWNERSHIPS ARE OTHERS!! WITHOUT diag_.
    delete [] diag_;
  }

  void factor()
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = prec_.xadj();
    const int *col = prec_.adjy();
    Range *coef = prec_.coef();

    prec_.setup();
    // L-part of coef: i < j, L( i, i ) = 1 is the implicit assumption.
    // L^T-part of coef: j < i
    // D-part of coef: i == j
    for ( int i = 0; i < A.m(); ++i )
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
        if ( col[ k ] == i ) { diag_[ i ] = k; break; }
    for ( int off = ind[ 0 ] + 1; off < ind[ 1 ]; ++off )
      coef[ off ] /= coef[ 0 ];
    for ( int i = 1; i < A.m(); ++i )
    {
      for ( int off = ind[ i ]; off < ind[ i + 1 ]; ++off )
      {
        int j = col[ off ];

        if ( i == j ) continue;
        else if ( i < j )
        {
          coef[ off ] /= coef[ diag_[ i ] ];
          continue;
        }
        for ( int off1 = off + 1, off2 = ind[ j ]
            ; off1 < ind[ i + 1 ] && off2 < ind[ j + 1 ]
            ;)
        {
          int j1 = col[ off1 ]; // j+1 < j1
          int j2 = col[ off2 ];

          if ( j1 < j2 ) { ++off1; continue; }
          else if ( j2 < j1 ) { ++off2; continue; }
          coef[ off1 ] -= coef[ off ] * coef[ off2 ];
          ++off1; ++off2;
        }
        coef[ off ] /= coef[ diag_[ j ] ];
      }
    }
  }
};

}

#endif//__ELAI_IC__

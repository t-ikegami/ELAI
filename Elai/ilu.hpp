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
#ifndef __ELAI_ILU__
#define __ELAI_ILU__

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
class ilu : public preconditioner< Range >
{
  fillin< Range > prec_;
  Range thr_;

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

        if ( fabs( coef[ k ] ) <= thr_ ) continue;
        if ( i <= j ) break;
        x( i ) -= coef[ k ] * x( j );
      }
      // L( i, i ) = 1.
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

        if ( fabs( coef[ k ] ) <= thr_ ) continue;
        if ( j < i ) continue;
        else if ( j == i ) diag = coef[ k ];
        else if ( i < j ) x( i ) -= coef[ k ] * x( j );
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

        if ( fabs( coef[ k ] ) <= thr_ ) continue;
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

        if ( fabs( coef[ k ] ) <= thr_ ) continue;
        if ( j < i ) continue;
        tmp( i ) += coef[ k ] * x( j );
      }
    }
    x = tmp;
  }

public:
  ilu
    ( const matrix< Range >& A
    , int lv = 0
    , Range thr = static_cast< Range >( 0e0 )
    , bool is_srule = false
    )
    : preconditioner< Range >( A ), prec_( A ), thr_( thr )
  {
    prec_( lv, thr_, is_srule );
  }
  ~ilu()
  {
    // DO NOTHING! OWNERSHIPS ARE OTHERS!!
  }

  void factor( Range thr = static_cast< Range >( -1e0 ) )
  {
    const matrix< Range >& A = preconditioner< Range >::A_;
    const int *ind = prec_.xadj();
    const int *col = prec_.adjy();
    Range *coef = prec_.coef();
    int diag[ A.m() ];

    if ( 0 < thr ) thr_ = thr;
    // Copy coefficients from A to prec
    prec_.setup();
    // L-part of coef: i < j, L( i, i ) = 1 is assumed implicitly.
    // U-part of coef: i <=j
    for ( int i = 0; i < A.m(); ++i )
    {
      diag[ i ] = -1;
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
        if ( col[ k ] == i ) { diag[ i ] = k; break; }
    }
    for ( int i = 1; i < A.m(); ++i )
    {
      for ( int off = ind[ i ]; off < ind[ i + 1 ]; ++off )
      {
        int j = col[ off ];

        if ( i <= j || diag[ j ] < 0 ) break;
        coef[ off ] /= coef[ diag[ j ] ];
        //if ( fabs( coef[ off ] ) <= thr ) continue;
        if ( fabs( coef[ off ] ) <= thr ) { coef[ off ] = static_cast< Range >( 0 ); continue; }
        for ( int off1 = off + 1, off2 = ind[ j ]
            ; off1 < ind[ i + 1 ] && off2 < ind[ j + 1 ]
            ;
            )
        {
          int j1 = col[ off1 ]; // j+1 < j1
          int j2 = col[ off2 ];

          if ( j1 < j2 ) { ++off1; continue; }
          else if ( j2 < j1 ) { ++off2; continue; }
          //if ( fabs( coef[ off2 ] ) <= thr ) { ++off1; ++off2; continue; }
          if ( fabs( coef[ off2 ] ) <= thr ) { coef[ off2 ] = static_cast< Range >( 0 ); ++off1; ++off2; continue; }
          coef[ off1 ] -= coef[ off ] * coef[ off2 ];
          ++off1; ++off2;
        }
      }
    }
  }
};

}

#endif//__ELAI_ILU__

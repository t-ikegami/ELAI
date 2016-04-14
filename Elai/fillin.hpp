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
#ifndef __ELAI_FILLIN__
#define __ELAI_FILLIN__

#include <algorithm>
#include <set>
#include <utility>
#include "def.hpp"
#include "matrix.hpp"

namespace elai
{

template< class Coef >
class fillin
{
  typedef matrix< Coef > Matrix;

  class fitem
  {
  public:
    fitem( int c, int l ) : col( c ), lvl( l ) {}
    fitem( const fitem& src ) : col( src.col ), lvl( src.lvl ) {}
    bool operator<( const fitem& lhs ) const { return col < lhs.col; }
    const int col;
    const int lvl;
  };

  const Matrix& A_;
  int nnz_;
  int *xadj_;
  int *adjy_;
  Coef *coef_;

  void sysfactor( int lv, Coef thres )
  {
    const int *ind = A_.ind();
    const int *col = A_.col();
    const Coef *val = A_.val();
    int m = A_.m(), n = A_.n();
    std::set < fitem > *adjs;

    adjs = new std::set< fitem >[ m ];

    for ( int i = 0; i < m; ++i )
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
        if ( ( col[ k ] == i ) || ( thres <= fabs( val[ k ] ) ) )
          adjs[ i ].insert( fitem( col[ k ], 0 ) );

    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];

      for ( typename std::set< fitem >::const_iterator it = adj.begin()
          ; it != adj.end(); ++it
          )
      {
        const fitem item( *it );

        if ( i <= item.col ) break;
        if ( lv <= item.lvl ) continue;
        for ( typename std::set< fitem >::const_iterator jt = adjs[ item.col ].begin()
            ; jt != adjs[ item.col ].end(); ++jt
            )
        {
          const fitem jtem( *jt );

          if ( i < jtem.col ) break;
          if ( ( jtem.col <= item.col ) || ( lv <= jtem.lvl ) ) continue;
          adj.insert( fitem( jtem.col, item.lvl + jtem.lvl + 1 ) );
          adjs[ jtem.col ].insert( fitem( i, item.lvl + jtem.lvl + 1 ) );
        }
      }
    }
#ifdef ELAI_DEBUG
    int acc = 0;
    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];
      acc += adj.size();
    }
    std::cerr << "N=" << A_.m() << ", NZ=" << A_.nnz() << ", NF=" << acc << std::endl;
#endif

    build( adjs );
    delete [] adjs;
  }

  void symfactor( int lv, Coef thres )
  {
    const int *ind = A_.ind();
    const int *col = A_.col();
    const Coef *val = A_.val();
    int m = A_.m(), n = A_.n();
    std::set < fitem > *adjs;

    adjs = new std::set< fitem >[ m ];

    for ( int i = 0; i < m; ++i )
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
        if ( ( col[ k ] == i ) || ( thres <= fabs( val[ k ] ) ) )
          adjs[ i ].insert( fitem( col[ k ], 0 ) );

    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];

      for ( typename std::set< fitem >::const_iterator it = adj.begin()
          ; it != adj.end(); ++it
          )
      {
        const fitem item( *it );

        if ( i <= item.col ) break;
        if ( lv <= item.lvl ) continue;
        for ( typename std::set< fitem >::const_iterator jt = adjs[ item.col ].begin()
            ; jt != adjs[ item.col ].end(); ++jt
            )
        {
          const fitem jtem( *jt );

          if ( i < jtem.col ) break;
          if ( ( jtem.col <= item.col ) || ( lv <= jtem.lvl ) ) continue;
          adj.insert( fitem( jtem.col, std::max( item.lvl, jtem.lvl ) + 1 ) );
          adjs[ jtem.col ].insert( fitem( i, std::max( item.lvl, jtem.lvl ) + 1 ) );
        }
      }
    }
#ifdef ELAI_DEBUG
    int acc = 0;
    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];
      acc += adj.size();
    }
    std::cerr << "N=" << A_.m() << ", NZ=" << A_.nnz() << ", NF=" << acc << std::endl;
#endif

    build( adjs );
    delete [] adjs;
  }

  void sfactor( int lv, Coef thres )
  {
    const int *ind = A_.ind();
    const int *col = A_.col();
    const Coef *val = A_.val();
    int m = A_.m(), n = A_.n();
    std::set < fitem > *adjs;

    adjs = new std::set< fitem >[ m ];

    for ( int i = 0; i < m; ++i )
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
        if ( ( col[ k ] == i ) || ( thres <= fabs( val[ k ] ) ) )
          adjs[ i ].insert( fitem( col[ k ], 0 ) );

    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];

      for ( typename std::set< fitem >::const_iterator it = adj.begin()
          ; it != adj.end(); ++it
          )
      {
        const fitem item( *it );

        if ( i <= item.col ) break;
        if ( lv <= item.lvl ) continue;
        for ( typename std::set< fitem >::const_iterator jt = adjs[ item.col ].begin()
            ; jt != adjs[ item.col ].end(); ++jt
            )
        {
          const fitem jtem( *jt );

          if ( ( jtem.col <= item.col ) || ( lv <= jtem.lvl ) ) continue;
          adj.insert( fitem( jtem.col, item.lvl + jtem.lvl + 1 ) );
        }
      }
    }
#ifdef ELAI_DEBUG
    int acc = 0;
    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];
      acc += adj.size();
    }
    std::cerr << "N=" << A_.m() << ", NZ=" << A_.nnz() << ", NF=" << acc << std::endl;
#endif

    build( adjs );
    delete [] adjs;
  }

  void mfactor( int lv, Coef thres )
  {
    const int *ind = A_.ind();
    const int *col = A_.col();
    const Coef *val = A_.val();
    int m = A_.m(), n = A_.n();
    std::set < fitem > *adjs;

    adjs = new std::set< fitem >[ m ];

    for ( int i = 0; i < m; ++i )
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
        if ( ( col[ k ] == i ) || ( thres <= fabs( val[ k ] ) ) )
          adjs[ i ].insert( fitem( col[ k ], 0 ) );

    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];

      for ( typename std::set< fitem >::const_iterator it = adj.begin()
          ; it != adj.end(); ++it
          )
      {
        const fitem item( *it );

        if ( i <= item.col ) break;
        if ( lv <= item.lvl ) continue;
        for ( typename std::set< fitem >::const_iterator jt = adjs[ item.col ].begin()
            ; jt != adjs[ item.col ].end(); ++jt
            )
        {
          const fitem jtem( *jt );

          if ( ( jtem.col <= item.col ) || ( lv <= jtem.lvl ) ) continue;
          adj.insert( fitem( jtem.col, std::max( item.lvl, jtem.lvl ) + 1 ) );
        }
      }
    }
#ifdef ELAI_DEBUG
    int acc = 0;
    for ( int i = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];
      acc += adj.size();
    }
    std::cerr << "N=" << A_.m() << ", NZ=" << A_.nnz() << ", NF=" << acc << std::endl;
#endif

    build( adjs );
    delete [] adjs;
  }

  void build( std::set< fitem > *adjs )
  {
    int m = A_.m();

    nnz_ = 0;
    xadj_ = new int[ m + 1 ];
    xadj_[ 0 ] = nnz_;
    for ( int i = 0; i < m; ++i )
    {
      nnz_ += adjs[ i ].size();
      xadj_[ i + 1 ] = nnz_;
    }
    adjy_ = new int[ nnz_ ];
    coef_ = new Coef[ nnz_ ];
    for ( int k = 0; k < nnz_; ++k ) coef_[ k ] = static_cast< Coef >( 0 );
    for ( int i = 0, k = 0; i < m; ++i )
    {
      std::set< fitem >& adj = adjs[ i ];

      for ( typename std::set< fitem >::const_iterator it = adj.begin()
          ; it != adj.end(); ++it
          ) adjy_[ k++ ] = fitem( *it ).col;
    }
    // Setup of coefficients moved to an independent method.
    // Because of scaling, coefficients are to be modified.
  }

  void destruct()
  {
    nnz_ = 0;
    if ( coef_ != NULL ) { delete [] coef_; coef_ = NULL; }
    if ( adjy_ != NULL ) { delete [] adjy_; adjy_ = NULL; }
    if ( xadj_ != NULL ) { delete [] xadj_; xadj_ = NULL; }
  }

public:
  fillin( const Matrix& A )
    : A_( A ), nnz_( 0 ), xadj_( NULL ), adjy_( NULL ), coef_( NULL )
  {}
  ~fillin() { destruct(); }

  void operator()
    ( int k = 0                 // Fill-in Level
    , Coef thres = static_cast< Coef >( 0e0 ) // Fill-in Threshold
    , bool is_srule = false     // Fill-in Rule true: Sum false: Max
    )
  {
    bool is_symm = A_.is_symmetric();

    destruct();
    if ( is_symm && is_srule ) sysfactor( k, thres );
    else if ( is_symm ) symfactor( k, thres );
    else if ( is_srule ) sfactor( k, thres );
    else mfactor( k, thres );
  }

  void setup()
  {
    if ( coef_ != NULL )
    {
      for ( int k = 0; k < nnz_; ++k ) coef_[ k ] = static_cast< Coef >( 0. );
      for ( int i = 0; i < A_.m(); ++i )
      {
        for ( int k0 = A_.ind( i ), k1 = xadj_[ i ]
            ; k0 < A_.ind( i + 1 ) && k1 < xadj_[ i + 1 ]
            ; )
        {
          int j0 = A_.col( k0 ), j1 = adjy_[ k1 ];

          if ( j0 < j1 ) { ++k0; continue; }
          else if ( j1 < j0 ) { ++k1; continue; }
          coef_[ k1 ] = A_.val( k0 );
          ++k0; ++k1;
        }
      }
    }
  }

  int m() const { return A_.m(); }
  int n() const { return A_.n(); }
  int nnz() const { return nnz_; }
  int *xadj() const { return xadj_; }
  int *adjy() const { return adjy_; }
  Coef *coef() { return coef_; }
  Coef *coef() const { return coef_; }
};

}

#endif//__ELAI_FILLIN__

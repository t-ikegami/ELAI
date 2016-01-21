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
#ifndef __ELAI_KSP__
#define __ELAI_KSP__

#include <iostream>
#include "def.hpp"
#include "coherence.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"
#include "preconditioner.hpp"
#include "util.hpp"

#ifdef ELAI_USE_MPI
#define ELAI_USE_KSP                \
public: /* FOR ODD COMPILERS */     \
  using ksp< Coef >::A_;            \
  using ksp< Coef >::b_;            \
  using ksp< Coef >::P_;            \
  using ksp< Coef >::res_;          \
  using ksp< Coef >::athres_;       \
  using ksp< Coef >::rthres_;       \
  using ksp< Coef >::elapsed_;      \
  using ksp< Coef >::prec_elapsed_; \
  using ksp< Coef >::sync;          \
  using ksp< Coef >::isOK;          \
  using ksp< Coef >::fix;           \
  using ksp< Coef >::fix_norm;      \
  using ksp< Coef >::sync_norm;     \
  using ksp< Coef >::is_trivial;    \
  using ksp< Coef >::abs_converged; \
  using ksp< Coef >::rel_converged; \
  using ksp< Coef >::iter_max;      \
  using ksp< Coef >::mem // ; is missed advisedly.
#else
#define ELAI_USE_KSP                \
public: /* FOR ODD COMPILERS */     \
  using ksp< Coef >::A_;            \
  using ksp< Coef >::b_;            \
  using ksp< Coef >::P_;            \
  using ksp< Coef >::res_;          \
  using ksp< Coef >::athres_;       \
  using ksp< Coef >::rthres_;       \
  using ksp< Coef >::elapsed_;      \
  using ksp< Coef >::prec_elapsed_; \
  using ksp< Coef >::isOK;          \
  using ksp< Coef >::fix_norm;      \
  using ksp< Coef >::sync_norm;     \
  using ksp< Coef >::is_trivial;    \
  using ksp< Coef >::abs_converged; \
  using ksp< Coef >::rel_converged; \
  using ksp< Coef >::iter_max;      \
  using ksp< Coef >::mem // ; is missed advisedly.
#endif

#ifdef ELAI_USE_MPI
#define ELAI_SYNC( v ) sync( ( v ) )
#define ELAI_PROD( acc, x, y )    \
  do {                            \
    ( acc ) = ( x ) * ( y );      \
    fix( ( acc ), ( x ), ( y ) ); \
  } while ( 0 )
#else
#define ELAI_SYNC( v )
#define ELAI_PROD( acc, x, y )    \
  do {                            \
    ( acc ) = ( x ) * ( y );      \
  } while ( 0 )
#endif

namespace elai
{

template< class Coef >
class ksp
{
protected:
  const matrix< Coef >& A_;
  const vector< Coef >& b_;
  const preconditioner< Coef > *P_;

  vector< Coef > res_;
  int iter_max_;
  Coef athres_, rthres_;
  double elapsed_, prec_elapsed_;

#ifdef ELAI_USE_MPI
  coherence *coherent_;

  inline void sync( vector< Coef >& u ) const
  { if ( coherent_ != NULL ) ( *coherent_ )( u.val() ); }

  inline void fix( Coef& acc, const vector< Coef >& u, const vector< Coef >& v ) const
  { if ( coherent_ != NULL ) coherent_->fix( &acc, u.val(), v.val() ); }

  inline bool isOK( const bool flg ) const
  {
    if ( coherent_ != NULL ) return coherent_->all_true( flg );
    else return flg;
  }
#else
  inline bool isOK( const bool flg ) const
  { return flg; }
#endif

  Coef fix_norm( const vector< Coef >& x ) const
  {
    Coef v;

    v = x * x;
#ifdef ELAI_USE_MPI
    fix( v, x, x );
#endif

    v = sqrt( v );

    return v;
  }

  Coef sync_norm( vector< Coef >& x ) const
  {
#ifdef ELAI_USE_MPI
    sync( x );
#endif

    return fix_norm( x );
  }

  bool is_trivial( Coef v ) const
  { return fabs( v ) <= static_cast< Coef >( 1e-50 ); }
  bool is_trivial( const vector< Coef >& x ) const
  { return fix_norm( x ) <= static_cast< Coef >( 1e-50 ); }

  bool abs_converged( const Coef v )
  { return fabs( v ) <= athres_; }
  bool abs_converged( const vector< Coef >& x )
  {
    res_ = A_ * x - b_;

    return sync_norm( res_ ) <= athres_;
  }

  bool rel_converged( const Coef r, const Coef r0 )
  { return fabs( r / r0 ) <= rthres_; }
  bool rel_converged( const vector< Coef >& x, const Coef r0 )
  {
    res_ = A_ * x - b_;

    return ( sync_norm( res_ ) / r0 ) <= rthres_;
  }

  virtual bool solve_( vector< Coef >& x ) = 0;
  virtual bool solveP_( vector< Coef >& x ) = 0;

public:
  ksp
    ( const matrix< Coef >& A
    , const vector< Coef >& b
    , const preconditioner< Coef > *P = NULL
#ifdef ELAI_USE_MPI
    , coherence *coherent = NULL
#endif
    )
    : A_( A ), b_( b ), P_( P ), res_( A_.m() )
    , iter_max_( A_.m() / 2 )
    , athres_( static_cast< Coef >( 1e-30 ) )
    , rthres_( static_cast< Coef >( 1e-12 ) )
    , elapsed_( 0. ), prec_elapsed_( 0. )
#ifdef ELAI_USE_MPI
    , coherent_( coherent )
#endif
  {}
  ~ksp() {}

  int iter_max() const { return iter_max_; }
  int iter_max( int max )
  {
    int old = iter_max_;

    iter_max_ = max;

    return old;
  }

  Coef abs_thres() const { return athres_; }
  Coef abs_thres( Coef thres )
  {
    Coef old = athres_;

    athres_ = thres;

    return old;
  }

  Coef rel_thres() const { return rthres_; }
  Coef rel_thres( Coef thres )
  {
    Coef old = rthres_;

    rthres_ = thres;

    return old;
  }

  double elapsed() const { return elapsed_; }
  double prec_elapsed() const { return prec_elapsed_; }

  bool solve( vector< Coef >& x )
  {
    bool flg;

    elapsed_ = 0;
    prec_elapsed_ = 0;
    ELAI_PROF_BEG( elapsed_ );
    if ( P_ == NULL ) flg = solve_( x );
    else flg = solveP_( x );
    ELAI_PROF_END( elapsed_ );

    return flg;
  }

  size_t mem() const
  {
    size_t sum = 0;

    sum += sizeof( P_ );
    sum += res_.mem();
    sum += sizeof( iter_max_ );
    sum += sizeof( elapsed_ );
    sum += sizeof( prec_elapsed_ );
#ifdef ELAI_USE_MPI
    sum += sizeof( coherent_ );
#endif

    return sum;
  }

};

}

#endif//__ELAI_KSP__

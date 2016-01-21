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
#ifndef __ELAI_GMRES__
#define __ELAI_GMRES__

#include "def.hpp"
#include "coherence.hpp"
#include "expression.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "blas.hpp"
#include "preconditioner.hpp"
#include "ksp.hpp"

namespace elai
{

template< class Coef >
class gmres : public ksp< Coef >
{
  ELAI_USE_KSP;

private:
  // WORKSPACE
  vector< Coef > *v_, y_, r_, c_, s_, e_, h_;
  int restart_;

  void setup_()
  {
    if ( v_ != NULL ) { delete [] v_; v_ = NULL; }
    v_ = new vector< Coef >[ restart_ + 1 ];
    for ( int i = 0; i <= restart_; ++i ) v_[ i ].setup( A_.m() );

    y_.setup( restart_ );
    c_.setup( restart_ );
    s_.setup( restart_ );
    e_.setup( restart_ + 1 );

    // Hessemberg matrix in vec.
    {
      int m = restart_ + 1;
      h_.setup( m * m - m - ( m - 2 ) * ( m - 1 ) / 2 );
    }
  }

  Coef& h( int i, int j )
  {
    if ( 1 < i )
    {
      int suffix = i - 1;
      int k = suffix * ( suffix - 1 ) / 2;

      return h_( i * restart_ - k + ( j - suffix ) );
    }
    else return h_( i * restart_ + j );
  }
  const Coef& h( int i, int j ) const { return h( i, j ); }

  void back_subst( vector< Coef >& y, int k, const vector< Coef >& x )
  {
    for ( int i = k; 0 <= i; --i )
    {
      y( i ) = x( i );
      for ( int j = k; i < j; --j ) y( i ) -= h( i, j ) * y( j );
      y( i ) /= h( i, i );
    }
  }

  bool solve_( vector< Coef >& x )
  {
    const int n = A_.m();
    int itr = 0;
    bool converged = false;
    Coef res, res0;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    if ( is_trivial( res ) ) return true;

    v_[ 0 ] = ( static_cast< Coef >( 1e0 ) / res ) * r_;
    e_.clear( static_cast< Coef >( 0e0 ) ); e_( 0 ) = res;
    while ( !converged )
    {
      int m;

      for ( m = 0; m < restart_; ++m )
      {
        v_[ m + 1 ] = A_ * v_[ m ];

        // Orthogonalization
        for ( int i = 0; i <= m; ++i ) ELAI_PROD( h( i, m ), v_[ m + 1 ], v_[ i ] );
        for ( int i = 0; i <= m; ++i ) v_[ m + 1 ] = v_[ m + 1 ] - h( i, m ) * v_[ i ];

        // Normalization
        h( m + 1, m ) = sync_norm( v_[ m + 1 ] );
        v_[ m + 1 ] = ( static_cast< Coef >( 1e0 ) / h( m + 1, m ) ) * v_[ m + 1 ];

        // Givens Transformation
        //for ( int i = 0; i < m - 1; ++i )
        for ( int i = 0; i < m; ++i )
        {
          Coef gamma;

          gamma = c_( i ) * h( i, m ) - conj_( s_( i ) ) * h( i + 1, m );
          h( i + 1, m ) = s_( i ) * h( i, m ) + c_( i ) * h( i + 1, m );
          h( i, m ) = gamma;
        }

        // Update Givens
        {
          Coef d;
          
          d = sqrt( conj_( h( m, m ) ) * h( m, m ) / ( conj_( h( m, m ) ) * h( m, m ) + conj_( h( m + 1, m ) ) * h( m + 1, m ) ) );
          c_( m ) = d;
          s_( m ) = - h( m + 1, m ) / h( m, m ) * d;
        }

        // Update error
        {
          Coef gamma;

          gamma = c_( m ) * e_( m ) - conj_( s_( m ) ) * e_( m + 1 );
          e_( m + 1 ) = s_( m ) * e_( m ) + c_( m ) * e_( m + 1 );
          e_( m ) = gamma;
        }

        // Update Hessemberg
        h( m, m ) = sqrt( conj_( h( m, m ) ) * h( m, m ) + conj_( h( m + 1, m ) ) * h( m + 1, m ) );
        h( m + 1, m ) = static_cast< Coef >( 0e0 );

#ifdef ELAI_DEBUG
        Coef norm = x * x;

        std::cerr << " ||x||=" << norm
                  << " ||Ax-b||=" << fabs( e_( m ) )
                  << "  " << fabs( e_( m ) ) / res0 << std::endl;
#endif
        // Convergence Check
        if ( rel_converged( e_( m ), res0 ) ) converged = true;
        if ( isOK( converged ) )
        {
          converged = true;
          break;
        }
        else converged = false;
      }

      // Solve via backward-substitute
      {
        int k = m < restart_ ? m : restart_ - 1;

        back_subst( y_, k, e_ );
        for ( int i = 0; i <= k; ++i ) x = x + y_( i ) * v_[ i ];
      }

      // DIVERGED
      if ( iter_max() <= ++itr ) break;

      r_ = b_ - A_ * x;
      res = sync_norm( r_ );
      v_[ 0 ] = ( static_cast< Coef >( 1e0 ) / res ) * r_;
      e_.clear( static_cast< Coef >( 0e0 ) ); e_( 0 ) = res;
    }

    return converged;
  }

  bool solveP_( vector< Coef >& x )
  {
    const int n = A_.m();
    int itr = 0;
    bool converged = false;
    Coef res, res0;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    if ( is_trivial( res ) ) return true;

    v_[ 0 ] = ( static_cast< Coef >( 1e0 ) / res ) * r_;
    e_.clear( static_cast< Coef >( 0e0 ) ); e_( 0 ) = res;
    while ( !converged )
    {
      int m;

      for ( m = 0; m < restart_; ++m )
      {
        // Preconditioning:
        ELAI_PROF_BEG( prec_elapsed_ );
        P_->forward( v_[ m ] );
        P_->backward( v_[ m ] );
        ELAI_SYNC( v_[ m ] );
        ELAI_PROF_END( prec_elapsed_ );

        v_[ m + 1 ] = A_ * v_[ m ];

        // Orthogonalization
        for ( int i = 0; i <= m; ++i ) ELAI_PROD( h( i, m ), v_[ m + 1 ], v_[ i ] );
        for ( int i = 0; i <= m; ++i ) v_[ m + 1 ] = v_[ m + 1 ] - h( i, m ) * v_[ i ];

        // Normalization
        h( m + 1, m ) = sync_norm( v_[ m + 1 ] );
        v_[ m + 1 ] = ( static_cast< Coef >( 1e0 ) / h( m + 1, m ) ) * v_[ m + 1 ];

        // Givens Transformation
        //for ( int i = 0; i < m - 1; ++i )
        for ( int i = 0; i < m; ++i )
        {
          Coef gamma;

          gamma = c_( i ) * h( i, m ) - conj_( s_( i ) ) * h( i + 1, m );
          h( i + 1, m ) = s_( i ) * h( i, m ) + c_( i ) * h( i + 1, m );
          h( i, m ) = gamma;
        }

        // Update Givens
        {
          Coef d;
          
          d = sqrt( conj_( h( m, m ) ) * h( m, m ) / ( conj_( h( m, m ) ) * h( m, m ) + conj_( h( m + 1, m ) ) * h( m + 1, m ) ) );
          c_( m ) = d;
          s_( m ) = - h( m + 1, m ) / h( m, m ) * d;
        }

        // Update error
        {
          Coef gamma;

          gamma = c_( m ) * e_( m ) - conj_( s_( m ) ) * e_( m + 1 );
          e_( m + 1 ) = s_( m ) * e_( m ) + c_( m ) * e_( m + 1 );
          e_( m ) = gamma;
        }

        // Update Hessemberg
        h( m, m ) = sqrt( conj_( h( m, m ) ) * h( m, m ) + conj_( h( m + 1, m ) ) * h( m + 1, m ) );
        h( m + 1, m ) = static_cast< Coef >( 0e0 );

#ifdef ELAI_DEBUG
        Coef norm = x * x;

        std::cerr << " ||x||=" << norm
                  << " ||Ax-b||=" << fabs( e_( m ) )
                  << "  " << fabs( e_( m ) ) / res0 << std::endl;
#endif
        // Convergence Check
        if ( rel_converged( e_( m ), res0 ) ) converged = true;
        if ( isOK( converged ) )
        {
          converged = true;
          break;
        }
        else converged = false;
      }

      // Solve via backward-substitute with preconditioner
      {
        int k = m < restart_ ? m : restart_ - 1;

        back_subst( y_, k, e_ );

        for ( int i = 0; i <= k; ++i )
        {
          v_[ i ] = y_( i ) * v_[ i ];

          // Preconditioning:
          ELAI_PROF_BEG( prec_elapsed_ );
          P_->forward( v_[ i ] );
          P_->backward( v_[ i ] );
          ELAI_SYNC( v_[ i ] );
          ELAI_PROF_END( prec_elapsed_ );

          x = x + v_[ i ];
        }

      }

      // DIVERGED
      if ( iter_max() <= ++itr ) break;

      r_ = b_ - A_ * x;
      res = sync_norm( r_ );
      v_[ 0 ] = ( static_cast< Coef >( 1e0 ) / res ) * r_;
      e_.clear( static_cast< Coef >( 0e0 ) ); e_( 0 ) = res;
    }

    return converged;
  }

public:
  gmres
    ( const matrix< Coef >& A
    , const vector< Coef >& b
    , preconditioner< Coef > *P = NULL
#ifdef ELAI_USE_MPI
    , coherence *coherent = NULL
#endif
    )
    : ksp< Coef >
      ( A, b, P
#ifdef ELAI_USE_MPI
      , coherent
#endif
      )
    , v_( NULL )
    , y_(), r_( A_.m() ), c_(), s_(), e_(), h_()
    , restart_( 50 )
  {
    iter_max( A_.m() / 2 );
    setup_();
  }
  ~gmres()
  {
    if ( v_ != NULL ) { delete [] v_; v_ = NULL; }
  }

  int restart() const { return restart_; }
  int restart( int restart )
  {
    int old = restart_;

    restart_ = restart;
    if ( old != restart_ ) setup_();

    return old;
  }

  size_t mem() const
  {
    size_t sum = ksp< Coef >::mem();

    for ( int i = 0; i < restart_; ++i ) sum += v_[ i ].mem();
    sum += r_.mem();
    sum += c_.mem();
    sum += s_.mem();
    sum += e_.mem();
    sum += h_.mem();
    sum += sizeof( restart_ );

    return sum;
  }
};

}

#endif//__ELAI_GMRES__

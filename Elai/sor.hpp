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
#ifndef __ELAI_SOR__
#define __ELAI_SOR__

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
class sor : public ksp< Coef >
{
  ELAI_USE_KSP;

private:
  // WORKSPACE
  vector< Coef > r_;
  Coef acc_;

  bool solve_( vector< Coef >& x )
  {
    Coef res, res0;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    if ( is_trivial( res0 ) ) return true;

    for ( int i = 0; i < iter_max(); ++i )
    {
      bool converged = false;
      Coef norm;

      norm = x * x;
#ifdef ELAI_DEBUG
      std::cerr << " ||x||=" << norm
                << " ||Ax-b||=" << res
                << "  " << res / res0 << std::endl;
#endif
      if ( std::isnan( res ) ) return false;
      else if ( rel_converged( res, res0 ) || abs_converged( res ) ) converged = true;
      else if ( !is_trivial( norm ) && rel_converged( res, norm ) ) converged = true;

      if ( isOK( converged ) ) return true;

      for ( int i = 0; i < A_.m(); ++i )
      {
        Coef acc = b_( i ), diag = static_cast< Coef >( 1. );

        for ( int k = A_.ind( i ); k < A_.ind( i + 1 ); ++k )
        {
          int j = A_.col( k );

          if ( i == j ) diag = A_.val( k );
          else acc -= A_.val( k ) * x( j );
        }
        x( i ) += acc_ * ( acc / diag - x( i ) );
      }
      ELAI_SYNC( x );

      r_ = b_ - A_ * x;
      res = sync_norm( r_ );
    }

    return false;
  }

  bool solveP_( vector< Coef >& x )
  {
    Coef res, res0;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    if ( is_trivial( res0 ) ) return true;

    for ( int i = 0; i < iter_max(); ++i )
    {
      bool converged = false;
      Coef norm;

      norm = x * x;
#ifdef ELAI_DEBUG
      std::cerr << " ||x||=" << norm
                << " ||Ax-b||=" << res
                << "  " << res / res0 << std::endl;
#endif
      if ( std::isnan( res ) ) return false;
      else if ( rel_converged( res, res0 ) || abs_converged( res ) ) converged = true;
      else if ( !is_trivial( norm ) && rel_converged( res, norm ) ) converged = true;

      if ( isOK( converged ) ) return true;

      // Preconditioning:
      ELAI_PROF_BEG( prec_elapsed_ );
      P_->forward( x );
      P_->backward( x );
      ELAI_SYNC( x );
      ELAI_PROF_END( prec_elapsed_ );

      for ( int i = 0; i < A_.m(); ++i )
      {
        Coef acc = b_( i ), diag = static_cast< Coef >( 1. );

        for ( int k = A_.ind( i ); k < A_.ind( i + 1 ); ++k )
        {
          int j = A_.col( k );

          if ( i == j ) diag = A_.val( k );
          else acc -= A_.val( k ) * x( j );
        }
        x( i ) += acc_ * ( acc / diag - x( i ) );
      }
      ELAI_SYNC( x );

      r_ = b_ - A_ * x;
      res = sync_norm( r_ );
    }

    return false;
  }

public:
  sor
    ( const matrix< Coef >& A
    , const vector< Coef >& b
    , const preconditioner< Coef > *P = NULL
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
    , r_( A_.m() )
    , acc_( static_cast< Coef >( 1. ) )
  {}
  ~sor() {}

  // acc = 1. => Gauss--Seidel
  Coef accel() const { return acc_; }
  Coef accel( Coef acc )
  {
    if ( 0 < acc && acc < 2 )
    {
      Coef old = acc_;
      acc_ = acc;
      return old;
    }
    return acc;
  }

  size_t mem() const
  {
    size_t sum = ksp< Coef >::mem();

    sum += r_.mem();

    return sum;
  }
};

}

#endif//__ELAI_SOR__

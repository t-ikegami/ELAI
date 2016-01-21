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
#ifndef __ELAI_CG__
#define __ELAI_CG__

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
class cg : public ksp< Coef >
{
  ELAI_USE_KSP;

  // WORKSPACE
  vector< Coef > p_, q_, r_, y_, z_;

protected:
  bool solve_( vector< Coef >& x )
  {
    Coef res, res0, rho0;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    rho0 = static_cast< Coef >( 1. );
    p_ = static_cast< Coef >( 0. );
    if ( is_trivial( res0 ) ) return true;

    for ( int i = 0; i < iter_max(); ++i )
    {
      bool converged = false;
      Coef alpha, beta, rho, pq;

      Coef norm = x * x;
#ifdef ELAI_DEBUG
      std::cerr << " ||x||=" << norm
                << " ||Ax-b||=" << res
                << "  " << res / res0 << std::endl;
#endif
      if ( std::isnan( res ) ) return false;
      else if ( rel_converged( res, res0 ) || abs_converged( res ) ) converged = true;
      else if ( !is_trivial( norm ) && rel_converged( res, norm ) ) converged = true;

      if ( isOK( converged ) ) return true;

      ELAI_PROD( rho, r_, r_ );

      beta = rho / rho0;
      rho0 = rho;
      q_ = r_ + beta * p_;
      p_ = q_;

      q_ = A_ * p_;
      ELAI_SYNC( q_ );

      ELAI_PROD( pq, p_, q_ );

      alpha = rho / pq;
      y_ = x + alpha * p_;
      x = y_;
      y_ = r_ - alpha * q_;
      r_ = y_;

      res = fix_norm( r_ );
    }

    return false;
  }

  bool solveP_( vector< Coef >& x )
  {
    Coef res, res0, rho0;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    rho0 = static_cast< Coef >( 1. );
    p_ = static_cast< Coef >( 0. );
    if ( is_trivial( res0 ) ) return true;

    for ( int i = 0; i < iter_max(); ++i )
    {
      bool converged = false;
      Coef alpha, beta, rho, pq;

      Coef norm = x * x;
#ifdef ELAI_DEBUG
      std::cerr << " ||x||=" << norm
                << " ||Ax-b||=" << res
                << "  " << res / res0 << std::endl;
#endif
      if ( std::isnan( res ) ) return false;
      else if ( rel_converged( res, res0 ) || abs_converged( res ) ) converged = true;
      else if ( !is_trivial( norm ) && rel_converged( res, norm ) ) converged = true;

      if ( isOK( converged ) ) return true;

      z_ = r_;

      // Preconditioning:
      ELAI_PROF_BEG( prec_elapsed_ );
      P_->forward( z_ );
      P_->backward( z_ );
      ELAI_SYNC( z_ );
      ELAI_PROF_END( prec_elapsed_ );

      ELAI_PROD( rho, r_, z_ );

      beta = rho / rho0;
      rho0 = rho;
      q_ = z_ + beta * p_;
      p_ = q_;

      q_ = A_ * p_;
      ELAI_SYNC( q_ );

      ELAI_PROD( pq, p_, q_ );

      alpha = rho / pq;
      y_ = x + alpha * p_;
      x = y_;
      y_ = r_ - alpha * q_;
      r_ = y_;

      res = fix_norm( r_ );
    }

    return false;
  }

public:
  cg
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
    , p_( A_.m() ), q_( A_.m() ), r_( A_.m() ), y_( A_.m() ), z_( A_.m() )
  {}
  ~cg() {}

  size_t mem() const
  {
    size_t sum = ksp< Coef >::mem();

    sum += p_.mem();
    sum += q_.mem();
    sum += r_.mem();
    sum += y_.mem();
    sum += z_.mem();

    return sum;
  }
};

}

#endif//__ELAI_CG__

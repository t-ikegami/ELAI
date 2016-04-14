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
#ifndef __ELAI_BICGSTAB__
#define __ELAI_BICGSTAB__

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
class bicgstab : public ksp< Coef >
{
  ELAI_USE_KSP;

private:
  // WORKSPACE
  vector< Coef > r_, r1_, r2_, rs0_, p_, p1_, Ap_, s_, s1_, s2_;
  Coef bthres_;

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

    rs0_ = r_;
    p_ = r_;

    for ( int i = 0; i < iter_max(); ++i )
    {
      bool converged = false;
      Coef alpha, beta, omega, tmp1, tmp2;

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

      Ap_ = A_ * p_;
      ELAI_SYNC( Ap_ );

      ELAI_PROD( tmp1, rs0_, r_ );
      ELAI_PROD( tmp2, rs0_, Ap_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      alpha = tmp1 / tmp2;
      s_ = r_ - alpha * Ap_;

      s1_ = A_ * s_;
      ELAI_SYNC( s1_ );

      ELAI_PROD( tmp1, s1_, s_ );
      ELAI_PROD( tmp2, s1_, s1_ );

      omega = tmp1 / tmp2;
      x = x + alpha * p_ + omega * s_;

      r1_ = s_ - omega * s1_;

      ELAI_PROD( tmp1, rs0_, r1_ );
      ELAI_PROD( tmp2, rs0_, r_ );

      beta = alpha / omega * tmp1 / tmp2;
      p_ = r1_ + beta * ( p_ - omega * Ap_ );

      r_ = r1_;
      res = fix_norm( r_ );
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

    rs0_ = r_;

    // Preconditioning:
    ELAI_PROF_BEG( prec_elapsed_ );
    P_->forward( rs0_ );
    P_->backward( rs0_ );
    ELAI_SYNC( rs0_ );
    ELAI_PROF_END( prec_elapsed_ );

    r1_ = rs0_;
    p_ = rs0_;

    for ( int i = 0; i < iter_max(); ++i )
    {
      bool converged = false;
      Coef alpha, beta, omega, tmp1, tmp2;

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

      Ap_ = A_ * p_;
      ELAI_SYNC( Ap_ );

      p1_ = Ap_;

      // Preconditioning:
      ELAI_PROF_BEG( prec_elapsed_ );
      P_->forward( p1_ );
      P_->backward( p1_ );
      ELAI_SYNC( p1_ );
      ELAI_PROF_END( prec_elapsed_ );

      ELAI_PROD( tmp1, rs0_, r1_ );
      ELAI_PROD( tmp2, rs0_, p1_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      alpha = tmp1 / tmp2;
      s_ = r_ - alpha * Ap_;
      s1_ = r1_ - alpha * p1_;

      s2_ = A_ * s1_;
      ELAI_SYNC( s2_ );

      ELAI_PROD( tmp1, s2_, s_ );
      ELAI_PROD( tmp2, s2_, s2_ );

      omega = tmp1 / tmp2;
      x = x + alpha * p_ + omega * s1_;

      r_ = s_ - omega * s2_;
      r2_ = r_;

      // Preconditioning:
      ELAI_PROF_BEG( prec_elapsed_ );
      P_->forward( r2_ );
      P_->backward( r2_ );
      ELAI_SYNC( r2_ );
      ELAI_PROF_END( prec_elapsed_ );

      ELAI_PROD( tmp1, rs0_, r2_ );
      ELAI_PROD( tmp2, rs0_, r1_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      beta = alpha / omega * tmp1 / tmp2;
      p_ = r2_ + beta * ( p_ - omega * p1_ );
      r1_ = r2_;

      res = fix_norm( r_ );
    }

    return false;
  }

public:
  bicgstab
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
    , r_( A_.m() ), r1_( A_.m() ), r2_( A_.m() ), rs0_( A_.m() )
    , p_( A_.m() ), p1_( A_.m() ), Ap_( A_.m() ), s_( A_.m() ), s1_( A_.m() ), s2_( A_.m() )
    , bthres_( static_cast< Coef >( 1e-16 ) )
  {
    iter_max( A_.m() );
  }
  ~bicgstab() {}

  Coef brk_thres() const { return bthres_; }
  Coef brk_thres( Coef thres )
  {
    Coef old = bthres_;

    bthres_ = thres;

    return old;
  }

  size_t mem() const
  {
    size_t sum = ksp< Coef >::mem();

    sum += r_.mem();
    sum += r1_.mem();
    sum += r2_.mem();
    sum += rs0_.mem();
    sum += p_.mem();
    sum += p1_.mem();
    sum += Ap_.mem();
    sum += s_.mem();
    sum += s1_.mem();
    sum += s2_.mem();
    sum += sizeof( bthres_ );

    return sum;
  }
};

}

#endif//__ELAI_BICGSTAB__

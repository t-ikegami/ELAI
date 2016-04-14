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
#ifndef __ELAI_BICGSAFE__
#define __ELAI_BICGSAFE__

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
class bicgsafe : public ksp< Coef >
{
  ELAI_USE_KSP;

private:
  // WORKSPACE
  vector< Coef > r_, r1_, rs0_, v_, u_, Au_, p_, Ap_, z_, y_, w_;
  Coef bthres_;

  bool solve_( vector< Coef >& x )
  {
    Coef res, res0, alpha, beta, eta, zeta, tmp1, tmp2, tmp3, tmp4, tmp5;

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

    v_ = A_ * r_;
    ELAI_SYNC( v_ );

    Ap_ = v_;

    ELAI_PROD( tmp1, rs0_, r_ );
    ELAI_PROD( tmp2, rs0_, Ap_ );

    alpha = tmp1 / tmp2;

    ELAI_PROD( tmp1, v_, r_ );
    ELAI_PROD( tmp2, v_, v_ );

    zeta = tmp1 / tmp2;
    eta = static_cast< Coef >( 0. );
    beta = static_cast< Coef >( 0. );
    u_ = r_; y_ = r_;

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

      w_ = zeta * Ap_ + eta * y_;
      u_ = w_ + eta * beta * u_;

      Au_ = A_ * u_;
      ELAI_SYNC( Au_ );

      z_ = zeta * r_ + eta * z_ - alpha * u_;
      y_ = zeta * v_ + eta * y_ - alpha * Au_;
      x = x + alpha * p_ + z_;
      r1_ = r_ - alpha * Ap_ - y_;

      ELAI_PROD( tmp1, rs0_, r1_ );
      ELAI_PROD( tmp2, rs0_, r_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      beta = ( alpha / zeta ) * ( tmp1 / tmp2 );
      r_ = r1_;
      p_ = r1_ + beta * ( p_ - u_ );

      v_ = A_ * r1_;
      ELAI_SYNC( v_ );

      Ap_ = v_ + beta * ( Ap_ - Au_ );

      ELAI_PROD( tmp1, rs0_, r_ );
      ELAI_PROD( tmp2, rs0_, Ap_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      alpha = tmp1 / tmp2;

      ELAI_PROD( tmp1, y_, y_ );
      ELAI_PROD( tmp2, v_, v_ );
      ELAI_PROD( tmp3, y_, v_ );
      ELAI_PROD( tmp4, y_, r_ );
      ELAI_PROD( tmp5, v_, r_ );

      zeta = ( tmp1 * tmp5 - tmp4 * tmp3 ) / ( tmp1 * tmp2 - tmp3 * tmp3 );
      eta = ( tmp2 * tmp4 - tmp3 * tmp5 ) / ( tmp1 * tmp2 - tmp3 * tmp3 );

      res = fix_norm( r_ );
    }

    return false;
  }

  bool solveP_( vector< Coef >& x )
  {
    Coef res, res0, alpha, beta, eta, zeta, tmp1, tmp2, tmp3, tmp4, tmp5;

    if ( is_trivial( b_ ) )
    {
      x = b_;

      return true;
    }

    r_ = b_ - A_ * x;
    res = res0 = sync_norm( r_ );
    if ( is_trivial( res0 ) ) return true;

    rs0_ = r_;
    r1_ = r_;

    // Preconditioning:
    ELAI_PROF_BEG( prec_elapsed_ );
    P_->forward( r1_ );
    P_->backward( r1_ );
    ELAI_SYNC( r1_ );
    ELAI_PROF_END( prec_elapsed_ );

    p_ = r1_;

    v_ = A_ * r1_;
    ELAI_SYNC( v_ );

    Ap_ = v_;

    ELAI_PROD( tmp1, rs0_, r_ );
    ELAI_PROD( tmp2, rs0_, Ap_ );

    alpha = tmp1 / tmp2;

    ELAI_PROD( tmp1, v_, r_ );
    ELAI_PROD( tmp2, v_, v_ );

    zeta = tmp1 / tmp2;
    eta = static_cast< Coef >( 0. );
    beta = static_cast< Coef >( 0. );
    u_ = r_; y_ = r_;

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

      w_ = zeta * Ap_ + eta * y_;

      // Preconditioning:
      ELAI_PROF_BEG( prec_elapsed_ );
      P_->forward( w_ );
      P_->backward( w_ );
      ELAI_SYNC( w_ );
      ELAI_PROF_END( prec_elapsed_ );

      u_ = w_ + eta * beta * u_;

      Au_ = A_ * u_;
      ELAI_SYNC( Au_ );

      z_ = zeta * r1_ + eta * z_ - alpha * u_;
      y_ = zeta * v_ + eta * y_ - alpha * Au_;
      x = x + alpha * p_ + z_;
      r1_ = r_ - alpha * Ap_ - y_;

      ELAI_PROD( tmp1, rs0_, r1_ );
      ELAI_PROD( tmp2, rs0_, r_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      beta = ( alpha / zeta ) * ( tmp1 / tmp2 );

      r_ = r1_;

      // Preconditioning:
      ELAI_PROF_BEG( prec_elapsed_ );
      P_->forward( r1_ );
      P_->backward( r1_ );
      ELAI_SYNC( r1_ );
      ELAI_PROF_END( prec_elapsed_ );

      p_ = r1_ + beta * ( p_ - u_ );

      v_ = A_ * r1_;
      ELAI_SYNC( v_ );

      Ap_ = v_ + beta * ( Ap_ - Au_ );

      ELAI_PROD( tmp1, rs0_, r_ );
      ELAI_PROD( tmp2, rs0_, Ap_ );

      // This avoidance for numerical breakdowns is NOT good.
      // However we have not implemented any other strategies.
      if ( sqrt( fabs( tmp2 ) ) / res0 <= bthres_ ) tmp2 = rs0_ * rs0_;

      alpha = tmp1 / tmp2;

      ELAI_PROD( tmp1, y_, y_ );
      ELAI_PROD( tmp2, v_, v_ );
      ELAI_PROD( tmp3, y_, v_ );
      ELAI_PROD( tmp4, y_, r_ );
      ELAI_PROD( tmp5, v_, r_ );

      zeta = ( tmp1 * tmp5 - tmp4 * tmp3 ) / ( tmp1 * tmp2 - tmp3 * tmp3 );
      eta = ( tmp2 * tmp4 - tmp3 * tmp5 ) / ( tmp1 * tmp2 - tmp3 * tmp3 );

      res = fix_norm( r_ );
    }

    return false;
  }

public:
  bicgsafe
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
    , r_( A_.m() ), r1_( A_.m() ), rs0_( A_.m() ), v_( A_.m() ), u_( A_.m() ), Au_( A_.m() )
    , p_( A_.m() ), Ap_( A_.m() ), z_( A_.m() ), y_( A_.m() ), w_( A_.m() )
    , bthres_( static_cast< Coef >( 1e-16 ) )
  {
    iter_max( A_.m() );
  }
  ~bicgsafe() {}

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
    sum += rs0_.mem();
    sum += v_.mem();
    sum += u_.mem();
    sum += Au_.mem();
    sum += p_.mem();
    sum += Ap_.mem();
    sum += z_.mem();
    sum += y_.mem();
    sum += w_.mem();
    sum += sizeof( bthres_ );

    return sum;
  }
};

}

#endif//__ELAI_BICGSAFE__

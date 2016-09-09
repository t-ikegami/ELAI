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
#ifndef __ELAI_EXPRESSION__
#define __ELAI_EXPRESSION__

#include <complex>

namespace elai
{

template< class Lhs, class Op, class Rhs >
class expression
{
  const Lhs& lhs_;
  const Rhs& rhs_;

public:
  typedef typename Op::range range;

  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}

  int m() const { return lhs_.m(); }
  int n() const { return lhs_.n(); }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()() const { return Op::apply( lhs_, rhs_ ); }
  range operator()( int i ) const { return Op::apply( lhs_( i ), rhs_( i ) ); }
  range operator()( int i, int j ) const { return Op::apply( lhs_( i, j ), rhs_( i, j ) ); }
};

template< class Coef >
struct expression_add
{
  typedef Coef range;
  static range apply( range lhs, range rhs ) { return lhs + rhs; }
};

template< class Coef >
struct expression_sub
{
  typedef Coef range;
  static range apply( range lhs, range rhs ) { return lhs - rhs; }
};

template< class Coef >
struct expression_mul
{
  typedef Coef range;
  static range apply( range lhs, range rhs ) { return lhs * rhs; }
};

template< class Coef >
struct expression_hat
{
  typedef Coef range;
  static range apply( range lhs, range rhs ) { return lhs * rhs; }
};

//
// operator+
//
template< class Lhs, class Rhs >
expression< Lhs, expression_add< typename Lhs::range >, Rhs >
operator+( const Lhs& lhs, const Rhs& rhs )
{ return expression< Lhs, expression_add< typename Lhs::range >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< typename Rhs::range, expression_add< typename Rhs::range>, Rhs >
operator+( const typename Rhs::range& lhs, const Rhs& rhs )
{ return expression< typename Rhs::range, expression_add< typename Rhs::range >, Rhs >(lhs, rhs); }


//
// operator-
//
template< class Lhs, class Rhs >
expression< Lhs, expression_sub< typename Lhs::range >, Rhs >
operator-( const Lhs& lhs, const Rhs& rhs )
{ return expression< Lhs, expression_sub< typename Lhs::range >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< typename Rhs::range, expression_sub< typename Rhs::range>, Rhs >
operator-( const typename Rhs::range& lhs, const Rhs& rhs )
{ return expression< typename Rhs::range, expression_sub< typename Rhs::range >, Rhs >(lhs, rhs); }


//
// operator*
//
template< class Lhs, class Rhs >
expression< Lhs, expression_mul< typename Lhs::range >, Rhs >
operator*( const Lhs& lhs, const Rhs& rhs )
{ return expression< Lhs, expression_mul< typename Lhs::range >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< typename Rhs::range, expression_mul< typename Rhs::range>, Rhs >
operator*( const typename Rhs::range& lhs, const Rhs& rhs )
{ return expression< typename Rhs::range, expression_mul< typename Rhs::range >, Rhs >(lhs, rhs); }


/*
//
// operator^
//
template< class Lhs, class Rhs >
expression< Lhs, expression_hat< typename Lhs::range >, Rhs >
operator^( const Lhs& lhs, const Rhs& rhs )
{ return expression< Lhs, expression_hat< typename Lhs::range >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< typename Rhs::range, expression_hat< typename Rhs::range>, Rhs >
operator^( const typename Rhs::range& lhs, const Rhs& rhs )
{ return expression< typename Rhs::range, expression_hat< typename Rhs::range >, Rhs >(lhs, rhs); }

*/

}

#endif//__ELAI_EXPRESSION__

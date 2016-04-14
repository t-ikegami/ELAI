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
expression< float, expression_add< float >, Rhs >
operator+( const float& lhs, const Rhs& rhs )
{ return expression< float, expression_add< float >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< double, expression_add< double >, Rhs >
operator+( const double& lhs, const Rhs& rhs )
{ return expression< double, expression_add< double >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< float >, expression_add< std::complex< float > >, Rhs >
operator+( const std::complex< float >& lhs, const Rhs& rhs )
{ return expression< std::complex< float >, expression_add< std::complex< float > >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< double >, expression_add< std::complex< double > >, Rhs >
operator+( const std::complex< double >& lhs, const Rhs& rhs )
{ return expression< std::complex< double >, expression_add< std::complex< double > >, Rhs >( lhs, rhs );
}

//
// operator-
//
template< class Lhs, class Rhs >
expression< Lhs, expression_sub< typename Lhs::range >, Rhs >
operator-( const Lhs& lhs, const Rhs& rhs )
{ return expression< Lhs, expression_sub< typename Lhs::range >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< float, expression_sub< float >, Rhs >
operator-( const float& lhs, const Rhs& rhs )
{ return expression< float, expression_sub< float >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< double, expression_sub< double >, Rhs >
operator-( const double& lhs, const Rhs& rhs )
{ return expression< double, expression_sub< double >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< float >, expression_sub< std::complex< float > >, Rhs >
operator-( const std::complex< float >& lhs, const Rhs& rhs )
{ return expression< std::complex< float >, expression_sub< std::complex< float > >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< double >, expression_sub< std::complex< double > >, Rhs >
operator-( const std::complex< double >& lhs, const Rhs& rhs )
{ return expression< std::complex< double >, expression_sub< std::complex< double > >, Rhs >( lhs, rhs );
}

//
// operator*
//
template< class Lhs, class Rhs >
expression< Lhs, expression_mul< typename Lhs::range >, Rhs >
operator*( const Lhs& lhs, const Rhs& rhs )
{ return expression< Lhs, expression_mul< typename Lhs::range >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< float, expression_mul< float >, Rhs >
operator*( const float& lhs, const Rhs& rhs )
{ return expression< float, expression_mul< float >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< double, expression_mul< double >, Rhs >
operator*( const double& lhs, const Rhs& rhs )
{ return expression< double, expression_mul< double >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< float >, expression_mul< std::complex< float > >, Rhs >
operator*( const std::complex< float >& lhs, const Rhs& rhs )
{ return expression< std::complex< float >, expression_mul< std::complex< float > >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< double >, expression_mul< std::complex< double > >, Rhs >
operator*( const std::complex< double >& lhs, const Rhs& rhs )
{ return expression< std::complex< double >, expression_mul< std::complex< double > >, Rhs >( lhs, rhs );
}

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
expression< float, expression_hat< float >, Rhs >
operator^( const float& lhs, const Rhs& rhs )
{ return expression< float, expression_hat< float >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< double, expression_hat< double >, Rhs >
operator^( const double& lhs, const Rhs& rhs )
{ return expression< double, expression_hat< double >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< float >, expression_hat< std::complex< float > >, Rhs >
operator^( const std::complex< float >& lhs, const Rhs& rhs )
{ return expression< std::complex< float >, expression_hat< std::complex< float > >, Rhs >( lhs, rhs );
}

template< class Rhs >
expression< std::complex< double >, expression_hat< std::complex< double > >, Rhs >
operator^( const std::complex< double >& lhs, const Rhs& rhs )
{ return expression< std::complex< double >, expression_hat< std::complex< double > >, Rhs >( lhs, rhs );
}
*/

}

#endif//__ELAI_EXPRESSION__

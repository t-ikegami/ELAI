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
#ifndef __ELAI_BLAS__
#define __ELAI_BLAS__

#include "expression.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace elai
{

template< class Op, class Rhs >
class expression< typename Rhs::range, Op, Rhs >
{
  typedef typename Rhs::range Lhs;

  const Lhs lhs_;
  const Rhs& rhs_;

public:
  typedef Lhs range;

  expression( const Lhs lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}

  int m() const { return rhs_.m(); }
  int n() const { return rhs_.n(); }
  int nnz() const { return rhs_.nnz(); }
  int ind( int i ) const { return rhs_.ind( i ); }
  int col( int k ) const { return rhs_.col( k ); }
  range operator()() const { return Op::apply( lhs_, rhs_ ); }
  range operator()( int i ) const { return Op::apply( lhs_, rhs_( i ) ); }
  range operator()( int i, int j ) const { return Op::apply( lhs_, rhs_( i, j ) ); }
};

template< class Lhs, class Op >
class expression< Lhs, Op, typename Lhs::range >
{
  typedef typename Lhs::range Rhs;

  const Lhs& lhs_;
  const Rhs rhs_;

public:
  typedef Rhs range;

  expression( const Lhs& lhs, const Rhs rhs ) : lhs_( lhs ), rhs_( rhs ) {}

  int m() const { return lhs_.m(); }
  int n() const { return lhs_.n(); }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()() const { return Op::apply( lhs_, rhs_ ); }
  range operator()( int i ) const { return Op::apply( lhs_(i), rhs_ ); }
  range operator()( int i, int j ) const { return Op::apply( lhs_(i, j), rhs_ ); }
};

// vector( i ) + scalar -> vector( i )
template< class Coef >
class expression
  < vector< Coef >
  , expression_add< Coef >
  , Coef
  >
{
  typedef vector< Coef > Lhs;
  typedef Coef Rhs;
  const Lhs& lhs_;
  const Rhs rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return lhs_.m(); }
  int n() const { return 1; }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i ) const { return lhs_( i ) + rhs_; }
};

// vector( i ) + vector( i ) -> vector( i )
template< class Coef >
class expression
  < vector< Coef >
  , expression_add< Coef >
  , vector< Coef >
  >
{
  typedef vector< Coef > Lhs;
  typedef vector< Coef> Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { assert( lhs_.m() == rhs_.m() ); return lhs_.m(); }
  int n() const { return 1; }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i ) const { return lhs_( i ) + rhs_( i ); }
};

// vector( i ) - scalar -> vector( i )
template< class Coef >
class expression
  < vector< Coef >
  , expression_sub< Coef >
  , Coef
  >
{
  typedef vector< Coef > Lhs;
  typedef Coef Rhs;
  const Lhs& lhs_;
  const Rhs rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return lhs_.m(); }
  int n() const { return 1; }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i ) const { return lhs_( i ) - rhs_; }
};

// vector( i ) - vector( i ) -> vector( i )
template< class Coef >
class expression
  < vector< Coef >
  , expression_sub< Coef >
  , vector< Coef >
  >
{
  typedef vector< Coef > Lhs;
  typedef vector< Coef> Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { assert( lhs_.m() == rhs_.m() ); return lhs_.m(); }
  int n() const { return 1; }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i ) const { return lhs_( i ) - rhs_( i ); }
};

// vector( i ) * scalar -> vector( i )
template< class Coef >
class expression
  < vector< Coef >
  , expression_mul< Coef >
  , Coef
  >
{
  typedef vector< Coef > Lhs;
  typedef Coef Rhs;
  const Lhs& lhs_;
  const Rhs rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return rhs_.m(); }
  int n() const { return 1; }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i ) const { return lhs_( i ) * rhs_; }
};

// vector( i ) * vector( i ) -> scalar( R )
template< class Coef >
class expression
  < vector< Coef >
  , expression_mul< Coef >
  , vector< Coef >
  >
{
  typedef vector< Coef > Lhs;
  typedef vector< Coef> Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return 1; }
  int n() const { return 1; }
  int nnz() const { return 1; }
  int ind( int i ) const { return 0; }
  int col( int k ) const { return 0; }
  template< class Return >
  operator Return() const { return static_cast< Return >( ( *this )() ); }
  range operator()() const
  {
    assert( lhs_.m() == rhs_.m() );
    range acc = static_cast< range >( 0 );
    for ( int i = 0; i < lhs_.m(); ++i ) acc += lhs_( i ) * rhs_( i );
    return acc;
  }
};

// vector( i ) * vector( i ) -> scalar( C )
template< class Coef >
class expression
  < vector< std::complex< Coef > >
  , expression_mul< std::complex< Coef > >
  , vector< std::complex< Coef > >
  >
{
  typedef vector< std::complex< Coef > > Lhs;
  typedef vector< std::complex< Coef > > Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef std::complex< Coef > range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return 1; }
  int n() const { return 1; }
  int nnz() const { return 1; }
  int ind( int i ) const { return 0; }
  int col( int k ) const { return 0; }
  template< class Return >
  operator Return() const { return static_cast< Return >( ( *this )() ); }
  range operator()() const
  {
    assert( lhs_.m() == rhs_.m() );
    range acc = static_cast< range >( 0 );
    for ( int i = 0; i < lhs_.m(); ++i ) acc += lhs_( i ) * std::conj( rhs_( i ) );
    return acc;
  }
};

/*
// vector( i ) ^ vector( j ) -> matrix( j, i )
template< class Coef >
class expression
  < vector< Coef >
  , expression_hat< Coef >
  , vector< Coef >
  >
{
  typedef vector< Coef > Lhs;
  typedef vector< Coef > Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return lhs_.m(); }
  int n() const { return rhs_.m(); }
  int nnz() const { return rhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return rhs_.col( k ); }
  range operator()( int i, int j ) const { return lhs_( j ) * rhs_( i ); }
};
*/

// matrix( i, j ) * scalar -> matrix( i, j )
template< class Coef >
class expression
  < matrix< Coef >
  , expression_mul< Coef >
  , Coef
  >
{
  typedef matrix< Coef > Lhs;
  typedef Coef Rhs;
  const Lhs& lhs_;
  const Rhs rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return lhs_.m(); }
  int n() const { return lhs_.n(); }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i, int j ) const { return lhs_( i, j ) * rhs_; }
};

// matrix( i, j ) + matrix( i, j ) -> matrix( i, j )
template< class Coef >
class expression
  < matrix< Coef >
  , expression_add< Coef >
  , matrix< Coef >
  >
{
  typedef matrix< Coef > Lhs;
  typedef matrix< Coef > Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return lhs_.m(); }
  int n() const { return lhs_.n(); }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i, int j ) const { return lhs_( i, j ) + rhs_( i, j ); }
};

// matrix( i, j ) - matrix( i, j ) -> matrix( i, j )
template< class Coef >
class expression
  < matrix< Coef >
  , expression_sub< Coef >
  , matrix< Coef >
  >
{
  typedef matrix< Coef > Lhs;
  typedef matrix< Coef > Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return lhs_.m(); }
  int n() const { return lhs_.n(); }
  int nnz() const { return lhs_.nnz(); }
  int ind( int i ) const { return lhs_.ind( i ); }
  int col( int k ) const { return lhs_.col( k ); }
  range operator()( int i, int j ) const { return lhs_( i, j ) - rhs_( i, j ); }
};

// matrix( i, j ) * vector( j ) -> vector( i )
template< class Coef >
class expression
  < matrix< Coef >
  , expression_mul< Coef >
  , vector< Coef >
  >
{
  typedef matrix< Coef > Lhs;
  typedef vector< Coef > Rhs;
  const Lhs& lhs_;
  const Rhs& rhs_;
public:
  typedef Coef range;
  expression( const Lhs& lhs, const Rhs& rhs ) : lhs_( lhs ), rhs_( rhs ) {}
  int m() const { return rhs_.m(); }
  int n() const { return 1; }
  int nnz() const { return rhs_.nnz(); }
  int ind( int i ) const { return rhs_.ind( i ); }
  int col( int k ) const { return rhs_.col( k ); }
  range operator()( int i ) const
  {
    range acc = static_cast< range >( 0 );

    //for ( int j = 0; j < lhs_.n(); ++j ) acc += lhs_( i, j ) * rhs_( j );
    for ( int k = lhs_.ind( i ); k < lhs_.ind( i + 1 ); ++k ) acc += lhs_.val( k ) * rhs_( lhs_.col( k ) );

    return acc;
  }
};

}

#endif//__ELAI_BLAS__

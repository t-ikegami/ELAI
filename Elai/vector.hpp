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
#ifndef __ELAI_VECTOR__
#define __ELAI_VECTOR__

#include <iostream>
#include <sstream>
#include <string>
#include "def.hpp"
#include "expression.hpp"

namespace elai
{

template< class Coef >
class vector
{
  int m_;
  Coef *f_;
  size_t mem_;

  void terminate()
  {
    m_ = 0;
    if ( f_ != NULL ) { delete [] f_; f_ = NULL; }
    mem_ = 0;
  }

  void init()
  {
    int padding = CACHE_LINE / sizeof( Coef );
    int dmy = ( padding - m_ % padding ) % padding;

    f_ = new Coef[ m_ + dmy ];

    for ( int i = 0; i < m_ + dmy; ++i ) f_[ i ] = static_cast< Coef >( 0 );
    mem_ = sizeof( Coef ) * ( m_ + dmy );
  }

  template< class Lhs, class Op, class Rhs >
  void eval( const expression< Lhs, Op, Rhs >& expr )
  {
#ifdef ELAI_USE_OPENMP
    #pragma omp parallel for
#endif
    for ( int i = 0; i < m_; ++i ) f_[ i ] = expr( i );
  }

public:
  typedef Coef range;

  vector() : m_( 0 ), f_( NULL ), mem_( 0 ) {}
  vector( const int m ) : m_( m ), f_( NULL ), mem_( 0 ) { init(); }
  vector( const int m, range *f ) : m_( m ), f_( NULL ), mem_( 0 )
  {
    init();
    for ( int i = 0; i < m_; ++i ) f_[ i ] = f[ i ];
  }
  vector( const vector< range >& src ) : m_( src.m_ ), f_( NULL ), mem_( 0 )
  {
    init();
    for ( int i = 0; i < m_; ++i ) f_[ i ] = src.f_[ i ];
  }
  vector( std::istream& is ) : m_( 0 ), f_( NULL ), mem_( 0 )
  {
    std::string header;

    std::getline( is, header );
    if ( header.find( "array" ) != std::string::npos )
    {
      int m, n = 0;

      while ( n == 0 )
      {
        std::string buf;

        std::getline( is, buf );
        if ( buf.find( "%" ) != std::string::npos ) continue;

        std::stringstream ss( buf );

        ss >> m >> n;
      }
      if ( n != 1 )
      {
        std::cerr << "YOUR MM-FILE CONTAINS MATRIX DATA." << std::endl;
        std::abort();
      }
      else
      {
        m_ = m;
        init();
        for ( int i = 0; i < m_; ++i ) is >> f_[ i ];
      }
    }
    else if ( header.find( "coordinate" ) != std::string::npos )
    {
      int m, n, nnz = 0;

      while ( nnz == 0 )
      {
        std::string buf;

        std::getline( is, buf );
        if ( buf.find( "%" ) != std::string::npos ) continue;

        std::stringstream ss( buf );

        ss >> m >> n >> nnz;
      }
      if ( n != 1 || m != nnz )
      {
        std::cerr << "YOUR MM-FILE CONTAINS MATRIX DATA." << std::endl;
        std::abort();
      }
      else
      {
        m_ = m;
        init();
        for ( int i = 0; i < m_; ++i ) is >> n >> m >> f_[ i ];
      }
    }
    else
    {
      std::cerr << "ELAI COULD NOT READ YOUR MM-FILE." << std::endl;
      std::abort();
    }
  }
  template< class Lhs, class Op, class Rhs >
  vector( const expression< Lhs, Op, Rhs >& expr ) : m_( expr.m() ), f_( NULL ), mem_( 0 )
  {
    init();
    eval( expr );
  }
  ~vector()
  {
    terminate();
  }

  void setup( int m )
  {
    terminate();
    m_ = m;
    init();
  }

  vector< range >& operator=( const range v )
  {
#ifdef ELAI_USE_OPENMP
    #pragma omp parallel for
#endif
    for ( int i = 0; i < m_; ++i ) f_[ i ] = v;
    return *this;
  }
  vector< range >& operator=( const vector< range >& rhs )
  {
    assert( m_ == rhs.m_ );
#ifdef ELAI_USE_OPENMP
    #pragma omp parallel for
#endif
    for ( int i = 0; i < m_; ++i ) f_[ i ] = rhs.f_[ i ];
    return *this;
  }
  template< class Lhs, class Op, class Rhs >
  vector< range >& operator=( const expression< Lhs, Op, Rhs >& expr )
  {
    eval( expr );
    return *this;
  }

  int m() const { return m_; }
  int n() const { return 1; }
  int nnz() const { return m_; }

  // FOR EXPRESSIONS, DO NOT TOUCH!
  int ind( int i ) const { return i; }
  int col( int k ) const { return 0; }
  range *val() { return f_; }
  const range *val() const { return f_; }

  range& operator()( int i )
  {
    assert( i < m_ );
    return f_[ i ];
  }
  range& operator()( int i ) const
  {
    assert( i < m_ );
    return f_[ i ];
  }

  vector< range >& clear( const range v )
  {
    for ( int i = 0; i < m_; ++i ) f_[ i ] = v;

    return *this;
  }

  vector< range >& reorder( const int *perm )
  {
    vector< range > tmp( *this );

    for ( int i = 0; i < m_; ++i ) f_[ i ] = tmp( perm[ i ] );

    return *this;
  }

  vector< range >& scale( const vector< range >& scal )
  {
    for ( int i = 0; i < m_; ++i ) f_[ i ] *= scal( i );

    return *this;
  }
  vector< range >& unscale( const vector< range >& scal )
  {
    for ( int i = 0; i < m_; ++i ) f_[ i ] /= scal( i );

    return *this;
  }

  vector< range >& operator<<( std::istream& is )
  {
    std::string header;

    terminate();
    std::getline( is, header );
    if ( header.find( "array" ) != std::string::npos )
    {
      int m, n = 0;

      while ( n == 0 )
      {
        std::string buf;

        std::getline( is, buf );
        if ( buf.find( "%" ) != std::string::npos ) continue;

        std::stringstream ss( buf );

        ss >> m >> n;
      }
      if ( m_ != m || 1 != n ) std::cerr << "YOUR MM-FILE HAS UNMATCHED SIZE." << std::endl;
      else for ( int i = 0; i < m_; ++i ) is >> f_[ i ];
    }
    else if ( header.find( "coordinate" ) != std::string::npos )
    {
      int m, n, nnz = 0;

      while ( nnz == 0 )
      {
        std::string buf;

        std::getline( is, buf );
        if ( buf.find( "%" ) != std::string::npos ) continue;

        std::stringstream ss( buf );

        ss >> m >> n >> nnz;
      }
      if ( m_ != m || 1 != n || m_ != nnz ) std::cerr << "YOUR MM-FILE HAS UNMATCHED SIZE." << std::endl;
      else for ( int i = 0; i < m_; ++i ) is >> m >> n >> f_[ i ];
    }
    else std::cerr << "ELAI COULD NOT READ YOUR MM-FILE." << std::endl;
    return *this;
  }
  std::ostream& operator>>( std::ostream& os ) const
  {
    os << "%%MatrixMarket matrix array real general" << std::endl;
    os << m_ << " 1" << std::endl;
    for ( int i = 0; i < m_; ++i ) os << "  " << f_[ i ] << std::endl;
    return os;
  }

  size_t mem() const { return mem_; }
};

template< class Coef >
std::ostream& operator<<( std::ostream& os, const vector< Coef >& v )
{
  v >> os;
  return os;
}

}

#endif//__ELAI_VECTOR__

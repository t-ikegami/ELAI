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
#ifndef __ELAI_MATRIX__
#define __ELAI_MATRIX__

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include "def.hpp"
#include "coherence.hpp"
#include "vector.hpp"
#include "expression.hpp"

namespace elai
{

template< class Coef >
class matrix
{
public:
  typedef Coef range;

private:
  int m_, n_, nnz_;
  int *ind_, *col_;
  range *c_, z_;
  size_t mem_;

  vector< range > scalR_;
  vector< range > scalC_;

#ifdef ELAI_USE_MPI
  inline void sync( coherence *coherent, vector< range >& u )
  {
    if ( coherent != NULL ) ( *coherent )( u.val() );
  }
#endif

  void init()
  {
    int padding, dmy;

    mem_ = 0;
    padding = CACHE_LINE / sizeof( int );
    dmy = ( padding - ( m_ + 1 ) % padding ) % padding;
    ind_ = new int[ ( m_ + 1 ) + dmy ];
    mem_ += sizeof( int ) * ( ( m_ + 1 ) + dmy );
    for ( int i = 0; i <= m_ + dmy; ++i ) ind_[ i ] = 0;
    dmy = ( padding - nnz_ % padding ) % padding;
    col_ = new int[ nnz_ + dmy ];
    mem_ += sizeof( int ) * ( nnz_ + dmy );
    for ( int i = 0; i < nnz_ + dmy; ++i ) col_[ i ] = 0;

    padding = CACHE_LINE / sizeof( range );
    dmy = ( padding - nnz_ % padding ) % padding;
    c_ = new range[ nnz_ + dmy ];
    mem_ += sizeof( range ) * ( nnz_ + dmy );
    z_ = static_cast< range >( 0 );
    for ( int i = 0; i < nnz_ + dmy; ++i ) c_[ i ] = z_;

    scalR_.setup( m_ );
    scalC_.setup( n_ );
  }

  void terminate()
  {
    if ( c_ != NULL ) { delete [] c_; c_ = NULL; }
    if ( col_ != NULL ) { delete [] col_; col_ = NULL; }
    if ( ind_ != NULL ) { delete [] ind_; ind_ = NULL; }
  }

  template< class Lhs, class Op, class Rhs >
  void eval( const expression< Lhs, Op, Rhs >& expr )
  {
    int k = 0;
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        c_[ k ] = expr( i, col_[ k ] );
        ++k;
      }
    }
  }

public:
#ifdef ELAI_USE_PYTHON
  void init(int m, int n, int nnz)
  {
    m_ = m;
    n_ = n;
    nnz_ = nnz;
    init();
  }
#endif

  matrix()
    : m_( 0 ), n_( 0 ), nnz_( 0 )
    , ind_( NULL ), col_( NULL ), c_( NULL ), z_( 0 )
    , mem_( 0 ), scalR_(), scalC_() {}
  matrix( int m, int n, int nnz, int *ind, int *col, range *c = NULL )
    : m_( m ), n_( n ), nnz_( nnz )
    , ind_( NULL ), col_( NULL ), c_( NULL ), z_( 0 )
    , mem_( 0 ), scalR_( m_ ), scalC_( n_ )
  {
    init();
    for ( int i = 0; i <= m_; ++i ) ind_[ i ] = ind[ i ];
    for ( int i = 0; i < nnz_; ++i ) col_[ i ] = col[ i ];
    if ( c != NULL ) for ( int i = 0; i < nnz_; ++i ) c_[ i ] = c[ i ];
    for ( int i = 0; i < m_; ++i ) scalR_( i ) = static_cast< range >( 1. );
    for ( int j = 0; j < n_; ++j ) scalC_( j ) = static_cast< range >( 1. );
  }
  matrix( const matrix< range >& src )
    : m_( src.m_ ), n_( src.n_ ), nnz_( src.nnz_ )
    , ind_( NULL ), col_( NULL ), c_( NULL ), z_( 0 )
    , mem_( src.mem_ ), scalR_( src.scalR_ ), scalC_( src.scalC_ )
  {
    init();
    for ( int i = 0; i <= m_; ++i ) ind_[ i ] = src.ind_[ i ];
    for ( int i = 0; i < nnz_; ++i ) col_[ i ] = src.col_[ i ];
    for ( int i = 0; i < nnz_; ++i ) c_[ i ] = src.c_[ i ];
  }
  matrix( std::istream& is )
    : m_( 0 ), n_( 0 ), nnz_( 0 )
    , ind_( NULL ), col_( NULL ), c_( NULL ), z_( 0 )
    , mem_( 0 ), scalR_(), scalC_()
  {
    std::string header;

    std::getline( is, header );
    if ( header.find( "coordinate" ) != std::string::npos )
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
      m_ = m; n_ = n;; nnz_ = nnz;
      init();
      int *cnt, *row, *col;
      range *val;
      cnt = new int[ m_ ];
      for ( int i = 0; i < m_; ++i ) cnt[ i ] = 0;
      row = new int[ nnz_ ];
      col = new int[ nnz_ ];
      val = new range[ nnz_ ];
      for ( int i = 0; i < nnz_; ++i )
      {
        is >> row[ i ] >> col[ i ] >> val[ i ];
        row[ i ] -= 1;
        col[ i ] -= 1;
        cnt[ row[ i ] ] += 1;
      }
      for ( int i = 0; i < m_; ++i ) ind_[ i + 1 ] = cnt[ i ] + ind_[ i ];
      for ( int i = 0; i < m_; ++i ) cnt[ i ] = 0;
      for ( int k = 0; k < nnz_; ++k )
      {
        // ROW/COL WISE IS NOT A MATTER, BUT ASCENDING ORDER IS NECESSARY.
        int i = row[ k ], j = col[ k ];
        col_[ ind_[ i ] + cnt[ i ] ] = j;
        c_[ ind_[ i ] + cnt[ i ] ] = val[ k ];
        cnt[ i ] += 1;
      }
      delete [] val;
      delete [] col;
      delete [] row;
      delete [] cnt;
    }
    else
    {
      std::cerr << "ELAI ALLOWS ONLY THE COORDINATE STYLE FOR MATRICES. SORRY." << std::endl;
      std::abort();
    }
  }
  template< class Lhs, class Op, class Rhs >
  matrix( const expression< Lhs, Op, Rhs >& expr ) : m_( expr.m() ), n_( expr.n() ), nnz_( expr.nnz() )
  {
    init();
    for ( int i = m_; 0 <= i; --i )
    {
      ind_[ i ] = expr.ind( i );
      if ( i < m_ ) for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) col_[ k ] = expr.col( k );
    }
    eval( expr );
  }
  ~matrix()
  {
    terminate();
  }

  matrix< range >& setup( int m, int n, int nnz, int *ind, int *col, range *c = NULL )
  {
    terminate();
    m_ = m; n_ = n; nnz_ = nnz;
    init();
    for ( int i = 0; i <= m_; ++i ) ind_[ i ] = ind[ i ];
    for ( int i = 0; i < nnz_; ++i ) col_[ i ] = col[ i ];
    if ( c != NULL ) for ( int i = 0; i < nnz_; ++i ) c_[ i ] = c[ i ];

    return *this;
  }

  matrix< range >& operator=( const range c )
  {
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) c_[ k ] = c;
    }
    return *this;
  }
  matrix< range >& operator=( const matrix< range >& rhs )
  {
    assert( m_ == rhs.m_ && n_ == rhs.n_ && nnz_ == rhs.nnz_ );
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) c_[ k ] = rhs( i, col_[ k ] );
    }
    return *this;
  }
  template< class Lhs, class Op, class Rhs >
  matrix< range >& operator=( const expression< Lhs, Op, Rhs >& expr )
  {
    // expr must be the same non-zero structure!
    eval( expr );
    return *this;
  }

  int m() const { return m_; }
  int n() const { return n_; }
  int nnz() const { return nnz_; }

  // FOR EXPRESSIONS, DO NOT TOUCH!
  int ind( int i ) const { return ind_[ i ]; }
  int col( int k ) const { return col_[ k ]; }
  range val( int k ) const { return c_[ k ]; }

#ifdef ELAI_USE_PYTHON
  void getInd(int *nelems, int **data)
  {
    *nelems = m_ + 1;
    *data = ind_;
  }

  void getCol(int *nelems, int **data)
  {
    *nelems = nnz_;
    *data = col_;
  }

  void getC(int *nelems, range **data)
  {
    *nelems = nnz_;
    *data = c_;
  }
#endif // ELAI_USE_PYTHON

  // FOR ONLY MUMPS, OTHERS DO NOT TOUCH!!
  int *ind() { return ind_; }
  const int *ind() const { return ind_; }
  int *col() { return col_; }
  const int *col() const { return col_; }
  range *val() { return c_; }
  const range *val() const { return c_; }

  range& operator()( int i, int j )
  {
    //std::cerr << i << " < " << m_ << ", " << j << " < " << n_ << std::endl;
    assert( 0 <= i && i < m_ && 0 <= j && j < n_ );
    for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
    {
      if ( col_[ k ] == j ) return c_[ k ];
    }
    //std::cerr << "HIT THE ZERO-REG" << std::endl;
    return z_;
  }
  const range& operator()( int i, int j ) const
  {
    //std::cerr << i << " < " << m_ << ", " << j << " < " << n_ << std::endl;
    assert( 0 <= i && i < m_ && 0 <= j && j < n_ );
    for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
    {
      if ( col_[ k ] == j ) return c_[ k ];
    }
    //std::cerr << "HIT THE ZERO-REG" << std::endl;
    return z_;
  }

  matrix< range >& clear( const range c )
  {
    for ( int k = 0; k < nnz_; ++k ) c_[ k ] = c;

    return *this;
  }

  bool is_symmetric( bool is_numeric = false ) const
  {
    bool is_symm = true;
    int *kind = new int[ m_ ];

    for ( int i = 0; i < m_; ++i ) kind[ i ] = ind_[ i ];
    for ( int i = 0; i < m_; ++i )
    {
      for ( int ki = ind_[ i ]; ki < ind_[ i + 1 ]; ++ki )
      {
        int j = col_[ ki ];
        int kj = kind[ j ];

        if ( j <= i ) continue;
        else if ( col_[ kj ] != i )
        {
          is_symm = false;

          return is_symm;
        }
        else if ( is_numeric )
        {
          if ( c_[ kj ] != c_[ ki ] )
          {
            is_symm = false;

            return is_symm;
          }
          else ++kj; // Numerical Matched
        }
        else ++kj; // Structural Matched

        kind[ j ] = kj;
      }
    }

    delete [] kind;

    return is_symm;
  }

  matrix< range >& transpose()
  {
    int *k0 = new int[ n_ ];
    int *ind = new int[ n_ + 1 ];
    int *col = new int[ nnz_ ];
    range *val = new range[ nnz_ ];

    for ( int i = 0; i < n_; ++i ) k0[ i ] = 0;
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) k0[ col_[ k ] ] += 1;
    ind[ 0 ] = 0;
    for ( int j = 1; j <= n_; ++j )
    {
      ind[ j ] = ind[ j - 1 ] + k0[ j - 1 ];
      k0[ j - 1 ] = 0;
    }
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        int j = col_[ k ];
        int k1 = ind[ j ] + k0[ j ];

        col[ k1 ] = i;
        val[ k1 ] = c_[ k ];
        k0[ j ] += 1;
      }
    }
    std::swap( m_, n_ );
    std::swap( ind, ind_ );
    std::swap( col, col_ );
    std::swap( val, c_ );
    delete [] val;
    delete [] col;
    delete [] ind;
    delete [] k0;

    return *this;
  }

  matrix< range >& conj()
  {
    int *k0 = new int[ n_ ];
    int *ind = new int[ n_ + 1 ];
    int *col = new int[ nnz_ ];
    range *val = new range[ nnz_ ];

    for ( int i = 0; i < n_; ++i ) k0[ i ] = 0;
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) k0[ col_[ k ] ] += 1;
    ind[ 0 ] = 0;
    for ( int j = 1; j <= n_; ++j )
    {
      ind[ j ] = ind[ j - 1 ] + k0[ j - 1 ];
      k0[ j - 1 ] = 0;
    }
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        int j = col_[ k ];
        int k1 = ind[ j ] + k0[ j ];

        col[ k1 ] = i;
        val[ k1 ] = conj_( c_[ k ] );
        k0[ j ] += 1;
      }
    }
    std::swap( m_, n_ );
    std::swap( ind, ind_ );
    std::swap( col, col_ );
    std::swap( val, c_ );
    delete [] val;
    delete [] col;
    delete [] ind;
    delete [] k0;

    return *this;
  }

  matrix< range >& perm_row( const int *perm )
  {
    int *ind = new int[ m_ + 1 ];
    int *col = new int[ nnz_ ];
    range *val = new range[ nnz_ ];

    ind[ 0 ] = 0;
    for ( int i = 0; i < m_; ++i ) ind[ i + 1 ] = ind_[ perm[ i ] + 1 ] - ind_[ perm[ i ] ];
    for ( int i = 0; i < m_; ++i ) ind[ i + 1 ] = ind[ i + 1 ] + ind[ i ];
    for ( int i = 0; i < m_; ++i )
    {
      int k1 = ind_[ perm[ i ] ];
      for ( int k = ind[ i ]; k < ind[ i + 1 ]; ++k )
      {
        col[ k ] = col_[ k1 ];
        val[ k ] = c_[ k1 ];
        ++k1;
      }
    }
    std::swap( ind, ind_ );
    std::swap( col, col_ );
    std::swap( val, c_ );
    delete [] val;
    delete [] col;
    delete [] ind;

    return *this;
  }

  matrix< range >& reorder( const int *perm )
  { // A' <- P A P^t = P ( P A^t )^t, P; perm-row
    transpose();
    perm_row( perm );
    transpose();
    perm_row( perm );

    return *this;
  }

  matrix< range >& normalize
    ( range thres = static_cast< range >( 1e-03 )
#ifdef ELAI_USE_MPI
    , coherence *coherent = NULL
#endif
    )
  {
    vector< range > diagR( m_ );
    vector< range > diagC( m_ );
    range *c = new range[ nnz_ ];

    for ( int i = 0; i < m_; ++i )
    {
      diagR( i ) = static_cast< range >( 0. );
      diagC( i ) = static_cast< range >( 0. );
      scalR_( i ) = static_cast< range >( 1. );
      scalC_( i ) = static_cast< range >( 1. );
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) c[ k ] = c_[ k ];
    }
    for ( int n = 0; n < m_; ++n )
    {
      range maxR = static_cast< range >( 0. );
      range maxC = static_cast< range >( 0. );

      for ( int i = 0; i < m_; ++i )
      {
        for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
        {
          int j = col_[ k ];
          range v = fabs( c[ k ] );

          if ( diagR( i ) < v ) diagR( i ) = v;
          if ( diagC( j ) < v ) diagC( j ) = v;
        }
      }
#ifdef ELAI_USE_MPI
      sync( coherent, diagR );
      sync( coherent, diagC );
#endif
      for ( int i = 0; i < m_; ++i )
      {
        // This iteration is capable of parallel executions excepted maxR.
        range tmpR = fabs( static_cast< range >( 1. ) - diagR( i ) );
        range tmpC = fabs( static_cast< range >( 1. ) - diagR( i ) );

        if ( maxR < tmpR ) maxR = tmpR;
        scalR_( i ) /= sqrt( diagR( i ) );
        diagR( i ) = static_cast< range >( 0. );
        if ( maxC < tmpC ) maxC = tmpC;
        scalC_( i ) /= sqrt( diagC( i ) );
        diagC( i ) = static_cast< range >( 0. );
      }
#ifdef ELAI_USE_MPI
      if ( coherent != NULL )
      {
        MPI_Allreduce( MPI_IN_PLACE, &maxR, 1, mpi_< range >().type, MPI_MAX, coherent->comm() );
        MPI_Allreduce( MPI_IN_PLACE, &maxC, 1, mpi_< range >().type, MPI_MAX, coherent->comm() );
      }
#endif
      for ( int i = 0; i < m_; ++i )
      {
        for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
        {
          int j = col_[ k ];

          c[ k ] = scalR_( i ) * c_[ k ] * scalC_( j );
        }
      }
      if ( maxR <= thres && maxC <= thres ) break;
      //if ( fabs( 1. - maxR ) <= thres && fabs( 1. - maxC ) <= thres ) break;
    }
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
        c_[ k ] = c[ k ];

    delete [] c;

    return *this;
  }
  matrix< range >& normalizeRow()
  {
    for ( int i = 0; i < m_; ++i )
    {
      scalR_( i ) = static_cast< range >( 0 );

      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        range s = fabs( c_[ k ] );

        if ( scalR_( i ) < s ) scalR_( i ) = s;
      }
      if ( 0 < scalR_( i ) ) scalR_( i ) = static_cast< range >( 1. ) / scalR_( i );
      else scalR_( i ) = static_cast< range >( 1. );
      for ( int k = ind_[ i ]; k != ind_[ i + 1 ]; ++k ) c_[ k ] *= scalR_( i );
    }

    return *this;
  }
  matrix< range >& normalizeCol()
  {
    for ( int j = 0; j < n_; ++j ) scalC_( j ) = static_cast< range >( 0 );
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        int j = col_[ k ];
        range s = fabs( c_[ k ] );

        if ( scalC_( j ) < s ) scalC_( j ) = s;
      }
    }
    for ( int j = 0; j < n_; ++j )
      if ( 0 < scalC_( j ) ) scalC_( j ) = static_cast< range >( 1. ) / scalC_( j );
      else scalC_( j ) = static_cast< range >( 1. );
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) c_[ k ] *= scalC_( col_[ k ] );

    return *this;
  }
  matrix< range >& unnormalize()
  {
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        int j = col_[ k ];
        
        c_[ k ] /= scalR_( i ) * scalC_( j );
      }

    return *this;
  }
  matrix< range >& unnormalizeRow()
  {
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) c_[ k ] /= scalR_( i );

    return *this;
  }
  matrix< range >& unnormalizeCol()
  {
    for ( int i = 0; i < m_; ++i )
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k ) c_[ k ] /= scalC_( col_[ k ] );

    return *this;
  }
  const vector< range >& scaleRow() const { return scalR_; }
  const vector< range >& scaleCol() const { return scalC_; }
  const range scaleRowNorm() const
  {
    range rnorm = static_cast< range >( 0. );

    for ( int i = 0; i < m_; ++i ) rnorm += scalR_( i ) * scalR_( i );

    return sqrt( rnorm );
  }
  const range scaleColNorm() const
  {
    range cnorm = static_cast< range >( 0. );

    for ( int i = 0; i < m_; ++i ) cnorm += scalC_( i ) * scalC_( i );

    return sqrt( cnorm );
  }
  const range scaleRatio() const
  {
    range rnorm = static_cast< range >( 0. );
    range cnorm = static_cast< range >( 0. );

    for ( int i = 0; i < m_; ++i ) rnorm += scalR_( i ) * scalR_( i );
    for ( int i = 0; i < m_; ++i ) cnorm += static_cast< range >( 1. ) / ( scalC_( i ) * scalC_( i ) );

    return sqrt( rnorm ) / sqrt( cnorm );
  }

  matrix< range >& operator<<( std::istream& is )
  {
    std::string header;

    std::getline( is, header );
    if ( header.find( "coordinate" ) != std::string::npos )
    {
      int m, n, nnz;

      while ( nnz == 0 )
      {
        std::string buf;

        std::getline( is, buf );
        if ( buf.find( "%" ) != std::string::npos ) continue;

        std::stringstream ss( buf );

        ss >> m >> n >> nnz;
      }
      if ( m_ != m || n_ != n || nnz_ != nnz )
      {
        std::cerr << "YOUR MM-FILE HAS UNMATCHED SIZE FOR THE MATRIX." << std::endl;
        std::abort();
      }
      for ( int k = 0; k < nnz_; ++k )
      {
        int i, j; range v;
        is >> i >> j >> v;
        ( *this )( i - 1, j - 1 ) = v;
      }
    }
    else
    {
      std::cerr << "ELAI ALLOWS ONLY THE COORDINATE STYLE FOR MATRICES. SORRY." << std::endl;
      std::abort();
    }

    return *this;
  }
  std::ostream& operator>>( std::ostream& os ) const
  {
    os << "%%MatrixMarket matrix coordinate real general" << std::endl; 
    os << m_ << " " << n_ << "  " << nnz_ << std::endl;
    for ( int i = 0; i < m_; ++i )
    {
      for ( int k = ind_[ i ]; k < ind_[ i + 1 ]; ++k )
      {
        os << " " << ( i + 1 ) << " " << ( col_[ k ] + 1 ) << "  " << c_[ k ] << std::endl;
      }
    }

    return os;
  }

  size_t mem() const { return mem_; }
};

template< class range >
std::ostream& operator<<( std::ostream& os, const matrix< range >& m )
{
  m >> os;

  return os;
}

}

#endif//__ELAI_MATRIX__
